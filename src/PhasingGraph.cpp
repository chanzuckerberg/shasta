// Shasta.
#include "PhasingGraph.hpp"
#include "AssemblyGraph2.hpp"
#include "deduplicate.hpp"
#include "diploidBayesianPhase.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "algorithm.hpp"
#include <cstdlib>
#include <queue>
#include "tuple.hpp"



void PhasingGraphEdge::runBayesianModel(double epsilon, bool allowRandomHypothesis)
{
    tie(logPin, logPout) = diploidBayesianPhase(matrix, epsilon);

    if(allowRandomHypothesis) {

        // Used for bubble removal.
        if(logPin >= logPout) {
            relativePhase = 0;
            logP = min(logPin-logPout, logPin);
        } else {
            relativePhase = 1;
            logP = min(logPout-logPin, logPout);
        }

    } else {

        // Used for phasing.
        logP = std::abs(logPin - logPout);
        if(logPin >= logPout) {
            relativePhase = 0;
        } else {
            relativePhase = 1;
        }

    }
}



PhasingGraph::PhasingGraph(
    const AssemblyGraph2& assemblyGraph2,
    uint64_t minConcordantReadCount,
    uint64_t maxDiscordantReadCount,
    double minLogP,
    double epsilon,
    size_t threadCount,
    bool allowRandomHypothesis) :
    MultithreadedObject<PhasingGraph>(*this)
{
    createVertices(assemblyGraph2);
    createOrientedReadsTable(assemblyGraph2.getReadCount());
    createEdges(
        minConcordantReadCount, maxDiscordantReadCount,
        minLogP, epsilon, threadCount, allowRandomHypothesis);
}



void PhasingGraph::createVertices(const AssemblyGraph2& assemblyGraph2)
{
    PhasingGraph& phasingGraph = *this;

    // Loop over assembly graph edges that describe diploid bubbles.
    BGL_FORALL_EDGES(e, assemblyGraph2, AssemblyGraph2) {
        const AssemblyGraph2Edge& edge = assemblyGraph2[e];

        // If not diploid, skip.
        if(edge.ploidy() != 2) {
            continue;
        }

        // If not phased, skip.
        if(not edge.isPhased()) {
            continue;
        }

        // Get the component this bubble belongs to.
        const uint64_t componentId = edge.componentId;

        // Get the vertex corresponding to this component, creating it if necessary.
        const vertex_descriptor v = getVertex(componentId);

        // Add this bubble to this vertex.
        phasingGraph[v].bubbles.push_back(make_pair(e, edge.phase));
    }


    // For each vertex of the phasing graph, find the oriented reads on each side.
    array< vector<OrientedReadId>, 2> orientedReadIds;
    BGL_FORALL_VERTICES(v, phasingGraph, PhasingGraph) {
        PhasingGraphVertex& vertex = phasingGraph[v];

        // Loop over the bubbles of this vertex.
        orientedReadIds[0].clear();
        orientedReadIds[1].clear();
        for(const auto& p: vertex.bubbles) {
            const AssemblyGraph2::edge_descriptor e = p.first;
            const uint64_t phase = p.second;
            const AssemblyGraph2Edge& edge = assemblyGraph2[e];
            SHASTA_ASSERT(edge.ploidy() == 2);

            // Loop over the two sides of the bubble.
            for(uint64_t bubbleSide=0; bubbleSide<2; bubbleSide++) {
                const AssemblyGraph2Edge::Branch& branch = edge.branches[bubbleSide];

                // Find the corresponding side on this vertex.
                uint64_t vertexSide = ((phase==0) ? bubbleSide : 1-bubbleSide);

                // Copy the oriented read ids of the bubble to the appropriate
                // side of the phasing graph vertex.
                copy(branch.orientedReadIds.begin(), branch.orientedReadIds.end(),
                    back_inserter(orientedReadIds[vertexSide]));
            }
        }

        // Sort and remove duplicates.
        for(uint64_t vertexSide=0; vertexSide<2; vertexSide++) {
            deduplicate(orientedReadIds[vertexSide]);
        }

        // Store the oriented read ids in the phasing graph vertex.
        // But don't store oriented read ids that appear on both sides.
        auto begin0 = orientedReadIds[0].begin();
        auto begin1 = orientedReadIds[1].begin();
        auto end0 = orientedReadIds[0].end();
        auto end1 = orientedReadIds[1].end();
        auto it0 = begin0;
        auto it1 = begin1;
        while(true) {

            if(it0 == end0) {
                for(; it1!=end1; ++it1) {
                    vertex.orientedReadIds[1].push_back(*it1);
                }
                break;
            }

            if(it1 == end1) {
                for(; it0!=end0; ++it0) {
                    vertex.orientedReadIds[0].push_back(*it0);
                }
                break;
            }

            if(*it0 < *it1) {
                vertex.orientedReadIds[0].push_back(*it0);
                ++it0;
            } else if(*it1 < *it0) {
                vertex.orientedReadIds[1].push_back(*it1);
                ++it1;
            } else {
                // This oriented read appears on both sides.
                // Don't store it.
                ++it0;
                ++it1;
            }

        }
    }
}



void PhasingGraph::createEdges(
    uint64_t minConcordantReadCount,
    uint64_t maxDiscordantReadCount,
    double minLogP,
    double epsilon,
    size_t threadCount,
    bool allowRandomHypothesis)
{
    performanceLog << timestamp << "AssemblyGraph2::PhasingGraph::createEdges begins." << endl;
    PhasingGraph& phasingGraph = *this;

    // Create a vector of all vertices, to be processed
    // one by one in parallel.
    createEdgesData.allVertices.clear();
    BGL_FORALL_VERTICES(v, phasingGraph, PhasingGraph) {
        createEdgesData.allVertices.push_back(v);
    }

    // Store the parameters so all threads can see them.
    createEdgesData.minConcordantReadCount = minConcordantReadCount;
    createEdgesData.maxDiscordantReadCount = maxDiscordantReadCount;
    createEdgesData.minLogP = minLogP;
    createEdgesData.epsilon = epsilon;
    createEdgesData.allowRandomHypothesis = allowRandomHypothesis;

    // Process all vertices in parallel.
    const uint64_t batchSize = 100;
    setupLoadBalancing(createEdgesData.allVertices.size(), batchSize);
    runThreads(&PhasingGraph::createEdgesThreadFunction, threadCount);

    performanceLog << timestamp << "AssemblyGraph2::PhasingGraph::createEdges ends." << endl;
}



void PhasingGraph::createEdgesThreadFunction(size_t threadId)
{
    PhasingGraph& phasingGraph = *this;

    const uint64_t minConcordantReadCount = createEdgesData.minConcordantReadCount;
    const uint64_t maxDiscordantReadCount = createEdgesData.maxDiscordantReadCount;
    const double minLogP = createEdgesData.minLogP;
    const double epsilon = createEdgesData.epsilon;
    const bool allowRandomHypothesis = createEdgesData.allowRandomHypothesis;
    vector<CreateEdgesData::EdgeData> edgeData;

    // Temporary storage of the edges found by this thread.
    vector< tuple<vertex_descriptor, vertex_descriptor, PhasingGraphEdge> > threadEdges;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all vertices assigned to this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            createEdges(createEdgesData.allVertices[i],
                minConcordantReadCount,
                maxDiscordantReadCount,
                minLogP, epsilon, edgeData, threadEdges,
                allowRandomHypothesis);
        }
    }

    std::lock_guard<std::mutex> lock(mutex);
    for(const auto& t: threadEdges) {
        const vertex_descriptor v0 = get<0>(t);
        const vertex_descriptor v1 = get<1>(t);
        const PhasingGraphEdge&  edge = get<2>(t);
        add_edge(v0, v1, edge, phasingGraph);
    }
}



// Create edges between vertex vA and vertices vB with id greater
// that the id of vA.
void PhasingGraph::createEdges(
    PhasingGraph::vertex_descriptor vA,
    uint64_t minConcordantReadCount,
    uint64_t maxDiscordantReadCount,
    double minLogP,
    double epsilon,
    vector<CreateEdgesData::EdgeData>& edgeData,
    vector< tuple<vertex_descriptor, vertex_descriptor, PhasingGraphEdge> >& threadEdges,
    bool allowRandomHypothesis)
{
    PhasingGraph& phasingGraph = *this;

    const PhasingGraphVertex& vertexA = phasingGraph[vA];

    // Gather EdgeData for this vA.
    edgeData.clear();
    for(uint64_t sideA=0; sideA<2; sideA++) {
        for(const OrientedReadId orientedReadId: vertexA.orientedReadIds[sideA]) {

            // Find all places where this read occurs.
            const auto& v = orientedReadsTable[orientedReadId.getValue()];

            for(const auto& p: v) {
                const vertex_descriptor vB = p.first;

                // Don't add it twice.
                if(vB <= vA) {
                    continue;
                }

                const uint64_t sideB = p.second;

                edgeData.push_back({vB, sideA, sideB});
            }
        }

    }

    // Sort the EdgeData by vB.
    sort(edgeData.begin(), edgeData.end());

    // Each streak with the same vB generates an edge, if there
    // is a sufficient number of reads.
    for(auto it=edgeData.begin(); it!=edgeData.end();  /* Increment later */) {

        auto streakBegin = it;
        auto streakEnd = streakBegin;
        const PhasingGraph::vertex_descriptor vB = streakBegin->vB;
        while((streakEnd != edgeData.end()) and (streakEnd->vB == vB)) {
            ++streakEnd;
        }

        // If this streak is sufficiently long, it can generate an edge.
        const uint64_t streakLength = uint64_t(streakEnd - streakBegin);
        if(streakLength >= minConcordantReadCount) {
            PhasingGraphEdge edge;
            for(auto jt=streakBegin; jt!=streakEnd; ++jt) {
                ++edge.matrix[jt->sideA][jt->sideB];
            }

            if( (edge.concordantCount() >= minConcordantReadCount) and
                (edge.discordantCount() <= maxDiscordantReadCount)) {

                edge.runBayesianModel(epsilon, allowRandomHypothesis);

                if(edge.logP > minLogP) {
                    threadEdges.push_back(make_tuple(vA, vB, edge));
                }
            }
        }

        // Prepare to process the next streak.
        it = streakEnd;
    }
}



void PhasingGraph::createOrientedReadsTable(uint64_t readCount)
{
    PhasingGraph& phasingGraph = *this;

    orientedReadsTable.clear();
    orientedReadsTable.resize(readCount * 2);
    BGL_FORALL_VERTICES(v, phasingGraph, PhasingGraph) {
        const PhasingGraphVertex& vertex = phasingGraph[v];
        for(uint64_t side=0; side<2; side++) {
            for(const OrientedReadId orientedReadId: vertex.orientedReadIds[side]) {
                orientedReadsTable[orientedReadId.getValue()].push_back(make_pair(v, side));
            }
        }
    }

}


// Find the optimal spanning tree using logFisher as the edge weight.
// Edges that are part of the optimal spanning tree get their
// isTreeEdge set.
void PhasingGraph::computeSpanningTree()
{
    PhasingGraph& phasingGraph = *this;

    // Create a vector of edges sorted by decreasing logFisher.
    vector<pair<PhasingGraph::edge_descriptor, double> > edgeTable;
    BGL_FORALL_EDGES(e, phasingGraph, PhasingGraph) {
        edgeTable.push_back(make_pair(e, phasingGraph[e].logP));
    }
    sort(edgeTable.begin(), edgeTable.end(),
        OrderPairsBySecondOnlyGreater<PhasingGraph::edge_descriptor, double>());


    // Computation of connected components and optimal spanning tree.

    // Initialize the disjoint sets data structure.
    const uint64_t n = num_vertices(phasingGraph);
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint32_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    // Process edges in order of decreasing logP.
    uint64_t treeEdgeCount = 0;
    for(const auto& p: edgeTable) {
        const PhasingGraph::edge_descriptor e = p.first;
        const PhasingGraph::vertex_descriptor v0 = source(e, phasingGraph);
        const PhasingGraph::vertex_descriptor v1 = target(e, phasingGraph);
        if(disjointSets.find_set(v0) != disjointSets.find_set(v1)) {
            disjointSets.union_set(v0, v1);
            phasingGraph[e].isTreeEdge = true;
            ++treeEdgeCount;
        }
    }
    cout << "Found " << treeEdgeCount << " edges of the optimal spanning tree." << endl;

}



// Phase vertices using the spanning tree.
void PhasingGraph::phase()
{
    PhasingGraph& phasingGraph = *this;

    // Loop over start vertices. One iteration for each connected component.
    uint64_t componentId = 0;
    BGL_FORALL_VERTICES(vStart, phasingGraph, PhasingGraph) {
        PhasingGraphVertex& vertexStart = phasingGraph[vStart];

        // If this vertex has already been assigned to a component,
        // don't use it as a start vertex.
        if(vertexStart.componentId != PhasingGraphVertex::invalidComponentId) {
            continue;
        }


        // BFS on the optimal spanning tree, starting at this vertex.
        std::queue<PhasingGraph::vertex_descriptor> q;
        q.push(vStart);
        vertexStart.componentId = componentId;
        vertexStart.phase = 0;
        while(not q.empty()) {

            // Dequeue a vertex.
            const PhasingGraph::vertex_descriptor v0 = q.front();
            q.pop();
            PhasingGraphVertex& vertex0 = phasingGraph[v0];
            SHASTA_ASSERT(vertex0.componentId == componentId);
            const uint64_t phase0 = vertex0.phase;

            // Loop over tree edges incident to this vertex.
            BGL_FORALL_OUTEDGES(v0, e, phasingGraph, PhasingGraph) {
                const PhasingGraphEdge& edge = phasingGraph[e];
                if(not edge.isTreeEdge) {
                    continue;
                }

                // If we already encountered the other vertex, skip.
                const PhasingGraph::vertex_descriptor v1 = target(e, phasingGraph);
                PhasingGraphVertex& vertex1 = phasingGraph[v1];
                if(vertex1.componentId != PhasingGraphVertex::invalidComponentId) {
                    SHASTA_ASSERT(vertex1.componentId == componentId);
                    continue;
                }

                // Add the other vertex to this component and phase it.
                q.push(v1);
                vertex1.componentId = componentId;

                if(edge.relativePhase == 0) {
                    vertex1.phase = phase0;
                } else {
                    vertex1.phase = 1 - phase0;
                }
            }
        }

        // Prepare to handle the next connected component.
        ++componentId;
    }

}



// Get the vertex corresponding to a componentId, creating it if necessary.
PhasingGraph::vertex_descriptor PhasingGraph::getVertex(uint64_t componentId)
{
    PhasingGraph& phasingGraph = *this;

    // Make sure the vertex exists.
    while(componentId >= num_vertices(phasingGraph)) {
        add_vertex(phasingGraph);
    }

    // The vertex descriptor is the same as the componentId.
    return componentId;
}



void PhasingGraph::writeEdgesCsv(
    const string& fileName,
    const AssemblyGraph2& assemblyGraph2) const
{
    const PhasingGraph& phasingGraph = *this;

    ofstream csv(fileName);
    csv << "v0,v1,m00,m01,m10,m11,logPin,logPout,logP,RelativePhase\n";

    BGL_FORALL_EDGES(e, phasingGraph,PhasingGraph) {
        const PhasingGraphEdge& edge = phasingGraph[e];
        csv << source(e, phasingGraph) << ",";
        csv << target(e, phasingGraph) << ",";
        csv << edge.matrix[0][0] << ",";
        csv << edge.matrix[0][1] << ",";
        csv << edge.matrix[1][0] << ",";
        csv << edge.matrix[1][1] << ",";
        csv << edge.logPin << ",";
        csv << edge.logPout << ",";
        csv << edge.logP << ",";
        csv << edge.relativePhase << "\n";
    }

}





// Store the phasing in the AssemblyGraph2.
void PhasingGraph::storePhasing(AssemblyGraph2& assemblyGraph2) const
{
    const PhasingGraph& phasingGraph = *this;

    // Remove all existing phasing information from the AssemblyGraph2.
    BGL_FORALL_EDGES(e, assemblyGraph2, AssemblyGraph2) {
        AssemblyGraph2Edge& edge = assemblyGraph2[e];
        edge.componentId = AssemblyGraph2Edge::invalidComponentId;
        edge.phase = AssemblyGraph2Edge::invalidPhase;
    }

    // Copy phasing information from the PhasingGraph to the AssemblyGraph2.
    BGL_FORALL_VERTICES(v, phasingGraph, PhasingGraph) {
        const PhasingGraphVertex& vertex = phasingGraph[v];
        SHASTA_ASSERT(vertex.isPhased());

        for(const auto& p: vertex.bubbles) {
            const AssemblyGraph2::edge_descriptor e = p.first;
            const uint64_t bubblePhase = p.second;

            AssemblyGraph2Edge& edge = assemblyGraph2[e];
            SHASTA_ASSERT(edge.ploidy() == 2);
            edge.componentId = vertex.componentId;

            edge.phase = vertex.phase;
            if(bubblePhase == 1) {
                edge.phase = 1 - edge.phase;
            }
        }
    }

}



void PhasingGraph::writeCsv(
    const string& baseName,
    const AssemblyGraph2& assemblyGraph2) const
{
    writeVerticesCsv(baseName + "-Vertices.csv");
    writeVerticesDetailsCsv(baseName + "-Vertices-Details.csv", assemblyGraph2);
    writeEdgesCsv(baseName + "-Edges.csv", assemblyGraph2);
}



void PhasingGraph::writeVerticesCsv(
    const string& fileName) const
{
    const PhasingGraph& phasingGraph = *this;

    ofstream csv(fileName);
    csv << "v,Bubbles,ComponentId,Phase\n";

    BGL_FORALL_VERTICES(v, phasingGraph,PhasingGraph) {
        const PhasingGraphVertex& vertex = phasingGraph[v];
        csv << v << ",";
        csv << vertex.bubbles.size() << ",";
        csv << vertex.componentId << ",";
        csv << vertex.phase << "\n";
    }

}



void PhasingGraph::writeVerticesDetailsCsv(
    const string& fileName,
    const AssemblyGraph2& assemblyGraph2) const
{
    const PhasingGraph& phasingGraph = *this;

    ofstream csv(fileName);
    csv << "v,Bubble,Phase\n";

    BGL_FORALL_VERTICES(v, phasingGraph,PhasingGraph) {
        const PhasingGraphVertex& vertex = phasingGraph[v];
        for(const auto& p: vertex.bubbles) {
            const AssemblyGraph2::edge_descriptor e = p.first;
            const uint64_t phase = p.second;
            csv << v << ",";
            csv << assemblyGraph2[e].id << ",";
            csv << phase << "\n";
        }
    }

}



void PhasingGraph::writeGraphviz(const string& fileName) const
{
    const PhasingGraph& phasingGraph = *this;

    ofstream out(fileName);
    out << "graph PhasingGraph {\n";



    BGL_FORALL_VERTICES(v, phasingGraph,PhasingGraph) {
        const PhasingGraphVertex& vertex = phasingGraph[v];
        out << v << " [tooltip=\"" << v << ", " << vertex.bubbles.size() << " bubbles\"];\n";
    }



    BGL_FORALL_EDGES(e, phasingGraph,PhasingGraph) {
        const PhasingGraphEdge& edge = phasingGraph[e];

        const PhasingGraph::vertex_descriptor v0 = source(e, phasingGraph);
        const PhasingGraph::vertex_descriptor v1 = target(e, phasingGraph);
        const uint64_t phase0 = phasingGraph[v0].phase;
        const uint64_t phase1 = phasingGraph[v1].phase;

        string color;
        if(edge.isTreeEdge) {
            color = "black";
        } else {
            if(edge.logPin >= edge.logPout) {
                if(phase0 == phase1) {
                    color = "green";
                } else {
                    color = "red";
                }
            } else {
                if(phase0 == phase1) {
                    color = "red";
                } else {
                    color = "green";
                }
            }
        }

        out << v0 << "--" << v1 <<
            " [tooltip=\"" << v0 << " " << v1 << " " << edge.logP <<
            "\" color=\"" << color << "\"];\n";
    }

    out << "}\n";
}


