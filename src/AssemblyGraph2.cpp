// Shasta.
#include "AssemblyGraph2.hpp"
#include "AssemblyGraph2Statistics.hpp"
#include "AssembledSegment.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "AssemblerOptions.hpp"
#include "copyNumber.hpp"
#include "deduplicate.hpp"
#include "dominatorTree.hpp"
#include "enumeratePaths.hpp"
#include "findLinearChains.hpp"
#include "findMarkerId.hpp"
#include "GfaAssemblyGraph.hpp"
#include "orderPairs.hpp"
#include "PhasingGraph.hpp"
#include "performanceLog.hpp"
#include "ReadFlags.hpp"
using namespace shasta;

// Boost libraries.
// boost/pending/disjoint_sets.hpp must be included first
// to avoid unexplained compilation problems.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/dominator_tree.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/reverse_graph.hpp>

// Standard library.
#include "fstream.hpp"
#include <limits>
#include <map>
#include <numeric>


// The constructor creates an edge for each linear path
// in the marker graph. Therefore, immediately after construction,
// each edge has a single MarkerGraphPath (no bubbles).
AssemblyGraph2::AssemblyGraph2(
    uint64_t readRepresentation,
    uint64_t k, // Marker length
    const MemoryMapped::Vector<ReadFlags>& readFlags,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph,
    uint64_t pruneLength,
    const Mode2AssemblyOptions& mode2Options,
    AssemblyGraph2Statistics& statistics,
    size_t threadCount
    ) :
    MultithreadedObject<AssemblyGraph2>(*this),
    readRepresentation(readRepresentation),
    k(k),
    readFlags(readFlags),
    markers(markers),
    markerGraph(markerGraph)
{



    // Threshold that defines a strong branch.
    // A branch is strong if it is supported by at least this number of
    // distinct oriented reads.
    // Weak branches are subject to removal by removeWeakBranches
    // (but at least one branch in each bubble will always be kept).
    const uint64_t strongBranchThreshold = mode2Options.strongBranchThreshold;

    // Epsilon for the Bayesian model used for phasing and for bubble removal.
    // This is the probability that a read appears on the wrong branch.
    const double epsilon = mode2Options.epsilon;

    // Parameters for bubble removal.
    const uint64_t minConcordantReadCountForBubbleRemoval = mode2Options.minConcordantReadCountForBubbleRemoval;
    const uint64_t maxDiscordantReadCountForBubbleRemoval = mode2Options.maxDiscordantReadCountForBubbleRemoval;
    const double minLogPForBubbleRemoval = mode2Options.minLogPForBubbleRemoval;
    const uint64_t componentSizeThresholdForBubbleRemoval = mode2Options.componentSizeThresholdForBubbleRemoval;

    // Parameters for phasing.
    const uint64_t minConcordantReadCountForPhasing = mode2Options.minConcordantReadCountForPhasing;
    const uint64_t maxDiscordantReadCountForPhasing = mode2Options.maxDiscordantReadCountForPhasing;
    const double minLogPForPhasing = mode2Options.minLogPForPhasing;

    // Parameters for superbubble removal.
    const uint64_t maxSuperbubbleSize = mode2Options.maxSuperbubbleSize;
    const uint64_t maxSuperbubbleChunkSize = mode2Options.maxSuperbubbleChunkSize;
    const uint64_t maxSuperbubbleChunkPathCount = mode2Options.maxSuperbubbleChunkPathCount;
    const uint64_t superbubbleRemovalEdgeLengthThreshold = mode2Options.superbubbleEdgeLengthThreshold;



    performanceLog << timestamp << "AssemblyGraph2 constructor begins." << endl;

    // Because of the way we write the GFA file (without overlaps),
    // k is required to be even.
    SHASTA_ASSERT((k % 2) == 0);

    // Create the assembly graph and do some initial simple transformations.
    create();
    prune(pruneLength);
    removeShortLoopbackEdges(superbubbleRemovalEdgeLengthThreshold);

    // Gather parallel edges into bubbles.
    gatherBubbles();

    // Handle superbubbles.
    handleSuperbubbles0(superbubbleRemovalEdgeLengthThreshold,
        maxSuperbubbleSize, maxSuperbubbleChunkSize, maxSuperbubbleChunkPathCount, false, false);
    merge(false, false);
    handleSuperbubbles1(
        maxSuperbubbleSize, maxSuperbubbleChunkSize, maxSuperbubbleChunkPathCount, false, false);
    merge(false, false);

    // Store the reads supporting each branch of each edge.
    storeReadInformationParallel(threadCount);

    // Remove weak branches.
    removeWeakBranches(strongBranchThreshold);
    merge(true, false);
    gatherBubbles();
    forceMaximumPloidy(2);

    // Assemble sequence.
    assembleParallel(threadCount);

    // Remove degenerate edges (both branches have the same sequence).
    removeDegenerateBranches();
    merge(true, true);
    prune(pruneLength);

    // Use the PhasingGraph to iteratively remove bad bubbles, then to phase.
    removeBadBubblesIterative(
        minConcordantReadCountForBubbleRemoval,
        maxDiscordantReadCountForBubbleRemoval,
        minLogPForBubbleRemoval,
        epsilon,
        superbubbleRemovalEdgeLengthThreshold,
        maxSuperbubbleSize,
        maxSuperbubbleChunkSize,
        maxSuperbubbleChunkPathCount,
        pruneLength,
        componentSizeThresholdForBubbleRemoval,
        threadCount);
    hierarchicalPhase(
        minConcordantReadCountForPhasing,
        maxDiscordantReadCountForPhasing,
        minLogPForPhasing,
        epsilon,
        threadCount);

    // Final pruning.
    prune(pruneLength);

    // Find chains of bubbles.
    // These are linear chains of edges of length at least 2.
    findBubbleChains();
    writeBubbleChains();
    findPhasingRegions();
    writePhasingRegions();

    // Write out what we have.
    storeGfaSequence();
    if(not mode2Options.suppressDetailedOutput) {
        writeDetailed("Assembly-Detailed", true, false, true,
            not mode2Options.suppressGfaOutput, not mode2Options.suppressFastaOutput);
        if(not mode2Options.suppressGfaOutput) {
            writeDetailed("Assembly-Detailed-NoSequence", false, false, false, true, false);
        }
    }
    if(not mode2Options.suppressHaploidOutput) {
        writeHaploid("Assembly-Haploid", true, true,
            not mode2Options.suppressGfaOutput, not mode2Options.suppressFastaOutput, &statistics);
        if(not mode2Options.suppressGfaOutput) {
            writeHaploid("Assembly-Haploid-NoSequence", false, false, true, false);
        }
    }
    if(not mode2Options.suppressPhasedOutput) {
        writePhased("Assembly-Phased", true, true,
            not mode2Options.suppressGfaOutput, not mode2Options.suppressFastaOutput, &statistics);
        if(not mode2Options.suppressGfaOutput) {
            writePhased("Assembly-Phased-NoSequence", false, false, true, false);
        }
        writePhasedDetails();
    }

    // Het snp statistics.
    uint64_t transitionCount, transversionCount, nonSnpCount;
    hetSnpStatistics(transitionCount, transversionCount, nonSnpCount);
    const uint64_t snpCount = transitionCount + transversionCount;
    cout <<
        "The following SNP statistics only count SNPs that are well separated "
        "from other heterozygous variants. \n"
        "There are " << snpCount <<
        " heterozygous SNPs (" <<
        transitionCount << " transitions, " <<
        transversionCount << " transversions).\n" <<
        "Transition/transversion ratio is " <<
        double(transitionCount) / double(transversionCount) << "\n"
        "There are " << nonSnpCount << " small bubbles which are not SNPs." << endl;
    statistics.simpleSnpBubbleTransitionCount = transitionCount;
    statistics.simpleSnpBubbleTransversionCount = transversionCount;
    statistics.nonSimpleSnpBubbleCount = nonSnpCount;

    performanceLog << timestamp << "AssemblyGraph2 constructor ends." << endl;
}



// Initial creation of vertices and edges.
void AssemblyGraph2::create()
{
    performanceLog << timestamp << "AssemblyGraph2::create begins." << endl;

    const bool debug = false;
    ofstream debugOut;
    if(debug) {
        debugOut.open("AssemblyGraph2-Constructor.txt");
    }

    const MarkerGraph::EdgeId edgeCount = markerGraph.edges.size();
    vector<bool> wasFound(edgeCount, false);

    MarkerGraphPath nextEdges;
    MarkerGraphPath previousEdges;
    MarkerGraphPath path;
    MarkerGraphPath reverseComplementedPath;



    // Main loop over all edges of the marker graph.
    // At each iteration we find a new linear path of edges.
    for(MarkerGraph::EdgeId startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {
        if(markerGraph.edges[startEdgeId].wasRemoved()) {
            continue;
        }
        if(debug) {
            const MarkerGraph::Edge& startEdge = markerGraph.edges[startEdgeId];
            debugOut << "Starting a new path at edge " << startEdgeId << " " <<
                startEdge.source << "->" << startEdge.target << endl;
        }

        // If we already found this edge, skip it.
        // It is part of a path we already found.
        if(wasFound[startEdgeId]) {
            continue;
        }

        // Follow the path forward.
        nextEdges.clear();
        MarkerGraph::EdgeId edgeId = startEdgeId;
        bool isCircular = false;
        while(true) {
            const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
            const MarkerGraph::VertexId v1 = edge.target;
            // const auto outEdges = markerGraph.edgesBySource[v1];
            if(markerGraph.outDegree(v1) != 1) {
                break;
            }
            // const auto inEdges = markerGraph.edgesByTarget[v1];
            if(markerGraph.inDegree(v1) != 1) {
                break;
            }
            edgeId = markerGraph.getFirstNonRemovedOutEdge(v1);
            if(debug) {
                debugOut << "Forward " << edgeId << " " <<
                    edge.source << "->" << edge.target << endl;
            }
            if(edgeId == startEdgeId) {
                isCircular = true;
                if(debug) {
                    cout << "Found a circular edge." << endl;
                    debugOut << "Found a circular edge." << endl;
                }
                break;
            }
            nextEdges.push_back(edgeId);
            SHASTA_ASSERT(not wasFound[edgeId]);
        }

        // Follow the path backward.
        previousEdges.clear();
        if(!isCircular) {
            edgeId = startEdgeId;
            while(true) {
                const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
                const MarkerGraph::VertexId v0 = edge.source;
                // const auto outEdges = markerGraph.edgesBySource[v0];
                if(markerGraph.outDegree(v0) != 1) {
                    break;
                }
                // const auto inEdges = markerGraph.edgesByTarget[v0];
                if(markerGraph.inDegree(v0) != 1) {
                    break;
                }
                edgeId = markerGraph.getFirstNonRemovedInEdge(v0);
                if(debug) {
                    debugOut << "Backward " << edgeId << " " <<
                        edge.source << "->" << edge.target << endl;
                }
                previousEdges.push_back(edgeId);
                SHASTA_ASSERT(not wasFound[edgeId]);
            }
        }

        // Gather the path.
        path.clear();
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter(path));
        path.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter(path));

        if(debug) {
            for(const MarkerGraph::EdgeId edgeId: path) {
                const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
                debugOut << "Path " << edgeId << " " <<
                    edge.source << "->" << edge.target << endl;
            }
        }

        // Mark all the edges in the path as found.
        for(const MarkerGraph::EdgeId edgeId: path) {
            if(wasFound[edgeId]) {
                cout << "Assertion failed at " << edgeId << endl;
                if(debug) {
                    debugOut << "Assertion failed at " << edgeId << endl;
                }
                SHASTA_ASSERT(0);
            }
            wasFound[edgeId] = true;
        }



        // Check if this path originated from one of the read graph
        // components that we want to assemble
        // (one of each reverse complemented pair).
        // If not, don't store it.
        // This way we do a single-stranded assembly.
        {
            const MarkerGraph::Edge& edge = markerGraph.edges[path.front()];
            const MarkerGraphVertexId vertexId = edge.source;
            const span<const MarkerId> markerIds = markerGraph.getVertexMarkerIds(vertexId);
            const MarkerId markerId = markerIds[0];
            const OrientedReadId orientedReadId = findMarkerId(markerId, markers).first;
            const ReadId readId = orientedReadId.getReadId();
            const Strand strand = orientedReadId.getStrand();
            if(readFlags[readId].strand != strand) {
                continue;
            }
        }

        // See if this path contains secondary edges.
        bool containsSecondaryEdges = false;
        for(const MarkerGraph::EdgeId edgeId: path) {
            const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
            if(edge.isSecondary) {
                containsSecondaryEdges = true;
                break;
            }
        }

        // Store this path as a new edge of the assembly graph.
        addEdge(path, containsSecondaryEdges, false, false);
    }



    // Check that all edges of the marker graph were found.
    for(MarkerGraph::EdgeId edgeId=0; edgeId<edgeCount; edgeId++) {
        SHASTA_ASSERT(wasFound[edgeId] or markerGraph.edges[edgeId].wasRemoved());
    }

    performanceLog << timestamp << "AssemblyGraph2::create ends." << endl;
}



// Prune leaves in which all branches are shorter than the specified length.
#if 0
void AssemblyGraph2::prune(uint64_t pruneLength)
{
    G& g = *this;

    const bool debug = false;
    performanceLog << timestamp << "AssemblyGraph2::prune begins" << endl;

    while(true) {
        performanceLog << timestamp << "Prune iteration begins." << endl;
        vector<edge_descriptor> edgesToBeRemoved;
        BGL_FORALL_EDGES(e, g, G) {
            const vertex_descriptor v0 = source(e, g);
            const vertex_descriptor v1 = target(e, g);
            const E& edge = g[e];

            // If not a leaf, skip it.
            const bool isLeaf = (in_degree(v0, g) == 0) or (out_degree(v1, g) == 0);
            if(not isLeaf) {
                if(debug) {
                    cout << g[e].id << " skipped because it is not a leaf." << endl;
                }
                continue;
            }

            // Find the length of the shortest branch, in markers.
            uint64_t minLength = std::numeric_limits<uint64_t>::max();
            for(const E::Branch& branch: edge.branches) {
                minLength = min(minLength, uint64_t (branch.path.size()));
            }
            if(debug) {
                cout << g[e].id << " minLength " << minLength << endl;
            }

            if(minLength < pruneLength) {
                edgesToBeRemoved.push_back(e);
                if(debug) {
                    cout << g[e].id << " flagged to be removed." << endl;
                }
            }
        }


        // Remove the edges we found.
        if(edgesToBeRemoved.empty()) {
            break;
        }
        for(const edge_descriptor e: edgesToBeRemoved) {
            boost::remove_edge(e, g);
        }

        cout << "Prune iteration removed " << edgesToBeRemoved.size() << " edges." << endl;
    }

    performanceLog << timestamp << "AssemblyGraph2::prune ends" << endl;
}
#else



// Prune leaves in which all branches are shorter than the specified length.
// More efficient version that uses a stack.
void AssemblyGraph2::prune(uint64_t pruneLength)
{
    G& g = *this;

    performanceLog << timestamp << "AssemblyGraph2::prune begins" << endl;

    // A stack to contain all of the short leaves currently in the graph.
    std::stack<edge_descriptor> s;

    // Add to the stack all the short leaves initially present.
    BGL_FORALL_EDGES(e, g, G) {
        const vertex_descriptor v0 = source(e, g);
        const vertex_descriptor v1 = target(e, g);
        const E& edge = g[e];

        // If not a leaf, skip it.
        const bool isLeaf = (in_degree(v0, g) == 0) or (out_degree(v1, g) == 0);
        if(not isLeaf) {
            continue;
        }

        // If it is short, add it to the stack.
        if(edge.minimumPathLength() < pruneLength) {
            s.push(e);
        }
    }



    // Remove the short leaves.
    // Every time a short leaf is removed, check if any nearby short leaves
    // appear and add them to the stack if they do.
    uint64_t pruneCount = 0;
    while(not s.empty()) {

        // Get the next short leaf edge from the stack.
        // This will be removed.
        const edge_descriptor eA = s.top();
        s.pop();

        // Get the two vertices.
        const vertex_descriptor vA0 = source(eA, g);
        const vertex_descriptor vA1 = target(eA, g);

        // Compute the degrees.
        const uint64_t inDegree0 = in_degree(vA0, g);
        const uint64_t outDegree1 = out_degree(vA1, g);

        // Sanity checks that this is indeed a short leaf.
        SHASTA_ASSERT((inDegree0==0) or (outDegree1==0));
        SHASTA_ASSERT(g[eA].minimumPathLength() < pruneLength);



        // Add to the stack any short leaves that will be created when we remove eA.
        if(inDegree0 ==0) {
            if(outDegree1 == 0) {

                // This edge is isolated, so removing it will not
                // create any new leafs.

            } else {

                // There are no parents.
                // If the inDegree of vA1 is 1, removing eA would turn
                // all of its children into leaves.
                if(in_degree(vA1, g) == 1) {

                    // Loop over the children edges to see if any should
                    // be added to the stack.
                    BGL_FORALL_OUTEDGES(vA1, eB, g, G) {
                        if(eB == eA) {
                            continue;
                        }
                        if(g[eB].minimumPathLength() >= pruneLength) {
                            continue;
                        }
                        const vertex_descriptor vB1 = target(eB, g);
                        if(out_degree(vB1, g) != 0) {
                            s.push(eB);
                        }
                    }
                }

            }
        } else {
            if(outDegree1 == 0) {

                // There are no children.
                // If the outDegree of vA0 is 1, removing eA would turn
                // all of its parents into leaves.
                if(out_degree(vA0, g) == 1) {

                    // Loop over the parent edges to see if any should
                    // be added to the stack.
                    BGL_FORALL_INEDGES(vA0, eB, g, G) {
                        if(eB == eA) {
                            continue;
                        }
                        if(g[eB].minimumPathLength() >= pruneLength) {
                            continue;
                        }
                        const vertex_descriptor vB1 = source(eB, g);
                        if(in_degree(vB1, g) != 0) {
                            s.push(eB);
                        }
                    }
                }

            } else {

                // This edge is not a leaf!
                SHASTA_ASSERT(0);
            }

        }



        // Now we can remove eA.
        boost::remove_edge(eA, g);
        ++pruneCount;
    }


    if(pruneCount > 0) {
        cout << "Pruned " << pruneCount << " edges." << endl;
    }
    performanceLog << timestamp << "AssemblyGraph2::prune ends" << endl;
}
#endif



// Get the vertex descriptor for the vertex corresponding to
// a given MarkerGraph::VertexId, creating the vertex if necessary.
AssemblyGraph2::vertex_descriptor AssemblyGraph2::getVertex(MarkerGraph::VertexId vertexId)
{
    const auto it = vertexMap.find(vertexId);
    if(it == vertexMap.end()) {
        const vertex_descriptor v = add_vertex(V(vertexId), *this);
        vertexMap.insert(make_pair(vertexId, v));
        return v;
    } else {
        return it->second;
    }
}



// Create a new edges corresponding to the given path.
// Also create the vertices if necessary.
AssemblyGraph2::edge_descriptor AssemblyGraph2::addEdge(
    const MarkerGraphPath& path,
    bool containsSecondaryEdges,
    bool storeReadInformation,  // If true, store read information for newly created edges.
    bool assemble               // If true, assemble sequence for newly created edges
    )
{
    // Get the first and last edge of the path.
    const MarkerGraph::EdgeId edgeId0 = path.front();
    const MarkerGraph::EdgeId edgeId1 = path.back();
    const MarkerGraph::Edge edge0 = markerGraph.edges[edgeId0];
    const MarkerGraph::Edge edge1 = markerGraph.edges[edgeId1];

    // Get the first and last vertex of the path.
    // This creates the vertices if necessary.
    const MarkerGraph::VertexId vertexId0 = edge0.source;
    const MarkerGraph::VertexId vertexId1 = edge1.target;
    const vertex_descriptor v0 = getVertex(vertexId0);
    const vertex_descriptor v1 = getVertex(vertexId1);

    // Create the edge.
    edge_descriptor e;
    bool edgeWasAdded = false;
    tie(e, edgeWasAdded) = add_edge(v0, v1,
        E(nextId++, path, containsSecondaryEdges), *this);
    SHASTA_ASSERT(edgeWasAdded);

    if(storeReadInformation) {
        (*this)[e].storeReadInformation(markerGraph);
    }
    if(assemble) {
        AssemblyGraph2::assemble(e);
    }

    return e;
}



void AssemblyGraph2::writeCsv(const string& baseName) const
{
    writeVerticesCsv(baseName + "-Vertices.csv");
    writeEdgesCsv(baseName + "-Edges.csv");
    writeEdgeDetailsCsv(baseName + "-EdgeDetails.csv");
}



void AssemblyGraph2::writeVerticesCsv(const string& fileName) const
{
    const AssemblyGraph2& graph = *this;

    ofstream csv(fileName);
    csv << "VertexId0\n";
    BGL_FORALL_VERTICES(v, graph, AssemblyGraph2) {
        csv << graph[v].markerGraphVertexId << "\n";
    }
}



void AssemblyGraph2::writeEdgesCsv(const string& fileName) const
{
    const AssemblyGraph2& graph = *this;

    ofstream csv(fileName);
    csv << "VertexId0,VertexId1\n";
    BGL_FORALL_EDGES(e, graph, AssemblyGraph2) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        csv << graph[v0].markerGraphVertexId << ",";
        csv << graph[v1].markerGraphVertexId << "\n";
    }

}



void AssemblyGraph2::writeEdgeDetailsCsv(const string& fileName) const
{
    const G& g = *this;

    ofstream csv(fileName);
    csv << "FirstVertexId,LastVertexId,Branch,Position,EdgeId,VertexId0,VertexId1\n";

    // Loop over edges.
    BGL_FORALL_EDGES(e, g, G) {
        const E& edge = g[e];
        const vertex_descriptor v0 = source(e, g);
        const vertex_descriptor v1 = target(e, g);
        const MarkerGraph::VertexId vertexId0 = g[v0].markerGraphVertexId;
        const MarkerGraph::VertexId vertexId1 = g[v1].markerGraphVertexId;

        // Loop over branches of this edge.
        for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
            const E::Branch& branch = edge.branches[branchId];
            const MarkerGraphPath& path = branch.path;
            for(uint64_t j=0; j<path.size(); j++) {
                const MarkerGraph::EdgeId edgeId = path[j];
                const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
                csv << vertexId0 << ",";
                csv << vertexId1 << ",";
                csv << branchId << ",";
                csv << j << ",";
                csv << edgeId << ",";
                csv << edge.source << ",";
                csv << edge.target << "\n";
            }
        }
    }

}



// Assemble sequence for every marker graph path of every edge.
void AssemblyGraph2::assemble()
{
    performanceLog << timestamp << "AssemblyGraph2::assemble begins." << endl;
    G& g = *this;

    cout << timestamp << "assemble begins." << endl;

    // Use assembled sequence from the marker graph to obtain
    // assembled sequence for all edges.
    BGL_FORALL_EDGES(e, g, G) {
        assemble(e);
    }

    performanceLog << timestamp << "AssemblyGraph2::assemble ends." << endl;
}



// Assemble sequence for every marker graph path of every edge. Multithreaded version.
void AssemblyGraph2::assembleParallel(uint64_t threadCount)
{
    performanceLog << timestamp << "AssemblyGraph2::assembleParallel begins." << endl;
    G& g = *this;

    // Store a vector of edge descriptors for all edges, to be processed in parallel.
    assembleParallelData.allEdges.clear();
    BGL_FORALL_EDGES(e, g, G) {
        assembleParallelData.allEdges.push_back(e);
    }

    // Process all edges in parallel.
    const uint64_t batchSize = 100;
    setupLoadBalancing(assembleParallelData.allEdges.size(), batchSize);
    runThreads(&AssemblyGraph2::assembleThreadFunction, threadCount);
    assembleParallelData.allEdges.clear();

    performanceLog << timestamp << "AssemblyGraph2::assembleParallel ends." << endl;
}



void AssemblyGraph2::assembleThreadFunction(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all edges in this batch.
        for(uint64_t i=begin; i!=end; i++) {
            const edge_descriptor e = assembleParallelData.allEdges[i];
            assemble(e);
        }
    }
}



// Assemble sequence for every marker graph path of a given edge.
void AssemblyGraph2::assemble(edge_descriptor e)
{
    G& g = *this;


    E& edge = g[e];
    for(E::Branch& branch: edge.branches) {
        const MarkerGraphPath& path = branch.path;

        AssembledSegment assembledSegment;
        MarkerGraph::EdgeId const * const begin = &path[0];
        MarkerGraph::EdgeId const * const end = begin + path.size();
        const span<const MarkerGraph::EdgeId> pathSpan(begin, end);
        assembleMarkerGraphPath(readRepresentation, k,
            markers, markerGraph, pathSpan, false, assembledSegment);



        // Store the sequence, excluding the first and last k/2 RLE bases.

        // Compute the number of raw bases to skip at the beginning.
        const uint64_t beginSkip = std::accumulate(
            assembledSegment.repeatCounts.begin(),
            assembledSegment.repeatCounts.begin() + k/2, 0ULL);

        // Compute the number of raw bases to skip at the end.
        const uint64_t endSkip = std::accumulate(
            assembledSegment.repeatCounts.end() - k/2,
            assembledSegment.repeatCounts.end(), 0ULL);

        // Copy the raw bases, excluding those corresponding to the first and last
        // k/2 RLE bases.
        branch.rawSequence.resize(assembledSegment.rawSequence.size() - beginSkip - endSkip);
        copy(
            assembledSegment.rawSequence.begin() + beginSkip,
            assembledSegment.rawSequence.end() - endSkip,
            branch.rawSequence.begin());
    }
}



// Finds edges that form bubbles, then combine
// each of them into a single edge with multiple paths.
void AssemblyGraph2::gatherBubbles()
{
    G& g = *this;

    performanceLog << timestamp << "AssemblyGraph2::gatherBubbles begins." << endl;


    // Look for sets of parallel edges v0->v1.
    BGL_FORALL_VERTICES(v0, g, G) {

        // Map keyed by child vertex v1, and containing for each v1
        // a vector of all the edges v0->v1.
        std::map< vertex_descriptor, vector<edge_descriptor> > successorMap;
        BGL_FORALL_OUTEDGES(v0, e01, g, G) {
            const vertex_descriptor v1 = target(e01, g);
            successorMap[v1].push_back(e01);
        }



        // Process each set with the same v1.
        for(const auto& p: successorMap) {

            // Get the edges v0->v1.
            const vector<edge_descriptor>& edges01= p.second;
            const uint64_t ploidy = edges01.size();

            // If less than two edges, it is not a bubble, so there is
            // nothing to do.
            if(ploidy < 2) {
                continue;
            }
            const vertex_descriptor v1 = p.first;

            // Create the bubble remove these edges.
            createBubble(v0, v1, edges01);
        }
    }


    performanceLog << timestamp << "AssemblyGraph2::gatherBubbles ends." << endl;


}



void AssemblyGraph2::writePloidyHistogram(ostream& s) const
{
    const G& g = *this;

    vector<uint64_t> ploidyHistogram;
    BGL_FORALL_EDGES(e, g, G) {
        const uint64_t ploidy = g[e].ploidy();
        if(ploidy >= ploidyHistogram.size()) {
            ploidyHistogram.resize(ploidy + 1);
        }
        ++ploidyHistogram[ploidy];
    }
    for(uint64_t ploidy=1; ploidy<ploidyHistogram.size(); ploidy++) {
        s << "Ploidy " << ploidy << ": " << ploidyHistogram[ploidy] << " edges." << endl;
    }

}



AssemblyGraph2::edge_descriptor AssemblyGraph2::createBubble(
    vertex_descriptor v0,
    vertex_descriptor v1,
    const vector<edge_descriptor>& edges01)
{
    G& g = *this;

    // Sanity check.
    for(const edge_descriptor e01: edges01) {
        SHASTA_ASSERT(source(e01, g) == v0);
        SHASTA_ASSERT(target(e01, g) == v1);
    }

    // Create the new edge to replace the old ones.
    edge_descriptor eNew;
    bool edgeWasAdded = false;
    tie(eNew, edgeWasAdded) = add_edge(v0, v1, E(nextId++), g);
    SHASTA_ASSERT(edgeWasAdded);
    E& edgeNew = g[eNew];

    // Copy the branches to the new edge and remove the old edges.
    for(const edge_descriptor e01: edges01) {
        const E& edge01 = g[e01];
        copy(edge01.branches.begin(), edge01.branches.end(),
           back_inserter(edgeNew.branches));
        boost::remove_edge(e01, g);
    }

    return eNew;
}



// The first two bool flags work as follows:

// - writeSequence=false, writeSequenceLengthInMarkers=false:
//   No sequence output (sequence is written as * instead),
//   and sequence length is expressed in bases.
//   This can only be called after storeGfaSequence has been called.

// - writeSequence=false, writeSequenceLengthInMarkers=true:
//   No sequence output (sequence is written as * instead),
//   and sequence length is expressed in markers.
//   This does not use the stored gfa sequences and
//   so can be called early, before storeGfaSequence has been called.

// - writeSequence=true, writeSequenceLengthInMarkers=false:
//   Complete output of sequence, with length in bases.
//   This can only be called after storeGfaSequence has been called.

// - writeSequence=true, writeSequenceLengthInMarkers=true:
//   This comnbination is illegal and causes an assertion.

void AssemblyGraph2::writeDetailed(
    const string& baseName,
    bool writeSequence,
    bool writeSequenceLengthInMarkers,
    bool writeCsv,
    bool writeGfa,
    bool writeFasta) const
{
    performanceLog << timestamp << "AssemblyGraph2::writeDetailed begins." << endl;

    // Check that we are not called with the forbidden combination
    // (see above comments).
    SHASTA_ASSERT(not(writeSequence and writeSequenceLengthInMarkers));

    const G& g = *this;


    // Open the accompanying csv file and write the header.
    ofstream csv;
    if(writeCsv) {
        csv.open(baseName + ".csv");
        csv << "Name,Component,Phase,Unphased strength,Color,"
            "First marker graph vertex,Last marker graph vertex,"
            "First marker graph edge,Last marker graph edge,"
            "Length in markers,"
            "Length in bases,"
            "Secondary,Period,"
            "Minimum marker graph edge coverage,Average marker graph edge coverage,Number of distinct oriented reads,";
        if(writeSequence) {
            csv << "Sequence,";
        }
        csv << "\n";
    }


    // Open the fasta file.
    ofstream fasta;
    if(writeFasta) {
        fasta.open(baseName + ".fasta");
    }

    // Create a GFA with a segment for each branch, then write it out.
    GfaAssemblyGraph<vertex_descriptor> gfa;
    BGL_FORALL_EDGES(e, g, G) {
        const E& edge = g[e];
        const vertex_descriptor v0 = source(e, g);
        const vertex_descriptor v1 = target(e, g);

        for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
            const E::Branch& branch = edge.branches[branchId];

            if(writeGfa) {
                if(writeSequence) {
                    gfa.addSegment(edge.pathId(branchId), v0, v1, branch.gfaSequence);
                } else {
                    if(writeSequenceLengthInMarkers) {
                        gfa.addSegment(edge.pathId(branchId), v0, v1, branch.path.size());
                    } else {
                        gfa.addSegment(edge.pathId(branchId), v0, v1, branch.gfaSequence.size());
                    }
                }
            }

            if(writeFasta) {
                fasta << ">" << edge.pathId(branchId) << " " << branch.gfaSequence.size() << "\n";
                copy(branch.gfaSequence.begin(), branch.gfaSequence.end(), ostream_iterator<Base>(fasta));
                fasta << "\n";
            }



            // Write a line for this segment to the csv file.
            if(writeCsv) {

                // Get some information we need below.
                const uint64_t lengthInMarkers = branch.path.size();
                SHASTA_ASSERT(lengthInMarkers > 0);
                const MarkerGraphEdgeId firstMarkerGraphEdgeId = branch.path.front();
                const MarkerGraphEdgeId lastMarkerGraphEdgeId = branch.path.back();
                const MarkerGraphVertexId firstMarkerGraphVertexId = markerGraph.edges[firstMarkerGraphEdgeId].source;
                const MarkerGraphVertexId lastMarkerGraphVertexId = markerGraph.edges[lastMarkerGraphEdgeId].target;
                const string color = edge.color(branchId);

                csv <<
                    edge.pathId(branchId) << ",";
                if(edge.componentId != std::numeric_limits<uint64_t>::max()) {
                    csv << edge.componentId;
                }
                csv << ",";

                if(edge.phase != std::numeric_limits<uint64_t>::max()) {
                    csv << (branchId == edge.phase ? 0 : 1);
                }
                csv << ",";

                if(edge.isBubble() and (edge.isBad or edge.phase == std::numeric_limits<uint64_t>::max())) {
                    if(branchId == edge.getStrongestBranchId()) {
                        csv << "Strong";
                    } else {
                        csv << "Weak";
                    }
                }
                csv << ",";

                csv <<
                    color << "," <<
                    firstMarkerGraphVertexId << "," <<
                    lastMarkerGraphVertexId << "," <<
                    firstMarkerGraphEdgeId << "," <<
                    lastMarkerGraphEdgeId << "," <<
                    lengthInMarkers << ",";

                if(writeSequence or (not writeSequenceLengthInMarkers)) {
                    csv << branch.gfaSequence.size();
                }
                csv << ",";

                csv <<
                    (branch.containsSecondaryEdges ? "S" : "") << "," <<
                    (edge.period ? to_string(edge.period) : string()) << "," <<
                    branch.minimumCoverage << "," <<
                    branch.averageCoverage() << "," <<
                    branch.orientedReadIds.size() << ",";
                if(writeSequence) {
                    if(branch.gfaSequence.size() == 0) {
                        csv << "-";
                    } else if(branch.gfaSequence.size() <= 6) {
                        copy(branch.gfaSequence.begin(), branch.gfaSequence.end(),
                            ostream_iterator<Base>(csv));
                    } else {
                        csv << "...";
                    }
                    csv << ",";
                }
                csv << "\n";
            }
        }
    }



    // Add paths.
    if(writeGfa) {
        for(uint64_t bubbleChainId=0; bubbleChainId<uint64_t(bubbleChains.size()); bubbleChainId++) {
            const BubbleChain& bubbleChain = bubbleChains[bubbleChainId];
            for(uint64_t phasingRegionId=0;
                phasingRegionId<uint64_t(bubbleChain.phasingRegions.size()); phasingRegionId++) {
                const auto& phasingRegion = bubbleChain.phasingRegions[phasingRegionId];

                vector<string> path0;
                vector<string> path1;

                for(uint64_t position=phasingRegion.firstPosition;
                    position<=phasingRegion.lastPosition; position++) {
                    const edge_descriptor e = bubbleChain.edges[position];
                    const E& edge = g[e];

                    if(edge.componentId == std::numeric_limits<uint64_t>::max()) {

                        // This edge is homozygous or unphased.
                        const string segmentName = edge.pathId(edge.getStrongestBranchId());
                        path0.push_back(segmentName);
                        path1.push_back(segmentName);

                    } else {

                        // This edge is diploid and phased.
                        SHASTA_ASSERT(edge.ploidy() == 2);
                        SHASTA_ASSERT(edge.componentId == phasingRegion.componentId);

                        string segmentName0 = edge.pathId(0);
                        string segmentName1 = edge.pathId(1);

                        if(edge.phase == 0) {
                            path0.push_back(segmentName0);
                            path1.push_back(segmentName1);
                        } else {
                            path0.push_back(segmentName1);
                            path1.push_back(segmentName0);
                        }

                    }
                }

                // Each phased (diploid) region generates two paths.
                // Each unphased (haploid) region generates one path.
                if(phasingRegion.isPhased) {
                    const string idPrefix =
                        "PR." +
                        to_string(bubbleChainId) + "." +
                        to_string(phasingRegionId) + "." +
                        to_string(phasingRegion.componentId) + ".";
                    gfa.addPath(idPrefix + "0", path0);
                    gfa.addPath(idPrefix + "1", path1);
                } else {
                    SHASTA_ASSERT(path0 == path1);
                    const string idString =
                        "UR." +
                        to_string(bubbleChainId) + "." +
                        to_string(phasingRegionId);
                    gfa.addPath(idString, path0);
                }
            }
        }
    }



    // Write out the GFA.
    if(writeGfa) {
        gfa.write(baseName + ".gfa");
    }

    performanceLog << timestamp << "AssemblyGraph2::writeDetailed ends." << endl;
}



void AssemblyGraph2::writeDetailedEarly(const string& baseName)
{
    const bool writeSequence = false;
    const bool writeSequenceLengthInMarkers = true;
    const bool writeCsv = true;
    const bool writeGfa = true;
    const bool writeFasta = false;

    writeDetailed(baseName,
        writeSequence, writeSequenceLengthInMarkers, writeCsv, writeGfa, writeFasta);
}



void AssemblyGraph2::writeHaploid(
    const string& baseName,
    bool writeSequence,
    bool writeCsv,
    bool writeGfa,
    bool writeFasta,
    AssemblyGraph2Statistics* statistics) const
{
    performanceLog << timestamp << "AssemblyGraph2::writeHaploid begins." << endl;
    const G& g = *this;

    vector<uint64_t> bubbleChainLengths;
    uint64_t totalNonBubbleChainLength = 0;

    // Open the fasta file.
    ofstream fasta;
    if(writeFasta) {
        fasta.open(baseName + ".fasta");
    }

    // Create a GFA and add a segment for each edge that is not part
    // of a bubble chain.
    GfaAssemblyGraph<vertex_descriptor> gfa;
    BGL_FORALL_EDGES(e, g, G) {
        const E& edge = g[e];
        if(edge.bubbleChain.first) {
            continue;
        }

        const vertex_descriptor v0 = source(e, g);
        const vertex_descriptor v1 = target(e, g);

        for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
            const E::Branch& branch = edge.branches[branchId];
            totalNonBubbleChainLength += branch.gfaSequence.size();

            if(writeGfa) {
                if(writeSequence) {
                    gfa.addSegment(edge.pathId(branchId), v0, v1, branch.gfaSequence);
                } else {
                    gfa.addSegment(edge.pathId(branchId), v0, v1, branch.gfaSequence.size());
                }
            }

            if(writeFasta) {
                fasta << ">" << edge.pathId(branchId) << " " << branch.gfaSequence.size() << "\n";
                copy(branch.gfaSequence.begin(), branch.gfaSequence.end(), ostream_iterator<Base>(fasta));
                fasta << "\n";
            }
        }
    }



    // Add a segment for each bubble chain.
    for(uint64_t bubbleChainId=0; bubbleChainId<uint64_t(bubbleChains.size()); bubbleChainId++) {
        const BubbleChain& bubbleChain = bubbleChains[bubbleChainId];
        const vertex_descriptor v0 = source(bubbleChain.edges.front(), g);
        const vertex_descriptor v1 = target(bubbleChain.edges.back(), g);

        vector<Base> sequence;
        computeBubbleChainGfaSequence(bubbleChain, sequence);
        bubbleChainLengths.push_back(uint64_t(sequence.size()));

        const string idString = "BC." + to_string(bubbleChainId);

        if(writeGfa) {
            if(writeSequence) {
                gfa.addSegment(idString, v0, v1, sequence);
            } else {
                gfa.addSegment(idString, v0, v1, sequence.size());
            }
        }

        if(writeFasta) {
            fasta << ">" << idString << " " << sequence.size() << "\n";
            copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
            fasta << "\n";
        }
    }



    // Write the GFA.
    if(writeGfa) {
        gfa.write(baseName + ".gfa");
    }



    // Also write a csv file that can be used in Bandage.
    if(writeCsv) {
        ofstream csv(baseName + ".csv");
        csv << "Name,ComponentId,Phase,Color,First marker graph edge,Last marker graph edge,"
            "Secondary,Period,"
            "Minimum edge coverage,Average edge coverage,Number of distinct oriented reads,\n";



        // Write a line to csv for each edge that is not part of a bubble chain.
        BGL_FORALL_EDGES(e, g, G) {
            const E& edge = g[e];
            if(edge.bubbleChain.first) {
                continue;
            }

            for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
                const E::Branch& branch = edge.branches[branchId];

                const string color = edge.color(branchId);
                csv <<
                    edge.pathId(branchId) << ",";
                if(edge.componentId != std::numeric_limits<uint64_t>::max()) {
                    csv << edge.componentId;
                }
                csv << ",";
                if(edge.phase != std::numeric_limits<uint64_t>::max()) {
                    csv << (branchId == edge.phase ? 0 : 1);
                }
                csv <<
                    "," <<
                    color << "," <<
                    branch.path.front() << "," << branch.path.back() << "," <<
                    (branch.containsSecondaryEdges ? "S" : "") << "," <<
                    (edge.period ? to_string(edge.period) : string()) << "," <<
                    branch.minimumCoverage << "," <<
                    branch.averageCoverage() << "," <<
                    branch.orientedReadIds.size() << "\n";
            }
        }



        // Write a line to csv for each bubble chain.
        for(uint64_t bubbleChainId=0; bubbleChainId<uint64_t(bubbleChains.size()); bubbleChainId++) {
            const string idString = "BC." + to_string(bubbleChainId);
            csv << idString << ",,,Cyan\n";
        }



        // Statistics.
        const uint64_t totalLength =
            accumulate(bubbleChainLengths.begin(), bubbleChainLengths.end(), 0ULL);
        sort(bubbleChainLengths.begin(), bubbleChainLengths.end(), std::greater<uint64_t>());
        uint64_t n50 = 0;
        uint64_t cumulativeLength = 0;
        for(const uint64_t length: bubbleChainLengths) {
            cumulativeLength += length;
            if(cumulativeLength >= totalLength/2) {
                n50 = length;
                break;
            }
        }
        cout << "Total length of bubble chains " << totalLength <<
            ", N50 " << n50 << endl;
        // cout << "Total length assembled outside of bubble chains " << totalNonBubbleChainLength << endl;
        if(statistics) {
            statistics->totalBubbleChainLength = totalLength;
            statistics->bubbleChainN50 = n50;
        }
    }

    performanceLog << timestamp << "AssemblyGraph2::writeHaploid ends." << endl;
}



void AssemblyGraph2::writePhased(
    const string& baseName,
    bool writeSequence,
    bool writeCsv,
    bool writeGfa,
    bool writeFasta,
    AssemblyGraph2Statistics* statistics) const
{
    performanceLog << timestamp << "AssemblyGraph2::writePhased begins." << endl;
    const G& g = *this;

    // Length statistics.
    uint64_t totalHaploidBases = 0;
    uint64_t totalDiploidBases = 0;
    uint64_t totalNonBubbleChainBases = 0;

    // Tables used to compute N50 for diploid and haploid segments
    // that are part of bubble chains.
    vector<uint64_t> haploidLengths;
    vector<uint64_t> diploidLengths;

    // Also write a csv file that can be used in Bandage.
    ofstream csv;
    if(writeCsv) {
        csv.open(baseName + ".csv");
        csv << "Name,Position in bubble chain,Ploidy,Bubble chain,Component,Haplotype,Length,Color\n";
    }

    // Open the fasta file.
    ofstream fasta;
    if(writeFasta) {
        fasta.open(baseName + ".fasta");
    }

    // Create a GFA and add a segment for each edge that is not part
    // of a bubble chain.
    GfaAssemblyGraph<vertex_descriptor> gfa;
    BGL_FORALL_EDGES(e, g, G) {
        const E& edge = g[e];
        if(edge.bubbleChain.first) {
            continue;
        }

        const vertex_descriptor v0 = source(e, g);
        const vertex_descriptor v1 = target(e, g);

        for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
            const E::Branch& branch = edge.branches[branchId];
            const string segmentId = edge.pathId(branchId);

            if(writeGfa) {
                if(writeSequence) {
                    gfa.addSegment(segmentId, v0, v1, branch.gfaSequence);
                } else {
                    gfa.addSegment(segmentId, v0, v1, branch.gfaSequence.size());
                }
            }

            if(writeFasta) {
                fasta << ">" << segmentId << " " << branch.gfaSequence.size() << "\n";
                copy(branch.gfaSequence.begin(), branch.gfaSequence.end(), ostream_iterator<Base>(fasta));
                fasta << "\n";
            }

            if(writeCsv) {
                csv << segmentId << ",,,,,,,#808080\n";
            }
            totalNonBubbleChainBases += uint64_t(branch.gfaSequence.size());
        }
    }



    // Add one or two segments, depending on ploidy, for each phasing region
    // of each bubble chain.
    vector<Base> sequence;
    for(uint64_t bubbleChainId=0; bubbleChainId<uint64_t(bubbleChains.size()); bubbleChainId++) {
        const BubbleChain& bubbleChain = bubbleChains[bubbleChainId];
        for(uint64_t phasingRegionId=0;
            phasingRegionId<uint64_t(bubbleChain.phasingRegions.size()); phasingRegionId++) {
            const auto& phasingRegion = bubbleChain.phasingRegions[phasingRegionId];

            const vertex_descriptor v0 = source(bubbleChain.edges[phasingRegion.firstPosition], g);
            const vertex_descriptor v1 = target(bubbleChain.edges[phasingRegion.lastPosition], g);

            if(phasingRegion.isPhased) {

                const string namePrefix =
                    "PR." +
                    to_string(bubbleChainId) + "." +
                    to_string(phasingRegionId) + "." +
                    to_string(phasingRegion.componentId) + ".";

                const string name0 = namePrefix + "0";
                computePhasedRegionGfaSequence(bubbleChain, phasingRegion, 0, sequence);

                if(writeGfa) {
                    if(writeSequence) {
                        gfa.addSegment(name0, v0, v1, sequence);
                    } else {
                        gfa.addSegment(name0, v0, v1, sequence.size());
                    }
                }

                if(writeFasta) {
                    fasta << ">" << name0 << " " << sequence.size() << "\n";
                    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
                    fasta << "\n";
                }

                totalDiploidBases += uint64_t(sequence.size());
                diploidLengths.push_back(uint64_t(sequence.size()));

                if(writeCsv) {
                    csv <<
                        name0 << "," <<
                        phasingRegionId << "," <<
                        "2," <<
                        bubbleChainId << "," <<
                        phasingRegion.componentId << "," <<
                        "0," <<
                        sequence.size() << ","
                        "Green\n";
                }

                const string name1 = namePrefix + "1";
                computePhasedRegionGfaSequence(bubbleChain, phasingRegion, 1, sequence);

                if(writeGfa) {
                    if(writeSequence) {
                        gfa.addSegment(name1, v0, v1, sequence);
                    } else {
                        gfa.addSegment(name1, v0, v1, sequence.size());
                    }
                }

                if(writeFasta) {
                    fasta << ">" << name1 << " " << sequence.size() << "\n";
                    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
                    fasta << "\n";
                }

                totalDiploidBases += uint64_t(sequence.size());
                diploidLengths.push_back(uint64_t(sequence.size()));

                if(writeCsv) {
                    csv <<
                        name1 << "," <<
                        phasingRegionId << "," <<
                        "2," <<
                        bubbleChainId << "," <<
                        phasingRegion.componentId << "," <<
                        "1," <<
                        sequence.size() << ","
                        "Green\n";
                }

            } else {

                computeUnphasedRegionGfaSequence(bubbleChain, phasingRegion, sequence);
                const string name = "UR." + to_string(bubbleChainId) + "." + to_string(phasingRegionId);

                if(writeGfa) {
                    if(writeSequence) {
                        gfa.addSegment(name, v0, v1, sequence);
                    } else {
                        gfa.addSegment(name, v0, v1, sequence.size());
                    }
                }

                totalHaploidBases += uint64_t(sequence.size());
                haploidLengths.push_back(uint64_t(sequence.size()));

                if(writeFasta) {
                    fasta << ">" << name << " " << sequence.size() << "\n";
                    copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fasta));
                    fasta << "\n";
                }

                if(writeCsv) {
                    csv <<
                        name << "," <<
                        phasingRegionId << "," <<
                        "1," <<
                        bubbleChainId << "," <<
                        "," <<
                        "," <<
                        sequence.size() << ","
                        "#eb4034\n";   // Near red.
                }

            }

        }
    }



    // Write the GFA.
    if(writeGfa) {
        gfa.write(baseName + ".gfa");
    }



    if(writeCsv) {
        // Compute N50 for regions assembled diploid and phased.
        sort(diploidLengths.begin(), diploidLengths.end(), std::greater<uint64_t>());
        /*
        cout << "Diploid lengths: ";
        copy(diploidLengths.begin(), diploidLengths.end(), ostream_iterator<uint64_t>(cout, " "));
        cout << endl;
        */
        uint64_t diploidN50 = 0;
        uint64_t cumulativeDiploidLength = 0;
        for(const uint64_t length: diploidLengths) {
            cumulativeDiploidLength += length;
            if(cumulativeDiploidLength >= totalDiploidBases/2) {
                diploidN50 = length;
                break;
            }
        }

        // Compute N50 for regions assembled haploid in bubble chains
        sort(haploidLengths.begin(), haploidLengths.end(), std::greater<uint64_t>());
        /*
        cout << "Haploid lengths: ";
        copy(haploidLengths.begin(), haploidLengths.end(), ostream_iterator<uint64_t>(cout, " "));
        cout << endl;
        */
        uint64_t haploidN50 = 0;
        uint64_t cumulativeHaploidLength = 0;
        for(const uint64_t length: haploidLengths) {
            cumulativeHaploidLength += length;
            if(cumulativeHaploidLength >= totalHaploidBases/2) {
                haploidN50 = length;
                break;
            }
        }

        cout << "Assembled diploid in bubble chains and phased: total " << totalDiploidBases <<
            " (" << totalDiploidBases/2 << " per haplotype), N50 " << diploidN50 << "."  << endl;
        cout << "Total length assembled haploid in bubble chains: " << totalHaploidBases <<
            ", N50 " << haploidN50 << "." << endl;
        cout << "Total genome length assembled in bubble chains, averaged over haplotypes: " <<
            totalDiploidBases/2 + totalHaploidBases << endl;
        cout << "Total length assembled outside bubble chains: " <<
            totalNonBubbleChainBases << endl;

        if(statistics) {
            statistics->totalDiploidLengthBothHaplotypes = totalDiploidBases;
            statistics->diploidN50 = diploidN50;
            statistics->totalHaploidLength = totalHaploidBases;
            statistics->haploidN50 = haploidN50;
            statistics->outsideBubbleChainsLength = totalNonBubbleChainBases;
        }
    }

    performanceLog << timestamp << "AssemblyGraph2::writePhased ends." << endl;
}



// This writes a csv file that relates coordinates in the phased
// assembly to coordinates in the detailed assembly.
void AssemblyGraph2::writePhasedDetails() const
{
    const AssemblyGraph2& g = *this;

    ofstream csv("Assembly-Phased-Details.csv");
    csv << "Segment,Detailed segment,Length,Begin,End\n";


    // Loop over all segments in bubble chains (segment names beginning
    // with "PR." or "UR.").
    for(uint64_t bubbleChainId=0; bubbleChainId<uint64_t(bubbleChains.size()); bubbleChainId++) {
        const BubbleChain& bubbleChain = bubbleChains[bubbleChainId];
        for(uint64_t phasingRegionId=0;
            phasingRegionId<uint64_t(bubbleChain.phasingRegions.size()); phasingRegionId++) {
            const auto& phasingRegion = bubbleChain.phasingRegions[phasingRegionId];

            if(phasingRegion.isPhased) {

                // Phased region. The code is similar to computePhasedRegionGfaSequence.
                const string namePrefix =
                    "PR." +
                    to_string(bubbleChainId) + "." +
                    to_string(phasingRegionId) + "." +
                    to_string(phasingRegion.componentId) + ".";

                for(uint64_t haplotype=0; haplotype<2; haplotype++) {
                    const string name = namePrefix + to_string(haplotype);

                    uint64_t n = 0; // The sequence accumulated so far.
                    for(uint64_t position=phasingRegion.firstPosition;
                        position<=phasingRegion.lastPosition; position++) {
                        const edge_descriptor e = bubbleChain.edges[position];
                        const E& edge = g[e];

                        if(edge.componentId == std::numeric_limits<uint64_t>::max()) {

                            // This edge is homozygous or unphased.
                            const uint64_t branchId = edge.getStrongestBranchId();
                            const E::Branch& branch = edge.branches[branchId];
                            const uint64_t length = branch.gfaSequence.size();
                            const uint64_t begin = n;
                            const uint64_t end = begin + length;

                            csv << name << ",";
                            csv << edge.pathId(branchId) << ",";
                            csv << length << ",";
                            csv << begin << ",";
                            csv << end << "\n";

                            n = end;

                        } else {

                            // This edge is diploid and phased.
                            SHASTA_ASSERT(edge.ploidy() == 2);
                            SHASTA_ASSERT(edge.componentId == phasingRegion.componentId);

                            // Figure out which branch we need.
                            uint64_t branchId = 0;
                            if(edge.phase != haplotype) {
                                branchId = 1;
                            }

                            const E::Branch& branch = edge.branches[branchId];
                            const uint64_t length = branch.gfaSequence.size();
                            const uint64_t begin = n;
                            const uint64_t end = begin + length;

                            csv << name << ",";
                            csv << edge.pathId(branchId) << ",";
                            csv << length << ",";
                            csv << begin << ",";
                            csv << end << "\n";

                            n = end;
                        }
                    }
                }


            } else {

                // Unphased region. The code is similar to computeUnphasedRegionGfaSequence.
                const string name = "UR." + to_string(bubbleChainId) + "." + to_string(phasingRegionId);

                uint64_t n = 0; // The sequence accumulated so far.
                for(uint64_t position=phasingRegion.firstPosition;
                    position<=phasingRegion.lastPosition; position++) {
                    const edge_descriptor e = bubbleChain.edges[position];
                    const E& edge = g[e];
                    const E::Branch& branch = edge.branches[edge.getStrongestBranchId()];
                    const uint64_t length = branch.gfaSequence.size();
                    const uint64_t begin = n;
                    const uint64_t end = begin + length;

                    csv << name << ",";
                    csv << edge.pathId(0) << ",";
                    csv << length << ",";
                    csv << begin << ",";
                    csv << end << "\n";

                    n = end;
                }


            }

        }
    }
}



// Compute the gfa sequence of a bubble chain
// by concatenating gfa sequence of the strongest branch of
// each of tis edges.
void AssemblyGraph2::computeBubbleChainGfaSequence(
    const BubbleChain& bubbleChain,
    vector<Base>& sequence
    ) const
{
    const G& g = *this;

    sequence.clear();
    for(const edge_descriptor e: bubbleChain.edges) {
        const E& edge = g[e];
        const E::Branch& branch = edge.branches[edge.getStrongestBranchId()];
        copy(branch.gfaSequence.begin(), branch.gfaSequence.end(),
            back_inserter(sequence));
    }
}



// Compute the gfa sequence of an unphased region
// by concatenating gfa sequence of the strongest branch of
// each of this edges.
void AssemblyGraph2::computeUnphasedRegionGfaSequence(
    const BubbleChain& bubbleChain,
    const BubbleChain::PhasingRegion& phasingRegion,
    vector<Base>& sequence
    ) const
{
    const G& g = *this;

    sequence.clear();
    for(uint64_t position=phasingRegion.firstPosition;
        position<=phasingRegion.lastPosition; position++) {
        const edge_descriptor e = bubbleChain.edges[position];
        const E& edge = g[e];
        const E::Branch& branch = edge.branches[edge.getStrongestBranchId()];
        copy(branch.gfaSequence.begin(), branch.gfaSequence.end(),
            back_inserter(sequence));
    }
}



// Compute the gfa sequence of an haplotype of a phased region.
void AssemblyGraph2::computePhasedRegionGfaSequence(
    const BubbleChain& bubbleChain,
    const BubbleChain::PhasingRegion& phasingRegion,
    uint64_t haplotype,
    vector<Base>& sequence
    ) const
{
    const G& g = *this;

    sequence.clear();
    for(uint64_t position=phasingRegion.firstPosition;
        position<=phasingRegion.lastPosition; position++) {
        const edge_descriptor e = bubbleChain.edges[position];
        const E& edge = g[e];

        if(edge.componentId == std::numeric_limits<uint64_t>::max()) {

            // This edge is homozygous or unphased.
            const E::Branch& branch = edge.branches[edge.getStrongestBranchId()];
            copy(branch.gfaSequence.begin(), branch.gfaSequence.end(),
                back_inserter(sequence));

        } else {

            // This edge is diploid and phased.
            SHASTA_ASSERT(edge.ploidy() == 2);
            SHASTA_ASSERT(edge.componentId == phasingRegion.componentId);

            // Figure out which branch we need.
            uint64_t branchId = 0;
            if(edge.phase != haplotype) {
                branchId = 1;
            }

            const E::Branch& branch = edge.branches[branchId];
            copy(branch.gfaSequence.begin(), branch.gfaSequence.end(),
                back_inserter(sequence));
        }
    }
}



string AssemblyGraph2Edge::color(uint64_t branchId) const
{
    if(isBubble()) {

        // A bad or unphased bubble is grey - darker on the strongest branch.
        if(isBad or phase == std::numeric_limits<uint64_t>::max()) {
            if(branchId == getStrongestBranchId()) {
                return "#808080";
            } else {
                return "#C0C0C0";
            }
        }

        if(branchId == phase) {
            return "#bf4040";   // HSL(0,50%,50%) Pink
        } else {
            return "#4040bf";   // HSL(240,50%,50%) Blue
        }

    } else {

        // It is not a bubble.
        return "#808080";
    }

}



// For each bubble edge, compute the number of raw sequence bases
// transferred in each direction for gfa output.
void AssemblyGraph2::countTransferredBases()
{
    G& g = *this;

    BGL_FORALL_EDGES(e, g, G) {
        E& edge = g[e];
        edge.backwardTransferCount = 0;
        edge.forwardTransferCount = 0;

        // To transfer any bases, the edge
        // must be a bubble preceded and followed by a single non-bubble.
        // If these conditions are not satisfied, leave the numbers
        // of transfered bases at 0.

        // The edge must be a bubble.
        if(not edge.isBubble()) {
            continue;
        }

        // v0 must have in-degree and out-degree 1.
        const vertex_descriptor v0 = source(e, g);
        if(in_degree(v0, g) != 1) {
            continue;
        }
        if(out_degree(v0, g) != 1) {
            continue;
        }

        // v1 must have in-degree and out-degree 1.
        const vertex_descriptor v1 = target(e, g);
        if(in_degree(v1, g) != 1) {
            continue;
        }
        if(out_degree(v1, g) != 1) {
            continue;
        }

        // The previous edge must not be a bubble.
        in_edge_iterator itPrevious;
        tie(itPrevious, ignore) = in_edges(v0, g);
        const E& previousEdge = g[*itPrevious];
        if(previousEdge.isBubble()) {
            continue;
        }

        // The next edge must not be a bubble.
        out_edge_iterator itNext;
        tie(itNext, ignore) = out_edges(v1, g);
        const E& nextEdge = g[*itNext];
        if(nextEdge.isBubble()) {
            continue;
        }

        // All conditions are satisfied.
        // Set the number of transfered bases equal to the number
        // of common identical prefix/suffix bases for all the
        // branches of this bubble.
        edge.backwardTransferCount = edge.countCommonPrefixBases();
        edge.forwardTransferCount = edge.countCommonSuffixBases();

        // Make sure we don't transfer more than the length of the
        // shortest branch of this edge.
        uint64_t shortestBranchLength = std::numeric_limits<uint64_t>::max();
        for(const E::Branch& branch:edge.branches) {
            shortestBranchLength = min(shortestBranchLength, uint64_t(branch.rawSequence.size()));
        }
        while(true) {
            if(edge.backwardTransferCount + edge.forwardTransferCount <= shortestBranchLength) {
                break;
            }
            --edge.backwardTransferCount;
            if(edge.backwardTransferCount + edge.forwardTransferCount <= shortestBranchLength) {
                break;
            }
            --edge.forwardTransferCount;
        }
    }
}



// Store GFA sequence in each edge.
// It is obtained from raw sequence by transferring identical bases
// of common sequence between all branches of a bubble to the
// preceding or following non-bubble edge.
void AssemblyGraph2::storeGfaSequence()
{
    performanceLog << timestamp << "storeGfaSequence begins." << endl;

    G& g = *this;

    // Count the number of sequence bases transferred forward/backward
    // from each bubble edge to adjacent non-bubble edges.
    countTransferredBases();



    BGL_FORALL_EDGES(e, g, G) {
        E& edge = g[e];
        const vertex_descriptor v0 = source(e, g);
        const vertex_descriptor v1 = target(e, g);

        for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
            E::Branch& branch = edge.branches[branchId];

            branch.gfaSequence.clear();

            // Add the sequence transferred forward by the preceding bubble, if appropriate.
            if(not edge.isBubble()) {
                if(in_degree(v0, g)==1 and out_degree(v0, g)==1) {
                    in_edge_iterator it;
                    tie(it, ignore) = in_edges(v0, g);
                    const E& previousEdge = g[*it];
                    if(previousEdge.isBubble()) {
                        const vector<Base>& s = previousEdge.branches.front().rawSequence;
                        copy(s.end() - previousEdge.forwardTransferCount, s.end(),
                            back_inserter(branch.gfaSequence));
                    }
                }
            }

            // Add the sequence of this branch, excluding sequence
            // that was transferred backward or forward.
            copy(
                branch.rawSequence.begin() + edge.backwardTransferCount,
                branch.rawSequence.end() - edge.forwardTransferCount,
                back_inserter(branch.gfaSequence));

            // Add the sequence transferred backward by the following bubble, if appropriate.
            if(not edge.isBubble()) {
                if(in_degree(v1, g)==1 and out_degree(v1, g)==1) {
                    out_edge_iterator it;
                    tie(it, ignore) = out_edges(v1, g);
                    const E& nextEdge = g[*it];
                    if(nextEdge.isBubble()) {
                        const vector<Base>& s = nextEdge.branches.front().rawSequence;
                        copy(s.begin(), s.begin() + nextEdge.backwardTransferCount,
                            back_inserter(branch.gfaSequence));
                    }
                }
            }
        }
    }

    performanceLog << timestamp << "storeGfaSequence ends." << endl;
}



// Return the number of raw bases of sequence identical between
// all branches at the beginning.
uint64_t AssemblyGraph2Edge::countCommonPrefixBases() const
{
    SHASTA_ASSERT(isBubble());

    const vector<Base>& firstRawSequence = branches.front().rawSequence;

    for(uint64_t position=0; position<firstRawSequence.size(); position++) {
        const Base base = firstRawSequence[position];

        for(uint64_t branchId=1; branchId<branches.size(); branchId++) {
            const vector<Base>& rawSequence = branches[branchId].rawSequence;
            if(position == rawSequence.size()) {
                return position;
            }
            if(rawSequence[position] != base) {
                return position;
            }
        }
    }

    return firstRawSequence.size();
}



// Return the number of raw bases of sequence identical between
// all branches at the end.
uint64_t AssemblyGraph2Edge::countCommonSuffixBases() const
{
    SHASTA_ASSERT(isBubble());

    const vector<Base>& firstRawSequence = branches.front().rawSequence;

    for(uint64_t position=0; position<firstRawSequence.size(); position++) {
        const Base base = firstRawSequence[firstRawSequence.size() - 1 - position];

        for(uint64_t branchId=1; branchId<branches.size(); branchId++) {
            const vector<Base>& rawSequence = branches[branchId].rawSequence;
            if(position == rawSequence.size()) {
                return position;
            }
            if(rawSequence[rawSequence.size() - 1 - position] != base) {
                return position;
            }
        }
    }

    return firstRawSequence.size();
}




// Figure out if this is a bubble is caused by copy number
// differences in repeats of period up to maxPeriod.
// If this is the case, returns the shortest period for which this is true.
// Otherwise, returns 0.
void AssemblyGraph2Edge::computeCopyNumberDifferencePeriod(uint64_t maxPeriod)
{
    if(not isBubble()) {
        period = 0;
    }

    // Check all pairs of branches.
    vector<uint64_t> periods;
    for(uint64_t i=0; i<branches.size()-1; i++) {
        const vector<Base>& iSequence = branches[i].rawSequence;
        for(uint64_t j=i+1; j<branches.size(); j++) {
            const vector<Base>& jSequence = branches[j].rawSequence;
            const uint64_t pairPeriod = shasta::isCopyNumberDifference(iSequence, jSequence, maxPeriod);
            if(pairPeriod == 0) {
                period = 0;
                return;
            }
            periods.push_back(pairPeriod);
        }
    }
    deduplicate(periods);


    if(periods.size() == 1) {
        period = periods.front();
    } else {
        period = 0;
    }
}




// Fill in orientedReads and average/minimum coverage.
void AssemblyGraph2Edge::Branch::storeReadInformation(const MarkerGraph& markerGraph)
{
    minimumCoverage = std::numeric_limits<uint64_t>::max();
    coverageSum = 0;
    orientedReadIds.clear();

    // Loop over the marker graph path of this branch.
    for(MarkerGraph::EdgeId edgeId: path) {

        // Loop over the marker intervals of this edge.
        const span<const MarkerInterval> markerIntervals =
            markerGraph.edgeMarkerIntervals[edgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            orientedReadIds.push_back(markerInterval.orientedReadId);
        }

        // Update coverage.
        minimumCoverage = min(minimumCoverage, uint64_t(markerIntervals.size()));
        coverageSum += markerIntervals.size();
    }

    deduplicate(orientedReadIds);
}



// Store read information on all edges.
void AssemblyGraph2::storeReadInformation()
{
    performanceLog << timestamp << "AssemblyGraph2::storeReadInformation begins." << endl;

    G& g = *this;

    BGL_FORALL_EDGES(e, g, G) {
        g[e].storeReadInformation(markerGraph);
    }
    performanceLog << timestamp << "AssemblyGraph2::storeReadInformation ends." << endl;
}



// Store read information on all edges. Multithreaded version.
void AssemblyGraph2::storeReadInformationParallel(uint64_t threadCount)
{
    performanceLog << timestamp << "AssemblyGraph2::storeReadInformationParallel begins." << endl;
    G& g = *this;

    // Store a vector of edge descriptors for all edges, to be processed in parallel.
    storeReadInformationParallelData.allEdges.clear();
    BGL_FORALL_EDGES(e, g, G) {
        storeReadInformationParallelData.allEdges.push_back(e);
    }

    // Process all edges in parallel.
    const uint64_t batchSize = 100;
    setupLoadBalancing(storeReadInformationParallelData.allEdges.size(), batchSize);
    runThreads(&AssemblyGraph2::storeReadInformationThreadFunction, threadCount);
    storeReadInformationParallelData.allEdges.clear();

    performanceLog << timestamp << "AssemblyGraph2::storeReadInformationParallel ends." << endl;
}



void AssemblyGraph2::storeReadInformationThreadFunction(size_t threadId)
{
    G& g = *this;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all edges in this batch.
        for(uint64_t i=begin; i!=end; i++) {
            const edge_descriptor e = storeReadInformationParallelData.allEdges[i];
            g[e].storeReadInformation(markerGraph);
        }
    }
}



// Store read information on all branches.
void AssemblyGraph2Edge::storeReadInformation(const MarkerGraph& markerGraph)
{
    for(Branch& branch: branches) {
        branch.storeReadInformation(markerGraph);
    }
}



uint64_t AssemblyGraph2Edge::getStrongestBranchId() const
{
    SHASTA_ASSERT(not branches.empty());
    uint64_t strongestBranchId = 0;

    uint64_t strongestBranchCoverage = branches.front().averageCoverage();

    for(uint64_t branchId=0; branchId<branches.size(); branchId++) {
        const uint64_t coverage = branches[branchId].averageCoverage();
        if (coverage > strongestBranchCoverage) {
            strongestBranchId = branchId;
            strongestBranchCoverage = coverage;
        }
    }

    return strongestBranchId;
}



// This is currently not used.
void AssemblyGraph2::removeSecondaryBubbles(uint64_t secondaryEdgeCleanupThreshold)
{
    G& g = *this;

    // Loop over edges.
    BGL_FORALL_EDGES(e, g, G) {

        // If not a bubble, do nothing.
        E& edge = g[e];
        if(not edge.isBubble()) {
            continue;
        }

        // See how many secondary branches we have.
        uint64_t secondaryCount = 0;
        for(const E::Branch& branch: edge.branches) {
            if(branch.containsSecondaryEdges) {
                ++secondaryCount;
            }
        }

        // If there are no secondary branches, do nothing.
        if(secondaryCount == 0) {
            continue;
        }
        const uint64_t primaryCount = edge.ploidy() - secondaryCount;

        // If any branches are long, skip it.
        bool isLong = false;
        for(const E::Branch& branch: edge.branches) {
            if(branch.path.size() > secondaryEdgeCleanupThreshold) {
                isLong = true;
                break;
            }
        }
        if(isLong) {
            continue;
        }


        // Remove secondary branches, keeping at most one.
        if(primaryCount > 0) {

            // There is at least one primary branch, so we can
            // remove all the secondary ones.
            edge.removeAllSecondaryBranches();
        } else {

            // There are no primary branches.
            // Remove all secondary branches except the strongest.
            edge.removeAllBranchesExceptStrongest();
        }
    }

}



void AssemblyGraph2::removeWeakBranches(uint64_t strongBranchThreshold)
{
    G& g = *this;

    // Loop over edges.
    BGL_FORALL_EDGES(e, g, G) {

        // If not a bubble, do nothing.
        E& edge = g[e];
        if(not edge.isBubble()) {
            continue;
        }

        // Find the weak branches.
        std::set<uint64_t> weakBranches;
        for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
            const E::Branch& branch = edge.branches[branchId];

            if( (branchId != edge.getStrongestBranchId()) and
                (uint64_t(branch.orientedReadIds.size()) < strongBranchThreshold)) {
                weakBranches.insert(branchId);
            }
        }

        // Remove them.
        if(not weakBranches.empty()) {
            vector<E::Branch> branches;
            for(uint64_t branchId=0; branchId<edge.ploidy(); branchId++) {
                if(weakBranches.find(branchId) == weakBranches.end()) {
                    const E::Branch& branch = edge.branches[branchId];
                    branches.push_back(branch);
                }
            }
            edge.branches.swap(branches);
        }

    }
}



void AssemblyGraph2Edge::removeAllSecondaryBranches()
{
    vector<Branch> newBranches;
    for(const Branch& branch: branches) {
        if(not branch.containsSecondaryEdges) {
            newBranches.push_back(branch);
        }
    }
    branches.swap(newBranches);
}



void AssemblyGraph2Edge::removeAllBranchesExceptStrongest()
{
    vector<Branch> newBranches(1, branches[getStrongestBranchId()]);
    branches.swap(newBranches);
}



void AssemblyGraph2Edge::forceMaximumPloidy(uint64_t maxPloidy)
{
    // If the ploidy is already not greater than maxPloidy, do nothing.
    if(ploidy() <= maxPloidy) {
        return;
    }

    // Sort the branches by decreasing average coverage.
    vector< pair<uint64_t, uint64_t> > v;  // Pairs(branchId, average coverage).
    for(uint64_t branchId=0; branchId<branches.size(); branchId++) {
        v.push_back(make_pair(branchId, branches[branchId].averageCoverage()));
    }
    sort(v.begin(), v.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());


    // Only keep the maxPloidy strongest branches.
    vector<Branch> newBranches;
    for(uint64_t i=0; i<maxPloidy; i++) {
        newBranches.push_back(branches[v[i].first]);
    }
    branches.swap(newBranches);
}



void AssemblyGraph2::forceMaximumPloidy(uint64_t maxPloidy)
{
    performanceLog << timestamp << "AssemblyGraph2::forceMaximumPloidy begins." << endl;

    G& g = *this;
    BGL_FORALL_EDGES(e, g, G) {
        g[e].forceMaximumPloidy(maxPloidy);
    }

    performanceLog << timestamp << "AssemblyGraph2::forceMaximumPloidy ends." << endl;
}



// Remove degenerate branches.
void AssemblyGraph2::removeDegenerateBranches()
{
    G& g = *this;

    BGL_FORALL_EDGES(e, g, G) {
        E& edge = g[e];

        const uint64_t ploidy = edge.ploidy();

        if(ploidy == 1) {
            continue;
        }

        if(ploidy == 2) {
            if(edge.isDegenerateBubble()) {
                edge.removeAllBranchesExceptStrongest();
            }
            continue;
        }

        // General case of ploidy 3 or more.

        // Gather branches that have the same sequence.
        // Map with:
        // Key = raw sequence.
        // Value = (average edge coverage, branchId).
        std::map<vector<shasta::Base>, vector<uint64_t> > branchMap;
        for(uint64_t branchId=0; branchId<edge.branches.size(); branchId++) {
            const E::Branch& branch = edge.branches[branchId];
            branchMap[branch.rawSequence].push_back(branchId);
        }



        // For each set of branches with identical sequence, only keep the strongest.
        vector<uint64_t> keep;
        for(const auto& p: branchMap) {
            const vector<uint64_t>& branchIds = p.second;
            if(branchIds.size() == 1) {
                // There is only one branch with this sequence. Keep it.
                keep.push_back(branchIds.front());
                continue;
            }

            // Find the one with the most coverage.
            uint64_t bestBranchId = branchIds.front();
            uint64_t bestCoverage = edge.branches[bestBranchId].averageCoverage();
            for(uint64_t branchId: branchIds) {
                const uint64_t branchCoverage = edge.branches[branchId].averageCoverage();
                if(branchCoverage > bestCoverage) {
                    bestBranchId = branchId;
                    bestCoverage = branchCoverage;
                }
            }
            SHASTA_ASSERT(bestBranchId < ploidy);
            keep.push_back(bestBranchId);
        }

        // Only keep the branches we marked to keep.
        vector<E::Branch> newBranches;
        for(uint64_t branchId: keep) {
            newBranches.push_back(edge.branches[branchId]);
        }
        edge.branches.swap(newBranches);
    }
}



void AssemblyGraph2::hetSnpStatistics(
    uint64_t& transitionCount,
    uint64_t& transversionCount,
    uint64_t& nonSnpCount
) const
{
    using shasta::Base;
    const G& g = *this;

    transitionCount = 0;
    transversionCount= 0;
    nonSnpCount = 0;
    BGL_FORALL_EDGES(e, g, G) {
        const E& edge = g[e];

        if(edge.ploidy() != 2) {
            continue;
        }

        if(edge.isBad) {
            continue;
        }

        const auto& s0 = edge.branches[0].gfaSequence;
        const auto& s1 = edge.branches[1].gfaSequence;

        if(s0.size() != 1) {
            ++nonSnpCount;
            continue;
        }
        if(s1.size() != 1) {
            ++nonSnpCount;
            continue;
        }

        const Base b0 = s0.front();
        const Base b1 = s1.front();

        const bool isPurine0 = (b0.value % 2) == 0;
        const bool isPurine1 = (b1.value % 2) == 0;

        const bool isTransition = (isPurine0 == isPurine1);
        if(isTransition) {
            ++transitionCount;
        } else {
            ++transversionCount;
        }
    }
}



// Merge consecutive non-bubbles, when possible.
void AssemblyGraph2::merge(
    bool storeReadInformation,  // If true, store read information for merged edges.
    bool assemble               // If true, assemble merged edges.
    )
{
    performanceLog << timestamp << "AssemblyGraph2::merge begins." << endl;

    // Find linear chains of non-bubbles.
    vector< vector<edge_descriptor> > chains;
    findNonBubbleLinearChains(chains);

    /*
    cout << "Found the following linear chains that will be merged:" << endl;
    for(const auto& chain: chains) {
        cout << "Chain";
        for(const edge_descriptor e: chain) {
            cout << " " << (*this)[e].id;
        }
        cout << endl;
    }
    */

    // Merge each chain.
    for(const vector<edge_descriptor>& chain: chains) {
        merge(chain, storeReadInformation, assemble);
    }


    performanceLog << timestamp << "AssemblyGraph2::merge ends." << endl;
}



AssemblyGraph2::edge_descriptor AssemblyGraph2::merge(
    const vector<edge_descriptor>& chain,
    bool storeReadInformation,  // If true, store read information for merged edges.
    bool assemble               // If true, assemble merged edges.
    )
{
    G& g = *this;

    /*
    cout << "Merging linear chain ";
    for(const edge_descriptor e: chain) {
        cout << " " << g[e].id;
    }
    cout << endl;
    */


    // Sanity checks.

    // Must be more than one edge.
    SHASTA_ASSERT(chain.size() > 1);

    // There must be no bubbles.
    for(const edge_descriptor e: chain) {
        SHASTA_ASSERT(not g[e].isBubble());
    }

    // It must be a path.
    for(uint64_t i=1; i<chain.size(); i++) {
        SHASTA_ASSERT(target(chain[i-1], g) == source(chain[i], g));
    }

    // Check the degrees of internal vertices.
    // It must be a linear chain.
    for(uint64_t i=1; i<chain.size(); i++) {
        const vertex_descriptor v = source(chain[i], g);
        SHASTA_ASSERT(in_degree(v, g) == 1);
        SHASTA_ASSERT(out_degree(v, g) == 1);
    }

    // Create the merged marker graph path.
    MarkerGraphPath newPath;
    bool containsSecondaryEdges = false;
    for(const edge_descriptor e: chain) {
        const AssemblyGraph2Edge::Branch& branch = g[e].branches.front();
        const MarkerGraphPath& path = g[e].branches.front().path;
        newPath.insert(newPath.end(), path.begin(), path.end());
        containsSecondaryEdges = containsSecondaryEdges or branch.containsSecondaryEdges;
    }

    // Add the new edge.
    const edge_descriptor eNew = addEdge(newPath, containsSecondaryEdges, storeReadInformation, assemble);

    // Gather the vertices in between, which will be removed.
    vector<vertex_descriptor> verticesToBeRemoved;
    for(uint64_t i=1; i<chain.size(); i++) {
        verticesToBeRemoved.push_back(source(chain[i], g));
    }

    // Remove the old edges.
    for(const edge_descriptor e: chain) {
        boost::remove_edge(e, g);
    }

    // Remove the vertices in between.
    for(const vertex_descriptor v: verticesToBeRemoved) {
        SHASTA_ASSERT(in_degree(v, g) == 0);
        SHASTA_ASSERT(out_degree(v, g) == 0);
        boost::remove_vertex(v, g);
    }

    return eNew;
}



// Merge an edge with the previous edge, if possible.
AssemblyGraph2::edge_descriptor AssemblyGraph2::mergeWithPreviousIfPossible(edge_descriptor e)
{
    AssemblyGraph2& g = *this;

    // This edge cannot be a bubble.
    if(g[e].isBubble()) {
        return e;
    }

    // The source vertex of e must have in-degree and out-degree 1.
    const vertex_descriptor v0 = source(e, g);
    if(in_degree(v0, g) != 1) {
        return e;
    }
    if(out_degree(v0, g) != 1) {
        return e;
    }

    // Locate the one and only previous edge.
    in_edge_iterator it;
    tie(it, ignore) = in_edges(v0, g);
    const edge_descriptor ePrevious = *it;

    // The previous edge cannot be a bubble.
    if(g[ePrevious].isBubble()) {
        return e;
    }


    // If getting here, we can merge.

    // Create the new edge.
    edge_descriptor eNew;
    tie(eNew, ignore) = add_edge(source(ePrevious, g), target(e, g), E(nextId++), g);
    g[eNew].branches.resize(1);

    // Access the branches we are working with.
    const AssemblyGraph2Edge::Branch& branch = g[e].branches.front();
    const AssemblyGraph2Edge::Branch& previousBranch = g[ePrevious].branches.front();
    AssemblyGraph2Edge::Branch& newBranch = g[eNew].branches.front();

    // Create the combined marker graph path.
    newBranch.path = previousBranch.path;
    copy(branch.path.begin(), branch.path.end(), back_inserter(newBranch.path));
    newBranch.containsSecondaryEdges = branch.containsSecondaryEdges or previousBranch.containsSecondaryEdges;

    // Recompute read support for the merged branch.
    newBranch.storeReadInformation(markerGraph);

    // Compute sequence for the updated edge.
    assemble(eNew);

    // Remove the edges we are merging.
    boost::remove_edge(e, g);
    boost::remove_edge(ePrevious, g);

    // Remove the vertex in between.
    SHASTA_ASSERT(in_degree(v0, g) == 0);
    SHASTA_ASSERT(out_degree(v0, g) == 0);
    remove_vertex(v0, g);

    // Done.
    return eNew;
}



// Merge an edge with the following edge, if possible.
AssemblyGraph2::edge_descriptor AssemblyGraph2::mergeWithFollowingIfPossible(edge_descriptor e)
{
    AssemblyGraph2& g = *this;

    // This edge cannot be a bubble.
    if(g[e].isBubble()) {
        return e;
    }

    // The target vertex of e must have in-degree and out-degree 1.
    const vertex_descriptor v1 = target(e, g);
    if(in_degree(v1, g) != 1) {
        return e;
    }
    if(out_degree(v1, g) != 1) {
        return e;
    }

    // Locate the one and only following edge.
    out_edge_iterator it;
    tie(it, ignore) = out_edges(v1, g);
    const edge_descriptor eFollowing = *it;

    // The following edge cannot be a bubble.
    if(g[eFollowing].isBubble()) {
        return e;
    }


    // If getting here, we can merge.

    // Create the new edge.
    edge_descriptor eNew;
    tie(eNew, ignore) = add_edge(source(e, g), target(eFollowing, g), E(nextId++), g);
    g[eNew].branches.resize(1);

    // Access the branches we are working with.
    const AssemblyGraph2Edge::Branch& branch = g[e].branches.front();
    const AssemblyGraph2Edge::Branch& followingBranch = g[eFollowing].branches.front();
    AssemblyGraph2Edge::Branch& newBranch = g[eNew].branches.front();

    // Create the combined marker graph path.
    newBranch.path = branch.path;
    copy(followingBranch.path.begin(), followingBranch.path.end(), back_inserter(newBranch.path));
    newBranch.containsSecondaryEdges = branch.containsSecondaryEdges or followingBranch.containsSecondaryEdges;

    // Recompute read support for the merged branch.
    newBranch.storeReadInformation(markerGraph);

    // Compute sequence for the updated edge.
    assemble(eNew);

    // Remove the edges we are merging.
    boost::remove_edge(e, g);
    boost::remove_edge(eFollowing, g);

    // Remove the vertex in between.
    SHASTA_ASSERT(in_degree(v1, g) == 0);
    SHASTA_ASSERT(out_degree(v1, g) == 0);
    remove_vertex(v1, g);

    // Done.
    return eNew;
}



// Find linear chains of adjacent non-bubbles.
// Used by merge.
void AssemblyGraph2::findNonBubbleLinearChains(
    vector< vector<edge_descriptor> >& chains) const
{
    const G& g = *this;

    // The edges we have already encountered.
    std::set<edge_descriptor> edgesFound;

    // Start with no chains.
    chains.clear();

    // Work vectors used in the main loop below.
    vector<edge_descriptor> forwardPath;
    vector<edge_descriptor> backwardPath;

    // Consider all possible start edges for the chain.
    BGL_FORALL_EDGES_T(eStart, g, G) {

        // If this is a bubble, skip it.
        if(g[eStart].isBubble()) {
            continue;
        }

        // If we already assigned this edge to a chain, skip it.
        if(edgesFound.find(eStart) != edgesFound.end()) {
            continue;
        }

        edgesFound.insert(eStart);

        // Extend forward.
        forwardPath.clear();
        bool isCircular = false;
        edge_descriptor e = eStart;
        while(true) {
            const vertex_descriptor v = target(e, g);
            if(in_degree(v, g) != 1) {
                break;
            }
            if(out_degree(v, g) != 1) {
                break;
            }
            BGL_FORALL_OUTEDGES_T(v, eNext, g, G) {
                e = eNext;
                break;
            }
            if(g[e].isBubble()) {
                break;
            }
            if(e == eStart) {
                isCircular = true;
                break;
            }
            forwardPath.push_back(e);
            SHASTA_ASSERT(edgesFound.find(e) == edgesFound.end());
            edgesFound.insert(e);
        }


        // Extend backward.
        backwardPath.clear();
        if(not isCircular) {
            edge_descriptor e = eStart;
            while(true) {
                const vertex_descriptor v = source(e, g);
                if(in_degree(v, g) != 1) {
                    break;
                }
                if(out_degree(v, g) != 1) {
                    break;
                }
                BGL_FORALL_INEDGES_T(v, ePrevious, g, G) {
                    e = ePrevious;
                    break;
                }
                if(g[e].isBubble()) {
                    break;
                }
                if(e == eStart) {
                    isCircular = true;
                    break;
                }
                backwardPath.push_back(e);
                SHASTA_ASSERT(edgesFound.find(e) == edgesFound.end());
                edgesFound.insert(e);
            }
        }

        // If the chain is too short, get rid of it.
        if(backwardPath.size() + 1 + forwardPath.size() < 2) {
            continue;
        }

        // Store it.
        chains.resize(chains.size() + 1);
        vector<edge_descriptor>& chain = chains.back();
        copy(backwardPath.rbegin(), backwardPath.rend(), back_inserter(chain));
        chain.push_back(eStart);
        copy(forwardPath.begin(), forwardPath.end(), back_inserter(chain));

    }

    // Check that all non-bubble edges were found.
    BGL_FORALL_EDGES_T(e, g, G) {
        if(not g[e].isBubble()) {
            SHASTA_ASSERT(edgesFound.find(e) != edgesFound.end());
        }
    }

}



void AssemblyGraph2::findBubbleChains()
{
    G& g = *this;

    vector< vector<edge_descriptor> > linearChains;
    findLinearChains(*this, 2, linearChains);
    bubbleChains.clear();
    bubbleChains.resize(linearChains.size());
    for(uint64_t i=0; i<linearChains.size(); i++) {
        BubbleChain& bubbleChain = bubbleChains[i];
        bubbleChain.edges.swap(linearChains[i]);
    }

    /*
    cout << "Found " << bubbleChains.size() << " bubble chains with the following numbers of edges:";
    for(const auto& bubbleChain: bubbleChains) {
        cout << " " << bubbleChain.edges.size();
    }
    cout << endl;
    */


    // Store pointers in the begin/end vertices.
    BGL_FORALL_VERTICES(v, g, G) {
        V& vertex = g[v];
        vertex.bubbleChainsBeginningHere.clear();
        vertex.bubbleChainsEndingHere.clear();
    }
    for(const BubbleChain& bubbleChain: bubbleChains) {
        SHASTA_ASSERT(not bubbleChain.edges.empty());
        const vertex_descriptor vBegin = source(bubbleChain.edges.front(), g);
        const vertex_descriptor vEnd = target(bubbleChain.edges.back(), g);
        g[vBegin].bubbleChainsBeginningHere.push_back(&bubbleChain);
        g[vEnd].bubbleChainsEndingHere.push_back(&bubbleChain);
    }


    // Stores pointers in edges.
    BGL_FORALL_EDGES(e, g, G) {
        g[e].bubbleChain = {0, 0};
    }
    for(const BubbleChain& bubbleChain: bubbleChains) {
        for(uint64_t position=0; position<bubbleChain.edges.size(); position++) {
            const edge_descriptor e = bubbleChain.edges[position];
            g[e].bubbleChain = make_pair(&bubbleChain, position);
        }
    }
}



void AssemblyGraph2::clearBubbleChains()
{
    G& g = *this;

    bubbleChains.clear();

    BGL_FORALL_VERTICES(v, g, G) {
        V& vertex = g[v];
        vertex.bubbleChainsBeginningHere.clear();
        vertex.bubbleChainsEndingHere.clear();
    }

    BGL_FORALL_EDGES(e, g, G) {
        g[e].bubbleChain = {0, 0};
    }
}



// Find PhasingRegions within BubbleChains.
void AssemblyGraph2::findPhasingRegions()
{
    for(BubbleChain& bubbleChain: bubbleChains) {
        findPhasingRegions(bubbleChain);
    }
}



void AssemblyGraph2::findPhasingRegions(BubbleChain& bubbleChain)
{
    const G& g = *this;
    const auto& edges = bubbleChain.edges;

    // cout << "findPhasingRegions begins for bubble chain " << bubbleChain.id << endl;

    // Gather all the positions that have a defined componentId.
    vector< pair<uint64_t, uint64_t> > bubbleChainTable;
    for(uint64_t position=0; position<edges.size(); position++) {
        const AssemblyGraph2Edge& edge = g[edges[position]];
        const uint64_t componentId = edge.componentId;
        if(componentId != std::numeric_limits<uint64_t>::max()) {
            bubbleChainTable.push_back(make_pair(position, componentId));
        }
    }

    // Use this table to find the boundaries of the phased regions.
    vector<uint64_t> firstPositions;
    vector<uint64_t> lastPositions;
    for(uint64_t i=0; i<bubbleChainTable.size(); i++) {
        const auto& p = bubbleChainTable[i];
        const uint64_t position = p.first;
        const uint64_t componentId = p.second;
        if(i==0 or componentId != bubbleChainTable[i-1].second) {
            firstPositions.push_back(position);
        }
        if(i==bubbleChainTable.size()-1 or componentId != bubbleChainTable[i+1].second) {
            lastPositions.push_back(position);
        }
    }



    // Now we can create the phased regions.
    bubbleChain.phasingRegions.clear();



    // If nothing was phased, generate a single phasing region for the entire bubble chain.
    SHASTA_ASSERT(firstPositions.size() == lastPositions.size());
    if(firstPositions.empty()) {
        BubbleChain::PhasingRegion unphasedRegion;
        unphasedRegion.firstPosition = 0;
        unphasedRegion.lastPosition = edges.size() - 1;
        unphasedRegion.isPhased = false;
        bubbleChain.phasingRegions.push_back(unphasedRegion);
        return;
    }



    // Create an initial unphased region, if necessary.
    if(firstPositions.front() != 0) {
        BubbleChain::PhasingRegion unphasedRegion;
        unphasedRegion.firstPosition = 0;
        unphasedRegion.lastPosition = firstPositions.front() - 1;
        unphasedRegion.isPhased = false;
        bubbleChain.phasingRegions.push_back(unphasedRegion);
    }

    for(uint64_t i=0; i<firstPositions.size(); i++) {

        // Add a phased region.
        BubbleChain::PhasingRegion phasedRegion;
        phasedRegion.firstPosition = firstPositions[i];
        phasedRegion.lastPosition = lastPositions[i];
        phasedRegion.isPhased = true;
        phasedRegion.componentId = g[edges[firstPositions[i]]].componentId;
        bubbleChain.phasingRegions.push_back(phasedRegion);

        // If necessary, add an unphased region to bridge the gap to the
        // next phased region.
        if(i != firstPositions.size() - 1 ) {
            if(firstPositions[i + 1] != lastPositions[i] + 1) {
                BubbleChain::PhasingRegion unphasedRegion;
                unphasedRegion.firstPosition = lastPositions[i] + 1;
                unphasedRegion.lastPosition = firstPositions[i + 1] - 1;
                unphasedRegion.isPhased = false;
                bubbleChain.phasingRegions.push_back(unphasedRegion);
            }
        }
    }

    // Create a final unphased region, if necessary.
    if(lastPositions.back() != edges.size()-1) {
        BubbleChain::PhasingRegion unphasedRegion;
        unphasedRegion.firstPosition = lastPositions.back() + 1;
        unphasedRegion.lastPosition = edges.size()-1;
        unphasedRegion.isPhased = false;
        bubbleChain.phasingRegions.push_back(unphasedRegion);
    }
}



void AssemblyGraph2::writePhasingRegions()
{
    performanceLog << timestamp << "AssemblyGraph2::writePhasingRegions begins." << endl;

    ofstream csv("PhasingRegions.csv");
    csv << "Bubble chain id,Phasing region id,First position,Last position,Phased,Component,\n";

    for(uint64_t bubbleChainId=0; bubbleChainId<uint64_t(bubbleChains.size()); bubbleChainId++) {
        const BubbleChain& bubbleChain = bubbleChains[bubbleChainId];
        for(uint64_t phasingRegionId=0; phasingRegionId<uint64_t(bubbleChain.phasingRegions.size()); phasingRegionId++) {
            const auto& phasingRegion = bubbleChain.phasingRegions[phasingRegionId];
            csv <<
                bubbleChainId << "," <<
                phasingRegionId << "," <<
                phasingRegion.firstPosition << "," <<
                phasingRegion.lastPosition << ",";
            if(phasingRegion.isPhased) {
                csv << "Yes," << phasingRegion.componentId << ",";
            } else {
                csv << "No,,";
            }
            csv << "\n";
        }
    }

    performanceLog << timestamp << "AssemblyGraph2::writePhasingRegions ends." << endl;
}



void AssemblyGraph2::writeBubbleChains()
{
    performanceLog << timestamp << "AssemblyGraph2::writeBubbleChains begins." << endl;
    G& g = *this;

    ofstream csv("BubbleChains.csv");
    csv << "Bubble chain,Position,Edge,Ploidy,Component,\n";

    for(uint64_t bubbleChainId=0; bubbleChainId<uint64_t(bubbleChains.size()); bubbleChainId++) {
        const BubbleChain& bubbleChain = bubbleChains[bubbleChainId];
        const vector<edge_descriptor>& edges = bubbleChain.edges;

        for(uint64_t position=0; position<uint64_t(edges.size()); position++) {
            const edge_descriptor e = edges[position];
            const AssemblyGraph2Edge& edge = g[e];

            csv << bubbleChainId << ",";
            csv << position << ",";
            csv << edge.id << ",";
            csv << edge.ploidy() << ",";
            if(edge.componentId != std::numeric_limits<uint64_t>::max()){
                csv << edge.componentId;
            }
            csv << ",";

            csv << "\n";

        }
    }

    performanceLog << timestamp << "AssemblyGraph2::writeBubbleChains ends." << endl;
}



void AssemblyGraph2::handleSuperbubbles0(
    uint64_t edgeLengthThreshold,
    uint64_t maxSuperbubbleSize,
    uint64_t maxSuperbubbleChunkSize,
    uint64_t maxSuperbubbleChunkPathCount,
    bool storeReadInformation,  // If true, store read information for newly created edges.
    bool assemble               // If true, assemble sequence for newly created edges
    )
{
    G& g = *this;
    performanceLog << timestamp << "AssemblyGraph2::handleSuperbubbles0 begins." << endl;

    // We cannot use boost::connected_components because it
    // only works for undirected graphs.

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, g, G) {
        vertexMap.insert(make_pair(v, vertexIndex++));
    }
    const uint64_t n = uint64_t(vertexMap.size());

    // Initialize the disjoint sets data structure.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint32_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }
    // Main loop over short edges.
    BGL_FORALL_EDGES(e, g, G) {
        const E& edge = g[e];
        if(edge.maximumPathLength() <= edgeLengthThreshold) {
            const vertex_descriptor v0 = source(e, g);
            const vertex_descriptor v1 = target(e, g);
            disjointSets.union_set(vertexMap[v0], vertexMap[v1]);
        }
    }

    // Gather the vertices in each connected component.
    std::map<uint64_t, vector<vertex_descriptor> > components;
    BGL_FORALL_VERTICES(v, g, G) {
        const uint64_t component = disjointSets.find_set(vertexMap[v]);
        components[component].push_back(v);
    }



    // Process one components one at a time.
    // Each component is used to create a superbubble.
    for(const auto& p: components) {
        const vector<vertex_descriptor>& componentVertices = p.second;

        // Create a superbubble with this component.
        Superbubble superbubble(g, componentVertices, edgeLengthThreshold);

        // Process it.
        handleSuperbubble1(superbubble,
            maxSuperbubbleSize, maxSuperbubbleChunkSize, maxSuperbubbleChunkPathCount,
            storeReadInformation, assemble);
    }
    performanceLog << timestamp << "AssemblyGraph2::handleSuperbubbles0 ends." << endl;
}



// This creates superbubbles using all edges not in bubble chains.
void AssemblyGraph2::handleSuperbubbles1(
    uint64_t maxSuperbubbleSize,
    uint64_t maxSuperbubbleChunkSize,
    uint64_t maxSuperbubbleChunkPathCount,
    bool storeReadInformation,  // If true, store read information for newly created edges.
    bool assemble               // If true, assemble sequence for newly created edges
    )
{
    G& g = *this;
    performanceLog << timestamp << "AssemblyGraph2::handleSuperbubbles1 begins." << endl;

    findBubbleChains();

    // We cannot use boost::connected_components because it
    // only works for undirected graphs.

    // Map vertices to integers.
    std::map<vertex_descriptor, uint64_t> vertexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, g, G) {
        vertexMap.insert(make_pair(v, vertexIndex++));
    }
    const uint64_t n = uint64_t(vertexMap.size());

    // Initialize the disjoint sets data structure.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint32_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    // Main loop over edges that don't belong to bubble chains.
    BGL_FORALL_EDGES(e, g, G) {
        const E& edge = g[e];
        if(edge.bubbleChain.first == 0) {
            const vertex_descriptor v0 = source(e, g);
            const vertex_descriptor v1 = target(e, g);
            disjointSets.union_set(vertexMap[v0], vertexMap[v1]);
        }
    }

    // Gather the vertices in each connected component.
    std::map<uint64_t, vector<vertex_descriptor> > components;
    BGL_FORALL_VERTICES(v, g, G) {
        const uint64_t component = disjointSets.find_set(vertexMap[v]);
        components[component].push_back(v);
    }

    // Process one components one at a time.
    // Each component is used to create a superbubble.
    for(const auto& p: components) {
        const vector<vertex_descriptor>& componentVertices = p.second;

        // Create a superbubble with this component.
        Superbubble superbubble(g, componentVertices);
        // superbubble.writeGraphviz(cout, g);

        // Process it.
        handleSuperbubble1(superbubble,
            maxSuperbubbleSize, maxSuperbubbleChunkSize, maxSuperbubbleChunkPathCount,
            storeReadInformation, assemble);
    }

    clearBubbleChains();
    performanceLog << timestamp << "AssemblyGraph2::handleSuperbubbles1 ends." << endl;
}



/*******************************************************************************

Version of superbubble removal that avoids enumerating paths over the entire superbubble.

For a superbubble with exactly one entrance and one exit, this uses
a dominator tree to divide the superbubble in chunks, and then
does path enumeration over each chunk separately.

Some nomenclature:

- On the dominator tree, there is only one path between the entrance
  and the exit. We call this the critical path. This does not necessarily
  corresponds to a path in the superbubble.

- The vertices on the critical path are called the choke points.
  They are numbered consecutively starting at 0.
  The entrance is choke point 0.

- Vertices that are not in the dominator tree are unreachable from the
  entrance. We call these the unreachable vertices.
  The remaining vertices are the reachable vertices.

- Superbubble edges for which the source vertex is reachable are called reachable.
  Superbubble edges for which the source vertex is unreachable are called unreachable.

- From each reachable vertex, we can follow the dominator tree up
  until we encounter the first choke point. This choke point
  is called the parent choke point of the vertex that we started from.

- For each reachable edge, the parent choke point of the source vertex of
  the edge is called the parent choke point of the edge.

- A chunk is the set of all reachable edges with the same parent choke point.
  That common parent choke point is the called the source of the chunk.
  The next choke point on the critical path is called the target of the chunk.

- A chunk is called trivial if all of its edges have as source the source
  of the chunk and as target the target of the chunk.

With these definitions, a superbubble with exactly one entrance and
one exit is processed as follows:

- The AssemblyGraph2 edges corresponding to all unreachable edges
  are removed from the AssemblyGraph2.

- For trivial chunks, no processing takes place.

- For non-trivial chunks:
  * We do path enumeration from the source to the target of the chunk
    and keep only the two strongest paths. These are used
    to generate a new bubble in the AssemblyGraph2.
  * All AssemblyGraph2 edges corresponding to edges in the chunk are removed.

*******************************************************************************/

void AssemblyGraph2::handleSuperbubble1(
    Superbubble& superbubble,
    uint64_t maxSuperbubbleSize,
    uint64_t maxSuperbubbleChunkSize,
    uint64_t maxSuperbubbleChunkPathCount,
    bool storeReadInformation,  // If true, store read information for newly created edges.
    bool assemble               // If true, assemble sequence for newly created edges
    )
{
    G& g = *this;
    const bool debug = false;

    // If there are no edges, don't do anything.
    if(num_edges(superbubble) == 0) {
        return;
    }

    // If just a simple linear chain, don't do anything.
    if(superbubble.isSimpleLinearChain()) {
        return;
    }

    if(debug) {
        cout << "Processing a non-trivial superbubble with " <<
            superbubble.entrances.size() << " entrances, " <<
            superbubble.exits.size() << " exits, " <<
            num_vertices(superbubble) << " vertices, and " << num_edges(superbubble) << " edges:\n";
        superbubble.writeGraphviz(cout, g);
    }



    // Ignore superbubbles that don't have exactly one entrance and one exit.
    if((superbubble.entrances.size() != 1) or (superbubble.exits.size() != 1)) {
        if(debug) {
            cout << "Superbubble ignored because does not have exactly one entrance and one exit." << endl;
        }
        return;
    }


    // If the superbubble is too big, ignore it.
    if(num_vertices(superbubble) > maxSuperbubbleSize) {
        if(debug) {
            cout << "Superbubble ignored because it is too big." << endl;
        }
        return;
    }

    const Superbubble::vertex_descriptor entrance = superbubble.entrances.front();
    const Superbubble::vertex_descriptor exit = superbubble.exits.front();

    // Compute the forward dominator tree.
    // Use the shasta version which includes a bug fix
    // (see dominatorTree.hpp).
    shasta::lengauer_tarjan_dominator_tree(
        superbubble,
        entrance,
        boost::get(&SuperbubbleVertex::immediateDominator0, superbubble));

    // Compute the backward dominator tree.
    shasta::lengauer_tarjan_dominator_tree(
        boost::reverse_graph<Superbubble>(superbubble),
        exit,
        boost::get(&SuperbubbleVertex::immediateDominator1, superbubble));



    if(debug) {
        cout << "Forward dominator tree:" << endl;
        BGL_FORALL_VERTICES (sv0, superbubble, Superbubble) {
            const AssemblyGraph2::vertex_descriptor av0 = superbubble[sv0].av;
            const AssemblyGraph2Vertex& aVertex0 = g[av0];

            cout << aVertex0.markerGraphVertexId;

            if(sv0 == entrance) {
                cout << " entrance" << endl;
            } else {

                const Superbubble::vertex_descriptor sv1 = superbubble[sv0].immediateDominator0;
                if(sv1 == Superbubble::null_vertex()) {
                    cout << " unreachable" << endl;
                } else {
                    const AssemblyGraph2::vertex_descriptor av1 = superbubble[sv1].av;
                    const AssemblyGraph2Vertex& aVertex1 = g[av1];
                    cout << " parent is " << aVertex1.markerGraphVertexId << endl;
                }
            }
        }



        cout << "Backward dominator tree:" << endl;
        BGL_FORALL_VERTICES (sv0, superbubble, Superbubble) {
            const AssemblyGraph2::vertex_descriptor av0 = superbubble[sv0].av;
            const AssemblyGraph2Vertex& aVertex0 = g[av0];

            cout << aVertex0.markerGraphVertexId;

            if(sv0 == exit) {
                cout << " exit" << endl;
            } else {

                const Superbubble::vertex_descriptor sv1 = superbubble[sv0].immediateDominator1;
                if(sv1 == Superbubble::null_vertex()) {
                    cout << " unreachable" << endl;
                } else {
                    const AssemblyGraph2::vertex_descriptor av1 = superbubble[sv1].av;
                    const AssemblyGraph2Vertex& aVertex1 = g[av1];
                    cout << " parent is " << aVertex1.markerGraphVertexId << endl;
                }
            }
        }
    }


    // In the exceptional case that the exit is unreachable from entrance,
    // do nothing.
    if(
        superbubble[exit].immediateDominator0 == Superbubble::null_vertex() or
        superbubble[entrance].immediateDominator1 == Superbubble::null_vertex()) {
        return;
    }


    // Construct the critical path.
    superbubble.computeCriticalPath();
    if(debug) {
        cout << "Critical path:" << endl;
        for(const auto sv: superbubble.criticalPath) {
            const AssemblyGraph2::vertex_descriptor av = superbubble[sv].av;
            const AssemblyGraph2Vertex& aVertex = g[av];
            cout << aVertex.markerGraphVertexId << endl;
        }
    }

    // Assign edges to chunks.
    superbubble.findChunks();
    if(debug) {
        for(uint64_t chunk=0; chunk<superbubble.chunkEdges.size(); chunk++) {
            cout << "Chunk " << chunk;
            for(const Superbubble::edge_descriptor se: superbubble.chunkEdges[chunk]) {
                const SuperbubbleEdge& sEdge = superbubble[se];
                cout << " " << g[sEdge.ae].pathId(sEdge.branchId);
            }
            cout << endl;
        }
    }



    // Remove edges not assigned to a chunk from both the Superbubble and the AssemblyGraph2.
    // These edges cannot belong to any path between the entrance and the exit.
    vector<Superbubble::edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, superbubble, Superbubble) {
        const SuperbubbleEdge& edge = superbubble[e];
        if(edge.chunk == std::numeric_limits<uint64_t>::max()) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const Superbubble::edge_descriptor se: edgesToBeRemoved) {
        const SuperbubbleEdge& sEdge = superbubble[se];
        if(debug) {
            cout << "Removing edge " << g[sEdge.ae].pathId(sEdge.branchId) << endl;
        }
        if(sEdge.branchId==0) {
            boost::remove_edge(sEdge.ae, g);
        }
        boost::remove_edge(se, superbubble);
    }



    // Loop over the chunks.
    // Chunk chunkId consists of all edges reachable forward from choke point chunkId
    // and backward from choke point chunkId+1.
    for(uint64_t chunkId=0; chunkId<superbubble.criticalPath.size()-1; chunkId++) {
        const Superbubble::vertex_descriptor chunkEntrance = superbubble.criticalPath[chunkId];
        const Superbubble::vertex_descriptor chunkExit = superbubble.criticalPath[chunkId+1];
        if(debug) {
            cout << "Working on chunk " << chunkId << endl;
            cout << "Chunk entrance " << g[superbubble[chunkEntrance].av].markerGraphVertexId << endl;
            cout << "Chunk exit " << g[superbubble[chunkExit].av].markerGraphVertexId << endl;
        }


        // If this is a trivial chunk, skip it.
        // A chunk is trivial if all out-edges of chunkEntrance have chunkExit
        // as their target vertex.
        bool isNonTrivial = false;
        BGL_FORALL_OUTEDGES(chunkEntrance, e, superbubble, Superbubble) {
            if(target(e, superbubble) != chunkExit) {
                isNonTrivial = true;
                break;
            }
        }
        if(not isNonTrivial) {
            if(debug) {
                cout << "This chunk is trivial. Nothing done." << endl;
            }
            continue;
        }



        // If getting here, we have a non-trivial chunk.

        // If the chunk is too big, ignore it.
        if(superbubble.chunkEdges[chunkId].size() > maxSuperbubbleChunkSize) {
            continue;
        }

        // At this stage, read support has not yet been computed.
        // So let's compute it for the edges in this chunk.
        for(const Superbubble::edge_descriptor se: superbubble.chunkEdges[chunkId]) {
            const SuperbubbleEdge& sEdge = superbubble[se];
            const AssemblyGraph2::edge_descriptor ae = sEdge.ae;
            AssemblyGraph2Edge& aEdge = g[ae];
            AssemblyGraph2Edge::Branch& branch = aEdge.branches[sEdge.branchId];
            branch.storeReadInformation(markerGraph);
        }

        // Enumerate paths between chunkEntrance and chunkExit.
        superbubble.enumeratePaths(chunkEntrance, chunkExit);

        // If we found too many paths, ignore this chunk.
        if(superbubble.paths.size() > maxSuperbubbleChunkPathCount) {
            if(debug) {
                cout << "Chunk ignored because it has too many paths." << endl;
            }
            continue;
        }

        if(debug) {
            cout << "Found " << superbubble.paths.size() << " paths for this chunk:" << endl;
            for(const vector<Superbubble::edge_descriptor>& path: superbubble.paths) {

                uint64_t coverageSum = 0;
                uint64_t lengthSum = 0;

                for(const Superbubble::edge_descriptor se: path) {
                    const SuperbubbleEdge& sEdge = superbubble[se];
                    const AssemblyGraph2::edge_descriptor ae = sEdge.ae;
                    const uint64_t branchId = sEdge.branchId;
                    const auto& branch = g[ae].branches[branchId];
                    cout << " " << g[ae].pathId(branchId);

                    coverageSum += branch.coverageSum;
                    lengthSum += branch.path.size();
                }

                cout << " average coverage " << double(coverageSum) / double(lengthSum);
                cout << endl;
            }
        }



        // Compute average coverage for each path.
        SHASTA_ASSERT(superbubble.paths.size() > 1);
        vector< pair<uint64_t, double> > pathCoverageTable;
        for(uint64_t i=0 ; i<superbubble.paths.size(); i++) {
            const auto& path = superbubble.paths[i];
            uint64_t coverageSum = 0;
            uint64_t lengthSum = 0;

            for(const Superbubble::edge_descriptor se: path) {
                const SuperbubbleEdge& sEdge = superbubble[se];
                const AssemblyGraph2::edge_descriptor ae = sEdge.ae;
                const uint64_t branchId = sEdge.branchId;
                const auto& branch = g[ae].branches[branchId];
                coverageSum += branch.coverageSum;
                lengthSum += branch.path.size();
            }
            const double averageCoverage = double(coverageSum) / double(lengthSum);
            pathCoverageTable.push_back(make_pair(i, averageCoverage));
        }
        sort(pathCoverageTable.begin(), pathCoverageTable.end(),
            OrderPairsBySecondOnlyGreater<uint64_t, double>());
        const array<vector<Superbubble::edge_descriptor>, 2> bestPaths = {
            superbubble.paths[pathCoverageTable[0].first],
            superbubble.paths[pathCoverageTable[1].first]};

        if(debug) {
            cout << "Best paths for this chunk:" << endl;
            for(const vector<Superbubble::edge_descriptor>& path: bestPaths) {

                for(const Superbubble::edge_descriptor se: path) {
                    const SuperbubbleEdge& sEdge = superbubble[se];
                    const AssemblyGraph2::edge_descriptor ae = sEdge.ae;
                    const uint64_t branchId = sEdge.branchId;
                    cout << " " << g[ae].pathId(branchId);

                }
                cout << endl;
            }
        }



        // The two best paths could have a common portion at their begin or end.
        // Find the length of the common portions.
        const uint64_t prefixLength = commonPrefixLength(bestPaths[0], bestPaths[1]);
        const uint64_t suffixLength = commonSuffixLength(bestPaths[0], bestPaths[1]);
        if(debug) {
            if(prefixLength) {
                cout << "The two best path have a common prefix of length " << prefixLength << endl;
            }
            if(suffixLength) {
                cout << "The two best path have a common suffix of length " << suffixLength << endl;
            }
        }



        // If there is a common prefix, generate a new haploid edge of the AssemblyGraph2.
        if(prefixLength) {
            const auto begin = bestPaths[0].begin();
            const auto end = begin + prefixLength;

            // Construct the marker graph path.
            MarkerGraphPath markerGraphPath;
            bool containsSecondaryEdges = false;
            for(auto it=begin; it!=end; ++it) {
                const Superbubble::edge_descriptor se = *it;
                const SuperbubbleEdge& sEdge = superbubble[se];
                const AssemblyGraph2::edge_descriptor ae = sEdge.ae;
                const AssemblyGraph2Edge& aEdge = g[ae];
                const AssemblyGraph2Edge::Branch& branch = aEdge.branches[sEdge.branchId];
                copy(branch.path.begin(), branch.path.end(), back_inserter(markerGraphPath));
                if(branch.containsSecondaryEdges) {
                    containsSecondaryEdges = true;
                }
            }

            // Create a new haploid edge with this path.
            addEdge(markerGraphPath, containsSecondaryEdges, storeReadInformation, assemble);
        }



        // Create a new AssemblyGraph2 edge to represent a bubble with the
        // two best paths, excluding their common prefix and suffix.
        if(
            (prefixLength + suffixLength < bestPaths[0].size()) and
            (prefixLength + suffixLength < bestPaths[1].size())) {
            const auto begin0 = bestPaths[0].begin() + prefixLength;
            const auto end0 = bestPaths[0].end() - suffixLength;
            const auto begin1 = bestPaths[1].begin() + prefixLength;
            const auto end1 = bestPaths[1].end() - suffixLength;

            const Superbubble::edge_descriptor first0 = *begin0;
            const Superbubble::edge_descriptor last0 = *(end0 - 1);
            const Superbubble::edge_descriptor first1 = *begin1;
            const Superbubble::edge_descriptor last1 = *(end1 - 1);

            // Find the source and target vertices.
            const Superbubble::vertex_descriptor sv0 = source(first0, superbubble);
            SHASTA_ASSERT(sv0 == source(first1, superbubble));
            const Superbubble::vertex_descriptor sv1 = target(last0, superbubble);
            SHASTA_ASSERT(sv1 == target(last1, superbubble));
            const AssemblyGraph2::vertex_descriptor av0 = superbubble[sv0].av;
            const AssemblyGraph2::vertex_descriptor av1 = superbubble[sv1].av;


            // Create the new edge with two branches, one for each of our best paths.
            MarkerGraphPath markerGraphPath0;
            bool containsSecondaryEdges0 = false;
            for(auto it=begin0; it!=end0; ++it) {
                const Superbubble::edge_descriptor se = *it;
                const SuperbubbleEdge& sEdge = superbubble[se];
                const AssemblyGraph2::edge_descriptor ae = sEdge.ae;
                const AssemblyGraph2Edge& aEdge = g[ae];
                const AssemblyGraph2Edge::Branch& branch = aEdge.branches[sEdge.branchId];
                copy(branch.path.begin(), branch.path.end(), back_inserter(markerGraphPath0));
                if(branch.containsSecondaryEdges) {
                    containsSecondaryEdges0 = true;
                }
            }
            MarkerGraphPath markerGraphPath1;
            bool containsSecondaryEdges1 = false;
            for(auto it=begin1; it!=end1; ++it) {
                const Superbubble::edge_descriptor se = *it;
                const SuperbubbleEdge& sEdge = superbubble[se];
                const AssemblyGraph2::edge_descriptor ae = sEdge.ae;
                const AssemblyGraph2Edge& aEdge = g[ae];
                const AssemblyGraph2Edge::Branch& branch = aEdge.branches[sEdge.branchId];
                copy(branch.path.begin(), branch.path.end(), back_inserter(markerGraphPath1));
                if(branch.containsSecondaryEdges) {
                    containsSecondaryEdges1 = true;
                }
            }
            edge_descriptor eNew;
            bool edgeWasAdded = false;
            tie(eNew, edgeWasAdded) = add_edge(av0, av1,
                E(nextId++,
                markerGraphPath0, containsSecondaryEdges0,
                markerGraphPath1, containsSecondaryEdges1),
                g);
            SHASTA_ASSERT(edgeWasAdded);
            if(storeReadInformation) {
                g[eNew].storeReadInformation(markerGraph);
            }
            if(assemble) {
                AssemblyGraph2::assemble(eNew);
            }
        }


        // If there is a common suffix, generate a new haploid edge of the AssemblyGraph2.
        if(suffixLength) {
            const auto begin = bestPaths[0].end() - suffixLength;
            const auto end = bestPaths[0].end();

            // Construct the marker graph path.
            MarkerGraphPath markerGraphPath;
            bool containsSecondaryEdges = false;
            for(auto it=begin; it!=end; ++it) {
                const Superbubble::edge_descriptor se = *it;
                const SuperbubbleEdge& sEdge = superbubble[se];
                const AssemblyGraph2::edge_descriptor ae = sEdge.ae;
                const AssemblyGraph2Edge& aEdge = g[ae];
                const AssemblyGraph2Edge::Branch& branch = aEdge.branches[sEdge.branchId];
                copy(branch.path.begin(), branch.path.end(), back_inserter(markerGraphPath));
                if(branch.containsSecondaryEdges) {
                    containsSecondaryEdges = true;
                }
            }

            // Create a new haploid edge with this path.
            addEdge(markerGraphPath, containsSecondaryEdges, storeReadInformation, assemble);
        }



        // Now we can remove all the chunk edges from the assembly graph.
        for(const Superbubble::edge_descriptor se: superbubble.chunkEdges[chunkId]) {
            const SuperbubbleEdge& sEdge = superbubble[se];
            if(sEdge.branchId == 0) {
                boost::remove_edge(sEdge.ae, g);
            }
        }
    }

    if(debug) {
        superbubble.writeGraphviz1(cout, g);
    }
}



// Find the chunk that each edge belongs to.
// This must be called after the dominator trees
// and the critical path are computed.
void AssemblyGraph2::Superbubble::findChunks()
{
    Superbubble& superbubble = *this;
    SHASTA_ASSERT(entrances.size() == 1);
    SHASTA_ASSERT(exits.size() == 1);

    chunkEdges.clear();
    const uint64_t chunkCount = superbubble[exits.front()].positionInCriticalPath;
    chunkEdges.resize(chunkCount);

    BGL_FORALL_EDGES(e, superbubble, Superbubble) {
        findChunk(e);

        const uint64_t chunk = superbubble[e].chunk;
        if(chunk != std::numeric_limits<uint64_t>::max()) {
            chunkEdges[chunk].push_back(e);
        }
    }
}



void AssemblyGraph2::Superbubble::findChunk(edge_descriptor e)
{
    Superbubble& superbubble = *this;

    vertex_descriptor v0 = source(e, superbubble);
    vertex_descriptor v1 = target(e, superbubble);

    // Walk the forward dominator tree up starting at v0
    // and until we get on the critical path.
    uint64_t chunk;
    while(true) {
        if(superbubble[v0].positionInCriticalPath != std::numeric_limits<uint64_t>::max()) {
            chunk = superbubble[v0].positionInCriticalPath;
            break;
        }
        v0 = superbubble[v0].immediateDominator0;
        if(v0 == null_vertex()) {
            return;
        }
    }

    // Do the same in the opposite direction.
    uint64_t nextChunk;
    while(true) {
        if(superbubble[v1].positionInCriticalPath != std::numeric_limits<uint64_t>::max()) {
            nextChunk = superbubble[v1].positionInCriticalPath;
            break;
        }
        v1 = superbubble[v1].immediateDominator1;
        if(v1 == null_vertex()) {
            return;
        }
    }

    if(nextChunk == chunk + 1) {
        superbubble[e].chunk = chunk;
    }
}



// This computes the critical path from the entrance to the exit
// for a superbubble with a single entrance and exit.
// This assumes that the dominator tree were already computed.
void AssemblyGraph2::Superbubble::computeCriticalPath()
{
    Superbubble& superbubble = *this;
    SHASTA_ASSERT(entrances.size() == 1);
    SHASTA_ASSERT(exits.size() == 1);
    const vertex_descriptor entrance = entrances.front();
    const vertex_descriptor exit = exits.front();

    // Compute the critical path using the forward dominator tree.
    criticalPath.clear();
    for(Superbubble::vertex_descriptor v = exit; ; ) {
        superbubble.criticalPath.push_back(v);
        if(v == entrance) {
            break;
        }
        v = superbubble[v].immediateDominator0;
    }
    reverse(superbubble.criticalPath.begin(), superbubble.criticalPath.end());

    // Also compute the critical path using the fbackward dominator tree
    // and check that we get the same result.
    vector<vertex_descriptor> criticalPathCheck;
    for(Superbubble::vertex_descriptor v = entrance; ; ) {
        criticalPathCheck.push_back(v);
        if(v == exit) {
            break;
        }
        v = superbubble[v].immediateDominator1;
    }
    SHASTA_ASSERT(criticalPathCheck == criticalPath);

    // Store positions in the critical path in the critical path vertices.
    for(uint64_t i=0; i<criticalPath.size(); i++) {
        const auto v = criticalPath[i];
        superbubble[v].positionInCriticalPath = i;
    }
}



AssemblyGraph2::Superbubble::Superbubble(
    const AssemblyGraph2& g,
    const vector<AssemblyGraph2::vertex_descriptor>& aVertices,
    uint64_t edgeLengthThreshold)
{
    Superbubble& superbubble = *this;

    // Create the vertices.
    std::map<AssemblyGraph2::vertex_descriptor, Superbubble::vertex_descriptor> vertexMap;
    for(const AssemblyGraph2::vertex_descriptor av: aVertices) {
        Superbubble::vertex_descriptor sv = boost::add_vertex(SuperbubbleVertex(av), superbubble);
        vertexMap.insert(make_pair(av, sv));
    }

    /*
    cout << "Vertex map:" << endl;
    for(const auto& p: vertexMap) {
        cout << p.second << " (" << p.first << ")" << endl;
    }
    */

    // Create the edges.
    // For an edge to be part of the superbubble:
    // - Both its source and target vertices must be part of the superbubble.
    // - It must be a short edge (maximumPathLength() <= edgeLengthThreshold).
    BGL_FORALL_VERTICES(sv0, superbubble, Superbubble) {
        const AssemblyGraph2::vertex_descriptor av0 = superbubble[sv0].av;
        BGL_FORALL_OUTEDGES(av0, ae, g, G) {
            const AssemblyGraph2::vertex_descriptor av1 = target(ae, g);
            auto it = vertexMap.find(av1);
            if(it != vertexMap.end()) {
                const E& aEdge = g[ae];
                if(aEdge.maximumPathLength() <= edgeLengthThreshold) {
                    const Superbubble::vertex_descriptor sv1 = it->second;
                    for(uint64_t branchId=0; branchId<aEdge.ploidy(); branchId++) {
                        add_edge(sv0, sv1, SuperbubbleEdge(ae, branchId), superbubble);
                    }
                }
            }
        }
    }



    // Find the entrances and exits.
    // An entrance has:
    // - At least one in-edge from a vertex outside the superbubble.
    // - At least one out-edge to a vertex inside the superbubble.
    // An exit has:
    // - At least one in-edge from a vertex inside the superbubble.
    // - At least one out-edge to a vertex outside the superbubble.
    BGL_FORALL_VERTICES(sv0, superbubble, Superbubble) {
        const AssemblyGraph2::vertex_descriptor av0 = superbubble[sv0].av;

        // Check in-edges.
        bool hasInedgesFromOutside = false;
        bool hasInedgesFromInside = false;
        BGL_FORALL_INEDGES(av0, ae, g, G) {
            const AssemblyGraph2::vertex_descriptor av1 = source(ae, g);
            if(av1 == av0) {
                continue;
            }
            if(vertexMap.find(av1) == vertexMap.end()) {
                hasInedgesFromOutside = true;
            } else {
                hasInedgesFromInside = true;
            }
        }


        // Check out-edges.
        bool hasOutedgesToOutside = false;
        bool hasOutedgesToInside = false;
        BGL_FORALL_OUTEDGES(av0, ae, g, G) {
            const AssemblyGraph2::vertex_descriptor av1 = target(ae, g);
            if(av1 == av0) {
                continue;
            }
            if(vertexMap.find(av1) == vertexMap.end()) {
                hasOutedgesToOutside = true;
            } else {
                hasOutedgesToInside = true;
            }
        }

        if(hasInedgesFromOutside and hasOutedgesToInside) {
            entrances.push_back(sv0);
        }
        if(hasInedgesFromInside and hasOutedgesToOutside) {
            exits.push_back(sv0);
        }
    }
}



AssemblyGraph2::Superbubble::Superbubble(
    const AssemblyGraph2& g,
    const vector<AssemblyGraph2::vertex_descriptor>& aVertices)
{
    Superbubble& superbubble = *this;

    // Create the vertices.
    std::map<AssemblyGraph2::vertex_descriptor, Superbubble::vertex_descriptor> vertexMap;
    for(const AssemblyGraph2::vertex_descriptor av: aVertices) {
        Superbubble::vertex_descriptor sv = boost::add_vertex(SuperbubbleVertex(av), superbubble);
        vertexMap.insert(make_pair(av, sv));
    }

    /*
    cout << "Vertex map:" << endl;
    for(const auto& p: vertexMap) {
        cout << p.second << " (" << p.first << ")" << endl;
    }
    */

    // Create the edges.
    // For an edge to be part of the superbubble:
    // - Both its source and target vertices must be part of the superbubble.
    BGL_FORALL_VERTICES(sv0, superbubble, Superbubble) {
        const AssemblyGraph2::vertex_descriptor av0 = superbubble[sv0].av;
        BGL_FORALL_OUTEDGES(av0, ae, g, G) {
            const AssemblyGraph2::vertex_descriptor av1 = target(ae, g);
            auto it = vertexMap.find(av1);
            if(it != vertexMap.end()) {
                const E& aEdge = g[ae];
                SHASTA_ASSERT(aEdge.bubbleChain.first == 0);
                const Superbubble::vertex_descriptor sv1 = it->second;
                for(uint64_t branchId=0; branchId<aEdge.ploidy(); branchId++) {
                    add_edge(sv0, sv1, SuperbubbleEdge(ae, branchId), superbubble);
                }
            }
        }
    }



    // Find the entrances and exits.
    // An entrance has:
    // - At least one in-edge from a vertex outside the superbubble.
    // - At least one out-edge to a vertex inside the superbubble.
    // An exit has:
    // - At least one in-edge from a vertex inside the superbubble.
    // - At least one out-edge to a vertex outside the superbubble.
    BGL_FORALL_VERTICES(sv0, superbubble, Superbubble) {
        const AssemblyGraph2::vertex_descriptor av0 = superbubble[sv0].av;

        // Check in-edges.
        bool hasInedgesFromOutside = false;
        bool hasInedgesFromInside = false;
        BGL_FORALL_INEDGES(av0, ae, g, G) {
            const AssemblyGraph2::vertex_descriptor av1 = source(ae, g);
            if(av1 == av0) {
                continue;
            }
            if(vertexMap.find(av1) == vertexMap.end()) {
                hasInedgesFromOutside = true;
            } else {
                hasInedgesFromInside = true;
            }
        }


        // Check out-edges.
        bool hasOutedgesToOutside = false;
        bool hasOutedgesToInside = false;
        BGL_FORALL_OUTEDGES(av0, ae, g, G) {
            const AssemblyGraph2::vertex_descriptor av1 = target(ae, g);
            if(av1 == av0) {
                continue;
            }
            if(vertexMap.find(av1) == vertexMap.end()) {
                hasOutedgesToOutside = true;
            } else {
                hasOutedgesToInside = true;
            }
        }

        if(hasInedgesFromOutside and hasOutedgesToInside) {
            entrances.push_back(sv0);
        }
        if(hasInedgesFromInside and hasOutedgesToOutside) {
            exits.push_back(sv0);
        }
    }
}



void AssemblyGraph2::Superbubble::writeGraphviz(
    ostream& out,
    const AssemblyGraph2& g) const
{
    const Superbubble& superbubble = *this;

    out << "digraph Superbubble {" << endl;
    BGL_FORALL_VERTICES(sv, superbubble, Superbubble) {
        const AssemblyGraph2::vertex_descriptor av = superbubble[sv].av;
        const AssemblyGraph2Vertex& aVertex = g[av];
        out << aVertex.markerGraphVertexId << ";\n";
    }

    BGL_FORALL_EDGES(se, superbubble, Superbubble) {
        const SuperbubbleEdge sEdge = superbubble[se];
        const AssemblyGraph2::edge_descriptor ae = sEdge.ae;
        const uint64_t branchId = sEdge.branchId;

        const Superbubble::vertex_descriptor sv0 = source(se, superbubble);
        const Superbubble::vertex_descriptor sv1 = target(se, superbubble);

        const AssemblyGraph2::vertex_descriptor av0 = superbubble[sv0].av;
        const AssemblyGraph2::vertex_descriptor av1 = superbubble[sv1].av;

        const AssemblyGraph2Vertex& aVertex0 = g[av0];
        const AssemblyGraph2Vertex& aVertex1 = g[av1];

        out <<
            aVertex0.markerGraphVertexId << "->" <<
            aVertex1.markerGraphVertexId <<
            " [label=\"" << g[ae].pathId(branchId) << "\"];\n";
    }

    out << "}" << endl;
}



// This assumes that there is exactly one entrance and one exit,
// and that the data structures created by handleSuperbubble1
// are available.
void AssemblyGraph2::Superbubble::writeGraphviz1(
    ostream& out,
    const AssemblyGraph2& g) const
{
    const Superbubble& superbubble = *this;
    SHASTA_ASSERT(entrances.size() == 1);
    SHASTA_ASSERT(exits.size() == 1);



    out << "digraph Superbubble {" << endl;



    BGL_FORALL_VERTICES(sv, superbubble, Superbubble) {
        const AssemblyGraph2::vertex_descriptor av = superbubble[sv].av;
        const AssemblyGraph2Vertex& aVertex = g[av];
        const uint64_t positionInCriticalPath = superbubble[sv].positionInCriticalPath;

        out << aVertex.markerGraphVertexId;
        out << "[";

        // Label.
        out << "label=\"" << aVertex.markerGraphVertexId;
        if(positionInCriticalPath != std::numeric_limits<uint64_t>::max()) {
            out << "\\n" << positionInCriticalPath;
        }
        out << "\"";


        // Color.
        if(sv == entrances.front()) {
            out << " style=filled fillcolor=green";
        } else if(sv == exits.front()) {
            out << " style=filled fillcolor=red";
        } else {
            if(positionInCriticalPath != std::numeric_limits<uint64_t>::max()) {
                out << " style=filled fillcolor=cyan";
            }
        }

        out << "];\n";
    }



    BGL_FORALL_EDGES(se, superbubble, Superbubble) {
        const SuperbubbleEdge sEdge = superbubble[se];
        const AssemblyGraph2::edge_descriptor ae = sEdge.ae;
        const uint64_t branchId = sEdge.branchId;

        const Superbubble::vertex_descriptor sv0 = source(se, superbubble);
        const Superbubble::vertex_descriptor sv1 = target(se, superbubble);

        const AssemblyGraph2::vertex_descriptor av0 = superbubble[sv0].av;
        const AssemblyGraph2::vertex_descriptor av1 = superbubble[sv1].av;

        const AssemblyGraph2Vertex& aVertex0 = g[av0];
        const AssemblyGraph2Vertex& aVertex1 = g[av1];

        out <<
            aVertex0.markerGraphVertexId << "->" <<
            aVertex1.markerGraphVertexId <<
            " [label=\"" << g[ae].pathId(branchId);
        if(sEdge.chunk != std::numeric_limits<uint64_t>::max()) {
            out << "\\n" << sEdge.chunk;
        }
        out << "\"];\n";
    }

    out << "}" << endl;
}



// Return true if the superbubble corresponds to a simple linear chain
// in the AssemblyGraph2.
bool AssemblyGraph2::Superbubble::isSimpleLinearChain() const
{
    const Superbubble& superbubble = *this;

    // A simple linear chain must have exactly one entrance and one exit.
    if(entrances.size() != 1) {
        return false;
    }
    if(exits.size() != 1) {
        return false;
    }

    // The entrance must have in_degree 0 and out_degree 1.
    const vertex_descriptor sEntrance = entrances.front();
    if(originalInDegree(sEntrance) != 0) {
        return false;
    }
    if(originalOutDegree(sEntrance) != 1) {
        return false;
    }

    // The exit must have in_degree 1 and out_degree 0.
    const vertex_descriptor sExit = exits.front();
    if(originalInDegree(sExit) != 1) {
        return false;
    }
    if(originalOutDegree(sExit) != 0) {
        return false;
    }

    // All other vertices must have in_degree and out_degree 1.
    BGL_FORALL_VERTICES(sv, superbubble, Superbubble) {
        if(sv == sEntrance) {
            continue;
        }
        if(sv == sExit) {
            continue;
        }
        if(originalInDegree(sv) != 1) {
            return false;
        }
        if(originalOutDegree(sv) != 1) {
            return false;
        }
    }

    // If getting here, all conditions for a simple linear chains are satisfied.
    return true;
}


// Return the number of distinct AssemblyGraph2 edges
// that begin/end at a given vertex.
uint64_t AssemblyGraph2::Superbubble::originalInDegree(vertex_descriptor v) const
{
    const Superbubble& superbubble = *this;
    uint64_t n = 0;
    BGL_FORALL_INEDGES(v, e, superbubble, Superbubble) {
        if(superbubble[e].branchId == 0) {
            ++n;
        }
    }
    return n;
}
uint64_t AssemblyGraph2::Superbubble::originalOutDegree(vertex_descriptor v) const
{
    const Superbubble& superbubble = *this;
    uint64_t n = 0;
    BGL_FORALL_OUTEDGES(v, e, superbubble, Superbubble) {
        if(superbubble[e].branchId == 0) {
            ++n;
        }
    }
    return n;
}



// Enumerate paths from the entrance to the exit.
// There must be exactly one entrance and one exit.
void AssemblyGraph2::Superbubble::enumeratePaths()
{
    SHASTA_ASSERT(entrances.size() == 1);
    SHASTA_ASSERT(exits.size() == 1);

    const vertex_descriptor entrance = entrances.front();
    const vertex_descriptor exit = exits.front();
    enumeratePaths(entrance, exit);
}



// Enumerate paths from an entrance to an exit.
void AssemblyGraph2::Superbubble::enumeratePaths(
    vertex_descriptor entrance,
    vertex_descriptor exit)
{
    enumerateSelfAvoidingPaths(*this, entrance, exit, paths);
}



// Iteratively remove bad bubbles using the PhasingGraph.
void AssemblyGraph2::removeBadBubblesIterative(
    uint64_t minConcordantReadCount,
    uint64_t maxDiscordantReadCount,
    double minLogP,
    double epsilon,
    uint64_t superbubbleRemovalEdgeLengthThreshold,
    uint64_t maxSuperbubbleSize,
    uint64_t maxSuperbubbleChunkSize,
    uint64_t maxSuperbubbleChunkPathCount,
    uint64_t pruneLength,
    uint64_t componentSizeThreshold,
    size_t threadCount)
{
    performanceLog << timestamp << "AssemblyGraph2::removeBadBubblesIterative begins." << endl;

    G& g = *this;
    const bool debug = false;

    for(uint64_t iteration=0; ; iteration++) {
        performanceLog << timestamp << "Removing bad bubbles: iteration " << iteration << " begins." << endl;

        // Assign each diploid bubble to its own component.
        // This way the PhasingGraph will have one bubble per vertex.
        uint64_t componentId = 0;
        BGL_FORALL_EDGES(e, g, G) {
            AssemblyGraph2Edge& edge = g[e];
            if(edge.ploidy() == 2) {
                edge.componentId = componentId++;
                edge.phase = 0;
            }
        }


        // Create the PhasingGraph.
        const bool allowRandomHypothesis = true;
        PhasingGraph phasingGraph(
            g,
            minConcordantReadCount, maxDiscordantReadCount, minLogP, epsilon,
            threadCount, allowRandomHypothesis);
        cout << "The number of diploid bubbles before iteration " <<
            iteration << " is " << num_vertices(phasingGraph) << endl;

        // Compute the optimal spanning tree.
        phasingGraph.computeSpanningTree();

        // Use the optimal spanning tree to phase the PhasingGraph,
        // but don't store the result in the AssemblyGraph2.
        phasingGraph.phase();
        if(debug) {
            const string s = "A-" + to_string(iteration);
            phasingGraph.writeCsv("PhasingGraph-" + s, g);
            phasingGraph.writeGraphviz("PhasingGraph-" + s + ".dot");
         }

        // Gather the bubbles in each connected component.
        vector< vector<PhasingGraph::vertex_descriptor> > components;
        BGL_FORALL_VERTICES(v, phasingGraph, PhasingGraph) {
            const uint64_t componentId = phasingGraph[v].componentId;

            // Make sure the components vector is long enough.
            if(componentId >= components.size()) {
                components.resize(componentId + 1);
            }

            // Store it.
            components[componentId].push_back(v);
        }

        // Write out a histogram of component sizes.
        if(debug) {
            std::map<uint64_t, uint64_t> histogram;
            for(const vector<PhasingGraph::vertex_descriptor>& component: components) {
                const uint64_t size = component.size();
                auto it = histogram.find(size);
                if(it == histogram.end()) {
                    tie(it, ignore) = histogram.insert(make_pair(size, 0));
                }
                ++(it->second);
            }

            /*
            cout << "Histogram of component sizes:\n";
            cout << "ComponentSize,Frequency,Total\n";
            for(const auto& p: histogram) {
                if(p.first > 0) {
                    cout << p.first << ",";
                    cout << p.second << ",";
                    cout << p.first * p.second<< "\n";
                }
            }
            */
        }


        // Gather the bubbles that are going to be removed.
        // Right now, these are the bubbles in small connected components.
        // Later, we will add bubbles with inconsistent edges.
        vector<PhasingGraph::vertex_descriptor> badBubbles;
        for(const vector<PhasingGraph::vertex_descriptor>& component: components) {
            if(component.size() >= componentSizeThreshold) {
                continue;
            }
            for(const PhasingGraph::vertex_descriptor v: component) {
                badBubbles.push_back(v);
            }
        }

        // If no bubbles were marked as bad, end the iteration.
        if(badBubbles.empty()) {
            break;
        }

        // Remove the bubbles we marked as bad.
        for(const PhasingGraph::vertex_descriptor v: badBubbles) {
            const PhasingGraphVertex& vertex = phasingGraph[v];
            SHASTA_ASSERT(vertex.bubbles.size() == 1);
            const AssemblyGraph2::edge_descriptor e = vertex.bubbles.front().first;
            g[e].removeAllBranchesExceptStrongest();
        }

        /*
        cout << "Removed " << badBubbles.size() << " bubbles." << endl;
        cout << num_vertices(phasingGraph) - badBubbles.size() <<
            " diploid bubbles are left." << endl;
        */


        // Some new bubbles may form after we merge.
        merge(true, true);
        gatherBubbles();
        forceMaximumPloidy(2);

        // Handle superbubbles that may have appeared as a result of removing bubbles.
        handleSuperbubbles0(superbubbleRemovalEdgeLengthThreshold,
            maxSuperbubbleSize, maxSuperbubbleChunkSize, maxSuperbubbleChunkPathCount, true, true);
        merge(true, true);
        handleSuperbubbles1(
            maxSuperbubbleSize, maxSuperbubbleChunkSize, maxSuperbubbleChunkPathCount, true, true);
        merge(true, true);
        prune(pruneLength);

        performanceLog << timestamp << "Removing bad bubbles: iteration " << iteration << " ends." << endl;
    }

    performanceLog << timestamp << "AssemblyGraph2::removeBadBubblesIterative ends." << endl;
}



void AssemblyGraph2::hierarchicalPhase(
    uint64_t minConcordantReadCount,
    uint64_t maxDiscordantReadCount,
    double minLogP,
    double epsilon,
    size_t threadCount)
{
    performanceLog << timestamp << "AssemblyGraph2::hierarchicalPhase begins." << endl;

    G& g = *this;
    const bool debug = false;

    // Start by assigning each diploid bubble to its own component.
    uint64_t componentId = 0;
    BGL_FORALL_EDGES(e, g, G) {
        AssemblyGraph2Edge& edge = g[e];
        if(edge.ploidy() == 2) {
            edge.componentId = componentId++;
            edge.phase = 0;
        }
    }

    // Main iteration loop.
    // At each iteration some phasing components can be combined into larger
    // phasing components.
    for(uint64_t iteration=0; ; iteration++) {
        performanceLog << timestamp << "Hierarchical phasing iteration " << iteration << " begins." << endl;


        // Create the PhasingGraph.
        const bool allowRandomHypothesis = false;
        PhasingGraph phasingGraph(
            g, minConcordantReadCount, maxDiscordantReadCount, minLogP, epsilon,
            threadCount, allowRandomHypothesis);
        cout << "The phasing graph has " << num_vertices(phasingGraph) <<
            " vertices and " << num_edges(phasingGraph) << " edges." << endl;

        // Compute the optimal spanning tree.
        phasingGraph.computeSpanningTree();

        // If the PhasingGraph has no edges, exit the inner iteration loop.
        // Phasing information in the AssemblyGraph2 is as stored in the previous inner iteration.
        if(num_edges(phasingGraph)== 0) {
            break;
        }

        // Use the optimal spanning tree to phase the PhasingGraph,
        // then store the result in the AssemblyGraph2.
        phasingGraph.phase();
        if(debug) {
            const string s = "B-" + to_string(iteration);
            phasingGraph.writeCsv("PhasingGraph-" + s, g);
            phasingGraph.writeGraphviz("PhasingGraph-" + s + ".dot");
        }

        phasingGraph.storePhasing(g);
        performanceLog << timestamp << "Hierarchical phasing iteration " << iteration << " ends." << endl;
    }

    // Create a final PhasinGraph with permissive criteria for edge creation.
    PhasingGraph phasingGraph(g, 0, 100, 0., epsilon, threadCount, false);
    if(true) {
        phasingGraph.writeCsv("PhasingGraph-Final", g);
        phasingGraph.writeGraphviz("PhasingGraph-Final.dot");
    }

    performanceLog << timestamp << "AssemblyGraph2::hierarchicalPhase ends." << endl;
}



// Renumber component to make them contiguous starting at 0.
void AssemblyGraph2::renumberComponents()
{
    G& g = *this;

    vector<uint64_t> componentIds;

    BGL_FORALL_EDGES(e, g, G) {
        const AssemblyGraph2Edge& edge = g[e];

        if(edge.ploidy() != 2) {
            continue;
        }

        const uint64_t componentId = edge.componentId;
        if(componentId == AssemblyGraph2Edge::invalidComponentId) {
            continue;
        }

        componentIds.push_back(componentId);
    }
    deduplicate(componentIds);



    // Replace component ids with the corresponding index in the componentIds vector.
    BGL_FORALL_EDGES(e, g, G) {
        AssemblyGraph2Edge& edge = g[e];

        if(edge.ploidy() != 2) {
            continue;
        }

        const uint64_t componentId = edge.componentId;
        if(componentId == AssemblyGraph2Edge::invalidComponentId) {
            continue;
        }

        auto it = std::lower_bound(componentIds.begin(), componentIds.end(), componentId);
        SHASTA_ASSERT(it != componentIds.end());
        SHASTA_ASSERT(*it == componentId);
        edge.componentId = it - componentIds.begin();
    }
}



// Remove short loop-back edges.
void AssemblyGraph2::removeShortLoopbackEdges(uint64_t edgeLengthThreshold)
{
    G& g = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, g, G) {
        const AssemblyGraph2Edge& edge = g[e];

        if(edge.ploidy() > 1) {
            continue;
        }

        if(edge.branches.front().path.size() >= edgeLengthThreshold) {
            continue;
        }

        if(source(e, g) != target(e, g)) {
            continue;
        }

        edgesToBeRemoved.push_back(e);
    }

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, g);
    }
}
