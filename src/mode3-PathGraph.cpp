// Shasta.
#include "mode3-PathGraph.hpp"
#include "findLinearChains.hpp"
#include "orderPairs.hpp"
#include "transitiveReduction.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>

// Standard library.
#include "iostream.hpp"
#include <queue>
#include <stack>



// Create the PathGraph from the AssemblyGraph.
// Start with a single segment for each vertex
// (that is, paths of length 1).
PathGraph::PathGraph(const AssemblyGraph& assemblyGraph) :
    MultithreadedObject<PathGraph>(*this),
    assemblyGraph(assemblyGraph)
{
    // HARDWIRED CONSTANTS TO BE EXPOSED WHEN CODE STABILIZES.
    const uint64_t minCoverage = 3;
    const uint64_t partitionMaxDistance = 10;
    const uint64_t minSubgraphSize = 8;

    PathGraph& pathGraph = *this;

    createVertices();
    createEdges(minCoverage);
    cout << "The initial path graph has " << num_vertices(pathGraph) <<
        " vertices and " << num_edges(pathGraph) << " edges." << endl;

    computeJourneys();
    writeJourneys("PathGraphJourneys.csv");

    // Partition the PathGraph into subgraphs.
    partition(partitionMaxDistance, minSubgraphSize);
    writeGfa("PathGraph");
    writeCsvDetailed("PathGraphDetailed.csv");

    uint64_t subgraphId;
    cout << "Enter a subgraph to detangle:" << endl;
    cin >> subgraphId;
    vector<PathGraphVertex> newVertices;
    detangleSubgraph(subgraphId, newVertices, true);
}



// Initial creation of the vertices.
// Start with a single segment for each vertex
// (that is, paths of length 1).
void PathGraph::createVertices() {

    PathGraph& pathGraph = *this;


    // Create a vertex for each segment in the AssemblyGraph.
    for(uint64_t segmentId=0; segmentId<assemblyGraph.paths.size(); segmentId++) {

        // Create the vertex.
        const vertex_descriptor v = add_vertex(pathGraph);
        PathGraphVertex& vertex = pathGraph[v];
        vertex.id = nextVertexId++;

        // Store the path.
        vertex.path.push_back(segmentId);

        // Store the AssemblyGraphJourneyInterval's.
        const span<const pair<OrientedReadId, uint64_t> > journeyInfos =
            assemblyGraph.assemblyGraphJourneyInfos[segmentId];
        for(const pair<OrientedReadId, uint64_t>& p: journeyInfos) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t position = p.second;
            AssemblyGraphJourneyInterval interval;
            interval.orientedReadId = orientedReadId;
            interval.first = position;
            interval.last = position;
            vertex.journeyIntervals.push_back(
                make_pair(interval, std::numeric_limits<uint64_t>::max()));
        }
    }

}



// Recreate all edges from scratch, using only the
// information stored in the vertices.
void PathGraph::createEdges(uint64_t minCoverage)
{
    PathGraph& pathGraph = *this;

    // Gather AssemblyGraphJourneyInterval's for all oriented reads.
    vector< vector<pair<AssemblyGraphJourneyInterval, vertex_descriptor> > >
        journeyIntervals(2 * assemblyGraph.readCount());
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        for(const auto& p: pathGraph[v].journeyIntervals) {
            const AssemblyGraphJourneyInterval& interval = p.first;
            journeyIntervals[interval.orientedReadId.getValue()].push_back(
                make_pair(interval, v));
        }
    }
    for(auto& v: journeyIntervals) {
        sort(v.begin(), v.end(),
            OrderPairsByFirstOnly<AssemblyGraphJourneyInterval, vertex_descriptor>());
    }

    // Create the edges.
    for(const auto& orientedReadJourneyIntervals: journeyIntervals) {

        for(uint64_t i=1; i<orientedReadJourneyIntervals.size(); i++) {
            const vertex_descriptor v0 = orientedReadJourneyIntervals[i-1].second;
            const vertex_descriptor v1 = orientedReadJourneyIntervals[i  ].second;

            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = edge(v0, v1, pathGraph);
            if(not edgeExists) {
                tie(e, edgeExists) = add_edge(v0, v1, pathGraph);
                SHASTA_ASSERT(edgeExists);
            }
            ++pathGraph[e].coverage;
        }
    }



    // Remove the low coverage edges.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        if(pathGraph[e].coverage < minCoverage) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, pathGraph);
    }
}



// Compute the journeys of all oriented reads in the PathGraph.
// The journey of an oriented read in the PathGraph is
// a sequence of vertex descriptors which is not necessarily a path.
// Indexed by OrientedReadId::getValue();
void PathGraph::computeJourneys()
{
    PathGraph& pathGraph = *this;
    const ReadId readCount = ReadId(assemblyGraph.readCount());

    // First create, for each oriented read, a vector
    // of pairs (AssemblyGraphJourneyInterval, vertex_descriptor).
    vector< vector< pair<AssemblyGraphJourneyInterval, vertex_descriptor> > >
        journeyTable(2 * readCount);
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        for(const auto& p: pathGraph[v].journeyIntervals) {
            const AssemblyGraphJourneyInterval& journeyInterval = p.first;
            journeyTable[journeyInterval.orientedReadId.getValue()].push_back(make_pair(journeyInterval, v));
        }
    }

    // Sort them and sanity check.
    for(vector< pair<AssemblyGraphJourneyInterval, vertex_descriptor> >& v: journeyTable) {
        sort(v.begin(), v.end());

        // Sanity check.
        if(v.size() > 1) {
            for(uint64_t i=1; i<v.size(); i++) {
                const AssemblyGraphJourneyInterval& previous = v[i-1].first;
                const AssemblyGraphJourneyInterval& current = v[i].first;
                SHASTA_ASSERT(previous.last < current.first);
            }
        }
    }


    // Store what we got.
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        pathGraph[v].journeyIntervals.clear();
    }
    journeys.clear();
    journeys.resize(2 * readCount);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const uint64_t index = orientedReadId.getValue();
            for(uint64_t position=0; position<journeyTable[index].size(); position++) {
                const auto& p = journeyTable[index][position];
                const AssemblyGraphJourneyInterval& interval = p.first;
                const vertex_descriptor v = p.second;
                journeys[index].push_back(v);
                pathGraph[v].journeyIntervals.push_back(make_pair(interval, position));
            }
        }
    }
}



void PathGraph::writeJourneys(const string& fileName) const
{
    const PathGraph& pathGraph = *this;
    ofstream csv(fileName);

    // Loop over all oriented reads.
    const ReadId readCount = ReadId(assemblyGraph.readCount());
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            csv << orientedReadId << ",";

            // Write the journey of this oriented read in the PathGraph.
            const auto journey = journeys[orientedReadId.getValue()];
            for(const vertex_descriptor v: journey) {
                csv << pathGraph[v].id << ",";
            }
            csv << "\n";
        }
    }
}



// Partition the PathGraph into subgraphs.
void PathGraph::partition(
    uint64_t maxDistance,
    uint64_t minSubgraphSize)
{
    PathGraph& pathGraph = *this;

    // Mark all vertices as not assigned to any partition.
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        pathGraph[v].subgraphId = noSubgraph;
    }

    // Start at all vertices with zero in-degree,
    // plus the boundary vertices we find that way.
    vector<vertex_descriptor> boundaryVertices;
    std::stack<vertex_descriptor> s;
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        if(in_degree(v, pathGraph) == 0) {
            s.push(v);
        }
    }
    uint64_t subgraphId = 0;
    while(not s.empty()) {
        const vertex_descriptor v = s.top();
        s.pop();

        if(pathGraph[v].subgraphId == noSubgraph) {
            partitionIteration(v, maxDistance, subgraphId++, boundaryVertices);
            for(const vertex_descriptor v: boundaryVertices) {
                s.push(v);
            }
        }
    }



    // In exceptional cases, the above procedure might not assign all
    // vertices to a subgraph.
    // This code takes care of that.
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        if(pathGraph[v].subgraphId == noSubgraph) {
            partitionIteration(v, maxDistance, subgraphId++, boundaryVertices);
        }
    }



    // Combine small subgraphs with adjacent subgraphs, if possible.
    // This can leave subgraphs with size 0, but we don't worry about that.
    while(true) {

        // Gather the subgraphs based on the current settings of
        // the vertices subgraphId.
        gatherSubgraphs();

        // Find the small subgraphs.
        std::set<uint64_t> smallSubgraphs;
        for(uint64_t subgraphId=0; subgraphId<subgraphs.size(); subgraphId++) {
            const vector<vertex_descriptor>& subgraph = subgraphs[subgraphId];
            const uint64_t subgraphSize = subgraph.size();
            if((subgraphSize != 0) and (subgraph.size() < minSubgraphSize)) {
                smallSubgraphs.insert(subgraphId);
            }
        }



        // Try and merge small subgraphs with adjacent subgraphs.

        // Loop over small subgraphs.
        bool changesWereMade = false;
        for(uint64_t subgraphId0: smallSubgraphs) {
            const vector<vertex_descriptor>& subgraph0 = subgraphs[subgraphId0];
            const uint64_t subgraph0Size = subgraph0.size();
            SHASTA_ASSERT(subgraph0Size < minSubgraphSize);

            // Find adjacent subgraphs and their sizes.
            vector< pair<uint64_t, uint64_t> > adjacentSubgraphsTable; // (size, subgraphId) of adjacent.
            for(const vertex_descriptor v0: subgraph0) {
                BGL_FORALL_OUTEDGES(v0, e, pathGraph, PathGraph) {
                    const vertex_descriptor v1 = target(e, pathGraph);
                    const uint64_t subgraphId1 = pathGraph[v1].subgraphId;
                    if(subgraphId1 != subgraphId0){
                        adjacentSubgraphsTable.push_back(make_pair(subgraphs[subgraphId1].size(), subgraphId1));
                    }
                }
                BGL_FORALL_INEDGES(v0, e, pathGraph, PathGraph) {
                    const vertex_descriptor v1 = source(e, pathGraph);
                    const uint64_t subgraphId1 = pathGraph[v1].subgraphId;
                    if(subgraphId1 != subgraphId0){
                        adjacentSubgraphsTable.push_back(make_pair(subgraphs[subgraphId1].size(), subgraphId1));
                    }
                }
            }
            sort(adjacentSubgraphsTable.begin(), adjacentSubgraphsTable.end());

            // Merge it with the smallest adjacent subgraph.
            const uint64_t subgraphId1 = adjacentSubgraphsTable.front().second;
            smallSubgraphs.erase(subgraphId1);
            for(const vertex_descriptor v0: subgraph0) {
                pathGraph[v0].subgraphId = subgraphId1;
            }
            changesWereMade = true;
        }

        if(not changesWereMade) {
            break;
        }
    }


    // Sort the vertex descriptors in each subgraph.
    for(vector<vertex_descriptor>& subgraph: subgraphs) {
        sort(subgraph.begin(), subgraph.end(), PathGraphOrderVerticesById(pathGraph));
    }



    // Subgraph statistics.
    cout << "Partitioned the path graph into " << subgraphs.size() << " subgraphs." << endl;
    histogramSubgraphs();

    // Count the edges across subgraphs.
    uint64_t crossEdgeCount = 0;
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        const vertex_descriptor v0 = source(e, pathGraph);
        const vertex_descriptor v1 = target(e, pathGraph);
        if(pathGraph[v0].subgraphId != pathGraph[v1].subgraphId) {
            ++crossEdgeCount;
        }
    }
    cout << "Number of edges across subgraphs is " << crossEdgeCount << endl;
}



// A partition iteration does a single BFS starting at v.
// It moves forward from v, avoiding vertices already
// assigned to a subgraph, and up to maxDistance from v.
// It also returns the boundaryVertices, that is the
// vertices found in the process that are at distance maxDistance+1
// from v and are not yet assigned to a subgraph.
// These can then used as starting points new partition iterations.
void PathGraph::partitionIteration(
    vertex_descriptor v,
    uint64_t maxDistance,
    uint64_t subgraphId,
    vector<vertex_descriptor>& boundaryVertices)
{
    PathGraph& pathGraph = *this;

    boundaryVertices.clear();

    // Initialize the BFS.
    std::queue<vertex_descriptor> q;
    q.push(v);
    PathGraphVertex& vertex = pathGraph[v];
    SHASTA_ASSERT(vertex.subgraphId == noSubgraph);
    vertex.subgraphId = subgraphId;
    vertex.distance = 0;

    // BFS loop.
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        const uint64_t distance0 = pathGraph[v0].distance;
        const uint64_t distance1 = distance0 + 1;
        SHASTA_ASSERT(distance0 <= maxDistance);

        // Loop over edges starting at v0.
        BGL_FORALL_OUTEDGES(v0, e01, pathGraph, PathGraph) {
            const vertex_descriptor v1 = target(e01, pathGraph);
            PathGraphVertex& vertex1 = pathGraph[v1];

            // If v1 is already in a subgraph, skip it.
            if(vertex1.subgraphId != noSubgraph) {
                continue;
            }

            // Assign v1 to this subgraph, if it is within maxDistance.
            if(distance1 <= maxDistance) {
                vertex1.subgraphId = subgraphId;
                vertex1.distance = distance1;
            }

            // Queue it or add it to the boundary vertices.
            if(distance1 <= maxDistance) {
                q.push(v1);
            } else {
                SHASTA_ASSERT(distance1 == maxDistance + 1);
                boundaryVertices.push_back(v1);
            }

        }

    }
}



// Gather subgraphs using the subgraphId stored in each vertex.
void PathGraph::gatherSubgraphs()
{
    PathGraph& pathGraph = *this;

    subgraphs.clear();
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        const uint64_t subgraphId = pathGraph[v].subgraphId;
        SHASTA_ASSERT(subgraphId != noSubgraph);

        if(subgraphId >= subgraphs.size()) {
            subgraphs.resize(subgraphId + 1);
        }

        subgraphs[subgraphId].push_back(v);
    }
}



void PathGraph::histogramSubgraphs()
{
    vector<uint64_t> histogram;
    for(const vector<vertex_descriptor>& subgraph: subgraphs) {
        const uint64_t subgraphSize = subgraph.size();
        if(subgraphSize >= histogram.size()) {
            histogram.resize(subgraphSize + 1, 0);
        }
        ++histogram[subgraphSize];
    }

    ofstream csv("PathGraphSubgraphHistogram.csv");
    csv << "Size,Frequency,Vertices\n";
    for(uint64_t subgraphSize=0; subgraphSize<histogram.size(); subgraphSize++) {
        const uint64_t frequency = histogram[subgraphSize];
        csv << subgraphSize << ",";
        csv << frequency << ",";
        csv << subgraphSize*frequency << "\n";
    }
}




void PathGraph::writeGfa(const string& baseName) const
{
    const PathGraph& pathGraph = *this;

    // Open the gfa and write the header.
    ofstream gfa(baseName + ".gfa");
    gfa << "H\tVN:Z:1.0\n";

    // Open the csv and write the header.
    ofstream csv(baseName + ".csv");
    csv << "PathGraph-VertexId,Color,SubgraphId\n";

    // Write each vertex as a segment in the gfa.
    // Note these segments are different from assembly graph segments:
    // here each segment represents a vertex of the path graph.
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        gfa <<
            "S\t" <<
            pathGraph[v].id << "\t" // Segment name
            "*"                     // Segment length
            "\n";


        // Color based on the subgraphId.
        const uint64_t subgraphId = pathGraph[v].subgraphId;
        string color = "LightGrey";
        if(subgraphId != noSubgraph) {
            const uint64_t r = MurmurHash2(&subgraphId,  sizeof(subgraphId),  231) &255;
            const uint64_t g = MurmurHash2(&subgraphId,  sizeof(subgraphId),  233) &255;
            const uint64_t b = MurmurHash2(&subgraphId,  sizeof(subgraphId),  235) &255;

            std::ostringstream s;
            s.fill('0');
            s << "#";
            s << hex << std::setw(2) << r;
            s << hex << std::setw(2) << g;
            s << hex << std::setw(2) << b;
            color = s.str();
        }

        csv << pathGraph[v].id << "," << color << "," << subgraphId << "\n";

    }

    // Write each edge as a link.
    BGL_FORALL_EDGES(e, pathGraph, PathGraph) {
        const vertex_descriptor v0 = source(e, pathGraph);
        const vertex_descriptor v1 = target(e, pathGraph);
        gfa <<
            "L\t" <<
            pathGraph[v0].id << "\t+\t" <<
            pathGraph[v1].id << "\t+\t0M\n";
    }

}



void PathGraph::writeCsvDetailed(const string& fileName) const
{
    const PathGraph& pathGraph = *this;
    ofstream csv(fileName);
    csv << "PathGraph-VertexId,SubgraphId,SegmentId\n";

    // Loop over vertices of the PathGraph.
    BGL_FORALL_VERTICES(v, pathGraph, PathGraph) {
        const PathGraphVertex& vertex = pathGraph[v];

        // Write the AssemblyGraph path corresponding to this vertex.
        for(const uint64_t segmentId: vertex.path) {
            csv << vertex.id << "," << vertex.subgraphId << "," << segmentId << "\n";
        }
    }
}



// Detangling of a subgraph.
// Returns new vertices for the next detangle iteration.
// The new vertices can only be used in a new PathGraph
// created from scratch.
// Only the path and journeyIntervals are filled in.
void PathGraph::detangleSubgraph(
    uint64_t subgraphId,
    vector<PathGraphVertex> newVertices,
    bool debug
)
{
    const vector<vertex_descriptor>& subgraph = subgraphs[subgraphId];

    // Call the templated function appropriate for the
    // size of this subgraph. This way we use the shortest possible
    // bitmap (with size multiple of 64).
    if(subgraph.size() <= 64) {
        detangleSubgraphTemplate<64>(subgraph, newVertices, debug);
    } else if(subgraph.size() <= 128) {
        detangleSubgraphTemplate<128>(subgraph, newVertices, debug);
    } else if(subgraph.size() <= 192) {
        detangleSubgraphTemplate<192>(subgraph, newVertices, debug);
    } else if(subgraph.size() <= 256) {
        detangleSubgraphTemplate<256>(subgraph, newVertices, debug);
    } else if(subgraph.size() <= 320) {
        detangleSubgraphTemplate<320>(subgraph, newVertices, debug);
    } else if(subgraph.size() <= 384) {
        detangleSubgraphTemplate<384>(subgraph, newVertices, debug);
    } else if(subgraph.size() <= 448) {
        detangleSubgraphTemplate<448>(subgraph, newVertices, debug);
    } else if(subgraph.size() <= 512) {
        detangleSubgraphTemplate<512>(subgraph, newVertices, debug);
    } else {
        SHASTA_ASSERT(0);
    }
}


// This code is similar to mode3::AssemblyGraph::analyzeSubgraphTemplate
// but it operates on a subgraph of the PathGraph, not of the AssemblyGraph.
template<uint64_t N> void PathGraph::detangleSubgraphTemplate(
    const vector<vertex_descriptor>& subgraph,
    vector<PathGraphVertex> newVertices,
    bool debug
)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double fractionThreshold = 0.05;
    const uint64_t minVertexCoverage = 6;
    const uint64_t minClusterCoverage = 6;

    const PathGraph& pathGraph = *this;

    // The bitmap type used to store which vertices are visited
    // by each journey snippet.
    using BitVector = std::bitset<N>;
    SHASTA_ASSERT(subgraph.size() <= N);

    if(debug) {
        cout << "Detangling a PathGraph subgraph consisting of the following " <<
            subgraph.size() << " vertices:" << endl;
        for(const vertex_descriptor v: subgraph) {
            cout << pathGraph[v].id << " ";
        }
        cout << endl;
    }

    // Sanity check: we expect the vertices in the subgraph to be sorted by vertex id.
    SHASTA_ASSERT(std::is_sorted(subgraph.begin(), subgraph.end(),
        PathGraphOrderVerticesById(pathGraph)));

    // For vertices in the subgraph, gather triplets
    // (orientedReadId, position in path graph journey, vertex_descriptor).
    using Triplet = tuple<OrientedReadId, uint64_t, vertex_descriptor>;
    vector<Triplet> triplets;
    for(const vertex_descriptor v: subgraph) {
        const PathGraphVertex& vertex = pathGraph[v];

        // Loop over oriented reads that visit this vertex.
        for(const pair<AssemblyGraphJourneyInterval, uint64_t>& p: vertex.journeyIntervals) {
            const AssemblyGraphJourneyInterval& assemblyGraphJourneyInterval = p.first;
            const uint64_t position = p.second;
            const OrientedReadId orientedReadId = assemblyGraphJourneyInterval.orientedReadId;
            triplets.push_back(Triplet(orientedReadId, position, v));
        }
    }
    sort(triplets.begin(), triplets.end());

    // Write the triplets.
    if(debug) {
        ofstream csv("Triplets.csv");
        for(const Triplet& triplet: triplets) {
            csv << get<0>(triplet) << ",";
            csv << get<1>(triplet) << ",";
            csv << pathGraph[get<2>(triplet)].id << "\n";
        }
    }



    // Find streaks for the same OrientedReadId where the position
    // increases by 1 each time.
    // Each streak generates a PathGraphJourneySnippet.
    vector<PathGraphJourneySnippet> snippets;
    for(uint64_t i=0; i<triplets.size(); /* Increment later */) {
        const OrientedReadId orientedReadId = get<0>(triplets[i]);

        // Find this streak.
        uint64_t streakBegin = i;
        uint64_t streakEnd = streakBegin + 1;
        for(; streakEnd<triplets.size(); streakEnd++) {
            if(get<0>(triplets[streakEnd]) != orientedReadId) {
                break;
            }
            if(get<1>(triplets[streakEnd]) != get<1>(triplets[streakEnd-1]) + 1) {
                break;
            }
        }

        // Add a snippet.
        PathGraphJourneySnippet snippet;
        snippet.orientedReadId = orientedReadId;
        snippet.firstPosition = get<1>(triplets[streakBegin]);
        for(uint64_t j=streakBegin; j!=streakEnd; ++j) {
            snippet.vertices.push_back(get<2>(triplets[j]));
        }
        snippets.push_back(snippet);

        // Prepare to process the next streak.
        i = streakEnd;
    }




    // Write the snippets.
    if(debug) {
        ofstream csv("PathGraphJourneySnippets.csv");
        csv << "SnippetIndex,OrientedReadId,First position,LastPosition,Vertices\n";
        for(uint64_t snippetIndex=0; snippetIndex<snippets.size(); snippetIndex++) {
            const PathGraphJourneySnippet& snippet = snippets[snippetIndex];
            csv << snippetIndex << ",";
            csv << snippet.orientedReadId << ",";
            csv << snippet.firstPosition << ",";
            csv << snippet.lastPosition() << ",";
            for(const vertex_descriptor v: snippet.vertices) {
                csv << pathGraph[v].id << ",";
            }
            csv << "\n";
        }
    }



    // For each snippet, create a BitVector that describes the segments
    // the snippet visits.
    const uint64_t snippetCount = snippets.size();
    vector<BitVector> bitVectors(snippetCount);
    vector<uint64_t> bitVectorsPopCount(snippetCount);  // The number of bits set in each of the bit vectors.
    for(uint64_t snippetIndex=0; snippetIndex<snippetCount; snippetIndex++) {
        const PathGraphJourneySnippet& snippet = snippets[snippetIndex];
        BitVector& bitVector = bitVectors[snippetIndex];

        for(const vertex_descriptor v: snippet.vertices) {
            auto it = lower_bound(subgraph.begin(), subgraph.end(), v, PathGraphOrderVerticesById(pathGraph));
            SHASTA_ASSERT(it != subgraph.end());
            SHASTA_ASSERT(*it == v);
            const uint64_t bitIndex = it - subgraph.begin();
            bitVector.set(bitIndex);
        }
        bitVectorsPopCount[snippetIndex] = bitVector.count();
    }



    // Create the SnippetGraph.
    // A vertex represents a set of snippets and stores
    // the corresponding snippet indexes.
    // An edge x->y is created if there is at least one snippet in y
    // that is an approximate subset of a snippet in x.
    // We express this condition as |y-x| < fractionThreshold * |y|
    // We start with one snippet per vertex.
    SnippetGraph snippetGraph;
    vector<SnippetGraph::vertex_descriptor> vertexTable;
    std::map<SnippetGraph::vertex_descriptor, uint64_t> vertexMap;
    for(uint64_t snippetIndex=0; snippetIndex<snippetCount; snippetIndex++) {
        const auto v = add_vertex(SnippetGraphVertex(snippetIndex), snippetGraph);
        vertexTable.push_back(v);
        vertexMap.insert(make_pair(v, snippetIndex));
    }
    for(uint64_t iy=0; iy<snippetCount; iy++) {
        const BitVector& y = bitVectors[iy];
        const uint64_t threshold = uint64_t(std::round(fractionThreshold * double(bitVectorsPopCount[iy])));
        const SnippetGraph::vertex_descriptor vy = vertexTable[iy];
        for(uint64_t ix=0; ix<snippetCount; ix++) {
            if(ix == iy) {
                continue;
            }
            const BitVector& x = bitVectors[ix];

            // Compute z = y-x.
            BitVector z = y;
            z &= (~x);

            if(z.count() <= threshold) {
                const SnippetGraph::vertex_descriptor vx = vertexTable[ix];
                add_edge(vx, vy, snippetGraph);
            }
        }
    }



    // Compute strongly connected components of the SnippetGraph.
    std::map<SnippetGraph::vertex_descriptor, uint64_t> componentMap;
    const uint64_t componentCount = boost::strong_components(
        snippetGraph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexMap)));
    // cout << "Found " << componentCount << " strongly connected components." << endl;

    // Gather the vertices of each strongly connected component.
    vector< vector<SnippetGraph::vertex_descriptor> > components(componentCount);
    BGL_FORALL_VERTICES_T(v, snippetGraph, SnippetGraph) {
        const uint64_t componentId = componentMap[v];
        SHASTA_ASSERT(componentId < componentCount);
        components[componentId].push_back(v);
    }
    if(false) {
        cout << "Strongly connected components:\n";
        for(uint64_t componentId=0; componentId<componentCount; componentId++) {
            cout << componentId << ": ";
            for(const SnippetGraph::vertex_descriptor v: components[componentId]) {
                cout << vertexMap[v] << " ";
            }
            cout << "\n";
        }
    }



    // Condense the strongly connected components.
    // After this, the SnippetGraph is guaranteed to be acyclic.
    for(const vector<SnippetGraph::vertex_descriptor>& component: components) {
        if(component.size() == 1) {
            continue;
        }

        // Create a new vertex to represent this component.
        const auto vNew = add_vertex(snippetGraph);
        vector<uint64_t>& snippetsNew = snippetGraph[vNew].snippetIndexes;
        for(const vertex_descriptor v: component) {
            const vector<uint64_t>& snippets = snippetGraph[v].snippetIndexes;
            SHASTA_ASSERT(snippets.size() == 1);
            snippetsNew.push_back(snippets.front());
        }

        // Create the new edges.
        for(const vertex_descriptor v0: component) {

            // Out-edges.
            BGL_FORALL_OUTEDGES_T(v0, e01, snippetGraph, SnippetGraph) {
                const vertex_descriptor v1 = target(e01, snippetGraph);
                if(v1 != vNew) {
                    add_edge(vNew, v1, snippetGraph);
                }
            }

            // In-edges.
            BGL_FORALL_INEDGES_T(v0, e10, snippetGraph, SnippetGraph) {
                const vertex_descriptor v1 = source(e10, snippetGraph);
                if(v1 != vNew) {
                    add_edge(v1, vNew, snippetGraph);
                }
            }
        }

        // Remove the old vertices and their edges.
        for(const vertex_descriptor v: component) {
            clear_vertex(v, snippetGraph);
            remove_vertex(v, snippetGraph);
        }
    }



    // Compute which maximal vertices each vertex is a descendant of.
    std::map<SnippetGraph::vertex_descriptor, vector<SnippetGraph::vertex_descriptor> > ancestorMap;
    BGL_FORALL_VERTICES_T(v0, snippetGraph, SnippetGraph) {
        if(in_degree(v0, snippetGraph) != 0) {
            continue;   // Not a maximal vertex.
        }

        // Find the descendants of this maximal vertex.
        vector<vertex_descriptor> descendants;
        snippetGraph.findDescendants(v0, descendants);

        // Update the ancestor map.
        for(const vertex_descriptor v1: descendants) {
            ancestorMap[v1].push_back(v0);
        }
    }



    // Each maximal vertex generates a cluster consisting of the vertices
    // that descend from it and from no other maximal vertex.
    // Gather the vertices in each cluster.
    std::map<SnippetGraph::vertex_descriptor, vector<SnippetGraph::vertex_descriptor> > clusterMap;
    uint64_t unclusterVertexCount = 0;
    BGL_FORALL_VERTICES_T(v1, snippetGraph, SnippetGraph) {
        const vector<SnippetGraph::vertex_descriptor>& ancestors = ancestorMap[v1];
        if(ancestors.size() == 1) {
            const vertex_descriptor v0 = ancestors.front();
            clusterMap[v0].push_back(v1);
        } else {
            ++unclusterVertexCount;
        }
    }
    cout << "Found " << unclusterVertexCount << " unclustered vertices." << endl;




    // Gather the snippets in each cluster.
    vector<PathGraphJourneySnippetCluster> clusters;
    for(const auto& p: clusterMap) {
        const vector<SnippetGraph::vertex_descriptor>& clusterVertices = p.second;
        clusters.resize(clusters.size() + 1);
        PathGraphJourneySnippetCluster& cluster = clusters.back();

        vector<uint64_t> clusterSnippetIndexes; // Only used for debug output.
        for(const SnippetGraph::vertex_descriptor v: clusterVertices) {
            const vector<uint64_t>& snippetIndexes = snippetGraph[v].snippetIndexes;
            for(const uint64_t snippetIndex: snippetIndexes) {
                cluster.snippets.push_back(snippets[snippetIndex]);
                clusterSnippetIndexes.push_back(snippetIndex);
            }
        }
        cluster.constructVertices(pathGraph);
        cluster.cleanupVertices(minVertexCoverage);
        cout << "Found a cluster candidate with " <<
            clusterVertices.size() << " vertices and " <<
            cluster.snippets.size() << " snippets:" << endl;
        for(const uint64_t snippetIndex: clusterSnippetIndexes) {
            cout << snippetIndex << " ";
        }
        cout << endl;

        // If coverage on this cluster is too low, discard it.
        if(cluster.coverage() < minClusterCoverage) {
            clusters.resize(clusters.size() - 1);
            cout << "This cluster candidate was discarded because of low coverage." << endl;
            continue;
        }

        // This cluster will be stored and is assigned this clusterId;
        const uint64_t clusterId = clusters.size() - 1;

        if(debug) {

            cout << "This cluster was stored as cluster " << clusterId << endl;
            cout << "Vertex(coverage) for this cluster:\n";
            for(const auto& p: cluster.vertices) {
                cout << pathGraph[p.first].id << "(" << p.second << ") ";
            }
            cout << endl;
        }

        // Mark the vertices of this cluster.
        for(const SnippetGraph::vertex_descriptor v: clusterVertices) {
            snippetGraph[v].clusterId = clusterId;
        }
    }
    snippetGraph.clusterCount = clusters.size();




    // Write out the SnippetGraph.
    if(debug) {
        snippetGraph.writeGraphviz("SnippetGraph.dot");
    }



    // Find the paths of each cluster.
    // Each of these paths generates a new vertex for the next detangle iteration.
    newVertices.clear();
    cout << "Kept " << clusters.size() << " clusters." << endl;
    for(uint64_t clusterId=0; clusterId<clusters.size(); clusterId++) {
        const PathGraphJourneySnippetCluster& cluster = clusters[clusterId];
        vector< vector<vertex_descriptor> > paths;
        ofstream graphOut("Cluster-" + to_string(clusterId) + ".dot");
        cout << "Finding paths generates by cluster " << clusterId << endl;
        findClusterPaths(cluster, paths, graphOut);

        // For each path, generate a new vertex for the next detangle iteration.
        for(const vector<vertex_descriptor>& path: paths) {
            newVertices.emplace_back();
            PathGraphVertex& newVertex = newVertices.back();

            // Construct the assembly graph path for the new vertex.
            for(const vertex_descriptor v: path) {
                const PathGraphVertex& vertex = pathGraph[v];
                copy(vertex.path.begin(), vertex.path.end(), back_inserter(newVertex.path));
            }

            // Construct the AssemblyGraphJourneyInterval for the new vertex.
        }
    }
}



void SnippetGraph::findDescendants(
    const vertex_descriptor vStart,
    vector<vertex_descriptor>& descendants) const
{
    const SnippetGraph& graph = *this;

    // Initialize the BFS.
    std::queue<vertex_descriptor> q;
    q.push(vStart);
    std::set<vertex_descriptor> descendantsSet;
    descendantsSet.insert(vStart);

    // BFS loop.
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        BGL_FORALL_OUTEDGES(v0, e01, graph, SnippetGraph) {
            const vertex_descriptor v1 = target(e01, graph);
            if(descendantsSet.find(v1) == descendantsSet.end()) {
                q.push(v1);
                descendantsSet.insert(v1);
            }
        }
    }

    descendants.clear();
    copy(descendantsSet.begin(), descendantsSet.end(), back_inserter(descendants));
}



void SnippetGraph::writeGraphviz(
    const string& fileName) const
{
    const SnippetGraph& graph = *this;

    ofstream dot(fileName);
    dot << "digraph SnippetGraph{\n"
        "node [shape=rectangle];\n";
    BGL_FORALL_VERTICES(v, graph, SnippetGraph) {
        dot << "\"" << v << "\" [label=\"";
        const vector<uint64_t>& snippetIndexes = graph[v].snippetIndexes;
        for(const uint64_t snippetIndex: snippetIndexes) {
            dot << snippetIndex;
            dot << "\\n";
        }
        dot << "\"";
        const uint64_t clusterId = graph[v].clusterId;
        if(clusterId != std::numeric_limits<uint64_t>::max()) {
            dot << " style=filled fillcolor=\"" <<
                float(clusterId)/float(clusterCount) <<
                ",0.3,1\"";
        }
        dot << "];\n";
    }
    BGL_FORALL_EDGES(e, graph, SnippetGraph) {
        const vertex_descriptor vx = source(e, graph);
        const vertex_descriptor vy = target(e, graph);
        dot << "\"" << vx << "\"->\"" << vy << "\";\n";
    }
    dot << "}\n";

}



vector<PathGraphBaseClass::vertex_descriptor> PathGraphJourneySnippetCluster::getVertices() const
{
    vector<PathGraphBaseClass::vertex_descriptor> v;
    for(const auto& p: vertices) {
        v.push_back(p.first);
    }
    return v;
}



void PathGraphJourneySnippetCluster::cleanupVertices(uint64_t minVertexCoverage)
{
    vector< pair<PathGraphBaseClass::vertex_descriptor, uint64_t > > newVertices;
    for(const auto& p: vertices) {
        if(p.second >= minVertexCoverage) {
            newVertices.push_back(p);
        }
    }
    vertices.swap(newVertices);
}



void PathGraphJourneySnippetCluster::constructVertices(const PathGraph& pathGraph)
{
    // A map with Key=vertex_descriptor, value = coverage.
    auto vertexMap = std::map<PathGraphBaseClass::vertex_descriptor, uint64_t, PathGraphOrderVerticesById>(
        PathGraphOrderVerticesById(pathGraph));

    for(const PathGraphJourneySnippet& snippet: snippets) {
        for(const PathGraphBaseClass::vertex_descriptor v: snippet.vertices) {
            auto it = vertexMap.find(v);
            if(it == vertexMap.end()) {
                vertexMap.insert(make_pair(v, 1));
            } else {
                ++(it->second);
            }
        }
    }

    vertices.clear();
    copy(vertexMap.begin(), vertexMap.end(), back_inserter(vertices));
}



// Given a PathGraphJourneySnippetCluster, find a plausible
// path for it in the PathGraph.
void PathGraph::findClusterPaths(
    const PathGraphJourneySnippetCluster& cluster,
    vector< vector<vertex_descriptor> >& paths,
    ostream& graphOut) const
{
    const PathGraph& pathGraph = *this;
    const bool debug = true;
    paths.clear();

    // Map vertex descriptors to indexes in cluster.vertices.
    std::map<vertex_descriptor, uint64_t> vertexMap;
    for(uint64_t i=0; i<cluster.vertices.size(); i++) {
        const vertex_descriptor v = cluster.vertices[i].first;
        vertexMap.insert(make_pair(v, i));
    }

    // Construct the subgraph induced by the vertices of the cluster.
    using Subgraph = boost::adjacency_list<boost::listS, boost::vecS, boost::bidirectionalS>;
    Subgraph subgraph(vertexMap.size());
    for(const auto& p: vertexMap) {
        const vertex_descriptor v0 = p.first;
        const uint64_t i0 = p.second;
        BGL_FORALL_OUTEDGES(v0, e, pathGraph, PathGraph) {
            const vertex_descriptor v1 = target(e, pathGraph);
            const auto it = vertexMap.find(v1);
            if(it == vertexMap.end()) {
                continue;
            }
            const uint64_t i1 = it->second;
            add_edge(i0, i1, subgraph);
        }
    }

    // Compute strong connected components of this subgraph.
    const auto indexMap = get(boost::vertex_index, subgraph);
    vector<uint64_t> strongComponent(num_vertices(subgraph));
    boost::strong_components(
        subgraph,
        boost::make_iterator_property_map(strongComponent.begin(), indexMap));

    // Remove edges internal to strong components.
    vector<Subgraph::edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, subgraph, Subgraph) {
        const uint64_t i0 = source(e, subgraph);
        const uint64_t i1 = target(e, subgraph);
        if(strongComponent[i0] == strongComponent[i1]) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const Subgraph::edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, subgraph);
    }

    // Transitive reduction.
    transitiveReduction(subgraph);


    // Write it out.
    if(debug) {
        graphOut << "digraph cluster {\n";
        for(uint64_t i=0; i<vertexMap.size(); i++) {
            const auto& p = cluster.vertices[i];
            const vertex_descriptor v = p.first;
            const uint64_t coverage = p.second;
            graphOut << pathGraph[v].id;
            graphOut << " [label=\"" << pathGraph[v].id << "\\n" << coverage << "\"]";
            graphOut << ";\n";
        }
        BGL_FORALL_EDGES(e, subgraph, Subgraph) {
            const uint64_t i0 = source(e, subgraph);
            const uint64_t i1 = target(e, subgraph);
            const vertex_descriptor v0 = cluster.vertices[i0].first;
            const vertex_descriptor v1 = cluster.vertices[i1].first;
            graphOut << pathGraph[v0].id << "->" << pathGraph[v1].id;
            graphOut << ";\n";
        }
        graphOut << "}\n";

    }


    // Find linear chains of vertices.
    vector< vector<Subgraph::vertex_descriptor> > chains;
    findLinearVertexChains(subgraph, chains);
    if(debug) {
        cout << "Found the following paths:" << endl;
        for(const vector<Subgraph::vertex_descriptor>& chain: chains) {
            for(const Subgraph::vertex_descriptor v: chain) {
                const PathGraph::vertex_descriptor u = cluster.vertices[v].first;
                cout << pathGraph[u].id << " ";
            }
            cout << endl;
        }
    }
}
