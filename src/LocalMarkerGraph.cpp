// Nanopore2
#include "LocalMarkerGraph.hpp"
#include "LongBaseSequence.hpp"
#include "Marker.hpp"
#include "orderPairs.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Boost libraries.
// Due to a bug in boost graph include files, disjoint_sets.hpp
// must be included before graphviz.hpp
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/topological_sort.hpp>

// Standard library.
#include "algorithm.hpp"
#include "fstream.hpp"
#include <set>
#include "tuple.hpp"



// Construct the graph given the markers of
// given oriented reads, sorted by position.
// Each pair in the second input vector is (position, KmerId).
// This constructs the initial graph without any
// alignment information, consisting of a
// separate linear sequence for each input oriented read.
LocalMarkerGraph::LocalMarkerGraph(
    size_t k,
    const vector<OrientedReadId>& orientedReadIds,
    const vector<LongBaseSequence>& sequences,
    const vector< vector<Marker0> >& markers0,
    size_t minCoverage,     // For a vertex to be considered strong.
    size_t minConsensus     // For an edge to be considered strong.
    ) :
    k(k),
    orientedReadIds(orientedReadIds),
    sequences(sequences),
    markers0(markers0),
    minCoverage(minCoverage),
    minConsensus(minConsensus)
{

    // Sanity check on the input.
    CZI_ASSERT(orientedReadIds.size() == markers0.size());

    // Shorthand for readability.
    LocalMarkerGraph& graph = *this;

    // Prepare the vertex map. It stores a vector of vertex descriptor
    // for each of the input reads.
    vertexMap.resize(markers0.size());

    // Loop over all input oriented reads to generate
    // a linear chain for each.
    for(size_t localOrientedReadId=0; localOrientedReadId<markers0.size(); localOrientedReadId++) {

        // Access information about this oriented read.
        const auto& orientedReadMarkers = markers0[localOrientedReadId];
        auto& orientedReadVertexMap = vertexMap[localOrientedReadId];
        orientedReadVertexMap.reserve(markers0.size());

        // Generate the vertices for this oriented read.
        for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
            const Marker0& marker = orientedReadMarkers[ordinal];
            CZI_ASSERT(marker.ordinal == ordinal);
            const KmerId kmerId = marker.kmerId;

            const vertex_descriptor v = boost::add_vertex(graph);

            orientedReadVertexMap.push_back(v);
            LocalMarkerGraphVertex& vertex = graph[v];
            vertex.kmerId = kmerId;
            vertex.markerIds.push_back(LocalMarkerGraphVertex::MarkerId(uint32_t(localOrientedReadId), ordinal));
            vertex.vertexId = nextVertexId++;
        }

        // Generate the edges.
        for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size()-1; ordinal++) {
            const vertex_descriptor v0 = orientedReadVertexMap[ordinal];
            const vertex_descriptor v1 = orientedReadVertexMap[ordinal + 1];
            boost::add_edge(v0, v1, graph);
        }
    }

}


// Merge two vertices.
// The vertices to be merged are specified by
// the local oriented read ids (indexes in the orientedReadIds vector)
// and the ordinals.
void LocalMarkerGraph::mergeVertices(
    uint32_t localOrientedReadId0,
    uint32_t ordinal0,
    uint32_t localOrientedReadId1,
    uint32_t ordinal1
)
{
    // Shorthand for readability.
    LocalMarkerGraph& graph = *this;

    // Get the vertices to be merged.
    // We will keep v0 and remove v1.
    const vertex_descriptor v0 = vertexMap[localOrientedReadId0][ordinal0];
    const vertex_descriptor v1 = vertexMap[localOrientedReadId1][ordinal1];
    if(v0 == v1) {
        return;
    }
    LocalMarkerGraphVertex& vertex0 = graph[v0];
    const LocalMarkerGraphVertex& vertex1 = graph[v1];

    // Copy the marker ids from v1 to v0.
    vertex0.markerIds.insert(
        vertex0.markerIds.end(),
        vertex1.markerIds.begin(),
        vertex1.markerIds.end());

    // Reroute the edges.
    std::set<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_OUTEDGES(v1, e12, graph, LocalMarkerGraph) {
        const vertex_descriptor v2 = target(e12, graph);
        boost::add_edge(v0, v2, graph);
        edgesToBeRemoved.insert(e12);
    }
    BGL_FORALL_INEDGES(v1, e21, graph, LocalMarkerGraph) {
        const vertex_descriptor v2 = source(e21, graph);
        boost::add_edge(v2, v0, graph);
        edgesToBeRemoved.insert(e21);
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }

    // Remove the vertex but keep track (using the vertex map)
    // of where it ended up.
    for(const auto& markerId: vertex1.markerIds) {
        vertexMap[markerId.localOrientedReadId][markerId.ordinal] = v0;
    }
    boost::clear_vertex(v1, graph);
    boost::remove_vertex(v1, graph);


}



// Remove vertices with coverage less than minCoverage.
void LocalMarkerGraph::removeWeakVertices()
{
    LocalMarkerGraph& graph = *this;
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        if(graph[v].markerIds.size() < minCoverage) {
            verticesToBeRemoved.push_back(v);
        }
    }

    for(const vertex_descriptor v: verticesToBeRemoved) {
        const LocalMarkerGraphVertex& vertex = graph[v];
        for(const auto& markerId: vertex.markerIds) {
            vertexMap[markerId.localOrientedReadId][markerId.ordinal] = null_vertex();
        }
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }

}



// Remove edges with consensus less than the minCconsensus
// (consensus is maximum number of oriented reads that agree on the sequence).
void LocalMarkerGraph::removeWeakEdges()
{
    LocalMarkerGraph& graph = *this;
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph) {
        if(edgeConsensus(e) < minConsensus) {
            edgesToBeRemoved.push_back(e);
        }
    }

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }

}



// Same as above, but keep all edges marked as spanning tree edges.
void LocalMarkerGraph::removeWeakNonSpanningTreeEdges()
{
    LocalMarkerGraph& graph = *this;
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph) {
        if( !graph[e].isInOptimalSpanningTree &&
            edgeConsensus(e) < minConsensus) {
            edgesToBeRemoved.push_back(e);
        }
    }

    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }

}



// Recursively prune weak leaves.
// This is a simple minded iterative implementation.
// It could be done faster recursively.
void LocalMarkerGraph::pruneWeakLeaves()
{
    LocalMarkerGraph& graph = *this;

    vector<vertex_descriptor> verticesToBeRemoved;
    while(true) {
        verticesToBeRemoved.clear();
        BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
            if(graph[v].markerIds.size() < minCoverage) {
                if(boost::in_degree(v, graph)==0 || boost::out_degree(v, graph)==0) {
                    verticesToBeRemoved.push_back(v);
                }
            }
        }

        if(verticesToBeRemoved.empty()) {
            break;
        }

        for(const vertex_descriptor v: verticesToBeRemoved) {
            boost::clear_vertex(v, graph);
            boost::remove_vertex(v, graph);
        }

    }


}



// Sort the marker ids in each vertex.
void LocalMarkerGraph::sort()
{
    LocalMarkerGraph& graph = *this;
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        auto& markerIds = graph[v].markerIds;
        std::sort(markerIds.begin(), markerIds.end());
    }
}



// Sort the occurrences in each vertex
// and remove ambiguous vertices.
// A vertex is ambiguous if it has more
// than one marker for at least one oriented read.
void LocalMarkerGraph::sortAndRemoveAmbiguousVertices()
{
    LocalMarkerGraph& graph = *this;

    // Make sure the marker ids in each vertex are sorted
    // by oriented read.
    sort();

    // Remove ambiguous vertices.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        auto& markerIds = graph[v].markerIds;
        if(markerIds.size() < 2) {
            continue;
        }
        for(auto it=markerIds.begin()+1; it!=markerIds.end(); ++it) {
            if(it->localOrientedReadId == (it-1)->localOrientedReadId) {
                verticesToBeRemoved.push_back(v);
                break;
            }
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        // Something is missing here.
        // We have to update the vertex map.
        // This function is currently not used.
        CZI_ASSERT(0);
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }
}



// Sort the occurrences in each vertex
// and split ambiguous vertices.
// A vertex is ambiguous if it has more
// than one marker for at least one oriented read.
// Ambiguous vertices can occur as a result of k-mers
// occurring in two or more neighboring copies.
void LocalMarkerGraph::sortAndSplitAmbiguousVertices()
{
    LocalMarkerGraph& graph = *this;

    // Make sure the occurrences in each vertex are sorted
    // by oriented read.
    sort();

    // Split ambiguous vertices.
    vector<vertex_descriptor> verticesToBeSplit;
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        auto& markerIds = graph[v].markerIds;
        if(markerIds.size() < 2) {
            continue;
        }
        for(auto it=markerIds.begin()+1; it!=markerIds.end(); ++it) {
            if(it->localOrientedReadId == (it-1)->localOrientedReadId) {
                verticesToBeSplit.push_back(v);
                break;
            }
        }
    }
    for(const vertex_descriptor v: verticesToBeSplit) {
        split(v);
    }
}



// Split a vertex, creating a separate new vertex for each marker.
void LocalMarkerGraph::split(vertex_descriptor vOld)
{

    LocalMarkerGraph& graph = *this;
    const LocalMarkerGraphVertex& oldVertex = graph[vOld];

    // Generate a new vertex for each marker in this vertex.
    const auto& markerIds = graph[vOld].markerIds;
    for(const auto& markerId: markerIds) {
        const uint32_t localOrientedReadId = markerId.localOrientedReadId;
        const uint32_t ordinal = markerId.ordinal;
        const vertex_descriptor vNew = add_vertex(graph);
        CZI_ASSERT(vertexMap[localOrientedReadId][ordinal] == vOld);
        vertexMap[localOrientedReadId][ordinal] = vNew;
        LocalMarkerGraphVertex& newVertex = graph[vNew];
        newVertex.kmerId = oldVertex.kmerId;
        newVertex.vertexId = nextVertexId++;
        newVertex.markerIds.push_back(markerId);

        // Add the edges.
        if(ordinal != markers0[localOrientedReadId].size()-1) {
            add_edge(vNew, vertexMap[localOrientedReadId][ordinal+1], graph);
        }
        if(ordinal != 0) {
            add_edge(vertexMap[localOrientedReadId][ordinal-1], vNew, graph);
        }

    }

    // Remove the vertex we split.
    boost::clear_vertex(vOld, graph);
    boost::remove_vertex(vOld, graph);
}



// Fill in edge data.
// This assumes that there are no ambiguous vertices
// and that the occurrences vectors in the vertices are sorted.
// This must be done before calling write.
// This only fills in edge data when the
// ordinal in the source and target vertices
// are contiguous.
void LocalMarkerGraph::fillEdgeData()
{
    BGL_FORALL_EDGES(e, *this, LocalMarkerGraph) {
        fillEdgeData(e);
    }
}
void LocalMarkerGraph::fillEdgeData(edge_descriptor e)
{
    // Access the graph and the edge.
    LocalMarkerGraph& graph = *this;
    LocalMarkerGraphEdge& edge = graph[e];

    // Access the vertices.
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);
    const LocalMarkerGraphVertex& vertex0 = graph[v0];
    const LocalMarkerGraphVertex& vertex1 = graph[v1];

    // Access the markers of each vertex.
    // They are assumed to be sorted and without duplicate
    // oriented read ids.
    const auto& markerIds0 = vertex0.markerIds;
    const auto& markerIds1 = vertex1.markerIds;

    // Joint loop over markers, looking for common oriented read ids.
    auto it0 = markerIds0.begin();
    auto it1 = markerIds1.begin();
    while(it0!=markerIds0.end() && it1!=markerIds0.end()) {
        if(it0->localOrientedReadId < it1->localOrientedReadId) {
            ++it0;
            continue;
        } else if(it1->localOrientedReadId < it0->localOrientedReadId) {
            ++it1;
            continue;
        }

        // Only add it if the ordinals are contiguous.
        if(it1->ordinal == it0->ordinal + 1) {

            // If getting here, we found a common oriented read id.
            const uint32_t localOrientedReadId = it0->localOrientedReadId;
            CZI_ASSERT(localOrientedReadId == it1->localOrientedReadId);
            edge.data.push_back(LocalMarkerGraphEdge::Data(localOrientedReadId, it0->ordinal, it1->ordinal));
        }

        // Continue the joint loop.
        // This assumes that there are no duplicate local oriented read ids
        // (that is, the vertices are not amboguous).
        ++it0;
        ++it1;
    }

    // Because the local oriented read ids in the vertices
    // are sorted, edge data is also sorted by localOrientedReadId
    // and there is no need to sort here.
}



// Compute the consensus of an edge.
// It is the maximum number of oriented read ids
// that agree on the sequence associated with the edge.
size_t LocalMarkerGraph::edgeConsensus(edge_descriptor e) const
{
    const LocalMarkerGraph& graph = *this;
    const LocalMarkerGraphEdge& edge = graph[e];

    // Get the sequence strings.
    std::map<string, vector<LocalMarkerGraphEdge::Data> > sequenceMap;
    for(const auto& data: edge.data) {
        sequenceMap[graph.getSequence(data)].push_back(data);
    }

    // Find the maximum number of oriented reads that agree.
    size_t consensus = 0;
    for(const auto& p: sequenceMap) {
        consensus = max(consensus, p.second.size());
    }
    return consensus;
}



// Add edges that were removed because of low coverage,
// skipping ordinals if necessary (we can only use surviving vertices).
void LocalMarkerGraph::addMissingEdges()
{
    LocalMarkerGraph& graph = *this;
    BGL_FORALL_VERTICES(v0, graph, LocalMarkerGraph) {
        const LocalMarkerGraphVertex& vertex0 = graph[v0];
        for(const auto& markerId: vertex0.markerIds) {
            const uint32_t localOrientedReadId = markerId.localOrientedReadId;
            const auto& orientedReadMarkers = markers0[localOrientedReadId];
            const auto& orientedReadVertexMap = vertexMap[localOrientedReadId];
            const uint32_t ordinal0 = markerId.ordinal;
            for(uint32_t ordinal1=ordinal0+1; ordinal1!=orientedReadMarkers.size(); ordinal1++) {

                // Get the vertex corresponding to this ordinal.
                const vertex_descriptor v1 = orientedReadVertexMap[ordinal1];
                if(v1 == null_vertex()) {
                    // The vertex does not exist. Try the next ordinal.
                    continue;
                }

                // See if we already have this edge.
                edge_descriptor e;
                bool edgeExists;
                tie(e, edgeExists) = boost::edge(v0, v1, graph);

                // If the edge does not exists, create it.
                // Otherwise just add a Data item for this ordinal0, ordinal1.
                const LocalMarkerGraphEdge::Data edgeData(localOrientedReadId, ordinal0, ordinal1);
                if(!edgeExists) {
                    tie(e, ignore) = boost::add_edge(v0, v1, graph);
                    graph[e].data.push_back(edgeData);
                } else {
                    auto& existingData = graph[e].data;
                    if(find(existingData.begin(), existingData.end(), edgeData) == existingData.end()) {
                        existingData.push_back(edgeData);
                    }
                }
                break;
            }
        }
    }

}



// Construct the sequence implied by an LocalMarkerGraphEdge::Data.
// If the source and target k-mers overlap,
// the string is a number representing the number of overlap bases.
// Otherwise it is a base sequence.
// We should phase this out in favor of getEdgeSequence.
string LocalMarkerGraph::getSequence(const LocalMarkerGraphEdge::Data& data) const
{
    const uint32_t localOrientedReadId = data.localOrientedReadId;
    const uint32_t ordinal0 = data.ordinals[0];
    const uint32_t ordinal1 = data.ordinals[1];

    const uint32_t beginPosition0 = markers0[localOrientedReadId][ordinal0].position;
    const uint32_t beginPosition1 = markers0[localOrientedReadId][ordinal1].position;
    const uint32_t endPosition0 = beginPosition0 + uint32_t(k);

    string s;
    if(beginPosition1 <= endPosition0) {
        s = to_string(endPosition0 - beginPosition1);
    } else {
        for(uint32_t position=endPosition0; position!=beginPosition1; position++) {
            s.push_back(sequences[localOrientedReadId][position].character());
        }
    }
    return s;
}



LocalMarkerGraph::EdgeSequence LocalMarkerGraph::getEdgeSequence(
    const LocalMarkerGraphEdge::Data& data) const
{
    const uint32_t localOrientedReadId = data.localOrientedReadId;
    const uint32_t ordinal0 = data.ordinals[0];
    const uint32_t ordinal1 = data.ordinals[1];

    const uint32_t beginPosition0 = markers0[localOrientedReadId][ordinal0].position;
    const uint32_t beginPosition1 = markers0[localOrientedReadId][ordinal1].position;
    const uint32_t endPosition0 = beginPosition0 + uint32_t(k);

    EdgeSequence s;
    if(beginPosition1 <= endPosition0) {
        s.overlap = endPosition0 - beginPosition1;
    } else {
        for(uint32_t position=endPosition0; position!=beginPosition1; position++) {
            s.sequence.push_back(sequences[localOrientedReadId][position]);
        }
    }
    return s;
}



pair<LocalMarkerGraph::EdgeSequence, int> LocalMarkerGraph::getDominantEdgeSequence(edge_descriptor e) const
{
    const LocalMarkerGraphEdge& edge = (*this)[e];

    // Gather the sequences.
    std::map<EdgeSequence, int> m;
    for(const auto& data: edge.data) {
        const EdgeSequence sequence = getEdgeSequence(data);
        const auto it = m.find(sequence);
        if(it == m.end()) {
            m.insert(make_pair(sequence, 1));
        } else {
            ++(it->second);
        }
    }
    CZI_ASSERT(!m.empty());

    // Find the dominant one.
    pair<EdgeSequence, int> p = *(m.begin());
    for(const auto& q: m) {
        if(q.second > p.second) {
            p = q;
        }
    }
    return p;
}



// Mark edges that are part of an optimal spanning tree.
// The optimal spanning tree is a minimum spanning tree
// (defined treating the graph as undirected)
// where the weight of each edge is the inverse of the consensus for the edge.
// We cannot use the minimum spanning tree algorithms
// from the boost graph library because
// they only work for undirected graph.
void LocalMarkerGraph::computeOptimalSpanningTree()
{
    LocalMarkerGraph& graph = *this;

    // Gather all the edges, sorted by decreasing consensus.
    // Also, mark all edge as initially not part of the optimal spanning tree.
    vector< pair<edge_descriptor, size_t> > edgeTable;
    BGL_FORALL_EDGES(e, graph, LocalMarkerGraph) {
        edgeTable.push_back(make_pair(e, edgeConsensus(e)));
        graph[e].isInOptimalSpanningTree = false;
    }
    std::sort(edgeTable.begin(), edgeTable.end(),
        OrderPairsBySecondOnlyGreater<edge_descriptor, size_t>());

    // Map the vertices to integers in [0, number of vertices).
    const size_t n = boost::num_vertices(graph);
    std::map<vertex_descriptor, uint32_t> vertexMap;
    uint32_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        vertexMap.insert(make_pair(v, vertexIndex++));
    }


    // Initialize the disjoint set data structures.
    vector<uint32_t> rank(n);
    vector<uint32_t> parent(n);
    boost::disjoint_sets<uint32_t*, uint32_t*> disjointSets(&rank[0], &parent[0]);
    for(size_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    // Process the edges in this order.
    // Only add each edge to the optimal spanning tree
    // if the two vertices are in two different connected components.
    for(const auto& p: edgeTable) {
        const edge_descriptor e = p.first;
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const uint32_t i0 = vertexMap[v0];
        const uint32_t i1 = vertexMap[v1];

        // If v0 and v1 are in separate components,
        // add this edge to the optimal spanning tree.
        const uint32_t component0 = disjointSets.find_set(i0);
        const uint32_t component1 = disjointSets.find_set(i1);
        if(component0 != component1) {
            graph[e].isInOptimalSpanningTree = true;
            disjointSets.union_set(i0, i1);
        }
    }

}



// Find the longest path in the graph.
// This assumes that the graph is acyclic.
// See here for a description of the algorithm:
// https://en.wikipedia.org/wiki/Longest_path_problem#Acyclic_graphs_and_critical_paths
vector<LocalMarkerGraph::edge_descriptor> LocalMarkerGraph::findLongestPath()
{
    LocalMarkerGraph& graph = *this;

    // Get the vertices in topological sort order.
    const vector<vertex_descriptor> topologicallySortedVertices = topologicalSort();

    // Construct the path. See the Wikipedia article linked above.
    // We use the vertex color to record distances.
    for(const vertex_descriptor v0: topologicallySortedVertices) {
        LocalMarkerGraphVertex& vertex0 = graph[v0];
        vertex0.color = 0;
        BGL_FORALL_INEDGES(v0, e10, graph, LocalMarkerGraph) {
            const vertex_descriptor v1 = source(e10, graph);
            const LocalMarkerGraphVertex& vertex1 = graph[v1];
            vertex0.color = max(vertex0.color, vertex1.color+1);
        }
    }

    // Find the vertex with the largest distance.
    // This will be the last vertex of the path.
    vertex_descriptor vLast = null_vertex();
    uint32_t maxDistance = 0;
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        if(vLast==null_vertex() || graph[v].color>maxDistance) {
            maxDistance = graph[v].color;
            vLast = v;
        }
    }


    // Now walk back to construct the path.
    vector<edge_descriptor> path;
    vertex_descriptor v = vLast;
    while(graph[v].color > 0) {
        int maxParentDistance = -1;
        edge_descriptor maxParentEdge;
        BGL_FORALL_INEDGES(v, e, graph, LocalMarkerGraph) {
            const vertex_descriptor parent = source(e, graph);
            if(int(graph[parent].color) > maxParentDistance) {
                maxParentDistance = int(graph[parent].color);
                maxParentEdge = e;
            }
        }
        path.push_back(maxParentEdge);
        v = source(maxParentEdge, graph);
    }
    reverse(path.begin(), path.end());

    return path;
}



// Topological sort of the graph.
// This assumes that the graph is acyclic.
vector<LocalMarkerGraph::vertex_descriptor> LocalMarkerGraph::topologicalSort()
{
    LocalMarkerGraph& graph = *this;

    vector<vertex_descriptor> topologicallySortedVertices;
    boost::topological_sort(graph, back_inserter(topologicallySortedVertices),
        boost::color_map(boost::get(&LocalMarkerGraphVertex::color, *this)));

    // Boost ::topological_sort returns the vertices in reverse topological.
    std::reverse(topologicallySortedVertices.begin(), topologicallySortedVertices.end());
    return topologicallySortedVertices;
}



// Extract the longest sequence present in the graph.
vector<pair<Nanopore2::Base, int> > LocalMarkerGraph::extractLongestSequence()
{

    // Find the longest path.
    const vector<LocalMarkerGraph::edge_descriptor> path = findLongestPath();
    CZI_ASSERT(!path.empty());

    // Write out the longest path.
    cout << "The longest path has " << path.size()+1 << " vertices." << endl;

    // Extract the sequence from this path.
    return getPathSequence(path);
}



// Extract the sequence corresponding to a path.
vector<pair<Nanopore2::Base, int> > LocalMarkerGraph::getPathSequence(
    const vector<edge_descriptor>& path)
{
    LocalMarkerGraph& graph = *this;
    using Nanopore2::Base;

    // Sanity check on the path.
    CZI_ASSERT(!path.empty());
    for(size_t i=1; i<path.size(); i++) {
        const edge_descriptor ePrevious = path[i-1];
        const edge_descriptor e = path[i];
        CZI_ASSERT(target(ePrevious, graph) == source(e, graph));
    }

    // Loop over vertices and edges of the path.
    // The number of vertices is one more than the number of edges.
    vector< pair<Base, int> > sequence;
    size_t position = 0;
    for(size_t i=0; i<=path.size(); i++) {

        // Get the i-th vertex.
        vertex_descriptor v;
        if(i < path.size()) {
            v = source(path[i], graph);
        } else {
            v = target(path.back(), graph);
        }
        const LocalMarkerGraphVertex& vertex = graph[v];
        const Kmer kmer(vertex.kmerId, k);

        // cout << "i=" << i << " vertex kmer: ";
        // kmer.write(cout, k);
        // cout << endl;

        // Update the sequence with the vertex sequence.
        const int vertexCoverage = int(vertex.markerIds.size());
        for(size_t j=0; j<k; j++) {
            const Base base = kmer[j];
            const size_t positionInSequence = position + j;
            if(positionInSequence < sequence.size()) {
                CZI_ASSERT(sequence[positionInSequence].first == base);
                sequence[positionInSequence].second = max(sequence[positionInSequence].second, vertexCoverage);
            } else {
                CZI_ASSERT(sequence.size() == positionInSequence);
                sequence.push_back(make_pair(base, vertexCoverage));
            }
        }

        // If this is the last vertex, we are done.
        if(i==path.size()) {
            break;
        }

        // Update the sequence with the dominant edge sequence.
        const edge_descriptor e = path[i];
        EdgeSequence edgeSequence;
        int edgeCoverage;
        tie(edgeSequence, edgeCoverage) = getDominantEdgeSequence(e);

        if(edgeSequence.sequence.empty()) {
            const int overlap = edgeSequence.overlap;
            position += (k - overlap);
            // cout << "overlap " << overlap << endl;

        } else {

            for(const Base base: edgeSequence.sequence) {
                sequence.push_back(make_pair(base, edgeCoverage));
            }
            position = sequence.size();
        }


    }

    return sequence;

}



// Write the graph in Graphviz format.
void LocalMarkerGraph::write(const string& fileName, bool detailed) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream, detailed);
}
void LocalMarkerGraph::write(ostream& s, bool detailed) const
{
    Writer writer(*this, detailed);
    boost::write_graphviz(s, *this, writer, writer, writer,
        boost::get(&LocalMarkerGraphVertex::vertexId, *this));
}

LocalMarkerGraph::Writer::Writer(const LocalMarkerGraph& graph, bool detailed) :
    graph(graph),
    detailed(detailed)
{
}



void LocalMarkerGraph::Writer::operator()(std::ostream& s) const
{
    if(detailed) {
        s << "layout=dot;\n";
        s << "ratio=expand;\n";
        s << "node [fontname = \"Courier New\" shape=rectangle];\n";
        s << "edge [fontname = \"Courier New\"];\n";
    } else {
        s << "layout=sfdp;\n";
        s << "ratio=expand;\n";
        s << "node [shape=point];\n";
    }
}


void LocalMarkerGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
    const LocalMarkerGraphVertex& vertex = graph[v];
    const auto coverage = vertex.markerIds.size();



    // For compact output, the node shape is already defaulted to point,
    // and we don't write a label. The tooltip contains the vertex id,
    // which can be used to create a local subgraph to be looked at
    // in detailed format (use scripts/CreateLocalSubgraph.py).
    if(!detailed) {

        // Begin vertex attributes.
        s << "[";

        // Tooltip.
        s << "tooltip=\"" << coverage << " " << vertex.vertexId << "\"";

        string color;
        if(coverage >= graph.minCoverage) {
            color = "black";
        } else if(coverage == 1) {
            color = "#ff000080";  // Red, half way transparent
        } else if(coverage == 2) {
            color = "#ff800080";  // Orange, half way transparent
        } else if(coverage < graph.minCoverage) {
            color = "#ff80ff80";  // Purple, half way transparent
        }
        if(!color.empty()) {
            s << " fillcolor=\"" << color << "\" color=\"" << color << "\"";
        }

        // Vertex size.
        s << " width=\"";
        const auto oldPrecision = s.precision(4);
        s << 0.05 * sqrt(double(coverage));
        s.precision(oldPrecision);
        s << "\"";

        // End vertex attributes.
        s << "]";
        return;
    }



    // If getting here, we are doing detailed output.
    const Kmer kmer(vertex.kmerId, graph.k);

    // Begin vertex attributes.
    s << "[";


    // Color.
    string color;
    if(coverage >= graph.minCoverage) {
        color = "green";
    } else if(coverage == 1) {
        color = "#ff0000";  // Red
    } else if(coverage == 2) {
        color = "#ff8000";  // Orange
    } else {
        color = "#ff80ff";  // Purple
    }
    s << " style=filled";
    s << " fillcolor=\"" << color << "\"";


    // Write the label using Graphviz html-like functionality.
    s << " label=<<font><table border=\"0\">";

    // Row with the vertex id.
    s << "<tr><td colspan=\"3\"><b>";
    s << vertex.vertexId;
    s << "</b></td></tr>";

    // Row with the k-mer sequence.
    s << "<tr><td colspan=\"3\"><b>";
    kmer.write(s, graph.k);
    s << "</b></td></tr>";

    // Row with the k-mer id.
    s << "<tr><td colspan=\"3\"><b>";
    s << vertex.kmerId;
    s << "</b></td></tr>";

    // Row with column headers.
    s << "<tr><td><b>Id</b></td><td><b>Ord</b></td><td><b>Pos</b></td></tr>";

    // A row for each k-mer occurrence.
    for(const auto& markerId: graph[v].markerIds) {
        // OrientedReadId
        s << "<tr><td align=\"right\"><b>" << graph.orientedReadIds[markerId.localOrientedReadId] << "</b></td>";
        // Ordinal.
        s << "<td align=\"right\"><b>" << markerId.ordinal << "</b></td>";
        // Position.
        const auto& markers0 = graph.markers0[markerId.localOrientedReadId];
        const int position = markers0[markerId.ordinal].position;
        s << "<td align=\"right\"><b>" << position << "</b></td></tr>";
    }
    s << "</table></font>>";

    // End vertex attributes.
    s << "]";
}



void LocalMarkerGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
    const LocalMarkerGraphEdge& edge = graph[e];

    // Find the maximum number of oriented reads that agree.
    const size_t consensus = graph.edgeConsensus(e);


    // Compact output.
    if(!detailed) {
        // Begin edge attributes.
        s << "[";

        s << "tooltip=\"" << consensus << "/" << edge.data.size() << "\"";

        s << " penwidth=";
        const auto oldPrecision = s.precision(4);
        s <<  0.5 * double(consensus);
        s.precision(oldPrecision);

        // Color.
        string color;
        if(consensus >= graph.minConsensus) {
            color = "black";
        } else if(consensus == 1) {
            color = "#ff000080";  // Red, half way transparent
        } else if(consensus == 2) {
            color = "#ff800080";  // Orange, half way transparent
        } else {
            color = "#ff80ff80";  // Purple, half way transparent
        }
        s << " color=\"" << color << "\"";

        // End edge attributes.
        s << "]";
        return;
    }



    // If getting here, we are doing detailed output.

    // Get the sequence strings.
    std::map<string, vector<LocalMarkerGraphEdge::Data> > sequenceMap;
    for(const auto& data: edge.data) {
        sequenceMap[graph.getSequence(data)].push_back(data);
    }


    // Begin edge attributes.
    s << "[";

    // Tooltip.
    s << "tooltip=\"" << consensus << "/" << edge.data.size() << "\"";
    s << " penwidth=" << consensus;

    // Color.
    string color;
    string fillColor;
    if(consensus >= graph.minConsensus) {
        color = "black";
        fillColor = "green";
    } else if(consensus == 1) {
        color = "red";
        fillColor = color;
    } else if(consensus == 2) {
        color = "#ff8000";  // Orange
        fillColor = color;
    } else {
        color = "#ff80ff";  // Purple
        fillColor = color;
    }
    s << " fillcolor=\"" << fillColor << "\"";
    s << " color=\"" << color << "\"";


    // Label.
    s << " label=<<font color=\"black\">";
    s << "<table";
    s << " color=\"black\"";
    s << " bgcolor=\"" << fillColor << "\"";
    s << " border=\"0\"";
    s << " cellborder=\"1\"";
    s << " cellspacing=\"0\"";
    s << ">";

    for(const auto& p: sequenceMap) {
        const string& sequenceString = p.first;
        const auto& data = p.second;
        for(auto it=data.begin(); it!=data.end(); ++it) {
            s << "<tr><td align=\"right\"><b>" << graph.orientedReadIds[it->localOrientedReadId] << "</b></td>";
            s << "<td align=\"right\"><b>" << it->shift() << "</b></td>";
            if(it!=data.begin()) {
                s << "<td><b>";
                s << "=";
            } else {
                s << "<td align=\"left\"><b>";
                s << sequenceString;
            }
            s << "</b></td></tr>";
        }

    }
    s << "</table></font>> decorate=true";


    // End edge attributes.
    s << "]";
}


// Check that the graph and the vertex map are consistent.
void LocalMarkerGraph::check()
{
    const LocalMarkerGraph& graph = *this;

    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        for(const auto& markerId: graph[v].markerIds) {
            CZI_ASSERT(vertexMap[markerId.localOrientedReadId][markerId.ordinal] == v);
        }
    }
}
