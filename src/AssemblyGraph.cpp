#include "AssemblyGraph.hpp"
#include "deduplicate.hpp"
using namespace shasta;

#include "fstream.hpp"
#include "iterator.hpp"



void AssemblyGraph::createMarkerToAssemblyTable(uint64_t markerGrapEdgeCount)
{
    markerToAssemblyTable.beginPass1(markerGrapEdgeCount);
    for(EdgeId assemblyGraphEdgeId=0; assemblyGraphEdgeId<edgeLists.size(); assemblyGraphEdgeId++) {
        const span<EdgeId> chain = edgeLists[assemblyGraphEdgeId];
        for(uint32_t position=0; position!=chain.size(); position++) {
            const EdgeId markerGraphEdgeId = chain[position];
            markerToAssemblyTable.incrementCount(markerGraphEdgeId);
        }
    }
    markerToAssemblyTable.beginPass2();
    for(EdgeId assemblyGraphEdgeId=0; assemblyGraphEdgeId<edgeLists.size(); assemblyGraphEdgeId++) {
        const span<EdgeId> chain = edgeLists[assemblyGraphEdgeId];
        for(uint32_t position=0; position!=chain.size(); position++) {
            const EdgeId markerGraphEdgeId = chain[position];
            markerToAssemblyTable.store(
                markerGraphEdgeId, make_pair(assemblyGraphEdgeId, position));
        }
    }
    markerToAssemblyTable.endPass2();

}



// Close all open data.
void AssemblyGraph::close()
{
    if(vertices.isOpen) {
        vertices.close();
    }

    if(reverseComplementVertex.isOpen) {
        reverseComplementVertex.close();
    }

    if(edges.isOpen) {
        edges.close();
    }

    if(reverseComplementEdge.isOpen) {
        reverseComplementEdge.close();
    }

    if(edgesBySource.isOpen()) {
        edgesBySource.close();
    }

    if(edgesByTarget.isOpen()) {
        edgesByTarget.close();
    }

    if(edgeLists.isOpen()) {
        edgeLists.close();
    }

    if(bubbles.isOpen) {
        bubbles.close();
    }

    if(markerToAssemblyTable.isOpen()) {
        markerToAssemblyTable.close();
    }

    if(sequences.isOpen()) {
        sequences.close();
    }

    if(repeatCounts.isOpen()) {
        repeatCounts.close();
    }

    if(orientedReadsByEdge.isOpen()) {
        orientedReadsByEdge.close();
    }
}


// Close and remove all open data.
void AssemblyGraph::remove()
{
    if(vertices.isOpen) {
        vertices.remove();
    }

    if(reverseComplementVertex.isOpen) {
    	reverseComplementVertex.remove();
    }

    if(edges.isOpen) {
        edges.remove();
    }

    if(reverseComplementEdge.isOpen) {
    	reverseComplementEdge.remove();
    }

    if(edgesBySource.isOpen()) {
        edgesBySource.remove();
    }

    if(edgesByTarget.isOpen()) {
        edgesByTarget.remove();
    }

    if(edgeLists.isOpen()) {
        edgeLists.remove();
    }

    if(bubbles.isOpen) {
        bubbles.remove();
    }

    if(markerToAssemblyTable.isOpen()) {
        markerToAssemblyTable.remove();
    }

    if(sequences.isOpen()) {
        sequences.remove();
    }

    if(repeatCounts.isOpen()) {
        repeatCounts.remove();
    }

    if(orientedReadsByEdge.isOpen()) {
        orientedReadsByEdge.remove();
    }
}


// Basic Graphviz output of the global assembly graph.
void AssemblyGraph::writeGraphviz(const string& fileName) const
{
    ofstream graphOut(fileName);
    graphOut << "digraph AssemblyGraph {\n";

    // Write the vertices.
    // The label contains the corresponding marker graph vertex id.
    for(VertexId vertexId=0; vertexId<vertices.size(); vertexId++) {
        graphOut <<
            vertexId <<
            " [label=\"" <<
            vertexId << "\\n" << vertices[vertexId] <<
             "\"];\n";
    }

    // Write the edges.
    // The label contains the edge id and the number of maker graph edges
    // that correspond to this assembly graph edge.
    for(EdgeId edgeId=0; edgeId<edges.size(); edgeId++) {
        const Edge& edge = edges[edgeId];
        graphOut <<
            edge.source << "->" << edge.target <<
            " [label=\"" << edgeId << "\\n" <<
            edgeLists.size(edgeId) <<
            "\"];\n";
    }

    graphOut << "}\n";
}



// Find the out-degree or in-degree of a vertex.
// This is not simply the same as counting edgesBySource
// and edgesByTarget, because we have to skip edges
// that were removed.
AssemblyGraph::VertexId AssemblyGraph::inDegree(VertexId vertexId) const
{
    const auto e = edgesByTarget[vertexId];
    VertexId inDegree = 0;
    for(const EdgeId edgeId: e) {
        if(!edges[edgeId].wasRemoved()) {
            ++inDegree;
        }
    }
    return inDegree;
}
AssemblyGraph::VertexId AssemblyGraph::outDegree(VertexId vertexId) const
{
    const auto e = edgesBySource[vertexId];
    VertexId outDegree = 0;
    for(const EdgeId edgeId: e) {
        if(!edges[edgeId].wasRemoved()) {
            ++outDegree;
        }
    }
    return outDegree;
}



// Fill in edgesBySource and edgesByTarget.
void AssemblyGraph::computeConnectivity()
{
    edgesBySource.beginPass1(vertices.size());
    edgesByTarget.beginPass1(vertices.size());
    for(const Edge& edge: edges) {
        edgesBySource.incrementCount(edge.source);
        edgesByTarget.incrementCount(edge.target);
    }
    edgesBySource.beginPass2();
    edgesByTarget.beginPass2();
    for(EdgeId edgeId=0; edgeId<edges.size(); edgeId++) {
        const Edge& edge = edges[edgeId];
        edgesBySource.store(edge.source, edgeId);
        edgesByTarget.store(edge.target, edgeId);
    }
    edgesBySource.endPass2();
    edgesByTarget.endPass2();

    // Make sure edges by source and by target are sorted.
    for(VertexId vertexId=0; vertexId<vertices.size(); vertexId++) {
        const auto es = edgesBySource[vertexId];
        const auto et = edgesByTarget[vertexId];
        sort(es.begin(), es.end());
        sort(et.begin(), et.end());
    }

}



// Find incoming/outgoing edges of a vertex
// that were not removed.
// They are returned sorted by edge id.
void AssemblyGraph::findInEdges(VertexId vertexId, vector<EdgeId>& edgeIds) const
{
    const auto e = edgesByTarget[vertexId];
    edgeIds.clear();
    for(const EdgeId edgeId: e) {
        if(!edges[edgeId].wasRemoved()) {
            edgeIds.push_back(edgeId);
        }
    }
}
void AssemblyGraph::findOutEdges(VertexId vertexId, vector<EdgeId>& edgeIds) const
{
    const auto e = edgesBySource[vertexId];
    edgeIds.clear();
    for(const EdgeId edgeId: e) {
        if(!edges[edgeId].wasRemoved()) {
            edgeIds.push_back(edgeId);
        }
    }
}


// Find bubbles in the assembly graph.
// A bubble is a set of two vertices v0, v1,
// such that the outgoing edges of v0 are the same
// as the outgoing edges of v1, and
// out-degree(v0) = in-degree(v1) > 1.
// v0 is called the bubble source.
// v1 is called the bubble target.
// In defining and detecting bubbles, edges
// that were removed are considered to not exist.
// This assumes that the bubble Vector was already initialized.
void AssemblyGraph::findBubbles()
{
    // Define vectors here to reduce memory allocation overhead.
    vector<EdgeId> v0OutEdges;
    vector<EdgeId> v1InEdges;

    // Histogram of number of bubbles by ploidy.
    vector<EdgeId> histogram;

    // Loop over possible source vertices for a bubble.
    bubbles.clear();
    for(VertexId v0=0; v0<vertices.size(); v0++) {

        // Find the out-edges of v0.
        findOutEdges(v0, v0OutEdges);

        // If less than 2 out-edges, v0 is not a bubble source.
        if(v0OutEdges.size() < 2) {
            continue;
        }

        // Find the ossible bubble target.
        const VertexId v1 = edges[v0OutEdges.front()].target;

        // Check if all out-edges of v0 have the same target.
        bool isBubble = true;
        for(const EdgeId e: v0OutEdges) {
            if(edges[e].target != v1) {
                isBubble = false;
                break;
            }
        }
        if(!isBubble) {
            continue;
        }

        // Find the in-edges of v1.
        findInEdges(v1, v1InEdges);

        // For this to be a bubble, the in-edges of v1
        // must be the same as the out-edges of v0.
        if(v1InEdges != v0OutEdges) {
            continue;
        }

        // Store this bubble.
        bubbles.push_back(Bubble(v0, v1));

        // Update the histogram by ploidy.
        const uint64_t ploidy = v0OutEdges.size();
        if(histogram.size() <= ploidy) {
            histogram.resize(ploidy+1, 0);
        }
        ++histogram[ploidy];

    }
    VertexId bubbleCount = 0;
    VertexId bubbleEdgeCount = 0;
    for(uint64_t ploidy=0; ploidy<histogram.size(); ploidy++) {
        const uint64_t frequency = histogram[ploidy];
        if(frequency) {
            bubbleCount += frequency;
            bubbleEdgeCount += frequency * ploidy;
            cout << "Found " << frequency << " bubbles with ploidy " << ploidy << endl;
        }
    }
    cout << "Total number of bubbles is " << bubbleCount << "." << endl;

    // Edge statistics.
    const AssemblyGraph::EdgeId totalEdgeCount = edgeCount();
    cout << "Assembly graph edge statistics:" << endl;
    cout << "Total number of edges is " << totalEdgeCount << "." << endl;
    cout << "Total number of bubble edges is " << bubbleEdgeCount << "." << endl;
    cout << "Total number of non-bubble edges is " << totalEdgeCount-bubbleEdgeCount << "." << endl;
}



// Return the number of edges that were not removed.
AssemblyGraph::EdgeId AssemblyGraph::edgeCount() const
{
    EdgeId count = 0;
    for(const Edge& edge: edges) {
        if(!edge.wasRemoved()) {
            ++count;
        }
    }
    return count;
}



// Find and store the assembly graph forks.
// The fork corresponding to a bubble is stored only once.
// This function assumes that the assembly graph has no removed edges.
void AssemblyGraph::createForks()
{
    // Simple-minded code to find forks.
    vector<ForkInfo> forkInfos;
    for(VertexId vertexId=0; vertexId<vertices.size(); vertexId++) {
        const span<EdgeId> outgoingEdges = edgesBySource[vertexId];
        if(outgoingEdges.size() > 1) {
            ForkInfo forkInfo;
            forkInfo.vertexId = vertexId;
            forkInfo.isForward = true;
            copy(outgoingEdges.begin(), outgoingEdges.end(), back_inserter(forkInfo.edgeIds));
            forkInfos.push_back(forkInfo);
        }
        const span<EdgeId> incomingEdges = edgesByTarget[vertexId];
        if(incomingEdges.size() > 1) {
            ForkInfo forkInfo;
            forkInfo.vertexId = vertexId;
            forkInfo.isForward = false;
            copy(incomingEdges.begin(), incomingEdges.end(), back_inserter(forkInfo.edgeIds));
            forkInfos.push_back(forkInfo);
        }
    }

    // Remove duplicates (the fork in a bubble is stored only once).
    const uint64_t forkCountBeforeDeduplication = forkInfos.size();
    deduplicate(forkInfos);

    // Store.
    forks.clear();
    for(const Fork& fork: forkInfos) {
        forks.push_back(fork);
    }
    const uint64_t forkCount = forks.size();
    const uint64_t bubbleCount = forkCountBeforeDeduplication - forkCount;
    cout << "Found " << forks.size() << " forks in the assembly graph, including " <<
        bubbleCount <<
        " bubbles." << endl;

    // Compute the size distribution of forks.
    vector<uint64_t> histogram;
    for(const Fork& fork: forks) {
        const uint64_t branchCount =
            fork.isForward ?
            edgesBySource.size(fork.vertexId) :
            edgesByTarget.size(fork.vertexId);
        if(branchCount >= histogram.size()) {
            histogram.resize(branchCount+1, 0);
        }
        ++histogram[branchCount];
    }
    cout << "Size distribution of forks:" << endl;
    for(uint64_t branchCount=0; branchCount<histogram.size(); branchCount++) {
        const uint64_t frequency = histogram[branchCount];
        if(frequency) {
            cout << branchCount << " " << frequency << endl;
        }
    }

    // Write details of each fork.
    for(size_t iFork=0; iFork<forks.size(); iFork++) {
        const Fork& fork = forks[iFork];
        const auto edges =
            fork.isForward ?
            edgesBySource[fork.vertexId] :
            edgesByTarget[fork.vertexId];
        cout << "Fork " << iFork << ":";
        for(const auto edge: edges) {
            cout << " " << edge;
        }
        cout << endl;
    }
}



span<AssemblyGraph::EdgeId> AssemblyGraph::getForkEdges(uint64_t forkId)
{
    const Fork& fork = forks[forkId];
    if(fork.isForward) {
        return edgesBySource[fork.vertexId];
    } else {
        return edgesByTarget[fork.vertexId];
    }
}
