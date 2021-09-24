// Shasta.
#include "AssemblyPathGraph2.hpp"
#include "deduplicate.hpp"
#include "html.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include <set>



AssemblyPathGraph2::AssemblyPathGraph2(
    const AssemblyGraph& assemblyGraph,
    uint64_t diagonalReadCountMin,
    uint64_t offDiagonalReadCountMax,
    double detangleOffDiagonalRatio) :
    diagonalReadCountMin(diagonalReadCountMin),
    offDiagonalReadCountMax(offDiagonalReadCountMax),
    detangleOffDiagonalRatio(detangleOffDiagonalRatio)
{
    AssemblyPathGraph2& graph = *this;

    // Create a vertex for each assembly graph vertex.
    vector<vertex_descriptor> vertexDescriptors;
    for(AssemblyGraph::VertexId vertexId=0; vertexId<assemblyGraph.vertices.size(); vertexId++) {
        const vertex_descriptor v = add_vertex(AssemblyPathGraph2Vertex(vertexId), graph);
        vertexDescriptors.push_back(v);
    }

    // Fill in the reverse complemented vertex.
    for(AssemblyGraph::VertexId vertexId=0; vertexId<assemblyGraph.vertices.size(); vertexId++) {
        vertex_descriptor v = vertexDescriptors[vertexId];
        graph[v].reverseComplementVertex = vertexDescriptors[assemblyGraph.reverseComplementVertex[vertexId]];
    }

    // Sanity check.
    BGL_FORALL_VERTICES(v, graph, AssemblyPathGraph2) {
        SHASTA_ASSERT(graph[graph[v].reverseComplementVertex].reverseComplementVertex == v);
    }



    // Create an edge for each assembly graph edge.
    vector<edge_descriptor> edgeDescriptors;
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[edgeId];
        const AssemblyGraph::VertexId vertexId0 = edge.source;
        const AssemblyGraph::VertexId vertexId1 = edge.target;
        const vertex_descriptor v0 = vertexDescriptors[vertexId0];
        const vertex_descriptor v1 = vertexDescriptors[vertexId1];
        edge_descriptor e;
        tie(e, ignore) = add_edge(v0, v1, AssemblyPathGraph2Edge(edgeId), graph);
        edgeDescriptors.push_back(e);
    }

    // Fill in the reverse complemented edge.
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
        edge_descriptor e = edgeDescriptors[edgeId];
        graph[e].reverseComplementEdge = edgeDescriptors[assemblyGraph.reverseComplementEdge[edgeId]];
    }

    // Sanity check.
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph2) {
        SHASTA_ASSERT(graph[graph[e].reverseComplementEdge].reverseComplementEdge == e);
    }
}



void AssemblyPathGraph2::writeGraphviz(const string& fileName) const
{
    ofstream file(fileName);
    writeGraphviz(file);
}



void AssemblyPathGraph2::writeGraphviz(ostream& s) const
{
    const AssemblyPathGraph2& graph = *this;

    s << "digraph G {\n";

    // Default attributes.
    // The font is necessary to avoid labels that go out
    // of their shape.
    s << "K=10;\n";
    s << "overlap=false;\n";
    s << "splines=true;\n";
    s << "smoothing=triangle;\n";
    s << "levels=10;\n";
    s << "node [shape=point fontname=\"Courier New\"];\n";

    // This turns off the tooltip on the graph and the edges.
    s << "tooltip = \" \";\n";
    s << "edge[tooltip = \" \"];\n";



    // Vertices.
    BGL_FORALL_VERTICES(v, graph, AssemblyPathGraph2) {
        const AssemblyPathGraph2Vertex& vertex = graph[v];
        s << vertex.vertexId;

        s << " [";

        s << "tooltip=\"" << vertex.vertexId << "\"";

        if(in_degree(v, graph)==0 or out_degree(v,graph)==0) {
            s << " color=red";
        }

        s << "]";

        s << "\n";
    }



    // Edges. We write each edge as an additional pseudovertex.
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph2) {
        const AssemblyPathGraph2Edge& edge = graph[e];

        // Get the vertices.
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const AssemblyPathGraph2Vertex& vertex0 = graph[v0];
        const AssemblyPathGraph2Vertex& vertex1 = graph[v1];


        // Write is as a pseudo vertex.
        const string pseudoVertexName = "edge" + string(edge);
        s << "\"" << pseudoVertexName << "\" [";
        s << "shape=rectangle";

        // Label.
        s << " label=\"" << edge << "\\n" << edge.pathLength << "";
        if(edge.tangle != invalidTangle2Id) {
            s << "\\n" << edge.tangle;
        }
        s << "\"";

        // Tooltip.
        s << " tooltip=\"Path " << edge << ", " << edge.pathLength << " markers";
        if(edge.tangle != invalidTangle2Id) {
            s << ", tangle " << edge.tangle;
        }
        s << "\"";



        // Color.
        if(edge.tangle != invalidTangle2Id) {
            // Tangle2 edge.
            SHASTA_ASSERT(edge.inTangle == invalidTangle2Id);
            SHASTA_ASSERT(edge.outTangle == invalidTangle2Id);
            const Tangle2& tangle = getTangle(edge.tangle);
            if(tangle.isSolvable) {
                s << " style=filled fillcolor=pink";
            } else {
                s << " style=filled fillcolor=orange";
            }
        } else if(edge.inTangle != invalidTangle2Id and edge.outTangle != invalidTangle2Id) {
            // The edge is an in-edge of a tangle and an out-edge of another tangle.
            s << " style=filled fillcolor=purple";
        } else if(edge.inTangle != invalidTangle2Id) {
            // The edge has an in-tangle, so it is an out-edge of a tangle.
            s << " style=filled fillcolor=red";
        } else if(edge.outTangle != invalidTangle2Id) {
            // The edge has an out-tangle, so it is an in-edge of a tangle.
            s << " style=filled fillcolor=green";
        }
        s << "];\n";



        // Write the arrows to/from the pseudovertex.
        s << vertex0.vertexId << "->\"" << pseudoVertexName << "\";\n";
        s << "\"" << pseudoVertexName << "\"->" << vertex1.vertexId << ";\n";
    }





    s << "}\n";
}



// Initial creation of all tangles.
void AssemblyPathGraph2::createTangles()
{
    AssemblyPathGraph2& graph = *this;

    // Just in case, clean up.
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph2) {
        graph[e].clearTangles();
    }
    tangles.clear();
    nextTangleId = 0;

    // Create the tangles.
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph2) {
        createTangleAtEdge(e);
    }
    cout << "Found " << tangles.size() << " tangles." << endl;
}



// Create a new tangle that has the specified edge
// as the tangle edge, if such a tangle is valid
// and does not already exist.
// Return true if the new tangle was created.
bool AssemblyPathGraph2::createTangleAtEdge(edge_descriptor e01)
{
    AssemblyPathGraph2& graph = *this;

    if(graph[e01].tangle != invalidTangle2Id) {
        return false;
    }

    const vertex_descriptor v0 = source(e01, graph);
    const vertex_descriptor v1 = target(e01, graph);

    // If this edge does not generate a tangle, return false.
    // See the top of AssemblyPathGraph2.hpp for details.
    if(out_degree(v0, graph) != 1) {
        return false;
    }
    if(in_degree(v1, graph) != 1) {
        return false;
    }
    if(in_degree(v0, graph) < 2) {
        return false;
    }
    if(out_degree(v1, graph) < 2) {
        return false;
    }

    // If there is an edge v1->v0, this is a reverse bubble, not a tangle.
    bool reverseEdgeExists = false;
    tie(ignore, reverseEdgeExists) = boost::edge(v1, v0, graph);
    if(reverseEdgeExists) {
        return false;
    }

    const auto inDegree = in_degree(v0, graph);
    const auto outDegree = out_degree(v1, graph);

    Tangle2 tangle;
    tangle.edge = e01;
    SHASTA_ASSERT(graph[e01].tangle == invalidTangle2Id);
    graph[e01].tangle = nextTangleId;

    // Gather the in-edges and out-edges.
    BGL_FORALL_INEDGES(v0, e, graph, AssemblyPathGraph2) {
        tangle.inEdges.push_back(e);
        SHASTA_ASSERT(graph[e].outTangle == invalidTangle2Id);
        graph[e].outTangle = nextTangleId;
    }
    BGL_FORALL_OUTEDGES(v1, e, graph, AssemblyPathGraph2) {
        tangle.outEdges.push_back(e);
        SHASTA_ASSERT(graph[e].inTangle == invalidTangle2Id);
        graph[e].inTangle = nextTangleId;
    }



    // Compute the tangle matrix, which contains the number of common oriented reads
    // for each pair of in-edges and out-edges.
    vector<OrientedReadId> commonOrientedReadIds;
    tangle.matrix.resize(inDegree, vector<uint64_t>(outDegree));
    for(uint64_t inEdgeIndex=0; inEdgeIndex<inDegree; inEdgeIndex++) {
        const AssemblyPathGraph2Edge& inEdge = graph[tangle.inEdges[inEdgeIndex]];
        for(uint64_t outEdgeIndex=0; outEdgeIndex<outDegree; outEdgeIndex++) {
            const AssemblyPathGraph2Edge& outEdge = graph[tangle.outEdges[outEdgeIndex]];
            commonOrientedReadIds.clear();
            std::set_intersection(
                inEdge.orientedReadIds.begin(), inEdge.orientedReadIds.end(),
                outEdge.orientedReadIds.begin(), outEdge.orientedReadIds.end(),
                back_inserter(commonOrientedReadIds));
            tangle.matrix[inEdgeIndex][outEdgeIndex] = commonOrientedReadIds.size();
        }
    }

    // Now that the tangle matrix is available we can find out if this tangle
    // is solvable, and if it is we can compute its priority.
    tangle.findIfSolvable(
        diagonalReadCountMin,
        offDiagonalReadCountMax,
        detangleOffDiagonalRatio);
    tangle.computePriority();

    tangle.tangleId = nextTangleId;
    tangles.insert(make_pair(nextTangleId++, tangle));
    cout << "Created tangle " << tangle.tangleId << " at " << graph[e01] << endl;

    return true;
}



// Create tangles involving a given edge.
// This can create up to two tangles involving
// the given edge as an in-edge, out-edge, or tangle edge.
// This is used for incrementally create new tangles as
// edges are created during detangling.
void AssemblyPathGraph2::createTanglesInvolvingEdge(edge_descriptor e)
{
    AssemblyPathGraph2& graph = *this;
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);

    createTangleAtEdge(e);

    BGL_FORALL_INEDGES(v0, e, graph, AssemblyPathGraph2) {
        createTangleAtEdge(e);
    }
    BGL_FORALL_OUTEDGES(v1, e, graph, AssemblyPathGraph2) {
        createTangleAtEdge(e);
    }
}



Tangle2& AssemblyPathGraph2::getTangle(Tangle2Id tangleId)
{
    auto it = tangles.find(tangleId);
    SHASTA_ASSERT(it != tangles.end());
    Tangle2& tangle = it->second;
    SHASTA_ASSERT(tangle.tangleId == tangleId);
    return tangle;
}

// Const version.
const Tangle2& AssemblyPathGraph2::getTangle(Tangle2Id tangleId) const
{
    auto it = tangles.find(tangleId);
    SHASTA_ASSERT(it != tangles.end());
    const Tangle2& tangle = it->second;
    SHASTA_ASSERT(tangle.tangleId == tangleId);
    return tangle;
}



Tangle2Id AssemblyPathGraph2::getReverseComplementTangle(
    Tangle2Id tangleId) const
{
    const AssemblyPathGraph2& graph = *this;

    // Get the edge of this tangle.
    const edge_descriptor e = getTangle(tangleId).edge;
    const AssemblyPathGraph2Edge& edge = graph[e];

    // Get the reverse complement edge.
    const edge_descriptor eReverseComplement = edge.reverseComplementEdge;
    const AssemblyPathGraph2Edge& reverseComplementEdge = graph[eReverseComplement];

    // Return its tangle.
    const Tangle2Id reverseComplementTangleId = reverseComplementEdge.tangle;
    SHASTA_ASSERT(reverseComplementTangleId != invalidTangle2Id);
    return reverseComplementTangleId;
}



// Detangle all we can.
// The average number of bases per marker is only used
// for GFA output.
void AssemblyPathGraph2::detangle(
    double basesPerMarker,
    const AssemblyGraph& assemblyGraph)
{
    AssemblyPathGraph2& graph = *this;
    const bool debug = false;

    // Detangle iteration.
    for(int iteration=0; ; ++iteration) {

        const Tangle2Id tangleId = findNextTangle();
        if(tangleId == invalidTangle2Id) {
            break;
        }

        if(debug) {
            cout << "Detangle iteration " << iteration <<
                " begins, working on tangle " << tangleId <<
                " and its reverse complement tangle " <<
                getReverseComplementTangle(tangleId) << endl;

            // Write the graph at the beginning of this iteration.
            graph.writeGraphviz("AssemblyPathGraph2-" + to_string(iteration) + ".dot");
            graph.writeHtml("AssemblyPathGraph2-" + to_string(iteration) + ".html");
            graph.writeGfa("AssemblyPathGraph2-" + to_string(iteration) + ".gfa",
                basesPerMarker);
        }

        // Detangle this tangle and its reverse complement.
        vector<edge_descriptor> newEdges;
        detangleComplementaryPair(tangleId, newEdges);

        // Fill in the reverseComplementEdge for the edges we just created.
        fillReverseComplementNewEdges(newEdges, assemblyGraph);

        // Create tangles involving the newly created edges.
        for(const edge_descriptor e: newEdges) {
            createTanglesInvolvingEdge(e);
        }

        // Remove any vertices that were left isolated.
        removeIsolatedVertices();
    }

    // graph.writeGraphviz("AssemblyPathGraph2-Final.dot");
    // graph.writeHtml("AssemblyPathGraph2-Final.html");
    // graph.writeGfa("AssemblyPathGraph2-Final.gfa", basesPerMarker);
}


void AssemblyPathGraph2::fillReverseComplementNewEdges(
    const vector<edge_descriptor>& newEdges,
    const AssemblyGraph& assemblyGraph)
{
    AssemblyPathGraph2& graph = *this;

    // Loop over the new edges.
    for(const edge_descriptor e: newEdges)  {

        // Extract the path for this edge.
        const vector <AssemblyGraph::EdgeId>& path = graph[e].path;

        // Use the AssemblyGraph to create the reverse complement path.
        vector <AssemblyGraph::EdgeId> pathRc;
        for(const AssemblyGraph::EdgeId edgeId: path) {
            const AssemblyGraph::EdgeId edgeIdRc = assemblyGraph.reverseComplementEdge[edgeId];
            pathRc.push_back(edgeIdRc);
        }
        std::reverse(pathRc.begin(), pathRc.end());

        // Now look for this in the new edges.
        bool wasFound = false;
        for(const edge_descriptor eRc: newEdges) {
            if(graph[eRc].path == pathRc) {
                graph[e].reverseComplementEdge = eRc;
                wasFound = true;
                break;
            }
        }
        SHASTA_ASSERT(wasFound);
    }

    // Sanity check.
    for(const edge_descriptor e: newEdges)  {
        const  edge_descriptor eRc = graph[e].reverseComplementEdge;
        SHASTA_ASSERT(graph[eRc].reverseComplementEdge == e);
    }
}



// Detangle a single tangle.
// This does not fill in the reverseComplementEdge of newly created edges,
// and does not create new tangles involving those edges.
void AssemblyPathGraph2::detangle(
    Tangle2Id tangleId,
    vector<edge_descriptor>& newEdges)
{
    cout << "Detangling tangle " << tangleId << endl;
    AssemblyPathGraph2& graph = *this;
    Tangle2& tangle = getTangle(tangleId);
    SHASTA_ASSERT(tangle.isSolvable);



    // Create the new edges.
    // We loop over all pairs of matching in-edges and out-edges.
    const AssemblyPathGraph2Edge& tangleEdge = graph[tangle.edge];
    for(uint64_t i=0; i<tangle.inEdges.size(); i++) {
        const edge_descriptor eIn = tangle.inEdges[i];
        const AssemblyPathGraph2Edge& inEdge = graph[eIn];
        const vertex_descriptor vA = source(eIn, graph);
        const uint64_t j = tangle.match[i];

        const edge_descriptor eOut = tangle.outEdges[j];
        const AssemblyPathGraph2Edge& outEdge = graph[eOut];
        const vertex_descriptor vB = target(eOut, graph);

        // Add the new edge and fill in what we can now.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(vA, vB, graph);
        newEdges.push_back(eNew);
        AssemblyPathGraph2Edge& newEdge = graph[eNew];
        newEdge.pathLength =
            inEdge.pathLength + tangleEdge.pathLength + outEdge.pathLength;
        // Don't include the reads of the tangle edge in the new edge!
        newEdge.mergeOrientedReadIds(
            inEdge.orientedReadIds,
            outEdge.orientedReadIds
        );
        newEdge.path = inEdge.path;
        copy(tangleEdge.path.begin(), tangleEdge.path.end(),
            back_inserter(newEdge.path));
        copy(outEdge.path.begin(), outEdge.path.end(),
            back_inserter(newEdge.path));
    }


    // Remove other tangles that the in-edges and out-edges
    // of this tangle are involved in.
    // Those will be recreated later, using the combined edges.
    vector<Tangle2Id> tanglesToBeRemoved;
    for(const edge_descriptor e: tangle.inEdges) {
        const AssemblyPathGraph2Edge& edge = graph[e];
        SHASTA_ASSERT(edge.outTangle == tangleId);
        SHASTA_ASSERT(edge.tangle == invalidTangle2Id);
        if(edge.inTangle != invalidTangle2Id) {
            tanglesToBeRemoved.push_back(edge.inTangle);
            cout << "Will remove preceding tangle " << edge.inTangle <<
                " due to tangle in-edge " << edge << endl;
        }
    }
    for(const edge_descriptor e: tangle.outEdges) {
        const AssemblyPathGraph2Edge& edge = graph[e];
        SHASTA_ASSERT(edge.tangle == invalidTangle2Id);
        SHASTA_ASSERT(edge.inTangle == tangleId);
        if(edge.outTangle != invalidTangle2Id) {
            tanglesToBeRemoved.push_back(edge.outTangle);
            cout << "Will remove following tangle " << edge.outTangle <<
                " due to tangle out-edge " << edge << endl;
        }
    }
    deduplicate(tanglesToBeRemoved);
    cout << "Removing " << tanglesToBeRemoved.size() << " adjacent tangles." << endl;
    for(const Tangle2Id tangleId: tanglesToBeRemoved) {
        removeTangle(tangleId);
        cout << "Removed adjacent tangle " << tangleId << endl;
    }


    // Remove all the edges involved in the tangle we are detangling,
    // as well as the source and target vertices of the tangle edge.
    for(const edge_descriptor e: tangle.inEdges) {
        boost::remove_edge(e, graph);
    }
    for(const edge_descriptor e: tangle.outEdges) {
        boost::remove_edge(e, graph);
    }
    const vertex_descriptor v0 = source(tangle.edge, graph);
    const vertex_descriptor v1 = target(tangle.edge, graph);
    boost::remove_edge(tangle.edge, graph);
    SHASTA_ASSERT(in_degree(v0, graph) == 0);
    SHASTA_ASSERT(out_degree(v0, graph) == 0);
    SHASTA_ASSERT(in_degree(v1, graph) == 0);
    SHASTA_ASSERT(out_degree(v1, graph) == 0);
    remove_vertex(v0, graph);
    remove_vertex(v1, graph);

    // Finally we can remove this tangle.
    tangles.erase(tangleId);

}



// Detangle a tangle and its reverse complement.
// This does not fill in the reverseComplementEdge of newly created edges,
// and does not create new tangles involving those edges.
// If the tangles in the pair don't collide, they are detangled separately
// using the above detangle function.
// Otherwise, they are detangled atthe same time using
// detangleCollidingComplementaryPair.
void AssemblyPathGraph2::detangleComplementaryPair(
    Tangle2Id tangleId,
    vector<edge_descriptor>& newEdges)
{
    if(collidesWithReverseComplement(tangleId)) {
        detangleCollidingComplementaryPair(tangleId, newEdges);
    } else {

        // Get the id of the reverse complement tangle.
        // This needs to be done before callign detangle, otherwise
        // our tangle goes away!
        const Tangle2Id  reverseComplementTangleId = getReverseComplementTangle(tangleId);

        // Detangle both of them separately.
        detangle(tangleId, newEdges);
        detangle(reverseComplementTangleId, newEdges);
    }
}



// Detangle a tangle and its reverse complement
// that collide with each other (that is, share edges).
// This does not fill in the reverseComplementEdge of newly created edges,
// and does not create new tangles involving those edges.
void AssemblyPathGraph2::detangleCollidingComplementaryPair(
    Tangle2Id tangleIdA,
    vector<edge_descriptor>& newEdges)
{

    // For now, call A the tangle passed in as an argument and B
    // its reverse complement.
    Tangle2& tangleA = getTangle(tangleIdA);
    const Tangle2Id tangleIdB = getReverseComplementTangle(tangleIdA);
    Tangle2& tangleB = getTangle(tangleIdB);

    cout << "Detangling colliding pair of reverse complement tangles " <<
        tangleIdA << " " << tangleIdB << endl;
    SHASTA_ASSERT(tangleA.isSolvable);
    SHASTA_ASSERT(tangleB.isSolvable);
    AssemblyPathGraph2& graph = *this;

    // Gather in-edges and out-edges and sort them.
    auto inEdgesA  = tangleA.inEdges;
    auto inEdgesB  = tangleB.inEdges;
    auto outEdgesA = tangleA.outEdges;
    auto outEdgesB = tangleB.outEdges;
    sort(inEdgesA.begin(), inEdgesA.end());
    sort(inEdgesB.begin(), inEdgesB.end());
    sort(outEdgesA.begin(), outEdgesA.end());
    sort(outEdgesB.begin(), outEdgesB.end());

    // Figure out which of the two tangles follows the other.
    const bool BFollowsA = (inEdgesB == outEdgesA);
    const bool AFollowsB = (inEdgesA == outEdgesB);
    if(not(BFollowsA or AFollowsB)) {
        // At first sight it seems that this is not possible,
        // but it can actually happen in tangles with in-degree/out-degree
        // greater than 2.
        // Just mark both of them as unsolvable.
        tangleA.isSolvable = false;
        tangleB.isSolvable = false;
        tangleA.priority = 0;
        tangleB.priority = 0;
        cout << "Unusual arrangement of colliding pair of reverse complement tangles " <<
            tangleIdA << " " << tangleIdB <<
            "  was marked as unsolvable." << endl;
        return;
    }
    if(BFollowsA and AFollowsB) {
        // This is a horrible mess where the two tangles follow each other.
        // Just mark both of them as unsolvable.
        tangleA.isSolvable = false;
        tangleB.isSolvable = false;
        tangleA.priority = 0;
        tangleB.priority = 0;
        cout << "Colliding pair of reverse complement tangles " <<
            tangleIdA << " " << tangleIdB <<
            " follow each other and were marked as unsolvable." << endl;
        return;
    }

    // Arrange the two tangles in the pair as tangle0 and tangle1,
    // such that the out-edges of tangle0 are the same as the in-edges of tangle1
    // (possibly in different order).
    // That is, tangle1 follows tangle 0.
    Tangle2Id tangleId0 = tangleIdA;
    Tangle2Id tangleId1 = tangleIdB;
    if(AFollowsB) {
        std::swap(tangleId0, tangleId1);
    } else {
        SHASTA_ASSERT(BFollowsA);
    }
    cout << "Tangle2 " << tangleId1 << " follows tangle " << tangleId0 << endl;
    const Tangle2& tangle0 = getTangle(tangleId0);
    const Tangle2& tangle1 = getTangle(tangleId1);


    // At this point, we know that tangle 1 follows tangle0.
    // Use this nomenclature:
    // inEdges: the in-edges of tangle0;
    // middleEdges: the out-edges of tangle 0 and in-edges of tangle1
    //     (but the two tangles could store them in different orders).
    // outEdges: the out-edges of tangle1.


    // Create the new edges.
    // We loop over triplets of matching (inEdge, middleEdge, outEdge).
    const AssemblyPathGraph2Edge& tangleEdge0 = graph[tangle0.edge];
    const AssemblyPathGraph2Edge& tangleEdge1 = graph[tangle1.edge];
    for(uint64_t i=0; i<tangle0.inEdges.size(); i++) {

        // Find the inEdge.
        const edge_descriptor eIn = tangle0.inEdges[i];
        const AssemblyPathGraph2Edge& inEdge = graph[eIn];
        const vertex_descriptor vA = source(eIn, graph);

        // Find the matching middleEdge.
        const uint64_t j0 = tangle0.match[i];
        const edge_descriptor eMiddle = tangle0.outEdges[j0];
        const AssemblyPathGraph2Edge& middleEdge = graph[eMiddle];

        // Locate this middleEdge in the in-edges of the second tangle.
        const auto it = find(tangle1.inEdges.begin(), tangle1.inEdges.end(), eMiddle);
        SHASTA_ASSERT(it != tangle1.inEdges.end());
        const uint64_t j1 = it - tangle1.inEdges.begin();

        // Find the matching outEdge.
        const uint64_t k = tangle1.match[j1];
        const edge_descriptor eOut = tangle1.outEdges[k];
        const AssemblyPathGraph2Edge& outEdge = graph[eOut];
        const vertex_descriptor vB = target(eOut, graph);

        // Create a new edge by combining the following:
        // (inEdge, tangleEdge0, middleEdge, tangleEdge1, outEdge).
        // Add the new edge.
        edge_descriptor eNew;
        tie(eNew, ignore) = add_edge(vA, vB, graph);
        newEdges.push_back(eNew);
        AssemblyPathGraph2Edge& newEdge = graph[eNew];

        // Fill in the path length.
        newEdge.pathLength =
            inEdge.pathLength +
            tangleEdge0.pathLength +
            middleEdge.pathLength +
            tangleEdge1.pathLength +
            outEdge.pathLength;

        // Fill in the reads.
        // Don't include the reads of the tangle edges in the new edge!
        newEdge.mergeOrientedReadIds(
            inEdge.orientedReadIds,
            middleEdge.orientedReadIds,
            outEdge.orientedReadIds
        );

        // Fill in the path.
        newEdge.path = inEdge.path;
        copy(tangleEdge0.path.begin(), tangleEdge0.path.end(),
            back_inserter(newEdge.path));
        copy(middleEdge.path.begin(), middleEdge.path.end(),
            back_inserter(newEdge.path));
        copy(tangleEdge1.path.begin(), tangleEdge1.path.end(),
            back_inserter(newEdge.path));
        copy(outEdge.path.begin(), outEdge.path.end(),
            back_inserter(newEdge.path));
    }



    // Remove other tangles that the in-edges of tangle0 and out-edges
    // of tangle1 are involved in.
    // Those will be recreated later, using the combined edges.
    vector<Tangle2Id> tanglesToBeRemoved;
    for(const edge_descriptor e: tangle0.inEdges) {
        const AssemblyPathGraph2Edge& edge = graph[e];
        SHASTA_ASSERT(edge.outTangle == tangleId0);
        SHASTA_ASSERT(edge.tangle == invalidTangle2Id);
        if(edge.inTangle != invalidTangle2Id) {
            tanglesToBeRemoved.push_back(edge.inTangle);
            cout << "Will remove preceding tangle " << edge.inTangle <<
                " due to tangle in-edge " << edge << endl;
        }
    }
    for(const edge_descriptor e: tangle1.outEdges) {
        const AssemblyPathGraph2Edge& edge = graph[e];
        SHASTA_ASSERT(edge.tangle == invalidTangle2Id);
        SHASTA_ASSERT(edge.inTangle == tangleId1);
        if(edge.outTangle != invalidTangle2Id) {
            tanglesToBeRemoved.push_back(edge.outTangle);
            cout << "Will remove following tangle " << edge.outTangle <<
                " due to tangle out-edge " << edge << endl;
        }
    }
    deduplicate(tanglesToBeRemoved);
    cout << "Removing " << tanglesToBeRemoved.size() << " adjacent tangles." << endl;
    for(const Tangle2Id tangleId: tanglesToBeRemoved) {
        removeTangle(tangleId);
        cout << "Removed adjacent tangle " << tangleId << endl;
    }



    // Gather the vertices that should be removed, but don't remove
    // them for now. We have to remove the edges first.
    vector<vertex_descriptor> verticesToBeRemoved;
    verticesToBeRemoved.push_back(source(tangle0.edge, graph));
    verticesToBeRemoved.push_back(target(tangle0.edge, graph));
    verticesToBeRemoved.push_back(source(tangle1.edge, graph));
    verticesToBeRemoved.push_back(target(tangle1.edge, graph));



    // Remove all the edges involved in the two tangles we are detangling,
    // as well as the source and target vertices of the tangle edge.
    for(const edge_descriptor e: tangle0.inEdges) {
        boost::remove_edge(e, graph);
    }
    for(const edge_descriptor e: tangle0.outEdges) {
        boost::remove_edge(e, graph);
    }
    // Tangle1 in-edges are the same as tangle1 out-edges.
    for(const edge_descriptor e: tangle1.outEdges) {
        boost::remove_edge(e, graph);
    }
    boost::remove_edge(tangle0.edge, graph);
    boost::remove_edge(tangle1.edge, graph);



    // Now we can remove the vertices.
    for(vertex_descriptor v: verticesToBeRemoved) {
        SHASTA_ASSERT(in_degree(v, graph) == 0);
        SHASTA_ASSERT(out_degree(v, graph) == 0);
        clear_vertex(v, graph);
        remove_vertex(v, graph);
    }



    // Finally we can remove these two tangles.
    tangles.erase(tangleId0);
    tangles.erase(tangleId1);
}



void AssemblyPathGraph2::removeTangle(Tangle2Id tangleId)
{
    AssemblyPathGraph2& graph = *this;
    const Tangle2& tangle = getTangle(tangleId);

    // Remove all references to this tangle.
    graph[tangle.edge].tangle = invalidTangle2Id;

    for(const edge_descriptor e: tangle.inEdges) {
        AssemblyPathGraph2Edge& edge = graph[e];
        SHASTA_ASSERT(edge.tangle == invalidTangle2Id);
        SHASTA_ASSERT(edge.outTangle == tangleId);
        edge.outTangle = invalidTangle2Id;
    }

    for(const edge_descriptor e: tangle.outEdges) {
        AssemblyPathGraph2Edge& edge = graph[e];
        SHASTA_ASSERT(edge.tangle == invalidTangle2Id);
        SHASTA_ASSERT(edge.inTangle == tangleId);
        edge.inTangle = invalidTangle2Id;
    }

    // Now we can remove the tangle.
    tangles.erase(tangleId);
}



void AssemblyPathGraph2Edge::mergeOrientedReadIds(
    const vector<OrientedReadId>& r0,
    const vector<OrientedReadId>& r1
    )
{
    orientedReadIds.clear();
    std::set_union(
        r0.begin(), r0.end(),
        r1.begin(), r1.end(),
        back_inserter(orientedReadIds));
}



void AssemblyPathGraph2Edge::mergeOrientedReadIds(
    const vector<OrientedReadId>& r0,
    const vector<OrientedReadId>& r1,
    const vector<OrientedReadId>& r2
    )
{
    // r01 = union(r0, r1)
    vector<OrientedReadId> r01;
    std::set_union(
        r0.begin(), r0.end(),
        r1.begin(), r1.end(),
        back_inserter(r01));


    // orientedReadIds = union(r01, r2)
    orientedReadIds.clear();
    std::set_union(
        r01.begin(), r01.end(),
        r2.begin(), r2.end(),
        back_inserter(orientedReadIds));
}



bool Tangle2::hasZeroMatrixElements() const
{
    for(const auto& v: matrix) {
        for(const auto x: v) {
            if(x == 0) {
                return true;
            }
        }
    }

    // If getting here, we did not find any non-zero matrix elements.
    return false;
}



bool Tangle2::hasNonZeroMatrixElements() const
{
    for(const auto& v: matrix) {
        for(const auto x: v) {
            if(x != 0) {
                return true;
            }
        }
    }

    // If getting here, we did not find any non-zero matrix elements.
    return false;
}



uint64_t Tangle2::countNonZeroElementsInRow(uint64_t i) const
{
    uint64_t n = 0;
    for(uint64_t j=0; j<outDegree(); j++) {
        if(matrix[i][j] != 0) {
            ++n;
        }
    }
    return n;
}

uint64_t Tangle2::countNonZeroElementsInColumn(uint64_t j) const
{
    uint64_t n = 0;
    for(uint64_t i=0; i<inDegree(); i++) {
        if(matrix[i][j] != 0) {
            ++n;
        }
    }
    return n;
}



// Inspect the tangle matrix to find out if the
// tangle is solvable.
// If it is, also create the match vector,
// which describes how in-edges should be matched with out-edges.
void Tangle2::findIfSolvable(
    uint64_t diagonalReadCountMin,
    uint64_t offDiagonalReadCountMax,
    double detangleOffDiagonalRatio)
{
    // If the in-degree and out-degree are not the same,
    // mark it as not solvable.
    const uint64_t n = inDegree();
    if(outDegree() != n) {
        isSolvable = false;
        return;
    }


    // Inspect the tangle matrix, one in-edge at a time.
    for(uint64_t i=0; i<n; i++) {

        // Find the out-edge with the largest tangle matrix element
        // (number of supporting reads) for this in-edge.
        const vector<uint64_t>& v = matrix[i];
        const uint64_t j = std::max_element(v.begin(), v.end()) - v.begin();

        // Tentatively match this in-edge with this out-edge,
        // subject to additional checks below.
        match.push_back(j);
    }




    // Check that diagonal matrix elements defined by the
    // match vector are greater than all other elements in the same row and column.
    for(uint64_t i=0; i<n; i++) {
        const uint64_t j = match[i];
        bool ok = true;
        for(uint64_t ii=0; ii<n; ii++) {
            if(ii != i) {
                if(matrix[i][j] <= matrix[ii][j]) {
                    ok =false;
                }
            }
        }
        for(uint64_t jj=0; jj<n; jj++) {
            if(jj != j) {
                if(matrix[i][j] <= matrix[i][jj]) {
                    ok =false;
                }
            }
        }
        if(not ok) {
            isSolvable = false;
            match.clear();
            return;
        }
    }



    // Verify that the match vector is a permutation
    // (each out-edge must appear once and only once).
    SHASTA_ASSERT(match.size() == n);
    vector<uint64_t> matchCountCheck(n, 0);
    for(const uint64_t j: match) {
        ++matchCountCheck[j];
    }
    for(const uint64_t count: matchCountCheck) {
        if(count != 1) {
            isSolvable = false;
            match.clear();
            return;
        }
    }



    // Construct the inverse permutation.
    inverseMatch.resize(n);
    for(uint64_t i=0; i<n; i++) {
        const uint64_t j = match[i];
        inverseMatch[j] = i;
    }



    // Now check that the tangle matrix elements satisfy our requirements
    // We call element ij of the tangle matrix "diagonal"
    // if match[i]=j. Otherwise we call it "off-diagonal".
    // For the tangle to be marked solvable, we require the following:
    // - All diagonal elements must be >= diagonalReadCountMin.
    // - All off-diagonal elements must be satisfy ONE OF THE TWO
    //   following conditions:
    //   * They are <= offDiagonalReadCountMax.
    //   * Their ratios to the two corresponding diagonal elements
    //     are less than detangleOffDiagonalRatio.
    // In other words, we require off-diagonal matrix elements to
    // be "small" either in absolute terms of relative to the corresponding
    // diagonal elements.
    for(uint64_t i=0; i<n; i++) {
        for(uint64_t j=0; j<n; j++) {
            bool ok = true;
            if(j == match[i]) {

                // Diagonal element.
                if(matrix[i][j] < diagonalReadCountMin) {
                    ok = false;
                }

            } else {

                // Off-diagonal element.
                if(matrix[i][j] > offDiagonalReadCountMax) {
                    // This element did not pass the absolute criterion,
                    // so we have to check the relative criterion.
                    if(double(matrix[i][j]) / double(matrix[i][match[i]]) >
                        detangleOffDiagonalRatio) {
                        ok = false;
                    } else {
                        // It was ok in relative terms versus the diagonal element with the same i.
                        // Now we have to check against the diagonal element with the same j.
                        if(double(matrix[i][j]) / double(matrix[inverseMatch[j]][j]) >
                            detangleOffDiagonalRatio) {
                            ok = false;
                        }
                    }

                }

            }

            if(not ok) {
                isSolvable = false;
                match.clear();
                return;
            }
        }
    }





    // If we get here, all is good, the tangle is solvable,
    // and the match and inverseMatch vectors tell us how to match
    // out-edges with in-edges.
    isSolvable = true;
}



// The tangle priority is the lowest diagonal element of the tangle
// matrix. Solvable tangles are processed in order of decreasing priority.
void Tangle2::computePriority()
{
    if(isSolvable) {
        priority = std::numeric_limits<uint64_t>::max();
        for(uint64_t i=0; i<match.size(); i++) {
            const uint64_t j = match[i];
            priority = min(priority, matrix[i][j]);
        }
    } else {
        priority = 0;
    }
}



// Return the next tangle to work on.
// This does a linear search, which could be eliminated
// with appropriated data structures if it becomes a
// performance problem.

Tangle2Id AssemblyPathGraph2::findNextTangle() const
{

    Tangle2Id bestTangleId = invalidTangle2Id;
    uint64_t bestPriority = 0;
    for(const auto& p: tangles) {
        const Tangle2& tangle = p.second;
        if(not tangle.isSolvable) {
            continue;
        }
        if(tangle.priority > bestPriority) {
            bestPriority = tangle.priority;
            bestTangleId = tangle.tangleId;
        }
    }
    return bestTangleId;
}



// Return true if a tangle collides with its reverse complement.
bool AssemblyPathGraph2::collidesWithReverseComplement(Tangle2Id tangleId) const
{
    const AssemblyPathGraph2& graph = *this;
    const Tangle2& tangle = getTangle(tangleId);
    const Tangle2Id reverseComplementTangleId = getReverseComplementTangle(tangleId);

    // If the tangle is the same as its reverse complement, we have a collision.
    // This is unusual but possible.
    if(reverseComplementTangleId == tangleId) {
        return true;
    }

    // Check the in-edges.
    for(const edge_descriptor e: tangle.inEdges) {
        if(graph[e].inTangle == reverseComplementTangleId) {
            return true;
        }
    }

    // Check the out-edges.
    for(const edge_descriptor e: tangle.outEdges) {
        if(graph[e].outTangle == reverseComplementTangleId) {
            return true;
        }
    }

    // If getting here, we did not find a collision between this tangle
    // and its reverse complement.
    return false;
}



void AssemblyPathGraph2::writeHtml(const string& fileName) const
{
    ofstream html(fileName);
    writeHtml(html);
}



void AssemblyPathGraph2::writeHtml(ostream& html) const
{
    writeHtmlBegin(html, "Assembly path graph");
    html << "<body>"
        "<h1>Assembly path graph</h1>";
    writeVerticesHtml(html);
    writeEdgesHtml(html);
    writeTanglesHtml(html);
    html << "</body>";
}



void AssemblyPathGraph2::writeVerticesHtml(ostream& html) const
{
    const AssemblyPathGraph2& graph = *this;

    html << "<h2>Vertices</h2>"
        "<p>Each vertex corresponds to a vertex of the assembly graph."
        "<p><table><tr>"
        "<th>Id"
        "<th>Id of<br>reverse<br>complement<br>vertex";

    BGL_FORALL_VERTICES(v, graph, AssemblyPathGraph2) {
        const AssemblyPathGraph2Vertex& vertex = graph[v];
        const vertex_descriptor vReverseComplement = vertex.reverseComplementVertex;
        html <<
            "<tr id=v" << vertex.vertexId << ">" <<
            "<td class=centered>" << vertex.vertexId <<
            "<td class=centered>" << graph[vReverseComplement].vertexId;
    }

    html << "</table>";
}



void AssemblyPathGraph2::writeEdgesHtml(ostream& html) const
{
    const AssemblyPathGraph2& graph = *this;

    html << "<h2>Edges</h2>"
        "<p>Each edge corresponds to a path in the assembly graph vertex."
        "<p><table><tr>"
        "<th>Path"
        "<th>Path of<br>reverse<br>complement<br>edge"
        "<th>Source<br>vertex"
        "<th>Target<br>vertex"
        "<th>Path<br>length<br>(markers)"
        "<th>In-tangle"
        "<th>Tangle"
        "<th>Out-tangle";



    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph2) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const AssemblyGraph::VertexId vertexId0 = graph[v0].vertexId;
        const AssemblyGraph::VertexId vertexId1 = graph[v1].vertexId;
        const AssemblyPathGraph2Edge& edge = graph[e];
        const edge_descriptor eReverseComplement = edge.reverseComplementEdge;
        html <<
            "<tr id='e" << edge << "'>"
            "<td class=centered>" << edge <<
            "<td class=centered>" << graph[eReverseComplement] <<
            "<td class=centered><a href='#v" << vertexId0 << "'>" << vertexId0 << "</a>" <<
            "<td class=centered><a href='#v" << vertexId1 << "'>" << vertexId1 << "</a>" <<
            "<td class=centered>" << edge.pathLength;

        html <<  "<td class=centered>";
        if(edge.inTangle != invalidTangle2Id) {
            html << "<a href='#t" << edge.inTangle << "'>" << edge.inTangle << "</a>";
        }

        html <<  "<td class=centered>";
        if(edge.tangle != invalidTangle2Id) {
            html << "<a href='#t" << edge.tangle << "'>" << edge.tangle << "</a>";
        }

        html <<  "<td class=centered>";
        if(edge.outTangle != invalidTangle2Id) {
            html << "<a href='#t" << edge.outTangle << "'>" << edge.outTangle << "</a>";
        }
    }

    html << "</table>";

}



void AssemblyPathGraph2::writeTanglesHtml(ostream& html) const
{
    const AssemblyPathGraph2& graph = *this;

    html << "<h2>Tangles</h2>"
        "A tangle is generated by each edge v<sub>0</sub>&rarr;v<sub>1</sub> "
        "for which the source vertex v<sub>0</sub> has in-degree greater than 1 "
        "and the target vertex v<sub>1</sub> has out-degree "
        "greater than 1."
        "<p><table><tr>"
        "<th>Id"
        "<th>In-edges"
        "<th>Tangle<br>edge"
        "<th>Out-edges"
        "<th>Tangle<br>matrix"
        "<th>Solvable?";

    for(const auto& p: tangles) {
        const Tangle2& tangle = p.second;
        html <<
            "<tr id=t" << tangle.tangleId << ">"
            "<td class=centered>" << tangle.tangleId;

        // In-edges.
        html << "<td class=centered>";
        for(const edge_descriptor e: tangle.inEdges) {
            html << "<a href='#e" << graph[e] << "'>" <<
                graph[e] << "</a>" << " ";
        }

        // Tangle2 edge.
        html << "<td class=centered><a href='#e" << graph[tangle.edge] << "'>" <<
            graph[tangle.edge] << "</a>";

        // Out-edges.
        html << "<td class=centered>";
        for(const edge_descriptor e: tangle.outEdges) {
            html << "<a href='#e" << graph[e] << "'>" <<
                graph[e] << "</a>" << " ";
        }



        // Tangle matrix.
        html << "<td class=centered>";
        html << "<table style='margin-left:auto;margin-right:auto;'>";
        html << "<tr><td class=centered>";
        for(uint64_t j=0; j<tangle.outEdges.size(); j++) {
            const edge_descriptor e = tangle.outEdges[j];
            html << "<td class=centered>" << graph[e];
        }
        for(uint64_t i=0; i<tangle.inEdges.size(); i++) {
            const edge_descriptor e = tangle.inEdges[i];
            html << "<tr><td class=centered>" << graph[e];
            for(uint64_t j=0; j<tangle.outEdges.size(); j++) {
                const uint64_t value = tangle.matrix[i][j];
                html << "<td class=centered>";
                if(value) {
                    html << value;
                }
            }
        }
        html << "</table>";



        html << "<td class=centered>";
        if(tangle.isSolvable) {
            html << "Yes";
        } else {
            html << "No";
        }

    }


    html << "</table>";

}



// GFA output (without sequence).
void AssemblyPathGraph2::writeGfa(
    const string& fileName,
    double basesPerMarker) const
{
    ofstream gfa(fileName);
    writeGfa(gfa, basesPerMarker);
}



void AssemblyPathGraph2::writeGfa(
    ostream& gfa,
    double basesPerMarker) const
{
    const AssemblyPathGraph2& graph = *this;

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";


    // Write a segment record for each edge.
    BGL_FORALL_EDGES(e, graph, AssemblyPathGraph2) {
        gfa <<
            "S\t" <<
            graph[e] << "\t" <<
            "*\t" <<
            "LN:i:" << uint64_t(basesPerMarker * double(graph[e].pathLength)) <<
            "\n";
    }



    // Write GFA links.
    // For each vertex in the assembly path graph there is a link for
    // each combination of in-edges and out-edges.
    // Therefore each vertex generates a number of
    // links equal to the product of its in-degree and out-degree.
    BGL_FORALL_VERTICES(v, graph, AssemblyPathGraph2) {
        BGL_FORALL_INEDGES(v, eIn, graph, AssemblyPathGraph2) {
            BGL_FORALL_OUTEDGES(v, eOut, graph, AssemblyPathGraph2) {
                gfa <<
                    "L\t" <<
                    graph[eIn] << "\t" <<
                    "+\t" <<
                    graph[eOut] << "\t" <<
                    "+\t" <<
                    "*\n";
            }
        }
    }
}


void AssemblyPathGraph2::removeIsolatedVertices()
{
    AssemblyPathGraph2& graph = *this;
    vector<vertex_descriptor> isolatedVertices;

    BGL_FORALL_VERTICES(v, graph, AssemblyPathGraph2) {
        if(in_degree(v, graph)==0 and out_degree(v, graph)==0) {
            isolatedVertices.push_back(v);
        }
    }

    for(const vertex_descriptor v: isolatedVertices) {
        clear_vertex(v, graph);
        remove_vertex(v, graph);
    }
}
