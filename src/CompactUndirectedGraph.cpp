// Unit test for class CompactUndirectedGraph.

#include "CompactUndirectedGraph.hpp"
#include "shortestPath.hpp"


namespace shasta {
    class TestCompactUndirectedGraphEdge;
    class TestCompactUndirectedGraphVertex {
    public:
        int id;
        uint64_t distance;
        uint8_t color;
        TestCompactUndirectedGraphVertex(int id = 0) : id(id) {}
        CompactUndirectedGraph<TestCompactUndirectedGraphVertex, TestCompactUndirectedGraphEdge>::vertex_descriptor predecessor;
    };
    class TestCompactUndirectedGraphEdge {
    public:
        uint64_t weight;
        TestCompactUndirectedGraphEdge(uint64_t weight) :
            weight(weight) {}
    };
}


void shasta::testCompactUndirectedGraph1()
{
    using G = CompactUndirectedGraph<double, double>;
    using vertex_descriptor = G::vertex_descriptor;
    using edge_descriptor = G::edge_descriptor;

    G g;

    const vertex_descriptor v0 = g.addVertex(0.);
    const vertex_descriptor v1 = g.addVertex(1.);
    const vertex_descriptor v2 = g.addVertex(2.);
    const vertex_descriptor v3 = g.addVertex(3.);
    SHASTA_ASSERT(g[v1] == 1.);
    g.doneAddingVertices();
    SHASTA_ASSERT(g[v2] == 2.);
    cout << "Done adding vertices." << endl;


    const edge_descriptor e01 = g.addEdge(v0, v1, 10.);
    g.addEdge(v1, v2, 20.);
    const edge_descriptor e23 = g.addEdge(v2, v3, 30.);
    SHASTA_ASSERT(g[e01] == 10.);
    g.doneAddingEdges();
    SHASTA_ASSERT(g[e23] == 30.);
    cout << "Done adding edges." << endl;

    // g.dump(cout);

    // Test iteration over all vertices.
    cout << "Vertices:" << endl;
    BGL_FORALL_VERTICES(v, g, G) {
        cout << g[v] << endl;
    }

    // Test iteration over all edges.
    cout << "Edges:" << endl;
    BGL_FORALL_EDGES(e, g, G) {
        const vertex_descriptor v0 = source(e, g);
        const vertex_descriptor v1 = target(e, g);
        cout << g[v0] << " " << g[v1] << " " << g[e] << endl;
    }

    // Test iteration over out-edges.
    BGL_FORALL_VERTICES(v0, g, G) {
        cout << "Out-edges of vertex " << g[v0] << ":" << endl;
        BGL_FORALL_OUTEDGES(v0, e01, g, G) {
            SHASTA_ASSERT(source(e01, g) == v0);
            const vertex_descriptor v1 = target(e01, g);
            cout << g[v0] << " " << g[v1] << " " << g[e01] << endl;
        }
    }
}




void shasta::testCompactUndirectedGraph2()
{
    using G = CompactUndirectedGraph<TestCompactUndirectedGraphVertex, TestCompactUndirectedGraphEdge>;
    using vertex_descriptor = G::vertex_descriptor;

    G g;

    const vertex_descriptor v0 = g.addVertex(0);
    const vertex_descriptor v1 = g.addVertex(1);
    const vertex_descriptor v2 = g.addVertex(2);
    const vertex_descriptor v3 = g.addVertex(3);
    const vertex_descriptor v4 = g.addVertex(4);
    const vertex_descriptor v5 = g.addVertex(5);
    g.doneAddingVertices();

    // Shortest path from 0 to 5 is 0-1-3-5.
    g.addEdge(v0, v1, TestCompactUndirectedGraphEdge(4));
    g.addEdge(v0, v2, TestCompactUndirectedGraphEdge(2));
    g.addEdge(v1, v2, TestCompactUndirectedGraphEdge(5));
    g.addEdge(v1, v3, TestCompactUndirectedGraphEdge(10));
    g.addEdge(v2, v4, TestCompactUndirectedGraphEdge(3));
    g.addEdge(v4, v3, TestCompactUndirectedGraphEdge(4));
    g.addEdge(v3, v5, TestCompactUndirectedGraphEdge(11));
    g.doneAddingEdges();

    vector<vertex_descriptor> path;
    FindShortestPathQueue<G> q;
    findShortestPath(g, v0, v5, path, q);

    cout << "Shortest path:";
    for(const vertex_descriptor v: path) {
        cout << " " << g[v].id;
    }
    cout << endl;
}
