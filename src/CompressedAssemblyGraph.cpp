#include "CompressedAssemblyGraph.hpp"
#include "findLinearChains.hpp"
using namespace shasta;

#include "vector.hpp"



// Create the CompressedAssemblyGraph from the AssemblyGraph.
CompressedAssemblyGraph::CompressedAssemblyGraph(
    const AssemblyGraph& assemblyGraph)
{
    CompressedAssemblyGraph& graph = *this;

    cout << "The assembly graph has " << assemblyGraph.vertices.size() <<
        " vertices and " << assemblyGraph.edges.size() << " edges." << endl;


    // Create a vertex for each vertex of the Assembly graph.
    vector<vertex_descriptor> vertexTable(assemblyGraph.vertices.size(), null_vertex());
    for(VertexId vertexId=0; vertexId<assemblyGraph.vertices.size(); vertexId++) {
        const vertex_descriptor v = add_vertex(CompressedAssemblyGraphVertex(vertexId), graph);
        vertexTable[vertexId] = v;
    }

    // Create an edge for each set of parallel edges of the assembly graph.
    for(const AssemblyGraph::Edge& edge: assemblyGraph.edges) {
        const vertex_descriptor v0 = vertexTable[edge.source];
        const vertex_descriptor v1 = vertexTable[edge.target];
        bool edgeExists = false;
        edge_descriptor e;
        tie(e, edgeExists) = boost::edge(v0, v1, graph);
        if(not edgeExists) {
            bool edgeWasAdded = false;
            tie(e, edgeWasAdded) = add_edge(v0, v1, graph);
            SHASTA_ASSERT(edgeWasAdded);
            graph[e].vertices.push_back(graph[v0].vertexId);
            graph[e].vertices.push_back(graph[v1].vertexId);
        }
    }

    // Find linear chains.
    vector< std::list<edge_descriptor> > chains;
    findLinearChains(graph, chains);



    // Replace each chain with a single edge.
    for(const std::list<edge_descriptor>& chain: chains) {

        // If the chain has length 1, leave it alone.
        if(chain.size() == 1) {
            continue;
        }

        // Add the new edge.
        const vertex_descriptor v0 = source(chain.front(), graph);
        const vertex_descriptor v1 = target(chain.back(), graph);
        edge_descriptor eNew;
        bool edgeWasAdded = false;
        tie(eNew, edgeWasAdded) = add_edge(v0, v1, graph);
        SHASTA_ASSERT(edgeWasAdded);
        CompressedAssemblyGraphEdge& newEdge = graph[eNew];

        // Fill in the assembly graph vertices corresponding to this new edge.
        newEdge.vertices.push_back(graph[v0].vertexId);
        for(const edge_descriptor e: chain) {
            const vertex_descriptor v = target(e, graph);
            newEdge.vertices.push_back(graph[v].vertexId);
        }

        // Remove the edges of the chain.
        // We will remove the vertices later.
        for(const edge_descriptor e: chain) {
            boost::remove_edge(e, graph);
        }
    }



    // Remove the vertices that have become isolated.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, CompressedAssemblyGraph) {
        if(in_degree(v, graph) == 0 and out_degree(v, graph) == 0) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        remove_vertex(v, graph);
    }
    cout << "The compressed assembly graph has " <<
        num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;

}
