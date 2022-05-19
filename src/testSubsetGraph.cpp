#include "SubsetGraph.hpp"
#include "testSubsetGraph.hpp"
using namespace shasta;

#include "iostream.hpp"
#include "fstream.hpp"
#include "utility.hpp"



void shasta::testSubsetGraph()
{
    using Graph = SubsetGraph<int>;
    Graph graph;

    vector<Graph::vertex_descriptor> v;
    v.push_back(graph.addVertex(vector<int>{10, 20, 30}));
    v.push_back(graph.addVertex(vector<int>{10, 20}));
    v.push_back(graph.addVertex(vector<int>{20, 30}));
    v.push_back(graph.addVertex(vector<int>{20}));
    v.push_back(graph.addVertex(vector<int>{20, 30, 40}));
    v.push_back(graph.addVertex(vector<int>{40, 50, 60}));
    v.push_back(graph.addVertex(vector<int>{30, 40}));
    v.push_back(graph.addVertex(vector<int>{40, 50}));
    v.push_back(graph.addVertex(vector<int>{40}));
    graph.createEdges();

    cout << "The subset graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;


    ofstream dot("SubsetGraph.dot");
    dot << "digraph SubsetGraph {\n";
    std::map<Graph::vertex_descriptor, int> indexMap;
    int index = 0;;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        dot << index << " [label=\"";
        for(int i=0; i<int(graph[v].size()); i++) {
            dot << graph[v][i];
            if(i != int(graph[v].size())-1) {
                dot << " ";
            }
        }
        dot << "\"];\n";
        indexMap.insert(make_pair(v, index++));
    }
    BGL_FORALL_EDGES(e, graph, Graph) {
        const auto v0 = source(e, graph);
        const auto v1 = target(e, graph);
        dot << indexMap[v0] << "->" << indexMap[v1] << ";\n";
    }
    dot << "}\n";
}
