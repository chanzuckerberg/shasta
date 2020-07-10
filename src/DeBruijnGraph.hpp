#ifndef SHASTA_DEBRUIJN_GRAPH_HPP
#define SHASTA_DEBRUIJN_GRAPH_HPP

// A general purpose templated class representing a De Bruijn graph.
// Each vertex represents a sequence of k symbols occurring in
// sequences of symbols.
// Template arguments:
// - Symbol: the symbols that the sequences are made of.
// - k: the number of symbols in the sequence represented by a vertex.
// - SequenceId: the type used to identify the sequences.

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"
#include "string.hpp"
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {

    template<class Symbol, uint64_t k, class SequenceId> class DeBruijnGraph;
    template<class Symbol, uint64_t k, class SequenceId> class DeBruijnGraphVertex;
    template<class Symbol, uint64_t k, class SequenceId> class DeBruijnGraphEdge;

    template<class Symbol, uint64_t k, class SequenceId> using DeBruijnGraphBaseClass =
        boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        DeBruijnGraphVertex<Symbol, k, SequenceId>,
        DeBruijnGraphEdge<Symbol, k, SequenceId>
        >;

}



template<class Symbol, uint64_t k, class SequenceId> class shasta::DeBruijnGraphVertex {
public:

    // The sequence of k symbols represented by this vertex.
    using VertexSequence = array<Symbol, k>;
    VertexSequence vertexSequence;

    // The occurrences of the k symbols represented by this vertex
    // in the input sequences.
    // The second element of each pair is the position in the sequence
    // of the first of the k symbols.
    vector< pair<SequenceId, uint64_t > > occurrences;

    uint64_t vertexId;

    DeBruijnGraphVertex(
        const VertexSequence& vertexSequence,
        uint64_t vertexId) :
        vertexSequence(vertexSequence),
        vertexId(vertexId) {}

};



template<class Symbol, uint64_t k, class SequenceId> class shasta::DeBruijnGraphEdge {
public:

    // The k-1 symbols associated with this edge.
    // This are the last k-1 symbols of the source vertex
    // and also the first k-1 symbols of the target vertex.
    using EdgeSequence = array<Symbol, k-1>;
    EdgeSequence edgeSequence;

    DeBruijnGraphEdge(const EdgeSequence& edgeSequence) :
        edgeSequence(edgeSequence) {}


};



template<class Symbol, uint64_t k, class SequenceId> class shasta::DeBruijnGraph :
    public DeBruijnGraphBaseClass<Symbol, k, SequenceId> {
public:

    using Graph = DeBruijnGraph;
    using Vertex = DeBruijnGraphVertex<Symbol, k, SequenceId>;
    using Edge = DeBruijnGraphEdge<Symbol, k, SequenceId>;
    using BaseClass = DeBruijnGraphBaseClass<Symbol, k, SequenceId>;
    using vertex_descriptor = typename BaseClass::vertex_descriptor;
    using edge_descriptor = typename BaseClass::edge_descriptor;
    using Sequence = vector<Symbol>;
    using VertexSequence = array<Symbol, k>;
    using EdgeSequence = array<Symbol, k-1>;

    // The vertices, keyed by the sequence they contain.
    std::map<VertexSequence, vertex_descriptor> vertexMap;
    uint64_t nextVertexId = 0;



    // Add a sequence to the graph.
    void addSequence(
        SequenceId sequenceId,
        const Sequence& sequence)
    {
        Graph& graph = *this;

        // Loop over possible starting positions.
        for(uint64_t startPosition=0; startPosition + k <= sequence.size(); startPosition++) {

            // Extract the k symbols starting here.
            VertexSequence vertexSequence;
            for(uint64_t i=0; i<k; i++) {
                vertexSequence[i] = sequence[startPosition + i];
            }

            // Get the corresponding vertex, creating it if necessary.
            vertex_descriptor v;
            auto it = vertexMap.find(vertexSequence);
            if(it == vertexMap.end()) {
                v = add_vertex(Vertex(vertexSequence, nextVertexId++), graph);
                vertexMap.insert(make_pair(vertexSequence, v));
            } else {
                v = it->second;
            }

            // Store this occurrence of the k symbols.
            graph[v].occurrences.push_back(make_pair(sequenceId, startPosition));
        }
    }



    void removeLowCoverageVertices(uint64_t minCoverage)
    {
        Graph& graph = *this;

        vector<vertex_descriptor> verticesTobeRemoved;
        BGL_FORALL_VERTICES_T(v, graph, Graph) {
            if(graph[v].occurrences.size() < minCoverage) {
                verticesTobeRemoved.push_back(v);
            }
        }

        for(const vertex_descriptor v: verticesTobeRemoved) {
            clear_vertex(v, graph);
            remove_vertex(v, graph);
        }

    }



    // Create the edges of the De Bruijn graph.
    // This should be called after all sequences are added.
    void createEdges()
    {
        Graph& graph = *this;

        // Index the vertices by their first k-1 symbols.
        std::map<EdgeSequence, vector<vertex_descriptor> > vertexIndex;
        BGL_FORALL_VERTICES_T(v, graph, Graph) {
            const VertexSequence& vertexSequence = graph[v].vertexSequence;
            EdgeSequence edgeSequence;
            const auto begin = vertexSequence.begin();
            const auto end = begin + (k-1);
            copy(begin, end, edgeSequence.begin());
            vertexIndex[edgeSequence].push_back(v);
        }

        // Use the index to create the edges.
        BGL_FORALL_VERTICES_T(v0, graph, Graph) {
            const VertexSequence& vertexSequence0 = graph[v0].vertexSequence;
            EdgeSequence edgeSequence;
            const auto begin = vertexSequence0.begin() + 1;
            const auto end = vertexSequence0.end();
            copy(begin, end, edgeSequence.begin());

            for(const vertex_descriptor v1: vertexIndex[edgeSequence]) {
                edge_descriptor e;
                tie(e, ignore) = add_edge(v0, v1, Edge(edgeSequence), graph);
            }
        }
    }


    void writeGraphviz(const string& fileName) const
    {
        const Graph& graph = *this;
        ofstream s(fileName);

        s << "digraph DeBruijnGraph {\n";

        BGL_FORALL_VERTICES_T(v, graph, Graph) {
            s << graph[v].vertexId <<
                "[label=\""  <<
                graph[v].occurrences.size() <<
                "\"];\n";
        }

        BGL_FORALL_EDGES_T(e, graph, Graph) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            s << graph[v0].vertexId << "->";
            s << graph[v1].vertexId << ";\n";
        }

        s << "}\n";

    }
};

#endif
