#ifndef SHASTA_COMPRESSED_ASSEMBLY_GRAPH_HPP
#define SHASTA_COMPRESSED_ASSEMBLY_GRAPH_HPP

/*******************************************************************************

The compressed assembly graph is a representation of the assembly graph
in which each linear sequence of bubbles is compressed to a single edge.

*******************************************************************************/

#include "AssemblyGraph.hpp"
#include <boost/graph/adjacency_list.hpp>

namespace shasta {
    class CompressedAssemblyGraph;
    class CompressedAssemblyGraphEdge;
    class CompressedAssemblyGraphVertex;

    using CompressedAssemblyGraphBaseClass =
        boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        CompressedAssemblyGraphVertex,
        CompressedAssemblyGraphEdge
        >;
}


class shasta::CompressedAssemblyGraphVertex {
public:
    AssemblyGraph::VertexId vertexId;

    CompressedAssemblyGraphVertex(AssemblyGraph::VertexId vertexId) :
        vertexId(vertexId) {}
};



class shasta::CompressedAssemblyGraphEdge {
public:

    // The chain of assembly graph vertices associated
    // with this edge.
    // This includes the assembly graph vertices
    // associated with the source and target of this edge.
    vector<AssemblyGraph::VertexId> vertices;
};



class shasta::CompressedAssemblyGraph :
    public CompressedAssemblyGraphBaseClass {
public:
    using VertexId = AssemblyGraph::VertexId;


    // Create the CompressedAssemblyGraph from the AssemblyGraph.
    CompressedAssemblyGraph(const AssemblyGraph&);

};


#endif
