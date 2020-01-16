#ifndef SHASTA_CONFLICT_READ_GRAPH_HPP
#define SHASTA_CONFLICT_READ_GRAPH_HPP


/*******************************************************************************

The conflict read graph is an undirected graph in which each vertex
corresponds to an oriented read. Therefore, a read is associated with
two vertices, one for each orientation (Strand).

We create an undirected edge between two vertices if we found a bad
induced alignment alignment between the corresponding oriented reads.
The existence of the bad induced alignment indicates that, even though
the two oriented reads are close in the read graph, they are likely
to originate from different regions of the genome.

We use approximate coloring of the conflict graph to separate reads
based on the regions of the genome they originate from.
This is used to prevent tangles in the marker graph and in the assembly,
where two or more similar region of the genome come together due to
sequence similarity.

*******************************************************************************/

// Shasta.
#include "MemoryMappedUndirectedGraph.hpp"
#include "ReadId.hpp"

namespace shasta {
    class Assembler;
    class ConflictReadGraph;
    class ConflictReadGraphEdge;
    class ConflictReadGraphVertex;


    using ConflictReadGraphBaseClass =
        MemoryMapped::UndirectedGraph<ConflictReadGraphVertex, ConflictReadGraphEdge>;
}



// A vertex of the conflict read graph.
class shasta::ConflictReadGraphVertex {
public:

    static const uint64_t invalid = std::numeric_limits<uint64_t>::max();

    // The color within this connected component, a positive number.
    // or invalid if the vertex is isolated and makes up its own trivial
    // connected component.
    uint64_t color = invalid;
};



// An edge of the conflict read graph.
class shasta::ConflictReadGraphEdge {
public:

};



class shasta::ConflictReadGraph :
    public ConflictReadGraphBaseClass {
public:
    using BaseClass = ConflictReadGraphBaseClass;
    using Vertex = ConflictReadGraphVertex;
    using Edge = ConflictReadGraphEdge;

    void createVertices(ReadId readCount)
    {
        vertices.resize(2 * readCount);
    }

    void colorConnectedComponent(const vector<VertexId>&);

    // Convert a VertexId to an OrientedReadId and vice versa.
    static OrientedReadId getOrientedReadId(VertexId vertexId)
    {
        return OrientedReadId(OrientedReadId::Int(vertexId));
    }
    static VertexId getVertexId(OrientedReadId orientedReadId)
    {
        return orientedReadId.getValue();
    }

private:


};

#endif

