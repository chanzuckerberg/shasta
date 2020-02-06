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
// Each vertex corresponds to an oriented read.
class shasta::ConflictReadGraphVertex {
public:

    static const uint32_t invalid = std::numeric_limits<uint32_t>::max();


    // The cluster id assigned during coloring of the conflict read graph.
    uint32_t clusterId = invalid;

    bool hasValidClusterId() const
    {
        return clusterId != invalid;
    }

    // Number of markers preceding the first marker with non-zero marker coverage
    // (that is, the first marker associated with a marker graph vertex).
    uint32_t leftTrim;

    // Number of markers following the last marker with non-zero marker coverage
    // (that is, the last marker associated with a marker graph vertex).
    uint32_t rightTrim;

    // The length of the longest internal streak of markers with zero marker coverage
    // (that is, not associated with a marker graph vertex).
    // This excludes the streaks at the beginning and end, which are
    // described by leftTrim and rightTrim.
    // If this is too long, the read is considered pathological and
    // excluded from the assembly.
    uint32_t longestGap;

    // This is set is longestGap>maxSkip.
    bool hasLongGap = false;

    // If set (by  cleanupConflictReadGraph), the corresponding read graph
    // edge is effectively excluded from assembly by marking all
    // edges incident to it as conflict edges.
    bool wasRemoved = false;
};



// An edge of the conflict read graph.
class shasta::ConflictReadGraphEdge {
public:

    // If set (by  cleanupConflictReadGraph), this edge is ignored
    // by markDirectedReadGraphConflictEdges2.
    bool wasRemoved = false;
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

    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;

private:


};

#endif

