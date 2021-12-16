#ifndef SHASTA_PHASING_GRAPH_HPP
#define SHASTA_PHASING_GRAPH_HPP

// The PhasingGraph is used for hierarchical phasing.
// Each vertex represents a set of bubbles already phased
// relative to each other.


// Shasta.
#include "MultithreadedObject.hpp"
#include "ReadId.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>

// Standard library.
#include "algorithm.hpp"
#include "array.hpp"
#include "cstdint.hpp"
#include <limits>
#include "tuple.hpp"
#include "utility.hpp"
#include "vector.hpp"

namespace shasta {

    class AssemblyGraph2Vertex;
    class AssemblyGraph2Edge;
    class AssemblyGraph2;

    using AssemblyGraph2BaseClass =
        boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
        AssemblyGraph2Vertex, AssemblyGraph2Edge>;

    class PhasingGraphVertex;
    class PhasingGraphEdge;
    class PhasingGraph;

    using PhasingGraphBaseClass =
        boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
        PhasingGraphVertex, PhasingGraphEdge>;

}




class shasta::PhasingGraphVertex {
public:

    // The bubbles of this vertex.
    // Each vertex is stored with a phase, which can be 0 or 1.
    // The phase is 0 if the bubble is in phase with the vertex
    // and 1 otherwise (the two branches need to be swapped for it to be in phase).
    vector<pair<AssemblyGraph2BaseClass::edge_descriptor, uint64_t> > bubbles;

    // The oriented reads on each side of the bubbles of this vertex
    // (after swapping sides for out of phase bubbles).
    array< vector<OrientedReadId>, 2> orientedReadIds;

    // The connected component this vertex belongs to.
    static const uint64_t invalidComponentId = std::numeric_limits<uint64_t>::max();
    uint64_t componentId = invalidComponentId;
    bool isPhased() const
    {
        return componentId != invalidComponentId;
    }

    // The phase assigned to this vertex.
    // It is only meaningful within each connected component.
    static const uint64_t invalidPhase = std::numeric_limits<uint64_t>::max();
    uint64_t phase = invalidPhase;
};



class shasta::PhasingGraphEdge {
public:
    // Store the number of common oriented reads.
    // matrix[sideA][sideB] stores the number of OrientedReadIds
    // that appear on sideA of the "first" vertex and on
    // sideB of the "second" vertex of this edge.
    // The "first" vertex of the edge is the lowered numbered.
    array<array<uint64_t, 2>, 2> matrix;

    PhasingGraphEdge()
    {
        for(uint64_t i=0; i<2; i++) {
            for(uint64_t j=0; j<2; j++) {
                matrix[i][j] = 0;
            }
        }
    }


    // Results of the Bayesian model computed using a call to diploidBayesianPhase.
    // Prandom = probability of the random hypothesis
    // Pin = probability of the in phase hypothesis
    // Pout = probability of the out of phase hypothesis
    // The meaning of logP is different for removal of bad bubbles
    // (allowRandomHypothesis=true) and for phasing (allowRandomHypothesis=false).
    // See PhasingGraphEdge::runBayesianModel for details.
    double logPin;  // log(Pin  / Prandom) in dB
    double logPout; // log(Pout / Prandom) in dB
    double logP;
    uint64_t relativePhase; // 0 = in phase, 1 = out of phase
    void runBayesianModel(double epsilon, bool allowRandomHypothesis);



    bool isTreeEdge = false;

    uint64_t diagonalCount() const
    {
        return matrix[0][0] + matrix[1][1];
    }
    uint64_t offDiagonalCount() const
    {
        return matrix[0][1] + matrix[1][0];
    }
    uint64_t totalCount() const
    {
        return diagonalCount() + offDiagonalCount();
    }
    uint64_t concordantCount() const
    {
        return max(diagonalCount(), offDiagonalCount());
    }
    uint64_t discordantCount() const
    {
        return min(diagonalCount(), offDiagonalCount());
    }
    int64_t delta() const
    {
        const int64_t d = int64_t(concordantCount()) - int64_t(discordantCount());
        SHASTA_ASSERT(d >= 0);
        return d;
    }

};



class shasta::PhasingGraph:
    public PhasingGraphBaseClass,
    public MultithreadedObject<PhasingGraph> {
public:
    PhasingGraph(
        const AssemblyGraph2&,
        uint64_t minConcordantReadCount,
        uint64_t maxDiscordantReadCount,
        double minLogP,
        double epsilon,
        size_t threadCount,
        bool allowRandomHypothesis);

    // Find the optimal spanning tree using logFisher as the edge weight.
    // Edges that are part of the optimal spanning tree get their
    // isTreeEdge set.
    void computeSpanningTree();

    // Phase vertices using the spanning tree.
    void phase();

    // Store the phasing in the AssemblyGraph2.
    void storePhasing(AssemblyGraph2&) const;

    void writeCsv(const string& baseName, const AssemblyGraph2&) const;
    void writeVerticesCsv(const string& fileName) const;
    void writeVerticesDetailsCsv(const string& fileName, const AssemblyGraph2&) const;
    void writeEdgesCsv(const string& fileName, const AssemblyGraph2&) const;
    void writeGraphviz(const string& fileName) const;

private:
    void createVertices(const AssemblyGraph2&);

    // Edge creation is expensive and runs in parallel.
    void createEdges(
        uint64_t minConcordantReadCount,
        uint64_t maxDiscordantReadCount,
        double minLogP,
        double epsilon,
        size_t threadCount,
        bool allowRandomHypothesis);
    void createEdgesThreadFunction(size_t threadId);
    class CreateEdgesData {
    public:
        uint64_t minConcordantReadCount;
        uint64_t maxDiscordantReadCount;
        double minLogP;
        double epsilon; // For Bayesian model.
        bool allowRandomHypothesis;
        vector<PhasingGraph::vertex_descriptor> allVertices;
        class EdgeData {
        public:
            PhasingGraph::vertex_descriptor vB;
            uint64_t sideA;
            uint64_t sideB;
            bool operator<(const EdgeData& that) const
            {
                return vB < that.vB;
            }
        };
    };
    CreateEdgesData createEdgesData;
    void createEdges(
        PhasingGraph::vertex_descriptor,
        uint64_t minConcordantReadCount,
        uint64_t maxDiscordantReadCount,
        double minLogP,
        double epsilon,
        vector<CreateEdgesData::EdgeData>&,
        vector< tuple<vertex_descriptor, vertex_descriptor, PhasingGraphEdge> >& threadEdges,
        bool allowRandomHypothesis);

    // Get the vertex corresponding to a component, creating it if necessary.
    PhasingGraph::vertex_descriptor getVertex(uint64_t componentId);

    // A table that, for each OrientedReadId, contains a list of
    // pairs(vertex, side) that the OrientedReadId appears on.
    // Indexed by OrientedReadId::getValue().
    vector< vector< pair<PhasingGraph::vertex_descriptor, uint64_t> > > orientedReadsTable;
    void createOrientedReadsTable(uint64_t readCount);

};

#endif
