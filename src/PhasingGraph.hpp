#ifndef SHASTA_PHASING_GRAPH_HPP
#define SHASTA_PHASING_GRAPH_HPP

/*******************************************************************************

We keep track of which oriented reads are internal to 
which assembly graph edges. We then define a phasing similarity
between oriented reads: two oriented reads have high phAsing similarity
if the sets of assembly graph edges they are internal to are similar.
Two reads with high similarity are likely to be phased together.

Using this similarity, we define the Phasing graph, an undirected
graph in which each edge represents an oriented read.
An edge is added between vertices corresponding to oriented reads
with high phasing similarity. 
 
*******************************************************************************/

#include "AssemblyGraph.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultitreadedObject.hpp"
#include "OrientedReadPair.hpp"
#include "ReadId.hpp"

namespace shasta {
    class PhasingGraph;
}



class shasta::PhasingGraph :
    public MultithreadedObject<PhasingGraph> {
public:    

    PhasingGraph();

    // The oriented reads internal to each assembly graph edge.
    // Indexed by assembly graph EdgeId.
    MemoryMapped::VectorOfVectors<OrientedReadId, uint64_t> orientedReads;

    // The assembly graph edges that each oriented read is internal to.
    // Indexed by OrientedReadId::getValue().
    MemoryMapped::VectorOfVectors<AssemblyGraph::EdgeId, uint64_t> assemblyGraphEdges;



#if 0
    // For an oriented read r, a turn in the assembly graph
    // is a path of length two (two consecutive edges e0, e1)
    // such that:
    // 1. v = target(e0) = source(e1), called the hinge of the turn.
    // 2. r is internal to both e0 and e1.
    // 3. r is not internal to any other incoming edges  of v.
    // 4. r is not internal to any other outgoing edges of v.
    // Note that, because of 3. and 4., r can have at most one
    // turn for each vertex v.
    class Turn {
    public:

        // The incoming edge.
        AssemblyGraph::EdgeId e0;

        // The outgoing edge.
        AssemblyGraph::EdgeId e1;

        // The hinge. v = target(e0) = source(e1).
        AssemblyGraph::VertexId v;

        Turn(
            AssemblyGraph::EdgeId e0 = AssemblyGraph::invalidEdgeId,
            AssemblyGraph::EdgeId e1 = AssemblyGraph::invalidEdgeId,
            AssemblyGraph::VertexId v = AssemblyGraph::invalidVertexId) :
            e0(e0), e1(e1), v(v) {}

        bool operator==(const Turn& that) const
        {
            // Can't compare turns on different vertices.
            SHASTA_ASSERT(v == that.v);

            return e0==that.e0 && e1==that.e1;
        }
    };

    // Store the turns for each oriented read.
    // This is indexed by OrientedReadId::getValue().
    // For a given oriented read, the turns are sorted in
    // order of increasing v, and there can be at most
    // one for each v.
    MemoryMapped::VectorOfVectors<Turn, uint64_t> turns;
    void findTurns();

    // Count concordant and discordant turns for two oriented reads.
    // Returns pair(concordantCount, discordantCount).
    pair<uint64_t, uint64_t> countCommonTurns(OrientedReadId, OrientedReadId);
#endif


    // Oriented read pairs with phasing similarity greater than the threshold used.
    // We only store the ones with readId0 < readId1.
    // Each pair is stored with its phasinbg similarity.
    MemoryMapped::Vector< pair<OrientedReadPair, float> > similarPairs;
    void findSimilarPairs(size_t threadCount, double phasingSimilarityThreshold);
    void keepBestSimilarPairs(int maxNeighborCount);

    // Same as above, for the similar pairs found by each read.
    // This is only used inside findSimilarPairs.
    vector<std::shared_ptr<MemoryMapped::Vector< pair<OrientedReadPair, float> > > > threadSimilarPairs;

    // File name prefix and page size for binary data.
    string dataFileNamePrefix;
    size_t dataPageSize;

    void writeGraphviz();

    double computePhasingSimilarity(OrientedReadId, OrientedReadId);

private:

    string dataName(const string& name) const;
    void findSimilarPairs(size_t threadId);
    double similarityThreshold;

    // Find similar pairs in which the first ReadId is readId0
    // and the second ReadId is readId1>readId0.
    // The last argument is a work area, passed in
    // to reduce memory allocation activity.
    void findSimilarPairs(
        ReadId readId0,
        MemoryMapped::Vector< pair<OrientedReadPair, float> >& pairVector,
        vector<OrientedReadId> orientedReadIds);



};

#endif
