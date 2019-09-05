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
// #include "OrientedReadPair.hpp"
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

    // Compute the phasing similarity of two assembly graph edges.
    // It is computed as the Jaccard similarity of the sets
    // of oriented reads internal to each of the two assembly graph edges.
    double computePhasingSimilarity(AssemblyGraph::EdgeId, AssemblyGraph::EdgeId);



#if 0
    // The assembly graph edges that each oriented read is internal to.
    // Indexed by OrientedReadId::getValue().
    MemoryMapped::VectorOfVectors<AssemblyGraph::EdgeId, uint64_t> assemblyGraphEdges;


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

    uint64_t countCommonInternalOrientedReads(AssemblyGraph::EdgeId, AssemblyGraph::EdgeId);

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
#endif


};

#endif
