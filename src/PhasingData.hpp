#ifndef SHASTA_PHASING_DATA_HPP
#define SHASTA_PHASING_DATA_HPP

/*******************************************************************************

Class PhasingData is used to keep track of which oriented reads are internal to
which assembly graph edges. There are data structures to find the oriented
reads internal to an assembly graph edge, and which assembly graph edges
each oriented read is internal to.

Given two assembly graph edges, there are functions that compute
the number of common internal reads, as well as the Jaccard similarity
of the two sets of internal reads.

*******************************************************************************/

#include "AssemblyGraph.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "MultitreadedObject.hpp"
// #include "OrientedReadPair.hpp"
#include "ReadId.hpp"

namespace shasta {
    class PhasingData;
}



class shasta::PhasingData :
    public MultithreadedObject<PhasingData> {
public:    

    PhasingData();

    // The oriented reads internal to each assembly graph edge.
    // Indexed by assembly graph EdgeId.
    MemoryMapped::VectorOfVectors<OrientedReadId, uint64_t> orientedReads;

    // Compute the phasing similarity of two assembly graph edges.
    // It is computed as the Jaccard similarity of the sets
    // of oriented reads internal to each of the two assembly graph edges.
    double computePhasingSimilarity(AssemblyGraph::EdgeId, AssemblyGraph::EdgeId);
    uint64_t countCommonInternalOrientedReads(AssemblyGraph::EdgeId, AssemblyGraph::EdgeId);



    // The assembly graph edges that each oriented read is internal to.
    // Indexed by OrientedReadId::getValue().
    // For each oriented read, they are stored sorted.
    MemoryMapped::VectorOfVectors<AssemblyGraph::EdgeId, uint64_t> assemblyGraphEdges;


#if 0
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
#endif


};

#endif
