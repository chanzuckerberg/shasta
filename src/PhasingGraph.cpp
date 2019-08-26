#include "PhasingGraph.hpp"
using namespace shasta;



PhasingGraph::PhasingGraph() :
    MultithreadedObject<PhasingGraph>(*this)
{
}



// Function to construct names for binary objects.
string PhasingGraph::dataName(const string& name) const
{
    if(dataFileNamePrefix.empty()) {
        return "";  // Anonymous;
    } else {
        return dataFileNamePrefix + name;
    }
}



void PhasingGraph::findSimilarPairs(
    size_t threadCount,
    double phasingSimilarityThreshold)
{
    // Store the tghreshold.
    similarityThreshold = phasingSimilarityThreshold;

    // Allocate space for the similar pairs found by each thread.
    threadSimilarPairs.resize(threadCount);

    // See how many reads we have.
    const uint64_t orientedReadCount = assemblyGraphEdges.size();
    const uint64_t readCount = orientedReadCount / 2;

    // Find the similar pairs with readId0<readId1 in multithreaded code.
    setupLoadBalancing(readCount, 1000);
    runThreads(&PhasingGraph::findSimilarPairs, threadCount);

    // Gather the pairs found by each thread.
    similarPairs.createNew(
        dataName("PhasingSimilarPairs"), dataPageSize);
    for(size_t threadId=0; threadId!=threadCount; threadId++) {
        auto& threadPairs = *threadSimilarPairs[threadId];
        for(const auto threadPair: threadPairs) {
            similarPairs.push_back(threadPair);
        }
        threadPairs.remove();
    }
    threadSimilarPairs.clear();
    cout << timestamp << "Found " << similarPairs.size() <<
        " oriented read pairs with phasing similarity " <<
        phasingSimilarityThreshold << " or more." << endl;



    // Write all the pairs in graphviz format.
    ofstream graphout("PhasingGraph.dot");
    graphout << "graph G {" << endl;
    for(const auto& p: similarPairs) {
        const auto& q = p.first;
        OrientedReadId orientedReadId0(q.readIds[0], 0);
        OrientedReadId orientedReadId1(q.readIds[1], q.isSameStrand ? 0 : 1);
        graphout << orientedReadId0.getValue();
        graphout << "--";
        graphout << orientedReadId1.getValue();
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        graphout << ";\n";
        graphout << orientedReadId0.getValue();
        graphout << "--";
        graphout << orientedReadId1.getValue();
        graphout << ";\n";
    }
    graphout << "}" << endl;

}



void PhasingGraph::findSimilarPairs(size_t threadId)
{
    // Allocate space for the similar pairs found by this thread.
    threadSimilarPairs[threadId] =
        std::make_shared<MemoryMapped::Vector< pair<OrientedReadPair, float> > >();
    MemoryMapped::Vector< pair<OrientedReadPair, float> >& thisThreadSimilarPairs =
        *threadSimilarPairs[threadId];
    thisThreadSimilarPairs.createNew(
        dataName("tmp-PhasingThreadSimilarPairs-" + to_string(threadId)), dataPageSize);

    // Work area, defined here to reduce memory allocation activity.
    vector<OrientedReadId> orientedReadIds;

    // Loop over batches allocated to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads allocated to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); ++readId) {
            findSimilarPairs(readId, thisThreadSimilarPairs, orientedReadIds);

        }
    }
}



double PhasingGraph::computePhasingSimilarity(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1)
{
    using EdgeId = AssemblyGraph::EdgeId;

    // Access the assembly graph edges that these two oriented reads
    // are internal to.
    const MemoryAsContainer<EdgeId> edges0 =
        assemblyGraphEdges[orientedReadId0.getValue()];
    const MemoryAsContainer<EdgeId> edges1 =
        assemblyGraphEdges[orientedReadId1.getValue()];

    // Count the number of common assembly graph edges.
    uint64_t intersectionSize = 0;
    auto begin0 = edges0.begin();
    auto begin1 = edges1.begin();
    auto end0   = edges0.end();
    auto end1   = edges1.end();
    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        const EdgeId edgeId0 = *it0;
        const EdgeId edgeId1 = *it1;
        if(edgeId0 < edgeId1) {
            ++it0;
        } else if(edgeId1 < edgeId0) {
            ++it1;
        } else {
            ++intersectionSize;
            ++it0;
            ++it1;
        }
    }

    const uint64_t n0 = edges0.size();
    const uint64_t n1 = edges1.size();
    SHASTA_ASSERT(n0 > 0);
    SHASTA_ASSERT(n1 > 0);
    return double(intersectionSize) / double(min(n0, n1));

}



// Find similar pairs in which the first ReadId is readId0
// and the second ReadId is readId1>readId0.
// The last argument is a work area, passed in
// to reduce memory allocation activity.
void PhasingGraph::findSimilarPairs(
    ReadId readId0,
    MemoryMapped::Vector< pair<OrientedReadPair, float> >& pairVector,
    vector<OrientedReadId> orientedReadIds
    )
{
    using EdgeId = AssemblyGraph::EdgeId;

    // Work with readId0 on strand 0.
    const OrientedReadId orientedReadId0(readId0, 0);
    // cout << "Working on " << readId0 << endl;

    // Store in the orientedReadIds vector
    // all oriented reads with readId>readId0
    // that share at least one assembly graph edge with orientedReadId0.
    const MemoryAsContainer<EdgeId> edges0 = assemblyGraphEdges[orientedReadId0.getValue()];
    orientedReadIds.clear();
    for(const EdgeId edge0: edges0) {

        // Loop over oriented read ids internal to this edge.
        const MemoryAsContainer<OrientedReadId> orientedReadIds1 =
            orientedReads[edge0];
        for(const OrientedReadId orientedReadId1: orientedReadIds1) {
            if(orientedReadId1.getReadId() > readId0) {
                orientedReadIds.push_back(orientedReadId1);
            }
        }
    }

    // Deduplicate.
    sort(orientedReadIds.begin(), orientedReadIds.end());
    orientedReadIds.resize(
        unique(orientedReadIds.begin(), orientedReadIds.end()) - orientedReadIds.begin());

    // Now compute phasing similarity for each.
    for(const OrientedReadId orientedReadId1: orientedReadIds) {
        const double similarity = computePhasingSimilarity(orientedReadId0, orientedReadId1);
        if(similarity >= similarityThreshold) {
            const bool isSameStrand = orientedReadId1.getStrand() == 0;
            pairVector.push_back(make_pair(
                OrientedReadPair(readId0, orientedReadId1.getReadId(), isSameStrand),
                float(similarity)
                ));
            // cout << readId0 << " " << orientedReadId1.getReadId() << " " << int(isSameStrand) << "\n";
        }
    }

}

