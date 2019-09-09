#include "PhasingData.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
#include "setOperations.hpp"
using namespace shasta;



PhasingData::PhasingData() :
    MultithreadedObject<PhasingData>(*this)
{
}



// Compute the phasing similarity of two assembly graph edges.
// It is computed as the Jaccard similarity of the sets
// of oriented reads internal to each of the two assembly graph edges.
double PhasingData::computePhasingSimilarity(
    AssemblyGraph::EdgeId edgeId0,
    AssemblyGraph::EdgeId edgeId1)
{
    // Access the oriented reads internals to these
    // assembly graph edges.
    const MemoryAsContainer<OrientedReadId> orientedReadIds0 =
        orientedReads[edgeId0];
    const MemoryAsContainer<OrientedReadId> orientedReadIds1 =
        orientedReads[edgeId1];

    // Compute the size of the intersection
    // (number of common oriented reads).
    const uint64_t intersectionSize = shasta::intersectionSize(
        orientedReadIds0.begin(),
        orientedReadIds0.end(),
        orientedReadIds1.begin(),
        orientedReadIds1.end()
    );

    // Compute the size of the union.
    const uint64_t n0 = orientedReadIds0.size();
    const uint64_t n1 = orientedReadIds1.size();
    SHASTA_ASSERT(n0 > 0);
    SHASTA_ASSERT(n1 > 0);
    const uint64_t unionSize = n0 + n1 - intersectionSize;

    // Return the Jaccard similarity.
    return double(intersectionSize) / double(unionSize);
}



uint64_t PhasingData::countCommonInternalOrientedReads(
    AssemblyGraph::EdgeId edgeId0,
    AssemblyGraph::EdgeId edgeId1)
{
    // Access the oriented reads internals to these
    // assembly graph edges
    const MemoryAsContainer<OrientedReadId> orientedReadIds0 =
        orientedReads[edgeId0];
    const MemoryAsContainer<OrientedReadId> orientedReadIds1 =
        orientedReads[edgeId1];

    // Count the number of common oriented reads.
    return intersectionSize(
        orientedReadIds0.begin(),
        orientedReadIds0.end(),
        orientedReadIds1.begin(),
        orientedReadIds1.end()
    );
}



// Find the assembly graph edges related to each assembly graph edge.
// Two assembly graph edges are related if they share at
// least one internal oriented read.
// For now this is single-threaded.
void PhasingData::gatherRelatedAssemblyGraphEdges()
{
    using EdgeId = AssemblyGraph::EdgeId;

    relatedAssemblyGraphEdges.createNew(
        dataName("PhasingRelatedAssemblyGraphEdges"), dataPageSize);

    // Work areas defined here to reduce memory allocation activity.
    vector<EdgeId> edges;
    vector<EdgeId> count;

    // Loop over all assembly graph edges.
    const EdgeId edgeCount = orientedReads.size();
    for(EdgeId e0=0; e0<edgeCount; e0++) {
        edges.clear();

        // Loop over the oriented reads internal to this edge.
        const MemoryAsContainer<OrientedReadId> orientedReadIds0 = orientedReads[e0];
        for(const OrientedReadId orientedReadId: orientedReadIds0) {

            // Loop over the assembly graph edges that this oriented
            // read is internal to.
            const MemoryAsContainer<EdgeId> edges1 = assemblyGraphEdges[orientedReadId.getValue()];
            for(const EdgeId e1: edges1) {
                if(e1 != e0) {
                    edges.push_back(e1);
                }
            }


        }

        // Count how many time we found each of the edges.
        deduplicateAndCount(edges, count);

        // Store.
        relatedAssemblyGraphEdges.appendVector();
        for(uint64_t i=0; i<edges.size(); i++) {
            if(count[i] >= 3) { //   ************** EXPOSE WHEN CODE IS STABLE  *****
                relatedAssemblyGraphEdges.append(make_pair(edges[i], count[i]));
                SHASTA_ASSERT(count[i] == countCommonInternalOrientedReads(e0, edges[i]));
            }
        }
        // cout << e0 << " " << edges.size() << endl;
    }
}



// Function to construct names for binary objects.
string PhasingData::dataName(const string& name) const
{
    if(dataFileNamePrefix.empty()) {
        return "";  // Anonymous;
    } else {
        return dataFileNamePrefix + name;
    }
}



#if 0
void PhasingData::findSimilarPairs(
    size_t threadCount,
    double phasingSimilarityThreshold)
{
    // Store the threshold.
    similarityThreshold = phasingSimilarityThreshold;

    // Allocate space for the similar pairs found by each thread.
    threadSimilarPairs.resize(threadCount);

    // See how many reads we have.
    const uint64_t orientedReadCount = assemblyGraphEdges.size();
    const uint64_t readCount = orientedReadCount / 2;

    // Find the similar pairs with readId0<readId1 in multithreaded code.
    setupLoadBalancing(readCount, 1000);
    runThreads(&PhasingData::findSimilarPairs, threadCount);

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

}



void PhasingData::findSimilarPairs(size_t threadId)
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



double PhasingData::computePhasingSimilarity(
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
    const uint64_t intersectionSize = shasta::intersectionSize(
        edges0.begin(), edges0.end(),
        edges1.begin(), edges1.end()
        );

    const uint64_t n0 = edges0.size();
    const uint64_t n1 = edges1.size();
    SHASTA_ASSERT(n0 > 0);
    SHASTA_ASSERT(n1 > 0);
    return double(intersectionSize) / double(min(n0, n1));

}



// Experimental version that uses turns.
double PhasingData::computePhasingSimilarity(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedparser.add_argument("edgeId0", type=int, help="The id of the firstassembly graph edge.")
ReadId1)
{
    uint64_t concordantCount;
    uint64_t discordantCount;
    tie(concordantCount, discordantCount) =
        countCommonTurns(orientedReadId0, orientedReadId1);

    return
        (double(concordantCount) - double(discordantCount)) /
        (double(concordantCount) + double(discordantCount));
}


// Find similar pairs in which the first ReadId is readId0
// and the second ReadId is readId1>readId0.
// The last argument is a work area, passed in
// to reduce memory allocation activity.
void PhasingData::findSimilarPairs(
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
        }
    }

}



// This is written single-threaded and without attention to performance, for now.
void PhasingData::keepBestSimilarPairs(int maxNeighborCount)
{
    // The indexes of the pairs that each oriented read is involved in.
    const uint64_t orientedReadCount = assemblyGraphEdges.size();
    vector< vector< pair<uint64_t, float> > > table(orientedReadCount);
    for(uint64_t i=0; i<similarPairs.size(); i++) {
        const auto& p = similarPairs[i];
        const auto& q = p.first;
        const float similarity = float(p.second);
        OrientedReadId orientedReadId0(q.readIds[0], 0);
        OrientedReadId orientedReadId1(q.readIds[1], q.isSameStrand ? 0 : 1);
        table[orientedReadId0.getValue()].push_back(make_pair(i, similarity));
        table[orientedReadId1.getValue()].push_back(make_pair(i, similarity));
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        table[orientedReadId0.getValue()].push_back(make_pair(i, similarity));
        table[orientedReadId1.getValue()].push_back(make_pair(i, similarity));
    }


    // For each oriented read, sort by decreasing similarity.
    // Then keep only the best maxNeighborCount.
    vector<bool> keep(similarPairs.size(), false);
    for(uint64_t orientedReadId=0; orientedReadId!=orientedReadCount; orientedReadId++) {
        auto& v = table[orientedReadId];
        if(v.size() > uint64_t(maxNeighborCount)) {
            sort(v.begin(), v.end(), OrderPairsBySecondOnlyGreater<uint64_t, float>());
        }
        for(uint64_t i=0; i<min(uint64_t(maxNeighborCount), uint64_t(v.size())); i++) {
            keep[v[i].first] = true;
        }
    }


    // Keep only the pairs we flagged.
    uint64_t i = 0;
    uint64_t j = 0;
    for(; i<similarPairs.size(); i++) {
        if(keep[i]) {
            similarPairs[j++] = similarPairs[i];
        }
    }
    similarPairs.resize(j);
    cout << "Kept " << similarPairs.size() << " best similar pairs." << endl;
}



// Count concordant and discordant turns for two oriented reads.
// Returns pair(concordantCount, discordantCount).
pair<uint64_t, uint64_t> PhasingData::countCommonTurns(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1)
{
    uint64_t concordantCount = 0;
    uint64_t discordantCount = 0;

    // Access the turns of these two oriented reads.
    const MemoryAsContainer<Turn> turns0 =
        turns[orientedReadId0.getValue()];
    const MemoryAsContainer<Turn> turns1 =
        turns[orientedReadId1.getValue()];

    // Iterate over common turns.
    // The turns are sorted by vertex id v.
    auto begin0 = turns0.begin();
    auto begin1 = turns1.begin();
    auto end0   = turns0.end();
    auto end1   = turns1.end();
    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        const Turn& turn0 = *it0;
        const Turn& turn1 = *it1;
        if(turn0.v < turn1.v) {
            ++it0;
        } else if(turn1.v < turn0.v) {
            ++it1;
        } else {
            if(turn0 == turn1) {
                ++concordantCount;
            } else {
                ++discordantCount;
            }
            ++it0;
            ++it1;
        }
    }

     return make_pair(concordantCount, discordantCount);
}



// Write all the pairs in graphviz format.
void PhasingData::writeGraphviz()
{
    ofstream graphout("PhasingData.dot");
    graphout << "graph G {" << endl;
    graphout  << "tooltip = \" \";\n";

    // Vertices.
    const ReadId orientedReadCount = ReadId(assemblyGraphEdges.size());
    const ReadId readCount = orientedReadCount / 2;
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            graphout << orientedReadId.getValue();
            graphout << "[";
            graphout << "tooltip=\"" << orientedReadId << "\"";
            graphout << "];\n";
        }
    }


    // Edges.
    for(const auto& p: similarPairs) {
        const auto& q = p.first;
        OrientedReadId orientedReadId0(q.readIds[0], 0);
        OrientedReadId orientedReadId1(q.readIds[1], q.isSameStrand ? 0 : 1);
        graphout << orientedReadId0.getValue();
        graphout << "--";
        graphout << orientedReadId1.getValue();
        graphout << "[";
        graphout << "tooltip=\"" << orientedReadId0 << " " << orientedReadId1 << " " <<
            p.second << "\"";
        graphout << "];\n";
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        graphout << orientedReadId0.getValue();
        graphout << "--";
        graphout << orientedReadId1.getValue();
        graphout << "[";
        graphout << "tooltip=\"" << orientedReadId0 << " " << orientedReadId1 << " " <<
            p.second << "\"";
        graphout << "];\n";
    }
    graphout << "}" << endl;

}
#endif

