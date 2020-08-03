// This file contains code for ReadGraph.creationMethod 2.

// Shasta.
#include "Assembler.hpp"
#include "MiniAssemblyMarkerGraph.hpp"
#include "orderPairs.hpp"
using namespace shasta;


void Assembler::createReadGraph2(size_t threadCount)
{
    // Parameters that control this function.
    // expose when code stabilizes.

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Mark all alignments as not to be kept.
    createReadGraph2Data.keepAlignment.resize(alignmentData.size());
    fill(
        createReadGraph2Data.keepAlignment.begin(),
        createReadGraph2Data.keepAlignment.end(),
        false);

    // Parallel loop over reads. For each read, flag the alignments we want to discard.
    const uint64_t batchSize = 100;
    setupLoadBalancing(reads.readCount(), batchSize);
    runThreads(&Assembler::createReadGraph2ThreadFunction, threadCount);;

    // Create the read graph using the alignments we selected.
    const size_t keepCount =
        count(
            createReadGraph2Data.keepAlignment.begin(),
            createReadGraph2Data.keepAlignment.end(),
            true);
    cout << "Keeping " << keepCount << " alignments of " << createReadGraph2Data.keepAlignment.size() << endl;
    createReadGraphUsingSelectedAlignments(createReadGraph2Data.keepAlignment);

}



void Assembler::createReadGraph2ThreadFunction(size_t threadId)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {
            createReadGraph2LowLevel(readId);
        }

    }
}



// For one read readId0, flag the alignments we want to discard.
// Like in analyzeAlignments3 (see AssemblerAlignments.cpp),
// we do a mini-assembly starting from readId0 on strand 0, that is
// orientedReadId0 = OrientedReadRead(readId0, 0).
// The mini-assembly uses orientedReadId0 plus
// the aligned portions of oriented reads for which
// we have an alignment with orientedReadId0.
// We then suppress alignments with oriented reads that tend to end up
// in different bubble branches from orientedReadId0.
void Assembler::createReadGraph2LowLevel(ReadId readId0)
{

    // Parameters controlling this function.
    // Expose when code stabilizes.
    const uint64_t minTotalEdgeCoverage = 4;
    const uint64_t minPerStrandEdgeCoverage = 1;
    const uint64_t neighborCount = 8;

    // Work with this read on strand 0.
    const OrientedReadId orientedReadId0(readId0, 0);

    // Get the alignments involving this oriented read.
    // This returns a vector alignments with swaps and/or
    // reverse complementing already done, as necessary.
    vector<StoredAlignmentInformation> alignments;
    getStoredAlignments(orientedReadId0, alignments);

    // Check that all alignments are strictly increasing.
    for(const auto& p: alignments) {
        p.alignment.checkStrictlyIncreasing();
    }



    // We will do a small assembly for the marker sequence of this oriented read
    // plus the aligned portions of the marker sequences of aligned reads.
    // Gather these sequences.
    // The marker sequence for this oriented read is stored
    // at the last position of this vector.
    using Sequence = vector<KmerId>;
    using SequenceId = uint64_t;
    vector<Sequence> sequences(alignments.size() + 1);
    vector<OrientedReadId> orientedReadIds(sequences.size());
    vector<uint32_t> firstOrdinals(sequences.size());
    vector<uint32_t> lastOrdinals(sequences.size());
    for(SequenceId sequenceId=0; sequenceId<alignments.size(); sequenceId++) {
        Sequence& sequence = sequences[sequenceId];
        const OrientedReadId orientedReadId1 = alignments[sequenceId].orientedReadId;
        orientedReadIds[sequenceId] = orientedReadId1;
        const span<CompressedMarker> markers1 = markers[orientedReadId1.getValue()];
        const Alignment& alignment = alignments[sequenceId].alignment;
        const uint32_t first1 = alignment.ordinals.front()[1];
        firstOrdinals[sequenceId] = first1;
        const uint32_t last1 = alignment.ordinals.back()[1];
        lastOrdinals[sequenceId] = last1;
        sequence.resize(last1 + 1 - first1);
        for(uint64_t i=0; i<sequence.size(); i++) {
            sequence[i] = markers1[first1 + i].kmerId;
        }
    }

    // Add the sequence of the oriented read we started from.
    Sequence& sequence0 = sequences.back();
    const SequenceId sequenceId0 = sequences.size() - 1;
    orientedReadIds.back() = orientedReadId0;
    const span<CompressedMarker> markers0 = markers[orientedReadId0.getValue()];
    const uint64_t markerCount0 = markers0.size();
    firstOrdinals.back() = 0;
    lastOrdinals.back() = uint32_t(markers0.size() - 1);
    sequence0.resize(markerCount0);
    for(uint32_t ordinal=0; ordinal!=markerCount0; ordinal++) {
        sequence0[ordinal] = markers0[ordinal].kmerId;
    }



    // Create a marker graph of these sequences.
    // Use as SequenceId the index into the sequences vector.
    using Graph = MiniAssemblyMarkerGraph;
    Graph graph(orientedReadIds);
    for(SequenceId sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        graph.addSequence(sequenceId, sequences[sequenceId]);
    }
    graph.doneAddingSequences();



    // Merge pairs of aligned markers.
    vector< pair<uint64_t, uint64_t> > v;
    for(SequenceId sequenceId1=0; sequenceId1<alignments.size(); sequenceId1++) {
        const Alignment& alignment = alignments[sequenceId1].alignment;
        v.clear();
        for(const auto& ordinals: alignment.ordinals) {
            // Merge ordinals relative to the start of the portion of
            // each sequence used in the mini-assembly.
            v.push_back({
                ordinals[0] - firstOrdinals[sequenceId0],
                ordinals[1] - firstOrdinals[sequenceId1]});
            SHASTA_ASSERT(
                markers[orientedReadIds[sequenceId0].getValue()][ordinals[0]].kmerId ==
                markers[orientedReadIds[sequenceId1].getValue()][ordinals[1]].kmerId
                );
        }
        graph.merge(sequenceId0, sequenceId1, v);
    }



    // We also need to merge vertices using alignments between the oriented reads
    // aligned with orientedReadId0.
    // Just for this portion of the code, take orientedReadId0 out of the orientedReadIds
    // vector.
    orientedReadIds.resize(alignments.size());
    for(SequenceId sequenceId1=0; sequenceId1<alignments.size(); sequenceId1++) {
        const OrientedReadId orientedReadId1 = orientedReadIds[sequenceId1];

        // Get alignments between orientedReadId1 and the other oriented reads in
        // orientedReadIds.
        getStoredAlignments(orientedReadId1, orientedReadIds, alignments);

        // Loop over the alignments we got.
        for(const auto& storedAlignment: alignments) {
            const OrientedReadId orientedReadId2 = storedAlignment.orientedReadId;

            // Look up the corresponding SequenceId.
            const auto it = std::lower_bound(orientedReadIds.begin(), orientedReadIds.end(),
                orientedReadId2);
            SHASTA_ASSERT(it != orientedReadIds.end());
            const SequenceId sequenceId2 = it - orientedReadIds.begin();

            // Merge vertices.
            const Alignment& alignment = storedAlignment.alignment;
            v.clear();
            for(const auto& ordinals: alignment.ordinals) {
                if(ordinals[0] < firstOrdinals[sequenceId1]) {
                    continue;
                }
                if(ordinals[1] < firstOrdinals[sequenceId2]) {
                    continue;
                }
                if(ordinals[0] > lastOrdinals[sequenceId1]) {
                    continue;
                }
                if(ordinals[1] > lastOrdinals[sequenceId2]) {
                    continue;
                }
                v.push_back({
                    ordinals[0] - firstOrdinals[sequenceId1],
                    ordinals[1] - firstOrdinals[sequenceId2]});
                SHASTA_ASSERT(
                    markers[orientedReadIds[sequenceId1].getValue()][ordinals[0]].kmerId ==
                    markers[orientedReadIds[sequenceId2].getValue()][ordinals[1]].kmerId
                    );
            }
            graph.merge(sequenceId1, sequenceId2, v);
        }
    }
    // Add orientedReadId0 back to our list.
    orientedReadIds.push_back(orientedReadId0);



    // Finish creation of the marker graph.
    graph.doneMerging();
    graph.removeSelfEdges();
    graph.removeLowCoverageEdges(minTotalEdgeCoverage, minPerStrandEdgeCoverage);
    graph.removeIsolatedVertices();
    graph.findBubbles();


    // For each of the aligned oriented reads, find the number of times
    // it is on the same or different branch from orientedReadId0.
    // Store tuples(sequenceId, sameBranchCount, differentBranchCount.
    using Tuple = tuple<SequenceId, uint64_t, uint64_t>;
    vector<Tuple> branchCounts;
    for(SequenceId sequenceId1=0; sequenceId1<orientedReadIds.size()-1; sequenceId1++) {
        uint64_t sameCount = 0;
        uint64_t differentCount = 0;
        for(const Graph::Bubble& bubble: graph.bubbles) {
            const int64_t branchId0 = bubble.branchTable[sequenceId0];
            if(branchId0 < 0) {
                continue;
            }
            const int64_t branchId1 = bubble.branchTable[sequenceId1];
            if(branchId1 < 0) {
                continue;
            }
            if(branchId0 == branchId1) {
                ++sameCount;
            } else {
                ++differentCount;
            }
        }
        branchCounts.push_back(make_tuple(sequenceId1, sameCount, differentCount));
    }

    // Sort so the "best" alignments come first.
    sort(branchCounts.begin(), branchCounts.end(),
        [](const Tuple& t0, const Tuple& t1)
        {
            const uint64_t different0 = get<2>(t0);
            const uint64_t different1 = get<2>(t1);
            const uint64_t same0 = get<1>(t0);
            const uint64_t same1 = get<1>(t1);
            const int64_t delta0 = int64_t(different0) - int64_t(same0);
            const int64_t delta1 = int64_t(different1) - int64_t(same1);
            return delta0 < delta1;
        });

    // Only keep up to neighborCount.
    if(branchCounts.size() > neighborCount) {
        branchCounts.resize(neighborCount);
    }


    // Flag the alignments we want to keep.
    const auto alignmentTable0 = alignmentTable[orientedReadId0.getValue()];
    for(const auto& t: branchCounts) {
        const SequenceId sequenceId = get<0>(t);
        const uint64_t alignmentId = alignmentTable0[sequenceId];
        createReadGraph2Data.keepAlignment[alignmentId] = true;
    }


}


