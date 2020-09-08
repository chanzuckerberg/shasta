
// Shasta.
#include "Assembler.hpp"
#include "MiniAssemblyMarkerGraph.hpp"
#include "orderPairs.hpp"
#include "seqan.hpp"
using namespace shasta;


# if 0
// Version that uses mini-assemblies to select alignments
// to be included in the read graph
// (see analyzeAlignments3 in AssemblerAlignments.cpp).
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
    setupLoadBalancing(reads->readCount(), batchSize);
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
#endif


// This versions use PseudoPaths to decide which alignments
// should be included in the read graph.
// See Assembler::alignPseudoPaths in AssemblerAnalyzePaths.cpp.
// This is a quick and dirty single threaded implementation for testing.
// If successful, a multithreaded version will be needed.
void Assembler::createReadGraphUsingPseudoPaths(
    int64_t matchScore,
    int64_t mismatchScore,
    int64_t gapScore,
    double mismatchSquareFactor,
    double minScore,
    uint64_t maxAlignmentCount,
    size_t threadCount)
{
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;



    // Compute the pseudo-path of each oriented read.
    // This vector is indexed by OrientedReadId::getValue().
    vector<MarkerGraph::EdgeId> path;
    vector< pair<uint32_t, uint32_t> > pathOrdinals;
    PseudoPath pseudoPath;
    using SegmentId = AssemblyGraph::EdgeId;
    const uint64_t readCount = reads->readCount();
    vector< vector<SegmentId> > pseudoPathSegments(2*readCount);
    cout << timestamp << "Computing pseudo-paths for " << readCount << " reads." << endl;
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            computePseudoPath(orientedReadId,
                path, pathOrdinals, pseudoPath);
            getPseudoPathSegments(pseudoPath,
                pseudoPathSegments[orientedReadId.getValue()]);
        }
    }



    // Vector to store the information we need for each alignment.
    vector<CreateReadGraphsingPseudoPathsAlignmentData> infos(alignmentData.size());



    // For each alignment we have, align the pseudo-paths
    // of the two oriented reads, putting the first read on strand 0.
    vector< pair<bool, bool> > alignment;
    cout << timestamp << "Computing pseudo-path alignments for " <<
        alignmentData.size() << " alignments." << endl;
    for(uint64_t alignmentId=0; alignmentId<alignmentData.size(); alignmentId++) {
        const AlignmentData& ad = alignmentData[alignmentId];

        // Gather the two oriented reads.
        const Strand strand0 = 0;
        const Strand strand1 = ad.isSameStrand ? 0 : 1;
        const OrientedReadId orientedReadId0(ad.readIds[0], strand0);
        const OrientedReadId orientedReadId1(ad.readIds[1], strand1);

        // Find the corresponding pseudo-paths.
        const vector<SegmentId>& pseudoPathSegments0 =
            pseudoPathSegments[orientedReadId0.getValue()];
        const vector<SegmentId>& pseudoPathSegments1 =
            pseudoPathSegments[orientedReadId1.getValue()];

        // Skip pathological case.
        if(pseudoPathSegments0.empty() or pseudoPathSegments1.empty()) {
            continue;
        }

        // Align them.
        /* const uint64_t alignmentScore =*/ shasta::seqanAlign(
            pseudoPathSegments0.begin(), pseudoPathSegments0.end(),
            pseudoPathSegments1.begin(), pseudoPathSegments1.end(),
            matchScore,
            mismatchScore,
            gapScore,
            true, true,
            alignment);

        // Analyze the alignment of the two pseudo-paths.
        uint64_t position0 = 0;
        uint64_t position1 = 0;
        uint64_t weakMatchCount =0;
        uint64_t strongMatchCount =0;
        uint64_t mismatchCount =0;
        uint64_t gapCount =0;
        uint64_t leftUnalignedCount =0;
        uint64_t rightUnalignedCount =0;
        for(const auto& p: alignment) {
            if(p.first and p.second) {
                if(pseudoPathSegments0[position0] != pseudoPathSegments1[position1]) {
                    ++mismatchCount;
                } else {
                    // Match. Figure out if it is a weak or strong match.
                    const SegmentId segmentId = pseudoPathSegments0[position0];
                    const AssemblyGraph::Edge& edge = assemblyGraph.edges[segmentId];
                    const AssemblyGraph::VertexId v0 = edge.source;
                    const AssemblyGraph::VertexId v1 = edge.target;
                    const auto out0 = assemblyGraph.outDegree(v0);
                    const auto in1 = assemblyGraph.inDegree(v1);
                    if(out0==1 and in1==1) {    // CONSIDER DOING OR INSTEAD?
                        ++weakMatchCount;
                    } else {
                        ++strongMatchCount;
                    }
                }
            } else if(position0 == 0 or position1==0) {
                ++leftUnalignedCount;
            } else if(
                position0 == pseudoPathSegments0.size() or
                position1 == pseudoPathSegments1.size()) {
                ++rightUnalignedCount;
            } else {
                ++gapCount;
            }

            if(p.first) {
                ++position0;
            }
            if(p.second) {
                ++position1;
            }
        }
        SHASTA_ASSERT(position0 == pseudoPathSegments0.size());
        SHASTA_ASSERT(position1 == pseudoPathSegments1.size());
        SHASTA_ASSERT(
            weakMatchCount + strongMatchCount + mismatchCount +
            gapCount + leftUnalignedCount + rightUnalignedCount ==
            alignment.size());

        // Store the information for this alignment.
        auto& info = infos[alignmentId];
        info.alignedMarkerCount = ad.info.markerCount;
        info.weakMatchCount = weakMatchCount;
        info.strongMatchCount = strongMatchCount;
        info.mismatchCount = mismatchCount;
    }



    // Write out this information, by read.
    ofstream csv("CreateReadGraph2.csv");
    csv << "ReadId,AlignmentId,ReadId0,ReadId1,SameStrand,AlignedMarkerCount,"
        "WeakMatchCount,StrongMatchCount,MismatchCount,Score\n";
    for(ReadId readId=0; readId<readCount; readId++) {

        // Put it on strand 0.
        const OrientedReadId orientedReadId(readId, 0);

        // Get the alignments it is involved in.
        const span<uint32_t> alignmentIds = alignmentTable[orientedReadId.getValue()];

        // Loop over those alignments.
        for(const uint32_t alignmentId: alignmentIds) {
            const AlignmentData& ad = alignmentData[alignmentId];
            const auto& info = infos[alignmentId];
            const double score = double(info.strongMatchCount) -
                mismatchSquareFactor * double(info.mismatchCount*info.mismatchCount);
            csv << readId << ",";
            csv << alignmentId << ",";
            csv << ad.readIds[0] << ",";
            csv << ad.readIds[1] << ",";
            csv << (ad.isSameStrand ? "Yes" : "No") << ",";
            csv << ad.info.markerCount << ",";
            csv << info.weakMatchCount << ",";
            csv << info.strongMatchCount << ",";
            csv << info.mismatchCount << ",";
            csv << score << "\n";
        }
    }



    // For each read, flag the alignments we want to keep.
    vector<bool> keepAlignment(alignmentData.size(), false);
    uint64_t tooFewCount = 0;
    for(ReadId readId=0; readId<readCount; readId++) {

        // Put it on strand 0.
        const OrientedReadId orientedReadId(readId, 0);

        // Get the alignments it is involved in.
        const span<uint32_t> alignmentIds = alignmentTable[orientedReadId.getValue()];

        // Sort them by score = segmentMatchCount - mismatchSquareFactor * segmentMismatchCount^2
        vector< pair<double, uint32_t> > table; // pair(score, alignmentId)
        for(const uint32_t alignmentId: alignmentIds) {
            const auto& info = infos[alignmentId];
            const double score = double(info.strongMatchCount) -
                mismatchSquareFactor * double(info.mismatchCount*info.mismatchCount);
            if(score >= minScore) {
                table.push_back(make_pair(score, alignmentId));
            }
        }
        sort(table.begin(), table.end(), OrderPairsByFirstOnlyGreater<double, uint32_t>());

        // Keep the best maxAlignmentCount.
        if(table.size() > maxAlignmentCount) {
            table.resize(maxAlignmentCount);
        }

        if(table.size() < maxAlignmentCount) {
            ++tooFewCount;
        }

        // Mark them to keep.
        for(const auto& p: table) {
            keepAlignment[p.second] = true;
        }
    }
    cout << "Too few: " << tooFewCount << endl;


    // Create the read graph using the alignments we selected.
    const size_t keepCount = count(keepAlignment.begin(), keepAlignment.end(), true);
    cout << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;
    readGraph.remove();
    createReadGraphUsingSelectedAlignments(keepAlignment);
}
