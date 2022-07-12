// Shasta.
#include "Assembler.hpp"
#include "AssemblyGraph.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
#include "seqan.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Standard library.
#include "fstream.hpp"



// This use PseudoPaths to decide which alignments
// should be included in the read graph.
// See Assembler::alignPseudoPaths in AssemblerAnalyzePaths.cpp.
void Assembler::createReadGraphUsingPseudoPaths(
    int64_t matchScore,
    int64_t mismatchScore,
    int64_t gapScore,
    double mismatchSquareFactor,
    double minScore,
    uint64_t maxAlignmentCount,
    size_t threadCount)
{
    using SegmentId = AssemblyGraphEdgeId;
    const bool debug = false;

    // Store the parameters so all threads can see them.
    createReadGraphUsingPseudoPathsData.matchScore = matchScore;
    createReadGraphUsingPseudoPathsData.mismatchScore = mismatchScore;
    createReadGraphUsingPseudoPathsData.gapScore = gapScore;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Compute the pseudo-path of each oriented read.
    // This vector is indexed by OrientedReadId::getValue().
    const uint64_t readCount = reads->readCount();
    auto& pseudoPathSegments = createReadGraphUsingPseudoPathsData.pseudoPaths;
    pseudoPathSegments.resize(2*readCount);
    size_t batchSize = 1000;
    cout << timestamp << "Computing pseudopaths for " << readCount << " reads." << endl;
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::createReadGraphUsingPseudoPathsThreadFunction1, threadCount);



    // Write a csv file with the pseudo-path of each oriented read.
    if(debug) {
        ofstream csv("PseudoPaths.csv");
        for(ReadId readId=0; readId<reads->readCount(); readId++) {
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);
                csv << orientedReadId << ",";

                const vector<SegmentId>& segments = pseudoPathSegments[orientedReadId.getValue()];
                for(const SegmentId segment: segments) {
                    csv << segment << ",";
                }
                csv << "\n";
            }
        }
    }




    // For each alignment we have, align the pseudo-paths
    // of the two oriented reads, putting the first read on strand 0.
    cout << timestamp << "Computing pseudopath alignments for " <<
        alignmentData.size() << " alignments." << endl;
    createReadGraphUsingPseudoPathsData.alignmentInfos.resize(alignmentData.size());
    setupLoadBalancing(alignmentData.size(), batchSize);
    runThreads(&Assembler::createReadGraphUsingPseudoPathsThreadFunction2, threadCount);


    // Write out this information, by read.
    if(debug) {
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
                const auto& info = createReadGraphUsingPseudoPathsData.alignmentInfos[alignmentId];
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
            const auto& info = createReadGraphUsingPseudoPathsData.alignmentInfos[alignmentId];
            const double score = double(info.strongMatchCount) -
                mismatchSquareFactor * double(info.mismatchCount*info.mismatchCount);
            if(score > minScore) {
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
    cout << timestamp << "Keeping " << keepCount << " alignments of " << keepAlignment.size() << endl;
    readGraph.remove();
    createReadGraphUsingSelectedAlignments(keepAlignment);
}



// Thread function used to compute pseudoPaths.
void Assembler::createReadGraphUsingPseudoPathsThreadFunction1(size_t threadId)
{
    vector<MarkerGraph::EdgeId> path;
    vector< pair<uint32_t, uint32_t> > pathOrdinals;
    PseudoPath pseudoPath;
    using SegmentId = AssemblyGraphEdgeId;
    vector< vector<SegmentId> >& pseudoPaths =
        createReadGraphUsingPseudoPathsData.pseudoPaths;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads in this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);
                computePseudoPath(orientedReadId, path, pathOrdinals, pseudoPath);
                getPseudoPathSegments(pseudoPath, pseudoPaths[orientedReadId.getValue()]);
            }
        }
    }
}



// Thread functions used to align pseudopaths.
void Assembler::createReadGraphUsingPseudoPathsThreadFunction2(size_t threadId)
{

    // Extract parameters.
    const auto matchScore = createReadGraphUsingPseudoPathsData.matchScore;
    const auto mismatchScore = createReadGraphUsingPseudoPathsData.mismatchScore;
    const auto gapScore = createReadGraphUsingPseudoPathsData.gapScore;

    // Access global objects.
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    using SegmentId = AssemblyGraphEdgeId;
    const vector< vector<SegmentId> >& pseudoPathSegments =
        createReadGraphUsingPseudoPathsData.pseudoPaths;
    auto& infos = createReadGraphUsingPseudoPathsData.alignmentInfos;

    vector< pair<bool, bool> > alignment;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all alignments in this batch.
        for(uint64_t alignmentId=begin; alignmentId!=end; alignmentId++) {
            const AlignmentData& ad = alignmentData[alignmentId];
            auto& info = infos[alignmentId];

            // Gather the two oriented reads.
            // Put the first one on strand 0.
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
                info.alignedMarkerCount = ad.info.markerCount;
                info.weakMatchCount = 0;
                info.strongMatchCount = 0;
                info.mismatchCount = 0;
                continue;
            }

            // Align them.
            shasta::seqanAlign(
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
            info.alignedMarkerCount = ad.info.markerCount;
            info.weakMatchCount = weakMatchCount;
            info.strongMatchCount = strongMatchCount;
            info.mismatchCount = mismatchCount;
        }
    }
}
