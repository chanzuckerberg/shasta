// This file contains code for ReadGraph.creationMethod 2.

// Shasta.
#include "Assembler.hpp"
#include "DeBruijnGraph.hpp"
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

    // Mark all alignments as to be kept.
    createReadGraph2Data.keepAlignment.resize(alignmentData.size());
    fill(
        createReadGraph2Data.keepAlignment.begin(),
        createReadGraph2Data.keepAlignment.end(),
        true);

    // Parallel loop over reads. For each read, flag the alignments we want to discard.
    const uint64_t batchSize = 100;
    setupLoadBalancing(readCount(), batchSize);
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
// Like in analyzeAlignments2 (see AssemblerAlignments.cpp),
// we do a mini-assembly startting from readId0 on strand 0, that is
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
    const uint64_t minTotalCoverage = 5;
    const uint64_t minSameStrandCoverage = 2;
    const uint64_t minOppositeStrandCoverage = 2;
    const uint64_t minDifferentBranchCount = 2;



    // Get the alignments of this oriented read, with the proper orientation,
    // and with this oriented read as the first oriented read in the alignment.
    const Strand strand0 = 0;
    const OrientedReadId orientedReadId0(readId0, strand0);
    const vector< pair<OrientedReadId, AlignmentInfo> > alignments =
        findOrientedAlignments(orientedReadId0);



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
    for(SequenceId sequenceId=0; sequenceId<alignments.size(); sequenceId++) {
        Sequence& sequence = sequences[sequenceId];
        const OrientedReadId orientedReadId1 = alignments[sequenceId].first;
        orientedReadIds[sequenceId] = orientedReadId1;
        const span<CompressedMarker> markers1 = markers[orientedReadId1.getValue()];
        const AlignmentInfo& alignmentInfo = alignments[sequenceId].second;
        const uint32_t first1 = alignmentInfo.data[1].firstOrdinal;
        firstOrdinals[sequenceId] = first1;
        const uint32_t last1 = alignmentInfo.data[1].lastOrdinal;
        sequence.resize(last1 + 1 - first1);
        for(uint64_t i=0; i<sequence.size(); i++) {
            sequence[i] = markers1[first1 + i].kmerId;
        }
    }
    Sequence& sequence0 = sequences.back();
    orientedReadIds.back() = orientedReadId0;
    firstOrdinals.back() = 0;
    const span<CompressedMarker> markers0 = markers[orientedReadId0.getValue()];
    const uint64_t markerCount0 = markers0.size();
    sequence0.resize(markerCount0);
    for(uint32_t ordinal=0; ordinal!=markerCount0; ordinal++) {
        sequence0[ordinal] = markers0[ordinal].kmerId;
    }



    // Create the De Bruijn graph.
    // Use as SequenceId the index into the above vector of sequences.
    using Graph = DeBruijnGraph<KmerId, 3, uint64_t>;
    using vertex_descriptor = Graph::vertex_descriptor;
    Graph graph;
    for(SequenceId sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        graph.addSequence(sequenceId, sequences[sequenceId]);
    }
    graph.removeAmbiguousVertices();



    // Remove low coverage vertices.
    vector<vertex_descriptor> verticesTobeRemoved;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {

        if(graph[v].occurrences.size() < minTotalCoverage) {

            // Total coverage is too low.
            verticesTobeRemoved.push_back(v);

        } else {

            // Total coverage is sufficient. Check coverage per strand.
            array<uint64_t, 2>coveragePerStrand = {0, 0};
            for(const auto& p: graph[v].occurrences) {
                const SequenceId sequenceId = p.first;
                const OrientedReadId orientedReadId = orientedReadIds[sequenceId];
                ++coveragePerStrand[orientedReadId.getStrand()];
            }

            if(
                (coveragePerStrand[strand0] < minSameStrandCoverage) or
                (coveragePerStrand[1 - strand0] < minOppositeStrandCoverage)) {
                verticesTobeRemoved.push_back(v);
            }
        }
    }
    for(const vertex_descriptor v: verticesTobeRemoved) {
        clear_vertex(v, graph);
        remove_vertex(v, graph);
    }



    // Create edges of the De Bruijn graph.
    graph.createEdges();

    // Find sets of incompatible vertices.
    std::set< std::set<Graph::vertex_descriptor> > incompatibleVertexSets;
    graph.findIncompatibleVertexSets(incompatibleVertexSets);



    // For each set of incompatible vertices,
    // construct a signature vector that tells us which of the incompatible vertices
    // each reads appears in, if any.
    // >=0: Gives the index of the vertex (in the incompatible set) in which the read appears.
    // -1 = Read does not appear in the incompatible vertex set.
    vector< vector<int64_t> > signatures(
        incompatibleVertexSets.size(), vector<int64_t>(sequences.size(), -1));

    uint64_t i = 0;
    for(const auto& incompatibleVertexSet : incompatibleVertexSets) {

        // Copy the set to a vector for ease in manipulating.
        vector<Graph::vertex_descriptor> incompatibleVertexVector(incompatibleVertexSet.size());
        copy(incompatibleVertexSet.begin(), incompatibleVertexSet.end(), incompatibleVertexVector.begin());

        // Find out in which branch each sequence appears.
        // -1 = does not appear.
        vector<int64_t>& signature = signatures[i];
        for(uint64_t branch=0; branch<incompatibleVertexVector.size(); branch++) {
            for(const auto& p: graph[incompatibleVertexVector[branch]].occurrences) {
                const SequenceId sequenceId = p.first;
                const int64_t oldValue = signature[sequenceId];
                if(oldValue == -1) {
                    signature[sequenceId] = branch; // This is the first time we see it.
                } else {
                    // We have already seen it.
                    SHASTA_ASSERT(0);   // findIncompatibleVertexSets should never generate this.
                }
            }
        }

        ++i;
    }



    // For each of the aligned reads, count how many times it appears
    // on a different branch than our orientedReadId0.
    const auto alignmentTable0 = alignmentTable[orientedReadId0.getValue()];
    for(SequenceId sequenceId=0; sequenceId<sequences.size()-1; sequenceId++) {
        uint64_t differentBranchCount = 0;
        for(const vector<int64_t>& signature: signatures) {
            const int64_t signature0 = signature.back();
            if(signature0 == -1) {
                continue;
            }
            const int64_t signature1 = signature[sequenceId];
            if(signature1 == -1) {
                continue;
            }
            if(signature0 != signature1) {
                ++differentBranchCount;
            }
        }

        // If on different branches enough times, discard the alignment
        // between orientedReadId0 and this oriented read.
        if(differentBranchCount >= minDifferentBranchCount) {
            const uint64_t alignmentId = alignmentTable0[sequenceId];
            createReadGraph2Data.keepAlignment[alignmentId] = false;
        }
    }

}


