// Nanopore2.
#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Standard libraries.
#include "chrono.hpp"
#include "tuple.hpp"



// Compute a marker alignment of two oriented reads.
void Assembler::alignOrientedReads(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    size_t maxSkip, // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer
)
{
    alignOrientedReads(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1),
        maxSkip, maxVertexCountPerKmer
        );
}
void Assembler::alignOrientedReads(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    size_t maxSkip, // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer
)
{
    checkReadsAreOpen();
    checkReadNamesAreOpen();
    checkMarkersAreOpen();


    // Get the markers sorted by kmerId.
    vector<MarkerWithOrdinal> markers0SortedByKmerId;
    vector<MarkerWithOrdinal> markers1SortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markers0SortedByKmerId);
    getMarkersSortedByKmerId(orientedReadId1, markers1SortedByKmerId);

    // Call the lower level function.
    AlignmentGraph graph;
    Alignment alignment;
    const bool debug = true;
    alignOrientedReads(
        markers0SortedByKmerId,
        markers1SortedByKmerId,
        maxSkip, maxVertexCountPerKmer, debug, graph, alignment);

    // Compute the AlignmentInfo.
    const AlignmentInfo alignmentInfo(alignment);
    uint32_t leftTrim;
    uint32_t rightTrim;
    tie(leftTrim, rightTrim) = computeTrim(
        orientedReadId0,
        orientedReadId1,
        alignmentInfo);
    cout << orientedReadId0 << " has " << reads[orientedReadId0.getReadId()].baseCount;
    cout << " bases and " << markers0SortedByKmerId.size() << " markers." << endl;
    cout << orientedReadId1 << " has " << reads[orientedReadId1.getReadId()].baseCount;
    cout << " bases and " << markers1SortedByKmerId.size() << " markers." << endl;
    cout << "The alignment has " << alignmentInfo.markerCount;
    cout << " markers. Left trim " << leftTrim;
    cout << " bases, right trim " << rightTrim << " bases." << endl;

    // For convenience, also write the two oriented reads.
    ofstream fasta("AlignedOrientedReads.fasta");
    writeOrientedRead(orientedReadId0, fasta);
    writeOrientedRead(orientedReadId1, fasta);
}



// This lower level version takes as input vectors of
// markers already sorted by kmerId.
void Assembler::alignOrientedReads(
    const vector<MarkerWithOrdinal>& markers0SortedByKmerId,
    const vector<MarkerWithOrdinal>& markers1SortedByKmerId,
    size_t maxSkip,             // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer
)
{
    // Compute the alignment.
    AlignmentGraph graph;
    Alignment alignment;
    const bool debug = true;
    alignOrientedReads(
        markers0SortedByKmerId,
        markers1SortedByKmerId,
        maxSkip, maxVertexCountPerKmer, debug, graph, alignment);
}



void Assembler::alignOrientedReads(
    const vector<MarkerWithOrdinal>& markers0SortedByKmerId,
    const vector<MarkerWithOrdinal>& markers1SortedByKmerId,
    size_t maxSkip,             // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer,
    bool debug,
    AlignmentGraph& graph,
    Alignment& alignment
)
{
    align(markers0SortedByKmerId, markers1SortedByKmerId,
        maxSkip, maxVertexCountPerKmer, graph, debug, alignment);
}



// Compute marker alignments of an oriented read with all reads
// for which we have an Overlap.
void Assembler::alignOverlappingOrientedReads(
    ReadId readId, Strand strand,
    size_t maxSkip,                 // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer,
    size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
    size_t maxTrim                  // Maximum trim allowed in an alignment.
)
{
    alignOverlappingOrientedReads(
        OrientedReadId(readId, strand),
        maxSkip, maxVertexCountPerKmer, minAlignedMarkerCount, maxTrim);
}



void Assembler::alignOverlappingOrientedReads(
    OrientedReadId orientedReadId0,
    size_t maxSkip,                 // Maximum ordinal skip allowed.
    size_t maxVertexCountPerKmer,
    size_t minAlignedMarkerCount,   // Minimum number of markers in an alignment.
    size_t maxTrim                  // Maximum trim allowed in an alignment.
)
{
    // Check that we have what we need.
    checkReadsAreOpen();
    checkMarkersAreOpen();
    checkOverlapsAreOpen();

    // Get the markers for orientedReadId0.
    vector<MarkerWithOrdinal> markers0SortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markers0SortedByKmerId);

    // Loop over all overlaps involving this oriented read.
    vector<MarkerWithOrdinal> markers1SortedByKmerId;
    size_t goodAlignmentCount = 0;
    for(const uint64_t i: overlapTable[orientedReadId0.getValue()]) {
        const Overlap& overlap = overlaps[i];

        // Get the other oriented read involved in this overlap.
        const OrientedReadId orientedReadId1 = overlap.getOther(orientedReadId0);

        // Get the markers for orientedReadId1.
        getMarkersSortedByKmerId(orientedReadId1, markers1SortedByKmerId);

        // Compute the alignment.
        AlignmentGraph graph;
        Alignment alignment;
        const bool debug = false;
        alignOrientedReads(
            markers0SortedByKmerId,
            markers1SortedByKmerId,
            maxSkip, maxVertexCountPerKmer, debug, graph, alignment);

        // Compute the AlignmentInfo.
        const AlignmentInfo alignmentInfo(alignment);
        uint32_t leftTrim;
        uint32_t rightTrim;
        tie(leftTrim, rightTrim) = computeTrim(
            orientedReadId0,
            orientedReadId1,
            alignmentInfo);

        cout << orientedReadId0 << " " << orientedReadId1 << " " << alignmentInfo.markerCount;
        if(alignmentInfo.markerCount) {
            cout << " " << leftTrim << " " << rightTrim;
            if(alignmentInfo.markerCount >= minAlignedMarkerCount && leftTrim<=maxTrim && rightTrim<=maxTrim) {
                cout << " good";
                ++goodAlignmentCount;
            }
        }
        cout << endl;

    }
    cout << "Found " << goodAlignmentCount << " good alignments among ";
    cout << overlapTable[orientedReadId0.getValue()].size() << " overlaps." << endl;

}



// Given two oriented reads and their computed AlignmentInfo,
// compute the left and right trim.
// This is the minimum number of bases (over the two reads)
// that are excluded from the alignment on each size.
// This takes as input markers sorted by position.
pair<uint32_t, uint32_t> Assembler::computeTrim(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const AlignmentInfo& alignmentInfo)
{
    const uint32_t baseCount0 = uint32_t(reads[orientedReadId0.getReadId()].baseCount);
    const uint32_t baseCount1 = uint32_t(reads[orientedReadId1.getReadId()].baseCount);

    const auto markers0 = markers[orientedReadId0.getValue()];
    const auto markers1 = markers[orientedReadId1.getValue()];

    if(alignmentInfo.markerCount) {
        const CompressedMarker& firstMarker0 = markers0[alignmentInfo.firstOrdinals.first];
        const CompressedMarker& firstMarker1 = markers1[alignmentInfo.firstOrdinals.second];
        const CompressedMarker& lastMarker0 = markers0[alignmentInfo.lastOrdinals.first];
        const CompressedMarker& lastMarker1 = markers1[alignmentInfo.lastOrdinals.second];
        return make_pair(
            min(firstMarker0.position, firstMarker1.position),
            min(baseCount0-1-lastMarker0.position, baseCount1-1-lastMarker1.position)
            );
    } else {
        // The alignment is empty.
        // The trim equal the entire lenght of the shortest read.
        const uint32_t trim = min(baseCount0, baseCount1);
        return make_pair(trim, trim);

    }

}



// Compute an Alignment for each Overlap, but  only store the AlignmentInfo.
// Optionally, the alignments are used for creation of the global marker graph.
void Assembler::computeAllAlignments(

    // Minimum number of MinHash hits for an alignment to be computed.
    size_t minFrequency,

    // The  maximum number of vertices in the alignment graph
    // that we allow a single k-mer to generate.
    size_t alignmentMaxVertexCountPerKmer,

    // The maximum ordinal skip to be tolerated between successive markers
    // in the alignment.
    size_t maxSkip,

    // Minimum number of alignment markers for an alignment to be used.
    size_t minAlignedMarkerCount,

    // Maximum left/right trim (in bases) for an alignment to be used.
    size_t maxTrim,

    // Number of threads. If zero, a number of threads equal to
    // the number of virtual processors is used.
    size_t threadCount,

    // Flag to control computation of the global marker graph.
    // This should normally be true.
    bool computeGlobalMarkerGraph
)
{
    const auto tBegin = steady_clock::now();
    cout << timestamp << "Begin computing alignments for ";
    cout << overlaps.size() << " overlaps." << endl;

    // Check that we have what we need.
    checkReadsAreOpen();
    checkKmersAreOpen();
    checkMarkersAreOpen();
    checkOverlapsAreOpen();

    // Store parameters so they are accessible to the threads.
    auto& data = computeAllAlignmentsData;
    data.maxSkip = maxSkip;
    data.maxVertexCountPerKmer = alignmentMaxVertexCountPerKmer;
    data.minAlignedMarkerCount = minAlignedMarkerCount;
    data.maxTrim = maxTrim;
    data.computeGlobalMarkerGraph = computeGlobalMarkerGraph;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Create the alignmentInfos.
    alignmentInfos.createNew(largeDataName("AlignmentInfos"), largeDataPageSize);
    alignmentInfos.resize(overlaps.size());

    // If requested, initialize computation of the global marker graph.
    if(computeGlobalMarkerGraph) {
        data.orientedMarkerCount = markers.totalSize();
        data.disjointSetsData.createNew(
            largeDataName("tmp-DisjointSetData"),
            largeDataPageSize);
        data.disjointSetsData.resize(data.orientedMarkerCount);
        data.disjointSetsPointer = std::make_shared<DisjointSets>(
            data.disjointSetsData.begin(),
            data.orientedMarkerCount
            );
    }

    // Do the computation.
    const size_t batchSize = 10000;
    setupLoadBalancing(overlaps.size(), batchSize);
    runThreads(&Assembler::computeAllAlignmentsThreadFunction,
        threadCount, "threadLogs/computeAllAlignments");
    cout << timestamp << "Multithreaded portion completed." << endl;




    // Use the disjoint sets to find the vertex of the global
    // marker graph that each oriented marker belongs to.
    if(computeGlobalMarkerGraph) {

        // Compute the global marker graph vertex corresponding
        // to each global OrientedMarkerId.
        globalMarkerGraphVertex.createNew(
            largeDataName("GlobalMarkerGraphVertex"),
            largeDataPageSize);
        globalMarkerGraphVertex.resize(markers.totalSize());
        DisjointSets& disjointSets = *data.disjointSetsPointer;
        for(MarkerId i=0; i<globalMarkerGraphVertex.size(); i++) {
            globalMarkerGraphVertex[i] = disjointSets.find(i);
        }

        // Clean up.
        data.disjointSetsPointer = 0;
        data.disjointSetsData.remove();



        // Now we want to replace the vertices with contiguous ids
        // beginning at 0 and sorted by decreasing size.

        // First, we have to count the number of oriented markers for each
        // vertex, with the current numbering.
        using VertexId = GlobalMarkerGraphVertexId;
        MemoryMapped::Vector<VertexId> count;
        count.createNew(
            largeDataName("tmp-Count"),
            largeDataPageSize);
        count.resize(globalMarkerGraphVertex.size());
        for(const auto i: globalMarkerGraphVertex) {
            ++count[i];
        }

        // Gather the ones with non-zero count.
        MemoryMapped::Vector< pair<VertexId, VertexId> > countTable;
        countTable.createNew(
            largeDataName("tmp-CountTable"),
            largeDataPageSize);
        for(VertexId i=0; i<count.size(); i++) {
            const VertexId n = count[i];
            if(n > 0) {
                countTable.push_back(make_pair(i, n));
            }
        }
        sort(countTable.begin(), countTable.end(),
            OrderPairsBySecondGreaterThenByFirstLess<VertexId, VertexId>());
        cout << "The global marker graph has " << countTable.size();
        cout << " vertices  for " << count.size() << " oriented markers." << endl;

        // Now we can remove the count vector.
        count.remove();

        // Create a histogram of number of markers per vertex.
        vector<uint64_t> histogram(countTable.front().second + 1, 0ULL);
        for(const auto& p: countTable) {
            ++histogram[p.second];
        }
        ofstream csv("GlobalMarkerGraphHistogram.csv");
        csv << "MarkerCount,Frequency\n";
        for(size_t i=0; i<histogram.size(); i++) {
            const auto n = histogram[i];
            if(n) {
                csv << i << "," << n << "\n";
            }
        }

        // In the count table, replace the count by the rank,
        // (the rank becomes the new vertex id),
        // then sort by first (the old vertex id).
        for(VertexId i=0; i<countTable.size(); i++) {
            countTable[i].second = i;
        }
        sort(countTable.begin(), countTable.end(),
            OrderPairsByFirstOnly<VertexId, VertexId>());

        // Redefine globalMarkerGraphVertex using this count table.
        for(VertexId& v: globalMarkerGraphVertex) {
            const auto it = std::lower_bound(countTable.begin(), countTable.end(),
                make_pair(v, 0),
                OrderPairsByFirstOnly<VertexId, VertexId>());
                CZI_ASSERT(it != countTable.end());
                CZI_ASSERT(it->first == v);
                v = it->second;
        }
        const VertexId vertexCount = countTable.size();
        countTable.remove();

        // Gather the oriented marker ids of each vertex of the global marker graph.
        globalMarkerGraphVertices.createNew(
            largeDataName("GlobalMarkerGraphVertices"),
            largeDataPageSize);
        globalMarkerGraphVertices.beginPass1(vertexCount);
        for(const VertexId& v: globalMarkerGraphVertex) {
            globalMarkerGraphVertices.incrementCount(v);
        }
        globalMarkerGraphVertices.beginPass2();
        for(VertexId i=0; i<globalMarkerGraphVertex.size(); i++) {
            globalMarkerGraphVertices.store(globalMarkerGraphVertex[i], i);
        }
        globalMarkerGraphVertices.endPass2();

        // Sort the markers in each vertex.
        for(VertexId i=0; i<globalMarkerGraphVertices.size(); i++) {
            sort(globalMarkerGraphVertices.begin(i), globalMarkerGraphVertices.end(i));
        }
    }


    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "Computation of alignments ";
    if(computeGlobalMarkerGraph) {
        cout << "and global marker graph ";
    }
    cout << "completed in " << tTotal << " s." << endl;
}



void Assembler::computeAllAlignmentsThreadFunction(size_t threadId)
{
    ostream& out = getLog(threadId);

    array<OrientedReadId, 2> orientedReadIds;
    array<OrientedReadId, 2> orientedReadIdsOppositeStrand;
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    AlignmentGraph graph;
    Alignment alignment;

    const bool debug = false;
    const auto& data = computeAllAlignmentsData;
    const size_t maxSkip = data.maxSkip;
    const size_t maxVertexCountPerKmer = data.maxVertexCountPerKmer;
    const size_t minAlignedMarkerCount = data.minAlignedMarkerCount;
    const size_t maxTrim = data.maxTrim;

    const std::shared_ptr<DisjointSets> disjointSetsPointer = data.disjointSetsPointer;

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << "Working on batch " << begin << " " << end << endl;

        for(size_t i=begin; i!=end; i++) {
            const Overlap& overlap = overlaps[i];
            AlignmentInfo& alignmentInfo = alignmentInfos[i];

            // Get the oriented read ids.
            orientedReadIds[0] = OrientedReadId(overlap.readIds[0], 0);
            orientedReadIds[1] = OrientedReadId(overlap.readIds[1], overlap.isSameStrand ? 0 : 1);

            // Get the oriented read ids for the opposite strand.
            orientedReadIdsOppositeStrand = orientedReadIds;
            orientedReadIdsOppositeStrand[0].flipStrand();
            orientedReadIdsOppositeStrand[1].flipStrand();


            // out << timestamp << "Working on " << i << " " << orientedReadIds[0] << " " << orientedReadIds[1] << endl;

            // Get the markers for the two oriented reads in this Overlap.
            for(size_t j=0; j<2; j++) {
                getMarkersSortedByKmerId(orientedReadIds[j], markersSortedByKmerId[j]);
            }

            // Compute the alignment.
            alignOrientedReads(
                markersSortedByKmerId[0],
                markersSortedByKmerId[1],
                maxSkip, maxVertexCountPerKmer, debug, graph, alignment);

            // Store the AlignmentInfo.
            alignmentInfo.create(alignment);

            // If not doing the global marker graph, we are done here.
            if(!data.computeGlobalMarkerGraph) {
                continue;
            }

            // If the alignment has too few markers, don't use it
            // for the global marker graph.
            if(alignment.ordinals.size() < minAlignedMarkerCount) {
                continue;
            }

            // If the alignment has too much trim, skip it.
            uint32_t leftTrim;
            uint32_t rightTrim;
            tie(leftTrim, rightTrim) = computeTrim(
                orientedReadIds[0],
                orientedReadIds[1],
                alignmentInfo);
            if(leftTrim>maxTrim || rightTrim>maxTrim) {
                continue;
            }



            // If getting here, this is a good alignment.
            // In the global marker graph, merge pairs
            // of aligned markers.
            for(const auto& p: alignment.ordinals) {
                const uint32_t ordinal0 = p.first;
                const uint32_t ordinal1 = p.second;
                MarkerId orientedMarkerId0 =
                    getGlobalOrientedMarkerId(orientedReadIds[0], ordinal0);
                MarkerId orientedMarkerId1 =
                    getGlobalOrientedMarkerId(orientedReadIds[1], ordinal1);
                disjointSetsPointer->unite(orientedMarkerId0, orientedMarkerId1);

                // Also do it for the corresponding oriented markers on
                // the opposite strand.
                const uint32_t ordinal0OppositeStrand = uint32_t(markersSortedByKmerId[0].size()) - 1 - ordinal0;
                const uint32_t ordinal1OppositeStrand = uint32_t(markersSortedByKmerId[1].size()) - 1 - ordinal1;
                const MarkerId orientedMarkerId0OppositeStrand =
                    getGlobalOrientedMarkerId(orientedReadIdsOppositeStrand[0], ordinal0OppositeStrand);
                const MarkerId orientedMarkerId1OppositeStrand =
                    getGlobalOrientedMarkerId(orientedReadIdsOppositeStrand[1], ordinal1OppositeStrand);
                disjointSetsPointer->unite(
                    orientedMarkerId0OppositeStrand,
                    orientedMarkerId1OppositeStrand);
            }
        }
    }

}



void Assembler::accessAlignmentInfos()
{
    alignmentInfos.accessExistingReadOnly(
        largeDataName("AlignmentInfos"));
}



void Assembler::checkAlignmentInfosAreOpen()
{
    if(!alignmentInfos.isOpen) {
        throw runtime_error("Alignment infos are not accessible.");
    }
}

