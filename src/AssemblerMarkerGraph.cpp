// Shasta.
#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "ConsensusCaller.hpp"
#include "LocalMarkerGraph.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// SeqAn.
#include <seqan/graph_msa.h>
#include <seqan/version.h>

// Spoa.
#include "spoa/spoa.hpp"

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "chrono.hpp"
#include <queue>



// Loop over all alignments in the read graph
// to create vertices of the global marker graph.
// Throw away vertices with coverage (number of markers)
// less than minCoverage or more than maxCoverage.
// Also throw away "bad" vertices - that is, vertices
// with more than one marker on the same oriented read.
void Assembler::createMarkerGraphVertices(

    // The  maximum number of vertices in the alignment graph
    // that we allow a single k-mer to generate.
    size_t alignmentMaxVertexCountPerKmer,

    // The maximum ordinal skip to be tolerated between successive markers
    // in the alignment.
    size_t maxSkip,

    // Minimum coverage (number of markers) for a vertex
    // of the marker graph to be kept.
    size_t minCoverage,

    // Maximum coverage (number of markers) for a vertex
    // of the marker graph to be kept.
    size_t maxCoverage,

    // Number of threads. If zero, a number of threads equal to
    // the number of virtual processors is used.
    size_t threadCount
)
{

    // Flag to control debug output.
    // Only turn on for a very small test run with just a few reads.
    const bool debug = false;

    // using VertexId = GlobalMarkerGraphVertexId;
    // using CompressedVertexId = CompressedGlobalMarkerGraphVertexId;

    const auto tBegin = steady_clock::now();
    cout << timestamp << "Begin computing marker graph vertices." << endl;

    // Check that we have what we need.
    checkReadsAreOpen();
    checkKmersAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();
    checkReadGraphIsOpen();

    // Store parameters so they are accessible to the threads.
    auto& data = createMarkerGraphVerticesData;
    data.maxSkip = maxSkip;
    data.maxVertexCountPerKmer = alignmentMaxVertexCountPerKmer;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Initialize computation of the global marker graph.
    data.orientedMarkerCount = markers.totalSize();
    data.disjointSetsData.createNew(
        largeDataName("tmp-DisjointSetData"),
        largeDataPageSize);
    data.disjointSetsData.reserveAndResize(data.orientedMarkerCount);
    data.disjointSetsPointer = std::make_shared<DisjointSets>(
        data.disjointSetsData.begin(),
        data.orientedMarkerCount
        );



    // Update the disjoint set data structure for each alignment
    // in the read graph.
    cout << "Begin processing " << readGraphEdges.size() << " alignments in the read graph." << endl;
    cout << timestamp << "Disjoint set computation begins." << endl;
    size_t batchSize = 10000;
    setupLoadBalancing(readGraphEdges.size(), batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction1,
        threadCount, "threadLogs/createMarkerGraphVertices1");
    cout << timestamp << "Disjoint set computation completed." << endl;



    // Find the disjoint set that each oriented marker was assigned to.
    cout << timestamp << "Finding the disjoint set that each oriented marker was assigned to." << endl;
    data.disjointSetTable.createNew(
        largeDataName("tmp-DisjointSetTable"),
        largeDataPageSize);
    data.disjointSetTable.reserveAndResize(data.orientedMarkerCount);
    batchSize = 1000000;
    cout << "Processing " << data.orientedMarkerCount << " oriented markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction2, threadCount,
        "threadLogs/createMarkerGraphVertices2");

    // Free the disjoint set data structure.
    data.disjointSetsPointer = 0;
    data.disjointSetsData.remove();


    // Debug output.
    if(debug) {
        ofstream out("DisjointSetTable-initial.csv");
        for(MarkerId markerId=0; markerId<data.orientedMarkerCount; markerId++) {
            out << markerId << "," << data.disjointSetTable[markerId] << "\n";
        }

    }



    // Count the number of markers in each disjoint set
    // and store it in data.workArea.
    // We don't want to combine this with the previous block
    // because it would significantly increase the peak memory usage.
    // This way, we allocate data.workArea only after freeing the
    // disjoint set data structure.
    cout << timestamp << "Counting the number of markers in each disjoint set." << endl;
    data.workArea.createNew(
        largeDataName("tmp-WorkArea"),
        largeDataPageSize);
    data.workArea.reserveAndResize(data.orientedMarkerCount);
    fill(data.workArea.begin(), data.workArea.end(), 0ULL);
    cout << "Processing " << data.orientedMarkerCount << " oriented markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction3, threadCount,
        "threadLogs/createMarkerGraphVertices3");



    // Debug output.
    if(debug) {
        ofstream out("WorkArea-initial-count.csv");
        for(MarkerId markerId=0; markerId<data.orientedMarkerCount; markerId++) {
            out << markerId << "," << data.workArea[markerId] << "\n";
        }

    }



    // At this point, data.workArea contains the number of oriented markers in
    // each disjoint set.
    // Replace it with a new numbering, counting only disjoint sets
    // with size not less than minCoverage and not greater than maxCoverage.
    // Note that this numbering is not yet the final vertex numbering,
    // as we will later remove "bad" vertices
    // (vertices with more than one marker on the same read).
    // This block is recursive and cannot be multithreaded.
    cout << timestamp << "Renumbering the disjoint sets." << endl;
    GlobalMarkerGraphVertexId newDisjointSetId = 0ULL;
    for(GlobalMarkerGraphVertexId oldDisjointSetId=0;
        oldDisjointSetId<data.orientedMarkerCount; ++oldDisjointSetId) {
        auto& w = data.workArea[oldDisjointSetId];
        const GlobalMarkerGraphVertexId markerCount = w;
        if(markerCount<minCoverage || markerCount>maxCoverage) {
            w = invalidGlobalMarkerGraphVertexId;
        } else {
            w = newDisjointSetId;
            ++newDisjointSetId;
        }
    }
    const auto disjointSetCount = newDisjointSetId;
    cout << "Kept " << disjointSetCount << " disjoint sets with coverage in the requested range." << endl;



    // Debug output.
    if(debug) {
        ofstream out("WorkArea-initial-renumbering.csv");
        for(MarkerId markerId=0; markerId<data.orientedMarkerCount; markerId++) {
            out << markerId << "," << data.workArea[markerId] << "\n";
        }
    }


    // Reassign vertices to disjoint sets using this new numbering.
    // Vertices assigned to no disjoint set will store invalidGlobalMarkerGraphVertexId.
    // This could be multithreaded if necessary.
    cout << timestamp << "Assigning vertices to renumbered disjoint sets." << endl;
    for(GlobalMarkerGraphVertexId markerId=0;
        markerId<data.orientedMarkerCount; ++markerId) {
        auto& d = data.disjointSetTable[markerId];
        const auto oldId = d;
        const auto newId = data.workArea[oldId];
        d = newId;
    }
    // We no longer need the workArea.
    data.workArea.remove();



    // Debug output.
    if(debug) {
        ofstream out("DisjointSetTable-renumbered.csv");
        for(MarkerId markerId=0; markerId<data.orientedMarkerCount; markerId++) {
            out << markerId << "," << data.disjointSetTable[markerId] << "\n";
        }
    }


    // At this point, data.disjointSetTable contains, for each oriented marker,
    // the disjoint set that the oriented marker is assigned to,
    // and all disjoint sets have a number of markers in the requested range.
    // Vertices not assigned to any disjoint set store a disjoint set
    // equal to invalidGlobalMarkerGraphVertexId.



    // Gather the markers in each disjoint set.
    data.disjointSetMarkers.createNew(
        largeDataName("tmp-DisjointSetMarkers"),
        largeDataPageSize);
    cout << timestamp << "Gathering markers in disjoint sets, pass1." << endl;
    data.disjointSetMarkers.beginPass1(disjointSetCount);
    cout << timestamp << "Processing " << data.orientedMarkerCount << " oriented markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction4, threadCount,
        "threadLogs/createMarkerGraphVertices4");
    cout << timestamp << "Gathering markers in disjoint sets, pass2." << endl;
    data.disjointSetMarkers.beginPass2();
    cout << timestamp << "Processing " << data.orientedMarkerCount << " oriented markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction5, threadCount,
        "threadLogs/createMarkerGraphVertices5");
    data.disjointSetMarkers.endPass2();



    // Sort the markers in each disjoint set.
    cout << timestamp << "Sorting the markers in each disjoint set." << endl;
    setupLoadBalancing(disjointSetCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction6, threadCount,
        "threadLogs/createMarkerGraphVertices6");



    // Flag disjoint sets that contain more than one marker on the same read.
    data.isBadDisjointSet.createNew(
        largeDataName("tmp-IsBadDisjointSet"),
        largeDataPageSize);
    data.isBadDisjointSet.reserveAndResize(disjointSetCount);
    cout << timestamp << "Flagging bad disjoint sets." << endl;
    setupLoadBalancing(disjointSetCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction7, threadCount,
        "threadLogs/createMarkerGraphVertices7");
    const size_t badDisjointSetCount = std::count(
        data.isBadDisjointSet.begin(), data.isBadDisjointSet.end(), true);
    cout << "Found " << badDisjointSetCount << " bad disjoint sets "
        "with more than one marker on a single read." << endl;



    // Renumber the disjoint sets again, this time without counting the ones marked as bad.
    cout << timestamp << "Renumbering disjoint sets to remove the bad ones." << endl;
    data.workArea.createNew(
        largeDataName("tmp-WorkArea"),
        largeDataPageSize);
    data.workArea.reserveAndResize(disjointSetCount);
    newDisjointSetId = 0ULL;
    for(GlobalMarkerGraphVertexId oldDisjointSetId=0;
        oldDisjointSetId<disjointSetCount; ++oldDisjointSetId) {
        auto& w = data.workArea[oldDisjointSetId];
        if(data.isBadDisjointSet[oldDisjointSetId]) {
            w = invalidGlobalMarkerGraphVertexId;
        } else {
            w = newDisjointSetId;
            ++newDisjointSetId;
        }
    }
    CZI_ASSERT(newDisjointSetId + badDisjointSetCount == disjointSetCount);



    // Debug output.
    if(debug) {
        ofstream out("WorkArea-final-renumbering.csv");
        for(MarkerId markerId=0; markerId<disjointSetCount; markerId++) {
            out << markerId << "," << data.workArea[markerId] << "\n";
        }
    }


    // Compute the final disjoint set number for each marker.
    // That becomes the vertex id assigned to that marker.
    // This could be multithreaded.
    cout << timestamp << "Assigning vertex ids to markers." << endl;
    globalMarkerGraphVertex.createNew(
        largeDataName("GlobalMarkerGraphVertex"),
        largeDataPageSize);
    globalMarkerGraphVertex.reserveAndResize(data.orientedMarkerCount);
    for(GlobalMarkerGraphVertexId markerId=0;
        markerId<data.orientedMarkerCount; ++markerId) {
        auto oldValue = data.disjointSetTable[markerId];
        if(oldValue == invalidGlobalMarkerGraphVertexId) {
            globalMarkerGraphVertex[markerId] = invalidCompressedGlobalMarkerGraphVertexId;
        } else {
            globalMarkerGraphVertex[markerId] = data.workArea[oldValue];
        }
    }



    // Store the disjoint sets that are not marker bad.
    // Each corresponds to a vertex of the global marker graph.
    // This could be multithreaded.
    cout << timestamp << "Gathering the markers of each vertex of the marker graph." << endl;
    globalMarkerGraphVertices.createNew(
        largeDataName("GlobalMarkerGraphVertices"),
        largeDataPageSize);
    for(GlobalMarkerGraphVertexId oldDisjointSetId=0;
        oldDisjointSetId<disjointSetCount; ++oldDisjointSetId) {
        if(data.isBadDisjointSet[oldDisjointSetId]) {
            continue;
        }
        globalMarkerGraphVertices.appendVector();
        const auto markers = data.disjointSetMarkers[oldDisjointSetId];
        for(const MarkerId markerId: markers) {
            globalMarkerGraphVertices.append(markerId);
        }
    }
    data.isBadDisjointSet.remove();
    data.workArea.remove();
    data.disjointSetMarkers.remove();
    data.disjointSetTable.remove();


    // Check that the data structures we created are consistent with each other.
    // This could be expensive. Remove when we know this code works.
    cout << timestamp << "Checking marker graph vertices." << endl;
    checkMarkerGraphVertices(minCoverage, maxCoverage);



    const auto tEnd = steady_clock::now();
    const double tTotal = seconds(tEnd - tBegin);
    cout << timestamp << "Computation of global marker graph vertices ";
    cout << "completed in " << tTotal << " s." << endl;

}



void Assembler::createMarkerGraphVerticesThreadFunction1(size_t threadId)
{
    ostream& out = getLog(threadId);

    array<OrientedReadId, 2> orientedReadIds;
    array<OrientedReadId, 2> orientedReadIdsOppositeStrand;
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    AlignmentGraph graph;
    Alignment alignment;

    const bool debug = false;
    auto& data = createMarkerGraphVerticesData;
    const size_t maxSkip = data.maxSkip;
    const size_t maxVertexCountPerKmer = data.maxVertexCountPerKmer;

    const std::shared_ptr<DisjointSets> disjointSetsPointer = data.disjointSetsPointer;

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << "Working on batch " << begin << " " << end << endl;

        for(size_t i=begin; i!=end; i++) {
            const ReadGraphEdge& readGraphEdge = readGraphEdges[i];
            const uint64_t alignmentId = readGraphEdge.alignmentId;
            const OrientedReadPair& candidate = alignmentData[alignmentId];
            CZI_ASSERT(candidate.readIds[0] < candidate.readIds[1]);

            // If either of the reads is flagged chimeric, skip it.
            if(isChimericRead[candidate.readIds[0]] || isChimericRead[candidate.readIds[1]]) {
                continue;
            }

            // Get the oriented read ids, with the first one on strand 0.
            orientedReadIds[0] = OrientedReadId(candidate.readIds[0], 0);
            orientedReadIds[1] = OrientedReadId(candidate.readIds[1], candidate.isSameStrand ? 0 : 1);

            // Get the oriented read ids for the opposite strand.
            orientedReadIdsOppositeStrand = orientedReadIds;
            orientedReadIdsOppositeStrand[0].flipStrand();
            orientedReadIdsOppositeStrand[1].flipStrand();


            // out << timestamp << "Working on " << i << " " << orientedReadIds[0] << " " << orientedReadIds[1] << endl;

            // Get the markers for the two oriented reads in this candidate.
            for(size_t j=0; j<2; j++) {
                getMarkersSortedByKmerId(orientedReadIds[j], markersSortedByKmerId[j]);
            }

            // Compute the Alignment.
            alignOrientedReads(
                markersSortedByKmerId[0],
                markersSortedByKmerId[1],
                maxSkip, maxVertexCountPerKmer, debug, graph, alignment);


            // In the global marker graph, merge pairs
            // of aligned markers.
            for(const auto& p: alignment.ordinals) {
                const uint32_t ordinal0 = p.first;
                const uint32_t ordinal1 = p.second;
                const MarkerId markerId0 = getMarkerId(orientedReadIds[0], ordinal0);
                const MarkerId markerId1 = getMarkerId(orientedReadIds[1], ordinal1);
                CZI_ASSERT(markers.begin()[markerId0].kmerId == markers.begin()[markerId1].kmerId);
                disjointSetsPointer->unite(markerId0, markerId1);

                // Also do it for the corresponding oriented markers on
                // the opposite strand.
                const uint32_t ordinal0OppositeStrand = uint32_t(markersSortedByKmerId[0].size()) - 1 - ordinal0;
                const uint32_t ordinal1OppositeStrand = uint32_t(markersSortedByKmerId[1].size()) - 1 - ordinal1;
                const MarkerId markerId0OppositeStrand =
                    getMarkerId(orientedReadIdsOppositeStrand[0], ordinal0OppositeStrand);
                const MarkerId markerId1OppositeStrand =
                    getMarkerId(orientedReadIdsOppositeStrand[1], ordinal1OppositeStrand);
                CZI_ASSERT(markers.begin()[markerId0OppositeStrand].kmerId == markers.begin()[markerId1OppositeStrand].kmerId);
                disjointSetsPointer->unite(
                    markerId0OppositeStrand,
                    markerId1OppositeStrand);
            }
        }
    }

}



void Assembler::createMarkerGraphVerticesThreadFunction2(size_t threadId)
{
    DisjointSets& disjointSets = *createMarkerGraphVerticesData.disjointSetsPointer;
    auto& disjointSetTable = createMarkerGraphVerticesData.disjointSetTable;

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            const uint64_t disjointSetId = disjointSets.find(i);
            disjointSetTable[i] = disjointSetId;
        }
    }
}



void Assembler::createMarkerGraphVerticesThreadFunction3(size_t threadId)
{
    const auto& disjointSetTable = createMarkerGraphVerticesData.disjointSetTable;
    auto& workArea = createMarkerGraphVerticesData.workArea;

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            const uint64_t disjointSetId = disjointSetTable[i];

            // Increment the set size in a thread-safe way.
            __sync_fetch_and_add(&workArea[disjointSetId], 1ULL);
        }
    }
}



void Assembler::createMarkerGraphVerticesThreadFunction4(size_t threadId)
{
    createMarkerGraphVerticesThreadFunction45(4);
}



void Assembler::createMarkerGraphVerticesThreadFunction5(size_t threadId)
{
    createMarkerGraphVerticesThreadFunction45(5);
}



void Assembler::createMarkerGraphVerticesThreadFunction6(size_t threadId)
{
    auto& disjointSetMarkers = createMarkerGraphVerticesData.disjointSetMarkers;

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(GlobalMarkerGraphVertexId i=begin; i!=end; ++i) {
            auto markers = disjointSetMarkers[i];
            sort(markers.begin(), markers.end());
        }
    }
}



void Assembler::createMarkerGraphVerticesThreadFunction7(size_t threadId)
{
    const auto& disjointSetMarkers = createMarkerGraphVerticesData.disjointSetMarkers;
    auto& isBadDisjointSet = createMarkerGraphVerticesData.isBadDisjointSet;

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(GlobalMarkerGraphVertexId disjointSetId=begin; disjointSetId!=end; ++disjointSetId) {
            auto markers = disjointSetMarkers[disjointSetId];
            const size_t markerCount = markers.size();
            CZI_ASSERT(markerCount > 0);
            isBadDisjointSet[disjointSetId] = false;
            if(markerCount == 1) {
                continue;
            }
            for(size_t j=1; j<markerCount; j++) {
                const MarkerId& previousMarkerId = markers[j-1];
                const MarkerId& markerId = markers[j];
                OrientedReadId previousOrientedReadId;
                OrientedReadId orientedReadId;
                tie(previousOrientedReadId, ignore) = findMarkerId(previousMarkerId);
                tie(orientedReadId, ignore) = findMarkerId(markerId);
                if(orientedReadId.getReadId() == previousOrientedReadId.getReadId()) {
                    isBadDisjointSet[disjointSetId] = true;
                    break;
                }
            }
        }
    }
}



void Assembler::createMarkerGraphVerticesThreadFunction45(int value)
{
    CZI_ASSERT(value==4 || value==5);
    const auto& disjointSetTable = createMarkerGraphVerticesData.disjointSetTable;
    auto& disjointSetMarkers = createMarkerGraphVerticesData.disjointSetMarkers;

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            const uint64_t disjointSetId = disjointSetTable[i];
            if(disjointSetId == invalidGlobalMarkerGraphVertexId) {
                continue;
            }
            if(value == 4) {
                disjointSetMarkers.incrementCountMultithreaded(disjointSetId);
            } else {
                disjointSetMarkers.storeMultithreaded(disjointSetId, i);
            }
        }
   }
}



// Check for consistency of globalMarkerGraphVertex and globalMarkerGraphVertices.
void Assembler::checkMarkerGraphVertices(
    size_t minCoverage,
    size_t maxCoverage)
{
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    CZI_ASSERT(markers.totalSize() == globalMarkerGraphVertex.size());
    const MarkerId markerCount = markers.totalSize();



    // Dump everything. Only turn on for a very small test case.
    if(false) {
        ofstream out1("globalMarkerGraphVertex.csv");
        out1 << "MarkerId,VertexId\n";
        for(MarkerId markerId=0; markerId<markerCount; markerId++) {
            out1 << markerId << "," << globalMarkerGraphVertex[markerId] << "\n";
        }
        ofstream out2("globalMarkerGraphVertices.csv");
        out1 << "VertexId,MarkerId\n";
        for(GlobalMarkerGraphVertexId vertexId=0;
            vertexId<globalMarkerGraphVertices.size(); vertexId++) {
            const auto markers = globalMarkerGraphVertices[vertexId];
            for(const MarkerId markerId: markers) {
                out2 << vertexId << "," << markerId << "\n";
            }
        }
    }



    for(GlobalMarkerGraphVertexId vertexId=0;
        vertexId!=globalMarkerGraphVertices.size(); vertexId++) {
        CZI_ASSERT(!isBadMarkerGraphVertex(vertexId));
        const auto markers = globalMarkerGraphVertices[vertexId];
        CZI_ASSERT(markers.size() >= minCoverage);
        CZI_ASSERT(markers.size() <= maxCoverage);
        for(const MarkerId markerId: markers) {
            if(globalMarkerGraphVertex[markerId] != vertexId) {
                cout << "Failure at vertex " << vertexId << " marker " << markerId << endl;
            }
            CZI_ASSERT(globalMarkerGraphVertex[markerId] == vertexId);
        }
    }
}



#if 0
void Assembler::createMarkerGraphVerticesThreadFunction3(size_t threadId)
{
    GlobalMarkerGraphVertexId* workArea =
        createMarkerGraphVerticesData.workArea.begin();

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            const MarkerId rawVertexId = globalMarkerGraphVertex[i];
            __sync_fetch_and_add(workArea + rawVertexId, 1ULL);
        }
    }
}



void Assembler::createMarkerGraphVerticesThreadFunction4(size_t threadId)
{
    GlobalMarkerGraphVertexId* workArea =
        createMarkerGraphVerticesData.workArea.begin();
    const uint64_t maxValue = std::numeric_limits<uint64_t>::max();
    const uint64_t maxValueMinus1 = maxValue - 1ULL;

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            const MarkerId rawVertexId = globalMarkerGraphVertex[i];
            CZI_ASSERT(rawVertexId != maxValueMinus1);
            const MarkerId finalVertexId = workArea[rawVertexId];
            globalMarkerGraphVertex[i] = finalVertexId;
        }
    }
}
#endif



void Assembler::accessMarkerGraphVertices()
{
    globalMarkerGraphVertex.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphVertex"));

    globalMarkerGraphVertices.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphVertices"));
}



void Assembler::checkMarkerGraphVerticesAreAvailable()
{
    if(!globalMarkerGraphVertices.isOpen() || !globalMarkerGraphVertex.isOpen) {
        throw runtime_error("Vertices of the marker graph are not accessible.");
    }
}



// Find the vertex of the global marker graph that contains a given marker.
GlobalMarkerGraphVertexId Assembler::getGlobalMarkerGraphVertex(
    ReadId readId,
    Strand strand,
    uint32_t ordinal) const
{
    return getGlobalMarkerGraphVertex(OrientedReadId(readId, strand), ordinal);

}
GlobalMarkerGraphVertexId Assembler::getGlobalMarkerGraphVertex(
    OrientedReadId orientedReadId,
    uint32_t ordinal) const
{
    const MarkerId markerId =  getMarkerId(orientedReadId, ordinal);
    return globalMarkerGraphVertex[markerId];
}


// Find the markers contained in a given vertex of the global marker graph.
// Returns the markers as tuples(read id, strand, ordinal).
vector< tuple<ReadId, Strand, uint32_t> >
    Assembler::getGlobalMarkerGraphVertexMarkers(
        GlobalMarkerGraphVertexId globalMarkerGraphVertexId) const
{
    // Call the lower level function.
    vector< pair<OrientedReadId, uint32_t> > markers;
    getGlobalMarkerGraphVertexMarkers(globalMarkerGraphVertexId, markers);

    // Create the return vector.
    vector< tuple<ReadId, Strand, uint32_t> > returnVector;
    for(const auto& marker: markers) {
        const OrientedReadId orientedReadId = marker.first;
        const uint32_t ordinal = marker.second;
        returnVector.push_back(make_tuple(orientedReadId.getReadId(), orientedReadId.getStrand(), ordinal));
    }
    return returnVector;
}
void Assembler::getGlobalMarkerGraphVertexMarkers(
    GlobalMarkerGraphVertexId vertexId,
    vector< pair<OrientedReadId, uint32_t> >& markers) const
{
    markers.clear();
    for(const MarkerId markerId: globalMarkerGraphVertices[vertexId]) {
        OrientedReadId orientedReadId;
        uint32_t ordinal;
        tie(orientedReadId, ordinal) = findMarkerId(markerId);
        markers.push_back(make_pair(orientedReadId, ordinal));
    }
}



// Find the children of a vertex of the global marker graph.
vector<GlobalMarkerGraphVertexId>
    Assembler::getGlobalMarkerGraphVertexChildren(
    GlobalMarkerGraphVertexId vertexId) const
{
    vector<GlobalMarkerGraphVertexId> children;
    getGlobalMarkerGraphVertexChildren(vertexId, children);
    return children;
}
void Assembler::getGlobalMarkerGraphVertexChildren(
    GlobalMarkerGraphVertexId vertexId,
    vector<GlobalMarkerGraphVertexId>& children,
    bool append
    ) const
{

    if(!append) {
        children.clear();
    }
    if(isBadMarkerGraphVertex(vertexId)) {
        return;
    }

    // Loop over the markers of this vertex.
    for(const MarkerId markerId: globalMarkerGraphVertices[vertexId]) {

        // Find the OrientedReadId and ordinal.
        OrientedReadId orientedReadId;
        uint32_t ordinal;
        tie(orientedReadId, ordinal) = findMarkerId(markerId);

        // Find the next marker that is contained in a vertex.
        ++ordinal;
        for(; ordinal<markers.size(orientedReadId.getValue()); ++ordinal) {

            // Find the vertex id.
            const MarkerId childMarkerId =  getMarkerId(orientedReadId, ordinal);
            const GlobalMarkerGraphVertexId childVertexId =
                globalMarkerGraphVertex[childMarkerId];

            // If this marker correspond to a vertex, add it to our list.
            if(childVertexId != invalidCompressedGlobalMarkerGraphVertexId &&
                !isBadMarkerGraphVertex(childVertexId)) {
                children.push_back(childVertexId);
                break;
            }
        }

    }



    // Deduplicate.
    sort(children.begin(), children.end());
    children.resize(std::unique(children.begin(), children.end()) - children.begin());
}



// This version also returns the oriented read ids and ordinals
// that caused a child to be marked as such.
void Assembler::getGlobalMarkerGraphVertexChildren(
    GlobalMarkerGraphVertexId vertexId,
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerInterval> > >& children,
    vector< pair<GlobalMarkerGraphVertexId, MarkerInterval> >& workArea
    ) const
{
    children.clear();
    workArea.clear();

    if(isBadMarkerGraphVertex(vertexId)) {
        return;
    }

    // Loop over the markers of this vertex.
    for(const MarkerId markerId: globalMarkerGraphVertices[vertexId]) {

        // Find the OrientedReadId and ordinal.
        MarkerInterval info;
        tie(info.orientedReadId, info.ordinals[0]) = findMarkerId(markerId);

        // Find the next marker that is contained in a vertex.
        const auto markerCount = markers.size(info.orientedReadId.getValue());
        for(info.ordinals[1]=info.ordinals[0]+1; info.ordinals[1]<markerCount; ++info.ordinals[1]) {

            // Find the vertex id.
            const MarkerId childMarkerId =  getMarkerId(info.orientedReadId, info.ordinals[1]);
            const GlobalMarkerGraphVertexId childVertexId =
                globalMarkerGraphVertex[childMarkerId];

            // If this marker correspond to a vertex, add it to our list.
            if( childVertexId!=invalidCompressedGlobalMarkerGraphVertexId &&
                !isBadMarkerGraphVertex(childVertexId)) {
                workArea.push_back(make_pair(childVertexId, info));
                break;
            }
        }

    }
    sort(workArea.begin(), workArea.end());



    // Now construct the children by gathering streaks of workArea entries
    // with the same child vertex id.
    for(auto streakBegin=workArea.begin(); streakBegin!=workArea.end(); ) {
        auto streakEnd = streakBegin + 1;
        for(;
            streakEnd!=workArea.end() && streakEnd->first==streakBegin->first;
            streakEnd++) {
        }
        children.resize(children.size() + 1);
        children.back().first = streakBegin->first;
        auto& v = children.back().second;
        for(auto it=streakBegin; it!=streakEnd; it++) {
            v.push_back(it->second);
        }

        // Process the next streak.
        streakBegin = streakEnd;
    }
}



// Find the parents of a vertex of the global marker graph.
vector<GlobalMarkerGraphVertexId>
    Assembler::getGlobalMarkerGraphVertexParents(
    GlobalMarkerGraphVertexId vertexId) const
{
    vector<GlobalMarkerGraphVertexId> parents;
    getGlobalMarkerGraphVertexParents(vertexId, parents);
    return parents;
}
void Assembler::getGlobalMarkerGraphVertexParents(
    GlobalMarkerGraphVertexId vertexId,
    vector<GlobalMarkerGraphVertexId>& parents,
    bool append
    ) const
{

    if(!append) {
        parents.clear();
    }
    if(isBadMarkerGraphVertex(vertexId)) {
        return;
    }

    // Loop over the markers of this vertex.
    for(const MarkerId markerId: globalMarkerGraphVertices[vertexId]) {

        // Find the OrientedReadId and ordinal.
        OrientedReadId orientedReadId;
        uint32_t ordinal;
        tie(orientedReadId, ordinal) = findMarkerId(markerId);

        // Find the previous marker that is contained in a vertex.
        if(ordinal == 0) {
            continue;
        }
        --ordinal;
        for(; ; --ordinal) {

            // Find the vertex id.
            const MarkerId parentMarkerId =  getMarkerId(orientedReadId, ordinal);
            const GlobalMarkerGraphVertexId parentVertexId =
                globalMarkerGraphVertex[parentMarkerId];

            // If this marker correspond to a vertex, add it to our list.
            if(parentVertexId != invalidCompressedGlobalMarkerGraphVertexId &&
                !isBadMarkerGraphVertex(parentVertexId)) {
                parents.push_back(parentVertexId);
                break;
            }

            if(ordinal == 0) {
                break;
            }
        }
    }

    // Deduplicate.
    sort(parents.begin(), parents.end());
    parents.resize(std::unique(parents.begin(), parents.end()) - parents.begin());
}



// This version also returns the oriented read ids and ordinals
// that caused a parent to be marked as such.
void Assembler::getGlobalMarkerGraphVertexParents(
    GlobalMarkerGraphVertexId vertexId,
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerInterval> > >& parents,
    vector< pair<GlobalMarkerGraphVertexId, MarkerInterval> >& workArea
    ) const
{
    parents.clear();
    workArea.clear();

    if(isBadMarkerGraphVertex(vertexId)) {
        return;
    }

    // Loop over the markers of this vertex.
    for(const MarkerId markerId: globalMarkerGraphVertices[vertexId]) {

        // Find the OrientedReadId and ordinal.
        MarkerInterval info;
        tie(info.orientedReadId, info.ordinals[0]) = findMarkerId(markerId);
        if(info.ordinals[0] == 0) {
            continue;
        }

        // Find the previous marker that is contained in a vertex.
        for(info.ordinals[1]=info.ordinals[0]-1; ; --info.ordinals[1]) {

            // Find the vertex id.
            const MarkerId parentMarkerId =  getMarkerId(info.orientedReadId, info.ordinals[1]);
            const GlobalMarkerGraphVertexId parentVertexId =
                globalMarkerGraphVertex[parentMarkerId];

            // If this marker correspond to a vertex, add it to our list.
            if( parentVertexId!=invalidCompressedGlobalMarkerGraphVertexId &&
                !isBadMarkerGraphVertex(parentVertexId)) {
                workArea.push_back(make_pair(parentVertexId, info));
                break;
            }

            if(info.ordinals[1]  == 0) {
                break;
            }
        }

    }
    sort(workArea.begin(), workArea.end());



    // Now construct the parents by gathering streaks of workArea entries
    // with the same child vertex id.
    for(auto streakBegin=workArea.begin(); streakBegin!=workArea.end(); ) {
        auto streakEnd = streakBegin + 1;
        for(;
            streakEnd!=workArea.end() && streakEnd->first==streakBegin->first;
            streakEnd++) {
        }
        parents.resize(parents.size() + 1);
        parents.back().first = streakBegin->first;
        auto& v = parents.back().second;
        for(auto it=streakBegin; it!=streakEnd; it++) {
            v.push_back(it->second);
        }

        // Process the next streak.
        streakBegin = streakEnd;
    }
}



// Python-callable function to get information about an edge of the
// global marker graph. Returns an empty vector if the specified
// edge does not exist.
vector<Assembler::GlobalMarkerGraphEdgeInformation> Assembler::getGlobalMarkerGraphEdgeInformation(
    GlobalMarkerGraphVertexId vertexId0,
    GlobalMarkerGraphVertexId vertexId1
    )
{
    const uint32_t k = uint32_t(assemblerInfo->k);

    // Find the children of vertexId0.
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerInterval> > > children;
    vector< pair<GlobalMarkerGraphVertexId, MarkerInterval> > workArea;
    getGlobalMarkerGraphVertexChildren(vertexId0, children, workArea);

    // Find vertexId1 in the children.
    vector<GlobalMarkerGraphEdgeInformation> v;
    for(const auto& child: children) {
        if(child.first != vertexId1) {
            continue;
        }

        // We found vertexId1. Loop over its MarkerGraphNeighborInfo.
        const auto& childInfos = child.second;
        v.resize(childInfos.size());
        for(size_t i=0; i<v.size(); i++) {
            const auto& childInfo = childInfos[i];
            auto& info = v[i];
            info.readId = childInfo.orientedReadId.getReadId();
            info.strand = childInfo.orientedReadId.getStrand();
            info.ordinal0 = childInfo.ordinals[0];
            info.ordinal1 = childInfo.ordinals[1];

            // Get the positions.
            const MarkerId markerId0 = getMarkerId(childInfo.orientedReadId, info.ordinal0);
            const MarkerId markerId1 = getMarkerId(childInfo.orientedReadId, info.ordinal1);
            const auto& marker0 = markers.begin()[markerId0];
            const auto& marker1 = markers.begin()[markerId1];
            info.position0 = marker0.position;
            info.position1 = marker1.position;

            // Construct the sequence.
            if(info.position1 <= info.position0+k) {
                // The marker overlap.
                info.overlappingBaseCount = info.position0+k - info.position1;
            } else {
                // The markers don't overlap.
                info.overlappingBaseCount = 0;
                for(uint32_t position=info.position0+k; position!=info.position1; position++) {
                    const Base base = getOrientedReadBase(childInfo.orientedReadId, position);
                    info.sequence.push_back(base.character());
                }
            }
        }
    }

    // If getting here, vertexId1 was not found, and we return
    // and empty vector.
    return v;
}



// Lower-level, more efficient version of the above
// (but it returns less information).
void Assembler::getGlobalMarkerGraphEdgeInfo(
    GlobalMarkerGraphVertexId vertexId0,
    GlobalMarkerGraphVertexId vertexId1,
    vector<MarkerInterval>& intervals
    )
{

    if(isBadMarkerGraphVertex(vertexId0)) {
        return;
    }
    if(isBadMarkerGraphVertex(vertexId1)) {
        return;
    }

    // Loop over markers of vertex0.
    for(const MarkerId markerId0: globalMarkerGraphVertices[vertexId0]) {
        OrientedReadId orientedReadId;
        uint32_t ordinal0;
        tie(orientedReadId, ordinal0) = findMarkerId(markerId0);

        // Find the next marker in orientedReadId that is contained in a vertex.
        uint32_t ordinal1 = ordinal0 + 1;;
        for(; ordinal1<markers.size(orientedReadId.getValue()); ++ordinal1) {

            // Find the vertex id.
            const MarkerId markerId1 =  getMarkerId(orientedReadId, ordinal1);
            const GlobalMarkerGraphVertexId vertexId1Candidate =
                globalMarkerGraphVertex[markerId1];

            // If this marker correspond to vertexId1, add it to our list.
            if(vertexId1Candidate != invalidCompressedGlobalMarkerGraphVertexId &&
                !isBadMarkerGraphVertex(vertexId1Candidate)) {
                if(vertexId1Candidate == vertexId1) {
                    intervals.push_back(MarkerInterval(orientedReadId, ordinal0, ordinal1));
                }
                break;
            }
        }
    }
}



// Return true if a vertex of the global marker graph has more than
// one marker for at least one oriented read id.
bool Assembler::isBadMarkerGraphVertex(GlobalMarkerGraphVertexId vertexId) const
{
    // Get the markers of this vertex.
    const auto& vertexMarkerIds = globalMarkerGraphVertices[vertexId];

    // The markers are sorted by OrientedReadId, so we can just check each
    // consecutive pairs.
    for(size_t i=1; i<vertexMarkerIds.size(); i++) {
        const MarkerId markerId0 = vertexMarkerIds[i-1];
        const MarkerId markerId1 = vertexMarkerIds[i];
        OrientedReadId orientedReadId0;
        OrientedReadId orientedReadId1;
        tie(orientedReadId0, ignore) = findMarkerId(markerId0);
        tie(orientedReadId1, ignore) = findMarkerId(markerId1);
        if(orientedReadId0 == orientedReadId1) {
            return true;
        }
    }
    return false;
}



void Assembler::extractLocalMarkerGraph(

    // The ReadId, Strand, and ordinal that identify the
    // marker corresponding to the start vertex
    // for the local marker graph to be created.
    ReadId readId,
    Strand strand,
    uint32_t ordinal,

    // Maximum distance from the start vertex (number of edges in the global marker graph).
    int distance,

    // Minimum coverage for a strong vertex or edge (affects coloring).
    size_t minCoverage

    )
{
    // Create the local marker graph.
    LocalMarkerGraph graph(
        uint32_t(assemblerInfo->k),
        reads,
        assemblerInfo->useRunLengthReads,
        readRepeatCounts,
        markers,
        globalMarkerGraphVertex,
        *consensusCaller);
    extractLocalMarkerGraph(OrientedReadId(readId, strand), ordinal, distance, 0., graph);

    cout << "The local marker graph has " << num_vertices(graph);
    cout << " vertices and " << num_edges(graph) << " edges." << endl;

    // Write it out.
    graph.write("MarkerGraph.dot", minCoverage, distance, false, true);
    graph.write("DetailedMarkerGraph.dot", minCoverage, distance, true, true);

}



bool Assembler::extractLocalMarkerGraph(
    OrientedReadId orientedReadId,
    uint32_t ordinal,
    int distance,
    double timeout,                 // Or 0 for no timeout.
    LocalMarkerGraph& graph
    )
{
    const GlobalMarkerGraphVertexId startVertexId =
        getGlobalMarkerGraphVertex(orientedReadId, ordinal);
    return extractLocalMarkerGraph(startVertexId, distance, timeout, graph);

}



bool Assembler::extractLocalMarkerGraph(
    GlobalMarkerGraphVertexId startVertexId,
    int distance,
    double timeout,                 // Or 0 for no timeout.
    LocalMarkerGraph& graph
    )
{


    using vertex_descriptor = LocalMarkerGraph::vertex_descriptor;
    using edge_descriptor = LocalMarkerGraph::edge_descriptor;
    const auto startTime = steady_clock::now();

    // Add the start vertex.
    if(startVertexId == invalidCompressedGlobalMarkerGraphVertexId) {
        return true;    // Because no timeout occurred.
    }
    const vertex_descriptor vStart = graph.addVertex(startVertexId, 0, globalMarkerGraphVertices[startVertexId]);

    // Some vectors used inside the BFS.
    // Define them here to reduce memory allocation activity.
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerInterval> > > children;
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerInterval> > > parents;
    vector< pair<GlobalMarkerGraphVertexId, MarkerInterval> > workArea;
    vector<MarkerInterval> markerIntervalVector;



    // Do the BFS.
    std::queue<vertex_descriptor> q;
    if(distance > 0) {
        q.push(vStart);
    }
    while(!q.empty()) {

        // See if we exceeded the timeout.
        if(timeout>0. && seconds(steady_clock::now() - startTime) > timeout) {
            graph.clear();
            return false;
        }

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const LocalMarkerGraphVertex& vertex0 = graph[v0];
        const GlobalMarkerGraphVertexId vertexId0 = vertex0.vertexId;
        const int distance0 = vertex0.distance;
        const int distance1 = distance0 + 1;

        // Get the children and parents.
        getGlobalMarkerGraphVertexChildren(vertexId0, children, workArea);
        getGlobalMarkerGraphVertexParents (vertexId0, parents, workArea);

        // Loop over the children.
        for(const auto& p: children) {
            const GlobalMarkerGraphVertexId vertexId1 = p.first;
            bool vertexExists;

            // Find the vertex corresponding to this child, creating it if necessary.
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);
            if(!vertexExists) {
                v1 = graph.addVertex(
                    vertexId1, distance1, globalMarkerGraphVertices[vertexId1]);
                if(distance1 < distance) {
                    q.push(v1);
                }
            }

            // Create the edge v0->v1, if it does not already exist.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            if(!edgeExists) {
                tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
                CZI_ASSERT(edgeExists);

                // Fill in edge information.
                markerIntervalVector.clear();
                const auto& v = p.second;
                for(const MarkerInterval& x: v) {
                    markerIntervalVector.push_back(MarkerInterval(
                        x.orientedReadId, x.ordinals[0], x.ordinals[1]));
                }
                graph.storeEdgeInfo(e, markerIntervalVector);
            }
        }


        // Loop over the parents.
        for(const auto& p: parents) {
            const GlobalMarkerGraphVertexId vertexId1 = p.first;
            bool vertexExists;

            // Find the vertex corresponding to this parent, creating it if necessary.
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);
            if(!vertexExists) {
                v1 = graph.addVertex(
                    vertexId1, distance1, globalMarkerGraphVertices[vertexId1]);
                if(distance1 < distance) {
                    q.push(v1);
                }
            }

            // Create the edge v1->v0, if it does not already exist.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v1, v0, graph);
            if(!edgeExists) {
                tie(e, edgeExists) = boost::add_edge(v1, v0, graph);
                CZI_ASSERT(edgeExists);

                // Fill in edge information.
                markerIntervalVector.clear();
                const auto& v = p.second;
                for(const MarkerInterval& x: v) {
                    markerIntervalVector.push_back(MarkerInterval(
                        x.orientedReadId, x.ordinals[1], x.ordinals[0]));
                }
                graph.storeEdgeInfo(e, markerIntervalVector);
            }
        }

    }



    // The BFS process did not create edges between vertices at maximum distance.
    // Do it now.
    // Loop over all vertices at maximum distance.
    BGL_FORALL_VERTICES(v0, graph, LocalMarkerGraph) {
        const LocalMarkerGraphVertex& vertex0 = graph[v0];
        if(vertex0.distance != distance) {
            continue;
        }

        // Loop over the children that exist in the local marker graph
        // and are also at maximum distance.
        getGlobalMarkerGraphVertexChildren(vertex0.vertexId, children, workArea);
        for(const auto& p: children) {
            const GlobalMarkerGraphVertexId vertexId1 = p.first;
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);

            // If it does not exist in the local marker graph, skip.
            if(!vertexExists) {
                continue;
            }

            // If it is not at maximum distance, skip.
            const LocalMarkerGraphVertex& vertex1 = graph[v1];
            if(vertex1.distance != distance) {
                continue;
            }

            // There is no way we already created this edge.
            // Check that this is the case.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            CZI_ASSERT(!edgeExists);

            // Add the edge.
            tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
            CZI_ASSERT(edgeExists);

            // Fill in edge information.
            markerIntervalVector.clear();
            const auto& v = p.second;
            for(const MarkerInterval& x: v) {
                markerIntervalVector.push_back(MarkerInterval(
                    x.orientedReadId, x.ordinals[0], x.ordinals[1]));
            }
            graph.storeEdgeInfo(e, markerIntervalVector);
        }
    }



    // Fill in the oriented read ids represented in the graph.
    graph.findOrientedReadIds();

    // If using run-length reads, also fill in the ConsensusInfo's
    // for each vertex.
    if(assemblerInfo->useRunLengthReads) {
        graph.computeVertexConsensusInfo();
    }

    return true;
}



bool Assembler::extractLocalMarkerGraphUsingStoredConnectivity(
    OrientedReadId orientedReadId,
    uint32_t ordinal,
    int distance,
    double timeout,                 // Or 0 for no timeout.
    bool useWeakEdges,
    bool usePrunedEdges,
    bool useBubbleEdges,
    bool useBubbleReplacementEdges,
    bool useShortCycleEdges,
    LocalMarkerGraph& graph
    )
{
    const GlobalMarkerGraphVertexId startVertexId =
        getGlobalMarkerGraphVertex(orientedReadId, ordinal);
    return extractLocalMarkerGraphUsingStoredConnectivity(
        startVertexId, distance, timeout,
        useWeakEdges, usePrunedEdges, useBubbleEdges, useBubbleReplacementEdges, useShortCycleEdges,
        graph);

}



bool Assembler::extractLocalMarkerGraphUsingStoredConnectivity(
    GlobalMarkerGraphVertexId startVertexId,
    int distance,
    double timeout,                 // Or 0 for no timeout.
    bool useWeakEdges,
    bool usePrunedEdges,
    bool useBubbleEdges,
    bool useBubbleReplacementEdges,
    bool useShortCycleEdges,
    LocalMarkerGraph& graph
    )
{
    // Sanity check.
    checkMarkerGraphEdgesIsOpen();

    // Some shorthands.
    using vertex_descriptor = LocalMarkerGraph::vertex_descriptor;
    using edge_descriptor = LocalMarkerGraph::edge_descriptor;

    // Start a timer.
    const auto startTime = steady_clock::now();

    // Add the start vertex.
    if(startVertexId == invalidCompressedGlobalMarkerGraphVertexId) {
        return true;    // Because no timeout occurred.
    }
    const vertex_descriptor vStart = graph.addVertex(startVertexId, 0, globalMarkerGraphVertices[startVertexId]);

    // Some vectors used inside the BFS.
    // Define them here to reduce memory allocation activity.
    vector<MarkerInterval> markerIntervals;


    // Do the BFS.
    std::queue<vertex_descriptor> q;
    if(distance > 0) {
        q.push(vStart);
    }
    while(!q.empty()) {

        // See if we exceeded the timeout.
        if(timeout>0. && seconds(steady_clock::now() - startTime) > timeout) {
            graph.clear();
            return false;
        }

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const LocalMarkerGraphVertex& vertex0 = graph[v0];
        const GlobalMarkerGraphVertexId vertexId0 = vertex0.vertexId;
        const int distance0 = vertex0.distance;
        const int distance1 = distance0 + 1;

        // Loop over the children.
        const auto childEdges = markerGraph.edgesBySource[vertexId0];
        for(uint64_t edgeId: childEdges) {
            const auto& edge = markerGraph.edges[edgeId];

            // Skip this edge if the arguments require it.
            if(edge.wasRemovedByTransitiveReduction && !useWeakEdges) {
                continue;
            }
            if(edge.wasPruned && !usePrunedEdges) {
                continue;
            }
            if(edge.isBubbleEdge && !useBubbleEdges) {
                continue;
            }
            if(edge.replacesBubbleEdges && !useBubbleReplacementEdges) {
                continue;
            }
            if(edge.isShortCycleEdge && !useShortCycleEdges) {
                continue;
            }

            const GlobalMarkerGraphVertexId vertexId1 = edge.target;
            CZI_ASSERT(edge.source == vertexId0);
            CZI_ASSERT(vertexId1 < globalMarkerGraphVertices.size());
            CZI_ASSERT(!isBadMarkerGraphVertex(vertexId1));

            // Find the vertex corresponding to this child, creating it if necessary.
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);
            if(!vertexExists) {
                v1 = graph.addVertex(
                    vertexId1, distance1, globalMarkerGraphVertices[vertexId1]);
                if(distance1 < distance) {
                    q.push(v1);
                }
            }

            // Create the edge v0->v1, if it does not already exist.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            if(!edgeExists) {
                tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
                CZI_ASSERT(edgeExists);

                // Fill in edge information.
                const auto storedMarkerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
                markerIntervals.resize(storedMarkerIntervals.size());
                copy(storedMarkerIntervals.begin(), storedMarkerIntervals.end(), markerIntervals.begin());
                graph.storeEdgeInfo(e, markerIntervals);
                graph[e].edgeId = edgeId;

                // Link to assembly graph edge.
                if(assemblyGraph.markerToAssemblyTable.isOpen) {
                    const auto& p = assemblyGraph.markerToAssemblyTable[edgeId];
                    graph[e].assemblyEdgeId = p.first;
                    graph[e].positionInAssemblyEdge = p.second;
                }
            }
        }

        // Loop over the parents.
        const auto parentEdges = markerGraph.edgesByTarget[vertexId0];
        for(uint64_t edgeId: parentEdges) {
            const auto& edge = markerGraph.edges[edgeId];

            // Skip this edge if the arguments require it.
            if(edge.wasRemovedByTransitiveReduction && !useWeakEdges) {
                continue;
            }
            if(edge.wasPruned && !usePrunedEdges) {
                continue;
            }
            if(edge.isBubbleEdge && !useBubbleEdges) {
                continue;
            }
            if(edge.replacesBubbleEdges && !useBubbleReplacementEdges) {
                continue;
            }
            if(edge.isShortCycleEdge && !useShortCycleEdges) {
                continue;
            }

            const GlobalMarkerGraphVertexId vertexId1 = edge.source;
            CZI_ASSERT(edge.target == vertexId0);
            CZI_ASSERT(vertexId1 < globalMarkerGraphVertices.size());

            // Find the vertex corresponding to this child, creating it if necessary.
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);
            if(!vertexExists) {
                v1 = graph.addVertex(
                    vertexId1, distance1, globalMarkerGraphVertices[vertexId1]);
                if(distance1 < distance) {
                    q.push(v1);
                }
            }

            // Create the edge v1->v0, if it does not already exist.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v1, v0, graph);
            if(!edgeExists) {
                tie(e, edgeExists) = boost::add_edge(v1, v0, graph);
                CZI_ASSERT(edgeExists);

                // Fill in edge information.
                const auto storedMarkerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
                markerIntervals.resize(storedMarkerIntervals.size());
                copy(storedMarkerIntervals.begin(), storedMarkerIntervals.end(), markerIntervals.begin());
                graph.storeEdgeInfo(e, markerIntervals);
                graph[e].edgeId = edgeId;

                // Link to assembly graph vertex.
                if(assemblyGraph.markerToAssemblyTable.isOpen) {
                    const auto& p = assemblyGraph.markerToAssemblyTable[edgeId];
                    graph[e].assemblyEdgeId = p.first;
                    graph[e].positionInAssemblyEdge = p.second;
                }
            }
        }

    }


    // The BFS process did not create edges between vertices at maximum distance.
    // Do it now.
    // Loop over all vertices at maximum distance.
    BGL_FORALL_VERTICES(v0, graph, LocalMarkerGraph) {
        const LocalMarkerGraphVertex& vertex0 = graph[v0];
        if(vertex0.distance != distance) {
            continue;
        }
        const GlobalMarkerGraphVertexId vertexId0 = vertex0.vertexId;

        // Loop over the children that exist in the local marker graph
        // and are also at maximum distance.
        const auto childEdges = markerGraph.edgesBySource[vertexId0];
        for(uint64_t edgeId: childEdges) {
            const auto& edge = markerGraph.edges[edgeId];

            // Skip this edge if the arguments require it.
            if(edge.wasRemovedByTransitiveReduction && !useWeakEdges) {
                continue;
            }
            if(edge.wasPruned && !usePrunedEdges) {
                continue;
            }
            if(edge.isBubbleEdge && !useBubbleEdges) {
                continue;
            }
            if(edge.replacesBubbleEdges && !useBubbleReplacementEdges) {
                continue;
            }
            if(edge.isShortCycleEdge && !useShortCycleEdges) {
                continue;
            }

            const GlobalMarkerGraphVertexId vertexId1 = edge.target;
            CZI_ASSERT(edge.source == vertexId0);
            CZI_ASSERT(vertexId1 < globalMarkerGraphVertices.size());

            // See if we have a vertex for this global vertex id.
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);

            // If it does not exist in the local marker graph, skip.
            if(!vertexExists) {
                continue;
            }

            // If it is not at maximum distance, skip.
            const LocalMarkerGraphVertex& vertex1 = graph[v1];
            if(vertex1.distance != distance) {
                continue;
            }

            // There is no way we already created this edge.
            // Check that this is the case.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            CZI_ASSERT(!edgeExists);

            // Add the edge.
            tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
            CZI_ASSERT(edgeExists);

            // Fill in edge information.
            const auto storedMarkerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
            markerIntervals.resize(storedMarkerIntervals.size());
            copy(storedMarkerIntervals.begin(), storedMarkerIntervals.end(), markerIntervals.begin());
            graph.storeEdgeInfo(e, markerIntervals);
            graph[e].edgeId = edgeId;

            // Link to assembly graph vertex.
            if(assemblyGraph.markerToAssemblyTable.isOpen) {
                const auto& p = assemblyGraph.markerToAssemblyTable[edgeId];
                graph[e].assemblyEdgeId = p.first;
                graph[e].positionInAssemblyEdge = p.second;
            }
        }
    }



    // Fill in the oriented read ids represented in the graph.
    graph.findOrientedReadIds();

    // If using run-length reads, also fill in the ConsensusInfo's
    // for each vertex.
    if(assemblerInfo->useRunLengthReads) {
        graph.computeVertexConsensusInfo();
    }

    return true;
}



// Create a local marker graph and return its local assembly path.
// The local marker graph is specified by its start vertex
// and maximum distance (number of edges) form the start vertex.
vector<GlobalMarkerGraphVertexId> Assembler::getLocalAssemblyPath(
    GlobalMarkerGraphVertexId startVertexId,
    int maxDistance
    )
{

    // Create the local marker graph.
    LocalMarkerGraph graph(
        uint32_t(assemblerInfo->k),
        reads,
        assemblerInfo->useRunLengthReads,
        readRepeatCounts,
        markers,
        globalMarkerGraphVertex,
        *consensusCaller);
    extractLocalMarkerGraph(startVertexId, maxDistance, 0., graph);

    // Construct the local assembly path.
    graph.approximateTopologicalSort();
    graph.computeOptimalSpanningTree();
    graph.computeOptimalSpanningTreeBestPath();
    graph.computeLocalAssemblyPath(maxDistance);

    // Get the vertex ids in the assembly path.
    vector<GlobalMarkerGraphVertexId> path;
    if(!graph.localAssemblyPath.empty()) {
        LocalMarkerGraph::edge_descriptor e = graph.localAssemblyPath.front();
        const LocalMarkerGraph::vertex_descriptor v = source(e, graph);
        path.push_back(graph[v].vertexId);
        for(LocalMarkerGraph::edge_descriptor e: graph.localAssemblyPath) {
            const LocalMarkerGraph::vertex_descriptor v = target(e, graph);
            path.push_back(graph[v].vertexId);
        }
    }
    return path;
}



// Compute edges of the global marker graph.
void Assembler::createMarkerGraphEdges(size_t threadCount)
{
    cout << timestamp << "createMarkerGraphEdges begins." << endl;

    // Check that we have what we need.
    checkMarkerGraphVerticesAreAvailable();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Each thread stores the edges it finds in a separate vector.
    markerGraph.threadEdges.resize(threadCount);
    markerGraph.threadEdgeMarkerIntervals.resize(threadCount);
    cout << timestamp << "Processing " << globalMarkerGraphVertices.size();
    cout << " marker graph vertices." << endl;
    setupLoadBalancing(globalMarkerGraphVertices.size(), 100000);
    runThreads(&Assembler::createMarkerGraphEdgesThreadFunction0, threadCount,
        "threadLogs/createMarkerGraphEdges0");

    // Combine the edges found by each thread.
    cout << timestamp << "Combining the edges found by each thread." << endl;
    markerGraph.edges.createNew(
            largeDataName("GlobalMarkerGraphEdges"),
            largeDataPageSize);
    markerGraph.edgeMarkerIntervals.createNew(
            largeDataName("GlobalMarkerGraphEdgeMarkerIntervals"),
            largeDataPageSize);
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        auto& thisThreadEdges = *markerGraph.threadEdges[threadId];
        auto& thisThreadEdgeMarkerIntervals = *markerGraph.threadEdgeMarkerIntervals[threadId];
        CZI_ASSERT(thisThreadEdges.size() == thisThreadEdgeMarkerIntervals.size());
        for(size_t i=0; i<thisThreadEdges.size(); i++) {
            const auto& edge = thisThreadEdges[i];
            const auto edgeMarkerIntervals = thisThreadEdgeMarkerIntervals[i];
            markerGraph.edges.push_back(edge);
            markerGraph.edgeMarkerIntervals.appendVector();
            for(auto edgeMarkerInterval: edgeMarkerIntervals) {
                markerGraph.edgeMarkerIntervals.append(edgeMarkerInterval);
            }
        }
        thisThreadEdges.remove();
        thisThreadEdgeMarkerIntervals.remove();
    }
    CZI_ASSERT(markerGraph.edges.size() == markerGraph.edgeMarkerIntervals.size());
    cout << timestamp << "Found " << markerGraph.edges.size();
    cout << " edges for " << globalMarkerGraphVertices.size() << " vertices." << endl;



    // Now we need to create edgesBySource and edgesByTarget.
    createMarkerGraphEdgesBySourceAndTarget(threadCount);
    cout << timestamp << "createMarkerGraphEdges ends." << endl;
}



void Assembler::createMarkerGraphEdgesBySourceAndTarget(size_t threadCount)
{
    markerGraph.edgesBySource.createNew(
        largeDataName("GlobalMarkerGraphEdgesBySource"),
        largeDataPageSize);
    markerGraph.edgesByTarget.createNew(
        largeDataName("GlobalMarkerGraphEdgesByTarget"),
        largeDataPageSize);

    cout << timestamp << "Create marker graph edges by source and target: pass 1 begins." << endl;
    markerGraph.edgesBySource.beginPass1(globalMarkerGraphVertices.size());
    markerGraph.edgesByTarget.beginPass1(globalMarkerGraphVertices.size());
    setupLoadBalancing(markerGraph.edges.size(), 100000);
    runThreads(&Assembler::createMarkerGraphEdgesThreadFunction1, threadCount);

    cout << timestamp << "Create marker graph edges by source and target: pass 2 begins." << endl;
    markerGraph.edgesBySource.beginPass2();
    markerGraph.edgesByTarget.beginPass2();
    setupLoadBalancing(markerGraph.edges.size(), 100000);
    runThreads(&Assembler::createMarkerGraphEdgesThreadFunction2, threadCount);
    markerGraph.edgesBySource.endPass2();
    markerGraph.edgesByTarget.endPass2();

}



void Assembler::createMarkerGraphEdgesThreadFunction0(size_t threadId)
{
    using std::shared_ptr;
    using std::make_shared;
    ostream& out = getLog(threadId);

    // Create the vector to contain the edges found by this thread.
    shared_ptr< MemoryMapped::Vector<MarkerGraph::Edge> > thisThreadEdgesPointer =
        make_shared< MemoryMapped::Vector<MarkerGraph::Edge> >();
    markerGraph.threadEdges[threadId] = thisThreadEdgesPointer;
    MemoryMapped::Vector<MarkerGraph::Edge>& thisThreadEdges = *thisThreadEdgesPointer;
    thisThreadEdges.createNew(
            largeDataName("tmp-ThreadGlobalMarkerGraphEdges-" + to_string(threadId)),
            largeDataPageSize);

    // Create the vector to contain the marker intervals for edges found by this thread.
    shared_ptr< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> >
        thisThreadEdgeMarkerIntervalsPointer =
        make_shared< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> >();
    markerGraph.threadEdgeMarkerIntervals[threadId] = thisThreadEdgeMarkerIntervalsPointer;
    MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t>&
        thisThreadEdgeMarkerIntervals = *thisThreadEdgeMarkerIntervalsPointer;
    thisThreadEdgeMarkerIntervals.createNew(
            largeDataName("tmp-ThreadGlobalMarkerGraphEdgeMarkerIntervals-" + to_string(threadId)),
            largeDataPageSize);

    // Some things used inside the loop but defined here for performance.
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerInterval> > > children;
    vector< pair<GlobalMarkerGraphVertexId, MarkerInterval> > workArea;
    MarkerGraph::Edge edge;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << begin << endl;

        // Loop over all marker graph vertices assigned to this batch.
        for(GlobalMarkerGraphVertexId vertex0=begin; vertex0!=end; ++vertex0) {
            // out << timestamp << vertex0 << " " << globalMarkerGraphVertices.size(vertex0) << endl;
            edge.source = vertex0;

            getGlobalMarkerGraphVertexChildren(vertex0, children, workArea);
            for(const auto& p: children) {
                const auto vertex1 = p.first;
                const auto& markerIntervals = p.second;
                edge.target = vertex1;
                size_t coverage = markerIntervals.size();
                if(coverage < 256) {
                    edge.coverage = uint8_t(coverage);
                } else {
                    edge.coverage = 255;
                }

                // Store the edge.
                thisThreadEdges.push_back(edge);

                // Store the marker intervals.
                thisThreadEdgeMarkerIntervals.appendVector();
                for(const MarkerInterval markerInterval: markerIntervals) {
                    thisThreadEdgeMarkerIntervals.append(markerInterval);
                }
            }
        }
    }

}



void Assembler::createMarkerGraphEdgesThreadFunction1(size_t threadId)
{
    createMarkerGraphEdgesThreadFunction12(threadId, 1);
}
void Assembler::createMarkerGraphEdgesThreadFunction2(size_t threadId)
{
    createMarkerGraphEdgesThreadFunction12(threadId, 2);
}
void Assembler::createMarkerGraphEdgesThreadFunction12(size_t threadId, size_t pass)
{
    CZI_ASSERT(pass==1 || pass==2);

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all marker graph edges assigned to this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            const auto& edge = markerGraph.edges[i];
            if(pass == 1) {
                markerGraph.edgesBySource.incrementCountMultithreaded(edge.source);
                markerGraph.edgesByTarget.incrementCountMultithreaded(edge.target);
            } else {
                markerGraph.edgesBySource.storeMultithreaded(edge.source, Uint40(i));
                markerGraph.edgesByTarget.storeMultithreaded(edge.target, Uint40(i));
            }
        }
    }

}


void Assembler::accessMarkerGraphEdges(bool accessEdgesReadWrite)
{
    if(accessEdgesReadWrite) {
        markerGraph.edges.accessExistingReadWrite(
            largeDataName("GlobalMarkerGraphEdges"));
        markerGraph.edgeMarkerIntervals.accessExistingReadWrite(
            largeDataName("GlobalMarkerGraphEdgeMarkerIntervals"));
    } else {
        markerGraph.edges.accessExistingReadOnly(
            largeDataName("GlobalMarkerGraphEdges"));
        markerGraph.edgeMarkerIntervals.accessExistingReadOnly(
            largeDataName("GlobalMarkerGraphEdgeMarkerIntervals"));
    }
    markerGraph.edgesBySource.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphEdgesBySource"));
    markerGraph.edgesByTarget.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphEdgesByTarget"));
}



void Assembler::checkMarkerGraphEdgesIsOpen()
{
    CZI_ASSERT(markerGraph.edges.isOpen);
    CZI_ASSERT(markerGraph.edgesBySource.isOpen());
    CZI_ASSERT(markerGraph.edgesByTarget.isOpen());
}



// Locate the edge given the vertices.
const Assembler::MarkerGraph::Edge*
    Assembler::MarkerGraph::findEdge(Uint40 source, Uint40 target) const
{
    const auto edgesWithThisSource = edgesBySource[source];
    for(const uint64_t i: edgesWithThisSource) {
        const Edge& edge = edges[i];
        if(edge.target == target) {
            return &edge;
        }
    }
    return 0;

}



#if 0
// (OLD VERSION)
// This old versdion works pretty well, but it causes contig breaks
// in some relatively rare end cases.
// The new version below uses approximate transitive reduction.

// Find weak edges in the marker graph.
// All edges with coverage up to lowCoverageThreshold are marked as weak.
// All edges with coverage at least highCoverageThreshold as marked as strong.
// The remaining edges are processed in order of increasing coverage,
// and they are flagged as weak if they don't disconnect the local subgraph
// up to maxDistance. Edges are processed in order of increasing coverage,
// and the local subgraph is created using only edges
// currently marked as strong. This way, low coverage edges
// are more likely to be marked as weak, and we never locally disconnect the graph.
void Assembler::flagMarkerGraphWeakEdges(
    size_t lowCoverageThreshold,
    size_t highCoverageThreshold,
    size_t maxDistance)
{
    cout << timestamp << "Flagging weak edges of the marker graph." << endl;
    cout << "The marker graph has " << globalMarkerGraphVertices.size() << " vertices and ";
    cout << markerGraph.edges.size() << " edges." << endl;

    // Initially flag all edges as strong.
    auto& edges = markerGraph.edges;
    for(auto& edge: edges) {
        edge.isWeak = 0;
    }

    // Gather edges for each coverage less than highCoverageThreshold.
    using EdgeId = GlobalMarkerGraphEdgeId;
    MemoryMapped::VectorOfVectors<EdgeId, EdgeId>  edgesByCoverage;
    edgesByCoverage.createNew(
            largeDataName("tmp-EdgesByCoverage"),
            largeDataPageSize);
    edgesByCoverage.beginPass1(highCoverageThreshold);
    for(EdgeId edgeId=0; edgeId!=edges.size(); edgeId++) {
        const MarkerGraph::Edge& edge = edges[edgeId];
        if(edge.coverage < highCoverageThreshold) {
            edgesByCoverage.incrementCount(edge.coverage);
        }
    }
    edgesByCoverage.beginPass2();
    for(EdgeId edgeId=0; edgeId!=edges.size(); edgeId++) {
        const MarkerGraph::Edge& edge = edges[edgeId];
        if(edge.coverage < highCoverageThreshold) {
            edgesByCoverage.store(edge.coverage, edgeId);
        }
    }
    edgesByCoverage.endPass2();

    // Work areas for markerGraphEdgeDisconnectsLocalStrongSubgraph.
    array<vector< vector<EdgeId> >, 2> verticesByDistance;
    array<vector<bool>, 2> vertexFlags;
    for(size_t i=0; i<2; i++) {
        verticesByDistance[i].resize(maxDistance+1);
        vertexFlags[i].resize(globalMarkerGraphVertices.size());
        fill(vertexFlags[i].begin(), vertexFlags[i].end(), false);
    }


    // Process the edges with coverage less than minCoverage
    // in order of increasing coverage.
    for(size_t coverage=0; coverage<highCoverageThreshold; coverage++) {
        const auto& edgesWithThisCoverage = edgesByCoverage[coverage];
        cout << timestamp << "Processing " << edgesWithThisCoverage.size() <<
            " edges with coverage " << coverage << endl;
        size_t count = 0;
        for(const EdgeId edgeId: edgesWithThisCoverage) {
            // CZI_ASSERT(find(vertexFlags[0].begin(), vertexFlags[0].end(), true) == vertexFlags[0].end());
            // CZI_ASSERT(find(vertexFlags[1].begin(), vertexFlags[1].end(), true) == vertexFlags[1].end());
            CZI_ASSERT(edges[edgeId].isWeak == 0);
            if(coverage<=lowCoverageThreshold || !markerGraphEdgeDisconnectsLocalStrongSubgraph(
                edgeId, maxDistance, verticesByDistance, vertexFlags)) {
                edges[edgeId].isWeak = 1;
                ++ count;
            }
        }
        cout << "Out of " << edgesWithThisCoverage.size() <<
            " edges with coverage " << coverage <<
            ", " << count << " were marked as weak." << endl;
    }
    edgesByCoverage.remove();


    // Count the number of edges that were flagged as weak.
    uint64_t weakEdgeCount = 0;;
    for(const auto& edge: markerGraph.edges) {
        if(edge.isWeak) {
            ++weakEdgeCount;
        }
    }
    cout << "Marked as weak " << weakEdgeCount << " marker graph edges out of ";
    cout << markerGraph.edges.size() << " total." << endl;

    cout << "The marker graph has " << globalMarkerGraphVertices.size() << " vertices and ";
    cout << markerGraph.edges.size()-weakEdgeCount << " strong edges." << endl;

    cout << timestamp << "Done flagging weak edges of the marker graph." << endl;
}
#endif



// Find weak edges in the marker graph.
// This new version uses approximate transitive reduction.
void Assembler::flagMarkerGraphWeakEdges(
    size_t lowCoverageThreshold,
    size_t highCoverageThreshold,
    size_t maxDistance)
{
    // Some shorthands for readability.
    auto& edges = markerGraph.edges;
    using VertexId = GlobalMarkerGraphVertexId;
    using EdgeId = GlobalMarkerGraphEdgeId;
    using Edge = MarkerGraph::Edge;

    // Initial message.
    cout << timestamp << "Flagging weak edges of the marker graph "
        "via approximate transitive reduction." << endl;
    cout << "The marker graph has " << globalMarkerGraphVertices.size() << " vertices and ";
    cout << edges.size() << " edges." << endl;

    // Initially flag all edges as not removed by transitive reduction.
    for(auto& edge: edges) {
        edge.wasRemovedByTransitiveReduction = 0;
    }

    // Gather edges for each coverage less than highCoverageThreshold.
    MemoryMapped::VectorOfVectors<EdgeId, EdgeId>  edgesByCoverage;
    edgesByCoverage.createNew(
            largeDataName("tmp-flagMarkerGraphWeakEdges-edgesByCoverage"),
            largeDataPageSize);
    edgesByCoverage.beginPass1(highCoverageThreshold);
    for(EdgeId edgeId=0; edgeId!=edges.size(); edgeId++) {
        const MarkerGraph::Edge& edge = edges[edgeId];
        if(edge.coverage < highCoverageThreshold) {
            edgesByCoverage.incrementCount(edge.coverage);
        }
    }
    edgesByCoverage.beginPass2();
    for(EdgeId edgeId=0; edgeId!=edges.size(); edgeId++) {
        const MarkerGraph::Edge& edge = edges[edgeId];
        if(edge.coverage < highCoverageThreshold) {
            edgesByCoverage.store(edge.coverage, edgeId);
        }
    }
    edgesByCoverage.endPass2();

    // Check that there are no edges with coverage 0.
    CZI_ASSERT(edgesByCoverage[0].size() == 0);

#if 0
    // Vector to store which edges should be flagged as
    // weak at each iteration over coverage.
    MemoryMapped::Vector<bool> edgeFlags;
    edgeFlags.createNew(
        largeDataName("tmp-flagMarkerGraphWeakEdges-edgeFlags"),
        largeDataPageSize);
    edgeFlags.resize(edges.size());
    fill(edgeFlags.begin(), edgeFlags.end(), false);
#endif

    // Vector to contain vertex distances during each BFS.
    // Is is set to -1 fore vertices nt reached by the BFS.
    MemoryMapped::Vector<int> vertexDistances;
    vertexDistances.createNew(
        largeDataName("tmp-flagMarkerGraphWeakEdges-vertexDistances"),
        largeDataPageSize);
    vertexDistances.resize(globalMarkerGraphVertices.size());
    fill(vertexDistances.begin(), vertexDistances.end(), -1);

    // Queue to be used for all BFSs.
    std::queue<VertexId> q;

    // Vector to store vertices encountered duirng a BFS.
    vector<VertexId> bfsVertices;



    // Flag as weak all edges with coverage <= lowCoverageThreshold
    for(size_t coverage=1; coverage<=lowCoverageThreshold; coverage++) {
        const auto& edgesWithThisCoverage = edgesByCoverage[coverage];
        cout << timestamp << "Flagging as weak " << edgesWithThisCoverage.size() <<
            " edges with coverage " << coverage << "." << endl;
        for(const EdgeId edgeId: edgesWithThisCoverage) {
            edges[edgeId].wasRemovedByTransitiveReduction = 1;
        }
    }



    // Process edges of intermediate coverage.
    for(size_t coverage=lowCoverageThreshold+1;
        coverage<highCoverageThreshold; coverage++) {
        const auto& edgesWithThisCoverage = edgesByCoverage[coverage];
        cout << timestamp << "Processing " << edgesWithThisCoverage.size() <<
            " edges with coverage " << coverage << "." << endl;
        size_t count = 0;

        // Loop over edges with this coverage.
        for(const EdgeId edgeId: edgesWithThisCoverage) {
            const Edge& edge = edges[edgeId];
            CZI_ASSERT(!edge.wasRemovedByTransitiveReduction);
            const VertexId u0 = edge.source;
            const VertexId u1 = edge.target;

            // Do a forward BFS starting at v0, up to distance maxDistance,
            // using only edges currently marked as strong
            // and without using this edge.
            // If we encounter v1, v1 is reachable from v0 without
            // using this edge, and so we can mark this edge as weak.
            q.push(u0);
            vertexDistances[u0] = 0;
            bfsVertices.push_back(u0);
            bool found = false;
            while(!q.empty()) {
                const VertexId v0 = q.front();
                q.pop();
                const int distance0 = vertexDistances[v0];
                const int distance1 = distance0 + 1;
                for(const auto edgeId01: markerGraph.edgesBySource[v0]) {
                    if(edgeId01 == edgeId) {
                        continue;
                    }
                    const Edge& edge01 = markerGraph.edges[edgeId01];
                    if(edge01.wasRemovedByTransitiveReduction) {
                        continue;
                    }
                    const VertexId v1 = edge01.target;
                    if(vertexDistances[v1] >= 0) {
                        continue;   // We already encountered this vertex.
                    }
                    if(v1 == u1) {
                        // We found it!
                        found = true;
                        break;
                    }
                    vertexDistances[v1] = distance1;
                    bfsVertices.push_back(v1);
                    if(distance1 < int(maxDistance)) {
                        q.push(v1);
                    }
                }
                if(found) {
                    break;
                }
            }

            if(found) {
                edges[edgeId].wasRemovedByTransitiveReduction = 1;
                ++count;
            }

            // Clean up to be ready to process the next edge.
            while(!q.empty()) {
                q.pop();
            }
            for(const VertexId v: bfsVertices) {
                vertexDistances[v] = -1;
            }
            bfsVertices.clear();
        }

#if 0
        // Actually flag these weak edges in the marker graph.
        size_t count = 0;
        for(EdgeId edgeId=0; edgeId<edgeFlags.size(); edgeId++) {
            if(edgeFlags[edgeId]) {
                edges[edgeId].isWeak = 1;
                edgeFlags[edgeId] = false;
                ++count;
            }
        }
#endif
        cout << "Flagged as weak " << count <<
            " edges with coverage " << coverage <<
            " out of "<< edgesWithThisCoverage.size() << " total." << endl;
    }


    // Clean up our work areas.
    edgesByCoverage.remove();
    // edgeFlags.remove();
    vertexDistances.remove();



    // Count the number of edges that were flagged as weak.
    uint64_t weakEdgeCount = 0;;
    for(const auto& edge: markerGraph.edges) {
        if(edge.wasRemovedByTransitiveReduction) {
            ++weakEdgeCount;
        }
    }
    cout << "Flagged as weak " << weakEdgeCount << " marker graph edges out of ";
    cout << markerGraph.edges.size() << " total." << endl;

    cout << "The marker graph has " << globalMarkerGraphVertices.size() << " vertices and ";
    cout << markerGraph.edges.size()-weakEdgeCount << " strong edges." << endl;

    cout << timestamp << "Done flagging weak edges of the marker graph." << endl;
}



// Return true if an edge disconnects the local subgraph.
bool Assembler::markerGraphEdgeDisconnectsLocalStrongSubgraph(
    GlobalMarkerGraphEdgeId startEdgeId,
    size_t maxDistance,

    // Each of these two must be sized maxDistance.
    array<vector< vector<GlobalMarkerGraphEdgeId> >, 2>& verticesByDistance,

    // Each of these two must be sized globalMarkerGraphVertices.size()
    // and set to all false on entry.
    // It is left set to all false on exit, so it can be reused.
    array<vector<bool>, 2>& vertexFlags
    ) const
{

    // Some shorthands for clarity.
    auto& edges = markerGraph.edges;
    using VertexId = GlobalMarkerGraphVertexId;
    using EdgeId = GlobalMarkerGraphEdgeId;
    using Edge = MarkerGraph::Edge;

    // Check that the work areas are sized as expected.
    for(size_t i=0; i<2; i++) {
        CZI_ASSERT(verticesByDistance[i].size() == maxDistance+1);
        CZI_ASSERT(vertexFlags[i].size() == globalMarkerGraphVertices.size());
    }

    // Find the two vertices of the starting edge.
    const Edge& startEdge = edges[startEdgeId];
    const array<VertexId, 2> startVertexIds = {startEdge.source, startEdge.target};

    // We want to the two BFS one step at a time,
    // going up by one in distance at each step,
    // and avoiding the start edge.
    // So, instead of the usual queue, we store the vertices found at each distance
    // for each of the two start vertices.
    // verticesByDistance is indexed by [0 or 1][distance].
    for(size_t i=0; i<2; i++) {
        CZI_ASSERT(verticesByDistance[i][0].size() == 0);
        verticesByDistance[i][0].clear();
        const VertexId startVertexId = startVertexIds[i];
        verticesByDistance[i][0].push_back(startVertexId);
        vertexFlags[i][startVertexId] = true;
    }
    // cout << "Working on edge " << startVertexIds[0] << " " << startVertexIds[1] << endl;



    // Do the two BFSs in order of increasing distance
    // and avoiding the start edge.
    // At each step we process vertices at distance
    // and find vertices at distance+1.
    bool disconnects = true;
    for(size_t distance=0; distance<maxDistance; distance++) {
        // cout << "Working on distance " << distance << endl;
        for(size_t i=0; i<2; i++) {
            // cout << "Working on i = " << i << endl;
            CZI_ASSERT(verticesByDistance[i][distance+1].size() == 0);
            for(const VertexId vertexId0: verticesByDistance[i][distance]) {
                // cout << "VertexId0 " << vertexId0 << endl;

                // Loop over children.
                auto childEdgeIds = markerGraph.edgesBySource[vertexId0];
                for(EdgeId edgeId: childEdgeIds) {
                    if(edgeId == startEdgeId) {
                        continue;
                    }
                    const Edge& edge = edges[edgeId];
                    if(edge.wasRemovedByTransitiveReduction) {
                        continue;
                    }
                    const VertexId vertexId1 = edge.target;
                    // cout << "Distance " << distance << " i " << i << " found child " << vertexId1 << endl;

                    // If we already found this vertex in this BFS, skip it.
                    if(vertexFlags[i][vertexId1]) {
                        // cout << "Already found." << endl;
                        continue;
                    }

                    // If we already found this vertex in the other BFS,
                    // we have found a path between the vertices of the starting edge
                    // that does not use the starting edge.
                    // Therefore, the starting edge does not disconnect the local subgraph
                    if(vertexFlags[1-i][vertexId1]) {
                        // cout << "Already found on other BFS" << endl;
                        disconnects = false;
                        break;
                    }

                    // This is the first time we find this vertex.
                    verticesByDistance[i][distance+1].push_back(vertexId1);
                    vertexFlags[i][vertexId1] = true;
                }
                if(!disconnects) {
                    break;
                }

                // Loop over parents.
                auto parentEdgeIds = markerGraph.edgesByTarget[vertexId0];
                for(EdgeId edgeId: parentEdgeIds) {
                    if(edgeId == startEdgeId) {
                        continue;
                    }
                    const Edge& edge = edges[edgeId];
                    if(edge.wasRemovedByTransitiveReduction) {
                        continue;
                    }
                    const VertexId vertexId1 = edge.source;
                    // cout << "Distance " << distance << " i " << i << " found parent " << vertexId1 << endl;

                    // If we already found this vertex in this BFS, skip it.
                    if(vertexFlags[i][vertexId1]) {
                        // cout << "Already found." << endl;
                        continue;
                    }

                    // If we already found this vertex in the other BFS,
                    // we have found a path between the vertices of the starting edge
                    // that does not use the starting edge.
                    // Therefore, the starting edge does not disconnect the local subgraph
                    if(vertexFlags[1-i][vertexId1]) {
                        disconnects = false;
                        // cout << "Already found on other BFS" << endl;
                        break;
                    }

                    // This is the first time we find this vertex.
                    verticesByDistance[i][distance+1].push_back(vertexId1);
                    vertexFlags[i][vertexId1] = true;
                }
                if(!disconnects) {
                    break;
                }

            }
            if(!disconnects) {
                break;
            }
        }
        if(!disconnects) {
            break;
        }
    }



    // Clean up.
    for(size_t distance=0; distance<=maxDistance; distance++) {
        for(size_t i=0; i<2; i++) {
            auto& v = verticesByDistance[i][distance];
            for(const VertexId vertexId: v) {
                vertexFlags[i][vertexId] = false;
            }
            v.clear();
        }
    }


    // cout << "Returning " << int(disconnects) << endl;
    // CZI_ASSERT(0);
    return disconnects;
}



// Prune leaves from the strong subgraph of the global marker graph.
void Assembler::pruneMarkerGraphStrongSubgraph(size_t iterationCount)
{
    // Some shorthands.
    using VertexId = GlobalMarkerGraphVertexId;
    using EdgeId = VertexId;

    // Check that we have what we need.
    checkMarkerGraphVerticesAreAvailable();
    checkMarkerGraphEdgesIsOpen();

    // Get the number of edges.
    auto& edges = markerGraph.edges;
    const EdgeId edgeCount = edges.size();

    // Flags to mark edges to prune at each iteration.
    MemoryMapped::Vector<bool> edgesToBePruned;
    edgesToBePruned.createNew(
        largeDataName("tmp-PruneMarkerGraphStrogngSubgraph"),
        largeDataPageSize);
    edgesToBePruned.resize(edgeCount);
    fill(edgesToBePruned.begin(), edgesToBePruned.end(), false);

    // Clear the wasPruned flag of all edges.
    for(MarkerGraph::Edge& edge: edges) {
        edge.wasPruned = 0;
    }



    // At each prune iteration we prune one layer of leaves.
    for(size_t iteration=0; iteration!=iterationCount; iteration++) {
        cout << timestamp << "Begin prune iteration " << iteration << endl;

        // Find the edges to be pruned at each iteration.
        for(EdgeId edgeId=0; edgeId<edgeCount; edgeId++) {
            MarkerGraph::Edge& edge = edges[edgeId];
            if(edge.wasRemovedByTransitiveReduction) {
                continue;
            }
            if(edge.wasPruned) {
                continue;
            }
            if(
                isForwardLeafOfMarkerGraphPrunedStrongSubgraph(edge.target) ||
                isBackwardLeafOfMarkerGraphPrunedStrongSubgraph(edge.source)
                ) {
                edgesToBePruned[edgeId] = true;
            }
        }



        // Flag the edges we found at this iteration.
        EdgeId count = 0;
        for(EdgeId edgeId=0; edgeId<edgeCount; edgeId++) {
            if(edgesToBePruned[edgeId]) {
                edges[edgeId].wasPruned = 1;
                ++count;
                edgesToBePruned[edgeId] = false;    // For next iteration.
            }
        }
        cout << "Pruned " << count << " edges at prune iteration " << iteration << "." << endl;
    }


    edgesToBePruned.remove();


    // Count the number of surviving edges in the pruned strong subgraph.
    size_t count = 0;
    for(MarkerGraph::Edge& edge: edges) {
        if(!edge.wasRemovedByTransitiveReduction && !edge.wasPruned) {
            ++count;
        }
    }
    cout << "The marker graph has " << globalMarkerGraphVertices.size();
    cout << " vertices and " << edgeCount << " edges." << endl;
    cout << "The pruned strong subgraph has " << globalMarkerGraphVertices.size();
    cout << " vertices and " << count << " edges." << endl;
}


// Find out if a vertex is a forward or backward leaf of the pruned
// strong subgraph of the marker graph.
// A forward leaf is a vertex with out-degree 0.
// A backward leaf is a vertex with in-degree 0.
bool Assembler::isForwardLeafOfMarkerGraphPrunedStrongSubgraph(GlobalMarkerGraphVertexId vertexId) const
{
    const auto& forwardEdges = markerGraph.edgesBySource[vertexId];
    for(const auto& edgeId: forwardEdges) {
        const auto& edge = markerGraph.edges[edgeId];
        if(!edge.wasRemovedByTransitiveReduction && !edge.wasPruned) {
            return false;   // We found a forward edge, so this is not a forward leaf.
        }

    }
    return true;    // We did not find any forward edges, so this is a forward leaf.
}
bool Assembler::isBackwardLeafOfMarkerGraphPrunedStrongSubgraph(GlobalMarkerGraphVertexId vertexId) const
{
    const auto& backwardEdges = markerGraph.edgesByTarget[vertexId];
    for(const auto& edgeId: backwardEdges) {
        const auto& edge = markerGraph.edges[edgeId];
        if(!edge.wasRemovedByTransitiveReduction && !edge.wasPruned) {
            return false;   // We found a backward edge, so this is not a backward leaf.
        }

    }
    return true;    // We did not find any backward edges, so this is a backward leaf.
}



// Given an edge of the pruned strong subgraph of the marker graph,
// return the next edge in the linear chain the edge belongs to.
// If the edge is the last edge in its linear chain, return invalidGlobalMarkerGraphEdgeId.
GlobalMarkerGraphEdgeId Assembler::nextEdgeInMarkerGraphPrunedStrongSubgraphChain(
    GlobalMarkerGraphEdgeId edgeId0) const
{
    // Some shorthands.
    using EdgeId = GlobalMarkerGraphEdgeId;
    using Edge = MarkerGraph::Edge;
    const auto& edges = markerGraph.edges;

    // Check that the edge we were passed belongs to the
    // pruned spanning subgraph of the marker graph.
    const Edge& edge0 = edges[edgeId0];
    CZI_ASSERT(!edge0.wasRemoved());

    // If the out-degree and in-degree of the target of this edge are not both 1,
    // this edge is the last of its chain.
    if(
        (markerGraphPrunedStrongSubgraphOutDegree(edge0.target) != 1) ||
        (markerGraphPrunedStrongSubgraphInDegree( edge0.target) != 1)
        ) {
        return invalidGlobalMarkerGraphEdgeId;
    }

    // Loop over all edges following it.
    EdgeId nextEdgeId = invalidGlobalMarkerGraphEdgeId;
    for(const EdgeId edgeId1: markerGraph.edgesBySource[edge0.target]) {
        const Edge& edge1 = edges[edgeId1];

        // Skip the edge if it is not part of the
        // pruned strong subgraph of the marker graph.
        if(edge1.wasRemoved()) {
            continue;
        }

        // Ok, this a possible next edge.
        if(nextEdgeId == invalidGlobalMarkerGraphEdgeId) {
            // This is the first one we find.
            nextEdgeId = edgeId1;
        } else {
            // This is not the first one we found, so the next edge is not unique.
            return invalidGlobalMarkerGraphEdgeId;
        }
    }

    return nextEdgeId;
}



// Given an edge of the pruned strong subgraph of the marker graph,
// return the previous edge in the linear chain the edge belongs to.
// If the edge is the first edge in its linear chain, return invalidGlobalMarkerGraphEdgeId.
GlobalMarkerGraphEdgeId Assembler::previousEdgeInMarkerGraphPrunedStrongSubgraphChain(
    GlobalMarkerGraphEdgeId edgeId0) const
{
    const bool debug = false;
    if(debug) {
        cout << "previousEdgeInMarkerGraphPrunedStrongSubgraphChain begins." << endl;
    }

    // Some shorthands.
    using EdgeId = GlobalMarkerGraphEdgeId;
    using Edge = MarkerGraph::Edge;
    const auto& edges = markerGraph.edges;

    // Check that the edge we were passed belongs to the
    // pruned spanning subgraph of the marker graph.
    const Edge& edge0 = edges[edgeId0];
    CZI_ASSERT(!edge0.wasRemoved());

    // If the out-degree and in-degree of the source of this edge are not both 1,
    // this edge is the last of its chain.
    if(
        (markerGraphPrunedStrongSubgraphOutDegree(edge0.source) != 1) ||
        (markerGraphPrunedStrongSubgraphInDegree( edge0.source) != 1)
        ) {
        return invalidGlobalMarkerGraphEdgeId;
    }

    // Loop over all edges preceding it.
    EdgeId previousEdgeId = invalidGlobalMarkerGraphEdgeId;
    for(const EdgeId edgeId1: markerGraph.edgesByTarget[edge0.source]) {
        const Edge& edge1 = edges[edgeId1];
        if(debug) {
            cout << "Found " << edgeId1 << " " << edge1.source << "->" << edge1.target << endl;
        }

        // Skip the edge if it is not part of the
        // pruned strong subgraph of the marker graph.
        if(edge1.wasRemoved()) {
            if(debug) {
                cout << "Edge was removed." << endl;
            }
            continue;
        }

        // Ok, this a possible previous edge.
        if(previousEdgeId == invalidGlobalMarkerGraphEdgeId) {
            // This is the first one we find.
            if(debug) {
                cout << "Tentative previous edge " << edgeId1 << " " << edge1.source << "->" << edge1.target << endl;
            }
            previousEdgeId = edgeId1;
        } else {
            // This is not the first one we found, so the previous edge is not unique.
            if(debug) {
                cout << "previousEdgeInMarkerGraphPrunedStrongSubgraphChain ends, case 1." << endl;
            }
            return invalidGlobalMarkerGraphEdgeId;
        }
    }
    if(debug) {
        cout << "previousEdgeInMarkerGraphPrunedStrongSubgraphChain ends, case 2 " << previousEdgeId << endl;
    }

    return previousEdgeId;
}



// Return the out-degree or in-degree (number of outgoing/incoming edges)
// of a vertex of the pruned spanning subgraph of the marker graph.
size_t Assembler::markerGraphPrunedStrongSubgraphOutDegree(
    GlobalMarkerGraphVertexId vertexId) const
{
    size_t outDegree = 0;
    for(const auto edgeId: markerGraph.edgesBySource[vertexId]) {
        const auto& edge = markerGraph.edges[edgeId];
        if(!edge.wasRemoved()) {
            ++outDegree;
        }
    }
    return outDegree;
}
size_t Assembler::markerGraphPrunedStrongSubgraphInDegree(
    GlobalMarkerGraphVertexId vertexId) const
{
    size_t inDegree = 0;
    for(const auto edgeId: markerGraph.edgesByTarget[vertexId]) {
        const auto& edge = markerGraph.edges[edgeId];
        if(!edge.wasRemoved()) {
            ++inDegree;
        }
    }
    return inDegree;
}



// Compute consensus sequence for a vertex of the marker graph.
void Assembler::computeMarkerGraphVertexConsensusSequence(
    GlobalMarkerGraphVertexId vertexId,
    vector<Base>& sequence,
    vector<uint32_t>& repeatCounts
    )
{
    // Access the markers of this vertex.
    const MemoryAsContainer<MarkerId> markerIds = globalMarkerGraphVertices[vertexId];
    const size_t markerCount = markerIds.size();
    CZI_ASSERT(markerCount > 0);

    // Find the corresponding oriented read ids, marker ordinals,
    // and marker positions in each oriented read id.
    vector< pair<OrientedReadId, uint32_t> > markerInfos;
    vector<uint32_t> markerPositions;
    markerInfos.reserve(markerIds.size());
    markerPositions.reserve(markerIds.size());
    for(const MarkerId markerId: markerIds) {
        markerInfos.push_back(findMarkerId(markerId));
        markerPositions.push_back(markers.begin()[markerId].position);
    }


    // Loop over all base positions of this marker.
    const size_t k = assemblerInfo->k;
    sequence.resize(k);
    repeatCounts.resize(k);
    for(uint32_t position=0; position<uint32_t(k); position++) {

        // Object to store base and repeat count information
        // at this position.
        Coverage coverage;

        // Loop over markers.
        for(size_t i=0; i<markerCount; i++) {

            // Get the base and repeat count.
            const OrientedReadId orientedReadId = markerInfos[i].first;
            const uint32_t markerPosition = markerPositions[i];
            Base base;
            uint8_t repeatCount;
            tie(base, repeatCount) = getOrientedReadBaseAndRepeatCount(orientedReadId, markerPosition + position);

            // Add it to the Coverage object.
            coverage.addRead(AlignedBase(base), orientedReadId.getStrand(), size_t(repeatCount));
        }

        // Sanity check that all the bases are the same.
        const vector<CoverageData>& coverageData = coverage.getReadCoverageData();
        CZI_ASSERT(coverageData.size() == markerCount);
        const Base firstBase = Base(coverageData.front().base);
        for(const CoverageData& c: coverageData) {
            CZI_ASSERT(Base(c.base) == firstBase);
        }

        // Compute the consensus.
        const Consensus consensus = (*consensusCaller)(coverage);
        sequence[position] = Base(consensus.base);
        repeatCounts[position] = uint32_t(consensus.repeatCount);
    }
}



// Compute consensus sequence for an edge of the marker graph.
// This includes the k bases corresponding to the flanking markers,
// but computed only using reads on this edge.
void Assembler::computeMarkerGraphEdgeConsensusSequenceUsingSeqan(
    GlobalMarkerGraphEdgeId edgeId,
    vector<Base>& sequence,
    vector<uint32_t>& repeatCounts
    )
{

    // Access the markerIntervals for this edge.
    // Each corresponds to an oriented read on this edge.
    const MemoryAsContainer<MarkerInterval> markerIntervals =
        markerGraph.edgeMarkerIntervals[edgeId];
    const size_t markerCount = markerIntervals.size();

    // Initialize a seqan alignment.
    seqan::Align<seqan::String<seqan::Dna> > seqanAlignment;
    static_assert(
        SEQAN_VERSION_MAJOR==2 &&
        SEQAN_VERSION_MINOR==4 &&
        SEQAN_VERSION_PATCH==0,
        "SeqAn version 2.4.0 is required.");
    resize(rows(seqanAlignment), markerCount);



    // Add all the sequences to the seqan alignment,
    // including the flanking markers.
    vector<uint32_t> positions(markerCount);
    for(size_t i=0; i!=markerCount; i++) {
        const MarkerInterval& markerInterval = markerIntervals[i];
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Get the two markers.
        const CompressedMarker& marker0 = orientedReadMarkers[markerInterval.ordinals[0]];
        const CompressedMarker& marker1 = orientedReadMarkers[markerInterval.ordinals[1]];

        // Get the position range, including the flanking markers.
        const uint32_t positionBegin = marker0.position;
        const uint32_t positionEnd = marker1.position + uint32_t(assemblerInfo->k);
        positions[i] = positionBegin;

        // Get the sequence
        string sequenceString;
        for(uint32_t position=positionBegin; position!=positionEnd; position++) {
            sequenceString.push_back(getOrientedReadBase(orientedReadId, position).character());
        }

        // Add it to the seqan alignment.
        seqan::assignSource(seqan::row(seqanAlignment, i), sequenceString);
    }

    // Use seqan to compute the multiple sequence alignment.
    seqan::globalMsaAlignment(seqanAlignment, seqan::Score<int, seqan::Simple>(0, -1, -1));
    // cout << seqanAlignment << endl;

    // The length of the alignment.
    // This includes gaps.
    const size_t alignmentLength = seqan::length(seqan::row(seqanAlignment, 0));



    // Loop over all positions in the seqan alignment.
    // At each position compute a consensus base and repeat count.
    // If the consensus base is not "-", store the base and repeat count.
    sequence.clear();
    repeatCounts.clear();
    for(size_t position=0; position<alignmentLength; position++) {
        // cout << "Alignment position " << position << endl;

        // Create and fill in a Coverage object for this position.
        Coverage coverage;
        for(size_t i=0; i!=markerCount; i++) {
            const MarkerInterval& markerInterval = markerIntervals[i];
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            if(seqan::isGap(seqan::row(seqanAlignment, i), position)) {
                coverage.addRead(AlignedBase::gap(), orientedReadId.getStrand(), 0);
                // cout << "Gap" << endl;
            } else {
                Base base;
                uint8_t repeatCount;
                tie(base, repeatCount) = getOrientedReadBaseAndRepeatCount(orientedReadId, positions[i]);
                ++positions[i];
                coverage.addRead(AlignedBase(base), orientedReadId.getStrand(), repeatCount);
                // cout << base << int(repeatCount) << endl;
            }
        }

        // Compute the consensus at this position.
        const Consensus consensus = (*consensusCaller)(coverage);

        // If not a gap, store the base and repeat count.
        if(!consensus.base.isGap()) {
            sequence.push_back(Base(consensus.base));
            repeatCounts.push_back(uint32_t(consensus.repeatCount));
        }
    }
}



// Compute consensus sequence for an edge of the marker graph.
// This includes the k bases corresponding to the flanking markers,
// but computed only using reads on this edge.
void Assembler::computeMarkerGraphEdgeConsensusSequenceUsingSpoa(
    GlobalMarkerGraphEdgeId edgeId,
    vector<Base>& sequence,
    vector<uint32_t>& repeatCounts
    )
{

    // Access the markerIntervals for this edge.
    // Each corresponds to an oriented read on this edge.
    const MemoryAsContainer<MarkerInterval> markerIntervals =
        markerGraph.edgeMarkerIntervals[edgeId];
    const size_t markerCount = markerIntervals.size();
    CZI_ASSERT(markerCount > 0);

#if 0
    {
        cout << "Working on edge " << edgeId << " ";
        const auto& edge = markerGraph.edges[edgeId];
        cout << int(edge.wasRemovedByTransitiveReduction);
        cout << int(edge.wasPruned);
        cout << int(edge.isBubbleEdge);
        cout << int(edge.replacesBubbleEdges);
        cout << " " << markerCount << endl;
    }
#endif

    // Initialize a spoa alignment.
    const spoa::AlignmentType alignmentType = spoa::AlignmentType::kSW;
    const int8_t match = 1;
    const int8_t mismatch = -1;
    const int8_t gap = -1;
    auto alignmentEngine = spoa::createAlignmentEngine(alignmentType, match, mismatch, gap);
    auto alignmentGraph = spoa::createGraph();


    // Add all the sequences to the alignment,
    // including the flanking markers.
    vector<uint32_t> positions(markerCount);
    string sequenceString;
    vector<string> msa;
    for(size_t i=0; i!=markerCount; i++) {
        const MarkerInterval& markerInterval = markerIntervals[i];
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        // cout << orientedReadId << endl;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Get the two markers.
        const CompressedMarker& marker0 = orientedReadMarkers[markerInterval.ordinals[0]];
        const CompressedMarker& marker1 = orientedReadMarkers[markerInterval.ordinals[1]];

        // Get the position range, including the flanking markers.
        const uint32_t positionBegin = marker0.position;
        const uint32_t positionEnd = marker1.position + uint32_t(assemblerInfo->k);
        // cout << "Position range " << positionBegin << " " << positionEnd << endl;
        positions[i] = positionBegin;

        // Get the sequence.
        sequenceString.clear();
        for(uint32_t position=positionBegin; position!=positionEnd; position++) {
            sequenceString.push_back(getOrientedReadBase(orientedReadId, position).character());
        }

        // Add it to the alignment.
        auto alignment = alignmentEngine->align_sequence_with_graph(sequenceString, alignmentGraph);
        alignmentGraph->add_alignment(alignment, sequenceString);
    }

    // Use spoa to compute the multiple sequence alignment.
    msa.clear();
    alignmentGraph->generate_multiple_sequence_alignment(msa);

    // The length of the alignment.
    // This includes gaps.
    const size_t alignmentLength = msa.front().size();



    // Loop over all positions in the alignment.
    // At each position compute a consensus base and repeat count.
    // If the consensus base is not "-", store the base and repeat count.
    sequence.clear();
    repeatCounts.clear();
    for(size_t position=0; position<alignmentLength; position++) {

        // Create and fill in a Coverage object for this position.
        Coverage coverage;
        for(size_t i=0; i!=markerCount; i++) {
            const MarkerInterval& markerInterval = markerIntervals[i];
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;
            if(msa[i][position] == '-') {
                coverage.addRead(AlignedBase::gap(), orientedReadId.getStrand(), 0);
            } else {
                Base base;
                uint8_t repeatCount;
                tie(base, repeatCount) = getOrientedReadBaseAndRepeatCount(orientedReadId, positions[i]);
                ++positions[i];
                coverage.addRead(AlignedBase(base), orientedReadId.getStrand(), repeatCount);
            }
        }

        // Compute the consensus at this position.
        const Consensus consensus = (*consensusCaller)(coverage);

        // If not a gap, store the base and repeat count.
        if(!consensus.base.isGap()) {
            sequence.push_back(Base(consensus.base));
            repeatCounts.push_back(uint32_t(consensus.repeatCount));
        }
    }
}



// Remove short cycles from the marker graph.
// The argument is the maximum length (number of edges)
// of a cycle path to be considered for removal.
// For now this only processes:
// - Self-edges of the assembly graph.
// - Reversed edges of the assembly graphs (pairs v0->v1, v1->v0).
void Assembler::removeShortMarkerGraphCycles(size_t maxLength)
{
    // To facilitate locating the cycles, create a temporary assembly graph.
    cout << timestamp << "Creating a temporary assembly graph for short cycle removal." << endl;
    createAssemblyGraphEdges();
    createAssemblyGraphVertices();
    assemblyGraph.writeGraphviz("AssemblyGraph-0.dot");
    cout << timestamp << "Done creating a temporary assembly graph for short cycle removal." << endl;

    // The assembly graph edges we want to remove.
    vector<AssemblyGraph::EdgeId> edgesToBeRemoved;



    // Flag self-edges of the assembly graph.
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[edgeId];

        // If not a self-edge, skip.
        if(edge.source != edge.target) {
            continue;
        }

        // If too long, skip.
        const MemoryAsContainer<GlobalMarkerGraphEdgeId> markerGrapEdgeIds =
            assemblyGraph.edgeLists[edgeId];
        if(markerGrapEdgeIds.size() > maxLength) {
            continue;
        }

        // Add it to our list.
        edgesToBeRemoved.push_back(edgeId);

    }



    // Flag reversed edges of the assembly graph.
    // We flag an edge v0->v1 if:
    // - It corresponds to maxLength or fewer edges of the marker graph.
    // - The edge v1->v0 also exists.
    // - v0 has no incoming edges other than v1->v0.
    // - v1 has not outgoing edges other than v1->v0.
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {

        // If it corresponds to maxLength or fewer edges of the marker graph, skip it.
        if(assemblyGraph.edgeLists[edgeId].size() > maxLength) {
            continue;
        }

        // Get the vertices.
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[edgeId];
        const AssemblyGraph::VertexId v0 = edge.source;
        const AssemblyGraph::VertexId v1 = edge.target;

        // If v0 has incoming edges other than v1->v0, skip it.
        if(assemblyGraph.edgesByTarget.size(v0) != 1) {
            continue;
        }

        // If v1 has outgoing edges other than v1->v0, skip it.
        if(assemblyGraph.edgesBySource.size(v1) != 1) {
            continue;
        }

        // Check that the only outgoing edge of v1 is to v0.
        const AssemblyGraph::EdgeId edgeId10 = assemblyGraph.edgesBySource[v1][0];
        const AssemblyGraph::Edge& edge10 = assemblyGraph.edges[edgeId10];
        if(edge10.target != v0) {
            continue;
        }

        // Add it to our list.
        edgesToBeRemoved.push_back(edgeId);
    }



    // Deduplicate the list of edges to be removed.
    sort(edgesToBeRemoved.begin(), edgesToBeRemoved.end());
    edgesToBeRemoved.resize(
        std::unique(edgesToBeRemoved.begin(), edgesToBeRemoved.end()) -
        edgesToBeRemoved.begin());



    // Mark edges of the marker graph.
    size_t markerGraphRemovedEdgeCount = 0;
    for(AssemblyGraph::EdgeId edgeId: edgesToBeRemoved) {
        const MemoryAsContainer<GlobalMarkerGraphEdgeId> markerGrapEdgeIds =
            assemblyGraph.edgeLists[edgeId];
        for(GlobalMarkerGraphEdgeId markerGraphEdgeId: markerGrapEdgeIds) {
            markerGraph.edges[markerGraphEdgeId].isShortCycleEdge = 1;
            ++markerGraphRemovedEdgeCount;
        }
    }
    cout << "Removed " << edgesToBeRemoved.size() <<
        " short cycle edges from the assembly graph\n"
        "and the corresponding " << markerGraphRemovedEdgeCount <<
        " marker graph edges." << endl;


    // Remove the temporary assembly graph.
    // We will create the final one later, after bubble removal.
    assemblyGraph.remove();

}



// Remove short bubbles from the marker graph.
// The argument is the maximum length (number of edges)
// of a bubble branch to be considered for removal.
void Assembler::removeMarkerGraphBubbles(size_t maxLength)
{

    // To facilitate locating the bubbles, create a temporary assembly graph.
    cout << timestamp << "Creating a temporary assembly graph for bubble removal." << endl;
    createAssemblyGraphEdges();
    createAssemblyGraphVertices();
    assemblyGraph.writeGraphviz("AssemblyGraph-1.dot");
    cout << timestamp << "Done creating a temporary assembly graph for bubble removal." << endl;



    // Vector to contain the new marker graph edges to be created.
    vector< pair<GlobalMarkerGraphVertexId, GlobalMarkerGraphVertexId> > newEdges;


    // To find bubbles, loop over assembly graph vertices with out-degree > 1.
    size_t bubbleCount = 0;
    size_t bubbleMarkerGraphEdgeCount = 0;
    for(AssemblyGraph::VertexId v0=0; v0<assemblyGraph.vertices.size(); v0++) {

        // Get the outgoing edges of v0.
        const MemoryAsContainer<AssemblyGraph::EdgeId> edgeIds01 =
            assemblyGraph.edgesBySource[v0];

        // If 0 or 1 outgoing edges, this vertex cannot be the start of a bubble. Skip it.
        if(edgeIds01.size() <= 1) {
            continue;
        }



        // All of the child edges must have the same child v1Common.
        bool isBubble = true;
        AssemblyGraph::VertexId v1Common = AssemblyGraph::invalidVertexId;
        for(AssemblyGraph::EdgeId edgeId01: edgeIds01) {
            const AssemblyGraph::Edge& edge01 = assemblyGraph.edges[edgeId01];
            const AssemblyGraph::VertexId v1 = edge01.target;

            // Check that the child vertex is the same.
            if(v1Common == AssemblyGraph::invalidVertexId) {
                // This is the first child. Store it as v1Common.
                v1Common = v1;
            } else {
                // This is not the first child.
                if(v1 != v1Common) {
                    isBubble = false;
                    break;
                }
            }

            // If this assembly graph edge corresponds to too many marker graph edges,
            // this bubble is too big and so we skip it.
            if(assemblyGraph.edgeLists[edgeId01].size() > maxLength) {
                isBubble = false;
                break;
            }
        }

        // If we found that a bubble does not start at v0, skip.
        if(!isBubble) {
            continue;
        }

        // We also have to check that the in-degree of v1Common is the
        // same as the out-degree of v0.
        if(edgeIds01.size() != assemblyGraph.edgesByTarget[v1Common].size()) {
            continue;
        }

        // Also skip if v0 and v1 coincide.
        if(v0 == v1Common) {
            continue;
        }


        // If getting here, we found a bubble that begins at v0 and ands at v2.
        CZI_ASSERT(v1Common != AssemblyGraph::invalidVertexId);
        ++bubbleCount;
        // const size_t degree = edgeIds01.size();
        // cout << "Bubble between assembly graph vertices " << v0 << " " << v2;
        // cout << " degree " << degree << endl;


        // Mark edges of the marker graph that are internal to this bubble.
        for(AssemblyGraph::EdgeId edgeId01: edgeIds01) {
            const MemoryAsContainer<GlobalMarkerGraphEdgeId> markerGraphEdgeIds =
                assemblyGraph.edgeLists[edgeId01];
            bubbleMarkerGraphEdgeCount += markerGraphEdgeIds.size();
            for(GlobalMarkerGraphEdgeId markerGraphEdgeId: markerGraphEdgeIds) {
                markerGraph.edges[markerGraphEdgeId].isBubbleEdge = 1;
            }
        }

        // The edges we just removed will be replaced with a single edge
        // joining the marker graph vertices corresponding to v0 and v1Common.
        newEdges.push_back(make_pair(
            assemblyGraph.vertices[v0],
            assemblyGraph.vertices[v1Common]));

        // cout << "Source and target marker graph vertices: " <<
        //     assemblyGraph.vertices[v0] << " " <<
        //     assemblyGraph.vertices[v1Common] << endl;
    }
    cout << "Found " << bubbleCount <<" bubbles to be removed." << endl;
    cout << "Found " << bubbleMarkerGraphEdgeCount <<
        " marker graph edges internals to bubbles to be removed." << endl;



    // Remove the temporary assembly graph.
    // We will create the final one later, after bubble removal.
    assemblyGraph.remove();



    // Add the new edges to the marker graph.
    vector<MarkerInterval> markerIntervals;
    for(const pair<GlobalMarkerGraphVertexId, GlobalMarkerGraphVertexId>& p: newEdges) {
        const GlobalMarkerGraphVertexId v0 = p.first;
        const GlobalMarkerGraphVertexId v1 = p.second;
        createBubbleReplacementEdge(v0, v1, false, markerIntervals);
    }
    markerGraph.edgesBySource.remove();
    markerGraph.edgesByTarget.remove();
    createMarkerGraphEdgesBySourceAndTarget(std::thread::hardware_concurrency());
}



// Remove short superbubbles from the marker graph.
// This uses the following process:
// - Create a temporary assembly graph.
// - Compute connected components of the temporary assembly graph,
//   considering only edges with length (number of corresponding
//   marker graph edges) up to maxLength.
// - Remove all edges internal to each connected component.
// - Create a new edge for each entry point/exit point
// conmbination for each connected component.
void Assembler::removeMarkerGraphSuperBubbles(size_t maxLength)
{
    // Create a temporary assembly graph.
    cout << timestamp << "Creating a temporary assembly graph for superbubble removal." << endl;
    createAssemblyGraphEdges();
    createAssemblyGraphVertices();
    assemblyGraph.writeGraphviz("AssemblyGraph-2.dot");
    cout << timestamp << "Done creating a temporary assembly graph for superbubble removal." << endl;


    // Compute connected components of the temporary assembly graph,
    // considering only edges with length (number of corresponding
    // marker graph edges) up to maxLength.
    // Initialize the disjoint set data structures.
    const size_t n = assemblyGraph.vertices.size();
    vector<AssemblyGraph::VertexId> rank(n);
    vector<AssemblyGraph::VertexId> parent(n);
    boost::disjoint_sets<AssemblyGraph::VertexId*, AssemblyGraph::VertexId*> disjointSets(&rank[0], &parent[0]);
    for(size_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }
    for(AssemblyGraph::EdgeId assemblyGraphEdgeId=0; assemblyGraphEdgeId<assemblyGraph.edges.size(); assemblyGraphEdgeId++) {
        if(assemblyGraph.edgeLists[assemblyGraphEdgeId].size() > maxLength) {
            continue;
        }
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[assemblyGraphEdgeId];
        disjointSets.union_set(edge.source, edge.target);
    }



    // Remove all marker graph edges corresponding to assembly graph edges
    // internal to these connected components.
    vector<bool> isEntry(n, false);
    vector<bool> isExit(n, false);
    for(AssemblyGraph::EdgeId assemblyGraphEdgeId=0; assemblyGraphEdgeId<assemblyGraph.edges.size(); assemblyGraphEdgeId++) {
        const MemoryAsContainer<GlobalMarkerGraphEdgeId> markerGraphEdgeIds =
            assemblyGraph.edgeLists[assemblyGraphEdgeId];

        const AssemblyGraph::Edge& assemblyGraphEdge = assemblyGraph.edges[assemblyGraphEdgeId];
        const AssemblyGraph::VertexId v0 = assemblyGraphEdge.source;
        const AssemblyGraph::VertexId v1 = assemblyGraphEdge.target;


        if(disjointSets.find_set(v0) == disjointSets.find_set(v1)) {
            if(markerGraphEdgeIds.size() <= maxLength) {
                // The edge is internal to a connected component.
                // Mark all of its marker graph edges as superbubble edges.
                cout << "Assembly graph edge " << v0 << "->" << v1 << " marked for removal " << markerGraphEdgeIds.size() << endl;
                for(GlobalMarkerGraphEdgeId markerGraphEdgeId: markerGraphEdgeIds) {
                    markerGraph.edges[markerGraphEdgeId].isSuperBubbleEdge = 1;
                }
            }
        } else {

            // The edge is across connected components.
            // Mark the two vertices as entries and exits.
            isExit[v0] = true;
            isEntry[v1] = true;
        }
    }


    // Create a new marker graph edge for each entry/exit combination
    // in the same connected component.
    // Vector to contain the new marker graph edges to be created.
    vector< pair<GlobalMarkerGraphVertexId, GlobalMarkerGraphVertexId> > newEdges;
    for(AssemblyGraph::VertexId v0=0; v0<n; v0++) {
        if(!isEntry[v0]) {
            continue;
        }
        const AssemblyGraph::VertexId component0 = disjointSets.find_set(v0);
        for(AssemblyGraph::VertexId v1=0; v1<n; v1++) {
            if(!isExit[v1]) {
                continue;
            }
            const AssemblyGraph::VertexId component1 = disjointSets.find_set(v1);
            if(component0 != component1) {
                continue;
            }
            cout << "Entry/exit pair " << v0 << " " << v1 << endl;
            newEdges.push_back(make_pair(
                assemblyGraph.vertices[v0],
                assemblyGraph.vertices[v1]));
        }
    }


    // Remove the temporary assembly graph.
    // We will create the final one later, after bubble removal.
    assemblyGraph.remove();



    // Add the new edges to the marker graph.
    vector<MarkerInterval> markerIntervals;
    for(const pair<GlobalMarkerGraphVertexId, GlobalMarkerGraphVertexId>& p: newEdges) {
        const GlobalMarkerGraphVertexId v0 = p.first;
        const GlobalMarkerGraphVertexId v1 = p.second;
        createBubbleReplacementEdge(v0, v1, true, markerIntervals);
    }
    markerGraph.edgesBySource.remove();
    markerGraph.edgesByTarget.remove();
    createMarkerGraphEdgesBySourceAndTarget(std::thread::hardware_concurrency());
}



// Used by removeMarkerGraphBubbles and removeMarkerGraphSuperBubbles.
void Assembler::createBubbleReplacementEdge(
    GlobalMarkerGraphVertexId v0,
    GlobalMarkerGraphVertexId v1,
    bool isSuperBubble,
    vector<MarkerInterval>& markerIntervals)
{
    // Access the markers of the two vertices.
    const MemoryAsContainer<MarkerId> markers0 = globalMarkerGraphVertices[v0];
    const MemoryAsContainer<MarkerId> markers1 = globalMarkerGraphVertices[v1];

    // Find the marker intervals for the new edge and store them.
    markerIntervals.clear();
    auto it0 = markers0.begin();
    auto it1 = markers1.begin();
    while(it0!=markers0.end() && it1!=markers1.end()) {
        const MarkerId markerId0 = *it0;
        const MarkerId markerId1 = *it1;
        OrientedReadId orientedReadId0;
        OrientedReadId orientedReadId1;
        uint32_t ordinal0;
        uint32_t ordinal1;
        tie(orientedReadId0, ordinal0) = findMarkerId(markerId0);
        tie(orientedReadId1, ordinal1) = findMarkerId(markerId1);
        if(orientedReadId0 < orientedReadId1) {
            ++it0;
        } else if (orientedReadId1 < orientedReadId0) {
            ++it1;
        } else {
            CZI_ASSERT(orientedReadId0 == orientedReadId1);
            if(ordinal1 > ordinal0) {
                markerIntervals.push_back(
                    MarkerInterval(orientedReadId0, ordinal0, ordinal1));
            }
            ++it0;
            ++it1;
        }

    }
    CZI_ASSERT(!markerIntervals.empty());
    markerGraph.edgeMarkerIntervals.appendVector(markerIntervals);



    // Store the new edge.
    MarkerGraph::Edge edge;
    edge.source = v0;
    edge.target = v1;
    edge.coverage = uint8_t(min(size_t(255), markerIntervals.size()));
    if(isSuperBubble) {
        edge.replacesSuperBubbleEdges = 1;
    } else {
        edge.replacesBubbleEdges = 1;
    }
    markerGraph.edges.push_back(edge);

}

