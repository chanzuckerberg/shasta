#include "Assembler.hpp"
using namespace shasta;


// "Strict" version of createMarkerGraphEdges.
// Differences from createMarkerGraphEdges:
// - Will only create edges in which all contributing oriented reads have
//   exactly the same RLE sequence. If more than one distinct RLE sequence
//   is present, the edge is split into two parallel edges.
// - Enforces minEdgeCoverage and minEdgeCoveragePerStrand.
//   An edge is not generated if the total number of oriented
//   reads on the edge is less than minEdgeCoverage,
//   of it the number of oriented reads on each strand is less
//   than minEdgeCoveragePerStrand.
// - The main loop is written differently - it loops over reads
//   rather than marker graph vertices.
// Because of these strict criteria, this version generates frequent breaks
// in contiguity that must later be fixed by other means.
void Assembler::createMarkerGraphEdgesStrict(
    uint64_t minEdgeCoverage,
    uint64_t minEdgeCoveragePerStrand,
    size_t threadCount)
{
    cout << timestamp << "createMarkerGraphEdgesStrict begins." << endl;

    // Check that we have what we need.
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store the arguments so the threads can see them.
    createMarkerGraphEdgesStrictData.minEdgeCoverage = minEdgeCoverage;
    createMarkerGraphEdgesStrictData.minEdgeCoveragePerStrand = minEdgeCoveragePerStrand;


    // Find marker intervals and gather them by source vertex id.
    createMarkerGraphEdgesStrictData.markerIntervalInfos.createNew(
        largeDataName("tmp-createMarkerGraphEdgesStrictData-MarkerIntervalInfos"),
        largeDataPageSize);
    createMarkerGraphEdgesStrictData.markerIntervalInfos.beginPass1(markerGraph.vertexCount());
    const uint64_t readCount = getReads().readCount();
    uint64_t batchSize = 10;
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::createMarkerGraphEdgesStrictPass1, threadCount);
    createMarkerGraphEdgesStrictData.markerIntervalInfos.beginPass2();
    setupLoadBalancing(readCount, batchSize);
    runThreads(&Assembler::createMarkerGraphEdgesStrictPass2, threadCount);
    const bool check = true;    // Eventually change to false for performance.
    const bool free = true;
    createMarkerGraphEdgesStrictData.markerIntervalInfos.endPass2(check, free);



    // In pass 3 we actually create the edges.
    // Each thread stores what it finds separately.
    createMarkerGraphEdgesStrictData.threadEdges.resize(threadCount);
    createMarkerGraphEdgesStrictData.threadEdgeMarkerIntervals.resize(threadCount);
    batchSize = 100;
    setupLoadBalancing(markerGraph.vertexCount(), batchSize);
    runThreads(&Assembler::createMarkerGraphEdgesStrictPass3, threadCount);
    createMarkerGraphEdgesStrictData.markerIntervalInfos.remove();



    // Combine the edges found by each thread.
    markerGraph.edges.createNew(
            largeDataName("GlobalMarkerGraphEdges"),
            largeDataPageSize);
    markerGraph.edgeMarkerIntervals.createNew(
            largeDataName("GlobalMarkerGraphEdgeMarkerIntervals"),
            largeDataPageSize);
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        auto& thisThreadEdges = *createMarkerGraphEdgesStrictData.threadEdges[threadId];
        auto& thisThreadEdgeMarkerIntervals = *createMarkerGraphEdgesStrictData.threadEdgeMarkerIntervals[threadId];
        SHASTA_ASSERT(thisThreadEdges.size() == thisThreadEdgeMarkerIntervals.size());
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
    markerGraph.edges.unreserve();
    markerGraph.edgeMarkerIntervals.unreserve();

    SHASTA_ASSERT(markerGraph.edges.size() == markerGraph.edgeMarkerIntervals.size());
    cout << timestamp << "Found " << markerGraph.edges.size();
    cout << " edges for " << markerGraph.vertexCount() << " vertices." << endl;



    // Now we need to create edgesBySource and edgesByTarget.
    createMarkerGraphEdgesBySourceAndTarget(threadCount);

    cout << timestamp << "createMarkerGraphEdgesStrict ends." << endl;
}



void Assembler::createMarkerGraphEdgesStrictPass1(size_t threadId)
{
    createMarkerGraphEdgesStrictPass12(threadId, 1);
}
void Assembler::createMarkerGraphEdgesStrictPass2(size_t threadId)
{
    createMarkerGraphEdgesStrictPass12(threadId, 2);
}



void Assembler::createMarkerGraphEdgesStrictPass12(size_t threadId, uint64_t pass)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.    }
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); ++readId) {

            // Loop over strands.
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);
                const OrientedReadId::Int orientedReadIdValue = orientedReadId.getValue();

                // The first MarkerId for this oriented read.
                const MarkerId startMarkerId = markers.begin(orientedReadIdValue) - markers.begin();

                // Loop over markers of this oriented read.
                const span<CompressedMarker>& orientedReadMarkers = markers[orientedReadIdValue];
                const uint32_t invalidOrdinal = std::numeric_limits<uint32_t>::max();
                uint32_t ordinal0 = invalidOrdinal;
                MarkerGraph::VertexId vertexId0 = MarkerGraph::invalidVertexId;
                for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {

                    // Find the vertex corresponding to this marker ordinal.
                    const MarkerId markerId = startMarkerId + ordinal;
                    const MarkerGraph::CompressedVertexId compressedVertexId = markerGraph.vertexTable[markerId];

                    // If not a valid vertex, skip.
                    if(compressedVertexId == MarkerGraph::invalidCompressedVertexId) {
                        continue;
                    }
                    const MarkerGraph::VertexId vertexId = MarkerGraph::VertexId(compressedVertexId);

                    // If this is the first one we found, store it as
                    // ordinal0 and vertexId0.
                    if(ordinal0 == invalidOrdinal) {
                        ordinal0 = ordinal;
                        vertexId0 = vertexId;
                        continue;
                    }

                    // We already have ordinal0. Generate a new marker interval.
                    if(pass == 1) {
                        createMarkerGraphEdgesStrictData.markerIntervalInfos.incrementCountMultithreaded(vertexId0);
                    } else {
                        CreateMarkerGraphEdgesStrictData::MarkerIntervalInfo markerIntervalInfo;
                        markerIntervalInfo.vertexId1 = vertexId;
                        markerIntervalInfo.orientedReadId = orientedReadId;
                        markerIntervalInfo.ordinal0 = ordinal0;
                        markerIntervalInfo.ordinal1 = ordinal;
                        createMarkerGraphEdgesStrictData.markerIntervalInfos.storeMultithreaded(
                            vertexId0, markerIntervalInfo);
                    }

                    // Store this as the previous ordinal and source vertex.
                    ordinal0 = ordinal;
                    vertexId0 = vertexId;

                }
            }
        }
    }
}



void Assembler::createMarkerGraphEdgesStrictPass3(size_t threadId)
{
    using std::shared_ptr;
    using std::make_shared;

    const uint32_t k = uint32_t(assemblerInfo->k);

    // Create the vector to contain the edges found by this thread.
    shared_ptr< MemoryMapped::Vector<MarkerGraph::Edge> > thisThreadEdgesPointer =
        make_shared< MemoryMapped::Vector<MarkerGraph::Edge> >();
    createMarkerGraphEdgesStrictData.threadEdges[threadId] = thisThreadEdgesPointer;
    MemoryMapped::Vector<MarkerGraph::Edge>& thisThreadEdges = *thisThreadEdgesPointer;
    thisThreadEdges.createNew(
            largeDataName("tmp-ThreadGlobalMarkerGraphEdges-" + to_string(threadId)),
            largeDataPageSize);

    // Create the vector to contain the marker intervals for edges found by this thread.
    shared_ptr< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> >
        thisThreadEdgeMarkerIntervalsPointer =
        make_shared< MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t> >();
    createMarkerGraphEdgesStrictData.threadEdgeMarkerIntervals[threadId] = thisThreadEdgeMarkerIntervalsPointer;
    MemoryMapped::VectorOfVectors<MarkerInterval, uint64_t>&
        thisThreadEdgeMarkerIntervals = *thisThreadEdgeMarkerIntervalsPointer;
    thisThreadEdgeMarkerIntervals.createNew(
            largeDataName("tmp-ThreadGlobalMarkerGraphEdgeMarkerIntervals-" + to_string(threadId)),
            largeDataPageSize);

    // Get coverage criteria.
    const uint64_t minEdgeCoverage = createMarkerGraphEdgesStrictData.minEdgeCoverage;
    const uint64_t minEdgeCoveragePerStrand = createMarkerGraphEdgesStrictData.minEdgeCoveragePerStrand;


    // Vector to store information for each streak of marker intervals
    // between the same two vertices.
    vector<CreateMarkerGraphEdgesStrictData::MarkerIntervalInfo3> markerIntervalInfos3;


    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all marker graph vertices assigned to this batch.
        for(MarkerGraph::VertexId vertexId0=begin; vertexId0!=end; ++vertexId0) {

            // Access the MarkerIntervalInfos for this vertex.
            const span<CreateMarkerGraphEdgesStrictData::MarkerIntervalInfo> markerIntervalInfos =
                createMarkerGraphEdgesStrictData.markerIntervalInfos[vertexId0];

            // Fast exit if there are not enough.
            if(markerIntervalInfos.size() < minEdgeCoverage) {
                continue;
            }
            if(markerIntervalInfos.size() == 0) {
                continue;
            }

            // Sort them by vertexId1, orientedReadId, ordinal0, ordinal1.
            sort(markerIntervalInfos.begin(), markerIntervalInfos.end());


            // Find streaks with the same vertexId1.
            for(uint64_t i0=0; i0!=markerIntervalInfos.size(); /* Incremented later */) {
                const MarkerGraph::VertexId vertexId1 = markerIntervalInfos[i0].vertexId1;

                // Find the end of the streak.
                uint64_t i1 = i0;
                while((i1!=markerIntervalInfos.size()) and (markerIntervalInfos[i1].vertexId1 == vertexId1)) {
                    ++i1;
                }
                const span<CreateMarkerGraphEdgesStrictData::MarkerIntervalInfo> streak(
                    markerIntervalInfos.begin() + i0,
                    markerIntervalInfos.begin() + i1);

                // If the streak is too short, skip.
                if(streak.size() < minEdgeCoverage) {
                    i0 = i1;
                    continue;
                }

                // Store information for this streak, including the RLE sequence.
                markerIntervalInfos3.clear();
                for(const CreateMarkerGraphEdgesStrictData::MarkerIntervalInfo& markerIntervalInfo : streak) {
                    CreateMarkerGraphEdgesStrictData::MarkerIntervalInfo3 markerIntervalInfo3;
                    markerIntervalInfo3.markerInterval.orientedReadId = markerIntervalInfo.orientedReadId;
                    const MarkerId startMarkerId = markers.begin(markerIntervalInfo.orientedReadId.getValue()) - markers.begin();
                    markerIntervalInfo3.markerInterval.ordinals[0] = markerIntervalInfo.ordinal0;
                    markerIntervalInfo3.markerInterval.ordinals[1] = markerIntervalInfo.ordinal1;
                    const MarkerId markerId0 = startMarkerId + markerIntervalInfo.ordinal0;
                    const MarkerId markerId1 = startMarkerId + markerIntervalInfo.ordinal1;
                    const uint32_t position0 = uint32_t(markers.begin()[markerId0].position);
                    const uint32_t position1 = uint32_t(markers.begin()[markerId1].position);
                    if(position1 <= position0 + k) {
                        // Store the overlap.
                        markerIntervalInfo3.overlap = (position0 + k) - position1;
                    } else {
                        // Store the sequence in between.
                        markerIntervalInfo3.overlap = 0;
                        for(uint32_t position = position0 + k; position < position1; position++) {
                            markerIntervalInfo3.sequence.push_back(getReads().getOrientedReadBase(
                                markerIntervalInfo.orientedReadId, position));
                        }
                    }
                    markerIntervalInfos3.push_back(markerIntervalInfo3);
                }

                // Sort again, this time sorting by sequence first.
                sort(markerIntervalInfos3.begin(), markerIntervalInfos3.end());


                // Now we look in the markerIntervalInfos3 vector for streaks
                // with the same sequence.
                for(uint64_t j0=0; j0!=markerIntervalInfos3.size(); /* Incremented later */) {
                    uint64_t j1 = j0;
                    while((j1 != markerIntervalInfos3.size()) and
                        (markerIntervalInfos3[j1].overlap == markerIntervalInfos3[j0].overlap) and
                        (markerIntervalInfos3[j1].sequence == markerIntervalInfos3[j0].sequence)) {
                        ++j1;
                    }
                    // This streak is a candidate edge.
                    const span<CreateMarkerGraphEdgesStrictData::MarkerIntervalInfo3> candidateEdge(
                        &markerIntervalInfos3[j0],
                        &markerIntervalInfos3[j1]);

                    // Only process it we have enough coverage.
                    if(candidateEdge.size() >= minEdgeCoverage) {

                        // Compute coverage by strand.
                        array<uint64_t, 2> strandCoverage = {0, 0};
                        for(const CreateMarkerGraphEdgesStrictData::MarkerIntervalInfo3& markerIntervalInfo3: candidateEdge) {
                            ++strandCoverage[markerIntervalInfo3.markerInterval.orientedReadId.getStrand()];
                        }

                        // Only process it we we have enough coverage on both strands.
                        if( (strandCoverage[0] >= minEdgeCoveragePerStrand) and
                            (strandCoverage[1] >= minEdgeCoveragePerStrand)) {

                            // If getting here, we actually generate an edge.
                            uint64_t coverage = candidateEdge.size();

                            // Store the edge.
                            MarkerGraph::Edge edge;
                            edge.clearFlags();
                            edge.source = vertexId0;
                            edge.target = vertexId1;
                            edge.coverage = (coverage > 255) ? 255 : uint8_t(coverage);
                            thisThreadEdges.push_back(edge);

                            // Store the marker intervals.
                            thisThreadEdgeMarkerIntervals.appendVector();
                            for(const CreateMarkerGraphEdgesStrictData::MarkerIntervalInfo3& markerIntervalInfo3: candidateEdge) {
                                thisThreadEdgeMarkerIntervals.append(markerIntervalInfo3.markerInterval);
                            }
                        }

                    }



                    // Get ready to process the next streak in markerIntervalInfos3.
                    j0 = j1;
                }



                // Get ready to process the next streak.
                i0 = i1;
            }

        }

    }
}


