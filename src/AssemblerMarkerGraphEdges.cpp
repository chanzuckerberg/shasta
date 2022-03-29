#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
#include "performanceLog.hpp"
#include "Reads.hpp"
using namespace shasta;

// Seqan
#include <seqan/align.h>

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>



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
    performanceLog << timestamp << "createMarkerGraphEdgesStrict begins." << endl;

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
    cout << "Found " << markerGraph.edges.size();
    cout << " edges for " << markerGraph.vertexCount() << " vertices." << endl;



    // Now we need to create edgesBySource and edgesByTarget.
    createMarkerGraphEdgesBySourceAndTarget(threadCount);

    performanceLog << timestamp << "createMarkerGraphEdgesStrict ends." << endl;
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



// Write out the sets of parallel marker graph edges.
// Only createMarkerGraphedgesStrict can create parallel edges.
void Assembler::writeParallelMarkerGraphEdges() const
{
    checkMarkerGraphEdgesIsOpen();
    ofstream csv("ParallelMarkerGraphEdges.csv");
    csv << "v0,v1\n";

    // Vector used to store pairs(v1, e01).
    vector< pair<uint64_t, uint64_t> > x;

    // Loop over source vertices v0.
    uint64_t count = 0;
    for(MarkerGraph::VertexId v0=0; v0!=markerGraph.edgesBySource.size(); v0++) {

        // Gather the target vertices and corresponding edge ids.
        const span<const Uint40> edgeIds01 = markerGraph.edgesBySource[v0];
        x.clear();
        for(const Uint40 e01: edgeIds01) {
            const MarkerGraph::Edge& edge01 = markerGraph.edges[e01];
            const uint64_t v1 = uint64_t(edge01.target);
            x.push_back(make_pair(v1, uint64_t(e01)));
        }

        // Sort so all the ones with the same v1 are together.
        sort(x.begin(), x.end());

        // Find streaks with the same v1.
        for(uint64_t i0=0; i0!=x.size(); /* Increment later */) {
            uint64_t i1 = i0;
            while(i1!= x.size() and x[i1].first==x[i0].first) {
                ++i1;
            }
            const uint64_t length = i1 - i0;
            if(length > 1) {
                ++count;
                csv << v0 << "," << x[0].first << ",";
                for(uint64_t i=i0; i!=i1; i++) {
                    csv << x[i].second << ",";
                }
                csv << "\n";
            }

            // Point to the next streak.
            i0 = i1;
        }
    }
    cout << "Found " << count << " sets of parallel marker graph edges." << endl;
}



// Function createMarkerGraphSecondaryEdges can be called after createMarkerGraphEdgesStrict
// to create a minimal amount of additional non-strict edges (secondary edges)
// sufficient to restore contiguity.
// This is a low performance version for testing:
// - It is not multithreaded.
// - It uses standard containers instead of the classes in
//   namespace shasta::MemoryMapped.
void Assembler::createMarkerGraphSecondaryEdges(
    uint32_t secondaryEdgeMaxSkip,
    size_t threadCount)
{
    createMarkerGraphSecondaryEdges(secondaryEdgeMaxSkip, false, threadCount);
    createMarkerGraphSecondaryEdges(secondaryEdgeMaxSkip, true, threadCount);
}
void Assembler::createMarkerGraphSecondaryEdges(
    uint32_t secondaryEdgeMaxSkip,
    bool aggressive,
    size_t threadCount)
{
    using VertexId = MarkerGraph::VertexId;

    // Check that we have what we need.
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    SHASTA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    SHASTA_ASSERT(markerGraph.edges.isOpenWithWriteAccess);
    SHASTA_ASSERT(markerGraph.edgesBySource.isOpen());
    SHASTA_ASSERT(markerGraph.edgesByTarget.isOpen());

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }


    const VertexId vertexCount = markerGraph.vertexCount();
    performanceLog << timestamp << "createMarkerGraphSecondaryEdges begins." << endl;
    cout << "The initial marker graph has " << vertexCount <<
        " vertices and " << markerGraph.edges.size() << " edges." << endl;



    // Vector to contain the dead end pairs.
    // A forward dead end is a vertex with out-degree 0.
    // A backward dead end is a vertex with in-degree 0.
    // Each dead end pair consists of a vertex with zero out-degree
    // (a forward dead end) and its reverse complement,
    // which has zero in-degree and is a backward dead end.
    vector< array<VertexId, 2> > deadEndPairs;

    // Also create a table to contains the forward/backward deadEndId for each vertex.
    // Index by VertexId.
    const uint64_t noDeadEnd = std::numeric_limits<uint64_t>::max();
    vector< array<uint64_t, 2> > deadEndTable(vertexCount, {noDeadEnd, noDeadEnd});

    // Gather the dead end pairs.
    for(VertexId v=0; v!=vertexCount; v++) {
        if(markerGraph.outDegree(v) == 0) {
            const VertexId vRc = markerGraph.reverseComplementVertex[v];
            SHASTA_ASSERT(markerGraph.inDegree(vRc) == 0);
            const uint64_t deadEndPairId = deadEndPairs.size();
            deadEndTable[v][0] = deadEndPairId;
            deadEndTable[vRc][1] = deadEndPairId;
            deadEndPairs.push_back({v, vRc});
        }
    }
    cout << "Found " << deadEndPairs.size() << " dead end pairs." << endl;



    // Look for secondary edges. Loop over forward dead ends.
    vector< array<VertexId, 2> > secondaryEdges;
    vector<MarkerGraph::VertexId> nextVertices;
    vector<MarkerGraph::VertexId> v1Candidates;
    vector<uint64_t> v1Count;
    for(uint64_t deadEndId0=0; deadEndId0<deadEndPairs.size(); deadEndId0++) {
        const VertexId v0 = deadEndPairs[deadEndId0][0];

        // Find the next vertex for each marker / oriented read on this vertex.
        findNextMarkerGraphVertices(v0, secondaryEdgeMaxSkip, nextVertices);

        // Loop over the markers in this vertex.
        // For each marker there is an OrientedReadId and a next vertex.
        v1Candidates.clear();
        for(uint64_t j=0; j<nextVertices.size(); j++) {
            const MarkerGraph::VertexId v1 = nextVertices[j];

            // Disregard secondary edges to self.
            if(v1 == v0) {
                continue;
            }

            // If there is no next vertex for this marker, skip.
            if(v1 == MarkerGraph::invalidVertexId) {
                continue;
            }

            // If not in aggressive mode and this is not a backward dead end, skip.
            // In aggressive mode we accept everything, which gives better
            // contiguity but can also create more artifacts.
            if(not aggressive) {
                const uint64_t deadEnd1 = deadEndTable[v1][1];
                if(deadEnd1 == noDeadEnd) {
                    continue;
                }
            }

            v1Candidates.push_back(v1);
        }
        if(v1Candidates.empty()) {
            continue;
        }

        // Count how many times each possible v1 appeared.
        deduplicateAndCount(v1Candidates, v1Count);
        const uint64_t v1CountMax = *std::max_element(v1Count.begin(), v1Count.end());

        for(uint64_t i=0; i<v1Candidates.size(); i++) {
            if(v1Count[i] == v1CountMax) {
                const VertexId v1 = v1Candidates[i];
                const VertexId v0Rc = markerGraph.reverseComplementVertex[v0];
                const VertexId v1Rc = markerGraph.reverseComplementVertex[v1];
                SHASTA_ASSERT(v0 != v1);
                SHASTA_ASSERT(v0Rc != v1Rc);
                secondaryEdges.push_back({v0, v1});
                secondaryEdges.push_back({v1Rc, v0Rc});
                break;
            }
        }
    }
    cout << "Found " << secondaryEdges.size() << " secondary edges." << endl;
    deduplicate(secondaryEdges);
    cout << "After deduplicating, there are " << secondaryEdges.size() << " secondary edges." << endl;



    // Add the edges we found.
    vector<MarkerInterval> markerIntervals;
    for(const array<VertexId, 2>& secondaryEdge: secondaryEdges) {
        const VertexId v0 = secondaryEdge[0];
        const VertexId v1 = secondaryEdge[1];
        SHASTA_ASSERT(v0 != v1);

        getMarkerIntervals(v0, v1, markerIntervals);

        // Add the edge.
        MarkerGraph::Edge edge;
        edge.source = v0;
        edge.target = v1;
        const uint64_t coverage = markerIntervals.size();
        if(coverage < 256) {
            edge.coverage = uint8_t(coverage);
        } else {
            edge.coverage = 255;
        }
        edge.isSecondary = 1;
        markerGraph.edges.push_back(edge);
        markerGraph.edgeMarkerIntervals.appendVector(markerIntervals);
    }



    // We need to recreate edgesBySource and edgesByTarget.
    if(markerGraph.edgesBySource.isOpen()) {
        markerGraph.edgesBySource.close();
    }
    if(markerGraph.edgesByTarget.isOpen()) {
        markerGraph.edgesByTarget.close();
    }
    createMarkerGraphEdgesBySourceAndTarget(threadCount);

    // We also need to recompute reverse complement marker graph edges.
    if(markerGraph.reverseComplementEdge.isOpen) {
        markerGraph.reverseComplementEdge.close();
    }
    findMarkerGraphReverseComplementEdges(threadCount);


    cout << "After adding secondary edges, the marker graph has " <<
        markerGraph.vertices().size() << " vertices and " <<
        markerGraph.edges.size() << " edges." << endl;


}



// Cluster the oriented reads on a marker graph edge based on their sequence.
// This returns a vector of connected components.
// Each connected component is an index into the marker intervals for the edge.
vector< vector<uint64_t> > Assembler::clusterMarkerGraphEdgeOrientedReads(
    MarkerGraphEdgeId edgeId,
    double errorRateThreshold,
    bool debug) const
{
    // Seqan types and functions used below.
    using TSequence = seqan::String<char>;
    using TAlign = seqan::Align<TSequence, seqan::ArrayGaps>;
    using seqan::assignSource;
    using seqan::globalAlignment;
    using seqan::row;
    using seqan::Score;
    using seqan::Simple;

    // The length of each marker sequence.
    const size_t k = assemblerInfo->k;

    const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
    const span<const MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
    const uint64_t n = markerIntervals.size();
    SHASTA_ASSERT(edge.coverage == n);



    // Get the sequence for each oriented read.
    // This includes sequence for the flanking markers.
    vector<TSequence> sequences;
    for(const MarkerInterval& markerInterval: markerIntervals) {
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Get the two markers.
        const CompressedMarker& marker0 = orientedReadMarkers[markerInterval.ordinals[0]];
        const CompressedMarker& marker1 = orientedReadMarkers[markerInterval.ordinals[1]];

        // Get the position range, including the flanking markers.
        const uint32_t positionBegin = marker0.position;
        const uint32_t positionEnd = marker1.position + uint32_t(k);

        // Get the sequence for this oriented read.
        TSequence sequence;
        for(uint32_t position=positionBegin; position!=positionEnd; ++position) {

            if(assemblerInfo->readRepresentation == 1) {
                Base base;
                uint8_t repeatCount;
                tie(base, repeatCount) = reads->getOrientedReadBaseAndRepeatCount(orientedReadId, position);
                const char c = base.character();
                for(uint64_t j=0; j<repeatCount; j++) {
                    appendValue(sequence, c);
                }
            } else {
                const Base base = reads->getOrientedReadBase(orientedReadId, position);
                const char c = base.character();
                appendValue(sequence, c);
            }
        }
        sequences.push_back(sequence);
    }


    ofstream csv;
    ofstream dot;
    if(debug) {
        csv.open("ClusterOrientedReads.csv");
        dot.open("ClusterOrientedReads.dot");
        dot << "graph G{\n";
        for(uint64_t i=0; i<n; i++) {
            dot << "\"" << markerIntervals[i].orientedReadId << "\";\n";
        }
    }


    // Implicitly create an undirected graph with a vertex
    // for each oriented read in this marker graph edge.
    // If the edit distance between the sequences of two
    // oriented reads is sufficiently small, add an edge
    // betwene the corresponding vertices.
    // Compute connected components of this undirected graph.
    // each connected component corresponds to a cluster.

    // Initialize incremental connected components.
    vector<uint64_t> rank(n);
    vector<uint64_t> parent(n);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    // Loop over pairs of oriented reads.
    for(uint64_t i0=0; i0<n-1; i0++) {
        const TSequence& sequence0 = sequences[i0];

        for(uint64_t i1=i0+1; i1<n; i1++) {
            const TSequence& sequence1 = sequences[i1];

            // Compute the error rate.
            TAlign align;
            resize(rows(align), 2);
            assignSource(row(align,0), sequence0);
            assignSource(row(align,1), sequence1);
            const int editDistance = -globalAlignment(align, Score<int, Simple>(0,-1,-1));
            const double errorRate = double(editDistance) / double(min(length(sequence0), length(sequence1)));

            // If the error rate is small, merge the connected components
            // corresponding to these two oriented reads.
            if(errorRate < errorRateThreshold) {
                disjointSets.union_set(i0, i1);
            }

            if(debug) {
                const OrientedReadId orientedReadId0 = markerIntervals[i0].orientedReadId;
                const OrientedReadId orientedReadId1 = markerIntervals[i1].orientedReadId;
                csv << orientedReadId0 << ",";
                csv << orientedReadId1 << ",";
                csv << errorRate << "\n";
                if((errorRate < errorRateThreshold)) {
                    dot << "\"" << orientedReadId0 << "\"--\"" << orientedReadId1 << "\";\n";
                }
            }
        }
    }

    if(debug) {
        dot << "}\n";
    }


    // Gather connected components.
    vector< vector<uint64_t> > components(n);
    for(uint64_t i=0; i<n; i++) {
        const uint64_t componentId = disjointSets.find_set(i);
        components[componentId].push_back(i);
    }

    // Gather a vector of pairs (componentId, size) for non-empty components.
    vector< pair<uint64_t, uint64_t> > componentTable;
    for(uint64_t componentId=0; componentId<n; componentId++) {
        const uint64_t componentSize = components[componentId].size();
        if(componentSize > 0) {
            componentTable.push_back(make_pair(componentId, componentSize));
        }
    }

    // Sort it by decreasing component size.
    sort(componentTable.begin(), componentTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    // Gather the connected components in order of drcreasing size.
    vector< vector<uint64_t> > sortedComponents;
    for(const auto& p: componentTable) {
        sortedComponents.push_back(components[p.first]);
    }

    return sortedComponents;
}



// Use clusterMarkerGraphEdgeOrientedReads to split secondary marker graph edges
// where necessary.
void Assembler::splitMarkerGraphSecondaryEdges(size_t threadCount)
{

    // CONSTANTS TO BE EXPOSED WHEN CODE STABILIZES.
    const double errorRateThreshold = 0.25;
    const uint64_t minCoverage = 4;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store parameters to all threads can see them.
    auto& data = splitMarkerGraphSecondaryEdgesData;
    data.errorRateThreshold = errorRateThreshold;
    data.minCoverage = minCoverage;

    // Initialize some counts and data structures.
    data.initialSecondaryCount = 0;
    data.splitCount = 0;
    data.createdCount = 0;
    data.threadEdges.resize(threadCount);

    cout << "Total number of marker graph edges before splitting of secondary edge " <<
        markerGraph.edges.size() << endl;

    // The edge splitting is done in multithreaded code.
    const uint64_t batchSize = 1000;
    setupLoadBalancing(markerGraph.edges.size(), batchSize);
    runThreads(&Assembler::splitMarkerGraphSecondaryEdgesThreadFunction, threadCount);

    // Add the new edges to the marker graph.
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        for(const auto& tmpNewEdge: data.threadEdges[threadId]) {
            markerGraph.edges.resize(markerGraph.edges.size() + 1);
            auto& newEdge = markerGraph.edges.back();
            newEdge.source = Uint40(tmpNewEdge.source);
            newEdge.target = Uint40(tmpNewEdge.target);
            newEdge.coverage = uint8_t(tmpNewEdge.markerIntervals.size());
            newEdge.isSecondary = 1;
            markerGraph.edgeMarkerIntervals.appendVector(tmpNewEdge.markerIntervals);
        }
    }
    data.threadEdges.clear();
    data.threadEdges.shrink_to_fit();

    cout << "Total number of marker graph edges after splitting of secondary edges " <<
        markerGraph.edges.size() << endl;
    cout << "Initial number of secondary edges " << data.initialSecondaryCount << endl;
    cout << "Number of secondary edges that were split " << data.splitCount << endl;
    cout << "Number of secondary edges that were created " << data.createdCount << endl;


    // We need to recreate edgesBySource and edgesByTarget.
    if(markerGraph.edgesBySource.isOpen()) {
        markerGraph.edgesBySource.close();
    }
    if(markerGraph.edgesByTarget.isOpen()) {
        markerGraph.edgesByTarget.close();
    }
    createMarkerGraphEdgesBySourceAndTarget(threadCount);

    // We also need to recompute reverse complement marker graph edges.
    if(markerGraph.reverseComplementEdge.isOpen) {
        markerGraph.reverseComplementEdge.close();
    }
    findMarkerGraphReverseComplementEdges(threadCount);

}



void Assembler::splitMarkerGraphSecondaryEdgesThreadFunction(size_t threadId)
{

    auto& data = splitMarkerGraphSecondaryEdgesData;
    const double errorRateThreshold = data.errorRateThreshold;
    const uint64_t minCoverage = data.minCoverage;
    auto& newEdges = data.threadEdges[threadId];

    uint64_t initialSecondaryCount = 0;
    uint64_t splitCount = 0;
    uint64_t createdCount = 0;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker graph edges assigned to this batch
        for(MarkerGraphEdgeId edgeId=begin; edgeId!=end; edgeId++) {
            MarkerGraph::Edge& edge = markerGraph.edges[edgeId];

            // Find the reverse complement edge.
            const MarkerGraphEdgeId edgeIdRc = markerGraph.reverseComplementEdge[edgeId];

            // Process each pair of reverse complemented edges only once.
            SHASTA_ASSERT(edgeIdRc != edgeId);
            if(edgeIdRc < edgeId) {
                continue;
            }
            MarkerGraph::Edge& edgeRc = markerGraph.edges[edgeIdRc];

            // Mark both of them as not removed due to splitting.
            edge.wasRemovedWhileSplittingSecondaryEdges = 0;
            edgeRc.wasRemovedWhileSplittingSecondaryEdges = 0;

            // If these are not secondary edges, do nothing.
            SHASTA_ASSERT(edge.isSecondary == edgeRc.isSecondary);
            if(not edge.isSecondary) {
                continue;
            }
            initialSecondaryCount += 2;

            // Sanity check on the marker intervals.
            const span<MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
            const span<MarkerInterval> markerIntervalsRc = markerGraph.edgeMarkerIntervals[edgeIdRc];
            SHASTA_ASSERT(markerIntervals.size() == markerIntervalsRc.size());
            for(uint64_t i=0; i<markerIntervals.size(); i++) {
                const OrientedReadId& orientedReadId = markerIntervals[i].orientedReadId;
                const OrientedReadId& orientedReadIdRc = markerIntervalsRc[i].orientedReadId;
                SHASTA_ASSERT(orientedReadId.getReadId() == orientedReadIdRc.getReadId());
                SHASTA_ASSERT(orientedReadId.getStrand() == 1 - orientedReadIdRc.getStrand());
            }

            // Cluster oriented reads on this edge based on their sequence.
            // The clusters are returned sorted by decreasing size.
            const vector< vector<uint64_t> > clusters =
                clusterMarkerGraphEdgeOrientedReads(edgeId, errorRateThreshold, false);

            // If we only found one cluster, we don't need to do anything.
            if(clusters.size() == 1) {
                continue;
            }

            // We will split this edge and its reverse complement.
            splitCount += 2;
            edge.wasRemovedWhileSplittingSecondaryEdges = 1;
            edgeRc.wasRemovedWhileSplittingSecondaryEdges = 1;

            // Each cluster generates a new edge,
            // if it is big enough. But always generate at least one.
            for(uint64_t i=0; i<clusters.size(); i++) {
                const vector<uint64_t>& cluster = clusters[i];

                // If this cluster is too small and it is not the first
                // cluster, stop here.
                if((cluster.size() < minCoverage) and (i > 0)) {
                    break;
                }

                // If getting here, generate a new edge for
                // each of edge and edgeRc.
                createdCount += 2;
                newEdges.resize(newEdges.size() + 2);
                auto& newEdge = newEdges[newEdges.size() - 2];
                auto& newEdgeRc = newEdges[newEdges.size() - 1];

                // Fill in the vertices.
                newEdge.source = edge.source;
                newEdge.target = edge.target;
                newEdgeRc.source = edgeRc.source;
                newEdgeRc.target = edgeRc.target;

                // Fill in the marker intervals.
                for(const uint64_t j: cluster) {
                    newEdge.markerIntervals.push_back(markerIntervals[j]);
                    newEdgeRc.markerIntervals.push_back(markerIntervalsRc[j]);
                }
            }
        }
    }

    // Increment global counts.
    __sync_fetch_and_add(&data.initialSecondaryCount, initialSecondaryCount);
    __sync_fetch_and_add(&data.splitCount, splitCount);
    __sync_fetch_and_add(&data.createdCount, createdCount);
}
