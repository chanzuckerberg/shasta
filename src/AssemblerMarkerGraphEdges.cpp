#include "Assembler.hpp"
#include "deduplicate.hpp"
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



// A "local" class used only in Assembler::createMarkerGraphSecondaryEdges below.
namespace shasta {
    namespace CreateMarkerGraphSecondaryEdges {

        class EdgeCandidate {
        public:

            // Information about the first vertex.
            uint32_t deadEndId0;
            uint32_t position0;

            // Information about the second vertex.
            uint32_t deadEndId1;
            uint32_t position1;

            bool operator==(const EdgeCandidate& that) const
            {
                return
                    tie(deadEndId0, position0, deadEndId1, position1) ==
                    tie(that.deadEndId0, that.position0, that.deadEndId1, that.position1);
            }

            bool operator<(const EdgeCandidate& that) const
            {
                return
                    tie(deadEndId0, position0, deadEndId1, position1) <
                    tie(that.deadEndId0, that.position0, that.deadEndId1, that.position1);
            }

            void reverseComplement()
            {
                swap(deadEndId0, deadEndId1);
                swap(position0, position1);
            }

        };
    }
}



// Function createMarkerGraphSecondaryEdges can be called after createMarkerGraphEdgesStrict
// to create a minimal amount of additional non-strict edges (secondary edges)
// sufficient to restore contiguity.
// This is a low performance version for testing:
// - It is not multithreaded.
// - It uses standard containers instead of trhe classes in
//   namespace shasta::MemoryMapped.
void Assembler::createMarkerGraphSecondaryEdges(
    uint64_t minEdgeCoverage,
    uint64_t minEdgeCoveragePerStrand,
    uint64_t neighborhoodSize,
    size_t threadCount)
{
    using VertexId = MarkerGraph::VertexId;
    using EdgeId = MarkerGraph::EdgeId;
    using Edge = MarkerGraph::Edge;
    using namespace CreateMarkerGraphSecondaryEdges;

    const bool debug = false;

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
    cout << timestamp << "createMarkerGraphSecondaryEdges begins." << endl;
    cout << "The initial marker graph has " << vertexCount <<
        " vertices and " << markerGraph.edges.size() << " edges." << endl;



    // Gather the dead ends (0-forward, 1=backward).
    // Go forward/backward to include up to neighborhoodSize
    // for each dead end.
    // A forward dead end has the same numbering as its
    // reverse complemented backward dead end.
    array< vector<vector<VertexId> >, 2> deadEnds;
    for(VertexId vertexId=0; vertexId!=vertexCount; vertexId++) {
        if(markerGraph.outDegree(vertexId) == 0) {
            SHASTA_ASSERT(markerGraph.inDegree(markerGraph.reverseComplementVertex[vertexId]) == 0);
            deadEnds[0].push_back(vector<VertexId>());
            deadEnds[1].push_back(vector<VertexId>());
            VertexId v = vertexId;
            for(uint64_t i=0; i<neighborhoodSize; i++) {
                deadEnds[0].back().push_back(v);
                deadEnds[1].back().push_back(markerGraph.reverseComplementVertex[v]);
                if(markerGraph.inDegree(v) != 1) {
                    break;
                }
                const EdgeId edgeId = EdgeId(markerGraph.edgesByTarget[v][0]);
                const Edge& edge = markerGraph.edges[edgeId];
                v = edge.source;
            }
        }
    }
    cout << "Found " << deadEnds[0].size() << " forward dead ends and " <<
        deadEnds[1].size() << " backward dead ends." << endl;

    if(debug) {
        ofstream csv("DeadEnds.csv");
        csv << "Direction,Id,Size,VertexIds\n";
        for(uint64_t direction=0; direction<2; direction++) {
            const vector<vector<VertexId> >& v = deadEnds[direction];
            for(uint64_t i=0; i<v.size(); i++) {
                const vector<VertexId>& neighborhood = v[i];
                csv << ((direction == 0) ? "Forward," : "Backward,");
                csv << i << ",";
                csv << neighborhood.size() << ",";
                for(const VertexId vertexId: neighborhood) {
                    csv << vertexId << ",";
                }
                csv << "\n";
            }
        }
    }



    // Construct two vectors that for each vertex tell us which dead ends
    // it belongs to, if any.
    // Each vector stores pairs(deadEndId, vertexPositionInDeadEnd).
    const uint32_t noDeadEnd = std::numeric_limits<uint32_t>::max();
    const pair<uint32_t, uint32_t> noDeadEndPair = make_pair(noDeadEnd, noDeadEnd);
    vector< array< pair<uint32_t, uint32_t>, 2> > vertexDeadEnd(vertexCount,
        array< pair<uint32_t, uint32_t>, 2>({noDeadEndPair, noDeadEndPair}));
    for(uint64_t direction=0; direction<2; direction++) {
        const vector<vector<VertexId> >& v = deadEnds[direction];
        for(uint32_t deadEndId=0; deadEndId<uint32_t(v.size()); deadEndId++) {
            const vector<VertexId>& deadEnd = v[deadEndId];
            for(uint32_t i=0; i<uint32_t(deadEnd.size()); i++) {
                const VertexId vertexId = deadEnd[i];
                vertexDeadEnd[vertexId][direction] = make_pair(deadEndId, i);
            }
        }
    }



    // Find edges to add.
    vector<MarkerGraph::VertexId> nextVertices;
    vector<EdgeCandidate> edgeCandidates;

    // Loop over forward dead ends.
    for(uint32_t deadEndId0=0; deadEndId0<uint32_t(deadEnds[0].size()); deadEndId0++) {
        if(debug) {
            cout << "Forward dead end " << deadEndId0 << endl;
        }

        // Coverage map for the candidate edges for this dead end.
        std::map<EdgeCandidate, array<uint64_t, 2> > coverageMap;

        // Loop over vertices of this forward dead end.
        const vector<VertexId>& neighborhood0 = deadEnds[0][deadEndId0];
        for(uint32_t position0=0; position0<uint32_t(neighborhood0.size()); position0++) {
            const VertexId vertexId0 = neighborhood0[position0];

            // Find the next vertex for each oriented read on this vertex.
            findNextMarkerGraphVertices(vertexId0, nextVertices);

            // Each distinct next vertex generates a possible edge to be added.
            // Compute coverage for each strand of this edge candidate.
            const span<MarkerId> markerIds = markerGraph.getVertexMarkerIds(vertexId0);
            SHASTA_ASSERT(nextVertices.size() == markerIds.size());
            for(uint64_t j=0; j<markerIds.size(); j++) {
                const MarkerId markerId = markerIds[j];
                const MarkerGraph::VertexId vertexId1 = nextVertices[j];
                if(vertexId1 == MarkerGraph::invalidVertexId) {
                    continue;
                }
                const uint32_t deadEndId1 = vertexDeadEnd[vertexId1][1].first;
                const uint32_t position1 = vertexDeadEnd[vertexId1][1].second;

                // The next vertex must belong to a backward dead end.
                if(deadEndId1 == noDeadEnd) {
                    continue;
                }

                // The next vertex must not belong to the same forward dead end.
                if(vertexDeadEnd[vertexId1][0].first == deadEndId0) {
                    continue;
                }

                EdgeCandidate edgeCandidate;
                edgeCandidate.deadEndId0 = deadEndId0;
                edgeCandidate.position0 = position0;
                edgeCandidate.deadEndId1 = deadEndId1;
                edgeCandidate.position1 = position1;

                // Increment the coverage map for this edge candidate.
                OrientedReadId orientedReadId;
                tie(orientedReadId, ignore) = findMarkerId(markerId);
                const uint64_t strand = orientedReadId.getStrand();
                auto it = coverageMap.find(edgeCandidate);
                if(it == coverageMap.end()) {
                    tie(it, ignore) = coverageMap.insert(make_pair(
                        edgeCandidate,
                        array<uint64_t, 2>({0, 0})));
                }
                ++it->second[strand];
            }
        }

#if 0
        // Write out the candidate edges we found for this forward dead end.
        if(debug) {
            for(const auto& p: coverageMap) {
                const VertexId vertexId = p.first.first;
                const VertexId nextVertexId = p.first.second;
                const uint64_t coverage0 = p.second[0];
                const uint64_t coverage1 = p.second[1];
                const uint64_t coverage = coverage0 + coverage1;
                cout << vertexId << "->" << nextVertexId << " " <<
                    coverage0 << "+" << coverage1 << "=" << coverage << "\n";
            }
        }
#endif

        if(coverageMap.empty()) {
            if(debug) {
                cout << "No candidate edges available for this forward dead end." << endl;
            }
            continue;
        }



        // In a first step, look for candidate edges that satisfy minEdgeCoverage
        // and minEdgeCoveragePerStrand, without requiring the RLE sequences to be
        // all identical. If there is more than one, choose the one with the best
        // coverage.
        auto itBest = coverageMap.end();
        uint64_t bestCoverage = 0;
        for(auto it=coverageMap.begin(); it!=coverageMap.end(); ++it) {

            // If coverage is insufficient, skip it.
            const array<uint64_t, 2>& strandCoverage = it->second;
            const uint64_t coverage = strandCoverage[0] + strandCoverage[1];
            if(coverage < minEdgeCoverage or
                strandCoverage[0] < minEdgeCoveragePerStrand or
                strandCoverage[1] < minEdgeCoveragePerStrand) {
                continue;
            }

            if(coverage > bestCoverage) {
                bestCoverage = coverage;
                itBest = it;
            }
        }



        if(itBest == coverageMap.end()) {

            // We did not find a candidate edge with sufficient coverage.
            // Pick the edge with best coverage.
            uint64_t bestCoverage = 0;
            for(auto it=coverageMap.begin(); it!=coverageMap.end(); ++it) {

                const array<uint64_t, 2>& strandCoverage = it->second;
                const uint64_t coverage = strandCoverage[0] + strandCoverage[1];

                if(coverage > bestCoverage) {
                    bestCoverage = coverage;
                    itBest = it;
                }
            }

        }

        // Add this edge candidate to our list.
        EdgeCandidate edgeCandidate = itBest->first;
        edgeCandidates.push_back(edgeCandidate);

        // Also add its reverse complement.
        edgeCandidate.reverseComplement();
        edgeCandidates.push_back(edgeCandidate);



#if 0
        if(debug) {
            const uint64_t coverage0 = itBest->second[0];
            const uint64_t coverage1 = itBest->second[1];
            const uint64_t coverage = coverage0 + coverage1;
            cout << "Chose candidate edge " << vertexId << "->" << nextVertexId << " " <<
                coverage0 << "+" << coverage1 << "=" << coverage << "\n";
        }
#endif
    }

    // Remove duplicates.
    deduplicate(edgeCandidates);



    // There is still the possibility of more than one secondary edge
    // between the same two dead ends.
    // Store them in a map, keyed by the two dead ends,
    // with the first numbered lower.
    std::map< pair<uint32_t, uint32_t>, vector<EdgeCandidate> > edgeCandidateMap;
    for(const EdgeCandidate& edgeCandidate: edgeCandidates) {
        uint32_t deadEndId0 = edgeCandidate.deadEndId0;
        uint32_t deadEndId1 = edgeCandidate.deadEndId1;
        if(deadEndId1 < deadEndId0) {
            swap(deadEndId0, deadEndId1);
        }
        edgeCandidateMap[make_pair(deadEndId0, deadEndId1)].push_back(edgeCandidate);
    }



    // Now recreate the edge candidates, keeping only a reverse complemented pair
    // for each pair of dead ends.
    edgeCandidates.clear();
    for(const auto& p: edgeCandidateMap) {
        EdgeCandidate edgeCandidate = p.second.front();
        edgeCandidates.push_back(edgeCandidate);
        edgeCandidate.reverseComplement();
        edgeCandidates.push_back(edgeCandidate);
    }



    // Create the secondary edges.
    vector<MarkerInterval> markerIntervals;
    for(const EdgeCandidate& edgeCandidate: edgeCandidates) {
        const uint32_t deadEndId0 = edgeCandidate.deadEndId0;
        const uint32_t position0 = edgeCandidate.position0;
        const uint32_t deadEndId1 = edgeCandidate.deadEndId1;
        const uint32_t position1 = edgeCandidate.position1;

        const VertexId v0 = deadEnds[0][deadEndId0][position0];
        const VertexId v1 = deadEnds[1][deadEndId1][position1];
        getMarkerIntervals(v0, v1, markerIntervals);

        if(debug) {
            cout << "Adding edge " << markerGraph.edges.size() <<
                " " << v0 << "->" << v1 << "\n";
        }

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


        // Prune edges closer to the dead end.
        // This can't be done because it can break contiguity.
        /*
        for(uint32_t position=0; position<position0; position++) {
            const VertexId vertexId0 = deadEnds[0][deadEndId0][position+1];
            const VertexId vertexId1 = deadEnds[0][deadEndId0][position];
            const EdgeId edgeId = markerGraph.findEdgeId(vertexId0,  vertexId1);
            cout << "Pruned " << vertexId0 << "->" << vertexId1 << endl;
            markerGraph.edges[edgeId].wasPruned = 1;
        }
        for(uint32_t position=0; position<position1; position++) {
            const VertexId vertexId0 = deadEnds[1][deadEndId1][position];
            const VertexId vertexId1 = deadEnds[1][deadEndId1][position+1];
            const EdgeId edgeId = markerGraph.findEdgeId(vertexId0,  vertexId1);
            cout << "Pruned " << vertexId0 << "->" << vertexId1 << endl;
            markerGraph.edges[edgeId].wasPruned = 1;
        }
        */
    }
    cout << "Created " << edgeCandidates.size() << " secondary edges." << endl;
    cout << "After adding secondary edges, the marker graph has " << vertexCount <<
        " vertices and " << markerGraph.edges.size() << " edges." << endl;
    SHASTA_ASSERT(markerGraph.edgeMarkerIntervals.size() == markerGraph.edges.size());


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

    cout << timestamp << "createMarkerGraphSecondaryEdges ends." << endl;

}


