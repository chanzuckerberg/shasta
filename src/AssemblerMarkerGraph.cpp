// Shasta.
#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "LocalMarkerGraph.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Standard library.
#include "chrono.hpp"
#include <queue>



// Loop over all alignments in the read graph
// to create vertices of the global marker graph.
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

    // Number of threads. If zero, a number of threads equal to
    // the number of virtual processors is used.
    size_t threadCount
)
{

    using VertexId = GlobalMarkerGraphVertexId;
    using CompressedVertexId = CompressedGlobalMarkerGraphVertexId;

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



    // Use the disjoint sets to find the vertex of the global
    // marker graph that each oriented marker belongs to.

    // Compute the global marker graph vertex corresponding
    // to each MarkerId.
    cout << timestamp << "Storing the global marker graph vertex each marker belongs to." << endl;
    globalMarkerGraphVertex.createNew(
        largeDataName("GlobalMarkerGraphVertex"),
        largeDataPageSize);
    globalMarkerGraphVertex.reserveAndResize(data.orientedMarkerCount);
    batchSize = 1000000;
    setupLoadBalancing(markers.totalSize(), batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction2, threadCount);

    // Clean up.
    data.disjointSetsPointer = 0;
    data.disjointSetsData.remove();


    // Using the disjoint sets data structure we computed a
    // "raw" vertex id that each marker belongs to.
    // These "raw" vertex ids are in [0, number of markers),
    // but are not contiguous, because the number of vertices is
    // less than the number of markers.
    // To compute final vertex ids that are contiguous
    // in [0, numberOfVertices), we begin by counting
    // the number of markers for each raw vertex id.
    // This is stored in workArea.
    cout << timestamp << "Counting the number of markers in each vertex." << endl;
    data.workArea.createNew(
        largeDataName("tmp-ComputeAllAlignmentsWorkArea"),
        largeDataPageSize);
    data.workArea.reserveAndResize(data.orientedMarkerCount);
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction3, threadCount);



    // Now we can loop over the work area, and replace each
    // entry at least equal to minCoverage with an increasing id
    // which will be the final vertex id.
    // Entries that are 0 are set to maxValueMinus1.
    // Entries that are greater than 0 and less than minCoverage
    // are set to maxValue, so markers that don't belong to
    // any vertex receive a vertex id of all one bits.
    const uint64_t maxValue = std::numeric_limits<uint64_t>::max();
    const uint64_t maxValueMinus1 = maxValue - 1ULL;
    cout << timestamp << "Generating contiguous vertex ids." << endl;
    GlobalMarkerGraphVertexId vertexId = 0;
    for(GlobalMarkerGraphVertexId& w: data.workArea) {
        if(w == 0ULL) {
            w = maxValueMinus1;
        } else if(w <minCoverage) {
            w = maxValue;
        } else {
            w = vertexId++;
        }
    }
    const GlobalMarkerGraphVertexId vertexCount = vertexId;
    cout << timestamp << "The global marker graph has " << vertexCount;
    cout << " vertices for " << data.orientedMarkerCount;
    cout << " oriented markers." << endl;

    // Now we can use the workArea to convert raw vertex ids to final vertex ids.
    cout << timestamp << "Storing final vertex ids for all markers." << endl;
    setupLoadBalancing(data.orientedMarkerCount, batchSize);
    runThreads(&Assembler::createMarkerGraphVerticesThreadFunction4, threadCount);
    data.workArea.remove();



    // Gather the oriented marker ids of each vertex of the global marker graph.
    cout << timestamp << "Gathering the oriented markers of each vertex." << endl;
    globalMarkerGraphVertices.createNew(
        largeDataName("GlobalMarkerGraphVertices"),
        largeDataPageSize);
    cout << timestamp << "... pass 1" << endl;
    globalMarkerGraphVertices.beginPass1(vertexCount);
    const CompressedVertexId maxValue40Bits = maxValue;
    for(const CompressedVertexId& v: globalMarkerGraphVertex) {
        if(v != maxValue40Bits) {
            globalMarkerGraphVertices.incrementCount(v);
        }
    }
    cout << timestamp << "... pass 2" << endl;
    globalMarkerGraphVertices.beginPass2();
    for(VertexId i=0; i<globalMarkerGraphVertex.size(); i++) {
        const CompressedVertexId v = globalMarkerGraphVertex[i];
        if(v != maxValue40Bits) {
            globalMarkerGraphVertices.store(globalMarkerGraphVertex[i], i);
        }
    }
    globalMarkerGraphVertices.endPass2();

    // Sort the markers in each vertex.
    cout << timestamp << "Sorting the oriented markers of each vertex." << endl;
    for(VertexId i=0; i<vertexCount; i++) {
        if(globalMarkerGraphVertices.size(i) > 1) {
            sort(globalMarkerGraphVertices.begin(i), globalMarkerGraphVertices.end(i));
        }
    }


    // Check that all the markers of a vertex have the same kmer id.
    if(true) {
        for(VertexId i=0; i<globalMarkerGraphVertices.size(); i++) {
            const MemoryAsContainer<MarkerId> vertexMarkers = globalMarkerGraphVertices[i];
            KmerId kmerId = 0;
            for(size_t j=0; j<vertexMarkers.size(); j++) {
                const MarkerId markerId = vertexMarkers[j];
                const CompressedMarker& marker = markers.begin()[markerId];
                if(j == 0) {
                    kmerId = marker.kmerId;
                } else {
                    CZI_ASSERT(kmerId == marker.kmerId);
                }
            }
        }
    }

    // Create a histogram of number of markers per vertex.
    cout << "Creating marker graph vertex coverage histogram." << endl;
    vector<uint64_t> histogram;
    for(size_t i=0; i<globalMarkerGraphVertices.size(); i++) {
        const size_t count = globalMarkerGraphVertices.size(i);
        if(histogram.size()<= count) {
            histogram.resize(count+1, 0ULL);
        }
        ++histogram[count];
    }
    ofstream csv("GlobalMarkerGraphHistogram.csv");
    csv << "MarkerCount,Frequency\n";
    for(size_t i=0; i<histogram.size(); i++) {
        const auto n = histogram[i];
        if(n) {
            csv << i << "," << n << "\n";
        }
    }


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

    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(MarkerId i=begin; i!=end; ++i) {
            globalMarkerGraphVertex[i] = disjointSets.find(i);
        }
    }
}



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
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > >& children,
    vector< pair<GlobalMarkerGraphVertexId, MarkerGraphNeighborInfo> >& workArea
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
        MarkerGraphNeighborInfo info;
        tie(info.orientedReadId, info.ordinal0) = findMarkerId(markerId);

        // Find the next marker that is contained in a vertex.
        const auto markerCount = markers.size(info.orientedReadId.getValue());
        for(info.ordinal1=info.ordinal0+1; info.ordinal1<markerCount; ++info.ordinal1) {

            // Find the vertex id.
            const MarkerId childMarkerId =  getMarkerId(info.orientedReadId, info.ordinal1);
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
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > >& parents,
    vector< pair<GlobalMarkerGraphVertexId, MarkerGraphNeighborInfo> >& workArea
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
        MarkerGraphNeighborInfo info;
        tie(info.orientedReadId, info.ordinal0) = findMarkerId(markerId);
        if(info.ordinal0 == 0) {
            continue;
        }

        // Find the previous marker that is contained in a vertex.
        for(info.ordinal1=info.ordinal0-1; ; --info.ordinal1) {

            // Find the vertex id.
            const MarkerId parentMarkerId =  getMarkerId(info.orientedReadId, info.ordinal1);
            const GlobalMarkerGraphVertexId parentVertexId =
                globalMarkerGraphVertex[parentMarkerId];

            // If this marker correspond to a vertex, add it to our list.
            if( parentVertexId!=invalidCompressedGlobalMarkerGraphVertexId &&
                !isBadMarkerGraphVertex(parentVertexId)) {
                workArea.push_back(make_pair(parentVertexId, info));
                break;
            }

            if(info.ordinal1  == 0) {
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
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > > children;
    vector< pair<GlobalMarkerGraphVertexId, MarkerGraphNeighborInfo> > workArea;
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
            info.ordinal0 = childInfo.ordinal0;
            info.ordinal1 = childInfo.ordinal1;

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
    // and empty veector.
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
    CZI_ASSERT(0);
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
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > > children;
    vector< pair<GlobalMarkerGraphVertexId, vector<MarkerGraphNeighborInfo> > > parents;
    vector< pair<GlobalMarkerGraphVertexId, MarkerGraphNeighborInfo> > workArea;
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
                for(const MarkerGraphNeighborInfo& x: v) {
                    markerIntervalVector.push_back(MarkerInterval(
                        x.orientedReadId, x.ordinal0, x.ordinal1));
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
                for(const MarkerGraphNeighborInfo& x: v) {
                    markerIntervalVector.push_back(MarkerInterval(
                        x.orientedReadId, x.ordinal1, x.ordinal0));
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
            for(const MarkerGraphNeighborInfo& x: v) {
                markerIntervalVector.push_back(MarkerInterval(
                    x.orientedReadId, x.ordinal0, x.ordinal1));
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
    LocalMarkerGraph& graph
    )
{
    const GlobalMarkerGraphVertexId startVertexId =
        getGlobalMarkerGraphVertex(orientedReadId, ordinal);
    return extractLocalMarkerGraphUsingStoredConnectivity(startVertexId, distance, timeout, graph);

}



bool Assembler::extractLocalMarkerGraphUsingStoredConnectivity(
    GlobalMarkerGraphVertexId startVertexId,
    int distance,
    double timeout,                 // Or 0 for no timeout.
    LocalMarkerGraph& graph
    )
{
    // Sanity check.
    checkMarkerGraphConnectivityIsOpen();

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
        const auto children = markerGraphConnectivity.edgesBySource[vertexId0];
        for(GlobalMarkerGraphVertexId vertexId1: children) {

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
                markerIntervals.clear();
                getGlobalMarkerGraphEdgeInfo(vertexId0, vertexId1, markerIntervals);
                graph.storeEdgeInfo(e, markerIntervals);
            }
        }

        CZI_ASSERT(0);
    }

    CZI_ASSERT(0);
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



// Compute connectivity of the global marker graph.
// Vertices with more than markerCountOverflow vertices are skipped.
void Assembler::createMarkerGraphConnectivity(
    size_t threadCount,
    size_t markerCountOverflow
    )
{
    cout << timestamp << "createMarkerGraphConnectivity begins." << endl;

    // Check that we have what we need.
    checkMarkerGraphVerticesAreAvailable();

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Store markerCountOverflow so all threads can see it.
    markerGraphConnectivity.markerCountOverflow = markerCountOverflow;

    // Each thread stores the edges it finds in a separate vector.
    markerGraphConnectivity.threadEdges.resize(threadCount);
    cout << timestamp << "Processing " << globalMarkerGraphVertices.size();
    cout << " marker graph vertices." << endl;
    setupLoadBalancing(globalMarkerGraphVertices.size(), 100000);
    runThreads(&Assembler::createMarkerGraphConnectivityThreadFunction0, threadCount,
        "threadLogs/createMarkerGraphConnectivity0");

    // Combine the edges found by each thread.
    cout << timestamp << "Combining the edges found by each thread." << endl;
    markerGraphConnectivity.edges.createNew(
            largeDataName("GlobalMarkerGraphEdges"),
            largeDataPageSize);
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        auto& thisThreadEdges = *markerGraphConnectivity.threadEdges[threadId];
        for(const auto& edge: thisThreadEdges) {
            markerGraphConnectivity.edges.push_back(edge);
        }
        thisThreadEdges.remove();
    }
    cout << timestamp << "Found " << markerGraphConnectivity.edges.size();
    cout << " edges for " << globalMarkerGraphVertices.size() << " vertices." << endl;



    // Now we need to create edgesBySource and edgesByTarget.
    markerGraphConnectivity.edgesBySource.createNew(
        largeDataName("GlobalMarkerGraphEdgesBySource"),
        largeDataPageSize);
    markerGraphConnectivity.edgesByTarget.createNew(
        largeDataName("GlobalMarkerGraphEdgesByTarget"),
        largeDataPageSize);

    cout << timestamp << "Creating connectivity: pass 1 begins." << endl;
    markerGraphConnectivity.edgesBySource.beginPass1(globalMarkerGraphVertices.size());
    markerGraphConnectivity.edgesByTarget.beginPass1(globalMarkerGraphVertices.size());
    setupLoadBalancing(markerGraphConnectivity.edges.size(), 100000);
    runThreads(&Assembler::createMarkerGraphConnectivityThreadFunction1, threadCount);

    cout << timestamp << "Creating connectivity: pass 2 begins." << endl;
    markerGraphConnectivity.edgesBySource.beginPass2();
    markerGraphConnectivity.edgesByTarget.beginPass2();
    setupLoadBalancing(markerGraphConnectivity.edges.size(), 100000);
    runThreads(&Assembler::createMarkerGraphConnectivityThreadFunction2, threadCount);
    markerGraphConnectivity.edgesBySource.endPass2();
    markerGraphConnectivity.edgesByTarget.endPass2();

    cout << timestamp << "createMarkerGraphConnectivity ends." << endl;
}



void Assembler::createMarkerGraphConnectivityThreadFunction0(size_t threadId)
{
    ostream& out = getLog(threadId);

    // Create the vector to contain the edges found by this thread.
    using std::shared_ptr;
    using std::make_shared;
    shared_ptr< MemoryMapped::Vector<MarkerGraphConnectivity::Edge> > thisThreadEdgesPointer =
        make_shared< MemoryMapped::Vector<MarkerGraphConnectivity::Edge> >();
    markerGraphConnectivity.threadEdges[threadId] = thisThreadEdgesPointer;
    MemoryMapped::Vector<MarkerGraphConnectivity::Edge>& thisThreadEdges = *thisThreadEdgesPointer;
    thisThreadEdges.createNew(
            largeDataName("tmp-ThreadGlobalMarkerGraphEdges-" + to_string(threadId)),
            largeDataPageSize);

    // Some things used inside the loop but defined here for performance.
    vector<GlobalMarkerGraphVertexId> children;
    MarkerGraphConnectivity::Edge edge;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << begin << endl;

        // Loop over all marker graph vertices assigned to this batch.
        for(GlobalMarkerGraphVertexId vertex0=begin; vertex0!=end; ++vertex0) {
            // out << timestamp << vertex0 << " " << globalMarkerGraphVertices.size(vertex0) << endl;
            edge.source = vertex0;
            const auto markerIds0 = globalMarkerGraphVertices[vertex0];

            // Skip it is it has too many markers.
            if(markerIds0.size() > markerGraphConnectivity.markerCountOverflow) {
                continue;
            }

            // Loop over children of this vertex.
            getGlobalMarkerGraphVertexChildren(vertex0, children);
            for(const GlobalMarkerGraphVertexId vertex1: children) {
                edge.target = vertex1;
                const auto markerIds1 = globalMarkerGraphVertices[vertex1];

                // Skip it is it has too many markers.
                if(markerIds1.size() > markerGraphConnectivity.markerCountOverflow) {
                    continue;
                }

                // Joint loop over markers to compute coverage.
                // This code is similar to LocalMarkerGraph::storeEdgeInfo.
                uint32_t coverage = 0;

                // Find pairs of markers for the same oriented read in the two vertices.
                // We exploit the fact that the markers in each
                // of the vertices are sorted.
                auto it0 = markerIds0.begin();
                auto it1 = markerIds1.begin();
                const auto end0 = markerIds0.end();
                const auto end1 = markerIds1.end();
                while(it0!=end0 && it1!=end1) {
                    const MarkerId markerId0 = *it0;
                    const MarkerId markerId1 = *it1;

                    // Find the oriented read ids and ordinals for these markers.
                    // This uses findMarkerId which requires binary searches
                    // and therefore could be expensive.
                    // If this becomes a performance problem we can
                    // store the OrientedReadId in each marker (4 extra bytes per marker).
                    OrientedReadId orientedReadId0;
                    uint32_t ordinal0;
                    tie(orientedReadId0, ordinal0) = findMarkerId(markerId0);
                    OrientedReadId orientedReadId1;
                    uint32_t ordinal1;
                    tie(orientedReadId1, ordinal1) = findMarkerId(markerId1);

                    if(orientedReadId0 < orientedReadId1) {
                        ++it0;
                        continue;
                    }
                    if(orientedReadId1 < orientedReadId0) {
                        ++it1;
                        continue;
                    }

                    // If getting here, the two oriented read ids are the same.
                    CZI_ASSERT(orientedReadId0 == orientedReadId1);
                    const OrientedReadId orientedReadId = orientedReadId0;

                    // Find the range of marker ids that correspond to this orientedReadId.
                    const auto thisOrientedReadMarkers = markers[orientedReadId.getValue()];
                    const MarkerId markerIdEnd   = thisOrientedReadMarkers.end()   - markers.begin();


                    // Find the streaks of markers for the same oriented readId.
                    auto it0StreakEnd = it0;
                    while(it0StreakEnd!=end0 && *it0StreakEnd<markerIdEnd) {
                        ++it0StreakEnd;
                    }
                    auto it1StreakEnd = it1;
                    while(it1StreakEnd!=end1 && *it1StreakEnd<markerIdEnd) {
                        ++it1StreakEnd;
                    }


                    // Only do it if both streaks contain one marker,
                    // the ordinal for the source vertex
                    // is less than the ordinal for the target vertex,
                    // and there are no intervening markers that also belong to a
                    // vertex of the marker graph.
                    if(it0StreakEnd-it0==1 && it1StreakEnd-it1==1 && ordinal0<ordinal1) {

                        // Check that there are no intervening markers that also belong to a
                        // vertex of the marker graph.
                        bool interveningVertexFound = false;
                        for(MarkerId markerId=markerId0+1; markerId!=markerId1; markerId++) {
                            if(globalMarkerGraphVertex[markerId] != invalidCompressedGlobalMarkerGraphVertexId) {
                                interveningVertexFound = true;
                                break;
                            }

                        }
                        if(!interveningVertexFound) {
                            ++coverage;
                        }
                    }

                    // Update the iterators to point to the end of the streaks.
                    it0 = it0StreakEnd;
                    it1 = it1StreakEnd;
                }


                // Store this edge.
                if(coverage >= 255) {
                    edge.coverage = 255;
                } else {
                    edge.coverage = uint8_t(coverage);
                }
                thisThreadEdges.push_back(edge);

            }

        }
    }

}



void Assembler::createMarkerGraphConnectivityThreadFunction1(size_t threadId)
{
    createMarkerGraphConnectivityThreadFunction12(threadId, 1);
}
void Assembler::createMarkerGraphConnectivityThreadFunction2(size_t threadId)
{
    createMarkerGraphConnectivityThreadFunction12(threadId, 2);
}
void Assembler::createMarkerGraphConnectivityThreadFunction12(size_t threadId, size_t pass)
{
    CZI_ASSERT(pass==1 || pass==2);

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all marker graph edges assigned to this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            const auto& edge = markerGraphConnectivity.edges[i];
            if(pass == 1) {
                markerGraphConnectivity.edgesBySource.incrementCountMultithreaded(edge.source);
                markerGraphConnectivity.edgesByTarget.incrementCountMultithreaded(edge.target);
            } else {
                markerGraphConnectivity.edgesBySource.storeMultithreaded(edge.source, Uint40(i));
                markerGraphConnectivity.edgesByTarget.storeMultithreaded(edge.target, Uint40(i));
            }
        }
    }

}


void Assembler::accessMarkerGraphConnectivity(bool accessEdgesReadWrite)
{
    if(accessEdgesReadWrite) {
        markerGraphConnectivity.edges.accessExistingReadWrite(
            largeDataName("GlobalMarkerGraphEdges"));
    } else {
        markerGraphConnectivity.edges.accessExistingReadOnly(
            largeDataName("GlobalMarkerGraphEdges"));
    }
    markerGraphConnectivity.edgesBySource.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphEdgesBySource"));
    markerGraphConnectivity.edgesByTarget.accessExistingReadOnly(
        largeDataName("GlobalMarkerGraphEdgesByTarget"));
}



void Assembler::checkMarkerGraphConnectivityIsOpen()
{
    CZI_ASSERT(markerGraphConnectivity.edges.isOpen);
    CZI_ASSERT(markerGraphConnectivity.edgesBySource.isOpen());
    CZI_ASSERT(markerGraphConnectivity.edgesByTarget.isOpen());
}



// Locate the edge given the vertices.
const Assembler::MarkerGraphConnectivity::Edge*
    Assembler::MarkerGraphConnectivity::findEdge(Uint40 source, Uint40 target) const
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



// Flag as not good a marker graph edge if:
// - It has coverage<minCoverage, AND
// - A path of length <= maxPathLength edges exists that:
//    * Starts at the source vertex of the edge.
//    * Ends at the target vertex of the edge.
//    * Only uses edges with coverage>=minCoverage.
void Assembler::flagMarkerGraphEdges(
    size_t threadCount,
    size_t minCoverage,
    size_t maxPathLength)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Store the parameters so all threads can see them.
    flagMarkerGraphEdgesData.minCoverage = minCoverage;
    flagMarkerGraphEdgesData.maxPathLength = maxPathLength;

    // Do it in parallel.
    const size_t vertexCount = markerGraphConnectivity.edgesBySource.size();
    setupLoadBalancing(vertexCount, 100000);
    runThreads(
        &Assembler::flagMarkerGraphEdgesThreadFunction,
        threadCount,
        "threadLogs/flagMarkerGraphEdges");
}



void Assembler::flagMarkerGraphEdgesThreadFunction(size_t threadId)
{
    ostream& out = getLog(threadId);

    const size_t minCoverage = flagMarkerGraphEdgesData.minCoverage;
    const size_t maxPathLength = flagMarkerGraphEdgesData.maxPathLength;

    vector< vector<GlobalMarkerGraphVertexId> > verticesByDistance(maxPathLength+1);
    verticesByDistance.front().resize(1);
    vector<GlobalMarkerGraphVertexId> vertices;

    // Loop over all batches assigned to this thread.
    size_t begin, end;
    while(getNextBatch(begin, end)) {
        out << timestamp << begin << endl;

        // Loop over vertices assigned to this batch.
        for(GlobalMarkerGraphVertexId startVertexId=begin; startVertexId!=end; ++startVertexId) {

            // Find all vertices within maxPathLength of this vertex,
            // computed using only edges with coverage >= minCoverage.
            // This uses a very simple algorithm that should work well
            // because the number of vertices involved is usually very small.
            verticesByDistance.front().front() = startVertexId;
            vertices.clear();
            vertices.push_back(startVertexId);
            for(size_t distance=1; distance<=maxPathLength; distance++) {
                auto& verticesAtThisDistance = verticesByDistance[distance];
                verticesAtThisDistance.clear();
                for(const GlobalMarkerGraphVertexId vertexId0: verticesByDistance[distance-1]) {
                    for(const uint64_t i: markerGraphConnectivity.edgesBySource[vertexId0]) {
                        const auto& edge = markerGraphConnectivity.edges[i];
                        if(edge.coverage < minCoverage) {
                            continue;
                        }
                        CZI_ASSERT(edge.source == vertexId0);
                        const GlobalMarkerGraphVertexId vertexId1 = edge.target;
                        verticesAtThisDistance.push_back(vertexId1);
                        vertices.push_back(vertexId1);
                    }
                }
            }
            sort(vertices.begin(), vertices.end());
            vertices.resize(unique(vertices.begin(), vertices.end()) - vertices.begin());

            // Loop over low coverage edges with source startVertexId.
            // If the target is one of the vertices we found, flag the edge as not good.
            for(const uint64_t i: markerGraphConnectivity.edgesBySource[startVertexId]) {
                auto& edge = markerGraphConnectivity.edges[i];
                if(edge.coverage >= minCoverage) {
                    continue;
                }
                CZI_ASSERT(edge.source == startVertexId);
                if(std::binary_search(vertices.begin(), vertices.end(), edge.target)) {
                    edge.flag0 = 1;
                }
            }

        }

    }

}

