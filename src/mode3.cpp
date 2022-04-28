
// Shasta
#include "mode3.hpp"
#include "findMarkerId.hpp"
#include "MarkerGraph.hpp"
#include "orderPairs.hpp"
#include "ReadFlags.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/right_open_interval.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include <map>
#include <queue>
#include <set>



// Each  linear chain of marker graph edges generates a segment.
void AssemblyGraph::createSegments()
{
    const bool debug = false;

    const MarkerGraph::EdgeId edgeCount = markerGraph.edges.size();
    vector<bool> wasFound(edgeCount, false);

    using MarkerGraphPath = vector<MarkerGraph::EdgeId>;
    MarkerGraphPath nextEdges;
    MarkerGraphPath previousEdges;
    MarkerGraphPath path;
    MarkerGraphPath reverseComplementedPath;

    // Main loop over all edges of the marker graph.
    // At each iteration we find a new linear path of edges.
    for(MarkerGraph::EdgeId startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {

        // If we already found this edge, skip it.
        // It is part of a path we already found.
        if(wasFound[startEdgeId]) {
            continue;
        }

        if(debug) {
            cout << "Starting a new path at edge " << startEdgeId << endl;
        }

        // Follow the path forward.
        nextEdges.clear();
        MarkerGraph::EdgeId edgeId = startEdgeId;
        bool isCircular = false;
        while(true) {
            const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
            const MarkerGraph::VertexId v1 = edge.target;
            const auto outEdges = markerGraph.edgesBySource[v1];
            if(outEdges.size() != 1) {
                break;
            }
            const auto inEdges = markerGraph.edgesByTarget[v1];
            if(inEdges.size() != 1) {
                break;
            }
            edgeId = outEdges[0];
            if(edgeId == startEdgeId) {
                isCircular = true;
                break;
            }
            nextEdges.push_back(edgeId);
            SHASTA_ASSERT(not wasFound[edgeId]);
            if(debug) {
                cout << "Moving forward: added " << edgeId << endl;
            }
        }

        // Follow the path backward.
        previousEdges.clear();
        if(!isCircular) {
            edgeId = startEdgeId;
            while(true) {
                const MarkerGraph::Edge edge = markerGraph.edges[edgeId];
                const MarkerGraph::VertexId v0 = edge.source;
                const auto outEdges = markerGraph.edgesBySource[v0];
                if(outEdges.size() != 1) {
                    break;
                }
                const auto inEdges = markerGraph.edgesByTarget[v0];
                if(inEdges.size() != 1) {
                    break;
                }
                edgeId = inEdges[0];
                previousEdges.push_back(edgeId);
                SHASTA_ASSERT(not wasFound[edgeId]);
                if(debug) {
                    cout << "Moving backward: added " << edgeId << endl;
                }
            }
        }

        // Gather the path.
        path.clear();
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter(path));
        path.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter(path));

        // Mark all the edges in the path as found.
        for(const MarkerGraph::EdgeId edgeId: path) {
            if(wasFound[edgeId]) {
                cout << "Assertion failed at " << edgeId << endl;
                SHASTA_ASSERT(0);
            }
            wasFound[edgeId] = true;
        }

        // Store this path as a new segment.
        paths.appendVector();
        for(const MarkerGraphEdgeId edgeId: path) {
            paths.append(edgeId);
        }
    }



    // Check that all edges of the marker graph were found.
    SHASTA_ASSERT(find(wasFound.begin(), wasFound.end(), false) == wasFound.end());


    // Debug output: write the paths.
    if(debug) {
        ofstream csv("Paths.csv");
        for(uint64_t segmentId=0; segmentId<paths.size(); segmentId++) {
            const auto path = paths[segmentId];
            for(const MarkerGraphEdgeId edgeId: path) {
                csv << segmentId << ",";
                csv << edgeId << "\n";
            }
        }
    }

}



// Compute coverage for all segments.
// It is computed as aompute average marker graph edge coverage.
void AssemblyGraph::computeSegmentCoverage()
{
    // Initialize segmentCoverage.
    segmentCoverage.createNew(
        largeDataName("Mode3-SegmentCoverage"), largeDataPageSize);
    const uint64_t segmentCount = paths.size();
    segmentCoverage.resize(segmentCount);

    // Loop over all segments.
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {

        // Access the marker graph path for this segment.
        const span<MarkerGraphEdgeId> path = paths[segmentId];


        // Loop over this path.
        uint64_t coverageSum = 0.;
        for(uint64_t position=0; position<path.size(); position++) {
            MarkerGraphEdgeId& edgeId = path[position];

            // Add the marker intervals on this marker graph edge.
            const span<const MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
            coverageSum += markerIntervals.size();
        }

        segmentCoverage[segmentId] = float(coverageSum) / float(path.size());

    }
}



// For each marker graph edge, store in the marker graph edge table
// the corresponding (segment)
// and position in the path, if any.
// This is needed when computing pseudopaths.
void AssemblyGraph::computeMarkerGraphEdgeTable(size_t threadCount)
{

    // Initialize the marker graph edge table.
    markerGraphEdgeTable.createNew(largeDataName("Mode3-MarkerGraphEdgeTable"), largeDataPageSize);
    markerGraphEdgeTable.resize(markerGraph.edges.size());
    fill(markerGraphEdgeTable.begin(), markerGraphEdgeTable.end(), make_pair(
        std::numeric_limits<uint64_t>::max(),
        std::numeric_limits<uint32_t>::max()
        ));

    // Fill in the marker graph edge table.
    const uint64_t batchSize = 100;
    setupLoadBalancing(paths.size(), batchSize);
    runThreads(&AssemblyGraph::computeMarkerGraphEdgeTableThreadFunction, threadCount);
}



void AssemblyGraph::computeMarkerGraphEdgeTableThreadFunction(size_t threadId)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all vertices assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; ++segmentId) {
            const span<MarkerGraphEdgeId> path = paths[segmentId];

            // Loop over the path of this segment.
            for(uint64_t position=0; position<path.size(); position++) {
                const MarkerGraphEdgeId edgeId = path[position];

                // Store the marker graph edge table entry for this edge.
                SHASTA_ASSERT(edgeId < markerGraphEdgeTable.size());
                markerGraphEdgeTable[edgeId] = make_pair(segmentId, position);
            }
        }

    }
}



void AssemblyGraph::computePseudoPaths(size_t threadCount)
{
    pseudoPaths.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "tmp-mode3-PseudoPaths-1"),
        largeDataPageSize);

    uint64_t batchSize = 1000;
    pseudoPaths.beginPass1(markers.size());
    setupLoadBalancing(markerGraphEdgeTable.size(), batchSize);
    runThreads(&AssemblyGraph::computePseudoPathsPass1, threadCount);
    pseudoPaths.beginPass2();
    setupLoadBalancing(markerGraphEdgeTable.size(), batchSize);
    runThreads(&AssemblyGraph::computePseudoPathsPass2, threadCount);
    pseudoPaths.endPass2();

    batchSize = 100;
    setupLoadBalancing(pseudoPaths.size(), batchSize);
    runThreads(&AssemblyGraph::sortPseudoPaths, threadCount);
}



void AssemblyGraph::computePseudoPathsPass1(size_t threadId)
{
    computePseudoPathsPass12(1);
}



void AssemblyGraph::computePseudoPathsPass2(size_t threadId)
{
    computePseudoPathsPass12(2);
}



void AssemblyGraph::computePseudoPathsPass12(uint64_t pass)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker graph edges assigned to this batch.
        for(MarkerGraph::EdgeId edgeId=begin; edgeId!=end; ++edgeId) {
            const auto& p = markerGraphEdgeTable[edgeId];
            const uint64_t segmentId = p.first;
            const uint32_t position = p.second;
            SHASTA_ASSERT(segmentId != std::numeric_limits<uint64_t>::max());
            SHASTA_ASSERT(position != std::numeric_limits<uint32_t>::max());

            // Loop over the marker intervals of this marker graph edge..
            const auto markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
            for(const MarkerInterval& markerInterval: markerIntervals) {
                const OrientedReadId orientedReadId = markerInterval.orientedReadId;

                if(pass == 1) {
                    pseudoPaths.incrementCountMultithreaded(orientedReadId.getValue());
                } else {
                    PseudoPathEntry pseudoPathEntry;
                    pseudoPathEntry.segmentId = segmentId;
                    pseudoPathEntry.position = position;
                    pseudoPathEntry.ordinals = markerInterval.ordinals;
                    pseudoPaths.storeMultithreaded(orientedReadId.getValue(), pseudoPathEntry);
                }
            }
        }
    }
}



void AssemblyGraph::sortPseudoPaths(size_t threadId)
{
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker graph edges assigned to this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            auto pseudoPath = pseudoPaths[i];
            sort(pseudoPath.begin(), pseudoPath.end());
        }
    }
}


// The compressed pseudopath of an oriented read
// is the sequence of segmentIds it encounters.
// For each entry, we also store the average offset
// of the beginning of the oriented read relative
// to the beginning of the segment, in markers.
void AssemblyGraph::computeCompressedPseudoPaths()
{
    // Initialize the compressed pseudopaths.
    compressedPseudoPaths.createNew(
        largeDataName("Mode3-CompressedPseudoPaths"), largeDataPageSize);

    // Work vector defined outside the loop to reduce memory allocation overhead.
    vector<uint64_t> compressedPseudoPath;

    // Loop over all oriented reads.
    for(uint64_t i=0; i<pseudoPaths.size(); i++) {

        // Access the pseudopath for this oriented read.
        const span<PseudoPathEntry> pseudoPath = pseudoPaths[i];

        // Compute the compressed pseudopath.
        computeCompressedPseudoPath(pseudoPath, compressedPseudoPath);

        // Store it.
        compressedPseudoPaths.appendVector(compressedPseudoPath);
    }



    // Write them out.
    ofstream csv("CompressedPseudoPaths.csv");
    for(uint64_t i=0; i<pseudoPaths.size(); i++) {
        const ReadId readId = ReadId(i >> 1);
        const Strand strand = i & 1;
        const OrientedReadId orientedReadId(readId, strand);
        const span<uint64_t> compressedPseudoPath = compressedPseudoPaths[i];

        csv << orientedReadId << ",";
        for(const uint64_t segmentId: compressedPseudoPath) {
            csv << segmentId << ",";
        }
        csv << endl;
    }
}



void AssemblyGraph::computeCompressedPseudoPath(
    const span<PseudoPathEntry> pseudoPath,
    vector<uint64_t>& compressedPseudoPath)
{
    // Start with an empty compressed pseudopath.
    compressedPseudoPath.clear();

    // Loop over the pseudopath.
    for(uint64_t i=0; i<pseudoPath.size(); /* Increment later */) {
        const PseudoPathEntry& pseudoPathEntry = pseudoPath[i];
        const uint64_t segmentId = pseudoPathEntry.segmentId;

        // Move to the end of the streak of pseudoPath entries with the same segmentId.
        for(i++; i<pseudoPath.size() and (pseudoPath[i].segmentId == segmentId); i++) {
        }

        // Store this segmentId in the compressed pseudopath.
        compressedPseudoPath.push_back(segmentId);

    }
}



void AssemblyGraph::findTransitions(std::map<SegmentPair, Transitions>& transitionMap)
{
    transitionMap.clear();

    for(ReadId readId=0; readId<pseudoPaths.size()/2; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto pseudoPath = pseudoPaths[orientedReadId.getValue()];

            if(pseudoPath.size() < 2) {
                continue;
            }

            for(uint64_t i=1; i<pseudoPath.size(); i++) {
                const auto& previous = pseudoPath[i-1];
                const auto& current = pseudoPath[i];
                if(previous.segmentId == current.segmentId) {
                    continue;
                }

                const SegmentPair segmentPair = make_pair(previous.segmentId, current.segmentId);
                transitionMap[segmentPair].push_back(
                    make_pair(orientedReadId, Transition({previous, current})));

            }
        }
    }
}



void AssemblyGraph::createLinks(
    const std::map<SegmentPair, Transitions>& transitionMap,
    uint64_t minCoverage)
{
    links.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "Mode3-Links"),
        largeDataPageSize);
    transitions.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "Mode3-Transitions"),
        largeDataPageSize);
    for(const auto& p: transitionMap) {
        const auto& transitionVector = p.second;
        const uint64_t coverage = transitionVector.size();
        if(coverage >= minCoverage) {
            const uint64_t segmentId0 = p.first.first;
            const uint64_t segmentId1 = p.first.second;
            links.push_back(Link(segmentId0, segmentId1, coverage));
            transitions.appendVector(transitionVector);
        }
    }
}



// Initial construction of the AssemblyGraph.
AssemblyGraph::AssemblyGraph(
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize,
    size_t threadCount,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph) :
    MultithreadedObject<AssemblyGraph>(*this),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize),
    markers(markers),
    markerGraph(markerGraph)
{
    // Minimum number of transitions (oriented reads) to create a link.
    const uint64_t minCoverage = 2; // EXPOSE WHEN CODE STABILIZES

    // Create a segment for each linear chain of marker graph edges.
    paths.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "Mode3-Paths"),
        largeDataPageSize);
    createSegments();
    computeSegmentCoverage();

    // Keep track of the segment and position each marker graph edge corresponds to.
    computeMarkerGraphEdgeTable(threadCount);

    // Compute pseudopaths of all oriented reads.
    computePseudoPaths(threadCount);
    computeCompressedPseudoPaths();

    // Find pseudopath transitions and store them keyed by the pair of segments.
    std::map<SegmentPair, Transitions> transitionMap;
    findTransitions(transitionMap);

    // Create a links between pairs of segments with a sufficient number of transitions.
    createLinks(transitionMap, minCoverage);
    pseudoPaths.remove();
    createConnectivity();

    cout << "The mode 3 assembly graph has " << paths.size() << " segments and " <<
        links.size() << " links." << endl;
}



string AssemblyGraph::largeDataName(const string& name) const
{
    if(largeDataFileNamePrefix.empty()) {
        return "";  // Anonymous;
    } else {
        return largeDataFileNamePrefix + name;
    }
}



// Constructor from binary data.
AssemblyGraph::AssemblyGraph(
    const string& largeDataFileNamePrefix,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph) :
    MultithreadedObject<AssemblyGraph>(*this),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    markers(markers),
    markerGraph(markerGraph)
{
    paths.accessExistingReadOnly(largeDataName("Mode3-Paths"));
    segmentCoverage.accessExistingReadOnly(
        largeDataName("Mode3-SegmentCoverage"));
    markerGraphEdgeTable.accessExistingReadOnly(largeDataName("Mode3-MarkerGraphEdgeTable"));
    compressedPseudoPaths.accessExistingReadOnly(
        largeDataName("Mode3-CompressedPseudoPaths"));
    links.accessExistingReadOnly(largeDataName("Mode3-Links"));
    transitions.accessExistingReadOnly(largeDataName("Mode3-Transitions"));
    linksBySource.accessExistingReadOnly(largeDataName("Mode3-LinksBySource"));
    linksByTarget.accessExistingReadOnly(largeDataName("Mode3-LinksByTarget"));
    clusterIds.accessExistingReadOnly(largeDataName("Mode3-ClusterIds"));
}



void AssemblyGraph::createConnectivity()
{
    linksBySource.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "Mode3-LinksBySource"),
        largeDataPageSize);
    linksByTarget.createNew(
        largeDataFileNamePrefix.empty() ? "" : (largeDataFileNamePrefix + "Mode3-LinksByTarget"),
        largeDataPageSize);

    linksBySource.beginPass1(links.size());
    linksByTarget.beginPass1(links.size());
    for(uint64_t linkId=0; linkId<links.size(); linkId++) {
        const Link& link = links[linkId];
        linksBySource.incrementCount(link.segmentId0);
        linksByTarget.incrementCount(link.segmentId1);
    }
    linksBySource.beginPass2();
    linksByTarget.beginPass2();
    for(uint64_t linkId=0; linkId<links.size(); linkId++) {
        const Link& link = links[linkId];
        linksBySource.store(link.segmentId0, linkId);
        linksByTarget.store(link.segmentId1, linkId);
    }
    linksBySource.endPass2();
    linksByTarget.endPass2();
}



void AssemblyGraph::writeGfa(const string& fileName) const
{
    ofstream gfa(fileName);
    writeGfa(gfa);
}



void AssemblyGraph::writeGfa(ostream& gfa) const
{
    // Write the header.
    gfa << "H\tVN:Z:1.0\n";

    // Write the segments.
    for(uint64_t segmentId=0; segmentId<paths.size(); segmentId++) {
        const auto path = paths[segmentId];
        gfa <<
            "S\t" << segmentId << "\t" <<
            "*\tLN:i:" << path.size() << "\n";
    }

    // Write the liks.
    for(const Link& link: links) {
        gfa << "L\t" <<
            link.segmentId0 << "\t+\t" <<
            link.segmentId1 << "\t+\t0M\n";
    }

}



// Find the distinct oriented reads that appear on the path
// of a segment. Also return the average edge coverage for the path.
double AssemblyGraph::findOrientedReadsOnSegment(
    uint64_t segmentId,
    vector<OrientedReadId>& orientedReadIdsArgument) const
{
    // Loop over the marker graph path corresponding to this segment.
    const span<const MarkerGraphEdgeId> path = paths[segmentId];
    double coverage = 0.;
    std::set<OrientedReadId> orientedReadIds;
    for(const MarkerGraphEdgeId& edgeId: path) {

        // Loop over the marker intervals for this marker graph edge.
        const span<const MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
        coverage += double(markerIntervals.size());
        for(const MarkerInterval& markerInterval: markerIntervals) {
            orientedReadIds.insert(markerInterval.orientedReadId);
        }
    }

    // Copy the oriented reads to the vector passed as an argument.
    orientedReadIdsArgument.clear();
    orientedReadIdsArgument.insert(orientedReadIdsArgument.end(),
        orientedReadIds.begin(), orientedReadIds.end());

    return coverage / double(path.size());
}



// Get information about the oriented reads that appear on the
// marker graph path of a segment.
void AssemblyGraph::getOrientedReadsOnSegment(
    uint64_t segmentId,
    SegmentOrientedReadInformation& information) const
{
    // A data structure that, for each oriented read we find,
    // contains a sum of offsets and the number of marker graph vertices
    // that contributed to the sum.
    std::map<OrientedReadId, pair<uint64_t, int64_t>  > table;

    // Loop over the marker graph path corresponding to this segment.
    const span<const MarkerGraphEdgeId> path = paths[segmentId];
    std::set<OrientedReadId> orientedReadIds;
    for(uint64_t position=0; position<path.size(); position++) {
        const MarkerGraphEdgeId& edgeId = path[position];

        // Loop over the marker intervals for this marker graph edge.
        const span<const MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
        for(const MarkerInterval& markerInterval: markerIntervals) {
            const OrientedReadId orientedReadId = markerInterval.orientedReadId;

            // Update our table for this oriented read.
            auto it = table.find(orientedReadId);
            if(it == table.end()) {
                tie(it, ignore) = table.insert(make_pair(orientedReadId, make_pair(0ULL, 0LL)));
            }
            auto& p = it->second;
            p.first += 2;
            p.second += int32_t(position) - int32_t(markerInterval.ordinals[0]);
            p.second += int32_t(position + 1) -int32_t(markerInterval.ordinals[1]);
        }
    }



    // Store what we found.
    information.infos.clear();
    for(const auto& p: table) {
        SegmentOrientedReadInformation::Info info;
        info.orientedReadId = p.first;
        const uint64_t n = p.second.first;
        const int64_t sum = p.second.second;
        info.averageOffset = int32_t(std::round(double(sum) / double(n)));
        information.infos.push_back(info);
    }
 }



// Estimate the offset between two segments.
// Takes as input SegmentOrientedReadInformation objects
// for the two segments.
// Common oriented reads between the two segments are used
// to estimate the average offset, in markers,
// between the beginning of the segments.
// The number of common oriented reads
// is computed and stored in the last argument.
// If that is zero, the computed offset is not valid.
void AssemblyGraph::estimateOffset(
    const SegmentOrientedReadInformation& info0,
    const SegmentOrientedReadInformation& info1,
    int64_t& offset,
    uint64_t& commonOrientedReadCount
    ) const
{
    offset = 0;
    commonOrientedReadCount = 0;

    // Joint loop over common oriented reads in the two segments.
    const auto begin0 = info0.infos.begin();
    const auto begin1 = info1.infos.begin();
    const auto end0 = info0.infos.end();
    const auto end1 = info1.infos.end();
    auto it0 = begin0;
    auto it1 = begin1;
    while((it0 != end0) and (it1 != end1)) {

        if(it0->orientedReadId < it1->orientedReadId) {
            ++it0;
        } else if(it1->orientedReadId < it0->orientedReadId) {
            ++it1;
        } else {
            SHASTA_ASSERT(it0->orientedReadId == it1->orientedReadId);

            commonOrientedReadCount++;
            offset += (int64_t(it0->averageOffset) - int64_t(it1->averageOffset));

            ++it0;
            ++it1;
        }
    }

    if(commonOrientedReadCount) {
        offset = int64_t(std::round(double(offset) / double(commonOrientedReadCount)));
    } else {
        offset = std::numeric_limits<uint64_t>::max();
    }

}



// Analyze a pair of segments for common oriented reads,
// offsets, missing reads, etc.
void AssemblyGraph::analyzeSegmentPair(
    uint64_t segmentId0,
    uint64_t segmentId1,
    const SegmentOrientedReadInformation& info0,
    const SegmentOrientedReadInformation& info1,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    SegmentPairInformation& info01
    ) const
{
    using boost::icl::discrete_interval;
    using boost::icl::intersects;

    // Store the number of oriented reads in each segment.
    info01.totalCount[0] = info0.infos.size();
    info01.totalCount[1] = info1.infos.size();

    // Use common oriented reads to estimate the offset between the two segments.
    // If there are no common oriented reads, stop here.
    estimateOffset(info0, info1, info01.offset, info01.commonCount);
    if(info01.commonCount == 0) {
        return;
    }


    // Count the oriented reads missing from each segment,
    // and which should have been present based on
    // the known relative offsets.
    info01.unexplainedCount = {0, 0};
    info01.shortCount = {0, 0};

    // Set up a joint loop over oriented reads in the two segments.
    const auto begin0 = info0.infos.begin();
    const auto begin1 = info1.infos.begin();
    const auto end0 = info0.infos.end();
    const auto end1 = info1.infos.end();
    auto it0 = begin0;
    auto it1 = begin1;

    const uint64_t length0 = paths.size(segmentId0);
    const uint64_t length1 = paths.size(segmentId1);
    while(true) {

        // At end of both segments.
        if((it0 == end0) and (it1 == end1)) {
            break;
        }



        // This read only appears on segment 0.
        if((it1 == end1) or ((it0 != end0) and (it0->orientedReadId < it1->orientedReadId))) {
            const int64_t orientedReadLength = markers.size(it0->orientedReadId.getValue());

            // Compute the hypothetical range of the oriented read relative
            // to the beginning of segment 1.
            const discrete_interval<int64_t> orientedReadRange1(
                it0->averageOffset - info01.offset,
                it0->averageOffset - info01.offset + orientedReadLength);
            const discrete_interval<int64_t> segment1Range(0, length1);

            // Figure out if it the oriented read would overlap segment 1.
            const bool wouldOverlap = intersects(orientedReadRange1, segment1Range);

            if(wouldOverlap) {
                ++info01.unexplainedCount[0];
            } else {
                ++info01.shortCount[0];
            }

            SHASTA_ASSERT(it0 != end0);
            ++it0;
        }



        // Only on segment 1
        else if((it0 == end0) or ((it1 != end1) and (it1->orientedReadId < it0->orientedReadId))) {
            const int64_t orientedReadLength = markers.size(it1->orientedReadId.getValue());

            // Compute the hypothetical range of the oriented read relative
            // to the beginning of segment 0.
            const discrete_interval<int64_t> orientedReadRange0(
                it1->averageOffset + info01.offset,
                it1->averageOffset + info01.offset + orientedReadLength);
            const discrete_interval<int64_t> segment0Range(0, length0);

            // Figure out if it the oriented read would overlap segment 0.
            const bool wouldOverlap = intersects(orientedReadRange0, segment0Range);

            if(wouldOverlap) {
                ++info01.unexplainedCount[1];
            } else {
                ++info01.shortCount[1];
            }

            SHASTA_ASSERT(it1 != end1);
            ++it1;
        }

        // On both segments.
        else {
            SHASTA_ASSERT(it0 != end0);
            SHASTA_ASSERT(it1 != end1);
            ++it0;
            ++it1;
        }
    }

    info01.check();

}



// Gather oriented read information for each segment.
void AssemblyGraph::storeSegmentOrientedReadInformation(size_t threadCount)
{
    const uint64_t segmentCount = paths.size();
    segmentOrientedReadInformation.resize(segmentCount);
    const uint64_t batchSize = 10;
    setupLoadBalancing(segmentCount, batchSize);
    runThreads(&AssemblyGraph::storeSegmentOrientedReadInformationThreadFunction, threadCount);
}




// Gather oriented read information for each segment.
void AssemblyGraph::storeSegmentOrientedReadInformationThreadFunction(size_t threadId)
{

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over segments assigned to this batch.
        for(uint64_t segmentId=begin; segmentId!=end; ++segmentId) {

            // Get oriented read information for this segment.
            getOrientedReadsOnSegment(segmentId, segmentOrientedReadInformation[segmentId]);
        }
    }
}

void AssemblyGraph::clusterSegments(size_t threadCount, uint64_t minClusterSize)
{
    // Gather oriented read information for all segments.
    storeSegmentOrientedReadInformation(threadCount);

    // Find the segment pairs.
    const uint64_t segmentCount = paths.size();
    const uint64_t batchSize = 10;
    setupLoadBalancing(segmentCount, batchSize);
    clusterSegmentsData.threadPairs.resize(threadCount);
    runThreads(&AssemblyGraph::clusterSegmentsThreadFunction1, threadCount);

    // For now, write a dot file with the pairs.
    ofstream dot("SegmentGraph.dot");
    dot << "graph segmentGraph {\n";
    for(const auto& threadPairs: clusterSegmentsData.threadPairs) {
        for(const auto& p: threadPairs) {
            dot << p.first << "--" << p.second << ";\n";
        }
    }
    dot << "}\n";



    // The segment pairs we found define a subgraph of the assembly graph.
    // Compute connected components of this subgraph.
    // The connected components of sufficient size become clusters.
    vector<uint64_t> rank(segmentCount);
    vector<uint64_t> parent(segmentCount);
    boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        disjointSets.make_set(segmentId);
    }
    for(const auto& threadPairs: clusterSegmentsData.threadPairs) {
        for(const auto& p: threadPairs) {
            disjointSets.union_set(p.first, p.second);
        }
    }

    // Gather the segments in each connected component.
    vector< vector<uint64_t> > components(segmentCount);
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        const uint64_t componentId = disjointSets.find_set(segmentId);
        components[componentId].push_back(segmentId);
    }

    // Each connected components of size at least minClusterSize
    // becomes a cluster.
    vector< pair<uint64_t, uint64_t> > clusterTable;
    for(uint64_t componentId=0; componentId<segmentCount; componentId++) {
        const vector<uint64_t>& component = components[componentId];
        const uint64_t componentSize = component.size();
        if(component.size() >= minClusterSize) {
            clusterTable.push_back(make_pair(componentId, componentSize));
        }
    }

    // Sort the clusters by decreasing size.
    sort(clusterTable.begin(), clusterTable.end(),
        OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    cout << "Found " << clusterTable.size() << " segment clusters with the following sizes:" << endl;
    uint64_t clusteredSegmentCount = 0;
    for(uint64_t clusterId=0; clusterId<clusterTable.size(); clusterId++) {
        const auto& p = clusterTable[clusterId];
        const uint64_t componentSize = p.second;
        cout << " " << componentSize;
        clusteredSegmentCount += componentSize;
    }
    cout << endl;
    cout << "Out of " << segmentCount << " segments, " <<
        clusteredSegmentCount << " were assigned to a cluster." << endl;



    // Store the cluster id of each segment.
    clusterIds.createNew(largeDataName("Mode3-ClusterIds"), largeDataPageSize);
    clusterIds.resize(segmentCount);
    fill(clusterIds.begin(), clusterIds.end(), std::numeric_limits<uint64_t>::max());
    for(uint64_t clusterId=0; clusterId<clusterTable.size(); clusterId++) {
        const auto& p = clusterTable[clusterId];
        const uint64_t componentId = p.first;
        const vector<uint64_t>& cluster = components[componentId];
        for(const uint64_t segmentId: cluster) {
            clusterIds[segmentId] = clusterId;
        }
    }



    // Clean up.
    clusterSegmentsData.threadPairs.clear();
    clusterSegmentsData.threadPairs.shrink_to_fit();
    segmentOrientedReadInformation.clear();
    segmentOrientedReadInformation.shrink_to_fit();
}



void AssemblyGraph::clusterSegmentsThreadFunction1(size_t threadId)
{
    // EXPOSE THESE CONSTANTS WHEN CODE STABILIZES.
    const uint64_t maxDistance = 20;
    const uint64_t minCommonReadCount = 4;
    const double maxUnexplainedFraction = 0.15;

    auto& threadPairs = clusterSegmentsData.threadPairs[threadId];
    threadPairs.clear();
    vector<uint64_t> descendants;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over segments assigned to this batch.
        for(uint64_t segmentId0=begin; segmentId0!=end; ++segmentId0) {
            findDescendants(segmentId0, maxDistance, descendants);

            // Check all the descendants.
            for(const uint64_t segmentId1: descendants) {
                SegmentPairInformation info;
                analyzeSegmentPair(segmentId0, segmentId1,
                    segmentOrientedReadInformation[segmentId0],
                    segmentOrientedReadInformation[segmentId1],
                    markers, info);
                if(info.commonCount < minCommonReadCount) {
                    continue;
                }
                if(info.maximumUnexplainedFraction() > maxUnexplainedFraction) {
                    continue;
                }
                pair<uint64_t, uint64_t> segmentPair(segmentId0, segmentId1);
                if(segmentPair.first > segmentPair.second) {
                    swap(segmentPair.first, segmentPair.second);
                }
                threadPairs.push_back(segmentPair);

            }
        }
    }
}



// Find descendants of a given segment, up to a given distance in the graph.
void AssemblyGraph::findDescendants(
    uint64_t startSegmentId,
    uint64_t maxDistance,
    vector<uint64_t>& descendants
    ) const
{
    // Initialize the BFS.
    descendants.clear();
    std::queue<uint64_t> q;
    q.push(startSegmentId);
    std::map<uint64_t, uint64_t> distanceMap;
    distanceMap.insert(make_pair(startSegmentId, 0));

    // BFS loop.
    while(not q.empty()) {
        const uint64_t segmentId0 = q.front();
        q.pop();

        const uint64_t distance0 = distanceMap[segmentId0];
        const uint64_t distance1 = distance0 + 1;

        // Loop over children of segmentId0.
        for(const uint64_t linkId01: linksBySource[segmentId0]) {
            const Link& link01 = links[linkId01];
            const uint64_t segmentId1 = link01.segmentId1;

            // If we already encountered segmentId1, skip it.
            if(distanceMap.find(segmentId1) != distanceMap.end()) {
                continue;
            }

            descendants.push_back(segmentId1);
            distanceMap.insert(make_pair(segmentId1, distance1));
            if(distance1 < maxDistance) {
                q.push(segmentId1);
            }
        }
    }
}




