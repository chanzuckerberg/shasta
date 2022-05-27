
// Shasta
#include "mode3.hpp"
#include "findMarkerId.hpp"
#include "MarkerGraph.hpp"
#include "orderPairs.hpp"
#include "shastaLapack.hpp"
#include "ReadFlags.hpp"
#include "SubsetGraph.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
// Include disjoint_sets.hpp first to avoid Boost problems.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/right_open_interval.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>

// Standard library.

#include <bitset>
#include <map>
#include <queue>
#include <set>
#include <unordered_set>



// Each  linear chain of marker graph edges generates a segment.
void AssemblyGraph::createSegmentPaths()
{
    const bool debug = false;

    createNew(paths, "Mode3-Paths");
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
// It is computed as average marker graph edge coverage
// over the marker graph edges in the path of each segment.
void AssemblyGraph::computeSegmentCoverage()
{
    // Initialize segmentCoverage.
    createNew(segmentCoverage, "Mode3-SegmentCoverage");
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


    // Write a histogram of segment coverage.
    vector<uint64_t> histogram;
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        const uint64_t coverage = uint64_t(std::round(segmentCoverage[segmentId]));
        if(coverage >= histogram.size()) {
            histogram.resize(coverage + 1, 0);
        }
        ++histogram[coverage];
    }
    ofstream csv("SegmentCoverageHistogram.csv");
    csv << "Coverage,Frequency\n";
    for(uint64_t coverage=0; coverage<histogram.size(); coverage++) {
        csv << coverage << "," << histogram[coverage] << "\n";
    }
}



// For each marker graph edge, store in the marker graph edge table
// the corresponding (segment)
// and position in the path, if any.
// This is needed when computing pseudopaths.
void AssemblyGraph::computeMarkerGraphEdgeTable(size_t threadCount)
{

    // Initialize the marker graph edge table.
    createNew(markerGraphEdgeTable, "Mode3-MarkerGraphEdgeTable");
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
    const bool debug = true;

    createNew(pseudoPaths, "tmp-mode3-PseudoPaths-1");

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

    if(debug) {
        ofstream csv("PseudoPaths.csv");
        csv << "OrientedReadId,SegmentId,Position,ordinal0,Ordinal1\n";
        for(uint64_t i=0; i<markers.size(); i++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
            const auto pseudoPath = pseudoPaths[i];
            for(uint64_t position=0; position<pseudoPath.size(); position++) {
                const auto& entry = pseudoPath[position];
                csv << orientedReadId << ",";
                csv << entry.segmentId << ",";
                csv << entry.position << ",";
                csv << entry.ordinals[0] << ",";
                csv << entry.ordinals[1] << "\n";
            }
        }

    }
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
    const bool debug = true;

    // Initialize the compressed pseudopaths.
    createNew(compressedPseudoPaths, "Mode3-CompressedPseudoPaths");

    // Work vector defined outside the loop to reduce memory allocation overhead.
    vector<CompressedPseudoPathEntry> compressedPseudoPath;

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
    if(debug) {
        ofstream csv("CompressedPseudoPaths.csv");
        for(uint64_t i=0; i<pseudoPaths.size(); i++) {
            const ReadId readId = ReadId(i >> 1);
            const Strand strand = i & 1;
            const OrientedReadId orientedReadId(readId, strand);
            const span<CompressedPseudoPathEntry> compressedPseudoPath = compressedPseudoPaths[i];

            csv << orientedReadId << ",";
            for(const CompressedPseudoPathEntry entry: compressedPseudoPath) {
                csv << entry.segmentId << ",";
            }
            csv << endl;
        }
    }



    // Write them out again, with more details.
    if(debug) {
        ofstream csv("CompressedPseudoPathsDetails.csv");
        csv << "OrientedReadId,Position,SegmentId,"
            "First position,First ordinal0,First ordinal1,"
            "Last position,Last ordinal0,Last ordinal1\n";
        for(uint64_t i=0; i<pseudoPaths.size(); i++) {
            const ReadId readId = ReadId(i >> 1);
            const Strand strand = i & 1;
            const OrientedReadId orientedReadId(readId, strand);
            const span<CompressedPseudoPathEntry> compressedPseudoPath = compressedPseudoPaths[i];

            for(uint64_t position=0; position<compressedPseudoPath.size(); position++) {
                const CompressedPseudoPathEntry& entry = compressedPseudoPath[position];
                const PseudoPathEntry& first = entry.pseudoPathEntries[0];
                const PseudoPathEntry& last = entry.pseudoPathEntries[1];
                csv << orientedReadId << ",";
                csv << position << ",";
                csv << entry.segmentId << ",";
                csv << first.position << ",";
                csv << first.ordinals[0] << ",";
                csv << first.ordinals[1] << ",";
                csv << last.position << ",";
                csv << last.ordinals[0] << ",";
                csv << last.ordinals[1] << "\n";
            }
        }
    }


}



void AssemblyGraph::computeCompressedPseudoPath(
    const span<PseudoPathEntry> pseudoPath,
    vector<CompressedPseudoPathEntry>& compressedPseudoPath)
{
    // Start with an empty compressed pseudopath.
    compressedPseudoPath.clear();

    // Loop over the pseudopath.
    for(uint32_t i=0; i<pseudoPath.size(); /* Increment later */) {
        const PseudoPathEntry& pseudoPathEntry = pseudoPath[i];
        const uint64_t segmentId = pseudoPathEntry.segmentId;

        // Move to the end of the streak of pseudoPath entries with the same segmentId.
        const uint32_t streakBegin = i;
        uint32_t streakEnd = streakBegin + 1;
        for(; streakEnd<pseudoPath.size() and (pseudoPath[streakEnd].segmentId == segmentId); streakEnd++) {
        }

        // Store this segmentId in the compressed pseudopath.
        CompressedPseudoPathEntry compressedPseudoPathEntry;
        compressedPseudoPathEntry.segmentId = segmentId;
        compressedPseudoPathEntry.pseudoPathEntries[0] =pseudoPath[streakBegin];
        compressedPseudoPathEntry.pseudoPathEntries[1] = pseudoPath[streakEnd - 1];
        compressedPseudoPath.push_back(compressedPseudoPathEntry);

        // Prepare to handle the next segment.
        i = streakEnd;
    }
}



void AssemblyGraph::computeSegmentCompressedPseudoPathInfo()
{
    const bool debug = true;

    const uint64_t segmentCount = paths.size();
    const uint64_t readCount = compressedPseudoPaths.size()/2;

    createNew(segmentCompressedPseudoPathInfo, "Mode3-SegmentCompressedPseudoPathInfo");

    // Pass 1.
    segmentCompressedPseudoPathInfo.beginPass1(segmentCount);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto compressedPseudoPath = compressedPseudoPaths[orientedReadId.getValue()];

            for(uint64_t position=0; position<compressedPseudoPath.size(); position++) {
                const CompressedPseudoPathEntry& compressedPseudoPathEntry =
                    compressedPseudoPath[position];
                const uint64_t segmentId = compressedPseudoPathEntry.segmentId;
                segmentCompressedPseudoPathInfo.incrementCount(segmentId);
            }
        }
    }

    // Pass 2.
    segmentCompressedPseudoPathInfo.beginPass2();
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto compressedPseudoPath = compressedPseudoPaths[orientedReadId.getValue()];

            for(uint64_t position=0; position<compressedPseudoPath.size(); position++) {
                const CompressedPseudoPathEntry& compressedPseudoPathEntry =
                    compressedPseudoPath[position];
                const uint64_t segmentId = compressedPseudoPathEntry.segmentId;
                segmentCompressedPseudoPathInfo.store(segmentId, make_pair(orientedReadId, position));
            }
        }
    }
    segmentCompressedPseudoPathInfo.endPass2();

    // Sort.
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        const auto v = segmentCompressedPseudoPathInfo[segmentId];
        sort(v.begin(), v.end());
    }


    if(debug) {
        ofstream csv("SegmentCompressedPseudoPathInfo.csv");
        csv << "SegmentId,OrientedReadId,Position in compressed pseudopath\n";
        for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
            const auto v = segmentCompressedPseudoPathInfo[segmentId];
            for(const auto& p: v) {
                csv << segmentId << ",";
                csv << p.first << ",";
                csv << p.second << "\n";
            }
        }
    }
}



void AssemblyGraph::findTransitions(std::map<SegmentPair, Transitions>& transitionMap)
{
    transitionMap.clear();

    for(ReadId readId=0; readId<compressedPseudoPaths.size()/2; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto compressedPseudoPath = compressedPseudoPaths[orientedReadId.getValue()];

            for(uint64_t i=1; i<compressedPseudoPath.size(); i++) {
                const auto& previous = compressedPseudoPath[i-1];
                const auto& current = compressedPseudoPath[i];
                SHASTA_ASSERT(previous.segmentId != current.segmentId);

                const SegmentPair segmentPair = make_pair(previous.segmentId, current.segmentId);
                transitionMap[segmentPair].push_back(
                    make_pair(orientedReadId, Transition({
                    previous.pseudoPathEntries[1],
                    current.pseudoPathEntries[0]})));

            }
        }
    }
}



void AssemblyGraph::createLinks(
    const std::map<SegmentPair, Transitions>& transitionMap,
    uint64_t minCoverage)
{
    createNew(links, "Mode3-Links");
    createNew(transitions, "Mode3-Transitions");
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
    // If this equals 1, then the sequence of segments visited by every
    // oriented read is a path in the graph.
    // But that is not desirable because of the extra edges it causes.
    const uint64_t minCoverage = 3; // EXPOSE WHEN CODE STABILIZES

    // Create a segment for each linear chain of marker graph edges.
    createSegmentPaths();
    computeSegmentCoverage();

    // Keep track of the segment and position each marker graph edge corresponds to.
    computeMarkerGraphEdgeTable(threadCount);

    // Compute pseudopaths of all oriented reads.
    // We permanently store only the compressed pseudopaths.
    computePseudoPaths(threadCount);
    computeCompressedPseudoPaths();
    pseudoPaths.remove();
    computeSegmentCompressedPseudoPathInfo();

    // Find pseudopath transitions and store them keyed by the pair of segments.
    std::map<SegmentPair, Transitions> transitionMap;
    findTransitions(transitionMap);

    // Create a links between pairs of segments with a sufficient number of transitions.
    createLinks(transitionMap, minCoverage);
    createConnectivity();
    flagBackSegments();

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
    accessExistingReadOnly(paths, "Mode3-Paths");
    accessExistingReadOnly(segmentCoverage, "Mode3-SegmentCoverage");
    accessExistingReadOnly(markerGraphEdgeTable, "Mode3-MarkerGraphEdgeTable");
    accessExistingReadOnly(compressedPseudoPaths, "Mode3-CompressedPseudoPaths");
    accessExistingReadOnly(segmentCompressedPseudoPathInfo, "Mode3-SegmentCompressedPseudoPathInfo");
    accessExistingReadOnly(links, "Mode3-Links");
    accessExistingReadOnly(transitions, "Mode3-Transitions");
    accessExistingReadOnly(linksBySource, "Mode3-LinksBySource");
    accessExistingReadOnly(linksByTarget, "Mode3-LinksByTarget");
    accessExistingReadOnly(isBackSegment, "Mode3-IsBackSegment");
    accessExistingReadOnly(clusterIds, "Mode3-ClusterIds");
}



void AssemblyGraph::createConnectivity()
{
    createNew(linksBySource, "Mode3-LinksBySource");
    createNew(linksByTarget, "Mode3-LinksByTarget");

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



// Flag back-segments.
// This does not do a full blown search for locally strongly connected components.
// A segment is marked as a back-segment if:
// - It has only a single incoming link.
// - It has a single outgoing link.
// - The incoming and outgoing links both connect to/from the same segment.
void AssemblyGraph::flagBackSegments()
{
    const uint64_t segmentCount = paths.size();
    createNew(isBackSegment, "Mode3-IsBackSegment");
    isBackSegment.resize(segmentCount);

    uint64_t backSegmentCount = 0;
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {

        // Initially flag it as not a back-segment.
        isBackSegment[segmentId] = false;

        // For a back-segment, there must be a single incoming link.
        const auto incomingLinks = linksByTarget[segmentId];
        if(incomingLinks.size() != 1) {
            continue;
        }

        // For a back-segment, there must be a single outgoing link.
        const auto outgoingLinks = linksBySource[segmentId];
        if(outgoingLinks.size() != 1) {
            continue;
        }

        // For a back-segment, the incoming and outgoing links
        // both connect to/from the same segment.
        const uint64_t incomingLinkId = incomingLinks[0];
        const uint64_t outgoingLinkId = outgoingLinks[0];
        const Link& incomingLink = links[incomingLinkId];
        const Link& outgoingLink = links[outgoingLinkId];
        if(incomingLink.segmentId0 != outgoingLink.segmentId1) {
            continue;
        }

        // Flag it as a back-segment.
        isBackSegment[segmentId] = true;
        ++backSegmentCount;
    }

    cout << "Found " << backSegmentCount << " back-segments." << endl;
}



// Get the children or parents of a given segment.
// Only use links with at least a specified coverage.
void AssemblyGraph::getChildrenOrParents(
    uint64_t segmentId,
    uint64_t direction, // 0=forward (children), 1=backward (parents).
    uint64_t minimumLinkCoverage,
    vector<uint64_t>& childrenOrParents) const
{
    switch(direction) {
    case 0:
        getChildren(segmentId, minimumLinkCoverage, childrenOrParents);
        break;
    case 1:
        getParents(segmentId, minimumLinkCoverage, childrenOrParents);
        break;
    default:
        SHASTA_ASSERT(0);
    }
}



void AssemblyGraph::getChildren(
    uint64_t segmentId,
    uint64_t minimumLinkCoverage,
    vector<uint64_t>& children) const
{
    children.clear();
    for(const auto linkId: linksBySource[segmentId]) {
        if(transitions.size(linkId) >= minimumLinkCoverage) {
            const Link& link = links[linkId];
            children.push_back(link.segmentId1);
        }
    }
}



void AssemblyGraph::getParents(
    uint64_t segmentId,
    uint64_t minimumLinkCoverage,
    vector<uint64_t>& parents) const
{
    parents.clear();
    for(const auto linkId: linksByTarget[segmentId]) {
        if(transitions.size(linkId) >= minimumLinkCoverage) {
            const Link& link = links[linkId];
            parents.push_back(link.segmentId0);
        }
    }
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
    createNew(clusterIds, "Mode3-ClusterIds");
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

    auto& threadPairs = clusterSegmentsData.threadPairs[threadId];
    threadPairs.clear();
    vector<uint64_t> descendants;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over segments assigned to this batch.
        for(uint64_t segmentId0=begin; segmentId0!=end; ++segmentId0) {

            // Add pairs for which the lowest numbered segment is segmentId0.
            addClusterPairs(threadId, segmentId0);
        }
    }
}



#if 0
void AssemblyGraph::addClusterPairs(size_t threadId, uint64_t segmentId0)
{
    // EXPOSE THESE CONSTANTS WHEN CODE STABILIZES.
    const uint64_t minCommonReadCount = 6;
    const double maxUnexplainedFraction = 0.1;

    // Loop over oriented reads in segmentId0.
    std::unordered_set<uint64_t> segmentId1s;
    for(const auto& info: segmentOrientedReadInformation[segmentId0].infos) {
        const OrientedReadId orientedReadId = info.orientedReadId;


        // Loop over the compressed pseudopath of this oriented read.
        const span<uint64_t> compressedPseudoPath = compressedPseudoPaths[orientedReadId.getValue()];
        for(const uint64_t segmentId1: compressedPseudoPath) {

            // Find each pair once.
            if(segmentId1 <= segmentId0) {
                continue;
            }
            segmentId1s.insert(segmentId1);
        }
    }



    // Check the pairs (segmentId0, segmentId1) we found.
    for(const uint64_t segmentId1: segmentId1s) {
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

        // Store it.
        clusterSegmentsData.threadPairs[threadId].push_back(make_pair(segmentId0, segmentId1));
    }
}
#else



void AssemblyGraph::addClusterPairs(size_t threadId, uint64_t startSegmentId)
{
    // EXPOSE THESE CONSTANTS WHEN CODE STABILIZES.
    const uint64_t minCommonReadCount = 6;
    const double maxUnexplainedFraction = 0.2;
    const uint64_t pairCountPerSegment = 3;
    const uint64_t maxDistance = 50;

    // std::lock_guard<std::mutex> lock(mutex);    // *********** TAKE OUT

    // Do a BFS and check each pair as we encounter it.
    // The BFS terminates when we found enough pairs.

    // Do the BFS in both directions.
    for(uint64_t direction=0; direction<1; direction++) { // ********* ONE DIRECTION ONLY
        // cout << startSegmentId << " direction " << direction << endl;

        // Initialize the BFS.
        std::queue<uint64_t> q;
        q.push(startSegmentId);
        std::map<uint64_t, uint64_t> distanceMap;
        distanceMap.insert(make_pair(startSegmentId, 0));
        uint64_t foundCount = 0;

        // BFS loop.
        while(not q.empty()) {
            const uint64_t segmentId0 = q.front();
            // cout << "Dequeued " << segmentId0 << endl;
            q.pop();

            const uint64_t distance0 = distanceMap[segmentId0];
            const uint64_t distance1 = distance0 + 1;

            // Loop over children or parents of segmentId0.
            const auto neighbors = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
            for(const uint64_t linkId01: neighbors) {
                const Link& link01 = links[linkId01];
                const uint64_t segmentId1 = (direction==0) ? link01.segmentId1 : link01.segmentId0;

                // If we already encountered segmentId1, skip it.
                if(distanceMap.find(segmentId1) != distanceMap.end()) {
                    continue;
                }

                // Enqueue it.
                if(distance1 < maxDistance) {
                    q.push(segmentId1);
                }
                distanceMap.insert(make_pair(segmentId1, distance1));

                // cout << "Found " << segmentId1 << endl;

                // Check the pair (startSegmentId, segmentId1).
                SegmentPairInformation info;
                analyzeSegmentPair(startSegmentId, segmentId1,
                    segmentOrientedReadInformation[startSegmentId],
                    segmentOrientedReadInformation[segmentId1],
                    markers, info);
                if(info.commonCount < minCommonReadCount) {
                    continue;
                }
                if(info.maximumUnexplainedFraction() > maxUnexplainedFraction) {
                    continue;
                }

                // Store it.
                // cout << "Stored " << segmentId1 << endl;
                clusterSegmentsData.threadPairs[threadId].push_back(make_pair(startSegmentId, segmentId1));
                ++foundCount;
                if(foundCount >= pairCountPerSegment) {
                    break;
                }
            }

            if(foundCount >= pairCountPerSegment) {
                break;
            }
        }
    }
}
#endif


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



// Analyze a subgraph of the assembly graph.
void AssemblyGraph::analyzeSubgraph1(
    const vector<uint64_t>& unsortedSegmentIds,
    vector<AnalyzeSubgraphClasses::Cluster>& clusters,
    bool debug) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minClusterCoverage = 3;

    using CompressedPseudoPathSnippet = AnalyzeSubgraphClasses::CompressedPseudoPathSnippet;
    using Cluster = AnalyzeSubgraphClasses::Cluster;

    // Create a sorted version of the segmentIds. We will need it later.
    vector<uint64_t> segmentIds = unsortedSegmentIds;
    sort(segmentIds.begin(), segmentIds.end());

    // Gather triplets (orientedReadId, position in compressedPseudoPath, segmentId).
    using Triplet = tuple<OrientedReadId, uint64_t, uint64_t>;
    vector<Triplet> triplets;
    for(const uint64_t segmentId: segmentIds) {
        const auto v = segmentCompressedPseudoPathInfo[segmentId];
        for(const auto& p: v) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t position = p.second;
            triplets.push_back(Triplet(orientedReadId, position, segmentId));
        }
    }
    sort(triplets.begin(), triplets.end());

    // Write the triplets.
    if(debug) {
        ofstream csv("Triplets.csv");
        for(const Triplet& triplet: triplets) {
            csv << get<0>(triplet) << ",";
            csv << get<1>(triplet) << ",";
            csv << get<2>(triplet) << "\n";
        }
    }



    // Find streaks for the same OrientedReadId where the position
    // increases by 1 each time.
    // Each streak generates a CompressedPseudoPathSnippet.
    vector<CompressedPseudoPathSnippet> snippets;
    for(uint64_t i=0; i<triplets.size(); /* Increment later */) {
        const OrientedReadId orientedReadId = get<0>(triplets[i]);

        // Find this streak.
        uint64_t streakBegin = i;
        uint64_t streakEnd = streakBegin + 1;
        for(; streakEnd<triplets.size(); streakEnd++) {
            if(get<0>(triplets[streakEnd]) != orientedReadId) {
                break;
            }
            if(get<1>(triplets[streakEnd]) != get<1>(triplets[streakEnd-1]) + 1) {
                break;
            }
        }

        // Add a snippet.
        CompressedPseudoPathSnippet snippet;
        snippet.orientedReadId = orientedReadId;
        snippet.firstPosition = get<1>(triplets[streakBegin]);
        for(uint64_t j=streakBegin; j!=streakEnd; ++j) {
            snippet.segmentIds.push_back(get<2>(triplets[j]));
        }
        snippets.push_back(snippet);

        // Prepare to process the next streak.
        i = streakEnd;
    }



    // Write the snippets.
    if(debug) {
        ofstream csv("CompressedPseudoPathSnippets.csv");
        csv << "SnippetIndex,OrientedReadId,First position,LastPosition,SegmentIds\n";
        for(uint64_t snippetIndex=0; snippetIndex<snippets.size(); snippetIndex++) {
            const CompressedPseudoPathSnippet& snippet = snippets[snippetIndex];
            csv << snippetIndex << ",";
            csv << snippet.orientedReadId << ",";
            csv << snippet.firstPosition << ",";
            csv << snippet.lastPosition() << ",";
            for(const uint64_t segmentId: snippet.segmentIds) {
                csv << segmentId << ",";
            }
            csv << "\n";
        }
    }



    // Gather the snippets into groups with the same segmentIds.
    // Key: segmentIds.
    // Value: vector of indexes into the snippets vector.
    std::map<vector<uint64_t>, vector<uint64_t> > snippetMap;
    for(uint64_t snippetIndex=0; snippetIndex<snippets.size(); snippetIndex++) {
        const CompressedPseudoPathSnippet& snippet = snippets[snippetIndex];
        snippetMap[snippet.segmentIds].push_back(snippetIndex);
    }
    vector< pair < vector<uint64_t>, vector<uint64_t> > > snippetGroups(snippetMap.size());
    copy(snippetMap.begin(), snippetMap.end(),snippetGroups.begin());
    for(uint64_t snippetGroupIndex=0; snippetGroupIndex<snippetGroups.size(); snippetGroupIndex++) {
        const auto& p = snippetGroups[snippetGroupIndex];

        if(debug) {
            const vector<uint64_t>& segmentIds = p.first;
            const vector<uint64_t>& snippetIndexes = p.second;
            cout << "Snippet group " << snippetGroupIndex <<
                " has " <<  snippetIndexes.size() << " snippets with " << segmentIds.size() << " segments." << endl;
            cout << "    Snippets: ";
            copy(snippetIndexes.begin(), snippetIndexes.end(), ostream_iterator<uint64_t>(cout, " "));
            cout << endl;
            cout << "    Segments: ";
            copy(segmentIds.begin(), segmentIds.end(), ostream_iterator<uint64_t>(cout, " "));
            cout << endl;
        }
    }



    // Create a subset graph for the groups.
    using SegmentSubsetGraph = SubsetGraph<uint64_t>;
    SegmentSubsetGraph subsetGraph;
    using vertex_descriptor = SegmentSubsetGraph::vertex_descriptor;
    vector<vertex_descriptor> subsetGraphVertices;
    std::map<vertex_descriptor, uint64_t> indexMap;
    for(uint64_t snippetGroupIndex=0; snippetGroupIndex<snippetGroups.size(); snippetGroupIndex++) {
        const auto& p = snippetGroups[snippetGroupIndex];
        const vector<uint64_t>& segmentIds = p.first;
        const SubsetGraph<uint64_t>::vertex_descriptor v = subsetGraph.addVertex(segmentIds);
        subsetGraphVertices.push_back(v);
        indexMap.insert(make_pair(v, snippetGroupIndex));
    }
    subsetGraph.createEdges();



    // Find the maximal snippet groups.
    vector<uint64_t> maximalSnippetGroups;
    for(uint64_t snippetGroupIndex=0; snippetGroupIndex<snippetGroups.size(); snippetGroupIndex++) {
        const auto v = subsetGraphVertices[snippetGroupIndex];
        if(in_degree(v, subsetGraph) == 0) {
            maximalSnippetGroups.push_back(snippetGroupIndex);
        }
    }
    if(debug) {
        cout << "Found " << maximalSnippetGroups.size() << " maximal snippet groups:" << endl;
        for(const uint64_t snippetGroupIndex: maximalSnippetGroups) {
            const auto& p = snippetGroups[snippetGroupIndex];
            const vector<uint64_t>& segmentIds = p.first;
            cout << "Group " << snippetGroupIndex <<
                " with " << p.second.size() <<
                " snippets and " << segmentIds.size() << " segments:" << endl;
            copy(segmentIds.begin(), segmentIds.end(), ostream_iterator<uint64_t>(cout, " "));
            cout << endl;
        }
    }



    // Find the descendants of each maximal snippet group.
    // Keep track of which maximal groups each vertex is a descendant of.
    std::map<vertex_descriptor, vector<vertex_descriptor> > ancestorMap;
    for(const uint64_t snippetGroupIndex: maximalSnippetGroups) {
        const vertex_descriptor vStart = subsetGraphVertices[snippetGroupIndex];

        // Do a BFS starting here.
        std::queue<vertex_descriptor> q;
        q.push(vStart);
        std::set<vertex_descriptor> descendants;
        descendants.insert(vStart);
        while(not q.empty()) {
            const vertex_descriptor v0 = q.front();
            q.pop();
            BGL_FORALL_OUTEDGES(v0, e01, subsetGraph, SegmentSubsetGraph) {
                const vertex_descriptor v1 = target(e01, subsetGraph);
                if(descendants.find(v1) == descendants.end()) {
                    q.push(v1);
                    descendants.insert(v1);
                }
            }
        }
        for(const vertex_descriptor v: descendants) {
            ancestorMap[v].push_back(vStart);
        }

        if(debug) {
            cout << "Descendants of " << snippetGroupIndex << ":" << endl;
            for(const vertex_descriptor v: descendants) {
                cout << indexMap[v] << " ";
            }
            cout << endl;
        }
    }

    if(debug) {
        cout << "Ancestor map:" << endl;
        for(uint64_t snippetGroupIndex=0; snippetGroupIndex<snippetGroups.size(); snippetGroupIndex++) {
            const vertex_descriptor v0 = subsetGraphVertices[snippetGroupIndex];
            cout << indexMap[v0] << ": ";
            for(const vertex_descriptor v1: ancestorMap[v0]) {
                cout << indexMap[v1] << " ";
            }
            cout << endl;
        }
    }



    // Each maximal snippet group, together with its descendants
    // that are not also descendants of another snippet group, generates a cluster.
    std::map<uint64_t, vector<uint64_t> > clusterMap;
    for(uint64_t snippetGroupIndex=0; snippetGroupIndex<snippetGroups.size(); snippetGroupIndex++) {
        const vertex_descriptor v0 = subsetGraphVertices[snippetGroupIndex];
        const vector<vertex_descriptor>& maximalAncestors = ancestorMap[v0];
        SHASTA_ASSERT(not maximalAncestors.empty());
        if(maximalAncestors.size() == 1) {
            const vertex_descriptor maximalAncestor = maximalAncestors.front();
            clusterMap[indexMap[maximalAncestor]].push_back(snippetGroupIndex);

        } else {
            if(debug) {
                cout << "Snippet group " << indexMap[v0] <<
                    " is ambiguous and was not clustered." << endl;
            }
        }
    }



    // Gather the snippets in each cluster.
    clusters.clear();
    for(const auto& p: clusterMap) {
        clusters.resize(clusters.size() + 1);
        Cluster& cluster = clusters.back();

        // Get the snippet groups in this cluster.
        const vector<uint64_t>& snippetGroupIndexes = p.second;
        SHASTA_ASSERT(not snippetGroupIndexes.empty());
        cluster.snippetGroupIndexes = snippetGroupIndexes;

        // Loop over the snippet groups in this cluster.
        for(const uint64_t snippetGroupIndex: snippetGroupIndexes) {
            const auto& q = snippetGroups[snippetGroupIndex];
            const vector<uint64_t>& segmentIds = q.first;
            const vector<uint64_t>& snippetIndexes = q.second;

            // Loop over the snippets in this snippet group.
            for(const uint64_t snippetIndex: snippetIndexes) {
                const CompressedPseudoPathSnippet& snippet = snippets[snippetIndex];
                SHASTA_ASSERT(snippet.segmentIds == segmentIds);
                cluster.snippets.push_back(snippet);
            }
        }

        cluster.constructSegments();
        cluster.cleanupSegments(minClusterCoverage);

        // If coverage on this cluster is too low, discard it.
        if(cluster.coverage() < minClusterCoverage) {
            clusters.resize(clusters.size() - 1);
            continue;
        }

        if(debug) {
            cout << "Cluster " << clusters.size() - 1 << " coverage " << cluster.snippets.size() << endl;
            cout << "    Snippet groups:\n    ";
            copy(snippetGroupIndexes.begin(), snippetGroupIndexes.end(),
                ostream_iterator<uint64_t>(cout, " "));
            cout << endl;
            cout << "    Segments with coverage:\n    ";
            for(const auto& p: cluster.segments) {
                cout << p.first << "(" << p.second << ") ";
            }
            cout << endl;
        }
    }



    // Due to errors in the reads, at this point we can have some
    // clusters very similar to each other.
    // Check each pair.
    // This could be made more efficient.
    vector<uint64_t> workArea;
    for(uint64_t i0=0; i0<clusters.size()-1; i0++) {
        const Cluster& cluster0 = clusters[i0];
        const vector<uint64_t> segments0 = cluster0.getSegments();
        for(uint64_t i1=i0+1; i1<clusters.size(); i1++) {
            const Cluster& cluster1 = clusters[i1];
            const vector<uint64_t> segments1 = cluster1.getSegments();

            // Intersection.
            workArea.clear();
            set_intersection(
                segments0.begin(), segments0.end(),
                segments1.begin(), segments1.end(),
                back_inserter(workArea));
            const uint64_t intersectionSize = workArea.size();

            // 0 - 1.
            workArea.clear();
            set_difference(
                segments0.begin(), segments0.end(),
                segments1.begin(), segments1.end(),
                back_inserter(workArea));
            const uint64_t difference01 = workArea.size();

            // 1 - 0.
            workArea.clear();
            set_difference(
                segments1.begin(), segments1.end(),
                segments0.begin(), segments0.end(),
                back_inserter(workArea));
            const uint64_t difference10 = workArea.size();

            const uint64_t unionSize = intersectionSize + difference01 + difference10;
            const double jaccard = double(intersectionSize) / double(unionSize);
            const double ratio01 = double(difference01) / double(unionSize);
            const double ratio10 = double(difference10) / double(unionSize);

            cout << i0 << " " << i1 << ": j=";
            cout << jaccard << ", ratio01=" << ratio01 << ", ratio10=" << ratio10 << endl;



        }
    }





    // Write the subset graph in Graphviz format.
    if(debug) {

        // A vector to contain the cluster that each snippet group
        // was assigned to, if any.
        vector<uint64_t> snippetGroupCluster(snippetGroups.size(), std::numeric_limits<uint64_t>::max());
        for(uint i=0; i<clusters.size(); i++) {
            const Cluster& cluster = clusters[i];
            for(const uint64_t snippetGroupIndex: cluster.snippetGroupIndexes) {
                snippetGroupCluster[snippetGroupIndex] = i;
            }
        }

        ofstream dot("SubsetGraph.dot");
        dot << "digraph SubsetGraph {\n"
            "node [shape=rectangle];\n";
        for(uint64_t snippetGroupIndex=0; snippetGroupIndex<snippetGroups.size(); snippetGroupIndex++) {
            const auto& p = snippetGroups[snippetGroupIndex];
            const uint64_t clusterId = snippetGroupCluster[snippetGroupIndex];
            const vector<uint64_t>& segmentIds = p.first;
            dot << snippetGroupIndex << " [label=\""
                "Group " << snippetGroupIndex << "\\n";
            if(clusterId != std::numeric_limits<uint64_t>::max()) {
                dot << "Cluster " << clusterId << "\\n";
            }
            dot << "Coverage " << p.second.size() << "\\n";
            for(uint64_t i=0; i<segmentIds.size(); i++) {
                dot << segmentIds[i];
                if(i != segmentIds.size()-1) {
                    dot << "\\n";
                }
            }
            dot << "\"";
            if(clusterId != std::numeric_limits<uint64_t>::max()) {
                dot << " style=filled fillcolor=\"" <<
                    float(clusterId)/float(clusters.size()) <<
                    ",0.3,1\"";
            }
            dot << "];\n";
        }
        BGL_FORALL_EDGES(e, subsetGraph, SegmentSubsetGraph) {
            const auto v0 = source(e, subsetGraph);
            const auto v1 = target(e, subsetGraph);
            dot << indexMap[v0] << "->" << indexMap[v1] << ";\n";
        }
        dot << "}\n";
    }




#if 0
    // Create a feature matrix that gives, for each snippet, a bit vector
    // of the segments it visits.
    vector< vector<bool> > featureMatrix(snippets.size(), vector<bool>(segmentIds.size(), false));
    for(uint64_t i=0; i<snippets.size(); i++) {
        const CompressedPseudoPathSnippet& snippet = snippets[i];
        for(const uint64_t segmentId: snippet.segmentIds) {
            const auto it = lower_bound(segmentIds.begin(), segmentIds.end(), segmentId);
            SHASTA_ASSERT(it != segmentIds.end());
            SHASTA_ASSERT(*it == segmentId);
            const uint64_t j = it - segmentIds.begin();
            featureMatrix[i][j] = true;
        }
    }

    // Write the feature matrix.
    {
        ofstream csv("FeatureMatrix.csv");
        csv << "OrientedReadId,";
        for(const uint64_t segmentId: segmentIds) {
            csv << "S" << segmentId << ",";
        }
        csv << "\n";
        for(uint64_t i=0; i<snippets.size(); i++) {
            csv << snippets[i].orientedReadId << ",";
            for(uint64_t j=0; j<segmentIds.size(); j++) {
                csv << int(featureMatrix[i][j]) << ",";
            }
            csv << "\n";
        }
    }


    // Compute the SVD of the feature matrix,
    const string JOBU = "A";
    const string JOBVT = "A";
    const int M = int(snippets.size());         // Number of rows. one row per snippet.
    const int N = int(segmentIds.size());       // Number of columns. One column pwer segment.
    vector<double> A(M*N);                      // Fortran storage (by column).
    for(uint64_t i=0; i<snippets.size(); i++) {
        for(uint64_t j=0; j<segmentIds.size(); j++) {
            A[j*M + i] = featureMatrix[i][j];
        }
    }
    const int LDA = M;
    vector<double> S(min(M, N));
    vector<double> U(M*M);
    const int LDU = M;
    vector<double> VT(N*N);
    const int LDVT = N;
    const int LWORK = 10 * max(M, N);
    vector<double> WORK(LWORK);
    int INFO = 0;

    dgesvd_(
        JOBU.data(), JOBVT.data(),
        M, N,
        &A[0], LDA, &S[0], &U[0], LDU, &VT[0], LDVT, &WORK[0], LWORK, INFO);
    cout << "dgesvd return code " << INFO << endl;

    if(INFO != 0) {
        return;
    }

    cout << "Singular values: " << endl;
    for(const double v: S) {
        cout << v << endl;
    }

    // Write the first 2 left eigenvectors (U, of length M).
    // There is one for each snippet.
    {
        ofstream csv("LeftEigenvectors.csv");
        for(uint64_t i=0; i<snippets.size(); i++) {
            csv << snippets[i].orientedReadId << ",";
            csv << U[i] << "," << U[i+M] << endl;
        }
    }
#endif
}


void AssemblyGraph::analyzeSubgraph2(
    const vector<uint64_t>& segmentIds,
    vector<AnalyzeSubgraphClasses::Cluster>& clusters,
    bool debug) const
{
    if(segmentIds.size() <= 64) {
        analyzeSubgraph2Template<64>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 128) {
        analyzeSubgraph2Template<128>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 192) {
        analyzeSubgraph2Template<192>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 256) {
        analyzeSubgraph2Template<256>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 320) {
        analyzeSubgraph2Template<320>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 384) {
        analyzeSubgraph2Template<384>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 448) {
        analyzeSubgraph2Template<448>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 512) {
        analyzeSubgraph2Template<512>(segmentIds, clusters, debug);
    } else {
        SHASTA_ASSERT(0);
    }
}



template<uint64_t N> void AssemblyGraph::analyzeSubgraph2Template(
    const vector<uint64_t>& unsortedSegmentIds,
    vector<AnalyzeSubgraphClasses::Cluster>& clusters,
    bool debug) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const double fractionThreshold = 0.05;
    const uint64_t minClusterCoverage = 6;
    const uint64_t minSegmentCoverage = 6;

    using BitVector = std::bitset<N>;
    using CompressedPseudoPathSnippet = AnalyzeSubgraphClasses::CompressedPseudoPathSnippet;
    using Cluster = AnalyzeSubgraphClasses::Cluster;
    using SnippetGraphVertex = AnalyzeSubgraphClasses::SnippetGraphVertex;
    using SnippetGraph = AnalyzeSubgraphClasses::SnippetGraph;
    using vertex_descriptor = SnippetGraph::vertex_descriptor;

    // Create a sorted version of the segmentIds. We will need it later.
    vector<uint64_t> segmentIds = unsortedSegmentIds;
    sort(segmentIds.begin(), segmentIds.end());

    // Gather triplets (orientedReadId, position in compressedPseudoPath, segmentId).
    using Triplet = tuple<OrientedReadId, uint64_t, uint64_t>;
    vector<Triplet> triplets;
    for(const uint64_t segmentId: segmentIds) {
        const auto v = segmentCompressedPseudoPathInfo[segmentId];
        for(const auto& p: v) {
            const OrientedReadId orientedReadId = p.first;
            const uint64_t position = p.second;
            triplets.push_back(Triplet(orientedReadId, position, segmentId));
        }
    }
    sort(triplets.begin(), triplets.end());

    // Write the triplets.
    if(debug) {
        ofstream csv("Triplets.csv");
        for(const Triplet& triplet: triplets) {
            csv << get<0>(triplet) << ",";
            csv << get<1>(triplet) << ",";
            csv << get<2>(triplet) << "\n";
        }
    }



    // Find streaks for the same OrientedReadId where the position
    // increases by 1 each time.
    // Each streak generates a CompressedPseudoPathSnippet.
    vector<CompressedPseudoPathSnippet> snippets;
    for(uint64_t i=0; i<triplets.size(); /* Increment later */) {
        const OrientedReadId orientedReadId = get<0>(triplets[i]);

        // Find this streak.
        uint64_t streakBegin = i;
        uint64_t streakEnd = streakBegin + 1;
        for(; streakEnd<triplets.size(); streakEnd++) {
            if(get<0>(triplets[streakEnd]) != orientedReadId) {
                break;
            }
            if(get<1>(triplets[streakEnd]) != get<1>(triplets[streakEnd-1]) + 1) {
                break;
            }
        }

        // Add a snippet.
        CompressedPseudoPathSnippet snippet;
        snippet.orientedReadId = orientedReadId;
        snippet.firstPosition = get<1>(triplets[streakBegin]);
        for(uint64_t j=streakBegin; j!=streakEnd; ++j) {
            snippet.segmentIds.push_back(get<2>(triplets[j]));
        }
        snippets.push_back(snippet);

        // Prepare to process the next streak.
        i = streakEnd;
    }



    // Write the snippets.
    if(debug) {
        ofstream csv("CompressedPseudoPathSnippets.csv");
        csv << "SnippetIndex,OrientedReadId,First position,LastPosition,SegmentIds\n";
        for(uint64_t snippetIndex=0; snippetIndex<snippets.size(); snippetIndex++) {
            const CompressedPseudoPathSnippet& snippet = snippets[snippetIndex];
            csv << snippetIndex << ",";
            csv << snippet.orientedReadId << ",";
            csv << snippet.firstPosition << ",";
            csv << snippet.lastPosition() << ",";
            for(const uint64_t segmentId: snippet.segmentIds) {
                csv << segmentId << ",";
            }
            csv << "\n";
        }
    }



    // For each snippet, create a BitVector that describes the segments
    // the snippet visits.
    const uint64_t snippetCount = snippets.size();
    vector<BitVector> bitVectors(snippetCount);
    vector<uint64_t> bitVectorsPopCount(snippetCount);
    for(uint64_t snippetIndex=0; snippetIndex<snippetCount; snippetIndex++) {
        const CompressedPseudoPathSnippet& snippet = snippets[snippetIndex];
        BitVector& bitVector = bitVectors[snippetIndex];

        for(const uint64_t segmentId: snippet.segmentIds) {
            auto it = lower_bound(segmentIds.begin(), segmentIds.end(), segmentId);
            SHASTA_ASSERT(it != segmentIds.end());
            SHASTA_ASSERT(*it == segmentId);
            const uint64_t bitIndex = it - segmentIds.begin();
            bitVector.set(bitIndex);
        }
        bitVectorsPopCount[snippetIndex] = bitVector.count();
    }



    // Create the SnippetGraph.
    // A vertex represents a set of snippets and stores
    // the corresponding snippet indexes.
    // An edge x->y is created if there is at least one snippet in y
    // that is an approximate subset of a snippet in x.
    // We express this condition as |y-x| < fractionThreshold * |y|
    // We start with one snippet per vertex.
    SnippetGraph graph;
    vector<vertex_descriptor> vertexTable;
    std::map<vertex_descriptor, uint64_t> vertexMap;
    for(uint64_t snippetIndex=0; snippetIndex<snippetCount; snippetIndex++) {
        const auto v = add_vertex(SnippetGraphVertex(snippetIndex), graph);
        vertexTable.push_back(v);
        vertexMap.insert(make_pair(v, snippetIndex));
    }
    for(uint64_t iy=0; iy<snippetCount; iy++) {
        const BitVector& y = bitVectors[iy];
        const uint64_t threshold = uint64_t(std::round(fractionThreshold * double(bitVectorsPopCount[iy])));
        const vertex_descriptor vy = vertexTable[iy];
        for(uint64_t ix=0; ix<snippetCount; ix++) {
            if(ix == iy) {
                continue;
            }
            const BitVector& x = bitVectors[ix];

            // Compute z = y-x.
            BitVector z = y;
            z &= (~x);

            if(z.count() <= threshold) {
                const vertex_descriptor vx = vertexTable[ix];
                add_edge(vx, vy, graph);
            }
        }
    }



    // Compute strongly connected components of the SnippetGraph.
    std::map<vertex_descriptor, uint64_t> componentMap;
    const uint64_t componentCount = boost::strong_components(
        graph,
        boost::make_assoc_property_map(componentMap),
        boost::vertex_index_map(boost::make_assoc_property_map(vertexMap)));
    // cout << "Found " << componentCount << " strongly connected components." << endl;

    // Gather the vertices of each strongly connected component.
    vector< vector<vertex_descriptor> > components(componentCount);
    BGL_FORALL_VERTICES_T(v, graph, SnippetGraph) {
        const uint64_t componentId = componentMap[v];
        SHASTA_ASSERT(componentId < componentCount);
        components[componentId].push_back(v);
    }
    if(false) {
        cout << "Strongly connected components:\n";
        for(uint64_t componentId=0; componentId<componentCount; componentId++) {
            cout << componentId << ": ";
            for(const vertex_descriptor v: components[componentId]) {
                cout << vertexMap[v] << " ";
            }
            cout << "\n";
        }
    }



    // Condense the strongly connected components.
    // After this, the SnippetGraph is guaranteed to be acyclic.
    for(const vector<vertex_descriptor>& component: components) {
        if(component.size() == 1) {
            continue;
        }

        // Create a new vertex to represent this component.
        const auto vNew = add_vertex(graph);
        vector<uint64_t>& snippetsNew = graph[vNew].snippetIndexes;
        for(const vertex_descriptor v: component) {
            const vector<uint64_t>& snippets = graph[v].snippetIndexes;
            SHASTA_ASSERT(snippets.size() == 1);
            snippetsNew.push_back(snippets.front());
        }

        // Create the new edges.
        for(const vertex_descriptor v0: component) {

            // Out-edges.
            BGL_FORALL_OUTEDGES_T(v0, e01, graph, SnippetGraph) {
                const vertex_descriptor v1 = target(e01, graph);
                if(v1 != vNew) {
                    add_edge(vNew, v1,graph);
                }
            }

            // In-edges.
            BGL_FORALL_INEDGES_T(v0, e10, graph, SnippetGraph) {
                const vertex_descriptor v1 = source(e10, graph);
                if(v1 != vNew) {
                    add_edge(v1, vNew, graph);
                }
            }
        }

        // Remove the old vertices and their edges.
        for(const vertex_descriptor v: component) {
            clear_vertex(v, graph);
            remove_vertex(v, graph);
        }
    }


    // Compute which maximal vertices each vertex is a descendant of.
    std::map<vertex_descriptor, vector<vertex_descriptor> > ancestorMap;
    BGL_FORALL_VERTICES_T(v0, graph, SnippetGraph) {
        if(in_degree(v0, graph) != 0) {
            continue;   // Not a maximal vertex.
        }

        // Find the descendants of this maximal vertex.
        vector<vertex_descriptor> descendants;
        graph.findDescendants(v0, descendants);

        // Update the ancestor map.
        for(const vertex_descriptor v1: descendants) {
            ancestorMap[v1].push_back(v0);
        }
    }



    // Each maximal vertex generates a cluster consisting of the vertices
    // that descend from it and from no other maximal vertex.
    // Gather the vertices in each cluster.
    std::map<vertex_descriptor, vector<vertex_descriptor> > clusterMap;
    uint64_t unclusterVertexCount = 0;
    BGL_FORALL_VERTICES_T(v1, graph, SnippetGraph) {
        const vector<vertex_descriptor>& ancestors = ancestorMap[v1];
        if(ancestors.size() == 1) {
            const vertex_descriptor v0 = ancestors.front();
            clusterMap[v0].push_back(v1);
        } else {
            ++unclusterVertexCount;
        }
    }
    cout << "Found " << unclusterVertexCount << " unclustered vertices." << endl;



    // Gather the snippets in each cluster.
    clusters.clear();
    for(const auto& p: clusterMap) {
        const vector<vertex_descriptor>& clusterVertices = p.second;
        clusters.resize(clusters.size() + 1);
        Cluster& cluster = clusters.back();

        vector<uint64_t> clusterSnippetIndexes; // Only used for debug output.
        for(const vertex_descriptor v: clusterVertices) {
            const vector<uint64_t>& snippetIndexes = graph[v].snippetIndexes;
            for(const uint64_t snippetIndex: snippetIndexes) {
                cluster.snippets.push_back(snippets[snippetIndex]);
                clusterSnippetIndexes.push_back(snippetIndex);
            }
        }
        cluster.constructSegments();
        cluster.cleanupSegments(minSegmentCoverage);
        cout << "Found a cluster candidate with " <<
            clusterVertices.size() << " vertices and " <<
            cluster.snippets.size() << " snippets:" << endl;
        for(const uint64_t snippetIndex: clusterSnippetIndexes) {
            cout << snippetIndex << " ";
        }
        cout << endl;

        // If coverage on this cluster is too low, discard it.
        if(cluster.coverage() < minClusterCoverage) {
            clusters.resize(clusters.size() - 1);
            cout << "This cluster candidate was discarded because of low coverage." << endl;
            continue;
        }

        // This cluster will be stored and is assigned this clusterId;
        const uint64_t clusterId = clusters.size() - 1;

        if(debug) {

            cout << "This cluster was stored as cluster " << clusterId << endl;
            cout << "Segment(coverage) for this cluster:\n";
            for(const auto& p: cluster.segments) {
                cout << p.first << "(" << p.second << ") ";
            }
            cout << endl;
        }

        // Mark the vertices of this cluster.
        for(const vertex_descriptor v: clusterVertices) {
            graph[v].clusterId = clusterId;
        }
    }
    graph.clusterCount = clusters.size();



    // Write out the SnippetGraph.
    if(debug) {
        graph.writeGraphviz("SnippetGraph.dot");
    }
}



void AssemblyGraph::AnalyzeSubgraphClasses::SnippetGraph::findDescendants(
    const vertex_descriptor vStart,
    vector<vertex_descriptor>& descendants) const
{
    const SnippetGraph& graph = *this;

    // Initialize the BFS.
    std::queue<vertex_descriptor> q;
    q.push(vStart);
    std::set<vertex_descriptor> descendantsSet;
    descendantsSet.insert(vStart);

    // BFS loop.
    while(not q.empty()) {
        const vertex_descriptor v0 = q.front();
        q.pop();

        BGL_FORALL_OUTEDGES(v0, e01, graph, SnippetGraph) {
            const vertex_descriptor v1 = target(e01, graph);
            if(descendantsSet.find(v1) == descendantsSet.end()) {
                q.push(v1);
                descendantsSet.insert(v1);
            }
        }
    }

    descendants.clear();
    copy(descendantsSet.begin(), descendantsSet.end(), back_inserter(descendants));
}



void AssemblyGraph::AnalyzeSubgraphClasses::SnippetGraph::writeGraphviz(
    const string& fileName) const
{
    const SnippetGraph& graph = *this;

    ofstream dot(fileName);
    dot << "digraph SnippetGraph{\n"
        "node [shape=rectangle];\n";
    BGL_FORALL_VERTICES(v, graph, SnippetGraph) {
        dot << "\"" << v << "\" [label=\"";
        const vector<uint64_t>& snippetIndexes = graph[v].snippetIndexes;
        for(const uint64_t snippetIndex: snippetIndexes) {
            dot << snippetIndex;
            dot << "\\n";
        }
        dot << "\"";
        const uint64_t clusterId = graph[v].clusterId;
        if(clusterId != std::numeric_limits<uint64_t>::max()) {
            dot << " style=filled fillcolor=\"" <<
                float(clusterId)/float(clusterCount) <<
                ",0.3,1\"";
        }
        dot << "];\n";
    }
    BGL_FORALL_EDGES(e, graph, SnippetGraph) {
        const vertex_descriptor vx = source(e, graph);
        const vertex_descriptor vy = target(e, graph);
        dot << "\"" << vx << "\"->\"" << vy << "\";\n";
    }
    dot << "}\n";

}



void AssemblyGraph::AnalyzeSubgraphClasses::Cluster::constructSegments()
{
    // A map with Key=segmentId, value = coverage.
    std::map<uint64_t, uint64_t> segmentMap;

    for(const CompressedPseudoPathSnippet& snippet: snippets) {
        for(const uint64_t segmentId: snippet.segmentIds) {
            auto it = segmentMap.find(segmentId);
            if(it == segmentMap.end()) {
                segmentMap.insert(make_pair(segmentId, 1));
            } else {
                ++(it->second);
            }
        }
    }

    segments.clear();
    copy(segmentMap.begin(), segmentMap.end(), back_inserter(segments));
}



void AssemblyGraph::AnalyzeSubgraphClasses::Cluster::cleanupSegments(uint64_t minSegmentCoverage)
{
    vector< pair<uint64_t, uint64_t > > newSegments;
    for(const auto& p: segments) {
        if(p.second >= minSegmentCoverage) {
            newSegments.push_back(p);
        }
    }
    segments.swap(newSegments);
}



vector<uint64_t> AssemblyGraph::AnalyzeSubgraphClasses::Cluster::getSegments() const
{
    vector<uint64_t> v;
    for(const auto& p: segments) {
        v.push_back(p.first);
    }
    return v;
}



// Create an assembly path starting at a given segment.
void AssemblyGraph::createAssemblyPath(
    uint64_t startSegmentId,
    uint64_t direction,    // 0 = forward, 1 = backward
    vector<uint64_t>& path // The segmentId's of the path.
    ) const
{
    // CONSTANTS TO EXPOSE WHEN CODE STABILIZES.
    const uint64_t minimumLinkCoverage = 4;
    const double minJaccardForReference = 0.8;

    const bool debug = true;
    if(debug) {
        cout << "Creating an assembly path starting at segment " << startSegmentId <<
            ", direction "<< direction << endl;
    }

    // Start with a path consisting only of startSegmentId.
    path.clear();
    if(isBackSegment[startSegmentId]) {
        return;
    }
    path.push_back(startSegmentId);

    // The last segment that was added to the path.
    uint64_t segmentId0 = startSegmentId;

    // Work areas used in the loop below and defined here
    // to reduce memory allocation activity.
    vector<uint64_t> childrenOrParents;
    SegmentOrientedReadInformation referenceOrientedReads;
    SegmentOrientedReadInformation orientedReads1;
    vector<SegmentPairInformation> segmentPairInformation;
    vector<uint64_t> newReferenceCandidates;
    vector<uint64_t> nextSegmentCandidates;

    // The current reference segment.
    uint64_t referenceSegmentId = startSegmentId;
    getOrientedReadsOnSegment(referenceSegmentId, referenceOrientedReads);

    // Main loop. At each step we append one segment to the path.
    while(true) {

        if(debug) {
            cout << "Loop iteration begins: "
                "reference segment " << referenceSegmentId <<
                ", segmentId0 " << segmentId0 << endl;
        }

        // Loop over children or parents of segmentId0.
        // For each, compute SegmentPairInformation relative to
        // our current reference segment.
        getChildrenOrParents(segmentId0, direction, minimumLinkCoverage, childrenOrParents);
        if(childrenOrParents.empty()) {
            if(debug) {
                cout << "There are no children or parents." << endl;
            }
            return;
        }
        segmentPairInformation.resize(childrenOrParents.size());
        for(uint64_t i=0; i<childrenOrParents.size(); i++) {
            const uint64_t segmentId1 = childrenOrParents[i];
            if(isBackSegment[segmentId1]) {
                continue;
            }
            getOrientedReadsOnSegment(segmentId1, orientedReads1);
            analyzeSegmentPair(referenceSegmentId, segmentId1,
                referenceOrientedReads, orientedReads1,
                markers, segmentPairInformation[i]);

            if(debug) {
                cout << "Segment " << segmentId1 << ": ";
                if(segmentPairInformation[i].commonCount == 0) {
                    cout << "no common reads";
                } else {
                    cout << "Jaccard " << segmentPairInformation[i].jaccard() <<
                    ", unexplained fractions " <<
                    segmentPairInformation[i].unexplainedFraction(0) << " " <<
                    segmentPairInformation[i].unexplainedFraction(1);
                }
                cout << endl;
            }
        }



        // Find the child or parent that has the best Jaccard similarity with
        // the current reference segment.
        uint64_t iBest = std::numeric_limits<uint64_t>::max();
        double jaccardBest = -1.;
        for(uint64_t i=0; i<childrenOrParents.size(); i++) {
            if(isBackSegment[childrenOrParents[i]]) {
                continue;
            }
            const double jaccard = segmentPairInformation[i].jaccard();
            if(jaccard > jaccardBest) {
                iBest = i;
                jaccardBest = jaccard;
            }
        }
        if(iBest == std::numeric_limits<uint64_t>::max()) { // Can happen due to back-segments.
            return;
        }
        if(debug) {
            cout << "Best Jaccard: segment " << childrenOrParents[iBest] <<
                ", jaccard " << jaccardBest << endl;
        }


        // Add the best segment to the path.
        segmentId0 = childrenOrParents[iBest];
        path.push_back(segmentId0);
        if(debug) {
            cout << "Segment " << segmentId0 << " added to the path." << endl;
        }

        // If good enough, make it a new reference segment.
        if(jaccardBest > minJaccardForReference) {
            referenceSegmentId = segmentId0;
            getOrientedReadsOnSegment(referenceSegmentId, referenceOrientedReads);
            if(debug) {
                cout << "Segment " << segmentId0 << " is the new reference segment." << endl;
            }
        }


#if 0
        // See if any of these children or parents are good enough
        // to become a new reference segment.
        // For this, both unexplained fractions must be low.
        // That is, the read composition of that segment
        // would be a subset of the read composition of the reference segment.
        newReferenceCandidates.clear();
        for(uint64_t i=0; i<childrenOrParents.size(); i++) {
            if(segmentPairInformation[i].maximumUnexplainedFraction() < unexplainedFractionThresholdForReference) {
                newReferenceCandidates.push_back(i);
            }
        }
        if(debug) {
            if(newReferenceCandidates.empty()) {
                cout << "No new reference candidates." << endl;
            } else {
                cout << "New reference candidates:"<< endl;
                for(const uint64_t i: newReferenceCandidates) {
                    cout << " " << childrenOrParents[i] << " " <<
                        segmentPairInformation[i].unexplainedFraction(0) << " " <<
                        segmentPairInformation[i].unexplainedFraction(1) << " " << endl;
                }
            }
        }


        // If there is a single new reference candidate, make it the new
        // reference and add it to the path.
        if(newReferenceCandidates.size() == 1) {
            referenceSegmentId = childrenOrParents[newReferenceCandidates.front()];
            getOrientedReadsOnSegment(referenceSegmentId, referenceOrientedReads);
            path.push_back(referenceSegmentId);

            // Keep going from here.
            segmentId0 = referenceSegmentId;
            continue;
        }



        if(newReferenceCandidates.size() > 1) {
            if(debug) {
                cout << "Multiple new reference candidates found." << endl;
                return;
            }
        }



        // See if any of these children or parents are good enough
        // to be added to the path, without becoming a new reference segment.
        // For this, only the first unexplained fractions must be low.
        // That is, the read composition of that segment
        // would be a supersetset of the read composition of the reference segment.
        nextSegmentCandidates.clear();
        for(uint64_t i=0; i<childrenOrParents.size(); i++) {
            if(segmentPairInformation[i].unexplainedFraction(0) < unexplainedFractionThreshold) {
                nextSegmentCandidates.push_back(i);
            }
        }
        if(debug) {
            if(nextSegmentCandidates.empty()) {
                cout << "No new segment candidates." << endl;
            } else {
                cout << "New segment candidates:";
                for(const uint64_t i: nextSegmentCandidates) {
                    cout << " " << childrenOrParents[i] << " " <<
                        segmentPairInformation[i].unexplainedFraction(0) << " " <<
                        segmentPairInformation[i].unexplainedFraction(1) << " " << endl;
                }
            }
        }



        // If there are no next segment candidates, stop here.
        // Later we will do better.
        if(nextSegmentCandidates.empty()) {
            if(debug) {
                cout << "No next segment candidates." << endl;
            }
            return;
        }


        // If there is a single next segment candidate, add it to the path.
        if(nextSegmentCandidates.size() == 1) {
            const uint64_t nextSegmentId = childrenOrParents[nextSegmentCandidates.front()];
            path.push_back(nextSegmentId);
            segmentId0 = nextSegmentId;
            continue;
        }



        // If there are multiple next segment candidates, stop here.
        // Later we will do better,
        if(nextSegmentCandidates.size() > 1) {
            if(debug) {
                cout << "Multiple next segment candidates." << endl;
            }
            return;
        }
#endif
    }
}
