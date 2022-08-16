
// Shasta
#include "mode3.hpp"
#include "assembleMarkerGraphPath.hpp"
#include "AssembledSegment.hpp"
#include "deduplicate.hpp"
#include "findMarkerId.hpp"
#include "html.hpp"
#include "MarkerGraph.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
#include "ReadFlags.hpp"
#include "SubsetGraph.hpp"
using namespace shasta;
using namespace mode3;

// Spoa.
#include "spoa/spoa.hpp"

// Boost libraries.
// Include disjoint_sets.hpp first to avoid Boost problems.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/icl/discrete_interval.hpp>
#include <boost/icl/right_open_interval.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>

// Standard library.
#include <bitset>
#include "fstream.hpp"
#include <map>
#include <queue>
#include <set>
#include <unordered_set>

#include "MultithreadedObject.tpp"
template class MultithreadedObject<mode3::AssemblyGraph>;


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



void AssemblyGraph::computeMarkerGraphJourneys(size_t threadCount)
{
    const bool debug = true;

    createNew(markerGraphJourneys, "tmp-mode3-MarkerGraphJourneys");

    uint64_t batchSize = 1000;
    markerGraphJourneys.beginPass1(markers.size());
    setupLoadBalancing(markerGraphEdgeTable.size(), batchSize);
    runThreads(&AssemblyGraph::computeMarkerGraphJourneysPass1, threadCount);
    markerGraphJourneys.beginPass2();
    setupLoadBalancing(markerGraphEdgeTable.size(), batchSize);
    runThreads(&AssemblyGraph::computeMarkerGraphJourneysPass2, threadCount);
    markerGraphJourneys.endPass2();

    batchSize = 100;
    setupLoadBalancing(markerGraphJourneys.size(), batchSize);
    runThreads(&AssemblyGraph::sortMarkerGraphJourneys, threadCount);

    if(debug) {
        ofstream csv("MarkerGraphJourneys.csv");
        csv << "OrientedReadId,SegmentId,Position,ordinal0,Ordinal1\n";
        for(uint64_t i=0; i<markers.size(); i++) {
            const OrientedReadId orientedReadId = OrientedReadId::fromValue(ReadId(i));
            const auto markerGraphJourney = markerGraphJourneys[i];
            for(uint64_t position=0; position<markerGraphJourney.size(); position++) {
                const MarkerGraphJourneyEntry& entry = markerGraphJourney[position];
                csv << orientedReadId << ",";
                csv << entry.segmentId << ",";
                csv << entry.position << ",";
                csv << entry.ordinals[0] << ",";
                csv << entry.ordinals[1] << "\n";
            }
        }

    }
}



void AssemblyGraph::computeMarkerGraphJourneysPass1(size_t threadId)
{
    computeMarkerGraphJourneysPass12(1);
}



void AssemblyGraph::computeMarkerGraphJourneysPass2(size_t threadId)
{
    computeMarkerGraphJourneysPass12(2);
}



void AssemblyGraph::computeMarkerGraphJourneysPass12(uint64_t pass)
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
                    markerGraphJourneys.incrementCountMultithreaded(orientedReadId.getValue());
                } else {
                    MarkerGraphJourneyEntry markerGraphJourneyEntry;
                    markerGraphJourneyEntry.segmentId = segmentId;
                    markerGraphJourneyEntry.position = position;
                    markerGraphJourneyEntry.ordinals = markerInterval.ordinals;
                    markerGraphJourneys.storeMultithreaded(orientedReadId.getValue(), markerGraphJourneyEntry);
                }
            }
        }
    }
}



void AssemblyGraph::sortMarkerGraphJourneys(size_t threadId)
{
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over marker graph edges assigned to this batch.
        for(uint64_t i=begin; i!=end; ++i) {
            auto markerGraphJourney = markerGraphJourneys[i];
            sort(markerGraphJourney.begin(), markerGraphJourney.end());
        }
    }
}


// The assembly graph journey of an oriented read
// is the sequence of segmentIds it encounters.
void AssemblyGraph::computeAssemblyGraphJourneys()
{
    const bool debug = true;

    // Initialize the assembly graph journeys.
    createNew(assemblyGraphJourneys, "Mode3-AssemblyGraphJourneys");

    // Work vector defined outside the loop to reduce memory allocation overhead.
    vector<AssemblyGraphJourneyEntry> assemblyGraphJourney;

    // Loop over all oriented reads.
    for(uint64_t i=0; i<markerGraphJourneys.size(); i++) {

        // Access the marker graph journey for this oriented read.
        const span<MarkerGraphJourneyEntry> markerGraphJourney = markerGraphJourneys[i];

        // Compute the assembly graph journey.
        computeAssemblyGraphJourney(markerGraphJourney, assemblyGraphJourney);

        // Store it.
        assemblyGraphJourneys.appendVector(assemblyGraphJourney);
    }



    // Write them out.
    if(debug) {
        ofstream csv("AssemblyGraphJourneys.csv");
        for(uint64_t i=0; i<assemblyGraphJourneys.size(); i++) {
            const ReadId readId = ReadId(i >> 1);
            const Strand strand = i & 1;
            const OrientedReadId orientedReadId(readId, strand);
            const span<AssemblyGraphJourneyEntry> assemblyGraphJourney = assemblyGraphJourneys[i];

            csv << orientedReadId << ",";
            for(const AssemblyGraphJourneyEntry entry: assemblyGraphJourney) {
                csv << entry.segmentId << ",";
            }
            csv << endl;
        }
    }



    // Write them out again, with more details.
    if(debug) {
        ofstream csv("AssemblyGraphJourneysDetails.csv");
        csv << "OrientedReadId,Position,SegmentId,"
            "First position,First ordinal0,First ordinal1,"
            "Last position,Last ordinal0,Last ordinal1\n";
        for(uint64_t i=0; i<assemblyGraphJourneys.size(); i++) {
            const ReadId readId = ReadId(i >> 1);
            const Strand strand = i & 1;
            const OrientedReadId orientedReadId(readId, strand);
            const span<AssemblyGraphJourneyEntry> assemblyGraphJourney = assemblyGraphJourneys[i];

            for(uint64_t position=0; position<assemblyGraphJourney.size(); position++) {
                const AssemblyGraphJourneyEntry& entry = assemblyGraphJourney[position];
                const MarkerGraphJourneyEntry& first = entry.markerGraphJourneyEntries[0];
                const MarkerGraphJourneyEntry& last = entry.markerGraphJourneyEntries[1];
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



// Given the marker graph journey of an oriented read,
// find the corresponding assembly graph journey.
void AssemblyGraph::computeAssemblyGraphJourney(
    const span<MarkerGraphJourneyEntry> markerGraphJourney,
    vector<AssemblyGraphJourneyEntry>& assemblyGraphJourney)
{
    // Start with an empty journey.
    assemblyGraphJourney.clear();

    // Loop over the marker graph journey, looking for places
    // where the segmentId changes.
    for(uint32_t i=0; i<markerGraphJourney.size(); /* Increment later */) {
        const MarkerGraphJourneyEntry& markerGraphJourneyEntry = markerGraphJourney[i];
        const uint64_t segmentId = markerGraphJourneyEntry.segmentId;

        // Move to the end of the streak with the same segmentId.
        const uint32_t streakBegin = i;
        uint32_t streakEnd = streakBegin + 1;
        for(;
            streakEnd<markerGraphJourney.size() and
            (markerGraphJourney[streakEnd].segmentId == segmentId);
            streakEnd++) {
        }

        // Store this segmentId in the assembly graph journey.
        AssemblyGraphJourneyEntry assemblyGraphJourneyEntry;
        assemblyGraphJourneyEntry.segmentId = segmentId;
        assemblyGraphJourneyEntry.markerGraphJourneyEntries[0] = markerGraphJourney[streakBegin];
        assemblyGraphJourneyEntry.markerGraphJourneyEntries[1] = markerGraphJourney[streakEnd - 1];
        assemblyGraphJourney.push_back(assemblyGraphJourneyEntry);

        // Prepare to handle the next segment.
        i = streakEnd;
    }
}



void AssemblyGraph::computeAssemblyGraphJourneyInfos()
{
    const bool debug = true;

    const uint64_t segmentCount = paths.size();
    const uint64_t readCount = assemblyGraphJourneys.size()/2;

    createNew(assemblyGraphJourneyInfos, "Mode3-AssemblyGraphJourneyInfos");

    // Pass 1.
    assemblyGraphJourneyInfos.beginPass1(segmentCount);
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto assemblyGraphJourney = assemblyGraphJourneys[orientedReadId.getValue()];

            for(uint64_t position=0; position<assemblyGraphJourney.size(); position++) {
                const AssemblyGraphJourneyEntry& entry = assemblyGraphJourney[position];
                assemblyGraphJourneyInfos.incrementCount(entry.segmentId);
            }
        }
    }

    // Pass 2.
    assemblyGraphJourneyInfos.beginPass2();
    for(ReadId readId=0; readId<readCount; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto assemblyGraphJourney = assemblyGraphJourneys[orientedReadId.getValue()];

            for(uint64_t position=0; position<assemblyGraphJourney.size(); position++) {
                const AssemblyGraphJourneyEntry& entry = assemblyGraphJourney[position];
                assemblyGraphJourneyInfos.store(entry.segmentId, make_pair(orientedReadId, position));
            }
        }
    }
    assemblyGraphJourneyInfos.endPass2();

    // Sort.
    for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
        const auto v = assemblyGraphJourneyInfos[segmentId];
        sort(v.begin(), v.end());
    }


    if(debug) {
        ofstream csv("SegmentJourneyInfo.csv");
        csv << "SegmentId,OrientedReadId,Position in assembly graph journey\n";
        for(uint64_t segmentId=0; segmentId<segmentCount; segmentId++) {
            const auto v = assemblyGraphJourneyInfos[segmentId];
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

    for(ReadId readId=0; readId<assemblyGraphJourneys.size()/2; readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const auto journey = assemblyGraphJourneys[orientedReadId.getValue()];

            for(uint64_t i=1; i<journey.size(); i++) {
                const auto& previous = journey[i-1];
                const auto& current = journey[i];
                SHASTA_ASSERT(previous.segmentId != current.segmentId);

                const SegmentPair segmentPair = make_pair(previous.segmentId, current.segmentId);
                transitionMap[segmentPair].push_back(
                    make_pair(orientedReadId, Transition({
                    previous.markerGraphJourneyEntries[1],
                    current.markerGraphJourneyEntries[0]})));

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
            links.push_back(Link(segmentId0, segmentId1));
            transitions.appendVector(transitionVector);
        }
    }

    // Store link separation.
    for(uint64_t linkId=0; linkId<links.size(); linkId++) {
        Link& link = links[linkId];

        // Check if these two segments are adjacent in the marker graph.
        const uint64_t segmentId0 = link.segmentId0;
        const uint64_t segmentId1 = link.segmentId1;
        const auto path0 = paths[segmentId0];
        const auto path1 = paths[segmentId1];
        const MarkerGraph::Edge lastEdge0 = markerGraph.edges[path0.back()];
        const MarkerGraph::Edge firstEdge1 = markerGraph.edges[path1.front()];
        if(lastEdge0.target == firstEdge1.source) {
            // The segments are adjacent. Set the link separation to 0.
            link.segmentsAreAdjacent = true;
            link.separation = 0;
        } else {
            // The segments are not adjacent.
            // Use the transitions to estimate the separation.
            const auto linkTransitions = transitions[linkId];
            const double separation = linkSeparation(linkTransitions, path0.size());

            link.segmentsAreAdjacent = false;
            link.separation = int32_t(std::round(separation));
        }
    }



    ofstream csv("Links.csv");
    csv << "LinkId,SegmentId0,SegmentId1,Coverage,Adjacent,Separation\n";
    for(uint64_t linkId=0; linkId<links.size(); linkId++) {
        Link& link = links[linkId];

        csv << linkId << ",";
        csv << link.segmentId0 << ",";
        csv << link.segmentId1 << ",";
        csv << transitions[linkId].size() << ",";
        csv << (link.segmentsAreAdjacent ? "Yes" : "No") << ",";
        csv << link.separation << "\n";
    }

}



// Initial construction of the AssemblyGraph.
AssemblyGraph::AssemblyGraph(
    const string& largeDataFileNamePrefix,
    size_t largeDataPageSize,
    size_t threadCount,
    uint64_t readRepresentation,
    uint64_t k, // Marker length
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph) :
    MultithreadedObject<AssemblyGraph>(*this),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    largeDataPageSize(largeDataPageSize),
    readRepresentation(readRepresentation),
    k(k),
    reads(reads),
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

    // Compute marker graph and assembly graph journeys of all oriented reads.
    // We permanently store only the assembly graph journeys.
    computeMarkerGraphJourneys(threadCount);
    computeAssemblyGraphJourneys();
    markerGraphJourneys.remove();
    computeAssemblyGraphJourneyInfos();

    // Find transitions from segment to segment in the marker graph
    // journeys of all oriented reads, and store them keyed by the pair of segments.
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
    uint64_t readRepresentation,
    uint64_t k, // Marker length
    const Reads& reads,
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph) :
    MultithreadedObject<AssemblyGraph>(*this),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    readRepresentation(readRepresentation),
    k(k),
    reads(reads),
    markers(markers),
    markerGraph(markerGraph)
{
    accessExistingReadOnly(paths, "Mode3-Paths");
    accessExistingReadOnly(segmentCoverage, "Mode3-SegmentCoverage");
    accessExistingReadOnly(markerGraphEdgeTable, "Mode3-MarkerGraphEdgeTable");
    accessExistingReadOnly(assemblyGraphJourneys, "Mode3-AssemblyGraphJourneys");
    accessExistingReadOnly(assemblyGraphJourneyInfos, "Mode3-AssemblyGraphJourneyInfos");
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



void AssemblyGraph::writeGfa(const string& baseName) const
{
    ofstream gfa(baseName + ".gfa");
    ofstream csv(baseName + ".csv");

    // Write the headers.
    gfa << "H\tVN:Z:1.0\n";
    csv << "Segment,Length,Average coverage,Read count\n";

    // Write the segments.
    for(uint64_t segmentId=0; segmentId<paths.size(); segmentId++) {
        const auto path = paths[segmentId];
        gfa <<
            "S\t" << segmentId << "\t" <<
            "*\tLN:i:" << path.size() << "\n";

        csv << segmentId << ",";
        csv << path.size() << ",";
        csv << segmentCoverage[segmentId] << ",";
        csv << assemblyGraphJourneyInfos[segmentId].size() << "\n";
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



void AssemblyGraph::addClusterPairs(size_t threadId, uint64_t startSegmentId)
{
    // EXPOSE THESE CONSTANTS WHEN CODE STABILIZES.
    const uint64_t minCommonReadCount = 10;
    const double maxUnexplainedFraction = 0.25;
    const double minJaccard = 0.7;
    const uint64_t pairCountPerSegment = 1;
    const uint64_t maxDistance = 200;

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
                if(info.jaccard() < minJaccard) {
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



void AssemblyGraph::analyzeSubgraph(
    const vector<uint64_t>& segmentIds,
    vector<AnalyzeSubgraphClasses::Cluster>& clusters,
    bool debug) const
{
    if(segmentIds.size() <= 64) {
        analyzeSubgraphTemplate<64>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 128) {
        analyzeSubgraphTemplate<128>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 192) {
        analyzeSubgraphTemplate<192>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 256) {
        analyzeSubgraphTemplate<256>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 320) {
        analyzeSubgraphTemplate<320>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 384) {
        analyzeSubgraphTemplate<384>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 448) {
        analyzeSubgraphTemplate<448>(segmentIds, clusters, debug);
    } else if(segmentIds.size() <= 512) {
        analyzeSubgraphTemplate<512>(segmentIds, clusters, debug);
    } else {
        SHASTA_ASSERT(0);
    }
}



template<uint64_t N> void AssemblyGraph::analyzeSubgraphTemplate(
    const vector<uint64_t>& unsortedSegmentIds,
    vector<AnalyzeSubgraphClasses::Cluster>& clusters,
    bool debug) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const double fractionThreshold = 0.05;
    const uint64_t minClusterCoverage = 6;
    const uint64_t minSegmentCoverage = 6;

    using BitVector = std::bitset<N>;
    using JourneySnippet = AnalyzeSubgraphClasses::JourneySnippet;
    using Cluster = AnalyzeSubgraphClasses::Cluster;
    using SnippetGraphVertex = AnalyzeSubgraphClasses::SnippetGraphVertex;
    using SnippetGraph = AnalyzeSubgraphClasses::SnippetGraph;
    using vertex_descriptor = SnippetGraph::vertex_descriptor;

    // Create a sorted version of the segmentIds. We will need it later.
    vector<uint64_t> segmentIds = unsortedSegmentIds;
    sort(segmentIds.begin(), segmentIds.end());

    // Gather triplets (orientedReadId, position in assembly graph journey, segmentId).
    using Triplet = tuple<OrientedReadId, uint64_t, uint64_t>;
    vector<Triplet> triplets;
    for(const uint64_t segmentId: segmentIds) {
        const auto v = assemblyGraphJourneyInfos[segmentId];
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
    // Each streak generates a JourneySnippet.
    vector<JourneySnippet> snippets;
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
        JourneySnippet snippet;
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
        ofstream csv("JourneySnippets.csv");
        csv << "SnippetIndex,OrientedReadId,First position,LastPosition,SegmentIds\n";
        for(uint64_t snippetIndex=0; snippetIndex<snippets.size(); snippetIndex++) {
            const JourneySnippet& snippet = snippets[snippetIndex];
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
        const JourneySnippet& snippet = snippets[snippetIndex];
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

    for(const JourneySnippet& snippet: snippets) {
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
    AssemblyPath& path // The segmentId's of the path.
    ) const
{
    createAssemblyPath3(startSegmentId, direction, path);
}



// Create an assembly path starting at a given segment.
void AssemblyGraph::createAssemblyPath1(
    uint64_t startSegmentId,
    uint64_t direction,    // 0 = forward, 1 = backward
    vector<uint64_t>& path // The segmentId's of the path.
    ) const
{
    // CONSTANTS TO EXPOSE WHEN CODE STABILIZES.
    const uint64_t minimumLinkCoverage = 4;
    const double minJaccardForReference = 0.8;
    const uint64_t minCommonForReference = 6;

    const bool debug = true;
    if(debug) {
        cout << "Creating an assembly path starting at segment " << startSegmentId <<
            ", direction "<< direction << endl;
    }

    // Start with a path consisting only of startSegmentId.
    path.clear();
    std::set<uint64_t> pathSet;
    if(isBackSegment[startSegmentId]) {
        return;
    }
    path.push_back(startSegmentId);
    pathSet.insert(startSegmentId);

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
        uint64_t commonBest = 0;
        double jaccardBest = -1.;
        for(uint64_t i=0; i<childrenOrParents.size(); i++) {
            if(isBackSegment[childrenOrParents[i]]) {
                continue;
            }
            const double jaccard = segmentPairInformation[i].jaccard();
            if(jaccard > jaccardBest) {
                iBest = i;
                jaccardBest = jaccard;
                commonBest = segmentPairInformation[i].commonCount;
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
        if(pathSet.contains(segmentId0)) {
            cout << "Circular path detected at segment " << segmentId0 << endl;
            return;
        }
        pathSet.insert(segmentId0);
        if(debug) {
            cout << "Segment " << segmentId0 << " added to the path." << endl;
        }

        // If good enough, make it a new reference segment.
        if(jaccardBest > minJaccardForReference and commonBest >= minCommonForReference) {
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



// Create an assembly path starting at a given segment.
void AssemblyGraph::createAssemblyPath2(
    uint64_t startSegmentId,
    uint64_t direction,    // 0 = forward, 1 = backward
    vector<uint64_t>& path // The segmentId's of the path.
    ) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxDistance = 20000;  // Markers.
    const uint64_t minLinkCoverage = 6;
    const uint64_t minCommon = 6;
    const double minJaccard = 0.7;
    const double maxUnexplainedFraction = 0.25;
    const int32_t minLinkSeparation = -20;


    const bool debug = true;
    if(debug) {
        cout << "Creating a " <<
            (direction==0 ? "forward" : "backward") <<
            " assembly path starting at segment " << startSegmentId << endl;
    }

    // Start with a path consisting only of startSegmentId.
    path.clear();
    path.push_back(startSegmentId);

    // Keep track of the milestones reached, to avoid cycles.
    std::set<uint64_t> milestones;
    milestones.insert(startSegmentId);



    // At each iteration, we start from segmentIdA and move
    // in the specified direction until we find segmentIdB with
    // sufficiently high Jaccard similarity and number of
    // common oriented reads with segmentIdA.
    // Then we fill in a path between segmentIdA and segmentIdB.
    // When the iteration begins, segmentIdA is the last segment
    // of the path.
    uint64_t segmentIdA = startSegmentId;
    vector<uint64_t> segments;
    while(true) {
        uint64_t segmentIdB = findSimilarSegment(
            segmentIdA, direction,
            maxDistance, minLinkCoverage, minLinkSeparation, minCommon, maxUnexplainedFraction, minJaccard, segments);
        SHASTA_ASSERT(not segments.empty());
        SHASTA_ASSERT(segments.front() == segmentIdA);
        SHASTA_ASSERT((segmentIdB == invalid<uint64_t>) or (segments.back() == segmentIdB));

        if(debug) {
            cout << "segmentIdA " << segmentIdA << ", segmentIdB " << segmentIdB << endl;
            cout << "Found the following segments:";
            for(const uint64_t segmentId: segments) {
                cout << " " << segmentId;
            }
            cout << endl;
        }

        // If we did not find a segment with high Jaccard similarity,
        // the path must end here.
        if(segmentIdB == invalid<uint64_t>) {
            break;
        }


        // If this is not the first time we reach this milestone, stop here
        // to avoid cycles.
        if(milestones.contains(segmentIdB)) {
            if(debug) {
                cout << "Previously reached milestone " << segmentIdB << endl;
            }
            break;
        }

        if(debug) {
            cout << "Milestone segment " << segmentIdB << endl;
        }

        // Add the segments to the path.
        for(uint64_t i=1; i<segments.size(); i++) {
            path.push_back(segments[i]);
        }

        // Prepare for the next iteration.
        milestones.insert(segmentIdA);
        segmentIdA = segmentIdB;
    }
}



void AssemblyGraph::createAssemblyPath3(
    uint64_t startSegmentId,
    uint64_t direction,    // 0 = forward, 1 = backward
    AssemblyPath& path
    ) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCommonForLink = 3;
    const uint64_t minCommonForReference = 3;
    const double minJaccard = 0.7;
    const int32_t minLinkSeparation = -20;

    const bool debug = false;
    if(debug) {
        cout << "createAssemblyPath3 begins at segment " << startSegmentId <<
            ", direction " << direction << endl;
    }



    // At each iteration, we start from segmentIdA (the current "reference segment")
    // and move in the specified direction until we find segmentIdB with
    // sufficiently high Jaccard similarity and number of
    // common oriented reads with segmentIdA.
    // At each step, we choose the links that has the most common oriented
    // reads with the current reference segment.
    uint64_t referenceSegmentId = startSegmentId;
    SegmentOrientedReadInformation infoReference;
    getOrientedReadsOnSegment(referenceSegmentId, infoReference);
    uint64_t segmentId0 = startSegmentId;
    path.clear();
    path.segments.push_back(make_pair(startSegmentId, true));
    vector<uint64_t> lastIterationSegments;
    while(true) {

        if(debug) {
            cout << "Reference segment " << referenceSegmentId <<
                ", segmentId0 " << segmentId0 << endl;
        }

        // Loop over outgoing or incoming links of the current segment.
        // Find the link with the most common reads with the reference segment.
        const auto linkIds = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
        if(linkIds.empty()) {
            if(debug) {
                cout << "No links in this direction." << endl;
            }
            break;
        }
        uint64_t linkIdBest = invalid<uint64_t>;
        uint64_t commonOrientedReadCountBest = 0;
        for(const uint64_t linkId: linkIds) {

            // If link separation is too negative, skip it.
            // The goal here is to avoid cycles in paths.
            const Link& link = links[linkId];
            if(link.separation < minLinkSeparation) {
                continue;
            }

            // Count the number of common oriented reads between the reference segment and this link.
            uint64_t commonOrientedReadCount;
            analyzeSegmentLinkPair(referenceSegmentId, linkId, commonOrientedReadCount);

            // If better than the one we have it, record it.
            if(commonOrientedReadCount > commonOrientedReadCountBest) {
                linkIdBest = linkId;
                commonOrientedReadCountBest = commonOrientedReadCount;
            }
        }
        if(commonOrientedReadCountBest < minCommonForLink) {
            if(debug) {
                cout << "No good links found." << endl;
            }
            break;
        }
        const uint64_t linkId = linkIdBest;
        if(debug) {
            cout << "Best link " << linkId <<
                ", " << commonOrientedReadCountBest <<
                " common oriented reads with the reference segment." << endl;
        }

        // Get the segment at the other side of this link.
        const Link& link = links[linkId];
        const uint64_t segmentId1 = (direction==0) ? link.segmentId1 : link.segmentId0;
        if(debug) {
            cout << "segmentId1 " << segmentId1 << endl;
        }
        lastIterationSegments.push_back(segmentId1);

        // Check segmentId1 against the reference segment.
        SegmentOrientedReadInformation info1;
        getOrientedReadsOnSegment(segmentId1, info1);
        SegmentPairInformation info;
        analyzeSegmentPair(
            referenceSegmentId, segmentId1,
            infoReference, info1,
            markers, info);
        if(debug) {
            cout << "Jaccard " << info.jaccard() << endl;
        }

        // If the Jaccard similarity is low, this becomes the new reference segment.
        if(info.commonCount >= minCommonForReference and info.jaccard() >= minJaccard) {
            referenceSegmentId = segmentId1;
            getOrientedReadsOnSegment(referenceSegmentId, infoReference);
            if(debug) {
                cout << "New reference segment is " << segmentId1 << endl;
            }
            for(const uint64_t segmentId: lastIterationSegments) {
                path.segments.push_back(make_pair(segmentId, false));
            }
            path.segments.back().second = true;
            lastIterationSegments.clear();
        }

        segmentId0 = segmentId1;
    }



    if(debug) {
        cout << "createAssemblyPath3 ends." << endl;
    }
}



// Count the number of common oriented reads between a segment and a link,
// without counting oriented reads that appear more than once on the
// segment or on the link.
void AssemblyGraph::analyzeSegmentLinkPair(
    uint64_t segmentId,
    uint64_t linkId,
    uint64_t& commonOrientedReadCount
) const
{
    // The oriented reads in this segment,
    // with some extra information that we don't care about here.
    const auto segmentOrientedReads = assemblyGraphJourneyInfos[segmentId];

    // The oriented reads in this link,
    // with some extra information that we don't care about here.
    const auto linkOrientedReads = transitions[linkId];

    // Joint loop over oriented reads.
    commonOrientedReadCount = 0;
    const auto segmentBegin = segmentOrientedReads.begin();
    const auto segmentEnd = segmentOrientedReads.end();
    const auto linkBegin = linkOrientedReads.begin();
    const auto linkEnd = linkOrientedReads.end();
    auto itSegment = segmentBegin;
    auto itLink = linkBegin;
    while(itSegment != segmentEnd and itLink != linkEnd) {

        if(itSegment->first < itLink->first) {
            ++itSegment;
            continue;
        }
        if(itLink->first < itSegment->first) {
            ++itLink;
            continue;
        }
        SHASTA_ASSERT(itSegment->first == itLink->first);

        // If it appears more than once in the segment, skip it.
        auto itSegmentNext = itSegment + 1;
        if(itSegmentNext != segmentEnd and  itSegmentNext->first == itSegment->first) {
            ++itSegment;
            ++itLink;
            continue;
        }
        if(itSegment != segmentBegin) {
            auto itSegmentPrevious = itSegment - 1;
            if(itSegmentPrevious->first == itSegment->first) {
                ++itSegment;
                ++itLink;
                continue;
            }
        }

        // If it appears more than once in the link, skip it.
        auto itLinkNext = itLink + 1;
        if(itLinkNext != linkEnd and  itLinkNext->first == itLink->first) {
            ++itSegment;
            ++itLink;
            continue;
        }
        if(itLink != linkBegin) {
             auto itLinkPrevious = itLink - 1;
            if(itLinkPrevious->first == itLink->first) {
                ++itSegment;
                ++itLink;
                continue;
            }
        }

        // Ok, this is a common oriented read that appears only once
        // in both the segment and the link.
        ++commonOrientedReadCount;
        ++itSegment;
        ++itLink;
    }

}



// Given a segment, use a BFS to move in the specified direction until
// we find a segment with sufficiently high Jaccard similarity
// and number of common reads.
// This returns invalid<uint64_t> if no such segment is found
// within the specified distance.
uint64_t AssemblyGraph::findSimilarSegmentBfs(
    uint64_t segmentIdA,
    uint64_t direction, // 0 = forward, 1 = backward
    uint64_t maxDistance,
    uint64_t minCommon,
    double minJaccard) const
{
    const bool debug = true;
    if(debug) {
        cout << "findSimilarSegmentBfs starts " << segmentIdA << " " << direction << endl;
    }

    // Sanity check.
    SHASTA_ASSERT(maxDistance > 0);

    // Get the oriented reads on segmentIdA.
    SegmentOrientedReadInformation infoA;
    getOrientedReadsOnSegment(segmentIdA, infoA);

    // Initialize a BFS.
    std::queue<uint64_t> q;
    q.push(segmentIdA);

    // Keep track of segments we already encountered and their distance.
    // Key = segmentId;
    // Value = distance.
    std::map<uint64_t, uint64_t> distanceMap;
    distanceMap.insert(make_pair(segmentIdA, 0));



    // BFS loop.
    while(not q.empty()) {

        // Dequeue a segment.
        const uint64_t segmentId0 = q.front();
        q.pop();
        const uint64_t distance0 = distanceMap[segmentId0];
        SHASTA_ASSERT(distance0 < maxDistance);
        const uint64_t distance1 = distance0 + 1;
        if(debug) {
            cout << "dequeued " << segmentId0 << " " << distance0 << endl;
        }

        // Loop over outgoing or incoming links.
        const auto linkIds = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
        for(const uint64_t linkId: linkIds) {
            const Link& link = links[linkId];

            // Get the segment at the other side of this link.
            const uint64_t segmentId1 = (direction==0) ? link.segmentId1 : link.segmentId0;

            // If we already found it, skip it.
            if(distanceMap.contains(segmentId1)) {
                continue;
            }

            if(debug) {
                cout << "found " << segmentId1 << " " << distance1 << endl;
            }

            // Get the oriented reads on segmentId1.
            SegmentOrientedReadInformation info1;
            getOrientedReadsOnSegment(segmentId1, info1);

            // See how similar this is to segmentIdA.
            SegmentPairInformation infoA1;
            analyzeSegmentPair(
                segmentIdA, segmentId1,
                infoA, info1,
                markers, infoA1);

            // If this satisfies our criteria, we are done.
            if(infoA1.commonCount >= minCommon and
                infoA1.jaccard() >= minJaccard) {
                if(debug) {
                    cout << "findSimilarSegmentBFS returns " << segmentId1 << " " << direction << endl;
                }
                return segmentId1;
            }

            // This segment did not satisfy our criteria, so we
            // have to continue the BFS.
            if(distance1 < maxDistance) {
                q.push(segmentId1);
                distanceMap.insert(make_pair(segmentId1, distance1));
                if(debug) {
                    cout << "enqueued " << segmentId1 << " " << distance1 << endl;
                }
            }
        }

    }



    // If getting here, we did not find a segment that satisfies
    // the requested criteria.
    if(debug) {
        cout << "findSimilarSegmentBfs returns invalid" << endl;
    }
    return invalid<uint64_t>;
}



// Given a segment, move in the specified direction,
// in order of increasing distance in markers, until
// we find a segment with sufficiently high Jaccard similarity
// and number of common reads.
// This returns invalid<uint64_t> if no such segment is found
// within the specified distance.
uint64_t AssemblyGraph::findSimilarSegment(
    uint64_t segmentIdA,
    uint64_t direction,     // 0 = forward, 1 = backward
    uint64_t maxDistance,   // In markers
    uint64_t minLinkCoverage,
    int32_t minLinkSeparation,
    uint64_t minCommon,
    double maxUnexplainedFraction,
    double minJaccard,
    vector<uint64_t>& segments) const
{
    const bool debug = false;
    if(debug) {
        cout << "findSimilarSegment begins, segmentIdA " << segmentIdA << endl;
    }
    // Sanity check.
    SHASTA_ASSERT(maxDistance > 0);

    segments.clear();

    // Get the oriented reads on segmentIdA.
    SegmentOrientedReadInformation infoA;
    getOrientedReadsOnSegment(segmentIdA, infoA);

    // (Offset, segmentId) for queued segments.
    std::multimap<uint64_t, uint64_t> q;
    q.insert(make_pair(0, segmentIdA));

    // The segments that we already encountered.
    std::set<uint64_t> visitedSegmentSet;
    visitedSegmentSet.insert(segmentIdA);



    // Search loop.
    while(not q.empty()) {

        // Dequeue the segment with the smallest offset.
        const auto it0 = q.begin();
        const uint64_t segmentId0 = it0->second;
        q.erase(it0);

        // Analyze against segmentIdA.
        SegmentOrientedReadInformation info0;
        getOrientedReadsOnSegment(segmentId0, info0);
        SegmentPairInformation infoA0;
        analyzeSegmentPair(
            segmentIdA, segmentId0,
            infoA, info0,
            markers, infoA0);

        // Add it to our list of segments, if possible.
        const double unexplainedFraction = infoA0.unexplainedFraction(0);
        if(unexplainedFraction < maxUnexplainedFraction) {
            segments.push_back(segmentId0);
        }

        // If unexplained fraction and Jaccard similarity are low, we are done.
        if(segmentId0 != segmentIdA) {
            if((unexplainedFraction < maxUnexplainedFraction) and (infoA0.jaccard() >= minJaccard)) {
                SHASTA_ASSERT(segments.back() == segmentId0);
                return segmentId0;
            }
        }

        if(debug) {
            cout << "Dequeued " << segmentId0 << endl;
        }

        // Loop over outgoing or incoming links.
        const auto linkIds = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
        for(const uint64_t linkId: linkIds) {
            // If link coverage is too low, skip.
            if(transitions.size(linkId) < minLinkCoverage) {
                continue;
            }

            // If link separation is too negative, skip it.
            // The goal here is to avoid cycles in paths.
            const Link& link = links[linkId];
            if(link.separation < minLinkSeparation) {
                continue;
            }

            // Get the segment at the other side of this link.
            const uint64_t segmentId1 = (direction==0) ? link.segmentId1 : link.segmentId0;
            if(debug) {
                cout << "Found " << segmentId1 << endl;
            }

            // If we already found it, skip it.
            if(visitedSegmentSet.contains(segmentId1)) {
                if(debug) {
                    cout << "Already found, skipping." << endl;
                }
                continue;
            }
            visitedSegmentSet.insert(segmentId1);

            // Get the oriented reads on segmentId1.
            SegmentOrientedReadInformation info1;
            getOrientedReadsOnSegment(segmentId1, info1);

            // Analyze similarity to segmentIdA.
            SegmentPairInformation infoA1;
            analyzeSegmentPair(
                segmentIdA, segmentId1,
                infoA, info1,
                markers, infoA1);

            // If not enough common segments, skip it.
            if(infoA1.commonCount < minCommon) {
                if(debug) {
                    cout << "Not enough common reads." << endl;
                }
                continue;
            }

            // Offset estimates are not reliable.
            // Don't use them to rule out segments.
#if 0
            // If not in the expected direction, skip it.
            uint64_t offset;
            if(direction == 0) {
                if(infoA1.offset < 0) {
                    if(debug) {
                        cout << "Not in the forward direction." << endl;
                    }
                    continue;
                } else {
                    offset = uint64_t(infoA1.offset);
                }
            } else {
                if(infoA1.offset > 0) {
                    if(debug) {
                        cout << "Not in the backward direction." << endl;
                    }
                    continue;
                } else {
                    offset = uint64_t(-infoA1.offset);
                }
            }
#endif

            // If we went too far, skip it.
            if(labs(infoA1.offset) > maxDistance) {
                if(debug) {
                    cout << "Too far." << endl;
                }
                continue;
            }

            // Queue it.
            q.insert(make_pair(labs(infoA1.offset), segmentId1));
            if(debug) {
                cout << "Queued " << segmentId1 << endl;
            }

        }

    }



    // If getting here, we did not find a segment that satisfies
    // the requested criteria.
    return invalid<uint64_t>;
}



// BFS with given begin/end.
// Does a BFS which starts at segmentIdA.
// and ends when segmentIdB is encountered.
// The BFS if forward if direction is 0
// and backward if direction is 1.
// Computes a vector of all the segments encountered,
// excluding segmentIdA and segmentIdB,
// in the order in which they are encountered in the BFS.
void AssemblyGraph::targetedBfs(
    uint64_t segmentIdA,
    uint64_t segmentIdB,
    uint64_t direction,
    vector<uint64_t>& segments
    ) const
{

    // Initialize the BFS.
    std::queue<uint64_t> q;
    q.push(segmentIdA);

    // Keep track of segments we already encountered.
    std::set<uint64_t> segmentSet;
    segmentSet.insert(segmentIdA);



    // BFS loop.
    segments.clear();
    while(not q.empty()) {

        // Dequeue a segment.
        const uint64_t segmentId0 = q.front();
        q.pop();

        // Loop over outgoing or incoming links.
        const auto linkIds = (direction == 0) ? linksBySource[segmentId0] : linksByTarget[segmentId0];
        for(const uint64_t linkId: linkIds) {
            const Link& link = links[linkId];

            // Get the segment at the other side of this link.
            const uint64_t segmentId1 = (direction==0) ? link.segmentId1 : link.segmentId0;

            // If we found segmentIdB, we are done.
            if(segmentId1 == segmentIdB) {
                break;
            }

            // If we already found it, skip it.
            if(segmentSet.contains(segmentId1)) {
                continue;
            }

            // Queue and store this segment.
            q.push(segmentId1);
            segments.push_back(segmentId1);
            segmentSet.insert(segmentId1);
        }
    }

}



// Assemble sequence for an AssemblyPath.
void AssemblyGraph::assemblePathSequence(const AssemblyPath& assemblyPath) const
{
    const bool debug = true;
    ofstream html;
    if(debug) {
        cout << timestamp << "AssemblyGraph::assemblePathSequence begins for a path of length " <<
            assemblyPath.segments.size() << endl;
        html.open("Msa.html");
        writeHtmlBegin(html, "Mode 3 assembly");
    }

    // The list of oriented reads that can be used to assemble links.
    // These are the oriented reads that appear in either the previous or next
    // reference segment in the path.
    // This list changes every time we encounter a new reference segment.
    // It is kept sorted for efficiency.
    vector<OrientedReadId> orientedReadsForLinks;


    // Assemble each segment on the path.
    vector<AssembledSegment> assembledSegments(assemblyPath.segments.size());
    for(uint64_t i=0; i<assemblyPath.segments.size(); i++) {
        const uint64_t segmentId = assemblyPath.segments[i].first;
        AssembledSegment& assembledSegment = assembledSegments[i];
        assembleMarkerGraphPath(readRepresentation, k,
            markers, markerGraph, paths[segmentId], false, assembledSegment);
    }



    // Write out assembled sequence of each segment.
    if(debug) {
        ofstream fasta("PathSegmentsSequence.fasta");
        ofstream txt("PathSegmentsRleSequence.txt");
        for(uint64_t i=0; i<assemblyPath.segments.size(); i++) {
            const uint64_t segmentId = assemblyPath.segments[i].first;
            AssembledSegment& assembledSegment = assembledSegments[i];

            fasta <<
                ">" << i <<
                " segment " << segmentId <<
                " length " << assembledSegment.rawSequence.size() << "\n";
            copy(assembledSegment.rawSequence.begin(), assembledSegment.rawSequence.end(),
                ostream_iterator<Base>(fasta));
            fasta << "\n";

            txt << i << " " << segmentId << endl;
            copy(assembledSegment.runLengthSequence.begin(), assembledSegment.runLengthSequence.end(),
                ostream_iterator<Base>(txt));
            txt << "\n";
            for(const uint32_t r: assembledSegment.repeatCounts) {
                if(r < 10) {
                    txt << r;
                } else {
                    txt << "*";
                }
            }
            txt << "\n";
        }
    }



    // Assemble segment and link sequence into path sequence.
    string sequence;
    for(uint64_t i=0; ; i++) {
        const auto& p = assemblyPath.segments[i];
        const uint64_t segmentId0 = p.first;
        if(debug) {
            cout << "Working on segment " << segmentId0 << " at position " << i <<
                " in assembly path." << endl;
        }
        const bool isReferenceSegment = p.second;
        AssembledSegment& assembledSegment0 = assembledSegments[i];

        // The first segment in the path must be a reference segment.
        if(i == 0) {
            SHASTA_ASSERT(isReferenceSegment);
        }



        // Add the sequence of this segment.
        for(const Base b: assembledSegment0.rawSequence) {
            sequence.push_back(b.character());
        }

        // If this is the last segment, we are done.
        if(i == assemblyPath.segments.size() - 1) {
            break;
        }



        // If this is a reference segment, update
        // the list of oriented reads that can be used to assemble links.
        // These are the oriented reads that appear in either the previous or next
        // reference segment in the path.
        if(isReferenceSegment) {
            orientedReadsForLinks.clear();

            // Add the oriented reads in segmentId0, our new
            // reference segment.
            for(const auto& p: assemblyGraphJourneyInfos[segmentId0]) {
                orientedReadsForLinks.push_back(p.first);
            }

            // Look for the next reference segment in the path.
            uint64_t nextReferenceSegmentId = invalid<uint64_t>;
            for(uint64_t j=i+1; j<assemblyPath.segments.size(); j++) {
                const auto& q = assemblyPath.segments[i];
                if(q.second) {
                    nextReferenceSegmentId = q.first;
                    break;
                }
            }
            SHASTA_ASSERT(nextReferenceSegmentId != invalid<uint64_t>);
            for(const auto& p: assemblyGraphJourneyInfos[nextReferenceSegmentId]) {
                orientedReadsForLinks.push_back(p.first);
            }

            // Deduplicate and sort.
            deduplicate(orientedReadsForLinks);
        }

        // This is not the last segment. Access the next segment in the path.
        const uint64_t segmentId1 = assemblyPath.segments[i+1].first;
        AssembledSegment& assembledSegment1 = assembledSegments[i+1];



        // If the segments are consecutive in the marker graph,
        // just take out the last k RLE bases from the sequence.
        // We will get them from the next segment.
        const auto path0 = paths[segmentId0];
        const auto path1 = paths[segmentId1];
        const MarkerGraph::Edge lastEdge0 = markerGraph.edges[path0.back()];
        const MarkerGraph::Edge firstEdge1 = markerGraph.edges[path1.front()];
        if(lastEdge0.target == firstEdge1.source) {

            // Compute the number of bases to be removed.
            uint64_t removeCount = 0;
            for(uint64_t j=assembledSegment0.repeatCounts.size() - k; j<assembledSegment0.repeatCounts.size(); j++) {
                removeCount += assembledSegment0.repeatCounts[j];
            }

            // Remove them from the sequence we already have.
            sequence.resize(sequence.size() - removeCount);

            // We are done.
            continue;
        }



        // If getting here, this segment and the next are not adjacent in the marker graph.
        // We need to assembly the link between them.
        // We can only use oriented reads that are in the link
        // and also in orientedReadsForLinks.

        // Locate the link between segmentId0 and segmentId1.
        uint64_t linkId01 = invalid<uint64_t>;
        for(uint64_t linkId: linksBySource[segmentId0]) {
            if(links[linkId].segmentId1 == segmentId1) {
                linkId01 = linkId;
                break;
            }
        }
        SHASTA_ASSERT(linkId01!= invalid<uint64_t>);

        // Access the oriented reads that transition over this link.
        const auto transitions01 = transitions[linkId01];

        if(debug) {
            cout << "Assembling link " << linkId01 << " " << segmentId0 << "->" << segmentId1 <<
                " with coverage " << transitions01.size() <<
                " and " << orientedReadsForLinks.size() <<
                " reference oriented reads." << endl;
            cout << segmentId0 << " path length " << path0.size() << endl;
            cout << segmentId1 << " path length " << path1.size() << endl;
        }



        // Consider all oriented reads that are both in transitions01 and in
        // orientedReadsForLinks. These are the oriented reads that we can
        // actually use to assemble the link.



        // First, find:
        // - The position in segmentId0 of the leftmost transition.
        // - The position in segmentId1 of the rightmost transition.
        uint64_t minEdgePosition0 = path0.size();
        uint64_t maxEdgePosition1 = 0;
        for(const auto& p: transitions[linkId01]) {
            const OrientedReadId orientedReadId = p.first;

            // If not in orientedReadsForLinks, skip it.
            if(not std::binary_search(
                orientedReadsForLinks.begin(), orientedReadsForLinks.end(),
                orientedReadId)) {
                continue;
            }

            // Access the transition from segmentId0 to segmentId1 for this oriented read.
            const Transition& transition = p.second;

            minEdgePosition0 = min(minEdgePosition0, uint64_t(transition[0].position));
            maxEdgePosition1 = max(maxEdgePosition1, uint64_t(transition[1].position));
        }
        if(debug) {
            cout << "minEdgePosition0 " << minEdgePosition0 << endl;
            cout << "maxEdgePosition1 " << maxEdgePosition1 << endl;
        }

        // When getting here:
        // - minPosition0 is the leftmost position of the transitions in path0.
        // - minPosition1 is the rightmost position of the transitions in path1.
        // These positions are edge positions in path0 and path1.
        // We will do a multiple sequence alignment of the oriented reads,
        // using the sequence of segmentId0 to extend to the left all reads to minPosition0,
        // and using the sequence of segmentId1 to extend to the right all reads to minPosition1,

        // Get the corresponding vertex positions in segmentId0 and segmentId1.
        const uint64_t minVertexPosition0 = minEdgePosition0 + 1;
        const uint64_t maxVertexPosition1 = maxEdgePosition1;
        if(debug) {
            cout << "minVertexPosition0 " << minVertexPosition0 << endl;
            cout << "maxVertexPosition1 " << maxVertexPosition1 << endl;
        }



        // Now extract the portion of each oriented read sequence that
        // will be used to assemble this link.
        uint64_t actualLinkCoverage = 0;
        vector< vector<Base> > orientedReadsSequencesForAssembly;
        vector< vector<uint64_t> > orientedReadsRepeatCountsForAssembly;
        for(const auto& p: transitions[linkId01]) {
            const OrientedReadId orientedReadId = p.first;

            // If not in orientedReadsForLinks, skip it.
            if(not std::binary_search(
                orientedReadsForLinks.begin(), orientedReadsForLinks.end(),
                orientedReadId)) {
                continue;
            }

            // Access the transition from segmentId0 to segmentId1 for this oriented read.
            const Transition& transition = p.second;
            ++actualLinkCoverage;

            // Get the ordinals of the last appearance of this oriented
            // read on segmentId0 and the first on segmentId1,
            // and the corresponding markers.
            const uint32_t ordinal0 = transition[0].ordinals[1];
            const uint32_t ordinal1 = transition[1].ordinals[0];
            if(debug) {
                cout << orientedReadId << " " << ordinal0 << " " << ordinal1 << endl;
            }
            const CompressedMarker& marker0 = markers[orientedReadId.getValue()][ordinal0];
            const CompressedMarker& marker1 = markers[orientedReadId.getValue()][ordinal1];

            // Get the positions of these markers on the oriented read.
            // If using RLE, these are RLE positions.
            const uint32_t position0 = marker0.position;
            const uint32_t position1 = marker1.position;

            // Extract the sequence between these markers (including the markers).
            vector<Base> orientedReadSequence;
            vector<uint8_t> orientedReadRepeatCounts;
            if(readRepresentation == 1) {
                // RLE.
                for(uint64_t position=position0; position<position1+k; position++) {
                    Base b;
                    uint8_t r;
                    tie(b, r) = reads.getOrientedReadBaseAndRepeatCount(orientedReadId, uint32_t(position));
                    orientedReadSequence.push_back(b);
                    orientedReadRepeatCounts.push_back(r);
                }
            } else {
                // Raw sequence.
                for(uint64_t position=position0; position<position1+k; position++) {
                    const Base b = reads.getOrientedReadBase(orientedReadId, uint32_t(position));
                    orientedReadSequence.push_back(b);
                    orientedReadRepeatCounts.push_back(uint8_t(1));
                }
            }

            if(debug) {
                copy(orientedReadSequence.begin(), orientedReadSequence.end(), ostream_iterator<Base>(cout));
                cout << endl;
                for(uint8_t r: orientedReadRepeatCounts) {
                    if(r < 10) {
                        cout << int(r);
                    } else {
                        cout << "*";
                    }
                }
                cout << endl;
            }

            // We need to extend the sequence of this read to the left,
            // using segmentId0 sequence, up to minVertexPosition0,
            // so the portions of all reads we will be using for the MSA
            // all begin in the same place.
            vector<Base> leftSequence;
            vector<uint32_t> leftRepeatCounts;
            const uint64_t vertexPosition0 = transition[0].position + 1;  // Add 1 to get vertex position.
            const uint64_t begin0 = assembledSegment0.vertexOffsets[minVertexPosition0];
            const uint64_t end0 = assembledSegment0.vertexOffsets[vertexPosition0];
            for(uint64_t position=begin0; position!=end0; position++) {
                leftSequence.push_back(assembledSegment0.runLengthSequence[position]);
                leftRepeatCounts.push_back(assembledSegment0.repeatCounts[position]);
            }
            if(debug) {
                cout << "Left extend sequence:" << endl;
                copy(leftSequence.begin(), leftSequence.end(), ostream_iterator<Base>(cout));
                cout << endl;
                for(uint32_t r: leftRepeatCounts) {
                    if(r < 10) {
                        cout << int(r);
                    } else {
                        cout << "*";
                    }
                }
                cout << endl;
            }

            vector<Base> rightSequence;
            vector<uint32_t> rightRepeatCounts;
            const uint64_t vertexPosition1 = transition[1].position;
            const uint64_t begin1 = assembledSegment1.vertexOffsets[vertexPosition1] + k;
            const uint64_t end1 = assembledSegment1.vertexOffsets[maxVertexPosition1] + k;
            // cout << "*** begin1 " << begin1 << endl;
            // cout << "*** end1 " << end1 << endl;
            for(uint64_t position=begin1; position!=end1; position++) {
                // cout << "*** " << position << " " << assembledSegment1.runLengthSequence[position] <<
                //     " " << int(assembledSegment1.repeatCounts[position]) << endl;
                rightSequence.push_back(assembledSegment1.runLengthSequence[position]);
                rightRepeatCounts.push_back(assembledSegment1.repeatCounts[position]);
            }
            if(debug) {
                cout << "Right extend sequence:" << endl;
                copy(rightSequence.begin(), rightSequence.end(), ostream_iterator<Base>(cout));
                cout << endl;
                for(uint32_t r: rightRepeatCounts) {
                    if(r < 10) {
                        cout << int(r);
                    } else {
                        cout << "*";
                    }
                }
                cout << endl;
            }

            // Construct the extended sequence for this oriented read,
            // to be used in the MSA.
            vector<Base> orientedReadExtendedSequence;
            vector<uint64_t> orientedReadExtendedRepeatCounts;
            copy(leftSequence.begin(), leftSequence.end(), back_inserter(orientedReadExtendedSequence));
            copy(leftRepeatCounts.begin(), leftRepeatCounts.end(), back_inserter(orientedReadExtendedRepeatCounts));
            copy(orientedReadSequence.begin(), orientedReadSequence.end(), back_inserter(orientedReadExtendedSequence));
            copy(orientedReadRepeatCounts.begin(), orientedReadRepeatCounts.end(), back_inserter(orientedReadExtendedRepeatCounts));
            copy(rightSequence.begin(), rightSequence.end(), back_inserter(orientedReadExtendedSequence));
            copy(rightRepeatCounts.begin(), rightRepeatCounts.end(), back_inserter(orientedReadExtendedRepeatCounts));

            orientedReadsSequencesForAssembly.push_back(orientedReadExtendedSequence);
            orientedReadsRepeatCountsForAssembly.push_back(orientedReadExtendedRepeatCounts);

            if(debug) {
                cout << "Sequence for this oriented read to be used for MSA:" << endl;
                copy(orientedReadExtendedSequence.begin(), orientedReadExtendedSequence.end(), ostream_iterator<Base>(cout));
                cout << endl;
                for(uint64_t r: orientedReadExtendedRepeatCounts) {
                    if(r < 10) {
                        cout << int(r);
                    } else {
                        cout << "*";
                    }
                }
                cout << endl;
            }
        }
        if(debug) {
            cout << "Link assembly will use " << actualLinkCoverage <<
                " oriented reads." << endl;
        }

        // Compute the consensus sequence for the link.
        vector<Base> consensusRleSequence;
        vector<uint64_t> consensusRepeatCounts;
        if(debug) {
            html<< "<h2>Link " << linkId01 << "</h2>\n";
        }
        computeLinkConsensusUsingSpoa(
            orientedReadsSequencesForAssembly,
            orientedReadsRepeatCountsForAssembly,
            readRepresentation,
            debug,
            html,
            consensusRleSequence,
            consensusRepeatCounts
            );
    }



    // Write out the path sequence.
    if(debug) {
        ofstream fasta("PathSequence.fasta");
        fasta << ">path\n" << sequence << "\n";
    }

    if(debug) {
        cout << timestamp << "AssemblyGraph::assemblePathSequence ends." << endl;
        writeHtmlEnd(html);
    }
}



// Compute consensus sequence for Link, given sequences of
// the oriented reads, which must all be anchored on both sides.
void AssemblyGraph::computeLinkConsensusUsingSpoa(
    const vector< vector<Base> > rleSequences,
    const vector< vector<uint64_t> > repeatCounts,
    uint64_t readRepresentation,
    bool debug,
    ostream& html,
    vector<Base>& consensusRleSequence,
    vector<uint64_t>& consensusRepeatCounts
    ) const
{
    // Create the spoa alignment engine and elignment graph.
    const spoa::AlignmentType alignmentType = spoa::AlignmentType::kNW;
    const int8_t match = 1;
    const int8_t mismatch = -1;
    const int8_t gap = -1;
    auto spoaAlignmentEngine = spoa::AlignmentEngine::Create(alignmentType, match, mismatch, gap);
    spoa::Graph spoaAlignmentGraph;

    // Add the oriented read sequences to the alignment.
    string sequenceString;
    for(const vector<Base>& sequence: rleSequences) {

        // Add it to the alignment.
        sequenceString.clear();
        for(const Base base: sequence) {
            sequenceString += base.character();
        }
        auto alignment = spoaAlignmentEngine->Align(sequenceString, spoaAlignmentGraph);
        spoaAlignmentGraph.AddAlignment(alignment, sequenceString);
    }

    // Compute the multiple sequence alignment.
    const vector<string> msa = spoaAlignmentGraph.GenerateMultipleSequenceAlignment();
    const string consensus = spoaAlignmentGraph.GenerateConsensus();
    const uint64_t msaLength = msa.front().size();
    if(debug) {
        cout << "Multiple sequence alignment has length " << msaLength << ":" << endl;
        for(const string& s: msa) {
            cout << s << endl;
        }
    }



    // Compute coverage for each base at each position of the MSA.
    // Use position 4 for gaps.
    vector < array < uint64_t, 5 > > coverage(msaLength, {0, 0, 0, 0, 0});
    for(const string& s: msa) {
        SHASTA_ASSERT(s.size() == msaLength);
        for(uint64_t i=0; i<msaLength; i++) {
            const char c = s[i];
            switch(c) {
            case 'A': ++coverage[i][0]; break;
            case 'C': ++coverage[i][1]; break;
            case 'G': ++coverage[i][2]; break;
            case 'T': ++coverage[i][3]; break;
            case '-': ++coverage[i][4]; break;
            default: SHASTA_ASSERT(0);
            }
        }
    }



    // Html output of the alignment.
    if(html.good()) {
        html << rleSequences.size() << " oriented reads<br>\n";
        html << "Alignment length " << msaLength << "<br>\n";
        html << "<div style='font-family:monospace;white-space:nowrap;'>\n";
        for(const string& s: msa) {
            for(const char c: s) {
                switch(c) {
                case 'A':
                    // html << "<span style='background-color:hsl(0,100%,70%)'>A</span>";
                    html << "<span style='background-color:#ff6666'>A</span>";
                    break;
                case 'C':
                    // html << "<span style='background-color:hsl(240,100%,70%)'>C</span>";
                    html << "<span style='background-color:#6666ff'>C</span>";
                    break;
                case 'G':
                    // html << "<span style='background-color:hsl(60,100%,70%)'>G</span>";
                    html << "<span style='background-color:#ffff66'>G</span>";
                    break;
                case 'T':
                    // html << "<span style='background-color:hsl(120,100%,70%)'>T</span>";
                    html << "<span style='background-color:#66ff66'>T</span>";
                    break;
                case '-':
                    html << "-";
                    break;
                default:
                    SHASTA_ASSERT(0);
                }
            }
            html << "<br>\n";
        }

        // Coverage.
        for(uint64_t k=0; k<=4; k++) {
            for(uint64_t i=0; i<msaLength; i++) {
                const uint64_t n = coverage[i][k];
                char c;
                if(n == 0) {
                    c = '.';
                } else if(n < 10) {
                    c = '0' + char(n);
                } else if(n < 100) {
                    c = 'A' + char(n/10);
                } else {
                    c = '*';
                }
                html << c;
            }
            html << " ";

            switch(k) {
            case 0: html << "A"; break;
            case 1: html << "C"; break;
            case 2: html << "G"; break;
            case 3: html << "T"; break;
            case 4: html << "-"; break;
            }
            html << "<br>\n";
        }
        html << "</div>\n";
#if 0
        for(uint64_t i=0; i<msaLength; i++) {
            html << i << " ";
            for(uint64_t k=0; k<=4; k++) {
                html << coverage[i][k] << " ";
            }
            html << std::accumulate(coverage[i].begin(), coverage[i].end(), 0);
            html << "<br>\n";
        }
#endif
    }
}
