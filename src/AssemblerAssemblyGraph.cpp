// Shasta.
#include "Assembler.hpp"
#include "LocalAssemblyGraph.hpp"
#include "orderPairs.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "chrono.hpp"
#include "iterator.hpp"
#include <numeric>
#include <queue>
#include <unordered_map>



// In the assembly graph, each edge corresponds to a linear chain
// of edges the marker graph.
// This code finds all linear chains of edges in the marker graph,
// and generates an assembly graph edge for each chain it finds.
// This creates:
// - assemblyGraph.edgeLists.
// - assemblyGraph.markerToAssemblyTable
void Assembler::createAssemblyGraphEdges()
{
    // Some shorthands.
    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::EdgeId;

    // Check that we have what we need.
    checkMarkerGraphVerticesAreAvailable();
    checkMarkerGraphEdgesIsOpen();
    const auto& edges = markerGraph.edges;

    // Flag to control debug output.
    const bool debug = false;

    // Vector used to keep track of marker graph edges that were already found.
    const EdgeId edgeCount = markerGraph.edges.size();
    MemoryMapped::Vector<bool> wasFound;
    wasFound.createNew(
        largeDataName("tmp-createAssemblyGraphVertices-wasFound"),
        largeDataPageSize);
    wasFound.resize(edgeCount);
    fill(wasFound.begin(), wasFound.end(), false);

    // Initialize the edge lists.
    assemblyGraph.edgeLists.createNew(
        largeDataName("AssemblyGraphEdgeLists"),
        largeDataPageSize);



    // Work vectors reused for each chain.
    vector<EdgeId> nextEdges;
    vector<EdgeId> previousEdges;
    vector<EdgeId> chain;



    // Main loop over all edges of the marker graph.
    // At each iteration we find a new linear chain of edges.
    for(EdgeId startEdgeId=0; startEdgeId<edgeCount; startEdgeId++) {
        const auto& startEdge = edges[startEdgeId];

        if(debug) {
            cout << "Working on start edge " << startEdgeId;
            cout << " " << startEdge.source << "->" << startEdge.target << endl;
        }

        // If this edge is not part of cleaned up marker graph, skip it.
        if(startEdge.wasRemoved()) {
            if(debug) {
                cout << "Start edge was removed." << endl;
            }
            continue;
        }

        // If we already found this edge, skip it.
        // It is part of a chain we already found.
        if(wasFound[startEdgeId]) {
            if(debug) {
                cout << "This edge is part of a chain we already found." << endl;
            }
            continue;
        }

        // Follow the chain forward.
        EdgeId edgeId = startEdgeId;
        bool isCircularChain = false;
        while(true) {
            edgeId = nextEdgeInMarkerGraphPrunedStrongSubgraphChain(edgeId);
            if(edgeId == invalidGlobalMarkerGraphEdgeId) {
                break;
            }
            if(edgeId == startEdgeId) {
                isCircularChain = true;
                break;
            }
            nextEdges.push_back(edgeId);
        }

        // Follow the chain backward.
        if(!isCircularChain) {
            edgeId = startEdgeId;
            while(true) {
                edgeId = previousEdgeInMarkerGraphPrunedStrongSubgraphChain(edgeId);
                if(edgeId == invalidGlobalMarkerGraphEdgeId) {
                    break;
                }
                previousEdges.push_back(edgeId);
            }
        }

        // Gather the chain.
        copy(previousEdges.rbegin(), previousEdges.rend(), back_inserter< vector<EdgeId> >(chain));
        chain.push_back(startEdgeId);
        copy(nextEdges.begin(), nextEdges.end(), back_inserter< vector<EdgeId> >(chain));

        // Mark all the edges in the chain as found.
        for(const EdgeId edgeId: chain) {
            wasFound[edgeId] = true;
        }

        // Store this chain as a new vertex of the assembly graph.
        assemblyGraph.edgeLists.appendVector(chain);

        // Write out the chain.
        if(debug) {
            cout << "Chain: ";
            cout << edges[chain.front()].source << " ";
            for(const EdgeId edgeId: chain) {
                cout << edges[edgeId].target << " ";
            }
            cout << endl;
        }

        // Cleanup.
        nextEdges.clear();
        previousEdges.clear();
        chain.clear();
    }



    // Check that only and all edges of the cleaned up marker graph
    // were found.
    for(EdgeId edgeId=0; edgeId<edgeCount; edgeId++) {
        const auto& edge = markerGraph.edges[edgeId];
        if(edge.wasRemoved()) {
            CZI_ASSERT(!wasFound[edgeId]);
        } else {
            CZI_ASSERT(wasFound[edgeId]);
        }
    }

    wasFound.remove();

    // cout << "The assembly graph has " << assemblyGraph.edgeLists.size() << " edges." << endl;



    // Create the markerToAssemblyTable.
    assemblyGraph.markerToAssemblyTable.createNew(
        largeDataName("MarkerToAssemblyTable"),
        largeDataPageSize);
    assemblyGraph.markerToAssemblyTable.resize(edges.size());
    fill(
        assemblyGraph.markerToAssemblyTable.begin(),
        assemblyGraph.markerToAssemblyTable.end(),
        make_pair(std::numeric_limits<VertexId>::max(), 0));
    for(EdgeId edgeId=0; edgeId<assemblyGraph.edgeLists.size(); edgeId++) {
        const MemoryAsContainer<EdgeId> chain = assemblyGraph.edgeLists[edgeId];
        for(uint32_t position=0; position!=chain.size(); position++) {
            const EdgeId edgeId = chain[position];
            assemblyGraph.markerToAssemblyTable[edgeId] = make_pair(edgeId, position);
        }
    }



    // Create a histogram of size (chain length) of assembly graph vertices.
    vector<size_t> histogram;
    for(EdgeId edgeId=0; edgeId<assemblyGraph.edgeLists.size(); edgeId++) {
        const size_t size = assemblyGraph.edgeLists.size(edgeId);
        if(histogram.size() <= size) {
            histogram.resize(size+1);
        }
        ++(histogram[size]);
    }
    ofstream csv("AssemblyGraphChainLengthHistogram.csv");
    csv << "ChainLength, Frequency\n";
    for(size_t size=0; size<histogram.size(); size++) {
        const size_t frequency = histogram[size];
        if(frequency) {
            csv << size << "," << frequency << "\n";
        }
    }
}



void Assembler::accessAssemblyGraphVertices()
{
    assemblyGraph.vertices.accessExistingReadOnly(
        largeDataName("AssemblyGraphVertices"));
    assemblyGraph.markerToAssemblyTable.accessExistingReadOnly(
        largeDataName("MarkerToAssemblyTable"));
}



// This uses assemblyGraph.edgeLists to create:
// - assemblyGraph.vertices
// - assemblyGraph.edges
// - assemblyGraph.edgesBySource
// - assemblyGraph.edgesByTarget

// In this function:
// mgv indicates a marker graph vertex id.
// mge indicates a marker graph edge id.
// agv indicates an assembly graph vertex id.
// age indicates an assembly graph edge id.
void Assembler::createAssemblyGraphVertices()
{

    // cout << timestamp << "Creating assembly graph vertices." << endl;

    // Check that we have what we need.
    CZI_ASSERT(assemblyGraph.edgeLists.isOpen());
    checkMarkerGraphEdgesIsOpen();

    // Shorthands for vertex and edge ids.
    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::EdgeId;



    // Create assemblyGraph.vertices.
    // Each marker graph vertex that is the first or last vertex
    // of a linear edge chain corresponding to an edge of the assembly graph.
    // generates an assembly graph vertex.
    assemblyGraph.vertices.createNew(
        largeDataName("AssemblyGraphVertices"),
        largeDataPageSize);
    for(EdgeId age=0; age<assemblyGraph.edgeLists.size(); age++) {
        const auto chain = assemblyGraph.edgeLists[age];
        CZI_ASSERT(chain.size() > 0);
        const EdgeId firstChainEdgeId = *(chain.begin());
        const EdgeId lastChainEdgeId = *(chain.end() - 1);
        const MarkerGraph::Edge& firstChainEdge = markerGraph.edges[firstChainEdgeId];
        const MarkerGraph::Edge& lastChainEdge = markerGraph.edges[lastChainEdgeId];
        const VertexId mgv0 = firstChainEdge.source;
        const VertexId mgv1 = lastChainEdge.target;
        assemblyGraph.vertices.push_back(mgv0);
        assemblyGraph.vertices.push_back(mgv1);
    }

    // Deduplicate.
    sort(assemblyGraph.vertices.begin(), assemblyGraph.vertices.end());
    assemblyGraph.vertices.resize(
        std::unique(assemblyGraph.vertices.begin(), assemblyGraph.vertices.end()) -
        assemblyGraph.vertices.begin());
    assemblyGraph.vertices.unreserve();
    // cout << "The assembly graph has " << assemblyGraph.vertices.size() <<
    //    " vertices and " << assemblyGraph.edgeLists.size() << " edges." << endl;

    // Create a map that gives the assembly graph vertex corresponding to a
    // marker graph vertex.
    std::unordered_map<VertexId, VertexId> vertexMap;
    for(VertexId agv=0; agv<assemblyGraph.vertices.size(); agv++) {
        const VertexId mgv = assemblyGraph.vertices[agv];
        vertexMap.insert(make_pair(mgv, agv));
    }



    // Create assemblyGraph edges.
    assemblyGraph.edges.createNew(
        largeDataName("AssemblyGraphEdges"),
        largeDataPageSize);
    assemblyGraph.edges.resize(assemblyGraph.edgeLists.size());
    for(EdgeId age=0; age<assemblyGraph.edgeLists.size(); age++) {
        const auto chain = assemblyGraph.edgeLists[age];
        CZI_ASSERT(chain.size() > 0);
        const EdgeId firstChainEdgeId = *(chain.begin());
        const EdgeId lastChainEdgeId = *(chain.end() - 1);
        const MarkerGraph::Edge& firstChainEdge = markerGraph.edges[firstChainEdgeId];
        const MarkerGraph::Edge& lastChainEdge = markerGraph.edges[lastChainEdgeId];
        const VertexId mgv0 = firstChainEdge.source;
        const VertexId mgv1 = lastChainEdge.target;
        const auto it0 = vertexMap.find(mgv0);
        const auto it1 = vertexMap.find(mgv1);
        CZI_ASSERT(it0 != vertexMap.end());
        CZI_ASSERT(it1 != vertexMap.end());
        const VertexId agv0 = it0->second;
        const VertexId agv1 = it1->second;
        AssemblyGraph::Edge& assemblyGraphEdge = assemblyGraph.edges[age];
        assemblyGraphEdge.source = agv0;
        assemblyGraphEdge.target = agv1;

        // Compute and store average coverage along the edges of this chain.
        size_t sum = 0;
        for(GlobalMarkerGraphEdgeId markerGraphEdgeId: chain) {
            const MarkerGraph::Edge& markerGraphEdge = markerGraph.edges[markerGraphEdgeId];
            sum += markerGraphEdge.coverage;
        }
        assemblyGraphEdge.averageCoverage = uint32_t(sum / chain.size());
    }



    // cout << timestamp << "Creating assembly graph edges by source and by target." << endl;

    assemblyGraph.edgesBySource.createNew(
        largeDataName("AssemblyGraphEdgesBySource"),
        largeDataPageSize);
    assemblyGraph.edgesByTarget.createNew(
        largeDataName("AssemblyGraphEdgesByTarget"),
        largeDataPageSize);
    assemblyGraph.edgesBySource.beginPass1(assemblyGraph.vertices.size());
    assemblyGraph.edgesByTarget.beginPass1(assemblyGraph.vertices.size());
    for(const AssemblyGraph::Edge& edge: assemblyGraph.edges) {
        assemblyGraph.edgesBySource.incrementCount(edge.source);
        assemblyGraph.edgesByTarget.incrementCount(edge.target);
    }
    assemblyGraph.edgesBySource.beginPass2();
    assemblyGraph.edgesByTarget.beginPass2();
    for(EdgeId edgeId=0; edgeId<assemblyGraph.edges.size(); edgeId++) {
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[edgeId];
        assemblyGraph.edgesBySource.store(edge.source, edgeId);
        assemblyGraph.edgesByTarget.store(edge.target, edgeId);
    }
    assemblyGraph.edgesBySource.endPass2();
    assemblyGraph.edgesByTarget.endPass2();

    // cout << timestamp << "Done creating assembly graph vertices." << endl;
}



void Assembler::accessAssemblyGraphEdgeLists()
{
    assemblyGraph.edgeLists.accessExistingReadOnly(
        largeDataName("AssemblyGraphEdgeLists"));
}



void Assembler::accessAssemblyGraphEdges()
{
    assemblyGraph.edges.accessExistingReadOnly(
        largeDataName("AssemblyGraphEdges"));
    assemblyGraph.edgesBySource.accessExistingReadOnly(
        largeDataName("AssemblyGraphEdgesBySource"));
    assemblyGraph.edgesByTarget.accessExistingReadOnly(
        largeDataName("AssemblyGraphEdgesByTarget"));
}

void Assembler::writeAssemblyGraph(const string& fileName) const
{
    assemblyGraph.writeGraphviz(fileName);
}


// Assemble sequence for all edges of the assembly graph.
void Assembler::assemble(
    size_t threadCount,

    uint32_t markerGraphEdgeLengthThresholdForConsensus,

    // Parameter to control whether we use spoa of marginPhase
    // to compute consensus sequence for marker graph edges.
    bool useMarginPhase
    )
{

    // Check that we have what we need.
    checkKmersAreOpen();
    checkReadsAreOpen();
    checkMarkersAreOpen();
    checkMarkerGraphVerticesAreAvailable();
    checkMarkerGraphEdgesIsOpen();
    CZI_ASSERT(assemblyGraph.edgeLists.isOpen());
    assembleData.markerGraphEdgeLengthThresholdForConsensus = markerGraphEdgeLengthThresholdForConsensus;
    if(useMarginPhase) {
        checkMarginPhaseWasSetup();
    }
    assembleData.useMarginPhase = useMarginPhase;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCount << " threads." << endl;

    // Allocate data structures to store assembly results for each thread.
    assembleData.allocate(threadCount);

    // Do all the assemblies.
    cout << "Assembly begins for " << assemblyGraph.edgeLists.size() <<
        " edges of the assembly graph." << endl;
    setupLoadBalancing(assemblyGraph.edgeLists.size(), 1);
    runThreads(&Assembler::assembleThreadFunction, threadCount,
        "threadLogs/assemble");

    // Find the pair(thread, index in thread) that the assembly for each edge is stored in.
    const auto uninitializedPair = make_pair(
        std::numeric_limits<size_t>::max(),
        std::numeric_limits<size_t>::max());
    vector< pair<size_t, size_t> > edgeTable(assemblyGraph.edgeLists.size(), uninitializedPair);
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        vector<AssemblyGraph::EdgeId>& edges = assembleData.edges[threadId];
        for(size_t i=0; i<edges.size(); i++) {
            edgeTable[edges[i]] = make_pair(threadId, i);
        }
    }


    // Store the assembly results found by each thread.
    assemblyGraph.sequences.createNew(largeDataName("AssembledSequences"), largeDataPageSize);
    assemblyGraph.repeatCounts.createNew(largeDataName("AssembledRepeatCounts"), largeDataPageSize);
    for(AssemblyGraph::EdgeId edgeId=0; edgeId<assemblyGraph.edgeLists.size(); edgeId++) {
        const auto& p = edgeTable[edgeId];
        CZI_ASSERT(p != uninitializedPair);
        const size_t threadId = p.first;
        const size_t i = p.second;

        // Store the sequence for this edge.
        LongBaseSequences& threadSequences =
            *(assembleData.sequences[threadId]);
        const auto vertexSequence = threadSequences[i];
        assemblyGraph.sequences.append(vertexSequence);

        // Store the repeat count for this edge.
        MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& threadRepeatCounts =
            *(assembleData.repeatCounts[threadId]);
        const MemoryAsContainer<uint8_t> vertexRepeatCounts = threadRepeatCounts[i];
        assemblyGraph.repeatCounts.appendVector();
        for(uint8_t r: vertexRepeatCounts) {
            assemblyGraph.repeatCounts.append(r);
        }
    }


    // Clean up the results stores by each thread.
    assembleData.free();

    // Compute the total number of bases assembled.
    size_t totalBaseCount = 0;
    for(size_t i=0; i<assemblyGraph.repeatCounts.totalSize(); i++) {
        totalBaseCount += assemblyGraph.repeatCounts.begin()[i];
    }
    cout << timestamp << "Assembled a total " << totalBaseCount <<
        " bases for " << assemblyGraph.edgeLists.size() << " assembly graph edges." << endl;
}



void Assembler::assembleThreadFunction(size_t threadId)
{
    ostream& out = getLog(threadId);
    const uint32_t markerGraphEdgeLengthThresholdForConsensus = assembleData.markerGraphEdgeLengthThresholdForConsensus;
    const bool useMarginPhase = assembleData.useMarginPhase;

    // Initialize data structures for this thread.
    vector<AssemblyGraph::EdgeId>& edges = assembleData.edges[threadId];

    assembleData.sequences[threadId] = make_shared<LongBaseSequences>();
    LongBaseSequences& sequences = *(assembleData.sequences[threadId]);
    sequences.createNew(largeDataName("tmp-Sequences-" + to_string(threadId)), largeDataPageSize);

    assembleData.repeatCounts[threadId] = make_shared<MemoryMapped::VectorOfVectors<uint8_t, uint64_t> >();
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& repeatCounts = *(assembleData.repeatCounts[threadId]);
    repeatCounts.createNew(largeDataName("tmp-RepeatCounts-" + to_string(threadId)), largeDataPageSize);

    vector<Base> edgeSequence;
    vector<uint32_t> edgeRepeatCounts;

    // Loop over batches allocated to this thread.
    size_t begin, end;
    while(getNextBatch(begin, end)) {
        for(AssemblyGraph::EdgeId edgeId=begin; edgeId!=end; edgeId++) {
            out << timestamp << "Assembling assembly graph edge " << edgeId <<
                " corresponding to " <<
                assemblyGraph.edgeLists[edgeId].size() <<
                " marker graph edges." << endl;
            assembleAssemblyGraphEdge(edgeId, markerGraphEdgeLengthThresholdForConsensus, useMarginPhase, edgeSequence, edgeRepeatCounts);

            // Store the edge id.
            edges.push_back(edgeId);

            // Store the sequence.
            sequences.append(edgeSequence);

            // Store the repeat counts.
            repeatCounts.appendVector();
            for(const uint32_t r: edgeRepeatCounts) {
                repeatCounts.append(uint8_t(min(uint32_t(255), r)));
            }
        }
    }

}



void Assembler::AssembleData::allocate(size_t threadCount)
{
    edges.resize(threadCount);
    sequences.resize(threadCount);
    repeatCounts.resize(threadCount);
}



void Assembler::AssembleData::free()
{
    edges.clear();
    for(auto& sequence: sequences) {
        sequence->remove();
    }
    sequences.clear();
    for(auto& r: repeatCounts) {
        r->remove();
    }
    repeatCounts.clear();
}



void Assembler::accessAssemblyGraphSequences()
{
    assemblyGraph.sequences.accessExistingReadOnly(
        largeDataName("AssembledSequences"));
    assemblyGraph.repeatCounts.accessExistingReadOnly(
        largeDataName("AssembledRepeatCounts"));
}



void Assembler::computeAssemblyStatistics()
{

    using EdgeId = AssemblyGraph::EdgeId;

    // Check that we have what we need.
    CZI_ASSERT(assemblyGraph.vertices.isOpen);
    const size_t vertexCount = assemblyGraph.vertices.size();
    CZI_ASSERT(assemblyGraph.edges.isOpen);
    const size_t edgeCount = assemblyGraph.edges.size();
    CZI_ASSERT(assemblyGraph.sequences.isOpen());
    CZI_ASSERT(assemblyGraph.sequences.size() == edgeCount);
    CZI_ASSERT(assemblyGraph.repeatCounts.isOpen());
    CZI_ASSERT(assemblyGraph.repeatCounts.size() == edgeCount);

    // Compute raw sequence length of each edge.
    vector< pair<EdgeId, size_t> > edgeTable(edgeCount);
    size_t totalLength = 0;
    for(EdgeId edgeId=0; edgeId<edgeCount; edgeId++) {
        const MemoryAsContainer<uint8_t> repeatCounts = assemblyGraph.repeatCounts[edgeId];
        size_t length = 0;
        for(uint8_t repeatCount: repeatCounts) {
            length += repeatCount;
        }
        edgeTable[edgeId] = make_pair(edgeId, length);
        totalLength += length;
    }

    cout << "The assembly graph has " << vertexCount << endl;
    cout << "vertices and " << edgeCount << " edges." << endl;
    cout << "Total length of assembled sequence is " << totalLength << endl;

    // Sort by decreasing length.
    sort(edgeTable.begin(), edgeTable.end(), OrderPairsBySecondOnlyGreater<EdgeId, size_t>());

    // Write a csv file.
    ofstream csv("AssemblySummary.csv");
    csv << "Rank,EdgeId,Length,CumulativeLength,LengthFraction,CumulativeFraction\n";
    size_t cumulativeLength = 0;
    bool n50MessageWritten = false;
    for(size_t rank=0; rank<edgeTable.size(); rank++) {
        const pair<EdgeId, size_t>& p = edgeTable[rank];
        const EdgeId edgeId = p.first;
        const size_t length = p.second;
        cumulativeLength += length;
        csv << rank << ",";
        csv << edgeId << ",";
        csv << length << ",";
        csv << cumulativeLength << ",";
        csv << double(length) / double(totalLength) << ",";
        csv << double(cumulativeLength) /double(totalLength) << "\n";
        if(!n50MessageWritten && cumulativeLength >= totalLength/2) {
            cout << "N50 for assembly segments is " << length << endl;
            n50MessageWritten = true;
        }
    }

}



// Write the assembly graph in GFA 1.0 format defined here:
// https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md
void Assembler::writeGfa1(const string& fileName)
{
    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::EdgeId;

    ofstream gfa(fileName);

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";

    // Write a segment record for each edge.
    for(EdgeId edgeId=0; edgeId<assemblyGraph.sequences.size(); edgeId++) {
        const auto sequence = assemblyGraph.sequences[edgeId];
        const auto repeatCounts = assemblyGraph.repeatCounts[edgeId];
        CZI_ASSERT(sequence.baseCount == repeatCounts.size());
        gfa << "S\t" << edgeId << "\t";
        for(size_t i=0; i<sequence.baseCount; i++) {
            const Base b = sequence[i];
            const uint8_t repeatCount = repeatCounts[i];
            for(size_t k=0; k<repeatCount; k++) {
                gfa << b;
            }
        }
        gfa << "\n";
    }



    // Write GFA links.
    // For each vertex in the assembly graph there is a link for
    // each combination of in-edges and out-edges.
    // Therefore each assembly graph vertex generates a number of
    // links equal to the product of its in-degree and out-degree.
    const size_t k = assemblerInfo->k;
    string cigarString;
    for(VertexId vertexId=0; vertexId<assemblyGraph.vertices.size(); vertexId++) {

        // In-edges.
        const MemoryAsContainer<EdgeId> edges0 = assemblyGraph.edgesByTarget[vertexId];

        // Out-edges.
        const MemoryAsContainer<EdgeId> edges1 = assemblyGraph.edgesBySource[vertexId];

        // Loop over combinations of in-edges and out-edges.
        for(const EdgeId edge0: edges0) {
            const MemoryAsContainer<uint8_t> repeatCounts0 = assemblyGraph.repeatCounts[edge0];
            for(const EdgeId edge1: edges1) {
                const MemoryAsContainer<uint8_t> repeatCounts1 = assemblyGraph.repeatCounts[edge1];

                // Locate the last k repeat counts of v0 and the first k of v1.
                const MemoryAsContainer<uint8_t> lastRepeatCounts0(
                    repeatCounts0.begin() + repeatCounts0.size() - k,
                    repeatCounts0.end());
                const MemoryAsContainer<uint8_t> firstRepeatCounts1(
                    repeatCounts1.begin(),
                    repeatCounts1.begin() + k);


                // Construct the cigar string.
                constructCigarString(lastRepeatCounts0, firstRepeatCounts1, cigarString);

                // Write out the link record for this edge.
                gfa << "L\t" <<
                    edge0 << "\t" <<
                    "+\t" <<
                    edge1 << "\t" <<
                    "+\t" <<
                    cigarString << "\n";
            }
        }
    }
}



// Write assembled sequences in FASTA format.
void Assembler::writeFasta(const string& fileName)
{
    using EdgeId = AssemblyGraph::EdgeId;

    ofstream fasta(fileName);

    // Write a sequence for each edge of the assembly graph.
    for(EdgeId edgeId=0; edgeId<assemblyGraph.sequences.size(); edgeId++) {
        const auto sequence = assemblyGraph.sequences[edgeId];
        const auto repeatCounts = assemblyGraph.repeatCounts[edgeId];
        CZI_ASSERT(sequence.baseCount == repeatCounts.size());

        // Compute the length so we can write it in the header.
        size_t length = 0;
        for(const uint8_t repeatCount: repeatCounts) {
            length += repeatCount;
        }

        fasta << ">" << edgeId << " length " << length << "\n";
        for(size_t i=0; i<sequence.baseCount; i++) {
            const Base b = sequence[i];
            const uint8_t repeatCount = repeatCounts[i];
            for(size_t k=0; k<repeatCount; k++) {
                fasta << b;
            }
        }
        fasta << "\n";
    }

}



// Construct the CIGAR string given two vectors of repeat counts.
// Used by writeGfa1.
void Assembler::constructCigarString(
    const MemoryAsContainer<uint8_t>& repeatCounts0,
    const MemoryAsContainer<uint8_t>& repeatCounts1,
    string& cigarString
    )
{
    // Check that the repeat counts have the same length.
    const size_t k = repeatCounts0.size();
    CZI_ASSERT(repeatCounts1.size() == k);

    if(std::equal(repeatCounts0.begin(), repeatCounts0.end(), repeatCounts1.begin())) {

        // Fast path when the repeat counts are identical.
        // This is the most common case.
        // The cigar string is all matches, for a number of bases
        // equal to the total base count (in raw representation).
        size_t totalBaseCount = 0;
        for(const uint8_t r: repeatCounts0) {
            totalBaseCount += r;
        }
        cigarString = to_string(totalBaseCount) + "M";

    } else {

        // General case.
        vector< pair<char, int> > cigar;
        for(size_t i=0; i<k; i++) {
            const uint8_t repeatCount0 = repeatCounts0[i];
            const uint8_t repeatCount1 = repeatCounts1[i];

            // Matching bases.
            const uint8_t matchingCount = min(repeatCount0, repeatCount1);
            if(matchingCount) {
                if(!cigar.empty() && cigar.back().first=='M') {
                    cigar.back().second += matchingCount;
                } else {
                    cigar.push_back(make_pair('M', matchingCount));
                }
            }

            // Inserted bases (bases present in repeatCounts1 but not in repeatCounts0).
            if(repeatCount1 > repeatCount0) {
                const uint8_t insertedCount = uint8_t(repeatCount1 - repeatCount0);
                if(!cigar.empty() && cigar.back().first=='I') {
                    cigar.back().second += insertedCount;
                } else {
                    cigar.push_back(make_pair('I', insertedCount));
                }
            }

            // Deleted bases (bases present in repeatCounts0 but not in repeatCounts1).
            if(repeatCount0 > repeatCount1) {
                const uint8_t deletedCount = uint8_t(repeatCount0 - repeatCount1);
                if(!cigar.empty() && cigar.back().first=='D') {
                    cigar.back().second += deletedCount;
                } else {
                    cigar.push_back(make_pair('D', deletedCount));
                }
            }
        }
    }

}



// Extract a local assembly graph from the global assembly graph.
// This returns false if the timeout was exceeded.
bool Assembler::extractLocalAssemblyGraph(
    AssemblyGraph::EdgeId startEdgeId,
    int distance,
    double timeout,
    LocalAssemblyGraph& graph) const
{
    using vertex_descriptor = LocalAssemblyGraph::vertex_descriptor;
    using edge_descriptor = LocalAssemblyGraph::edge_descriptor;
    using VertexId = AssemblyGraph::VertexId;
    using EdgeId = AssemblyGraph::EdgeId;

    const bool debug = false;
    if(debug) {
        cout << "Begin extractLocalAssemblyGraph for edge "
            << startEdgeId << " distance " << distance << endl;
    }


    const auto startTime = steady_clock::now();

    // Add the start vertices.
    const AssemblyGraph::Edge startEdge = assemblyGraph.edges[startEdgeId];
    if(startEdge.source == startEdge.target) {
        throw runtime_error("The local assembly graph starts at a self-edge. Display for this situation is not implemented.");
    }
    const vertex_descriptor vStart0 = graph.addVertex(
        startEdge.source,
        assemblyGraph.vertices[startEdge.source],
        0);
    const vertex_descriptor vStart1 = graph.addVertex(
        startEdge.target,
        assemblyGraph.vertices[startEdge.target],
        0);


    // Do the BFS.
    std::queue<vertex_descriptor> q;
    if(distance > 0) {
        q.push(vStart0);
        q.push(vStart1);
    }
    std::set<vertex_descriptor> processedVertices;
    while(!q.empty()) {

        // See if we exceeded the timeout.
        if(timeout>0. && seconds(steady_clock::now() - startTime) > timeout) {
            graph.clear();
            return false;
        }

        // Dequeue a vertex.
        const vertex_descriptor v0 = q.front();
        q.pop();
        const LocalAssemblyGraphVertex& vertex0 = graph[v0];
        const VertexId assemblyGraphVertexId0 = vertex0.assemblyGraphVertexId;
        const int distance0 = vertex0.distance;
        const int distance1 = distance0 + 1;

        if(debug) {
            cout << "Dequeued " << assemblyGraphVertexId0 << " at distance " << distance0 << endl;
        }



        // Loop over children.
        const auto childEdges = assemblyGraph.edgesBySource[assemblyGraphVertexId0];
        for(const EdgeId edgeId: childEdges) {
            const AssemblyGraph::Edge& globalEdge = assemblyGraph.edges[edgeId];
            const VertexId assemblyGraphVertexId1 = globalEdge.target;

            if(debug) {
                cout << "Found child " << assemblyGraphVertexId1 << endl;
            }

            // Find the vertex corresponding to this child, creating it if necessary.
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(assemblyGraphVertexId1);
            if(!vertexExists) {
                v1 = graph.addVertex(
                    assemblyGraphVertexId1,
                    assemblyGraph.vertices[assemblyGraphVertexId1],
                    distance1);
                if(distance1 < distance) {
                    q.push(v1);
                }
                if(debug) {
                    cout << "Vertex added " << assemblyGraphVertexId1 << endl;
                }
            }

            // Create the edge v0->v1.
            if(processedVertices.find(v1) == processedVertices.end()) {
                edge_descriptor e;
                bool edgeExists;
                tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
                CZI_ASSERT(edgeExists);
                graph[e].edgeId = edgeId;
                if(debug) {
                    cout << "Edge added " << assemblyGraphVertexId0 << "->" << assemblyGraphVertexId1 << endl;
                }
            }
        }



        // Loop over parents.
        const auto parentEdges = assemblyGraph.edgesByTarget[assemblyGraphVertexId0];
        for(const EdgeId edgeId: parentEdges) {
            const AssemblyGraph::Edge& globalEdge = assemblyGraph.edges[edgeId];
            const VertexId assemblyGraphVertexId1 = globalEdge.source;

            if(debug) {
                cout << "Found parent " << assemblyGraphVertexId1 << endl;
            }

            // Find the vertex corresponding to this child, creating it if necessary.
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(assemblyGraphVertexId1);
            if(!vertexExists) {
                v1 = graph.addVertex(
                    assemblyGraphVertexId1,
                    assemblyGraph.vertices[assemblyGraphVertexId1],
                    distance1);
                if(distance1 < distance) {
                    q.push(v1);
                }
                if(debug) {
                    cout << "Vertex added " << assemblyGraphVertexId1 << endl;
                }
            }

            // Create the edge v1->v0.
            if(processedVertices.find(v1) == processedVertices.end()) {
                edge_descriptor e;
                bool edgeExists;
                    tie(e, edgeExists) = boost::add_edge(v1, v0, graph);
                    CZI_ASSERT(edgeExists);
                    graph[e].edgeId = edgeId;
                    if(debug) {
                        cout << "Edge added " << assemblyGraphVertexId1 << "->"
                            << assemblyGraphVertexId0 << endl;
                }
            }
        }

        processedVertices.insert(v0);
    }


#if 0
    // The BFS process did not create edges between vertices at maximum distance.
    // Do it now.
    // Loop over all vertices at maximum distance.
    BGL_FORALL_VERTICES(v0, graph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex0 = graph[v0];
        if(vertex0.distance != distance) {
            continue;
        }
        const AssemblyGraph::VertexId vertexId0 = vertex0.vertexId;

        // Loop over the children that exist in the local assembly graph
        // and are also at maximum distance.
        const auto childEdges = assemblyGraph.edgesBySource[vertexId0];
        for(uint64_t edgeId: childEdges) {
            const auto& edge = assemblyGraph.edges[edgeId];

            const AssemblyGraph::VertexId vertexId1 = edge.target;
            CZI_ASSERT(edge.source == vertexId0);
            CZI_ASSERT(vertexId1 < assemblyGraph.vertices.size());

            // See if we have a vertex for this global vertex id.
            bool vertexExists;
            vertex_descriptor v1;
            tie(vertexExists, v1) = graph.findVertex(vertexId1);

            // If it does not exist in the local marker graph, skip.
            if(!vertexExists) {
                continue;
            }

            // If it is not at maximum distance, skip.
            const LocalAssemblyGraphVertex& vertex1 = graph[v1];
            if(vertex1.distance != distance) {
                continue;
            }

            // Add the edge.
            edge_descriptor e;
            bool edgeExists;
            tie(e, edgeExists) = boost::add_edge(v0, v1, graph);
            CZI_ASSERT(edgeExists);

        }
    }
#endif

    if(debug) {
        cout << "Vertices:" << endl;
        BGL_FORALL_VERTICES(v, graph, LocalAssemblyGraph) {
            cout << graph[v].assemblyGraphVertexId << endl;
        }
        cout << "Edges:" << endl;
        BGL_FORALL_EDGES(e, graph, LocalAssemblyGraph) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            cout << graph[v0].assemblyGraphVertexId << "->";
            cout << graph[v1].assemblyGraphVertexId << endl;
        }

    }


    return true;
}



// Assemble sequence for an edge of the assembly graph.
// Optionally outputs detailed assembly information
// in html (skipped if the html pointer is 0).
void Assembler::assembleAssemblyGraphEdge(
    AssemblyGraph::EdgeId edgeId,
    uint32_t markerGraphEdgeLengthThresholdForConsensus,
    bool useMarginPhase,
    vector<Base>& assembledRunLengthSequence,
    vector<uint32_t>& assembledRepeatCounts,
    ostream* htmlPointer)
{
    const auto k = assemblerInfo->k;

    // The edges of this chain in the marker graph.
    const MemoryAsContainer<GlobalMarkerGraphEdgeId> edgeIds = assemblyGraph.edgeLists[edgeId];
    const size_t edgeCount = edgeIds.size();
    const size_t vertexCount = edgeCount + 1;

    // The vertices of this chain in the marker graph.
    vector<GlobalMarkerGraphVertexId> vertexIds;
    vertexIds.reserve(vertexCount);
    for(const GlobalMarkerGraphEdgeId edgeId: edgeIds) {
        const MarkerGraph::Edge& edge =
            markerGraph.edges[edgeId];
        vertexIds.push_back(edge.source);
    }
    const MarkerGraph::Edge& lastEdge =
        markerGraph.edges[edgeIds[edgeIds.size()-1]];
    vertexIds.push_back(lastEdge.target);

    // Vertex coverage.
    vector<uint32_t> vertexCoverage(vertexCount);
    for(size_t i=0; i<vertexCount; i++) {
        vertexCoverage[i] = uint32_t(globalMarkerGraphVertices.size(vertexIds[i]));
    }

    // Edge coverage.
    vector<uint32_t> edgeCoverage(edgeCount);
    for(size_t i=0; i<edgeCount; i++) {
        edgeCoverage[i] = uint32_t(markerGraph.edgeMarkerIntervals.size(edgeIds[i]));
    }




    // Compute consensus sequence for the vertices of the chain.
    vector< vector<Base> > vertexSequences(vertexCount);
    vector< vector<uint32_t> > vertexRepeatCounts(vertexCount);
    for(size_t i=0; i<vertexCount; i++) {
        computeMarkerGraphVertexConsensusSequence(
            vertexIds[i], vertexSequences[i], vertexRepeatCounts[i]);
    }

    // Compute consensus sequence for the edges of the chain.
    vector< vector<Base> > edgeSequences(edgeCount);
    vector< vector<uint32_t> > edgeRepeatCounts(edgeCount);
    for(size_t i=0; i<edgeCount; i++) {
        /*
        if((i%1000) == 0) {
            cout << timestamp << i << "/" << edgeCount << endl;
        }
        */
        if(useMarginPhase) {
            computeMarkerGraphEdgeConsensusSequenceUsingMarginPhase(
                edgeIds[i], edgeSequences[i], edgeRepeatCounts[i]);
        } else {
            computeMarkerGraphEdgeConsensusSequenceUsingSpoa(
                edgeIds[i], markerGraphEdgeLengthThresholdForConsensus, edgeSequences[i], edgeRepeatCounts[i]);
        }
    }


    // Compute offsets for edges and vertices.
    // A vertex/edge offset is the position of the first base
    // of the vertex/edge consensus sequence (run-length)
    // relative to the first base of assembled sequence (run-length).
    vector<uint32_t> vertexOffsets(vertexCount);
    vector<uint32_t> edgeOffsets(edgeCount);
    vertexOffsets[0] = 0;
    for(size_t i=0; i<edgeCount; i++) {
        edgeOffsets[i] = vertexOffsets[i];
        vertexOffsets[i+1] = edgeOffsets[i] + uint32_t(edgeSequences[i].size()) - uint32_t(k);
    }



    // Compute, for each vertex, the portion of vertex sequence that contributes
    // to the assembly. This is the portion that does not overlap a vertex with greater coverage.
    // (Break ties using vertex ids).
    // An edge always contributes to the assembly all but its first and last k bases
    // (this can mean no bases, for edges between overlapping vertices).
    vector< pair<uint32_t, uint32_t> > vertexAssembledPortion(vertexCount);
    for(int i=0; i<int(vertexCount); i++) {

        // Check previous vertices.
        vertexAssembledPortion[i].first = 0;
        for(int j=i-1; j>=0; j--) {
            if(vertexOffsets[j]+k < vertexOffsets[i]) {
                break;
            }
            if(vertexCoverage[j]>vertexCoverage[i] ||
                (vertexCoverage[j]==vertexCoverage[i] && vertexIds[j]<vertexIds[i])) {
                vertexAssembledPortion[i].first = vertexOffsets[j] + uint32_t(k) - vertexOffsets[i];
                break;
            }
        }

        // Check following vertices.
        vertexAssembledPortion[i].second = uint32_t(k);
        for(int j=i+1; j<int(vertexCount); j++) {
            if(vertexOffsets[i]+k < vertexOffsets[j]) {
                break;
            }
            if(vertexCoverage[j]>vertexCoverage[i] ||
                (vertexCoverage[j]==vertexCoverage[i] && vertexIds[j]<vertexIds[i])) {
                vertexAssembledPortion[i].second = vertexOffsets[j] - vertexOffsets[i];
                break;
            }
        }

        // Handle the case of a vertex that contributes nothing.
        if(vertexAssembledPortion[i].second <= vertexAssembledPortion[i].first) {
            vertexAssembledPortion[i].first = 0;
            vertexAssembledPortion[i].second = 0;
        }
    }



    // Assemble run-length sequence and raw sequence.
    // Keep track of the range each vertex and edge contributes.
    assembledRunLengthSequence.clear();
    assembledRepeatCounts.clear();
    vector<Base> assembledRawSequence;
    vector< pair<uint32_t, uint32_t> > vertexRunLengthRange(vertexCount);
    vector< pair<uint32_t, uint32_t> > vertexRawRange(vertexCount);
    vector< pair<uint32_t, uint32_t> > edgeRunLengthRange(edgeCount);
    vector< pair<uint32_t, uint32_t> > edgeRawRange(edgeCount);
    for(size_t i=0; ; i++) {

        // Vertex.
        vertexRunLengthRange[i].first = uint32_t(assembledRunLengthSequence.size());
        vertexRawRange[i].first = uint32_t(assembledRawSequence.size());
        for(uint32_t j=vertexAssembledPortion[i].first; j!=vertexAssembledPortion[i].second; j++) {
            const Base base = vertexSequences[i][j];
            const uint32_t repeatCount = vertexRepeatCounts[i][j];
            assembledRunLengthSequence.push_back(base);
            assembledRepeatCounts.push_back(repeatCount);
            for(uint32_t k=0; k!=repeatCount; k++) {
                assembledRawSequence.push_back(base);
            }
        }
        vertexRunLengthRange[i].second = uint32_t(assembledRunLengthSequence.size());
        vertexRawRange[i].second = uint32_t(assembledRawSequence.size());

        // This was the last vertex.
        if(i == edgeCount) {
            break;
        }

        // Edge.
        edgeRunLengthRange[i].first = uint32_t(assembledRunLengthSequence.size());
        edgeRawRange[i].first = uint32_t(assembledRawSequence.size());
        if(edgeSequences[i].size() > 2*k) {
            for(uint32_t j=uint32_t(k); j!=uint32_t(edgeSequences[i].size()-k); j++) {
                const Base base = edgeSequences[i][j];
                const uint32_t repeatCount = edgeRepeatCounts[i][j];
                assembledRunLengthSequence.push_back(base);
                assembledRepeatCounts.push_back(repeatCount);
                for(uint32_t k=0; k!=repeatCount; k++) {
                    assembledRawSequence.push_back(base);
                }
            }
        }
        edgeRunLengthRange[i].second = uint32_t(assembledRunLengthSequence.size());
        edgeRawRange[i].second = uint32_t(assembledRawSequence.size());
    }



    // If requested, write out details in html format.
    if(htmlPointer) {
        ostream& html = *htmlPointer;

        // Write a title.
        html <<
            "<h1>Assembly graph edge <a href="
            "'exploreAssemblyGraph?edgeId=" << edgeId <<
            "&maxDistance=6&detailed=on&sizePixels=1600&timeout=30'>" <<
            edgeId << "</a></h1>";

#if 0
        // Write parent and child vertices in the assembly graph.
        html << "<p>Parent vertices in the assembly graph:";
        for(const auto parentEdge: assemblyGraph.edgesByTarget[vertexId]) {
            const AssemblyGraph::VertexId parent = assemblyGraph.edges[parentEdge].source;
            html <<
                " <a href='exploreAssemblyGraphVertex?vertexId=" << parent << "'>"
                << parent << "</a>";
        }

        html << "<p>Child vertices in the assembly graph:";
        for(const auto childEdge: assemblyGraph.edgesBySource[vertexId]) {
            const AssemblyGraph::VertexId child = assemblyGraph.edges[childEdge].target;
            html <<
                " <a href='exploreAssemblyGraphVertex?vertexId=" << child << "'>"
                << child << "</a>";
        }
#endif
        // Assembled run-length sequence.
        html << "<p>Assembled run-length sequence (" << assembledRunLengthSequence.size() <<
            " bases):<br><span style='font-family:courier'>";
        copy(assembledRunLengthSequence.begin(), assembledRunLengthSequence.end(),
            ostream_iterator<Base>(html));
        html << "<br>";
        const uint32_t maxRepeatCount =
            *max_element(assembledRepeatCounts.begin(), assembledRepeatCounts.end());
        for(size_t j=0; j<assembledRepeatCounts.size(); j++) {
            const uint32_t repeatCount = assembledRepeatCounts[j];
            html << repeatCount%10;
        }
        if(maxRepeatCount >= 10) {
            html << "<br>";
            for(size_t j=0; j<assembledRepeatCounts.size(); j++) {
                const uint32_t repeatCount = assembledRepeatCounts[j];
                const uint32_t digit = (repeatCount/10) % 10;
                if(digit == 0) {
                    html << "&nbsp;";
                } else {
                    html << digit;
                }
            }
        }
        if(maxRepeatCount >= 100) {
            html << "<br>";
            for(size_t j=0; j<assembledRepeatCounts.size(); j++) {
                const uint32_t repeatCount = assembledRepeatCounts[j];
                if(repeatCount >= 1000) {
                    html << "*";
                } else {
                    const uint32_t digit = (repeatCount/100) % 10;
                    if(digit == 0) {
                        html << "&nbsp;";
                    } else {
                        html << digit;
                    }
                }
            }
        }
        html << "</span>";

        // Assembled raw sequence.
        html << "<p>Assembled raw sequence (" << assembledRawSequence.size() <<
            " bases):<br><span style='font-family:courier'>";
        copy(assembledRawSequence.begin(), assembledRawSequence.end(),
            ostream_iterator<Base>(html));
        html << "</span>";




        // Write a table with a row for each marker graph vertex or edge
        // in the marker graph chain.
        html <<
            "<p>This vertex of the assembly graph corresponds to a chain of " <<
            edgeIds.size() << " edges in the marker graph. "
            "The table below shows consensus sequences "
            "for the vertices and edges of this chain of the marker graph. "
            "All vertex and edge ids in the table refer to the marker graph."
            "<p><table><tr>"
            "<th rowspan=2>Vertex<br>or<br>edge"
            "<th rowspan=2>Id"
            "<th rowspan=2>Coverage"
            "<th colspan=4>Run-length"
            "<th colspan=3>Raw"
            "<tr>"
            "<th>Offset"
            "<th>Begin"
            "<th>End"
            "<th>Sequence"
            "<th>Begin"
            "<th>End"
            "<th>Sequence";
        const string urlPrefix = "exploreMarkerGraph?vertexId=";
        const string urlSuffix =
            "&maxDistance=5"
            "&detailed=on"
            "&minCoverage=3"
            "&minConsensus=3"
            "&sizePixels=600&timeout=10"
            "&useStoredConnectivity=on"
            "&showVertexId=on"
            "&showAssembledSequence=on";
        for(size_t i=0; ; i++) {

            // Vertex.
            const GlobalMarkerGraphVertexId vertexId = vertexIds[i];
            const string url = urlPrefix + to_string(vertexId) + urlSuffix;
            const vector<Base>& vertexSequence = vertexSequences[i];
            const vector<uint32_t>& vertexRepeatCount = vertexRepeatCounts[i];
            // const uint32_t maxVertexRepeatCount =
            //     *std::max_element(vertexRepeatCount.begin(), vertexRepeatCount.end());
            html <<
                 "<tr><td>Vertex" <<
                "<td class=centered><a href='" << url << "'>" << vertexId << "</a>"
                "<td class=centered>" << vertexCoverage[i] <<
                "<td class=centered>" << vertexOffsets[i] <<
                "<td class=centered>" << vertexRunLengthRange[i].first <<
                "<td class=centered>" << vertexRunLengthRange[i].second <<
                "<td style='font-family:courier'>";
            for(size_t j=0; j<vertexSequence.size(); j++) {
                if(j==vertexAssembledPortion[i].first && vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "<span style='background-color:LightGreen'>";
                }
                html << vertexSequence[j];
                if(j==vertexAssembledPortion[i].second-1  && vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "</span>";
                }
            }
            html << "<br>";
            for(size_t j=0; j<vertexSequence.size(); j++) {
                const uint32_t repeatCount = vertexRepeatCount[j];
                if(j==vertexAssembledPortion[i].first && vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "<span style='background-color:LightGreen'>";
                }
                if(repeatCount < 10) {
                    html << repeatCount % 10;
                } else {
                    html << "*";
                }
                if(j==vertexAssembledPortion[i].second-1 && vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "</span>";
                }
            }
            html <<
                "<td class=centered>" << vertexRawRange[i].first <<
                "<td class=centered>" << vertexRawRange[i].second <<
                "<td style='font-family:courier'>";
            for(size_t j=0; j<vertexSequence.size(); j++) {
                if(j==vertexAssembledPortion[i].first && vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "<span style='background-color:LightGreen'>";
                }
                const Base b = vertexSequence[j];
                const uint32_t repeatCount = vertexRepeatCount[j];
                for(uint32_t k=0; k<repeatCount; k++) {
                    html << b;
                }
                if(j==vertexAssembledPortion[i].second-1 && vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "</span>";
                }
            }

            // This was the last vertex.
            if(i == edgeCount) {
                break;
            }



            // Edge.
            const GlobalMarkerGraphEdgeId edgeId = edgeIds[i];
            const MarkerGraph::Edge& edge =
                markerGraph.edges[edgeId];
            const string sourceUrl = urlPrefix + to_string(edge.source) + urlSuffix;
            const string targetUrl = urlPrefix + to_string(edge.target) + urlSuffix;
            const vector<Base>& edgeSequence = edgeSequences[i];
            const vector<uint32_t>& edgeRepeatCount = edgeRepeatCounts[i];
            const size_t edgeSequenceLength = edgeSequence.size();
            CZI_ASSERT(edgeRepeatCount.size() == edgeSequenceLength);
            // const uint32_t maxEdgeRepeatCount =
            //    *std::max_element(edgeRepeatCount.begin(), edgeRepeatCount.end());
            // CZI_ASSERT(maxEdgeRepeatCount < 10);  // For now. Add additional code when this fails.
            html <<
                "<tr><td>Edge<td class=centered>" << edgeId <<
                "<td class=centered>" << edgeCoverage[i] <<
                "<td class=centered>" << edgeOffsets[i] <<
                "<td class=centered>" << edgeRunLengthRange[i].first <<
                "<td class=centered>" << edgeRunLengthRange[i].second <<
                "<td style='font-family:courier'>";
            for(size_t j=0; j<edgeSequenceLength; j++) {
                if(edgeSequenceLength>2*k && j==k) {
                    html << "<span style='background-color:LightGreen'>";
                }
                html << edgeSequence[j];
                if(edgeSequenceLength>2*k && j==edgeSequenceLength-k-1) {
                    html << "</span>";
                }
            }
            html << "<br>";
            for(size_t j=0; j<edgeSequenceLength; j++) {
                const uint32_t repeatCount = edgeRepeatCount[j];
                if(edgeSequenceLength>2*k && j==k) {
                    html << "<span style='background-color:LightGreen'>";
                }
                if(repeatCount < 10) {
                    html << repeatCount % 10;
                } else {
                    html << "*";
                }
                if(edgeSequenceLength>2*k && j==edgeSequenceLength-k-1) {
                    html << "</span>";
                }
            }
            html <<
                "<td class=centered>" << edgeRawRange[i].first <<
                "<td class=centered>" << edgeRawRange[i].second <<
                "<td style='font-family:courier'>";
            for(size_t j=0; j<edgeSequenceLength; j++) {
                const Base b = edgeSequence[j];
                const uint32_t repeatCount = edgeRepeatCount[j];
                if(edgeSequenceLength>2*k && j==k) {
                    html << "<span style='background-color:LightGreen'>";
                }
                for(uint32_t k=0; k<repeatCount; k++) {
                    html << b;
                }
                if(edgeSequenceLength>2*k && j==edgeSequenceLength-k-1) {
                    html << "</span>";
                }
            }
         }
    }

}
