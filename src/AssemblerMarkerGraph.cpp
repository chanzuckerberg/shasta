#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "LocalMarkerGraph.hpp"
#include "LocalReadGraph.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;



// Compute a local marker graph for a set of oriented reads.
void Assembler::createLocalMarkerGraph(
    const vector< pair<ReadId, Strand> >& readIdsAndStrands,
    bool alignAllPairs,
    size_t alignmentMaxSkip,
    size_t alignmentMaxVertexCountPerKmer,
    size_t minAlignedMarkerCount,
    size_t minCoverage,
    size_t minConsensus)
{
    vector<OrientedReadId> orientedReadIds;
    for(const auto& p: readIdsAndStrands) {
        checkReadId(p.first);
        orientedReadIds.push_back(OrientedReadId(p.first, p.second));
    }
    createLocalMarkerGraph(orientedReadIds, alignAllPairs,
        alignmentMaxSkip, alignmentMaxVertexCountPerKmer,
        minAlignedMarkerCount, minCoverage, minConsensus);
}
void Assembler::createLocalMarkerGraph(
    const vector<OrientedReadId>& orientedReadIds,
    bool alignAllPairs,
    size_t alignmentMaxSkip,
    size_t alignmentMaxVertexCountPerKmer,
    size_t minAlignedMarkerCount,
    size_t minCoverage,
    size_t minConsensus)
{
    checkReadsAreOpen();
    checkKmersAreOpen();
    checkMarkersAreOpen();
    if(!alignAllPairs) {
        checkOverlapsAreOpen();
    }

    // Flag to control debug output.
    const bool debug = true;


    if(debug) {
        cout << "Creating a local marker graph the following ";
        cout << orientedReadIds.size() << " oriented reads:\n";
        for(size_t i=0; i<orientedReadIds.size(); i++) {
            const OrientedReadId orientedReadId = orientedReadIds[i];
            cout << i << " " << orientedReadId << "\n";
        }
        cout << flush;
    }

    // Extract the sequences of the oriented reads.
    vector<LongBaseSequence> sequences;
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        sequences.push_back(reads[orientedReadId.getReadId()]);
        if(orientedReadId.getStrand() == 1) {
            sequences.back().reverseComplement();
        }
    }

    // Extract the k-mer occurrences sorted by position
    // for these oriented reads.
    vector< vector<Marker> > markersInGraphSortedByPosition(orientedReadIds.size());
    vector< vector<Marker> > markersInGraphSortedByKmerId(orientedReadIds.size());
    for(size_t localOrientedReadId=0; localOrientedReadId!=orientedReadIds.size(); ++localOrientedReadId) {
        const OrientedReadId orientedReadId = orientedReadIds[localOrientedReadId];
        getMarkers(orientedReadId, markersInGraphSortedByPosition[localOrientedReadId]);
        markersInGraphSortedByKmerId[localOrientedReadId] = markersInGraphSortedByPosition[localOrientedReadId];
        sort(
            markersInGraphSortedByKmerId[localOrientedReadId].begin(),
            markersInGraphSortedByKmerId[localOrientedReadId].end(),
            OrderMarkersByKmerId());
    }

    // Construct the initial local marker graph.
    LocalMarkerGraph graph(assemblerInfo->k, orientedReadIds, sequences, markersInGraphSortedByPosition,
        minCoverage, minConsensus);



    // Add the alignments. This merges vertices whose markers are aligned.
    Alignment alignment;
    AlignmentGraph alignmentGraph;
    const bool alignDebug = false;
    if(alignAllPairs) {

        // Align all pairs of oriented reads in the graph.
        for(uint32_t i0=0; i0<uint32_t(orientedReadIds.size()-1); i0++) {
            for(uint32_t i1=i0+1; i1<uint32_t(orientedReadIds.size()); i1++) {

                // Compute the alignment.
                align(
                    markersInGraphSortedByKmerId[i0],
                    markersInGraphSortedByKmerId[i1],
                    int(alignmentMaxSkip),
                    alignmentMaxVertexCountPerKmer,
                    alignmentGraph,
                    alignDebug,
                    alignment);

                // If the alignment is too short, skip.
                if(alignment.ordinals.size() < minAlignedMarkerCount) {
                    continue;
                }

                // Merge alignment vertices.
                cout << "Alignment of " << orientedReadIds[i0] << " " << orientedReadIds[i1];
                cout << " of length " << alignment.ordinals.size() << endl;
                for(const auto& ordinals: alignment.ordinals) {
                    graph.mergeVertices(
                        i0, ordinals.first,
                        i1, ordinals.second);
                    // cout << ordinals.first << " " << ordinals.second << endl;
                }
            }
        }

    } else {
        // Only align pairs of oriented reads in the graph
        // for which we have an overlap.
        CZI_ASSERT(0);
    }
    if(debug) {
        cout << "The initial local marker graph graph  has ";
        cout << boost::num_vertices(graph) << " vertices and ";
        cout << boost::num_edges(graph) << " edges." << endl;
    }



    graph.sortAndSplitAmbiguousVertices();
    graph.fillEdgeData();
    // graph.computeOptimalSpanningTree();
    // graph.removeWeakNonSpanningTreeEdges();
    graph.pruneWeakLeaves();
    if(debug) {
        cout << "The local marker graph graph  has ";
        cout << boost::num_vertices(graph) << " vertices and ";
        cout << boost::num_edges(graph) << " edges." << endl;
        const bool addEdgeLabels = true;
        graph.write("LocalMarkerGraph.dot", addEdgeLabels);
    }

#if 0
    const vector< pair<Base, int> > longestSequence = graph.extractLongestSequence();
    // Write out the sequence.
    ofstream fastaOut("LongestSequence.txt");
    for(const auto& p: longestSequence) {
        fastaOut << p.first;
    }
    fastaOut << endl;
    for(const auto& p: longestSequence) {
        const int coverage = p.second;
        if(coverage < 10) {
            fastaOut << coverage;
        } else {
            fastaOut << "*";
        }
    }
    fastaOut << endl;
#endif
}



// Create the local marker graph that corresponds to a local read graph
// constructed starting at a given oriented read and extending out
// up to a specified distance.
void Assembler::createLocalMarkerGraph(
    ReadId readId,
    Strand strand,
    size_t minFrequency,            // Minimum number of minHash hits to generate an edge.
    size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
    size_t maxTrim,                 // Maximum left/right trim to generate an edge.
    size_t distance,                // How far to go from starting oriented read.
    size_t alignmentMaxSkip,
    size_t alignmentMaxVertexCountPerKmer,
    size_t minCoverage,
    size_t minConsensus)
{
    // Check that we have what we need.
    checkReadsAreOpen();
    checkReadNamesAreOpen();
    checkKmersAreOpen();
    checkMarkersAreOpen();
    checkOverlapsAreOpen();
    checkAlignmentInfosAreOpen();
    CZI_ASSERT(overlaps.size() == alignmentInfos.size());

    // Create the local read graph.
    LocalReadGraph localReadGraph;
    createLocalReadGraph(OrientedReadId(readId, strand),
        minFrequency, minAlignedMarkerCount, maxTrim, distance,
        localReadGraph);
    const size_t orientedReadIdCount = num_vertices(localReadGraph);
    cout << "The local read graph has " << num_vertices(localReadGraph);
    cout << " vertices and " << num_edges(localReadGraph) << " edges." << endl;
    localReadGraph.write("LocalReadGraph.dot");
    writeLocalReadGraphToFasta(localReadGraph, "LocalReadGraph.fasta");



    // Gather the oriented read ids, sequences, and markers
    // (sorted by position) of the oriented reads in the local read graph.
    // These vectors are indexed by the local oriented read id,
    // an index that runs from 0 through orientedReadIdCount (excluded).
    vector<OrientedReadId> orientedReadIds(orientedReadIdCount);
    vector<LongBaseSequence> sequences(orientedReadIdCount);
    vector< vector<Marker> > markersInGraphSortedByPosition(orientedReadIdCount);
    vector< vector<Marker> > markersInGraphSortedByKmerId(orientedReadIdCount);
    std::map<OrientedReadId, uint32_t> orientedReadIdMap;    // Mapo global to local oriented read id.
    uint32_t localOrientedReadId = 0;
    BGL_FORALL_VERTICES(v, localReadGraph, LocalReadGraph) {

        // Oriented read id.
        const OrientedReadId orientedReadId = OrientedReadId(localReadGraph[v].orientedReadId);
        orientedReadIds[localOrientedReadId] = orientedReadId;
        orientedReadIdMap.insert(make_pair(orientedReadId, localOrientedReadId));

        // Sequence.
        sequences[localOrientedReadId] = reads[orientedReadId.getReadId()];
        if(orientedReadId.getStrand() == 1) {
            sequences[localOrientedReadId].reverseComplement();
        }

        // Markers sorted by position.
        getMarkers(orientedReadId, markersInGraphSortedByPosition[localOrientedReadId]);

        // Markers sorted by kmer id.
        markersInGraphSortedByKmerId[localOrientedReadId] = markersInGraphSortedByPosition[localOrientedReadId];
        sort(
            markersInGraphSortedByKmerId[localOrientedReadId].begin(),
            markersInGraphSortedByKmerId[localOrientedReadId].end(),
            OrderMarkersByKmerId());

        ++localOrientedReadId;
    }



    // Construct the initial local marker graph.
    // The initial local marker graph has a linear chain corresponding
    // to each oriented read.
    LocalMarkerGraph localMarkerGraph(assemblerInfo->k,
        orientedReadIds, sequences, markersInGraphSortedByPosition,
        minCoverage, minConsensus);
    cout << "The initial local marker graph has " << num_vertices(localMarkerGraph);
    cout << " vertices and " << num_edges(localMarkerGraph) << " edges." << endl;



    // Merge vertices of the marker graph corresponding to aligned markers.
    // We loop over all edges of the local read graph.
    // For each edge we compute the alignment between the two oriented read,
    // then merge pairs of marker graph vertices corresponding to aligned markers.
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    AlignmentGraph alignmentGraph;
    const bool alignDebug = false;
    BGL_FORALL_EDGES(e, localReadGraph, LocalReadGraph) {
        CZI_ASSERT(localReadGraph[e].overlap.minHashFrequency >= minFrequency);

        // Get the global and local ids of the oriented reads to align.
        const LocalReadGraph::vertex_descriptor v0 = source(e, localReadGraph);
        const LocalReadGraph::vertex_descriptor v1 = target(e, localReadGraph);
        const OrientedReadId orientedReadId0 = OrientedReadId(localReadGraph[v0].orientedReadId);
        const OrientedReadId orientedReadId1 = OrientedReadId(localReadGraph[v1].orientedReadId);
        const uint32_t localOrientedReadId0 = orientedReadIdMap[orientedReadId0];
        const uint32_t localOrientedReadId1 = orientedReadIdMap[orientedReadId1];

        // Compute the alignment.
        align(
            markersInGraphSortedByKmerId[localOrientedReadId0],
            markersInGraphSortedByKmerId[localOrientedReadId1],
            int(alignmentMaxSkip),
            alignmentMaxVertexCountPerKmer,
            alignmentGraph,
            alignDebug,
            alignment);
        alignmentInfo.create(alignment);

        // If the alignment has too few markers, skip it.
        if(alignment.ordinals.size() < minAlignedMarkerCount) {
            continue;
        }

        // If the alignment has too much trim, skip it.
        uint32_t leftTrim;
        uint32_t rightTrim;
        tie(leftTrim, rightTrim) = computeTrim(
            orientedReadId0,
            orientedReadId1,
            markersInGraphSortedByPosition[localOrientedReadId0],
            markersInGraphSortedByPosition[localOrientedReadId1],
            alignmentInfo);
        if(leftTrim>maxTrim || rightTrim>maxTrim) {
            continue;
        }

        // Merge aligned vertices.
        for(const auto& ordinals: alignment.ordinals) {
            localMarkerGraph.mergeVertices(
                localOrientedReadId0, ordinals.first,
                localOrientedReadId1, ordinals.second);
        }
    }

    // Clean up the local marker graph.
    localMarkerGraph.sortAndSplitAmbiguousVertices();
    localMarkerGraph.fillEdgeData();
    // localMarkerGraph.computeOptimalSpanningTree();
    // localMarkerGraph.removeWeakNonSpanningTreeEdges();
    localMarkerGraph.pruneWeakLeaves();

    cout << "The local marker graph has " << num_vertices(localMarkerGraph);
    cout << " vertices and " << num_edges(localMarkerGraph) << " edges." << endl;
    localMarkerGraph.write("MarkerGraph.dot", false);
    localMarkerGraph.write("DetailedMarkerGraph.dot", true);

}
