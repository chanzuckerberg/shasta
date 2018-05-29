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
    LocalReadGraph graph;
    createLocalReadGraph(OrientedReadId(readId, strand),
        minFrequency, minAlignedMarkerCount, maxTrim, distance,
        graph);
    cout << "The local read graph has " << num_vertices(graph);
    cout << " vertices and " << num_edges(graph) << " edges." << endl;
    graph.write("LocalReadGraph.dot");
    writeLocalReadGraphToFasta(graph, "LocalReadGraph.fasta");

    CZI_ASSERT(0);
}
