// Shasta.
#include "Assembler.hpp"
#include "deduplicate.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Seqan.
#include <seqan/align.h>

// Standard library.
#include <map>



// Analyze oriented read paths in the marker graph and in the assembly graph.

// An oriented read always corresponds to a path in the marker graph
// if all edges are used. However, because of edges removed
// during transitive reduction, an oriented read does not correspond
// to a path in the assembly graph.
// The sequence of assembly graph edges (segments) encountered
// by an oriented read is referred to below as the pseudo-path
// of that oriented read. It is a sequence of assembly graph
// edges, but not necessarily a path in the assembly graph.

// In the code and comments below, "segment" is synonym for
// "assembly graph edge".

void Assembler::analyzeOrientedReadPaths(int readGraphCreationMethod) const
{
    using SegmentId = AssemblyGraph::EdgeId;
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    const uint64_t segmentCount = assemblyGraph.edges.size();



    // Compute the pseudo-path of each oriented read.
    // This vector is indexed by OrientedReadId::getValue().
    vector< vector<SegmentId> > pseudoPaths(2*readCount());
    vector<MarkerGraph::EdgeId> markerGraphPath;
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            vector<SegmentId>& pseudoPath =
                pseudoPaths[orientedReadId.getValue()];

            // Find the marker graph path of this oriented read.
            computeOrientedReadMarkerGraphPath(
                orientedReadId,
                0, uint32_t(markers.size(orientedReadId.getValue())-1),
                markerGraphPath);

            // Loop over the marker graph path.
            SegmentId previousSegmentId =
                std::numeric_limits<SegmentId>::max();
            for(const MarkerGraph::EdgeId markerGraphEdgeId: markerGraphPath) {

                // Get the corresponding segments.
                const span<const pair<SegmentId, uint32_t> > v =
                    assemblyGraph.markerToAssemblyTable[markerGraphEdgeId];

                // If no segments, skip.
                if(v.size() == 0) {
                    continue;
                }

                // If detangling was used, there can be more than one,
                // and we don't want this here.
                SHASTA_ASSERT(v.size() == 1);

                // There is only one segment.
                const SegmentId segmentId = v.front().first;

                // If same as the previous, slip.
                if(segmentId == previousSegmentId) {
                    continue;
                }

                // This is the next segment edge encountered
                // by this oriented read along its marker graph path.
                // Add it to the pseudo-path.
                pseudoPath.push_back(segmentId);
                previousSegmentId = segmentId;
            }
        }
    }



    // Write a csv file with the pseudo-path of each oriented read.
    {
        ofstream csv("PseudoPaths.csv");
        for(ReadId readId=0; readId<readCount(); readId++) {
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);
                csv << orientedReadId << ",";

                vector<SegmentId>& pseudoPath =
                    pseudoPaths[orientedReadId.getValue()];
                for(const SegmentId segmentId: pseudoPath) {
                    csv << segmentId << ",";
                }
                csv << "\n";
            }
        }
    }



    // Create the pseudo-path table which contains, for each segment,
    // its occurrences in oriented read pseudo-paths.
    // For each segmentId, we store a vector of pairs (orientedReadId, index) such that
    // pseudoPaths[orientedReadId.getValue()][index] == segmentId
    vector< vector< pair<OrientedReadId, uint64_t> > >  pseudoPathTable(segmentCount);
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const vector<SegmentId>& pseudoPath =
                pseudoPaths[orientedReadId.getValue()];

            for(uint64_t index=0; index<pseudoPath.size(); index++) {
                const SegmentId segmentId = pseudoPath[index];
                pseudoPathTable[segmentId].push_back(make_pair(orientedReadId, index));
            }
        }
    }



    // Write out the pseudo-path table.
    {
        ofstream csv("PseudoPathTable.csv");
        csv << "Segment,OrientedRead,Position\n";
        for(SegmentId segmentId=0; segmentId<segmentCount; segmentId++) {
            const vector< pair<OrientedReadId, uint64_t> >& v = pseudoPathTable[segmentId];
            for(const auto& p: v) {
                csv << segmentId << ",";
                csv << p.first << ",";
                csv << p.second << "\n";
            }

        }
    }



    // Gather all pairs of oriented reads that occur in non-branching segments.
    // A segment (assembly graph edge) v0->v1 is no branching
    // if in-degree(v0)<2 and out_degree(v1)<2.
    vector< pair<OrientedReadId, OrientedReadId> > orientedReadPairs;
    for(SegmentId segmentId=0; segmentId<segmentCount; segmentId++) {

        // If this is a branching segment, skip it.
        const AssemblyGraph::Edge& segment = assemblyGraph.edges[segmentId];
        const AssemblyGraph::VertexId v0 = segment.source;
        if(assemblyGraph.inDegree(v0) > 1) {
            continue;
        }
        const AssemblyGraph::VertexId v1 = segment.target;
        if(assemblyGraph.outDegree(v1) > 1) {
            continue;
        }

        // Loop over pairs of  oriented reads that have this segment on their pseudo-path.
        const vector< pair<OrientedReadId, uint64_t> >& v = pseudoPathTable[segmentId];
        for(uint64_t i0=0; i0<v.size(); i0++) {
            const OrientedReadId orientedReadId0 = v[i0].first;
            for(uint64_t i1=i0+1; i1<v.size(); i1++) {
                const OrientedReadId orientedReadId1 = v[i1].first;
                orientedReadPairs.push_back(make_pair(orientedReadId0, orientedReadId1));
            }
        }
    }
    deduplicate(orientedReadPairs);
    cout << "Found " << orientedReadPairs.size() <<
        " oriented read pairs in non-branching segments." << endl;



    // For each such pair of oriented reads, compute an alignment between their pseudo-paths.
    ofstream alignmentsCsv("Alignments.csv");
    for(const auto& p: orientedReadPairs) {
        const OrientedReadId orientedReadId0 = p.first;
        const OrientedReadId orientedReadId1 = p.second;

        // Get the pseudo-paths of these two oriented reads.
        const vector<AssemblyGraph::EdgeId>& pseudoPath0 = pseudoPaths[orientedReadId0.getValue()];
        const vector<AssemblyGraph::EdgeId>& pseudoPath1 = pseudoPaths[orientedReadId1.getValue()];

        // Use SeqAn to compute an alignment free at both ends.
        // https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html
        using namespace seqan;

        // Hide shasta::Alignment.
        using seqan::Alignment;

        // An oriented read is represented by its pseudo-path.
        // We want to align a pair of such sequences.
        using TSequence = String<AssemblyGraph::EdgeId>;

        // Other SeqAn types we need.
        using TStringSet = StringSet<TSequence>;
        using TDepStringSet = StringSet<TSequence, Dependent<> >;
        using TAlignGraph = Graph<Alignment<TDepStringSet> >;

        // Construct the sequences we want to pass to SeqAn.
        // Add 100 to all segment ids to avoid collision with the
        // value 45, used by SeqAn to represent gaps.
        TSequence seq0;
        for(const AssemblyGraph::EdgeId segmentId: pseudoPath0) {
            appendValue(seq0, segmentId+100);
        }
        TSequence seq1;
        for(const AssemblyGraph::EdgeId segmentId: pseudoPath1) {
            appendValue(seq1, segmentId+100);
        }

        // Store them in a SeqAn string set.
        TStringSet sequences;
        appendValue(sequences, seq0);
        appendValue(sequences, seq1);

        // Compute the alignment.
        TAlignGraph graph(sequences);
        const int matchScore = 1;
        const int mismatchScore = -1;
        const int gapScore = -1;
        globalAlignment(
                graph,
                Score<int, Simple>(matchScore, mismatchScore, gapScore),
                AlignConfig<true, true, true, true>(),
                LinearGaps());

        // Extract the alignment from the graph.
        // This creates a single sequence consisting of the two rows
        // of the alignment, concatenated.
        TSequence align;
        convertAlignment(graph, align);
        const uint64_t totalAlignmentLength = seqan::length(align);
        SHASTA_ASSERT((totalAlignmentLength % 2) == 0);    // Because we are aligning two sequences.
        const uint64_t alignmentLength = totalAlignmentLength / 2;

        // Write out the alignment.
        uint64_t index = 0;
        for(uint64_t i=0; i<2; i++) {
            alignmentsCsv << (i==0 ? orientedReadId0 : orientedReadId1) << ",";
            for(uint64_t j=0; j<alignmentLength; j++, index++) {
                const uint64_t value = align[index];
                if(value == 45) {
                    alignmentsCsv << "-";
                } else {
                    alignmentsCsv << value - 100;
                }
                alignmentsCsv << ",";
            }
            alignmentsCsv << "\n";
        }
        alignmentsCsv << "Alignment,";
        for(uint64_t j=0; j<alignmentLength; j++) {
            const uint64_t value0 = align[j];
            const uint64_t value1 = align[j+alignmentLength];
            if(value0==45 and value1==45) {
                alignmentsCsv << "?";    // This should never happen.
            } else if(value0==45 or value1==45) {
                // Gap on one of the two.
                alignmentsCsv << "-";
            } else if(value0 == value1) {
                // Match.
                alignmentsCsv << ".";
            } else {
                // Mismatch.
                alignmentsCsv << "*";
            }
            alignmentsCsv << ",";
        }
        alignmentsCsv << "\n";
    }



#if 0
    // Loop over read graph edges and the corresponding alignments.
    // We process read graph edges in pairs.
    // In each pair, the second edge is the reverse complement of the first.
    const uint64_t readGraphEdGeCount =
        (readGraphCreationMethod==0) ? readGraph.edges.size() : directedReadGraph.edges.size();
    for(uint64_t readGraphEdgeId=0; readGraphEdgeId!=readGraphEdGeCount; readGraphEdgeId+=2) {

        // Get the oriented read ids
        array<OrientedReadId, 2> orientedReadIds;
        if(readGraphCreationMethod == 0) {

            // We use the undirected read graph.
            const ReadGraphEdge& readGraphEdge = readGraph.edges[readGraphEdgeId];

            // Check that the next edge is the reverse complement of
            // this edge.
            {
                const ReadGraphEdge& readGraphNextEdge = readGraph.edges[readGraphEdgeId + 1];
                array<OrientedReadId, 2> nextEdgeOrientedReadIds = readGraphNextEdge.orientedReadIds;
                nextEdgeOrientedReadIds[0].flipStrand();
                nextEdgeOrientedReadIds[1].flipStrand();
                SHASTA_ASSERT(nextEdgeOrientedReadIds == readGraphEdge.orientedReadIds);
            }


            if(readGraphEdge.crossesStrands) {
                continue;
            }
            orientedReadIds = readGraphEdge.orientedReadIds;
            SHASTA_ASSERT(orientedReadIds[0] < orientedReadIds[1]);

            // If either of the reads is flagged chimeric, skip it.
            if( readFlags[orientedReadIds[0].getReadId()].isChimeric ||
                readFlags[orientedReadIds[1].getReadId()].isChimeric) {
                continue;
            }
        } else if(readGraphCreationMethod == 1) {

            // We use the directed read graph.
            const DirectedReadGraphEdge& edge = directedReadGraph.getEdge(readGraphEdgeId);

            // Sanity checks.
            // Pairs of reverse complemented adges are stored consecutively.
            SHASTA_ASSERT(edge.reverseComplementedEdgeId == readGraphEdgeId+1);
            const DirectedReadGraphEdge& nextEdge = directedReadGraph.getEdge(readGraphEdgeId+1);
            SHASTA_ASSERT(nextEdge.reverseComplementedEdgeId == readGraphEdgeId);
            SHASTA_ASSERT(nextEdge.keep == edge.keep);
            SHASTA_ASSERT(nextEdge.isConflict == edge.isConflict);

            // Skip if not marked as "keep".
            if(edge.keep == 0) {
                continue;
            }

            // Skip if marked as "conflict".
            if(edge.isConflict == 1) {
                continue;
            }


            // Get the oriented read ids.
            const DirectedReadGraph::VertexId v0 = directedReadGraph.source(readGraphEdgeId);
            const DirectedReadGraph::VertexId v1 = directedReadGraph.target(readGraphEdgeId);
            orientedReadIds[0] = OrientedReadId(OrientedReadId::Int(v0));
            orientedReadIds[1] = OrientedReadId(OrientedReadId::Int(v1));

        } else {
            throw runtime_error("Invalid read graph creation method " + to_string(readGraphCreationMethod));
        }


        // Get the pseudo-paths of these two oriented reads.
        const vector<AssemblyGraph::EdgeId>& pseudoPath0 = pseudoPaths[orientedReadIds[0].getValue()];
        const vector<AssemblyGraph::EdgeId>& pseudoPath1 = pseudoPaths[orientedReadIds[1].getValue()];

        cout << "\n";
        cout << orientedReadIds[0];
        for(const AssemblyGraph::EdgeId segmentId: pseudoPath0) {
            cout << " " << segmentId;
        }
        cout << "\n";
        cout << orientedReadIds[1];
        for(const AssemblyGraph::EdgeId segmentId: pseudoPath1) {
            cout << " " << segmentId;
        }
        cout << "\n";



        // Use SeqAn to compute an alignment free at both ends.
        // https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html
        using namespace seqan;

        // Hide shasta::Alignment.
        using seqan::Alignment;

        // An oriented read is represented by its pseudo-path.
        // We want to align a pair of such sequences.
        using TSequence = String<AssemblyGraph::EdgeId>;

        // Other SeqAn types we need.
        using TStringSet = StringSet<TSequence>;
        using TDepStringSet = StringSet<TSequence, Dependent<> >;
        using TAlignGraph = Graph<Alignment<TDepStringSet> >;

        // Construct the sequences we want to pass to SeqAn.
        // Add 100 to all segment ids to avoid collision with the
        // value 45, used by SeqAn to represent gaps.
        TSequence seq0;
        for(const AssemblyGraph::EdgeId segmentId: pseudoPath0) {
            appendValue(seq0, segmentId+100);
        }
        TSequence seq1;
        for(const AssemblyGraph::EdgeId segmentId: pseudoPath1) {
            appendValue(seq1, segmentId+100);
        }

        // Store them in a SeqAn string set.
        TStringSet sequences;
        appendValue(sequences, seq0);
        appendValue(sequences, seq1);

        // Compute the alignment.
        TAlignGraph graph(sequences);
        const int matchScore = 3;
        const int mismatchScore = -3;
        const int gapScore = -1;
        globalAlignment(
                graph,
                Score<int, Simple>(matchScore, mismatchScore, gapScore),
                AlignConfig<true, true, true, true>(),
                LinearGaps());

        // Extract the alignment from the graph.
        // This creates a single sequence consisting of the two rows
        // of the alignment, concatenated.
        TSequence align;
        convertAlignment(graph, align);
        const uint64_t totalAlignmentLength = seqan::length(align);
        SHASTA_ASSERT((totalAlignmentLength % 2) == 0);    // Because we are aligning two sequences.
        const uint64_t alignmentLength = totalAlignmentLength / 2;
        cout << "Alignment length " << alignmentLength << endl;

        // Write out the alignment.
        uint64_t index = 0;
        for(uint64_t i=0; i<2; i++) {
            for(uint64_t j=0; j<alignmentLength; j++, index++) {
                const uint64_t value = align[index];
                if(value == 45) {
                    cout << "-";
                } else {
                    cout << value - 100;
                }
                cout << ",";
            }
            cout << endl;
        }
        for(uint64_t j=0; j<alignmentLength; j++) {
            const uint64_t value0 = align[j];
            const uint64_t value1 = align[j+alignmentLength];
            if(value0==45 and value1==45) {
                cout << "?";    // This should never happen.
            } else if(value0==45 or value1==45) {
                // Gap on one of the two.
                cout << "-";
            } else if(value0 == value1) {
                // Match.
                cout << "|";
            } else {
                // Mismatch.
                cout << "*";
            }
            cout << ",";
        }
        cout << endl;
    }
#endif



#if 0
    // Find segments that are encountered more than once in the sequence of any
    // oriented reads.
    vector<bool> isDuplicate(assemblyGraph.edges.size(), false);
    vector<AssemblyGraph::EdgeId> deduplicatedSegments;
    vector<uint64_t> frequency;
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const vector<AssemblyGraph::EdgeId>& segments =
                orientedReadSegments[orientedReadId.getValue()];

            // See if there are any duplicates.
            deduplicatedSegments = segments;
            deduplicateAndCount(deduplicatedSegments, frequency);
            for(uint64_t i=0; i<deduplicatedSegments.size(); i++) {
                if(frequency[i] > 1) {
                    isDuplicate[deduplicatedSegments[i]] = true;
                    //
                    cout << "Segment " << deduplicatedSegments[i] << " encountered " <<
                        frequency[i] << " times by oriented read " <<
                        orientedReadId << endl;
                    //
                }
            }
        }
    }
    cout << "The following assembly graph edges (segments) are encountered "
        "more than once by one or more oriented reads:" << endl;
    uint64_t duplicateCount = 0;
    for(AssemblyGraph::EdgeId segmentId=0; segmentId<segmentCount; segmentId++) {
        if(isDuplicate[segmentId]) {
            cout << segmentId << " ";
            ++duplicateCount;
        }
    }
    cout << endl;
    cout << duplicateCount << " such segments found." << endl;



    // Table to contain, for each segment, its occurrences in
    // orientedReadSegments.
    // For each segmentId, we store a vector of pairs (orientedReadId, index)
    // such that
    // orientedReadSegments[orientedReadId.getValue()][index] == segmentId
    vector< vector< pair<OrientedReadId, uint64_t> > >
        segmentTable(segmentCount);
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const vector<AssemblyGraph::EdgeId>& segments =
                orientedReadSegments[orientedReadId.getValue()];

            for(uint64_t index=0; index<segments.size(); index++) {
                const AssemblyGraph::EdgeId segmentId = segments[index];
                segmentTable[segmentId].push_back(make_pair(orientedReadId, index));
            }
        }
    }



    // Write out the segment table.
    {
        ofstream csv("SegmentTable.csv");
        csv << "Segment,OrientedRead,Position\n";
        for(AssemblyGraph::EdgeId segmentId=0; segmentId<segmentCount; segmentId++) {
            const vector< pair<OrientedReadId, uint64_t> >& v = segmentTable[segmentId];
            for(const auto& p: v) {
                csv << segmentId << ",";
                csv << p.first << ",";
                csv << p.second << "\n";
            }

        }
    }



    // For each segment, analyze the paths of oriented reads that
    // encounter that segment.
    for(AssemblyGraph::EdgeId segmentId=0; segmentId<segmentCount; segmentId++) {
        const vector< pair<OrientedReadId, uint64_t> >& v = segmentTable[segmentId];

        if(segmentId != 1489) {
            continue;       // For debugging.
        }

        // Create a graph in which each edge corresponds of
        // successive segments encountered by an oriented read.
        using boost::adjacency_list;
        using boost::setS;  // No parallel edges.
        using boost::vecS;
        using boost::bidirectionalS;
        using Graph = adjacency_list<setS, vecS, bidirectionalS, AssemblyGraph::EdgeId>;
        Graph graph;
        std::map<AssemblyGraph::EdgeId, Graph::vertex_descriptor> vertexMap;

        for(const auto& p: v) {
            const OrientedReadId orientedReadId = p.first;

            // Access the sequence of segments encountered by this read.
            const vector<AssemblyGraph::EdgeId>& segments =
                orientedReadSegments[orientedReadId.getValue()];
            for(uint64_t i=1; i<segments.size(); i++) {
                const AssemblyGraph::EdgeId segmentId0 = segments[i-1];
                const AssemblyGraph::EdgeId segmentId1 = segments[i];

                // Locate the vertices corresponding to these segments,
                // creating them if they don't exist.
                auto it0 = vertexMap.find(segmentId0);
                if(it0 == vertexMap.end()) {
                    bool wasInserted = false;
                    tie(it0, wasInserted) = vertexMap.insert(
                        make_pair(segmentId0, add_vertex(segmentId0, graph)));
                    SHASTA_ASSERT(wasInserted);
                    SHASTA_ASSERT(it0->first == segmentId0);
                }
                const Graph::vertex_descriptor v0 = it0->second;
                auto it1 = vertexMap.find(segmentId1);
                if(it1 == vertexMap.end()) {
                    bool wasInserted = false;
                    tie(it1, wasInserted) = vertexMap.insert(
                        make_pair(segmentId1, add_vertex(segmentId1, graph)));
                    SHASTA_ASSERT(wasInserted);
                    SHASTA_ASSERT(it1->first == segmentId1);
                }
                const Graph::vertex_descriptor v1 = it1->second;

                // Add the edge.
                add_edge(v0, v1, graph);
            }
        }

        // Write out the graph.
        ofstream graphOut("ReadPaths-" + to_string(segmentId) + ".dot");
        graphOut << "digraph G {\n";
        BGL_FORALL_EDGES(e, graph, Graph) {
            graphOut << graph[source(e, graph)] << "->";
            graphOut << graph[target(e, graph)] << ";\n";
        }
        graphOut << "}\n";
    }
#endif
}
