// Shasta.
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "MetaMarkerGraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Seqan.
#include <seqan/align.h>

// Standard library.
#include "array.hpp"
#include <map>



// An experimental De Bruijn graph used below to align pseudo-paths.
// Each vertex stores a sequence of k segment ids that occur
// consecutively in one or more pseudo-paths.
// An directed edge is v0->v1 is created between vertices v0 and v1
// if the last k-1 segments of v0 are the same as the first
// k-1 segment of v1.
namespace shasta {
    namespace pseudoPaths {
        template<int k> class DeBruijnGraph;
        template<int k> class DeBruijnGraphVertex;
        template<int k> class DeBruijnGraphEdge;
        template<int k> using DeBruijnGraphBaseClass =
            boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            DeBruijnGraphVertex<k>,
            DeBruijnGraphEdge<k> >;
    }
}
using namespace shasta::pseudoPaths;

template<int k> class shasta::pseudoPaths::DeBruijnGraphVertex {
public:

    // The k segment ids associated with this vertex.
    array<uint64_t, k> segmentIds;

    // The oriented reads that have this sequence of segment ids
    // on their pseudo-path, each stored with the starting position
    // in the pseudo-path.
    vector< pair<OrientedReadId, uint64_t > > orientedReadIds;

    void writeGraphviz(ostream& s) const
    {
        s << "\"";
        for(uint64_t i=0; i<k; i++) {
            s << segmentIds[i];
            if(i != k-1) {
                s << "\\n";
            }
        }
        s << "\"";
    }
    void write(ostream& s) const
    {
        for(uint64_t i=0; i<k; i++) {
            s << segmentIds[i];
            if(i != k-1) {
                s << ",";
            }
        }
    }
};



template<int k> class shasta::pseudoPaths::DeBruijnGraphEdge {
public:

    // The k-1 segment ids associated with this edge.
    array<uint64_t, k-1> segmentIds;

    // The oriented reads that have this sequence of segment ids
    // on their pseudo-path, each stored with the starting position
    // in the pseudo-path.
    vector< pair<OrientedReadId, uint64_t > > orientedReadIds;

    DeBruijnGraphEdge(const array<uint64_t, k-1>& segmentIds) :
        segmentIds(segmentIds) {}
};



template<int k> class shasta::pseudoPaths::DeBruijnGraph :
    public DeBruijnGraphBaseClass<k> {
public:
    using Graph = DeBruijnGraph;
    using SegmentId = AssemblyGraph::EdgeId;
    using BaseClass = DeBruijnGraphBaseClass<k>;
    using vertex_descriptor = typename BaseClass::vertex_descriptor;
    using edge_descriptor = typename BaseClass::edge_descriptor;

    // The vertices, keyed by the k segment ids.
    std::map< array<uint64_t, k>, vertex_descriptor> vertexMap;



    // Add vertices for an oriented read.
    void addOrientedReadVertices(
        OrientedReadId orientedReadId,
        const vector<SegmentId>& pseudoPath)
    {
        Graph& graph = *this;

        // Loop over possible starting positions such that the k segments
        // starting there are on the pseudo-path.
        for(uint64_t startPosition=0; startPosition+k<=pseudoPath.size(); startPosition++) {

            // Extract the k segments.
            array<uint64_t, k> segmentIds;
            const auto begin = pseudoPath.begin() + startPosition;
            const auto end = begin + k;
            copy(begin, end, segmentIds.begin());

            // Get the vertex corresponding to these k segments, creating it
            // if necessary.
            vertex_descriptor v;
            auto it = vertexMap.find(segmentIds);
            if(it == vertexMap.end()) {
                v = add_vertex(graph);
                graph[v].segmentIds = segmentIds;
                vertexMap.insert(make_pair(segmentIds, v));
            } else {
                v = it->second;
            }
            graph[v].orientedReadIds.push_back(make_pair(orientedReadId, startPosition));
        }
    }


    void createEdges()
    {
        Graph& graph = *this;

        // Index the vertices by their first k-1 segment ids.
        std::map< array<uint64_t, k-1>, vector<vertex_descriptor> > vertexIndex;
        BGL_FORALL_VERTICES_T(v, graph, Graph) {
            const array<uint64_t, k>& segmentIds = graph[v].segmentIds;
            array<uint64_t, k-1> firstSegmentIds;
            const auto begin = segmentIds.begin();
            const auto end = begin + (k-1);
            copy(begin, end, firstSegmentIds.begin());
            vertexIndex[firstSegmentIds].push_back(v);
        }

        // Use the index to create the edges.
        BGL_FORALL_VERTICES_T(v0, graph, Graph) {
            const array<uint64_t, k>& segmentIds0 = graph[v0].segmentIds;
            array<uint64_t, k-1> lastSegmentIds0;
            const auto begin = segmentIds0.begin() + 1;
            const auto end = segmentIds0.end();
            copy(begin, end, lastSegmentIds0.begin());

            for(const vertex_descriptor v1: vertexIndex[lastSegmentIds0]) {
                edge_descriptor e;
                tie(e, ignore) = add_edge(v0, v1, DeBruijnGraphEdge<k>(lastSegmentIds0), graph);
                findOrientedReadIds(e);
            }
        }
    }



    void findOrientedReadIds(edge_descriptor e)
    {
        Graph& graph = *this;
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const auto& orientedReadIds0 = graph[v0].orientedReadIds;
        const auto& orientedReadIds1 = graph[v1].orientedReadIds;

        // This could be done better.
        for(const auto& p0: orientedReadIds0) {
            auto p1 = p0;
            ++p1.second;
            if(binary_search(orientedReadIds1.begin(), orientedReadIds1.end(), p1)) {
                graph[e].orientedReadIds.push_back(p0);
            }
        }
    }


    void removeLowCoverageEdges(uint64_t minCoverage)
    {
        Graph& graph = *this;
        vector<edge_descriptor> edgesToBeRemoved;
        BGL_FORALL_EDGES_T(e, graph, Graph) {
            if(graph[e].orientedReadIds.size() < minCoverage) {
                edgesToBeRemoved.push_back(e);
            }
        }

        for(const edge_descriptor e: edgesToBeRemoved) {
            remove_edge(e, graph);
        }
    }



    void writeGraphviz() const
    {
        const Graph& graph = *this;

        ofstream out("DeBruijnGraph.dot");
        out << "digraph DeBruijnGraph {\n";
        BGL_FORALL_EDGES_T(e, graph, Graph) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            const auto coverage = graph[e].orientedReadIds.size();
            graph[v0].writeGraphviz(out);
            out << "->";
            graph[v1].writeGraphviz(out);
            out << "[";
            out << "penwidth=\"" << sqrt(double(coverage)) << "\"";
            out << " tooltip=\"(";
            graph[v0].write(out);
            out << ")->(";
            graph[v1].write(out);
            out << ") ";
            out << coverage << "\"";
            out << "]";
            out << ";\n";
        }
        out << "}\n";

    }
};



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



    // Parameters that control the process below. EXPOSE WHEN CODE STABILIZES. *********

    // The minimum number of markers for a segment to be used.
    const uint64_t minMarkerCount = 0;

    const uint64_t minEdgeCoverage = 2;

#if 0
    // The minimum length of a pseudo-path for a read t be used.
    const uint64_t minPseudoPathLength = 3;

    // The minimum number of aligned meta-markers for an alignment to be used.
    const uint64_t minAlignedMetaMarkerCount = 3;
#endif



    // Only consider segments that are sufficiently long.
    vector<bool> useSegment(segmentCount);
    for(SegmentId segmentId=0; segmentId<segmentCount; segmentId++) {
        useSegment[segmentId] = assemblyGraph.edgeLists.size(segmentId) >= minMarkerCount;
    }



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

                // If this segment is not being used, skip.
                if(not useSegment[segmentId]) {
                    continue;
                }

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



    // Use these pseudo-paths to create a de Bruijn graph.
    DeBruijnGraph<3> graph;
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const vector<SegmentId>& pseudoPath =
                pseudoPaths[orientedReadId.getValue()];
            graph.addOrientedReadVertices(orientedReadId, pseudoPath);
        }
    }
    graph.createEdges();
    graph.removeLowCoverageEdges(minEdgeCoverage);
    graph.writeGraphviz();
    cout << "The De Bruijn graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;


#if 0
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



    // Gather all pairs of oriented reads that occur in non-branching segments that are in use.
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
            if(pseudoPaths[orientedReadId0.getValue()].size() < minPseudoPathLength) {
                continue;
            }
            for(uint64_t i1=i0+1; i1<v.size(); i1++) {
                const OrientedReadId orientedReadId1 = v[i1].first;
                if(pseudoPaths[orientedReadId1.getValue()].size() < minPseudoPathLength) {
                    continue;
                }
                orientedReadPairs.push_back(make_pair(orientedReadId0, orientedReadId1));
            }
        }
    }
    deduplicate(orientedReadPairs);
    cout << "Found " << orientedReadPairs.size() <<
        " oriented read pairs." << endl;



    // The following process is similar to the one used to create the marker graph.
    // The difference is that, when creating the marker graph, each read is
    // represented by sequence of markers. But here, each read is represented
    // by the sequence of segments it encounters - its pseudo-path.

    // We use a disjoint set data structure where each entry represents one
    // segment of a pseudo-path.

    // Vector to contain, for each oriented read, the starting point of
    // its pseudo-path in the disjoint set data structure.
    vector<uint64_t> start(2*readCount(), std::numeric_limits<uint64_t>::max());
    uint64_t n = 0;
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const vector<SegmentId>& pseudoPath =
                pseudoPaths[orientedReadId.getValue()];
            if(pseudoPath.size() < minPseudoPathLength) {
                continue;
            }
            start[orientedReadId.getValue()] = n;
            n += pseudoPath.size();
        }
    }



    // Initialize the disjoint set data structure.
    vector<SegmentId> rank(n);
    vector<SegmentId> parent(n);
    boost::disjoint_sets<SegmentId*, SegmentId*> disjointSets(&rank[0], &parent[0]);
    for(SegmentId i=0; i<n; i++) {
        disjointSets.make_set(i);
    }
    cout << "The disjoint set data structure has size " << n << endl;


    // For each such pair of oriented reads, compute an alignment between their pseudo-paths.
    const bool writeAlignments = false;
    ofstream alignmentsCsv;
    if(writeAlignments) {
        alignmentsCsv.open("Alignments.csv");
    }
    uint64_t alignmentsDone = 0;
    for(const auto& p: orientedReadPairs) {
        if((alignmentsDone % 10000) == 0) {
            cout << alignmentsDone << "/" << orientedReadPairs.size() << endl;
        }
        ++alignmentsDone;
        const OrientedReadId orientedReadId0 = p.first;
        const OrientedReadId orientedReadId1 = p.second;

        // Get the pseudo-paths of these two oriented reads.
        const vector<AssemblyGraph::EdgeId>& pseudoPath0 = pseudoPaths[orientedReadId0.getValue()];
        const vector<AssemblyGraph::EdgeId>& pseudoPath1 = pseudoPaths[orientedReadId1.getValue()];
        SHASTA_ASSERT(pseudoPath0.size() >= minPseudoPathLength);
        SHASTA_ASSERT(pseudoPath1.size() >= minPseudoPathLength);

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
        if(writeAlignments) {
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



        // Find pairs of disjoint sets to be merged, based on this alignment.
        vector< pair<uint64_t, uint64_t> > toBeMerged;
        const uint64_t start0 = start[orientedReadId0.getValue()];
        const uint64_t start1 = start[orientedReadId1.getValue()];
        uint64_t i0 = 0;
        uint64_t i1 = 0;
        for(uint64_t j=0; j<alignmentLength; j++) {
            const uint64_t value0 = align[j];
            const uint64_t value1 = align[j+alignmentLength];
            // cout << j << " " << i0 << " " << i1 << " " << value0 << " " << value1 << endl;
            if(value0!=45 and value1!=45 and value0 == value1) {
                SHASTA_ASSERT(value0 == pseudoPath0[i0] + 100);
                SHASTA_ASSERT(value1 == pseudoPath1[i1] + 100);
                toBeMerged.push_back(make_pair(start0+i0, start1+i1));
            }
            if(value0 != 45) {
                ++i0;
            }
            if(value1 != 45) {
                ++i1;
            }
        }
        SHASTA_ASSERT(i0 == pseudoPath0.size());
        SHASTA_ASSERT(i1 == pseudoPath1.size());

        // If the alignment is too short, skip it.
        if(toBeMerged.size() < minAlignedMetaMarkerCount) {
            continue;
        }



        // In the disjoint set data structure, merge entries corresponding to
        // aligned segments.
        for(const auto& p: toBeMerged) {
            disjointSets.union_set(p.first, p.second);
        }
    }



    // Each of the disjoint data sets becomes a vertex of the MetaMarkerGraph.
    // Store the disjoint set for each position of the pseudo-path of each oriented read.
    vector< vector<uint64_t> > vertexTable(2*readCount());
    vector< vector< pair<OrientedReadId, uint64_t > > > vertices(n);
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const vector<SegmentId>& pseudoPath =
                pseudoPaths[orientedReadId.getValue()];
            if(pseudoPath.size() < minPseudoPathLength) {
                continue;
            }
            vector<uint64_t>& v = vertexTable[orientedReadId.getValue()];
            const uint64_t firstMetaMarkerId = start[orientedReadId.getValue()];

            for(uint64_t index=0; index<pseudoPath.size(); index++) {
                const uint64_t metaMarkerId = firstMetaMarkerId + index;
                const uint64_t disjointSetId = disjointSets.find_set(metaMarkerId);
                v.push_back(disjointSetId);
                vertices[disjointSetId].push_back(make_pair(orientedReadId, index));
            }
        }
    }


    vector<uint64_t> histogram;
    for(uint64_t i=0; i<n; i++) {
        const uint64_t size = vertices[i].size();
        if(size) {
            if(histogram.size() <= size) {
                histogram.resize(size+1, 0);
            }
            ++histogram[size];
        }
    }
    ofstream histogramCsv("Histogram.csv");
    histogramCsv << "Size,Frequency\n";
    for(size_t i=1; i<histogram.size(); i++) {
        const uint64_t frequency = histogram[i];
        if(frequency) {
            histogramCsv << i << "," << frequency << "\n";
        }
    }



    // Create vertices of the MetaMarkerGraph.
    MetaMarkerGraph graph;
    uint64_t vertexId = 0;
    for(uint64_t i=0; i<n; i++) {
        const vector< pair<OrientedReadId, uint64_t > >& v = vertices[i];
        if(v.empty()) {
            continue;
        }

        // Sanity check: they must all correspond to the same SegmentId.
        SegmentId firstSegmentId = std::numeric_limits<SegmentId>::max();
        for(uint64_t i=0; i<v.size(); i++) {
            const auto& p = v[i];
            const OrientedReadId orientedReadId = p.first;
            const uint64_t metaOrdinal = p.second;
            const vector<SegmentId>& pseudoPath = pseudoPaths[orientedReadId.getValue()];
            const SegmentId segmentId = pseudoPath[metaOrdinal];
            if(i == 0) {
                firstSegmentId = segmentId;
            } else {
                SHASTA_ASSERT(firstSegmentId == segmentId);
            }
        }

        // Create a vertex with these oriented reads and meta ordinals.
        add_vertex(MetaMarkerGraphVertex(
            vertexId++,
            firstSegmentId,
            assemblyGraph.edgeLists.size(firstSegmentId),
            v),
            graph);
    }
    graph.createEdges();
    graph.writeGfa("MetaMarkerGraph.gfa");
    graph.writeVerticesCsv("MetaMarkerGraphVertices.csv");
    graph.writeEdgesCsv("MetaMarkerGraphEdges.csv");
    cout << "The MetaMarkerGraph has " << num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges." << endl;
#endif



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
