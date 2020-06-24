// Shasta.
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "MarkerGraph2.hpp"
#include "MetaMarkerGraph.hpp"
#include "orderPairs.hpp"
#include "seqan.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Seqan.
#include <seqan/align.h>

// Standard library.
#include "array.hpp"
#include "fstream.hpp"
#include <map>
#include <set>



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
        const Assembler::PseudoPath& pseudoPath)
    {
        Graph& graph = *this;

        // Loop over possible starting positions such that the k segments
        // starting there are on the pseudo-path.
        for(uint64_t startPosition=0; startPosition+k<=pseudoPath.size(); startPosition++) {

            // Extract the k segments from the pseudo-path
            array<uint64_t, k> segmentIds;
            for(uint64_t i=0; i<k; i++) {
                segmentIds[i] = pseudoPath[startPosition+i].segmentId;
            }

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
        out << "digraph DeBruijnGraph {\n"
            "node[shape=rectangle];\n";
        BGL_FORALL_EDGES_T(e, graph, Graph) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            const auto coverage = graph[e].orientedReadIds.size();
            graph[v0].writeGraphviz(out);
            out << "->";
            graph[v1].writeGraphviz(out);
            out << "[";
            out << "penwidth=\"" << sqrt(double(coverage + 1)) << "\"";
            out << " tooltip=\"(";
            graph[v0].write(out);
            out << ")->(";
            graph[v1].write(out);
            out << ") ";
            out << coverage << "\"";
            if(coverage == 0) {
                out << " style=dotted";
            }
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

    // The minimum length of a pseudo-path for a read to be used.
    const uint64_t minPseudoPathLength = 3;

    // The minimum number of aligned meta-markers for an alignment to be used.
    const uint64_t minAlignedMetaMarkerCount = 3;

    // Alignment scores.
    const int matchScore = 1;
    const int mismatchScore = -1;
    const int gapScore = -2;

    // The minimum score for an alignment to be used.
    const int minScore = 3;



    // Compute the pseudo-path of each oriented read.
    // This vector is indexed by OrientedReadId::getValue().
    vector<MarkerGraph::EdgeId> path;
    vector< pair<uint32_t, uint32_t> > pathOrdinals;
    vector<PseudoPath> pseudoPaths(2*readCount());
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            computePseudoPath(orientedReadId, path, pathOrdinals,
                pseudoPaths[orientedReadId.getValue()]);
        }
    }



    // Write a csv file with the pseudo-path of each oriented read.
    {
        ofstream csv("PseudoPaths.csv");
        for(ReadId readId=0; readId<readCount(); readId++) {
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);
                csv << orientedReadId << ",";

                const PseudoPath& pseudoPath = pseudoPaths[orientedReadId.getValue()];
                for(const auto& pseudoPathEntry: pseudoPath) {
                    csv << pseudoPathEntry.segmentId << ",";
                }
                csv << "\n";
            }
        }
    }



    // Create the pseudo-path table which contains, for each segment,
    // its occurrences in oriented read pseudo-paths.
    // For each segmentId, we store a vector of pairs (orientedReadId, ordinal) such that
    // pseudoPaths[orientedReadId.getValue()][ordinal] == segmentId
    vector< vector< pair<OrientedReadId, uint64_t> > >  pseudoPathTable(segmentCount);
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const PseudoPath& pseudoPath = pseudoPaths[orientedReadId.getValue()];

            for(uint64_t ordinal=0; ordinal<pseudoPath.size(); ordinal++) {
                const SegmentId segmentId = pseudoPath[ordinal].segmentId;
                pseudoPathTable[segmentId].push_back(make_pair(orientedReadId, ordinal));
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
    // A segment (assembly graph edge) v0->v1 is non-branching
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

                // Store the pair with the lowest oriented read id first.
                if(orientedReadId0 < orientedReadId1) {
                    orientedReadPairs.push_back(make_pair(orientedReadId0, orientedReadId1));
                } else {
                    orientedReadPairs.push_back(make_pair(orientedReadId1, orientedReadId0));
                }
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
            const PseudoPath& pseudoPath = pseudoPaths[orientedReadId.getValue()];
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
        const PseudoPath& pseudoPath0 = pseudoPaths[orientedReadId0.getValue()];
        const PseudoPath& pseudoPath1 = pseudoPaths[orientedReadId1.getValue()];
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
        for(const auto& pseudoPathEntry: pseudoPath0) {
            appendValue(seq0, pseudoPathEntry.segmentId + 100);
        }
        TSequence seq1;
        for(const auto& pseudoPathEntry: pseudoPath1) {
            appendValue(seq1, pseudoPathEntry.segmentId + 100);
        }

        // Store them in a SeqAn string set.
        TStringSet sequences;
        appendValue(sequences, seq0);
        appendValue(sequences, seq1);

        // Compute the alignment.
        TAlignGraph graph(sequences);
        const int score = globalAlignment(
                graph,
                Score<int, Simple>(matchScore, mismatchScore, gapScore),
                AlignConfig<true, true, true, true>(),
                LinearGaps());
        if(score < minScore) {
            continue;
        }

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

        // If the alignment contains any mismatches, discard it.
        uint64_t i0 = 0;
        uint64_t i1 = 0;
        bool mismatchFound =  false;
        for(uint64_t j=0; j<alignmentLength; j++) {
            const uint64_t value0 = align[j];
            const uint64_t value1 = align[j+alignmentLength];
            // cout << j << " " << i0 << " " << i1 << " " << value0 << " " << value1 << endl;
            if(value0!=45 and value1!=45 and value0 != value1) {
                mismatchFound = true;
                break;
            }
            if(value0 != 45) {
                ++i0;
            }
            if(value1 != 45) {
                ++i1;
            }
        }
        if(mismatchFound) {
            continue;
        }
        SHASTA_ASSERT(i0 == pseudoPath0.size());
        SHASTA_ASSERT(i1 == pseudoPath1.size());



        // Find pairs of disjoint sets to be merged, based on this alignment.
        vector< pair<uint64_t, uint64_t> > toBeMerged;
        const uint64_t start0 = start[orientedReadId0.getValue()];
        const uint64_t start1 = start[orientedReadId1.getValue()];
        i0 = 0;
        i1 = 0;
        for(uint64_t j=0; j<alignmentLength; j++) {
            const uint64_t value0 = align[j];
            const uint64_t value1 = align[j+alignmentLength];
            // cout << j << " " << i0 << " " << i1 << " " << value0 << " " << value1 << endl;
            if(value0!=45 and value1!=45 and value0 == value1) {
                SHASTA_ASSERT(value0 == pseudoPath0[i0].segmentId + 100);
                SHASTA_ASSERT(value1 == pseudoPath1[i1].segmentId + 100);
                toBeMerged.push_back(make_pair(start0 + i0, start1 + i1));
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
            const PseudoPath& pseudoPath = pseudoPaths[orientedReadId.getValue()];
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
            const PseudoPath& pseudoPath = pseudoPaths[orientedReadId.getValue()];
            const SegmentId segmentId = pseudoPath[metaOrdinal].segmentId;
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
    graph.writeGraphviz("MetaMarkerGraph.dot");
    graph.writeGfa("MetaMarkerGraph.gfa");
    graph.writeVerticesCsv("MetaMarkerGraphVertices.csv");
    graph.writeEdgesCsv("MetaMarkerGraphEdges.csv");
    cout << "The MetaMarkerGraph has " << num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges." << endl;
}



// Analyze paths of oriented reads that go through a given assembly graph edge (segment).
void Assembler::analyzeOrientedReadPathsThroughSegment(
    AssemblyGraph::EdgeId startSegmentId)
{
    vector<AssemblyGraph::EdgeId> forwardChokePoints;
    vector<AssemblyGraph::EdgeId> backwardChokePoints;
    const bool debug = true;
    analyzeOrientedReadPathsThroughSegment(
        startSegmentId,
        forwardChokePoints,
        backwardChokePoints,
        debug
        );
}



void Assembler::findOrientedReadsOnAssemblyGraphEdge(
    AssemblyGraph::EdgeId segmentId,
    vector<OrientedReadId>& orientedReadIds
) const
{
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    std::set<OrientedReadId> orientedReadIdsSet;

    // Find the marker graph edges that are on this assembly graph edge (segment).
    const span< const MarkerGraph::EdgeId> markerGraphEdges =
        assemblyGraph.edgeLists[segmentId];

    // Loop over these marker graph edges.
    for(const MarkerGraph::EdgeId markerGraphEdgeId: markerGraphEdges) {

        // Loop over oriented read ids on this marker graph edge.
        for(const MarkerInterval& interval: markerGraph.edgeMarkerIntervals[markerGraphEdgeId]) {
            orientedReadIdsSet.insert(interval.orientedReadId);
        }
    }

    orientedReadIds.clear();
    copy(orientedReadIdsSet.begin(), orientedReadIdsSet.end(),
        back_inserter(orientedReadIds));
}



void Assembler::followOrientedReadPaths(
    AssemblyGraph::EdgeId startSegmentId,
    bool forward)
{
    const bool debug = true;
    using SegmentId = AssemblyGraph::EdgeId;

    // PARAMETERS TO EXPOSE WHEN THE CODE STABILIZES.
    const int64_t matchScore = 1;
    const int64_t mismatchScore = -1;
    const int64_t gapScore = -2;
    const double scoreFraction = 0.5;

    cout << "Starting at segment " << startSegmentId <<
        " and moving " << (forward ? "forward" : "backward") << "." << endl;

    // Start with consensus equal to just our start segment.
    vector<SegmentId> consensus = {startSegmentId};

    // The position in the consensus of the segment we are currently working on.
    uint64_t consensusPosition = 0;

    // A map to contain oriented reads we already encountered and their pseudo-paths.
    std::map<OrientedReadId, PseudoPath> orientedReadIdsAlreadyFound;



    // Main iteration loop.
    for(int iteration=0; ; ++iteration) {
        cout << "Begin iteration " << iteration <<
            " with current consensus position " << consensusPosition << endl;

        // The segment we process at this iteration.
        const SegmentId segmentId0 = consensus[consensusPosition];

        // Find all oriented reads on segmentId0.
        vector<OrientedReadId> orientedReadIds0;
        findOrientedReadsOnAssemblyGraphEdge(segmentId0, orientedReadIds0);



        // Of the oriented reads on segmentId0, find the ones
        // we did not already encounter.
        // For each one compute and store the pseudo-path.
        // The newly found oriented reads are stored with the length
        // of their pseudo-path.
        vector< pair<OrientedReadId, uint64_t> > newOrientedReadIds0;
        for(const OrientedReadId orientedReadId: orientedReadIds0) {

            if(orientedReadIdsAlreadyFound.find(orientedReadId) !=
                orientedReadIdsAlreadyFound.end()) {
                continue;
            }

            vector<MarkerGraph::EdgeId> path;
            vector< pair<uint32_t, uint32_t> > pathOrdinals;
            PseudoPath pseudoPath;
            computePseudoPath(orientedReadId, path, pathOrdinals, pseudoPath);

            orientedReadIdsAlreadyFound.insert(make_pair(orientedReadId, pseudoPath));
            newOrientedReadIds0.push_back(make_pair(orientedReadId, pseudoPath.size()));
        }

        // Order them by decreasing pseudo-path length.
        sort(newOrientedReadIds0.begin(), newOrientedReadIds0.end(),
            OrderPairsBySecondOnlyGreater<OrientedReadId, uint64_t>());

        if(debug) {
            cout << "Found " << orientedReadIds0.size() <<
                " oriented reads on segment " << segmentId0 <<
                " of which " << newOrientedReadIds0.size() <<
                " are new." << endl;
        }



        // Loop over the new oriented reads in order of decreasing pseudo-path length.
        for(const auto& p: newOrientedReadIds0) {
            const OrientedReadId orientedReadId = p.first;
            if(debug) {
                cout << "Current consensus: ";
                copy(consensus.begin(), consensus.end(), ostream_iterator<SegmentId>(cout, " "));
                cout << "\nCurrent consensus position: " << consensusPosition << endl;
                cout << "Processing oriented read " << orientedReadId << endl;
            }
            const PseudoPath& pseudoPath = orientedReadIdsAlreadyFound[orientedReadId];
            if(debug) {
                cout << orientedReadId << " pseudo-path: ";
                for(auto& pseudoPathEntry: pseudoPath) {
                    cout << " " << pseudoPathEntry.segmentId;
                }
                cout << endl;
            }

            // Align this pseudo-path with the current consensus.
            vector<SegmentId> pseudoPathSegments;
            for(auto& pseudoPathEntry: pseudoPath) {
                pseudoPathSegments.push_back(pseudoPathEntry.segmentId);
            }
            if(consensus.size()==1) {
                consensus = pseudoPathSegments;
                consensusPosition = std::find(consensus.begin(), consensus.end(), segmentId0) -
                    consensus.begin();
            } else {
                vector< pair<bool, bool> > alignment;
                const int64_t score = seqanAlign(
                    consensus.begin(), consensus.end(),
                    pseudoPathSegments.begin(), pseudoPathSegments.end(),
                    matchScore, mismatchScore, gapScore,
                    true, true,
                    alignment);

                if(debug) {
                    cout << "Alignment score " << score << endl;
                    for(const auto& p: alignment) {
                        cout << (p.first ? '.' : '-');
                    }
                    cout << endl;
                    for(const auto& p: alignment) {
                        cout << (p.second ? '.' : '-');
                    }
                    cout << endl;
                }


                // If the alignment score is too low, ignore it.
                // This check will need some refinement.
                if(score < int64_t(double(pseudoPathSegments.size())*scoreFraction)) {
                    if(debug) {
                        cout << "Ignored due to low alignment score." << endl;
                    }
                    continue;
                }

                // If the alignment has any mismatches, ignore it.
                uint64_t position0 = 0;
                uint64_t position1 = 0;
                bool mismatchFound = false;
                for(const auto& p: alignment) {
                    if(p.first and p.second) {
                        if(consensus[position0] != pseudoPathSegments[position1]) {
                            mismatchFound = true;
                            break;
                        }
                        ++position0;
                        ++position1;
                    } else if(p.first) {
                        ++position0;
                    } else if(p.second) {
                        ++position1;
                    }
                }
                if(mismatchFound) {
                    if(debug) {
                        cout << "Ignored due to mismatch." << endl;
                    }
                    continue;
                }
                SHASTA_ASSERT(position0 == consensus.size());
                SHASTA_ASSERT(position1 == pseudoPathSegments.size());



                // If getting here, we know that the alignment has only
                // matches and gaps - no mismatches.
                // Update the consensus.
                vector<SegmentId> newConsensus;
                position0 = 0;
                position1 = 0;
                for(const auto& p: alignment) {
                    if(p.first and p.second) {
                        SHASTA_ASSERT(consensus[position0] == pseudoPathSegments[position1]);
                        if(position0 == consensusPosition) {
                            consensusPosition = newConsensus.size();
                        }
                        newConsensus.push_back(consensus[position0]);
                        ++position0;
                        ++position1;
                    } else if(p.first) {
                        if(position0 == consensusPosition) {
                            consensusPosition = newConsensus.size();
                        }
                        newConsensus.push_back(consensus[position0]);
                        ++position0;
                    } else if(p.second) {
                        newConsensus.push_back(pseudoPathSegments[position1]);
                        ++position1;
                    }
                }
                SHASTA_ASSERT(position0 == consensus.size());
                SHASTA_ASSERT(position1 == pseudoPathSegments.size());
                if(debug) {
                    if(newConsensus == consensus) {
                        cout << "No change in consensus." << endl;
                    } else {
                        cout << "Consensus updated." << endl;
                    }
                    cout << "Current consensus position is now " << consensusPosition << endl;
                }
                consensus = newConsensus;
            }

        }

        // Update position in consensus for the next iteration.
        if(forward) {
            ++consensusPosition;
            if(consensusPosition == consensus.size()) {
                break;
            }
        } else {
            if(consensusPosition == 0) {
                break;
            }
            --consensusPosition;
        }
        if(debug) {
            cout << "Consensus position updated to " << consensusPosition <<
                " for next iteration." << endl;
        }
    }
}



#if 0
void Assembler::followOrientedReadPaths(
    AssemblyGraph::EdgeId startSegmentId,
    bool forward)
{
    const bool debug = true;
    using SegmentId = AssemblyGraph::EdgeId;

    cout << "Starting at segment " << startSegmentId <<
        " and moving " << (forward ? "forward" : "backward") << "." << endl;

    vector<SegmentId> forwardChokePoints;
    vector<SegmentId> backwardChokePoints;
    analyzeOrientedReadPathsThroughSegment(
        startSegmentId,
        forwardChokePoints,
        backwardChokePoints,
        false
        );


    // Begin with consensus equal to the choke points of the start segment.
    vector<SegmentId> consensus = {startSegmentId};
    if(forward) {
        copy(forwardChokePoints.begin(), forwardChokePoints.end(),
            back_inserter(consensus));
    } else {
        copy(backwardChokePoints.begin(), backwardChokePoints.end(),
            back_inserter(consensus));
    }



    // Main iteration loop.
    // At each iteration we move forward (or backward) by one segment.
    // We compute the choke points of the segment at consensus[index].
    uint64_t index = 1;
    while(index < consensus.size()) {

        if(debug) {
            cout << "Consensus at beginning of iteration " << index-1 << ":" << endl;
            copy(consensus.begin(), consensus.end(), ostream_iterator<SegmentId>(cout, " "));
            cout << endl;
        }


        // Compute choke points for the current segment in the consensus.
        const SegmentId currentSegmentId = consensus[index];
        analyzeOrientedReadPathsThroughSegment(
            currentSegmentId,
            forwardChokePoints,
            backwardChokePoints,
            false
            );
        const vector<SegmentId>& chokePoints = (forward ? forwardChokePoints : backwardChokePoints);

        if(debug) {
            cout << "Choke points of segment " << currentSegmentId << ":" << endl;
            copy(chokePoints.begin(), chokePoints.end(), ostream_iterator<SegmentId>(cout, " "));
            cout << endl;
        }

        // If there are no choke points, don't do anything.
        if(chokePoints.empty()) {
            if(debug) {
                cout << "Alignment skipped because there are no choke points." << endl;
            }
            ++index;
            continue;
        }

        // If this was the last segment in the consensus, just append the choke points
        // to the consensus.
        if(index == consensus.size() - 1) {
            copy(chokePoints.begin(), chokePoints.end(), back_inserter(consensus));
            if(debug) {
                cout << "Choke points appended at the end of consensus." << endl;
            }
            ++index;
            continue;
        }



        // Use SeqAn to compute an alignment of these choke points with the portion
        // of the current consensus following the current segment.
        const auto& constConsensus = consensus;
        vector< pair<bool, bool> > alignment;
        const int64_t matchScore = 1;
        const int64_t mismatchScore = -1;
        const int64_t gapScore = -1;
        const int64_t score = seqanAlign(
            constConsensus.begin() + index + 1, constConsensus.end(),
            chokePoints.begin(), chokePoints.end(),
            matchScore, mismatchScore, gapScore,
            false, true,
            alignment);
        cout << "Alignment score is " << score << endl;
        for(uint64_t i=0; i<alignment.size(); i++) {
            cout << (alignment[i].first ? '.' : '-');
        }
        cout << endl;
        for(uint64_t i=0; i<alignment.size(); i++) {
            cout << (alignment[i].second ? '.' : '-');
        }
        cout << endl;



        // Use the alignment to update the consensus.
        vector<SegmentId> newConsensus;
        copy(consensus.begin(), consensus.begin() + index + 1, back_inserter(newConsensus));
        uint64_t consensusPosition = newConsensus.size();
        uint64_t chokePointsPosition = 0;
        for(const auto& p: alignment) {
            if(p.first and p.second) {
                const SegmentId oldConsensusValue = consensus[consensusPosition];
                SHASTA_ASSERT(oldConsensusValue == chokePoints[chokePointsPosition]);
                newConsensus.push_back(oldConsensusValue);
                ++consensusPosition;
                ++chokePointsPosition;
            } else if(p.first) {
                const SegmentId oldConsensusValue = consensus[consensusPosition];
                newConsensus.push_back(oldConsensusValue);
                ++consensusPosition;
            } else if(p.second) {
                newConsensus.push_back(chokePoints[chokePointsPosition]);
                ++chokePointsPosition;
            } else {
                SHASTA_ASSERT(0);
            }
         }
        consensus.swap(newConsensus);



        // The next iteration will process the next segment in the current consensus,
        // unless we reached the end of the consensus.
        ++index;
    }


    if(debug) {
        cout << "Final consensus:" << endl;
        copy(consensus.begin(), consensus.end(), ostream_iterator<SegmentId>(cout, " "));
        cout << endl;
    }
}
#endif



void Assembler::analyzeOrientedReadPathsThroughSegment(
    AssemblyGraph::EdgeId startSegmentId,
    vector<AssemblyGraph::EdgeId>& forwardChokePoints,
    vector<AssemblyGraph::EdgeId>& backwardChokePoints,
    bool debug)
{

    using SegmentId = AssemblyGraph::EdgeId;
    AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    const uint64_t segmentCount = assemblyGraph.edges.size();
    SHASTA_ASSERT(startSegmentId < segmentCount);



    // ********************************* PARAMETERS TO EXPOSE WHEN CODE STABILIZES.

    // The maximum ordinal skip between adjacent segments on the
    // pseudo-path of an oriented read. If an oriented read has a skip larger than
    // that, it is not used.
    const uint32_t maxOrdinalSkip = 500;
    const uint64_t minEdgeCount = 2;
    const int matchScore = 1;
    const int mismatchScore = -1;
    const int gapScore = -1;



    // Find the oriented reads that have edges on this assembly graph edge (segment).
    std::map<OrientedReadId, uint64_t> orientedReadIdsThroughSegment;
    // Loop over the marker graph edges that are on this assembly graph edge (segment).
    const span<MarkerGraph::EdgeId> markerGraphEdges =
        assemblyGraph.edgeLists[startSegmentId];
    for(const MarkerGraph::EdgeId markerGraphEdgeId: markerGraphEdges) {

        // Loop over oriented read ids on this marker graph edge.
        for(const MarkerInterval& interval: markerGraph.edgeMarkerIntervals[markerGraphEdgeId]) {
            ++orientedReadIdsThroughSegment[interval.orientedReadId];
        }
    }
    if(debug) {
        cout << "Found " << orientedReadIdsThroughSegment.size() <<
            " oriented reads on segment " << startSegmentId << endl;
    }



    // Compute pseudo-paths for these oriented reads.
    // Only keep the oriented reads without a large marker skip.
    vector<OrientedReadId> orientedReadIds;
    vector<PseudoPath> pseudoPaths;
    vector<MarkerGraph::EdgeId> markerGraphPath;
    vector< pair<uint32_t, uint32_t> > pathOrdinals;
    for(const auto& p: orientedReadIdsThroughSegment) {
        if(p.second < minEdgeCount) {
            continue;
        }
        const OrientedReadId orientedReadId = p.first;

        // Compute the pseudo-path.
        PseudoPath pseudoPath;
        computePseudoPath(orientedReadId, markerGraphPath, pathOrdinals, pseudoPath);



        // If the pseudo-path has a large ordinal skip,
        // disregard this oriented read.
        bool disregard = false;
        for(uint64_t j=1; j<pseudoPath.size(); j++) {
            const uint32_t ordinalSkip = pseudoPath[j].firstOrdinal - pseudoPath[j-1].lastOrdinal;
            if(ordinalSkip > maxOrdinalSkip) {
                disregard = true;
                break;
            }
        }
        // Also check the begin and end.
        if(not disregard) {
            if(pseudoPath.front().firstOrdinal > maxOrdinalSkip) {
                disregard = true;
            }
        }
        if(not disregard) {
            if(markers.size(orientedReadId.getValue()) - pseudoPath.back().lastOrdinal > maxOrdinalSkip) {
                disregard = true;
            }
        }
        if(disregard) {
            if(debug) {
                cout << orientedReadId << " disregarded because of a large ordinal skip." << endl;
            }
            continue;
        }



        // Store this oriented read and its pseudo-path.
        orientedReadIds.push_back(orientedReadId);
        pseudoPaths.push_back(pseudoPath);
    }
    if(debug) {
        cout << "Continuing with " << orientedReadIds.size() <<
            " oriented reads." << endl;
    }



    // The following process is similar to the one used to create the marker graph.
    // The difference is that, when creating the marker graph, each oriented read is
    // represented by sequence of markers. But here, each oriented read is represented
    // by the sequence of segments it encounters - its pseudo-path.

    // We use a disjoint set data structure where each entry represents one
    // segment of a pseudo-path.

    // Vector to contain, for each oriented read, the starting point of
    // its pseudo-path in the disjoint set data structure.
    vector<uint64_t> start(orientedReadIds.size(), std::numeric_limits<uint64_t>::max());
    uint64_t n = 0;
    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        start[i] = n;
        n += pseudoPaths[i].size();
    }



    // Initialize the disjoint set data structure.
    vector<SegmentId> rank(n);
    vector<SegmentId> parent(n);
    boost::disjoint_sets<SegmentId*, SegmentId*> disjointSets(&rank[0], &parent[0]);
    for(SegmentId segmentId=0; segmentId<n; segmentId++) {
        disjointSets.make_set(segmentId);
    }
    if(debug) {
        cout << "The disjoint set data structure has size " << n << endl;
    }



    // For each pair of the oriented reads we kept, compute an alignment
    // between their pseudo-paths.
    const bool writeAlignments = debug;
    ofstream alignmentsCsv;
    if(writeAlignments) {
        alignmentsCsv.open("Alignments.csv");
    }
    for(uint64_t i0=0; i0<orientedReadIds.size(); i0++) {
        const OrientedReadId orientedReadId0 = orientedReadIds[i0];
        const PseudoPath& pseudoPath0 = pseudoPaths[i0];
        for(uint64_t i1=i0+1; i1<orientedReadIds.size(); i1++) {
            const OrientedReadId orientedReadId1 = orientedReadIds[i1];
            const PseudoPath& pseudoPath1 = pseudoPaths[i1];

            // Use SeqAn to compute an alignment free at both ends.
            // https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html
            using namespace seqan;

            // Hide shasta::Alignment.
            using seqan::Alignment;

            // An oriented read is represented by the segment ids in its pseudo-path.
            // We want to align a pair of such sequences.
            using TSequence = String<SegmentId>;

            // Other SeqAn types we need.
            using TStringSet = StringSet<TSequence>;
            using TDepStringSet = StringSet<TSequence, Dependent<> >;
            using TAlignGraph = Graph<Alignment<TDepStringSet> >;

            // Construct the sequences we want to pass to SeqAn.
            // Add 100 to all segment ids to avoid collision with the
            // value 45, used by SeqAn to represent gaps.
            TSequence seq0;
            for(const auto& pseudoPathEntry: pseudoPath0) {
                appendValue(seq0, pseudoPathEntry.segmentId+100);
            }
            TSequence seq1;
            for(const auto& pseudoPathEntry: pseudoPath1) {
                appendValue(seq1, pseudoPathEntry.segmentId+100);
            }

            // Store them in a SeqAn string set.
            TStringSet sequences;
            appendValue(sequences, seq0);
            appendValue(sequences, seq1);

            // Compute the alignment.
            TAlignGraph graph(sequences);
            const auto alignmentScore = globalAlignment(
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
                alignmentsCsv << "Score," << alignmentScore << "\n\n";
            }



            // Find pairs of disjoint sets to be merged, based on this alignment.
            vector< pair<uint64_t, uint64_t> > toBeMerged;
            const uint64_t start0 = start[i0];
            const uint64_t start1 = start[i1];
            uint64_t position0 = 0;
            uint64_t position1 = 0;
            for(uint64_t j=0; j<alignmentLength; j++) {
                const uint64_t value0 = align[j];
                const uint64_t value1 = align[j+alignmentLength];
                // cout << j << " " << i0 << " " << i1 << " " << value0 << " " << value1 << endl;
                if(value0!=45 and value1!=45 and value0 == value1) {
                    SHASTA_ASSERT(value0 == pseudoPath0[position0].segmentId + 100);
                    SHASTA_ASSERT(value1 == pseudoPath1[position1].segmentId + 100);
                    toBeMerged.push_back(make_pair(start0+position0, start1+position1));
                }
                if(value0 != 45) {
                    ++position0;
                }
                if(value1 != 45) {
                    ++position1;
                }
            }
            SHASTA_ASSERT(position0 == pseudoPath0.size());
            SHASTA_ASSERT(position1 == pseudoPath1.size());

            // In the disjoint set data structure, merge entries corresponding to
            // aligned segments.
            for(const auto& p: toBeMerged) {
                disjointSets.union_set(p.first, p.second);
            }
        }
    }



    // Each of the disjoint data sets becomes a vertex of the MetaMarkerGraph.
    // Store them in data structures that allow us to easily go from
    // vertices to oriented reads and vice versa.
    vector< vector<uint64_t> > vertexTable(orientedReadIds.size());
    vector< vector< pair<uint64_t, uint64_t > > > vertices(n);
    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const PseudoPath& pseudoPath = pseudoPaths[i];

        for(uint64_t position=0; position<pseudoPath.size(); position++) {
            const uint64_t metaMarkerId = start[i] + position;
            const uint64_t disjointSetId = disjointSets.find_set(metaMarkerId);
            vertexTable[i].push_back(disjointSetId);
            vertices[disjointSetId].push_back(make_pair(i, position));
        }
    }


    if(debug) {
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
    }



    // Create the MetaMarkerGraph, with one vertex for each of
    // the disjoint sets we found.
    MetaMarkerGraph graph;
    uint64_t vertexId = 0;
    for(uint64_t i=0; i<n; i++) {
        const vector< pair<uint64_t, uint64_t > >& v = vertices[i];
        if(not v.empty()) {
            const SegmentId segmentId = pseudoPaths[v.front().first][v.front().second].segmentId;
            vector< pair<OrientedReadId, uint64_t> > u;
            for(const auto& p: v) {
                u.push_back(make_pair(orientedReadIds[p.first], p.second));
            }
            add_vertex(MetaMarkerGraphVertex(
                vertexId++,
                segmentId,
                assemblyGraph.edgeLists.size(segmentId),
                u),
                graph);
        }
    }
    graph.createEdges();
    if(debug) {
        cout << "Before transitive reduction, the MetaMarkerGraph has " <<
            num_vertices(graph) << " vertices and " <<
            num_edges(graph) << " edges." << endl;
    }
    graph.transitiveReduction();

    if(debug) {
        graph.writeGraphviz("MetaMarkerGraph.dot", startSegmentId);
        graph.writeGfa("MetaMarkerGraph.gfa");
        graph.writeVerticesCsv("MetaMarkerGraphVertices.csv");
        graph.writeEdgesCsv("MetaMarkerGraphEdges.csv");
        cout << "The final MetaMarkerGraph has " << num_vertices(graph) << " vertices and " <<
            num_edges(graph) << " edges." << endl;
    }



#if 0
    // We want to construct the linear chain that includes our start segment.
    // First, locate the vertices corresponding to the start segment.
    // If there is more than one (possible but unusual), give up.
    vector<SegmentId> chain;
    graph.findLinearChain(startSegmentId, chain);
    cout << "Found a linear chain with " << chain.size() << " segments:" << endl;
    copy(chain.begin(), chain.end(), ostream_iterator<SegmentId>(cout, " "));
    cout << endl;
#endif


    // Find choke points.
    vector<MetaMarkerGraph::vertex_descriptor> forwardChokePointVertices;
    vector<MetaMarkerGraph::vertex_descriptor> backwardChokePointVertices;
    graph.findForwardChokePoints(startSegmentId, forwardChokePointVertices);
    graph.findBackwardChokePoints(startSegmentId, backwardChokePointVertices);
    forwardChokePoints.clear();
    backwardChokePoints.clear();
    for(const auto v: forwardChokePointVertices) {
        forwardChokePoints.push_back(graph[v].segmentId);
    }
    for(const auto v: backwardChokePointVertices) {
        backwardChokePoints.push_back(graph[v].segmentId);
    }

    if(debug) {
        cout << "Forward choke points: ";
        copy(forwardChokePoints.begin(), forwardChokePoints.end(),
            ostream_iterator<SegmentId>(cout, " "));
        cout << endl;
        cout << "Backward choke points: ";
        copy(backwardChokePoints.begin(), backwardChokePoints.end(),
            ostream_iterator<SegmentId>(cout, " "));
        cout << endl;
    }


#if 0
    // Write a fasta file with the sequence of the segments in the chain.
    if(debug) {
        ofstream fasta("Chain.fasta");
        fasta << ">Chain\n";
        for(uint64_t i=0; i<chain.size(); i++) {
            const SegmentId segmentId = chain[i];

            SHASTA_ASSERT(not assemblyGraph.edges[segmentId].wasRemoved());

            // Get the id of the reverse complemented edge.
            const SegmentId segmentIdRc = assemblyGraph.reverseComplementEdge[segmentId];


            // Write the sequence.
            if(assemblyGraph.isAssembledEdge(segmentId)) {

                // This segment was assembled. We can just write the stored sequence.
                auto sequence = assemblyGraph.sequences[segmentId];
                auto repeatCounts = assemblyGraph.repeatCounts[segmentId];
                SHASTA_ASSERT(sequence.baseCount == repeatCounts.size());
                for(size_t i=0; i<sequence.baseCount; i++) {
                    const Base b = sequence[i];
                    const uint8_t repeatCount = repeatCounts[i];
                    for(size_t k=0; k<repeatCount; k++) {
                        fasta << b;
                    }
                }
            } else {

                // This segment was not assembled. We write out the reverse
                // complemented sequence of the reverse complemented edge.
                SHASTA_ASSERT(assemblyGraph.isAssembledEdge(segmentIdRc));
                const auto sequence = assemblyGraph.sequences[segmentIdRc];
                const auto repeatCounts = assemblyGraph.repeatCounts[segmentIdRc];
                SHASTA_ASSERT(sequence.baseCount == repeatCounts.size());
                for(size_t i=0; i<sequence.baseCount; i++) {
                    const size_t j = sequence.baseCount - 1 - i;
                    const Base b = sequence[j].complement();
                    const uint8_t repeatCount = repeatCounts[j];
                    for(size_t k=0; k<repeatCount; k++) {
                        fasta << b;
                    }
                }
            }
        }
    }
#endif


#if 0
    // Now create a De Bruijn graph using these pseudo-paths.
    ofstream csv("PseudoPathsThroughSegment.csv");
    DeBruijnGraph<2> graph;
    for(uint64_t i=0; i<orientedReadIds.size(); i++) {
        const OrientedReadId orientedReadId = orientedReadIds[i];
        const PseudoPath& pseudoPath = pseudoPaths[i];

        // If the pseudo-path has a large ordinal skip,
        // disregard this read.
        bool disregard = false;
        for(uint64_t j=1; j<pseudoPath.size(); j++) {
            const uint32_t ordinalSkip = pseudoPath[j].firstOrdinal - pseudoPath[j-1].lastOrdinal;
            if(ordinalSkip > maxOrdinalSkip) {
                disregard = true;
                break;
            }
        }
        if(disregard) {
            cout << orientedReadId << " disregarded because of a large ordinal skip." << endl;
            continue;
        }

        // Write a line to the csv file.
        csv << orientedReadId << ",";
        for(const PseudoPathEntry& pseudoPathEntry: pseudoPath) {
            csv << pseudoPathEntry.segmentId << ",";
        }
        csv << "\n";

        // Update the De Bruijn graph for this oriented read.
        graph.addOrientedReadVertices(orientedReadId, pseudoPath);
    }
    graph.createEdges();
    graph.writeGraphviz();
    cout << "The De Bruijn graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;
#endif

}



void Assembler::computePseudoPath(
    OrientedReadId orientedReadId,

    // The marker graph path computed using computeOrientedReadMarkerGraphPath.
    // This is computed by this function - it does not neet to be filled in
    // in advance.
    vector<MarkerGraph::EdgeId>& path,
    vector< pair<uint32_t, uint32_t> >& pathOrdinals,

    // The pseudo-path computed by this function.
    PseudoPath& pseudoPath) const
{
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    using SegmentId = AssemblyGraph::EdgeId;

    // Compute the marker graph path.
    computeOrientedReadMarkerGraphPath(
        orientedReadId,
        0, uint32_t(markers.size(orientedReadId.getValue()) - 1),
        path, pathOrdinals);
    SHASTA_ASSERT(path.size() == pathOrdinals.size());



    // Now compute the pseudo-path.
    pseudoPath.clear();
    pseudoPath.clear();
    PseudoPathEntry pseudoPathEntry;
    pseudoPathEntry.segmentId = std::numeric_limits<SegmentId>::max();
    for(uint64_t i=0; i<path.size(); i++) {
        const MarkerGraph::EdgeId markerGraphEdgeId = path[i];
        const pair<uint32_t, uint32_t>& ordinals = pathOrdinals[i];

        // Get the corresponding assembly graph segments.
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
        const uint32_t positionInSegment = v.front().second;

        // Same as the previous.
        if(segmentId == pseudoPathEntry.segmentId) {
            pseudoPathEntry.lastOrdinal = ordinals.second;
            pseudoPathEntry.lastPosition = positionInSegment;
            ++pseudoPathEntry.markerGraphEdgeCount;
            continue;
        }

        // This is the next segment edge encountered
        // by this oriented read along its marker graph path.
        if(pseudoPathEntry.segmentId != std::numeric_limits<SegmentId>::max()) {
            pseudoPath.push_back(pseudoPathEntry);
        }
        pseudoPathEntry.segmentId = segmentId;
        pseudoPathEntry.firstOrdinal = ordinals.first;
        pseudoPathEntry.lastOrdinal = ordinals.second;
        pseudoPathEntry.firstPosition = positionInSegment;
        pseudoPathEntry.lastPosition = positionInSegment;
        pseudoPathEntry.markerGraphEdgeCount = 1;
    }

    // Add the last entry.
    if(pseudoPathEntry.segmentId != std::numeric_limits<SegmentId>::max()) {
        pseudoPath.push_back(pseudoPathEntry);
    }
}



// Write the pseudo-path of an oriented read to a csv file.
// The pseudo-path is the sequence of assembly graph edges
// (not necsssarily all adjacent, so not necessatily a path)
// encountered by the oriented read.
void Assembler::writePseudoPath(ReadId readId, Strand strand) const
{
    // Compute the pseudo path.
    const OrientedReadId orientedReadId(readId, strand);
    vector<MarkerGraph::EdgeId> path;
    vector< pair<uint32_t, uint32_t> > pathOrdinals;
    PseudoPath pseudoPath;
    computePseudoPath(orientedReadId, path, pathOrdinals, pseudoPath);

    // Write it out.
    ofstream csv("PseudoPath.csv");
    csv << "Segment id,First ordinal,Last ordinal,"
        "First position in segment,Last position in segment, Marker graph edge count\n";
    for(const PseudoPathEntry& pseudoPathEntry: pseudoPath) {
        csv << pseudoPathEntry.segmentId << ",";
        csv << pseudoPathEntry.firstOrdinal << ",";
        csv << pseudoPathEntry.lastOrdinal << ",";
        csv << pseudoPathEntry.firstPosition << ",";
        csv << pseudoPathEntry.lastPosition << ",";
        csv << pseudoPathEntry.markerGraphEdgeCount << "\n";
    }
}



# if 0
// This version uses choke points.
void Assembler::analyzeOrientedReadPaths()
{
    using SegmentId = AssemblyGraph::EdgeId;
    AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    const uint64_t segmentCount = assemblyGraph.edges.size();
    cout << "Found " << segmentCount << " segments." << endl;

    ofstream csv1("ChokePoints.csv");
    csv1 << "SegmentId,Direction\n";
    ofstream csv2("ChokePaths.csv");
    csv1 << "SegmentId\n";

    vector<SegmentId> forwardChokePoints;
    vector<SegmentId> backwardChokePoints;
    vector< vector<SegmentId> > paths(segmentCount);
    for(SegmentId segmentId=0; segmentId!=segmentCount; segmentId++) {
        analyzeOrientedReadPathsThroughSegment(
            segmentId,
            forwardChokePoints,
            backwardChokePoints,
            false);

        // Reverse the backward choke points so everything
        // is written in the forward direction.
        reverse(backwardChokePoints.begin(), backwardChokePoints.end());

        // Store the path.
        vector<SegmentId>& path = paths[segmentId];
        copy(backwardChokePoints.begin(), backwardChokePoints.end(), back_inserter(path));
        path.push_back(segmentId);
        copy(forwardChokePoints.begin(), forwardChokePoints.end(), back_inserter(path));

        // Write to ChokePoints.csv on two lines.
        csv1 << segmentId << ",Forward,";
        copy(forwardChokePoints.begin(), forwardChokePoints.end(),
            ostream_iterator<SegmentId>(csv1, ","));
        csv1 << "\n";
        csv1 << segmentId << ",Backward,";
        copy(backwardChokePoints.begin(), backwardChokePoints.end(),
            ostream_iterator<SegmentId>(csv1, ","));
        csv1 << "\n";

        // Write to ChokePaths.csv on one line.
        csv2 << segmentId << ",";
        copy(path.begin(), path.end(),
            ostream_iterator<SegmentId>(csv2, ","));
        csv2 << "\n";

    }


    // Find non-branching segments.
    vector<bool> isBranchingSegment(segmentCount, false);
    for(SegmentId segmentId=0; segmentId!=segmentCount; segmentId++) {
        const AssemblyGraph::Edge& edge = assemblyGraph.edges[segmentId];
        const AssemblyGraph::VertexId v0 = edge.source;
        const AssemblyGraph::VertexId v1 = edge.target;
        isBranchingSegment[segmentId] =
            (assemblyGraph.edgesBySource[v1].size() >= 2) or
            (assemblyGraph.edgesByTarget[v0].size() >= 2);
    }
    cout << "Found " << std::count(isBranchingSegment.begin(), isBranchingSegment.end(), false) <<
        " non-branching segments." << endl;


    // Use the marker graph pattern to align these sequences.
    using Graph = MarkerGraph2<SegmentId, SegmentId>;
    Graph graph;
    for(SegmentId segmentId=0; segmentId!=segmentCount; segmentId++) {
        if(isBranchingSegment[segmentId]) {
            continue;
        }
        vector<SegmentId> sequenceWithoutBranchingSegments;
        for(const SegmentId segmentId: paths[segmentId]) {
            if(not isBranchingSegment[segmentId]) {
                sequenceWithoutBranchingSegments.push_back(segmentId);
            }
        }
        graph.addSequence(segmentId, sequenceWithoutBranchingSegments);
    }
    const uint64_t disjointSetsSize = graph.doneAddingSequences();
    cout << "The disjoint sets data structure has size " << disjointSetsSize << endl;

    // Align and merge.
    // Just test for now.
    std::set< pair<SegmentId, SegmentId> > pairs;
    for(SegmentId segmentId0=0; segmentId0!=segmentCount; segmentId0++) {
        if(isBranchingSegment[segmentId0]) {
            continue;
        }
        for(const SegmentId segmentId1: paths[segmentId0]) {
            if(isBranchingSegment[segmentId1]) {
                continue;
            }
            if(segmentId1 < segmentId0) {
                pairs.insert(make_pair(segmentId1, segmentId0));
            } else if(segmentId0 < segmentId1) {
                pairs.insert(make_pair(segmentId0, segmentId1));
            }
        }
    }
    cout << "Found " << pairs.size() << " pairs to align." << endl;
    for(const auto& p: pairs) {
        graph.alignAndMerge(p.first, p.second);
    }

    // Finish creation of the marker graph.
    graph.doneMerging();
    cout << "The marker graph for assembly segments has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;
    graph.writeGraphviz("SegmentMarkerGraph.dot");

}
#endif

