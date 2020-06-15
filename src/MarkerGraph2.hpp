#ifndef SHASTA_MARKER_GRAPH2_HPP
#define SHASTA_MARKER_GRAPH2_HPP

// A general purpose class that implements the marker graph
// pattern for sequences of arbitrary symbols.
// Templatized by the type of the symbol in the
// sequences and the type used to identify the sequences
// to be aligned.

// Shasta
#include "SHASTA_ASSERT.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Seqan.
#include <seqan/align.h>

// Standard library.
#include "fstream.hpp"
#include <map>
#include "string"
#include "utility.hpp"
#include "vector.hpp"



namespace shasta {

    template<class Symbol, class SequenceId> class MarkerGraph2;
    template<class Symbol, class SequenceId> class MarkerGraph2Vertex;
    template<class Symbol, class SequenceId> class MarkerGraph2Edge;

    template<class Symbol, class SequenceId> using MarkerGraph2BaseClass =
        boost::adjacency_list<
        boost::listS,
        boost::vecS,
        boost::bidirectionalS,
        MarkerGraph2Vertex<Symbol, SequenceId>,
        MarkerGraph2Edge<Symbol, SequenceId>
        >;

}



template<class Symbol, class SequenceId> class shasta::MarkerGraph2Vertex {
public:

    uint64_t vertexId;

    // The markers that were merged into this vertex.
    // For each, we store the id of the sequence in which
    // the marker appeared and the position in that sequence.
    vector< pair<SequenceId, uint64_t> > markers;

    MarkerGraph2Vertex(uint64_t vertexId) :
        vertexId(vertexId)
    {}

    // Unclear why the default constructor is needed.
    MarkerGraph2Vertex()
    {
        SHASTA_ASSERT(0);
    }

};



template<class Symbol, class SequenceId> class shasta::MarkerGraph2Edge {
public:

    // Sequences and positions for the markers of this edge.
    // For each one we store the SequenceId and a pair
    // containing positions in the source and target vertex.
    vector< pair<SequenceId, pair<uint64_t, uint64_t> > > markers;
};



template<class Symbol, class SequenceId>
    class shasta::MarkerGraph2 : public MarkerGraph2BaseClass<Symbol, SequenceId> {
public:

    // A Sequence is a vector of Symbol's.
    using Sequence = vector<Symbol>;

    // Add one of the input sequences.
    void addSequence(SequenceId sequenceId, const Sequence&);

    // Call this to indicate that we are done adding sequences.
    // Returns the size of the disjoint sets data structure.
    uint64_t doneAddingSequences();

    // Align two sequences, then merge pairs of aligned markers.
    void alignAndMerge(SequenceId, SequenceId);

    // Call this to indicate that we are done calling alignAndMerge.
    // This creates vertices and edges.
    void doneMerging();

    void writeGraphviz(const string& fileName) const;

    using vertex_descriptor = typename MarkerGraph2BaseClass<Symbol, SequenceId>::vertex_descriptor;
    using edge_descriptor = typename MarkerGraph2BaseClass<Symbol, SequenceId>::edge_descriptor;

    // Get the Symbol corresponding ot all markers in a given vertex.
    Symbol getVertexSymbol(vertex_descriptor) const;

private:

    // Our input sequences and their ids.
    vector<Sequence> sequences;
    vector<SequenceId> sequenceIds;

    // A map that gives the index in sequences and sequenceIds
    // corresponding to each SequenceId.
    std::map<SequenceId, uint64_t> idMap;



    // Keep track of what phase we are in.
    enum class Phase {

        // In this phase we are just storing sequences.
        addingSequences,

        // In this phase we are merging markers
        // using the disjoint set structure.
        merging,

        // In this phase the marker graph has been created.
        // Its vertices and edges are available.
        done
    };
    Phase phase = Phase::addingSequences;


    // The disjoint sets data structure.
    vector<uint64_t> rank;
    vector<uint64_t> parent;
    shared_ptr< boost::disjoint_sets<uint64_t*, uint64_t*> > disjointSetsPointer;

    // The position of the first marker of each sequence
    // in the disjoint data set structure.
    vector<uint64_t> start;

    // The vertex corresponding to each marker of each sequence.
    vector< vector<vertex_descriptor> > vertexTable;

    // Align two sequences and return a vector of identical aligned markers
    // that should be merged.
    // The compute alignment vector contains pairs of marker indexes in
    // the two sequences (marker ordinals) that are aligned and should we merged.
    void align(
        const Sequence&,
        const Sequence&,
        vector< pair<uint64_t, uint64_t> >& alignment
    ) const;


    // Scores used to compute SeqAn alignments.
    int matchScore = 1;
    int mismatchScore = -1;
    int gapScore = -1;

    void createVertices();
    void createEdges();


};



// Add one of the input sequences.
template<class Symbol, class SequenceId>
    void shasta::MarkerGraph2<Symbol, SequenceId>::addSequence(
        SequenceId sequenceId,
        const Sequence& sequence)
{
    SHASTA_ASSERT(phase == Phase::addingSequences);
    SHASTA_ASSERT(idMap.find(sequenceId) == idMap.end());
    idMap.insert(make_pair(sequenceId, sequenceIds.size()));
    sequenceIds.push_back(sequenceId);
    sequences.push_back(sequence);
}



// Call this to indicate that we are done adding sequences.
// Returns the size of the disjoint sets data structure.
template<class Symbol, class SequenceId>
    uint64_t shasta::MarkerGraph2<Symbol, SequenceId>::doneAddingSequences()
{
    SHASTA_ASSERT(sequenceIds.size() == sequences.size());
    SHASTA_ASSERT(phase == Phase::addingSequences);
    phase = Phase::merging;

    // Compute the position of the first marker of each sequence
    // in the disjoint data set structure.
    uint64_t n = 0;
    start.resize(sequences.size());
    for(uint64_t i=0; i<sequences.size(); i++) {
        start[i] = n;
        n += sequences[i].size();
    }

    // Initialize the disjoint sets data structure.
    rank.resize(n);
    parent.resize(n);
    disjointSetsPointer =
        make_shared<boost::disjoint_sets<uint64_t*, uint64_t*> >(&rank[0], &parent[0]);
    auto& disjointSets = *disjointSetsPointer;
    for(uint64_t i=0; i<n; i++) {
        disjointSets.make_set(i);
    }

    return n;
}



// Align two sequences, then merge pairs of aligned markers.
template<class Symbol, class SequenceId>
    void shasta::MarkerGraph2<Symbol, SequenceId>::alignAndMerge(
        SequenceId sequenceId0,
        SequenceId sequenceId1)
{
    SHASTA_ASSERT(phase == Phase::merging);

    // Locate the two sequences.
    const auto it0 = idMap.find(sequenceId0);
    const auto it1 = idMap.find(sequenceId1);
    SHASTA_ASSERT(it0 != idMap.end());
    SHASTA_ASSERT(it1 != idMap.end());
    const Sequence& sequence0 = sequences[it0->second];
    const Sequence& sequence1 = sequences[it1->second];

    // Align them.
    vector< pair<uint64_t, uint64_t> > alignment;
    align(sequence0, sequence1, alignment);
    /*
    cout << "Alignment length " << alignment.size() << endl;
    for(const auto& p: alignment) {
        cout << sequence0[p.first] << " " << sequence1[p.second] << endl;
    }
    */

    // Merge pairs of aligned markers.
    auto& disjointSets = *disjointSetsPointer;
    const uint64_t start0 = start[idMap[sequenceId0]];
    const uint64_t start1 = start[idMap[sequenceId1]];
    for(const auto& p: alignment) {
        const uint64_t i0 = p.first;
        const uint64_t i1 = p.second;
        disjointSets.union_set(start0 + i0, start1 + i1);
    }

}



// Align two sequences and return a vector of identical aligned markers
// that should be merged.
// The compute alignment vector contains pairs of marker indexes in
// the two sequences (marker ordinals) that are aligned and should we merged.
template<class Symbol, class SequenceId>
    void shasta::MarkerGraph2<Symbol, SequenceId>::align(
    const Sequence& sequence0,
    const Sequence& sequence1,
    vector< pair<uint64_t, uint64_t> >& alignment
) const
{
    // Use SeqAn to compute an alignment free at both ends.
    // https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html
    using namespace seqan;

    // Hide shasta::Alignment.
    using seqan::Alignment;

    // An oriented read is represented by its pseudo-path.
    // We want to align a pair of such sequences.
    using TSequence = String<Symbol>;

    // Other SeqAn types we need.
    using TStringSet = StringSet<TSequence>;
    using TDepStringSet = StringSet<TSequence, Dependent<> >;
    using TAlignGraph = Graph<Alignment<TDepStringSet> >;

    // Construct the sequences we want to pass to SeqAn.
    // Add 100 to all Symbols to avoid collisions with the
    // value 45, used by SeqAn to represent gaps.
    TSequence seq0;
    for(const Symbol symbol: sequence0) {
        appendValue(seq0, symbol + 100);
    }
    TSequence seq1;
    for(const Symbol symbol: sequence1) {
        appendValue(seq1, symbol + 100);
    }

    // Store them in a SeqAn string set.
    TStringSet sequences;
    appendValue(sequences, seq0);
    appendValue(sequences, seq1);

    // Compute the alignment.
    TAlignGraph graph(sequences);
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




    // Find pairs of aligned markers.
    alignment.clear();
    uint64_t i0 = 0;
    uint64_t i1 = 0;
    for(uint64_t j=0; j<alignmentLength; j++) {
        const uint64_t value0 = align[j];
        const uint64_t value1 = align[j+alignmentLength];
        if(value0!=45 and value1!=45 and value0 == value1) {
            SHASTA_ASSERT(value0 == sequence0[i0] + 100);
            SHASTA_ASSERT(value1 == sequence1[i1] + 100);
            alignment.push_back(make_pair(i0, i1));
        }
        if(value0 != 45) {
            ++i0;
        }
        if(value1 != 45) {
            ++i1;
        }
    }
    SHASTA_ASSERT(i0 == sequence0.size());
    SHASTA_ASSERT(i1 == sequence1.size());

}


// Call this to indicate that we are done calling alignAndMerge.
// This creates vertices and edges.
template<class Symbol, class SequenceId>
    void shasta::MarkerGraph2<Symbol, SequenceId>::doneMerging()
{
    SHASTA_ASSERT(phase == Phase::merging);
    phase = Phase::done;


    // Create vertices.
    createVertices();

    // Clean up the disjoint sets data structure.
    disjointSetsPointer = 0;
    rank.clear();
    parent.clear();
    rank.shrink_to_fit();
    parent.shrink_to_fit();

    // Create edges.
    createEdges();

}



template<class Symbol, class SequenceId>
    void shasta::MarkerGraph2<Symbol, SequenceId>::createVertices()
{
    MarkerGraph2<Symbol, SequenceId>& graph = *this;

    // Map that gives the vertex corresponding to each
    // disjoint set.
    std::map<uint64_t, vertex_descriptor> vertexMap;

    auto& disjointSets = *disjointSetsPointer;
    uint64_t vertexId = 0;

    // Loop over all sequences.
    vertexTable.resize(sequences.size());
    for(uint64_t i=0; i<sequences.size(); i++) {
        const SequenceId sequenceId = sequenceIds[i];
        const uint64_t sequenceLength = sequences[i].size();
        const uint64_t k = start[i];

        // Loop over markers of this sequence.
        for(uint64_t j=0; j<sequenceLength; j++) {

            // Get the disjoint set.
            const uint64_t disjointSet = disjointSets.find_set(k + j);

            // Get the vertex, creating it if necessary.
            vertex_descriptor v;
            const auto it = vertexMap.find(disjointSet);
            if(it == vertexMap.end()) {
                v = add_vertex(
                    MarkerGraph2Vertex<Symbol, SequenceId>(vertexId++), graph);
                vertexMap.insert(make_pair(disjointSet, v));
            } else {
                v = it->second;
            }

            // Add this marker to the vertex.
            graph[v].markers.push_back(make_pair(sequenceId, j));
            vertexTable[i].push_back(v);
        }
    }
}



template<class Symbol, class SequenceId>
    void shasta::MarkerGraph2<Symbol, SequenceId>::createEdges()
{
    MarkerGraph2<Symbol, SequenceId>& graph = *this;

    for(uint64_t i=0; i<sequences.size(); i++) {
        const SequenceId sequenceId = sequenceIds[i];

        for(uint64_t j=1; j<sequences[i].size(); j++) {
            const vertex_descriptor v0 = vertexTable[i][j-1];
            const vertex_descriptor v1 = vertexTable[i][j];

            edge_descriptor e;
            bool edgeExists = false;
            tie(e, edgeExists) = boost::edge(v0, v1, graph);
            if(not edgeExists) {
                tie(e, ignore) = boost::add_edge(v0, v1, graph);
            }
            graph[e].markers.push_back(make_pair(sequenceId, make_pair(j-1, j)));
        }
    }

}



template<class Symbol, class SequenceId>
    void shasta::MarkerGraph2<Symbol, SequenceId>::writeGraphviz(
    const string& fileName) const
{
    using Graph = MarkerGraph2<Symbol, SequenceId>;
    const Graph& graph = *this;

    ofstream out(fileName);
    out << "digraph MarkerGraph2 {\n";
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        out << graph[v].vertexId << "[label=\"" << getVertexSymbol(v) << "\"];\n";
    }
    BGL_FORALL_EDGES_T(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        out << graph[v0].vertexId << "->" <<
            graph[v1].vertexId << ";\n";
    }
    out << "}\n";

}



// Get the Symbol corresponding ot all markers in a given vertex.
template<class Symbol, class SequenceId>
    Symbol shasta::MarkerGraph2<Symbol, SequenceId>::getVertexSymbol(
    vertex_descriptor v) const
{
    const MarkerGraph2Vertex<Symbol, SequenceId>& vertex = (*this)[v];
    SHASTA_ASSERT(not vertex.markers.empty());
    const auto& p = vertex.markers.front();
    const SequenceId sequenceId = p.first;
    const uint64_t position = p.second;
    const auto it = idMap.find(sequenceId);
    SHASTA_ASSERT(it != idMap.end());
    const uint64_t sequenceIndex = it->second;
    return sequences[sequenceIndex][position];
}

#endif
