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
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include <map>
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
    phase = Phase::merging;

    // Compute the position of the first marker of each sequence
    // in the disjoint data set structure.
    uint64_t n = 0;
    start.resize(sequences.size());
    for(uint64_t i=0; i<sequences.size(); i++) {
        start[i] = n;
        n += sequences[i].size();
    }

    // Initialize the disjoint sets data structure..
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
    // Locate the two sequences.
    const auto it0 = idMap.find(sequenceId0);
    const auto it1 = idMap.find(sequenceId1);
    SHASTA_ASSERT(it0 != idMap.end());
    SHASTA_ASSERT(it1 != idMap.end());
    const Sequence& sequence0 = sequences[it0->second];
    const Sequence& sequence1 = sequences[it1->second];

}


#endif
