#ifndef SHASTA_PHASING_GRAPH_HPP
#define SHASTA_PHASING_GRAPH_HPP

/*******************************************************************************

We keep track of which oriented reads are internal to 
which assembly graph edges. We then define a phasing similarity
between oriented reads: two oriented reads have high phAsing similarity
if the sets of assembly graph edges they are internal to are similar.
Two reads with high similarity are likely to be phased together.

Using this similarity, we define the Phasing graph, an undirected
graph in which each edge represents an oriented read.
An edge is added between vertices corresponding to oriented reads
with high phasing similarity. 
 
*******************************************************************************/

#include "AssemblyGraph.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
#include "ReadId.hpp"

namespace shasta {
    class PhasingGraph;
}



class shasta::PhasingGraph {
public:    

    // The oriented reads internal to each assembly graph edge.
    // Indexed by assembly graph EdgeId.
    MemoryMapped::VectorOfVectors<OrientedReadId, uint64_t> orientedReads;

    // The assembly graph edges that each oriented read is internal to.
    // Indexed by OrientedReadId::getValue().
    MemoryMapped::VectorOfVectors<AssemblyGraph::EdgeId, uint64_t> assemblyGraphEdges;
};

#endif
