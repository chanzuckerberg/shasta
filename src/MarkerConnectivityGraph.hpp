#ifndef SHASTA_MARKER_CONNECTIVITY_GRAPH_HPP
#define SHASTA_MARKER_CONNECTIVITY_GRAPH_HPP

/*******************************************************************************

In the MarkerConnectivityGraph, each vertex represents a marker
of an oriented read, identified by a pair(OrientedRead, marker ordinal)
(a MarkerPair).

Two vertices are joined by an undirected edge if there is an
alignment in which the two markers corresponding to the two vertices
are aligned. We can create a MarkerConnectivityGraph using all alignments,
or just alignments that are in the ReadGraph.

*******************************************************************************/

#include "ReadId.hpp"
#include <boost/graph/adjacency_list.hpp>
#include "utility.hpp"

namespace shasta {
    class MarkerConnectivityGraph;

    using MarkerPair = pair<OrientedReadId, uint32_t>;
    using MarkerConnectivityGraphBaseClass = boost::adjacency_list<
        boost::setS,
        boost::vecS,
        boost::undirectedS,
        MarkerPair
        >;
}



class shasta::MarkerConnectivityGraph : public MarkerConnectivityGraphBaseClass {
public:
};


#endif
