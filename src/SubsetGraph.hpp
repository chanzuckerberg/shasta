#ifndef SHASTA_SUBSET_GRAPH_HPP
#define SHASTA_SUBSET_GRAPH_HPP

/******************************************************************************

A subset graph is a directed graph in which each vertex represents a set
of items.

A directed edge x->y is created if y is a strict subset of x.
That is, is x is a strict superset of y.

Vertices that are not a strict subset of any other vertices are the
"maximal sets" (or "maximal vertices"), and have in-degree 0 in the graph.

Vertices that are not a strict superset of any other vertices are the
"minimal sets" (or "minimal vertices"), and have out-degree 0 in the graph.

Sets that are maximal or minimal are also called "extremal sets"
(or "extremal vertices").

There is a literature for computing the extremal sets efficiently.
For example:
Daniel M. Yellin, Charanjit S. Jutla: Finding Extremal Sets in Less than Quadratic Time.
Inf. Process. Lett. 48(1): 29-34 (1993).

The code below implements the simple algorithm described here:
https://stackoverflow.com/questions/14106121/efficient-algorithm-for-finding-all-maximal-subsets
It is not optimal but sufficient for our purposes.
It can be improving by using bitmaps.

******************************************************************************/

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

//Standard library.
#include "vector.hpp"
#include <map>

namespace shasta {

    // The SubsetGraph is a directed graph in which each vertex represents a set of
    // items stored in a vector (which gets sorted at the beginning of edge creation).
    template<class Item> class SubsetGraph;
    template<class Item> using SubsetGraphBaseClass =
        boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, vector<Item> >;

}


template<class Item> class shasta::SubsetGraph :
    public SubsetGraphBaseClass<Item> {
    public:

    using vertex_descriptor = typename SubsetGraphBaseClass<Item>::vertex_descriptor;

    vertex_descriptor addVertex(const vector<Item>& v)
    {
        return add_vertex(v, *this);
    }

    void createEdges()
    {
        sortSets();
        createIndex();



        // Now the main loop of the algorithm linked above.
        vector<vertex_descriptor> work;

        // Look over all vertices in the graph.
        BGL_FORALL_VERTICES_T(v0, *this, SubsetGraphBaseClass<Item> ) {

            // The intersection set required by the algorithm.
            vector<vertex_descriptor> intersection;
            bool firstTime = true;

            // Loop over all items in the set corresponding to this vertex.
            auto& items0 = (*this)[v0];
            for(const Item& item0: items0) {

                const vector<vertex_descriptor>& v1s = index[item0];
                if(firstTime) {
                    intersection = v1s;
                    firstTime = false;
                } else {
                    work.clear();
                    set_intersection(
                        intersection.begin(), intersection.end(),
                        v1s.begin(), v1s.end(),
                        back_inserter(work));
                    intersection = work;
                }
            }

            // Add the edges.
            for(const vertex_descriptor v1: intersection) {
                if(v1 != v0) {
                    add_edge(v1, v0, *this);
                }
            }
        }
    }



    private:
    void sortSets()
    {
        BGL_FORALL_VERTICES_T(v, *this, SubsetGraphBaseClass<Item> ) {
            auto& items = (*this)[v];
            sort(items.begin(), items.end());
        }
    }

    // And index that, for each Item, it tells us which sets it appears in.
    std::map<Item, vector<vertex_descriptor> > index;
    void createIndex()
    {
        BGL_FORALL_VERTICES_T(v, *this, SubsetGraphBaseClass<Item> ) {
            auto& items = (*this)[v];
            for(const Item& item: items) {
                index[item].push_back(v);
            }
        }

        for(auto& p: index) {
            auto& v = p.second;
            sort(v.begin(), v.end());
        }
    }

};


#endif
