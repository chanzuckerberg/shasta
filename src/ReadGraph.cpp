// Shasta.
#include "ReadGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/graphviz.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>

const uint32_t ReadGraph::infiniteDistance = std::numeric_limits<uint32_t>::max();


void RawReadGraph::Visitor::examine_edge(edge_descriptor e, const RawReadGraph& constGraph)
{
    RawReadGraph& graph = const_cast<RawReadGraph&>(constGraph);

    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);
    const bool isSameStrand = graph[e].isSameStrand;

    const uint8_t strand0 = graph[v0].strand;
    SHASTA_ASSERT(strand0 != undiscovered);

    const uint8_t strand1 = isSameStrand ? strand0 : ((~strand0) & 1);
    if(graph[v1].strand == undiscovered) {
        graph[v1].strand = strand1;
    } else {
        if(graph[v1].strand != strand1) {
            crossStrandEdges.push_back(e);
            // cout << "Edge " << v0 << " " << v1 << " " << int(isSameStrand) <<
            //     " marked as cross-strand." << endl;
        }
    }

}



// Write the graph in Graphviz format.
void RawReadGraph::write(const string& fileName) const
{
    ofstream outputFileStream(fileName);
    if(!outputFileStream) {
        throw runtime_error("Error opening " + fileName);
    }
    write(outputFileStream);
}
void RawReadGraph::write(ostream& s) const
{
    Writer writer(*this);
    boost::write_graphviz(s, *this, writer, writer, writer);
}

RawReadGraph::Writer::Writer(const RawReadGraph& graph) :
    graph(graph)
{
}



void RawReadGraph::Writer::operator()(std::ostream& s) const
{
    s << "layout=sfdp;\n";
    s << "ratio=expand;\n";
}


void RawReadGraph::Writer::operator()(std::ostream& s, vertex_descriptor v) const
{
}



void RawReadGraph::Writer::operator()(std::ostream& s, edge_descriptor e) const
{
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);
    if(graph[e].isSameStrand != (graph[v0].strand==graph[v1].strand)) {
        s << "[color=red]";
    }
}



// Compute a shortest path, disregarding edges flagged as cross-strand edges.
void ReadGraph::computeShortPath(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    size_t maxDistance,

    // Edge ids of the shortest path starting at orientedReadId0 and
    // ending at orientedReadId1.
    vector<uint32_t>& path,

    // Work areas.
    vector<uint32_t>& distance, // One per vertex, equals infiniteDistance before and after.
    vector<OrientedReadId>& reachedVertices,   // For which distance is not infiniteDistance.
    vector<uint32_t>& parentEdges  // One per vertex
    )
{
    const bool debug = false;
    if(debug) {
        cout << "Looking for a shortest path between " << orientedReadId0 <<
            " " << orientedReadId1 << endl;
    }

    path.clear();

    // Initialize the BFS.
    std::queue<OrientedReadId> queuedVertices;
    queuedVertices.push(orientedReadId0);
    distance[orientedReadId0.getValue()] = 0;
    reachedVertices.clear();
    reachedVertices.push_back(orientedReadId0);


    // Do the BFS.
    while(!queuedVertices.empty()) {

        // Dequeue a vertex.
        const OrientedReadId vertex0 = queuedVertices.front();
        queuedVertices.pop();
        const uint32_t distance0 = distance[vertex0.getValue()];
        const uint32_t distance1 = distance0 + 1;
        if(debug) {
            cout << "Dequeued " << vertex0 << " at distance " << distance0 << endl;
        }

        // Loop over adjacent vertices.
        bool pathFound = false;
        for(const uint32_t edgeId: connectivity[vertex0.getValue()]) {
            const Edge& edge = edges[edgeId];
            if(edge.crossesStrands) {
                continue;
            }
            const OrientedReadId vertex1 = edge.getOther(vertex0);

            // If we did not encounter this vertex before, process it.
            if(distance[vertex1.getValue()] == infiniteDistance) {
                distance[vertex1.getValue()] = distance1;
                reachedVertices.push_back(vertex1);
                parentEdges[vertex1.getValue()] = edgeId;
                if(distance1 < maxDistance) {
                    if(debug) {
                        cout << "Enqueued " << vertex1 << endl;
                    }
                    queuedVertices.push(vertex1);
                }
            }

            // If we reached orientedReadId1, construct the path.
            if(vertex1 == orientedReadId1) {
                pathFound = true;
                // We have already cleared the path above.
                OrientedReadId vertex = vertex1;
                while(vertex != orientedReadId0) {
                    const uint32_t edgeId = parentEdges[vertex.getValue()];
                    path.push_back(edgeId);
                    vertex = edges[edgeId].getOther(vertex);
                }
                std::reverse(path.begin(), path.end());
                break;
            }

        }
        if(pathFound) {
            break;
        }

    }



    // Clean up.
    for(const OrientedReadId orientedReadId: reachedVertices) {
        distance[orientedReadId.getValue()] = infiniteDistance;
    }
    reachedVertices.clear();
}
