// Shasta.
#include "CompressedAssemblyGraph.hpp"
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "findLinearChains.hpp"
#include "html.hpp"
#include "platformDependent.hpp"
#include "runCommandWithTimeout.hpp"
#include "subgraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"
#include "vector.hpp"



// Create the CompressedAssemblyGraph from the AssemblyGraph.
CompressedAssemblyGraph::CompressedAssemblyGraph(
    const Assembler& assembler)
{
    CompressedAssemblyGraph& graph = *this;
    const AssemblyGraph& assemblyGraph = *(assembler.assemblyGraphPointer);

    cout << "The assembly graph has " << assemblyGraph.vertices.size() <<
        " vertices and " << assemblyGraph.edges.size() << " edges." << endl;

    // Create a vertex for each vertex of the assembly graph.
    vector<vertex_descriptor> vertexTable;
    createVertices(assemblyGraph.vertices.size(), vertexTable);

    // Create an edge for each set of parallel edges of the assembly graph.
    createEdges(assemblyGraph, vertexTable);
    removeReverseBubbles();

    // Merge linear chains of edges.
    mergeLinearChains();

    cout << "The compressed assembly graph has " <<
        num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." << endl;

    // Assign an id to each edge.
    assignEdgeIds();

    // Fill in the assembly graph edges that go into each
    // edge of the compressed assembly graph.
    fillContributingEdges(assemblyGraph);

    // Fill in minimum and maximum marker counts for each edge.
    fillMarkerCounts(assemblyGraph);

    // Find the oriented reads that appear in marker graph vertices
    // internal to each edge of the compressed assembly graph.
    findOrientedReads(assembler);
    fillOrientedReadTable(assembler);

    // Find edges that have at least one common oriented read
    // which each edge.
    findRelatedEdges();
}



// Create a vertex for each vertex of the assembly graph.
void CompressedAssemblyGraph::createVertices(
    uint64_t vertexCount,
    vector<vertex_descriptor>& vertexTable)
{
    CompressedAssemblyGraph& graph = *this;

    vertexTable.resize(vertexCount, null_vertex());

    for(VertexId vertexId=0; vertexId<vertexCount; vertexId++) {
        const vertex_descriptor v = add_vertex(CompressedAssemblyGraphVertex(vertexId), graph);
        vertexTable[vertexId] = v;
    }
}



// Create an edge for each set of parallel edges of the assembly graph.
void CompressedAssemblyGraph::createEdges(
    const AssemblyGraph& assemblyGraph,
    const vector<vertex_descriptor>& vertexTable)
{
    CompressedAssemblyGraph& graph = *this;

    // Loop over assembly graph edges.
    for(const AssemblyGraph::Edge& edge: assemblyGraph.edges) {
        const vertex_descriptor v0 = vertexTable[edge.source];
        const vertex_descriptor v1 = vertexTable[edge.target];

        // Do we already have an edge between these two vertices?
        bool edgeExists = false;
        tie(ignore, edgeExists) = boost::edge(v0, v1, graph);

        // If we don't already have an edge between these two vertices, add it.
        // We only create one of each set of parallel edges.
        if(not edgeExists) {
            edge_descriptor e;
            bool edgeWasAdded = false;
            tie(e, edgeWasAdded) = add_edge(v0, v1, graph);
            SHASTA_ASSERT(edgeWasAdded);

            // Store the assembly graph vertices.
            CompressedAssemblyGraphEdge& edge = graph[e];
            edge.vertices.push_back(graph[v0].vertexId);
            edge.vertices.push_back(graph[v1].vertexId);
        }
    }

}


// Remove back edges that create reverse bubbles.
// A reverse bubble is defined by two vertices v0 and v1 such that:
// - Edge v0->v1 exists.
// - Edge v1->v0 exists.
// - Out-degree(v0) = 1
// - In-degree(v1) = 1.
// Under these conditions, we remove vertex v1->v0.
void CompressedAssemblyGraph::removeReverseBubbles()
{
    CompressedAssemblyGraph& graph = *this;

    // Vector to contain the edges that we wil remove.
    vector<edge_descriptor> edgesToBeRemoved;

    // Try all edges.
    BGL_FORALL_EDGES(e01, graph, CompressedAssemblyGraph) {

        // Check that v0 has out-degree 1.
        const vertex_descriptor v0 = source(e01, graph);
        if(out_degree(v0, graph) != 1) {
            continue;
        }

        // Check that v1 has in-degree 1.
        const vertex_descriptor v1 = target(e01, graph);
        if(in_degree(v1, graph) != 1) {
            continue;
        }

        // Look for edges v1->v0.
        // Try all edges v1->v2. Fkag to be renmoved if v2=v0.
        BGL_FORALL_OUTEDGES(v1, e12, graph, CompressedAssemblyGraph) {
            const vertex_descriptor v2 = target(e12, graph);
            if(v2 == v0) {
                edgesToBeRemoved.push_back(e12);
            }
        }
    }

    // Remove the edges we flagged.
    deduplicate(edgesToBeRemoved);
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }

}



// Merge linear chains of edges.
void CompressedAssemblyGraph::mergeLinearChains()
{
    CompressedAssemblyGraph& graph = *this;

    // Find linear chains.
    vector< std::list<edge_descriptor> > chains;
    findLinearChains(graph, chains);



    // Replace each chain with a single edge.
    for(const std::list<edge_descriptor>& chain: chains) {

        // If the chain has length 1, leave it alone.
        if(chain.size() == 1) {
            continue;
        }

        // Add the new edge.
        const vertex_descriptor v0 = source(chain.front(), graph);
        const vertex_descriptor v1 = target(chain.back(), graph);
        edge_descriptor eNew;
        bool edgeWasAdded = false;
        tie(eNew, edgeWasAdded) = add_edge(v0, v1, graph);
        SHASTA_ASSERT(edgeWasAdded);
        CompressedAssemblyGraphEdge& newEdge = graph[eNew];

        // Fill in the assembly graph vertices corresponding to this new edge.
        newEdge.vertices.push_back(graph[v0].vertexId);
        for(const edge_descriptor e: chain) {
            const vertex_descriptor v = target(e, graph);
            newEdge.vertices.push_back(graph[v].vertexId);
        }

        // Remove the edges of the chain.
        // We will remove the vertices later.
        for(const edge_descriptor e: chain) {
            boost::remove_edge(e, graph);
        }
    }



    // Remove the vertices that have become isolated.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, CompressedAssemblyGraph) {
        if(in_degree(v, graph) == 0 and out_degree(v, graph) == 0) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        remove_vertex(v, graph);
    }

}



// Assign an id to each edge.
void CompressedAssemblyGraph::assignEdgeIds()
{
    CompressedAssemblyGraph& graph = *this;

    uint64_t edgeId = 0;
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        graph[e].id = edgeId++;
    }

}



// Fill in the assembly graph edges that go into each
// edge of the compressed assembly graph.
void CompressedAssemblyGraph::fillContributingEdges(
    const AssemblyGraph& assemblyGraph)
{
    CompressedAssemblyGraph& graph = *this;

    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        CompressedAssemblyGraphEdge& edge = graph[e];
        edge.edges.resize(edge.vertices.size() - 1);
        for(uint64_t i=0; i<edge.edges.size(); i++) {
            const VertexId vertexId0 = edge.vertices[i];
            const VertexId vertexId1 = edge.vertices[i+1];
            const span<const EdgeId> edges0 = assemblyGraph.edgesBySource[vertexId0];
            for(const EdgeId edge01: edges0) {
                if(assemblyGraph.edges[edge01].target == vertexId1) {
                    edge.edges[i].push_back(edge01);
                }
            }
        }

    }
}



// Find the oriented reads that appear in marker graph vertices
// internal to each edge of the compressed assembly graph.
void CompressedAssemblyGraph::findOrientedReads(
    const Assembler& assembler)
{
    CompressedAssemblyGraph& graph = *this;

    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        graph[e].findOrientedReads(assembler);
    }

    // Fill in the oriented read table, which tells us
    // which edges each read appears in.
    orientedReadTable.resize(2 * assembler.readCount());
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        for(const OrientedReadId orientedReadId: graph[e].orientedReadIds) {
            orientedReadTable[orientedReadId.getValue()].push_back(e);
        }
    }
}



void CompressedAssemblyGraph::fillOrientedReadTable(
    const Assembler& assembler)
{
    CompressedAssemblyGraph& graph = *this;

    orientedReadTable.clear();
    orientedReadTable.resize(2 * assembler.readCount());
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        for(const OrientedReadId orientedReadId: graph[e].orientedReadIds) {
            orientedReadTable[orientedReadId.getValue()].push_back(e);
        }
    }
}




// Find the oriented reads that appear in marker graph vertices
// internal to an edge of the compressed assembly graph.
void CompressedAssemblyGraphEdge::findOrientedReads(
    const Assembler& assembler)
{
    const AssemblyGraph& assemblyGraph = *assembler.assemblyGraphPointer;

    // Loop over assembly graph edges.
    for(const vector<AssemblyGraph::EdgeId>& edgesHere: edges) {
        for(const AssemblyGraph::EdgeId assemblyGraphEdgeId: edgesHere) {

            // Loop over marker graph edges corresponding to this
            // assembly graph edge.
            const span<const MarkerGraph::EdgeId> markerGraphEdgeIds =
                assemblyGraph.edgeLists[assemblyGraphEdgeId];
            for(const MarkerGraph::EdgeId markerGraphEdgeId: markerGraphEdgeIds) {
                findOrientedReads(assembler, markerGraphEdgeId);
            }
        }
    }


    // Deduplicate oriented reads and count their occurrences.
    deduplicateAndCount(orientedReadIds, orientedReadIdsFrequency);
}



// Append to orientedReadIds the oriented reads that
// appear in a given marker graph edge.
void CompressedAssemblyGraphEdge::findOrientedReads(
    const Assembler& assembler,
    const MarkerGraph::EdgeId& markerGraphEdgeId)
{
    const span<const MarkerInterval> markerIntervals =
        assembler.markerGraph.edgeMarkerIntervals[markerGraphEdgeId];
    for(const MarkerInterval markerInterval: markerIntervals) {
        orientedReadIds.push_back(markerInterval.orientedReadId);
    }
}



// Find edges that have at least one common oriented read
// which each edge.
void CompressedAssemblyGraph::findRelatedEdges()
{
    CompressedAssemblyGraph& graph = *this;
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        findRelatedEdges(e);
    }
}
void CompressedAssemblyGraph::findRelatedEdges(edge_descriptor e0)
{
    CompressedAssemblyGraph& graph = *this;
    CompressedAssemblyGraphEdge& edge0 = graph[e0];
    for(const OrientedReadId orientedReadId: edge0.orientedReadIds) {
        const vector<edge_descriptor>& edges = orientedReadTable[orientedReadId.getValue()];
        for(const edge_descriptor e1: edges) {
            if(e1 != e0) {
                edge0.relatedEdges.push_back(e1);
            }
        }
    }
    deduplicate(edge0.relatedEdges);
    edge0.relatedEdges.shrink_to_fit();

    /*
    cout << edge0.gfaId() << ":";
    for(const edge_descriptor e1: edge0.relatedEdges) {
        cout << " " << graph[e1].gfaId();
    }
    cout << endl;
    */
}


string CompressedAssemblyGraphEdge::gfaId() const
{
    if(edges.size()==1 and edges.front().size()==1) {
        // Return the one and only assembly graph edge associated with this
        // compressed assembly graph edge.
        return to_string(edges.front().front());
    } else {
        // Return the id of this compressed assembly graph edge,
        // prefixed with "C".
        return "C" + to_string(id);
    }
}



// Return the edge with a given GFA id.
pair<CompressedAssemblyGraph::edge_descriptor, bool>
    CompressedAssemblyGraph::getEdgeFromGfaId(
    const string& s) const
{
    const CompressedAssemblyGraph& graph = *this;
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        if(graph[e].gfaId() == s) {
            return make_pair(e, true);
        }
    }
    return make_pair(edge_descriptor(), false);
}



uint64_t CompressedAssemblyGraph::maxPloidy() const
{
    const CompressedAssemblyGraph& graph = *this;
    uint64_t returnValue = 0;
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        returnValue = max(returnValue, graph[e].maxPloidy());
    }
    return returnValue;
}



uint64_t CompressedAssemblyGraphEdge::maxPloidy() const
{
    uint64_t returnValue = 0;
    for(const auto& v: edges) {
        returnValue = max(returnValue, uint64_t(v.size()));
    }
    return returnValue;
}



// GFA output (without sequence).
void CompressedAssemblyGraph::writeGfa(const string& fileName, double basesPerMarker) const
{
    ofstream gfa(fileName);
    writeGfa(gfa, basesPerMarker);
}
void CompressedAssemblyGraph::writeGfa(ostream& gfa, double basesPerMarker) const
{
    const CompressedAssemblyGraph& graph = *this;

    // Write the header line.
    gfa << "H\tVN:Z:1.0\n";


    // Write a segment record for each edge.
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        const CompressedAssemblyGraphEdge& edge = graph[e];

        gfa <<
            "S\t" <<
            edge.gfaId() << "\t" <<
            "*\t" <<
            "LN:i:" << uint64_t(basesPerMarker * 0.5 *
                double(edge.minMarkerCount + edge.maxMarkerCount)) <<
            "\n";
    }


    // Write GFA links.
    // For each vertex in the compressed assembly graph there is a link for
    // each combination of in-edges and out-edges.
    // Therefore each vertex generates a number of
    // links equal to the product of its in-degree and out-degree.
    BGL_FORALL_VERTICES(v, graph, CompressedAssemblyGraph) {
        BGL_FORALL_INEDGES(v, eIn, graph, CompressedAssemblyGraph) {
            BGL_FORALL_OUTEDGES(v, eOut, graph, CompressedAssemblyGraph) {
                gfa <<
                    "L\t" <<
                    graph[eIn].gfaId() << "\t" <<
                    "+\t" <<
                    graph[eOut].gfaId() << "\t" <<
                    "+\t" <<
                    "*\n";
            }
        }
    }
}



void CompressedAssemblyGraph::writeCsv() const
{
    writeCsvEdges();
    writeCsvBubbleChains();
    writeCsvOrientedReadsByEdge();
    writeCsvOrientedReads();
}



void CompressedAssemblyGraph::writeCsvEdges() const
{
    const CompressedAssemblyGraph& graph = *this;

    ofstream csv("CompressedGraph-Edges.csv");
    csv << "Id,GFA id,Source,Target,MinMarkerCount,MaxMarkerCount,OrientedReadsCount,RelatedEdgesCount,\n";
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        const CompressedAssemblyGraphEdge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        csv << edge.id << ",";
        csv << edge.gfaId() << ",";
        csv << graph[v0].vertexId << ",";
        csv << graph[v1].vertexId << ",";
        csv << edge.minMarkerCount << ",";
        csv << edge.maxMarkerCount << ",";
        csv << edge.orientedReadIds.size() << ",";
        csv << edge.relatedEdges.size() << ",";
        csv << "\n";
    }

}



void CompressedAssemblyGraph::writeCsvOrientedReadsByEdge() const
{
    const CompressedAssemblyGraph& graph = *this;

    ofstream csv("CompressedGraph-OrientedReadsByEdge.csv");
    csv << "Id,GFA id,OrientedRead,Frequency\n";
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        const CompressedAssemblyGraphEdge& edge = graph[e];
        SHASTA_ASSERT(edge.orientedReadIds.size() == edge.orientedReadIdsFrequency.size());
        for(uint64_t i=0; i<edge.orientedReadIds.size(); i++) {
            const OrientedReadId orientedReadId = edge.orientedReadIds[i];
            const uint64_t frequency = edge.orientedReadIdsFrequency[i];
            csv << edge.id << ",";
            csv << edge.gfaId() << ",";
            csv << orientedReadId << ",";
            csv << frequency << "\n";
        }
    }

}



void CompressedAssemblyGraph::writeCsvBubbleChains() const
{
    const CompressedAssemblyGraph& graph = *this;
    const uint64_t maxPloidy = graph.maxPloidy();

    ofstream csv("CompressedGraph-BubbleChains.csv");
    csv << "Id,GFA id,Position,";
    for(uint64_t i=0; i<maxPloidy; i++) {
        csv << "Edge" << i << ",";
    }
    csv << "\n";

    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        const CompressedAssemblyGraphEdge& edge = graph[e];

        for(uint64_t position=0; position<edge.edges.size(); position++) {
            const vector<AssemblyGraph::EdgeId>& edgesAtPosition = edge.edges[position];

            csv << edge.id << ",";
            csv << edge.gfaId() << ",";
            csv << position << ",";
            for(const AssemblyGraph::EdgeId edgeId: edgesAtPosition) {
                csv << edgeId << ",";
            }
            csv << "\n";
        }
    }
}



void CompressedAssemblyGraph::writeCsvOrientedReads() const
{
    const CompressedAssemblyGraph& graph = *this;

    ofstream csv("CompressedGraph-OrientedReads.csv");
    csv << "OrientedReadId,Id,GFA id,\n";
    for(OrientedReadId::Int orientedReadId=0; orientedReadId<orientedReadTable.size();
        orientedReadId++) {
        const vector<edge_descriptor>& edges = orientedReadTable[orientedReadId];
        for(const edge_descriptor e: edges) {
            const CompressedAssemblyGraphEdge& edge = graph[e];
            csv << OrientedReadId(orientedReadId) << ",";
            csv << edge.id << ",";
            csv << edge.gfaId() << "\n";
        }
    }
}



// Fill in minimum and maximum marker counts for each edge.
void CompressedAssemblyGraph::fillMarkerCounts(const AssemblyGraph& assemblyGraph)
{
    CompressedAssemblyGraph& graph = *this;

    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        graph[e].fillMarkerCounts(assemblyGraph);
    }
}
void CompressedAssemblyGraphEdge::fillMarkerCounts(const AssemblyGraph& assemblyGraph)
{
    minMarkerCount = 0;
    maxMarkerCount = 0;
    for(const vector<AssemblyGraph::EdgeId>& parallelEdges: edges) {
        SHASTA_ASSERT(not parallelEdges.empty());

        // Compute the minimum and maximum number of markers
        // over this set of parallel edges.
        uint64_t minMarkerCountHere = std::numeric_limits<uint64_t>::max();
        uint64_t maxMarkerCountHere = 0;
        for(const AssemblyGraph::EdgeId edgeId: parallelEdges) {
            const uint64_t markerCount = assemblyGraph.edgeLists.size(edgeId);
            minMarkerCountHere = min(minMarkerCountHere, markerCount);
            maxMarkerCountHere = max(maxMarkerCountHere, markerCount);
        }

        // Update the totals.
        minMarkerCount += minMarkerCountHere;
        maxMarkerCount += maxMarkerCountHere;
    }
}


// Create a local subgraph.
// See createLocalSubgraph for argument explanation.
CompressedAssemblyGraph::CompressedAssemblyGraph(
    const CompressedAssemblyGraph& graph,
    const Assembler& assembler,
    const vector<vertex_descriptor>& startVertices,
    uint64_t maxDistance,
    boost::bimap<vertex_descriptor, vertex_descriptor>& vertexMap,
    boost::bimap<edge_descriptor, edge_descriptor>& edgeMap,
    std::map<vertex_descriptor, uint64_t>& distanceMap
    )
{
    CompressedAssemblyGraph& subgraph = *this;
    createLocalSubgraph(
        graph,
        startVertices,
        maxDistance,
        subgraph,
        vertexMap,
        edgeMap,
        distanceMap);

    // Make sure the relatedEdges of each edge contain
    // edge descriptors in the subgraph.
    BGL_FORALL_EDGES(e0, subgraph, CompressedAssemblyGraph) {
        vector<edge_descriptor> relatedEdges;
        for(const edge_descriptor e1: subgraph[e0].relatedEdges) {
            const auto it = edgeMap.right.find(e1);
            if(it != edgeMap.right.end()) {
                relatedEdges.push_back(it->second);
            }
        }
        subgraph[e0].relatedEdges.swap(relatedEdges);
    }

    subgraph.fillOrientedReadTable(assembler);
}



// Graphviz output.
void CompressedAssemblyGraph::writeGraphviz(
    const string& fileName,
    uint64_t sizePixels,
    double vertexScalingFactor,
    double edgeLengthScalingFactor,
    double edgeThicknessScalingFactor,
    double edgeArrowScalingFactor,
    std::map<vertex_descriptor, array<double, 2 > >& vertexPositions) const
{
    ofstream s(fileName);
    writeGraphviz(s,
        sizePixels,
        vertexScalingFactor,
        edgeLengthScalingFactor,
        edgeThicknessScalingFactor,
        edgeArrowScalingFactor,
        vertexPositions) ;
}

void CompressedAssemblyGraph::writeGraphviz(
    ostream& s,
    uint64_t sizePixels,
    double vertexScalingFactor,
    double edgeLengthScalingFactor,
    double edgeThicknessScalingFactor,
    double edgeArrowScalingFactor,
    std::map<vertex_descriptor, array<double, 2 > >& vertexPositions) const
{
    const CompressedAssemblyGraph& graph = *this;

    s << "digraph CompressedAssemblyGraph {\n"
        "layout=neato;\n"
        "size=" << uint64_t(double(sizePixels)/72.) << ";\n"
        "ratio=expand;\n"
        "splines=true;\n"
        "node [shape=point];\n"
        "node [width=" << vertexScalingFactor << "];\n"
        "edge [penwidth=" << edgeThicknessScalingFactor << "];\n"
        "edge [arrowsize=" << edgeArrowScalingFactor << "];\n";



    // Write the vertices.
    BGL_FORALL_VERTICES(v, graph, CompressedAssemblyGraph) {
        const auto it = vertexPositions.find(v);
        SHASTA_ASSERT(it != vertexPositions.end());
        const auto& x = it->second;

        s << graph[v].vertexId <<
            " [pos=\"" << x[0] << "," << x[1] << "\"];\n";
    }



    // Write the edges.
    // Each edge is written as a number of dummy edges,
    // to make it look longer, proportionally to its number of markers.
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        const CompressedAssemblyGraphEdge& edge = graph[e];
        const string gfaId = edge.gfaId();
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        // To color the edge, use a hash function of the edge id,
        // so the same edge always gets colored the same way.
        const uint32_t hashValue = MurmurHash2(&edge.id, sizeof(edge.id), 757);
        const double H = double(hashValue) / double(std::numeric_limits<uint32_t>::max());
        const double S = 0.7;
        const double V = 0.7;

        s <<
            graph[v0].vertexId << "->" <<
            graph[v1].vertexId <<
            "[tooltip=\"" << gfaId << "\" "
            "color = \"" << H << "," << S << "," << "," << V << "\"" <<
            "];\n";
    }

    s << "}";
}



#if 0
void CompressedAssemblyGraph::writeGraphviz(
    ostream& s,
    uint64_t sizePixels,
    double vertexScalingFactor,
    double edgeLengthScalingFactor,
    double edgeThicknessScalingFactor,
    double edgeArrowScalingFactor,
    std::map<vertex_descriptor, array<double, 2 > >& vertexPositions) const
{
    const CompressedAssemblyGraph& graph = *this;

    s << "digraph CompressedAssemblyGraph {\n"
        "layout=sfdp;\n"
        "size=" << uint64_t(double(sizePixels)/72.) << ";\n"
        "ratio=expand;\n"
        "rankdir=LR;\n"
        "node [shape=point];\n"
        "node [width=" << vertexScalingFactor << "];\n"
        "edge [penwidth=" << edgeThicknessScalingFactor << "];\n"
        "edge [arrowsize=" << edgeArrowScalingFactor << "];\n"
        "edge [arrowhead=none];\n";

    // Write the vertices.
    BGL_FORALL_VERTICES(v, graph, CompressedAssemblyGraph) {
        s << graph[v].vertexId << ";\n";
    }



    // Write the edges.
    // Each edge is written as a number of dummy edges,
    // to make it look longer, proportionally to its number of markers.
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        const CompressedAssemblyGraphEdge& edge = graph[e];
        const string gfaId = edge.gfaId();
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        const uint64_t dummyEdgeCount =
            max(uint64_t(1), uint64_t(0.5 + edgeLengthScalingFactor * edge.averageMarkerCount()));

        // Write the dummy vertices.
        for(uint64_t i=1; i<dummyEdgeCount; i++) {
            s << "\"" << gfaId << "-dummy" << i << "\" [width=" <<
                edgeThicknessScalingFactor/72 << "];\n";
        }


        // Write the dummy edges.
        for(uint64_t i=0; i<dummyEdgeCount; i++) {

            // First vertex - either v0 or a dummy vertex.
            s << "\"";
            if(i == 0) {
                s << graph[v0].vertexId;
            } else {
                s << gfaId << "-dummy" << i;
            }

            s << "\"->\"";

            // Second vertex - either v1 or a dummy vertex.
            if(i == dummyEdgeCount-1) {
                s << graph[v1].vertexId;
            } else {
                s << gfaId << "-dummy" << i+1;
            }
            s << "\"";

            s << " [";
            if(i == dummyEdgeCount-1) {
                s << "style=tapered";
            }
            /*
            if(i == dummyEdgeCount/2) {
                s << " label=" << edge.gfaId();
            }
            */
            s << "]";

            s << ";\n";
        }
    }



    s << "}";
}
#endif


// Use sfdp to compute a vertex layout.
// But sfdp does not support edges of different length
// (it tries to make all edges the same length).
// Here, we expand each edge into a variable number
// of dummy edges. Longer edges are expanded into more
// dummy edges. Then we use sfdp to compute a layout
// including all these dummy edges (and their vertices).
// Finally, we extract the coordinates for the original
// vertices only.
void CompressedAssemblyGraph::computeVertexLayout(
    uint64_t sizePixels,
    double vertexScalingFactor,
    double edgeLengthPower,
    double edgeLengthScalingFactor,
    double timeout,
    std::map<vertex_descriptor, array<double, 2 > >& vertexPositions
) const
{
    const CompressedAssemblyGraph& graph = *this;

    // Create a dot file to contain the graph including the dummy edges.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    ofstream graphOut(dotFileName);
    graphOut <<
        "digraph G{\n"
        "layout=sfdp;\n"
        "smoothing=triangle;\n"
        "overlap=false;\n" <<
        "size=" << uint64_t(double(sizePixels)/72.) << ";\n"
        "ratio=expand;\n"
        "node [shape=point];\n"
        "node [width=" << vertexScalingFactor << "];\n";

    // Write the vertices.
    BGL_FORALL_VERTICES(v, graph, CompressedAssemblyGraph) {
        graphOut << graph[v].vertexId << ";\n";
    }



    // Write the dummy edges.
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        const CompressedAssemblyGraphEdge& edge = graph[e];
        const string gfaId = edge.gfaId();
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        // Figure out the number of dummy edges, based on the number of
        // markers on this edge.
        const double desiredDummyEdgeCount =
            edgeLengthScalingFactor * std::pow(edge.averageMarkerCount(), edgeLengthPower);
        const uint64_t dummyEdgeCount = max(
            uint64_t(1),
            uint64_t(0.5 + desiredDummyEdgeCount));

        // Write the dummy edges.
        for(uint64_t i=0; i<dummyEdgeCount; i++) {

            // First vertex - either v0 or a dummy vertex.
            if(i == 0) {
                graphOut << graph[v0].vertexId;
            } else {
                graphOut << "D" << i << "_" << gfaId;
            }

            graphOut << "->";

            // Second vertex - either v1 or a dummy vertex.
            if(i == dummyEdgeCount-1) {
                graphOut << graph[v1].vertexId;
            } else {
                graphOut << "D" << i+1 << "_" << gfaId;
            }

            graphOut << ";\n";
        }
    }
    graphOut << "}";
    graphOut.close();



    // Now use sfdp to compute the layout.
    // Use plain format output described here
    // https://www.graphviz.org/doc/info/output.html#d:plain
    const string plainFileName = dotFileName + ".txt";
    const string command = "sfdp -T plain " + dotFileName + " -o " + plainFileName;
    bool timeoutTriggered = false;
    bool signalOccurred = false;
    int returnCode = 0;
    runCommandWithTimeout(command, timeout,
        timeoutTriggered, signalOccurred, returnCode);
    if(signalOccurred) {
        throw runtime_error("Unable to compute graph layout: terminated by a signal. "
            "The failing command was: " + command);
    }
    if(timeoutTriggered) {
        throw runtime_error("Timeout exceeded during graph layout computation. "
            "Increase the timeout or decrease the maximum distance to simplify the graph");
    }
    if(returnCode!=0 ) {
        throw runtime_error("Unable to compute graph layout: return code " +
            to_string(returnCode) +
            ". The failing command was: " + command);
    }
    filesystem::remove(dotFileName);

    // Map vertex ids to vertex descriptors.
    std::map<uint64_t, vertex_descriptor> vertexMap;
    BGL_FORALL_VERTICES(v, graph, CompressedAssemblyGraph) {
        vertexMap.insert(make_pair(graph[v].vertexId, v));
    }



    // Extract vertex coordinates.
    ifstream plainFile(plainFileName);
    string line;
    vector<string> tokens;
    while(true) {

        // Read the next line.
        std::getline(plainFile, line);
        if( not plainFile) {
            break;
        }

        // Parse it.
        boost::algorithm::split(tokens, line, boost::algorithm::is_any_of(" "));
        // SHASTA_ASSERT(not tokens.empty());

        // If not a line describing a vertex, skip it.
        if(tokens.front() != "node") {
            continue;
        }
        SHASTA_ASSERT(tokens.size() >= 4);

        // If a dummy vertex, skip it.
        const string& vertexName = tokens[1];
        SHASTA_ASSERT(not vertexName.empty());
        if(vertexName[0] == 'D') {
            continue;
        }
        // Get the vertex id.
        const uint64_t vertexId = boost::lexical_cast<uint64_t>(vertexName);

        // Get the corresponding vertex descriptor.
        const auto it = vertexMap.find(vertexId);
        const vertex_descriptor v = it->second;

        // Store it in the layout.
        array<double, 2> x;
        x[0] = boost::lexical_cast<double>(tokens[2]);
        x[1] = boost::lexical_cast<double>(tokens[3]);
        vertexPositions.insert(make_pair(v, x));
    }
    plainFile.close();
    filesystem::remove(plainFileName);
}



// Create a csv file with coloring.
// If the string passed in is an oriented read,
// it colors all edges that have that read.
// If it is the gfaId of an edge, it colors that edge in red
// and all related edges in green.
// This can be loaded in Bandage to color the edges.
void CompressedAssemblyGraph::color(
    const string& s,
    const string& fileName) const
{
    ofstream csv(fileName);
    color(s, csv);
}
void CompressedAssemblyGraph::color(
    const string& s,
    ostream& csv) const
{
    const CompressedAssemblyGraph& graph = *this;
    std::map<edge_descriptor, string> colorMap;



    // If it is the gfaId of an edge, color that edge in red
    // and all related edges in green.
    edge_descriptor e0;
    bool edgeWasFound = false;
    tie(e0, edgeWasFound) = getEdgeFromGfaId(s);
    if(edgeWasFound) {
        colorMap.insert(make_pair(e0, "Red"));
        for(const edge_descriptor e1: graph[e0].relatedEdges) {
            colorMap.insert(make_pair(e1, "Green"));
        }
    }


    // If it is an oriented read id, color in green all segments
    // that have that oriented read.
    else {
        try {
            const OrientedReadId orientedReadId = OrientedReadId(s);
            for(const edge_descriptor e1: orientedReadTable[orientedReadId.getValue()]) {
                colorMap.insert(make_pair(e1, "Green"));
            }
        } catch(...) {
            cout << "The string to color by does not correspond to a valid "
                "GFA id or a valid oriented read id.\n";
        }
    }



    csv << "Segment,Color\n";
    BGL_FORALL_EDGES(e, graph, CompressedAssemblyGraph) {
        csv << graph[e].gfaId() << ",";
        const auto it = colorMap.find(e);
        if(it == colorMap.end()) {
            csv << "Grey";
        } else {
            csv << it->second;
        }
        csv << "\n";
    }

}
