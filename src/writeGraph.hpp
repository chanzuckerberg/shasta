#ifndef SHASTA_WRITE_GRAPH_HPP
#define SHASTA_WRITE_GRAPH_HPP

// Code to write a Boost graph directly, without using Graphviz rendering.

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "iostream.hpp"
#include <functional>
#include <limits>
#include <map>
#include "string.hpp"

namespace shasta {

    namespace WriteGraph {

        class VertexAttributes {
        public:
            double radius = 1.;
            string id;
            string color = "black";
            string tooltip;
            string url;
        };

        class EdgeAttributes {
        public:
            double thickness = 1.;
            string id;
            string color = "black";
            string tooltip;
            string url;
        };

        template<class Graph> void writeSvg(
            const Graph&,
            const string& svgId,
            uint64_t width,
            uint64_t height,
            const std::map<typename Graph::vertex_descriptor, VertexAttributes>&,
            const std::map<typename Graph::edge_descriptor, EdgeAttributes>&,
            ostream&);

        // This method expands on the writeSvg method to allow the Graph Vertex/Edge classes to dictate a numeric
        // ordering with a getSvgOrdering() member function. Since SVG uses the written order of objects to determine
        // the render order of shapes, Edges/Vertexes with lower ordinals will end up underneath those with higher
        // ordinals.
        // Graph - expected to be a graph with edges and vertexes that are accessible by their descriptor objects
        //         using the [] operator. e.g. the boost::adjacency_list
        // svgId - labels the DOM object
        // width and height - the size of the figure in pixels
        template<class Graph> void writeOrderedSvg(
                const Graph& graph,
                const string& svgId,
                uint64_t width,
                uint64_t height,
                const std::map<typename Graph::vertex_descriptor, VertexAttributes>& vertexAttributes,
                const std::map<typename Graph::edge_descriptor, EdgeAttributes>& edgeAttributes,
                ostream& svg);
    }
}


template<class Graph> void shasta::WriteGraph::writeSvg(
    const Graph& graph,
    const string& svgId,
    uint64_t width,
    uint64_t height,
    const std::map<typename Graph::vertex_descriptor, VertexAttributes>& vertexAttributes,
    const std::map<typename Graph::edge_descriptor, EdgeAttributes>& edgeAttributes,
    ostream& svg)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;
    // using edge_descriptor = typename Graph::edge_descriptor;

    // Compute the view box.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::min();
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        const auto& position = graph[v].position;
        VertexAttributes attributes;
        auto it = vertexAttributes.find(v);
        if(it != vertexAttributes.end()) {
            attributes = it->second;
        }
        const double radius = attributes.radius;

        // Update the view box to include this vertex.
        xMin = min(xMin, position[0] - radius);
        xMax = max(xMax, position[0] + radius);
        yMin = min(yMin, position[1] - radius);
        yMax = max(yMax, position[1] + radius);
    }



    // Begin the svg.
    svg << "<svg id='" << svgId << "' width='" << width << "' height='" << height <<
        "' viewbox='" << xMin << " " << yMin << " " << xMax-xMin << " " << yMax-yMin <<
        "'>\n";



    // Write the edges first, so they don't cover the vertices.
    svg << "<g id='" << svgId << "-edges'>\n";
    BGL_FORALL_EDGES_T(e, graph, Graph) {

        // Get the attributes for this vertex.
        EdgeAttributes attributes;
        auto it = edgeAttributes.find(e);
        if(it != edgeAttributes.end()) {
            attributes = it->second;
        }

        // Get vertex positions.
        const vertex_descriptor v1 = source(e, graph);
        const vertex_descriptor v2 = target(e, graph);
        const auto& position1 = graph[v1].position;
        const auto& position2 = graph[v2].position;

        svg << "<line x1='" << position1[0] << "' y1='" << position1[1] <<
            "' x2='" << position2[0] << "' y2='" << position2[1];

        if(not attributes.id.empty()) {
            svg << " id='" << attributes.id << "'";
        }

        svg << "' stroke='" << attributes.color <<
            "' stroke-width='" << attributes.thickness <<
            "'>";

        if(not attributes.tooltip.empty()) {
            svg << "<title>" << attributes.tooltip << "</title>";
        }

        svg << "</line>\n";
    }
    svg << "</g>\n";



    // Write the vertices.
    svg << "<g id='" << svgId << "-vertices' stroke='none'>\n";
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        VertexAttributes attributes;
        auto it = vertexAttributes.find(v);
        if(it != vertexAttributes.end()) {
            attributes = it->second;
        }
        const auto& position = graph[v].position;

        if(not attributes.url.empty()) {
            svg << "<a href='" << attributes.url << "'>";
        }

        svg << "<circle cx='" << position[0] << "' cy='" << position[1] <<
            "' r='" << attributes.radius << "'";

        if(not attributes.id.empty()) {
            svg << " id='" << attributes.id << "'";
        }

        svg << " fill='" << attributes.color << "'>";

        if(not attributes.tooltip.empty()) {
            svg << "<title>" << attributes.tooltip << "</title>";
        }

        svg << "</circle>";

        if(not attributes.url.empty()) {
            svg << "</a>";
        }
        svg << "\n";
    }
    svg << "</g>\n";



    // End the svg.
    svg << "</svg>\n";
}


template<class Graph> void shasta::WriteGraph::writeOrderedSvg(
        const Graph& graph,
        const string& svgId,
        uint64_t width,
        uint64_t height,
        const std::map<typename Graph::vertex_descriptor, VertexAttributes>& vertexAttributes,
        const std::map<typename Graph::edge_descriptor, EdgeAttributes>& edgeAttributes,
        ostream& svg)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;
    using edge_descriptor = typename Graph::edge_descriptor;

    // Compute the view box.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::min();
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        const auto& position = graph[v].position;
        VertexAttributes attributes;
        auto it = vertexAttributes.find(v);
        if(it != vertexAttributes.end()) {
            attributes = it->second;
        }
        const double radius = attributes.radius;

        // Update the view box to include this vertex.
        xMin = min(xMin, position[0] - radius);
        xMax = max(xMax, position[0] + radius);
        yMin = min(yMin, position[1] - radius);
        yMax = max(yMax, position[1] + radius);
    }

    // Begin the svg.
    svg << "<svg id='" << svgId << "' width='" << width << "' height='" << height <<
        "' viewbox='" << xMin << " " << yMin << " " << xMax-xMin << " " << yMax-yMin <<
        "'>\n";

    // Accumulate all the edge data
    std::multimap <uint8_t, pair<edge_descriptor, EdgeAttributes> > orderedEdgeAttributes;
    BGL_FORALL_EDGES_T(e, graph, Graph) {
        // Get the attributes for this edge.
        EdgeAttributes attributes;
        auto it = edgeAttributes.find(e);
        if(it != edgeAttributes.end()) {
            attributes = it->second;
        }

        // Find an ordering specified by the Edge class itself
        auto edge = graph[e];
        auto info = std::make_pair(e, attributes);

        orderedEdgeAttributes.emplace(edge.getSvgOrdering(), info);
    }

    // Write the edges first, so they don't cover the vertices.
    svg << "<g id='" << svgId << "-edges'>\n";
    for(const auto& item: orderedEdgeAttributes) {
        const edge_descriptor& e = item.second.first;
        const EdgeAttributes& attributes = item.second.second;

        // Get vertex positions.
        const vertex_descriptor v1 = source(e, graph);
        const vertex_descriptor v2 = target(e, graph);
        const auto& position1 = graph[v1].position;
        const auto& position2 = graph[v2].position;

        auto edge = graph[e];
        string svgClassName = edge.getSvgClassName();

        svg << "<line class= '" << svgClassName << "'x1='" << position1[0] << "' y1='" << position1[1] <<
            "' x2='" << position2[0] << "' y2='" << position2[1];

        if(not attributes.id.empty()) {
            svg << " id='" << attributes.id << "'";
        }

        svg << "' stroke='" << attributes.color <<
            "' stroke-width='" << attributes.thickness <<
            "'>";

        if(not attributes.tooltip.empty()) {
            svg << "<title>" << attributes.tooltip << "</title>";
        }

        svg << "</line>\n";
    }
    svg << "</g>\n";

    // Accumulate all the vertex data
    std::multimap <uint8_t, pair<vertex_descriptor, VertexAttributes> > orderedVertexAttributes;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        VertexAttributes attributes;
        auto it = vertexAttributes.find(v);
        if (it != vertexAttributes.end()) {
            attributes = it->second;
        }

        // Find an ordering specified by the Vertex class itself
        auto vertex = graph[v];
        auto info = std::make_pair(v, attributes);

        orderedVertexAttributes.emplace(vertex.getSvgOrdering(), info);
    }

    // Write the vertices.
    svg << "<g id='" << svgId << "-vertices' stroke='none'>\n";
    for (const auto& item: orderedVertexAttributes){
        const vertex_descriptor& v = item.second.first;
        const VertexAttributes& attributes = item.second.second;

        const auto& position = graph[v].position;

        if(not attributes.url.empty()) {
            svg << "<a href='" << attributes.url << "'>";
        }

        svg << "<circle cx='" << position[0] << "' cy='" << position[1] <<
            "' r='" << attributes.radius << "'";

        if(not attributes.id.empty()) {
            svg << " id='" << attributes.id << "'";
        }

        svg << " fill='" << attributes.color << "'>";

        if(not attributes.tooltip.empty()) {
            svg << "<title>" << attributes.tooltip << "</title>";
        }

        svg << "</circle>";

        if(not attributes.url.empty()) {
            svg << "</a>";
        }
        svg << "\n";
    }
    svg << "</g>\n";


    // End the svg.
    svg << "</svg>\n";
}


#endif
