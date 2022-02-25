#ifdef SHASTA_HTTP_SERVER



// Shasta.
#include "mode3-LocalAssemblyGraph.hpp"
#include "computeLayout.hpp"
#include "HttpServer.hpp"
#include "MarkerGraph.hpp"
#include "MurmurHash2.hpp"
#include "writeGraph.hpp"
using namespace shasta;
using namespace mode3;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/topology.hpp>

// Standard library.
#include <map>
#include <queue>
#include "tuple.hpp"



// Create the LocalAssemblyGraph using a BFS
// that starts at the specified vertex and moves away
// (in both directions) up to the specified distance
mode3::LocalAssemblyGraph::LocalAssemblyGraph(
    const MarkerGraph& markerGraph,
    const AssemblyGraph& assemblyGraph,
    uint64_t startSegmentId,
    uint64_t maxDistance) :
    markerGraph(markerGraph),
    assemblyGraph(assemblyGraph),
    maxDistance(maxDistance)
{
    LocalAssemblyGraph& localAssemblyGraph= *this;

    // The BFS queue.
    std::queue<uint64_t> q;

    // Map segments in the AssemblyGraph to vertices in
    // the LocalAssemblyGraph.
    std::map<uint64_t, vertex_descriptor> segmentMap;

    // Initialize the BFS.
    if(maxDistance > 0) {
        q.push(startSegmentId);
    }
    const vertex_descriptor vStart = addVertex(startSegmentId, 0);
    segmentMap.insert(make_pair(startSegmentId, vStart));



    // BFS.
    while(not q.empty()) {

        // Dequeue a segment.
        const uint64_t segmentId0 = q.front();
        q.pop();
        const vertex_descriptor v0 = segmentMap[segmentId0];
        const uint64_t distance0 = localAssemblyGraph[v0].distance;
        const uint64_t distance1 = distance0 + 1;

        // Loop over children.
        for(const uint64_t linkId: assemblyGraph.linksBySource[segmentId0]) {
            const Link& link = assemblyGraph.links[linkId];
            const uint64_t segmentId1 = link.segmentId1;
            if(segmentMap.find(segmentId1) != segmentMap.end()) {
                // We already encountered this segment.
                continue;
            }
            const vertex_descriptor v1 = addVertex(segmentId1, distance1);
            segmentMap.insert(make_pair(segmentId1, v1));
            if(distance1 < maxDistance) {
                q.push(segmentId1);
            }
        }

        // Loop over parents.
        for(const uint64_t linkId: assemblyGraph.linksByTarget[segmentId0]) {
            const Link& link = assemblyGraph.links[linkId];
            const uint64_t segmentId1 = link.segmentId0;
            if(segmentMap.find(segmentId1) != segmentMap.end()) {
                // We already encountered this segment.
                continue;
            }
            const vertex_descriptor v1 = addVertex(segmentId1, distance1);
            segmentMap.insert(make_pair(segmentId1, v1));
            if(distance1 < maxDistance) {
                q.push(segmentId1);
            }
        }
    }



    // Add the edges.
    for(const auto& p: segmentMap) {
        const uint64_t segmentId0 = p.first;
        const vertex_descriptor v0 = p.second;

        for(const uint64_t linkId: assemblyGraph.linksBySource[segmentId0]) {
            const Link& link = assemblyGraph.links[linkId];
            const uint64_t segmentId1 = link.segmentId1;
            const auto it1 = segmentMap.find(segmentId1);
            if(it1 == segmentMap.end()) {
                continue;
            }
            const vertex_descriptor v1 = it1->second;
            boost::add_edge(v0, v1, LocalAssemblyGraphEdge(linkId), localAssemblyGraph);
        }
    }

}



mode3::LocalAssemblyGraphVertex::LocalAssemblyGraphVertex(
    uint64_t segmentId,
    uint64_t distance) :
    segmentId(segmentId),
    distance(distance)
{
}



mode3::LocalAssemblyGraphVertex::LocalAssemblyGraphVertex() :
    segmentId(0),
    distance(0)
{
}



mode3::LocalAssemblyGraph::vertex_descriptor mode3::LocalAssemblyGraph::addVertex(
    uint64_t segmentId,
    uint64_t distance)
{
    return add_vertex(LocalAssemblyGraphVertex(segmentId, distance), *this);
}






void mode3::LocalAssemblyGraph::writeSvg(
    const string& fileName,
    const SvgOptions& options) const
{
    ofstream svg(fileName);
    writeSvg(svg, options);
}
void mode3::LocalAssemblyGraph::writeSvg(
    ostream& svg,
    const SvgOptions& options) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;

    // Compute the view box.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::min();
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex = localAssemblyGraph[v];
        SHASTA_ASSERT(vertex.position.size() >= 2);
        const array<double, 2>& p1 = vertex.position.front();
        const array<double, 2>& p2 = vertex.position.back();

        // Update the view box to include these end points.
        xMin = min(xMin, p1[0]);
        xMax = max(xMax, p1[0]);
        yMin = min(yMin, p1[1]);
        yMax = max(yMax, p1[1]);
        xMin = min(xMin, p2[0]);
        xMax = max(xMax, p2[0]);
        yMin = min(yMin, p2[1]);
        yMax = max(yMax, p2[1]);
    }
    // Add a bit of extra space.
    const double dx = 0.05 * max(xMax-xMin, yMax-yMin);
    const double dy = 0.05 * max(xMax-xMin, yMax-yMin);
    xMin -= dx;
    xMax += dx;
    yMin -= dy;
    yMax += dy;

    // Figure out the required size of the svg.
    const uint64_t svgSizeX = uint64_t(std::round(options.pixelsPerUnitLength * (xMax - xMin)));
    const uint64_t svgSizeY = uint64_t(std::round(options.pixelsPerUnitLength * (yMax - yMin)));

    // Begin the svg.
    const string svgId = "LocalAssemblyGraph";
    svg << "\n<svg id='" << svgId <<
        "' width='" <<  svgSizeX <<
        "' height='" << svgSizeY <<
        "' viewbox='" << xMin << " " << yMin << " " <<
        max(xMax-xMin, yMax-yMin) << " " <<
        max(xMax-xMin, yMax-yMin) <<
        "'>\n";



    // Write the links first, so they don't overwrite the segments.
    svg << "<g id='" << svgId << "-links'>\n";
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const uint64_t linkId = localAssemblyGraph[e].linkId;

        // Access the LocalAssemblyGraph vertices corresponding to
        // the two segments of this Link and extract some information
        // from them.
        const vertex_descriptor v1 = source(e, localAssemblyGraph);
        const vertex_descriptor v2 = target(e, localAssemblyGraph);
        const LocalAssemblyGraphVertex& vertex1 = localAssemblyGraph[v1];
        const LocalAssemblyGraphVertex& vertex2 = localAssemblyGraph[v2];
        const uint64_t segmentId1 = vertex1.segmentId;
        const uint64_t segmentId2 = vertex2.segmentId;
        const bool areConsecutivePaths = haveConsecutivePaths(v1, v2);

        // Get the positions of the ends of this link.
        SHASTA_ASSERT(vertex1.position.size() >= 2);
        SHASTA_ASSERT(vertex2.position.size() >= 2);
        const array<double, 2>& p1 = vertex1.position.back();
        const array<double, 2>& p2 = vertex2.position.front();

        const double linkThickness =
            min(options.segmentThickness,
            options.minimumLinkThickness +  options.additionalLinkThicknessPerRead * double(assemblyGraph.linkCoverage(linkId) - 1));

        const string dash =
            areConsecutivePaths ? "" :
            " stroke-dasharray='0 " + to_string(1.5 * linkThickness) + "'";


        svg <<
            "<g><title>"
            "Link " << linkId <<
            " from segment " << segmentId1 <<
            " to segment " << segmentId2 <<
            ", coverage " << assemblyGraph.linkCoverage(linkId) <<
            "</title>"
            "<path d='M " << p1[0] << " " << p1[1] << " L " << p2[0] << " " << p2[1] << "'"
            // " stroke='" << (areConsecutivePaths ?  options.linkColor : options.nonConsecutiveLinkColor ) << "'" <<
            " stroke='" << options.linkColor << "'" <<
            dash <<
            " stroke-width='" << linkThickness << "'"
            " stroke-linecap='round'"
            " fill='transparent'"
            // " vector-effect='non-scaling-stroke'"
            "/></g>\n";

    }
    svg << "</g>\n";



    // Write the segments.
    svg << "<g id='" << svgId << "-segments'>\n";
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex = localAssemblyGraph[v];
        const uint64_t distance = localAssemblyGraph[v].distance;

        // Get the positions of the ends of this segment.
        SHASTA_ASSERT(vertex.position.size() >= 2);
        const array<double, 2>& p1 = vertex.position.front();
        const array<double, 2>& p2 = vertex.position.back();
        const array<double, 2>& q1 = vertex.q1;
        const array<double, 2>& q2 = vertex.q2;

        const uint64_t segmentId = localAssemblyGraph[v].segmentId;
        const string color =
            (distance == maxDistance) ?
            options.segmentAtMaxDistanceColor :
            segmentColor(segmentId);

        // Create a marker to show the arrow for this segment.
        const string arrowMarkerName = "arrow" + to_string(segmentId);
        svg <<
            "<defs>\n"
            "<marker id='" << arrowMarkerName <<
            "' viewBox='0 0 0.51 1'\n"
            "refX='0.01' refY='0.5'\n"
            "markerUnits='strokeWidth'\n"
            "markerWidth='0.51' markerHeight='1'\n"
            "orient='auto'>\n"
            "<path d='M 0 0 L 0.01 0 L 0.51 0.5 L 0.01 1 L 0 1 z' "
            "fill='" << color << "' "
            "/>\n"
            "</marker>\n"
            "</defs>\n";

        // Add this segment to the svg.
        svg <<
            "<g><title>"
            "Segment " << segmentId <<
            ", path length " << assemblyGraph.paths.size(segmentId) <<
            ", distance " << distance <<
            "</title>"
            "<path id='Segment-" << segmentId <<
            // Straight line.
            // "' d='M " <<
            // p1[0] << " " << p1[1] << " L " <<
            // p2[0] << " " << p2[1] << "'" <<
            // Cubic using the control points.
            "' d='M " <<
            p1[0] << " " << p1[1] << " C " <<
            q1[0] << " " << q1[1] << ", " <<
            q2[0] << " " << q2[1] << ", " <<
            p2[0] << " " << p2[1] << "'" <<
            " stroke='" << color << "'"
            " stroke-width='" <<  options.segmentThickness << "'"
            " fill='transparent'"
            // " vector-effect='non-scaling-stroke'"
            " marker-end='url(#" <<
            arrowMarkerName <<
            ")'"
            "/></g>\n";
    }
    svg << "</g>\n";



    // End the svg.
    svg << "</svg>\n";
}



void mode3::LocalAssemblyGraph::computeLayout(
    const SvgOptions& options)
{
    LocalAssemblyGraph& localAssemblyGraph = *this;


    // Create an auxiliary graph with two vertices for each segment.
    using G = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
    G g;
    std::map<vertex_descriptor, array<G::vertex_descriptor, 2> > vertexMap;
    std::map<G::edge_descriptor, double> edgeLengthMap;
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const uint64_t segmentId = localAssemblyGraph[v].segmentId;

        const uint64_t pathLength = assemblyGraph.paths.size(segmentId);
        const double displayLength =
            options.minimumSegmentLength +
            double(pathLength - 1) * options.additionalSegmentLengthPerMarker;

        // Add the auxiliary vertices.
        array<G::vertex_descriptor, 2>& auxiliaryVertices = vertexMap[v];
        for(uint64_t i=0; i<2; i++) {
            auxiliaryVertices[i] = boost::add_vertex(g);
        }

        // Add the edge between these auxiliary vertices.
        G::edge_descriptor e;
        tie(e, ignore) = boost::add_edge(auxiliaryVertices[0], auxiliaryVertices[1], g);
        edgeLengthMap.insert(make_pair(e, displayLength));
    }



    // Add auxiliary graph edges between vertices corresponding to different
    // LocalAssemblyGraph vertices.
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const vertex_descriptor v1 = source(e, localAssemblyGraph);
        const vertex_descriptor v2 = target(e, localAssemblyGraph);

        double edgeLength;
        if(haveConsecutivePaths(v1, v2)) {
            edgeLength = options.minimumLinkLength;
        } else {
            const double linkSeparation = max(this->linkSeparation(e), 0.);
            edgeLength = 3. * options.minimumLinkLength + linkSeparation * options.additionalLinkLengthPerMarker;
        }
        G::edge_descriptor eAuxiliary;
        tie(eAuxiliary, ignore) = add_edge(
            vertexMap[v1].back(),
            vertexMap[v2].front(),
            g);
        edgeLengthMap.insert(make_pair(eAuxiliary, edgeLength));
    }



    // Compute the layout of the auxiliary graph.
    std::map<G::vertex_descriptor, array<double, 2> > positionMap;
    if(options.layoutMethod == "neato") {
        if(shasta::computeLayoutGraphviz(g, "neato", 30., positionMap, "", &edgeLengthMap) !=
            ComputeLayoutReturnCode::Success) {
            throw runtime_error("Graph layout failed.");
        }
    } else if(options.layoutMethod == "custom") {
        if(shasta::computeLayoutCustom(g, edgeLengthMap, positionMap, 30.) !=
            ComputeLayoutReturnCode::Success) {
            throw runtime_error("Graph layout failed.");
        }
    } else {
        throw runtime_error("Invalid layout method specified: " + options.layoutMethod);
    }



    // Store the layout in the vertices of the localAssemblyGraph.
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        LocalAssemblyGraphVertex& vertex = localAssemblyGraph[v];
        vertex.position.clear();

        // Locate the auxiliary vertices corresponding to this segment.
        auto it = vertexMap.find(v);
        SHASTA_ASSERT(it != vertexMap.end());
        const array<G::vertex_descriptor, 2>& auxiliaryVertices = it->second;

        // Loop over the auxiliary vertices.
        for(const G::vertex_descriptor u: auxiliaryVertices) {
            auto jt = positionMap.find(u);
            SHASTA_ASSERT(jt != positionMap.end());
            vertex.position.push_back(jt->second);
        }
    }
}



void LocalAssemblyGraph::computeControlPoints()
{
    LocalAssemblyGraph& localAssemblyGraph = *this;
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        computeControlPoints(v);
    }
}




void LocalAssemblyGraph::computeControlPoints(vertex_descriptor v0)
{
    LocalAssemblyGraph& localAssemblyGraph = *this;
    LocalAssemblyGraphVertex& vertex0 = localAssemblyGraph[v0];
    SHASTA_ASSERT(vertex0.position.size() >= 2);
    const array<double, 2> vertex0Start = vertex0.position.front();
    const array<double, 2> vertex0End = vertex0.position.back();

    // Compute the "length" of this segment.
    const double dx0 = vertex0End[0] - vertex0Start[0];
    const double dy0 = vertex0End[1] - vertex0Start[1];
    const double d0 = sqrt(dx0 * dx0 + dy0 *dy0);

    // The distance of the control points from the ends of the segment.
    const double dq = 0.25 * d0;



    // To compute the first control point, q1,
    // average the unit vectors of the backward links.
    array<double, 2> direction = {0., 0.};
    uint64_t n = 0;
    BGL_FORALL_INEDGES(v0, e, localAssemblyGraph, LocalAssemblyGraph) {
        const vertex_descriptor v1 = source(e, localAssemblyGraph);
        LocalAssemblyGraphVertex& vertex1 = localAssemblyGraph[v1];
        SHASTA_ASSERT(vertex1.position.size() >= 2);
        const array<double, 2> vertex1Start = vertex1.position.front();

        const double dx = vertex1Start[0] - vertex0End[0];
        const double dy = vertex1Start[1] - vertex0End[1];
        const double d = sqrt(dx * dx + dy * dy);
        if(d == 0.) {
            continue;
        }

        // Accumulate the unit vector.
        ++n;
        direction[0] += dx / d;
        direction[1] += dy / d;
    }
    // Compute the average,normalized direction.
    double dLength = sqrt(direction[0] * direction[0] + direction[1] * direction[1]);
    direction[0] /= dLength;
    direction[1] /= dLength;

    // The first control point is in the opposite direction.
    vertex0.q1[0] = vertex0Start[0] - dq * direction[0];
    vertex0.q1[1] = vertex0Start[1] - dq * direction[1];



    // To compute the second control point, q2,
    // average the unit vectors of the forward links.
    direction = {0., 0.};
    n = 0;
    BGL_FORALL_OUTEDGES(v0, e, localAssemblyGraph, LocalAssemblyGraph) {
        const vertex_descriptor v1 = target(e, localAssemblyGraph);
        LocalAssemblyGraphVertex& vertex1 = localAssemblyGraph[v1];
        SHASTA_ASSERT(vertex1.position.size() >= 2);
        const array<double, 2> vertex1Start = vertex1.position.front();

        const double dx = vertex1Start[0] - vertex0End[0];
        const double dy = vertex1Start[1] - vertex0End[1];
        const double d = sqrt(dx * dx + dy * dy);
        if(d == 0.) {
            continue;
        }

        // Accumulate the unit vector.
        ++n;
        direction[0] += dx / d;
        direction[1] += dy / d;
    }
    // Compute the average,normalized direction.
    dLength = sqrt(direction[0] * direction[0] + direction[1] * direction[1]);
    direction[0] /= dLength;
    direction[1] /= dLength;

    // The second control point is in the oppotite direction.
    vertex0.q2[0] = vertex0End[0] - dq * direction[0];
    vertex0.q2[1] = vertex0End[1] - dq * direction[1];
}



// Return the svg color for a segment.
string LocalAssemblyGraph::segmentColor(uint64_t segmentId)
{
    const uint32_t hue = MurmurHash2(&segmentId, sizeof(segmentId), 231) % 360;
    return "hsl(" + to_string(hue) + ",50%,50%)";
}



// Find out if the paths of two segments are consecutive.
bool LocalAssemblyGraph::haveConsecutivePaths(
    vertex_descriptor v0,
    vertex_descriptor v1
) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;

    const LocalAssemblyGraphVertex& vertex0 = localAssemblyGraph[v0];
    const LocalAssemblyGraphVertex& vertex1 = localAssemblyGraph[v1];

    const uint64_t segmentId0 = vertex0.segmentId;
    const uint64_t segmentId1 = vertex1.segmentId;

    const auto path0 = assemblyGraph.paths[segmentId0];
    const auto path1 = assemblyGraph.paths[segmentId1];

    const MarkerGraphEdgeInfo& info0 = path0.back();
    const MarkerGraphEdgeInfo& info1 = path1.front();
    SHASTA_ASSERT(not info0.isVirtual);
    SHASTA_ASSERT(not info1.isVirtual);

    const MarkerGraph::EdgeId edgeId0 = info0.edgeId;
    const MarkerGraph::EdgeId edgeId1 = info1.edgeId;

    const MarkerGraph::Edge& edge0 = markerGraph.edges[edgeId0];
    const MarkerGraph::Edge& edge1 = markerGraph.edges[edgeId1];

    return edge0.target == edge1.source;
}



// Return the average link separation for the Link
// described by an edge.
double LocalAssemblyGraph::linkSeparation(edge_descriptor e) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;

    // Get the path length of the first segment.
    const vertex_descriptor v0 = source(e, localAssemblyGraph);
    const LocalAssemblyGraphVertex& vertex0 = localAssemblyGraph[v0];
    const uint64_t segmentId0 = vertex0.segmentId;
    const auto path0 = assemblyGraph.paths[segmentId0];
    const uint64_t pathLength0 = path0.size();

    // Now we can compute the link separation.
    const uint64_t linkId = localAssemblyGraph[e].linkId;
    return mode3::linkSeparation(assemblyGraph.transitions[linkId], pathLength0);

}



// Construct the svg options from an html request.
LocalAssemblyGraph::SvgOptions::SvgOptions(const vector<string>& request)
{
    HttpServer::getParameterValue(request, "pixelsPerUnitLength", pixelsPerUnitLength);
    HttpServer::getParameterValue(request, "layoutMethod", layoutMethod);

    // Segment length and thickness.
    HttpServer::getParameterValue(request, "minimumSegmentLength", minimumSegmentLength);
    HttpServer::getParameterValue(request, "additionalSegmentLengthPerMarker", additionalSegmentLengthPerMarker);
    HttpServer::getParameterValue(request, "segmentThickness", segmentThickness);

    // Link length and thickness.
    HttpServer::getParameterValue(request, "minimumLinkLength", minimumLinkLength);
    HttpServer::getParameterValue(request, "additionalLinkLengthPerMarker", additionalLinkLengthPerMarker);
    HttpServer::getParameterValue(request, "minimumLinkThickness", minimumLinkThickness);
    HttpServer::getParameterValue(request, "additionalLinkThicknessPerRead", additionalLinkThicknessPerRead);
}



// Add rows to the html request form.
void LocalAssemblyGraph::SvgOptions::addFormRows(ostream& html)
{
    html <<
        "<tr>"
        "<td>Graphics scaling factor"
        "<td class=centered><input type=text name=pixelsPerUnitLength size=8 style='text-align:center'"
        " value='" << pixelsPerUnitLength <<
        "'>"

        "<tr>"
        "<td>Graph layout method"
        "<td class=centered>"
        "<input type=radio name=layoutMethod value=neato"
        << (layoutMethod=="neato" ? " checked=checked" : "") <<
        ">neato<br>"
        "<input type=radio name=layoutMethod value=custom"
        << (layoutMethod=="custom" ? " checked=checked" : "") <<
        ">custom<br>"

        "<tr>"
        "<td>Segments"
        "<td class=centered>"
        "<table>"
        "<tr><td class=left>"
        "Minimum display length "
        "<td><input type=text name=minimumSegmentLength size=8 style='text-align:center'"
        " value='" << minimumSegmentLength << "'>"
        "<tr><td class=left>"
        "Additional display length per marker"
        "<td><input type=text name=additionalSegmentLengthPerMarker size=8 style='text-align:center'"
        " value='" << additionalSegmentLengthPerMarker << "'>"
        "<tr>"
        "<td class=left>Thickness"
        "<td class=centered><input type=text name=segmentThickness size=8 style='text-align:center'"
        " value='" << segmentThickness <<
        "'>"
        "</table>"



        "<tr>"
        "<td>Links"
        "<td class=centered>"
        "<table>"
        "<tr><td class=left>"
        "Minimum display length "
        "<td><input type=text name=minimumLinkLength size=8 style='text-align:center'"
        " value='" << minimumLinkLength << "'>"
        "<tr><td class=left>"
        "Additional display length per marker"
        "<td><input type=text name=additionalLinkLengthPerMarker size=8 style='text-align:center'"
        " value='" << additionalLinkLengthPerMarker << "'>"
        "<tr>"
        "<td class=left>Minimum thickness"
        "<td class=centered><input type=text name=minimumLinkThickness size=8 style='text-align:center'"
        " value='" << minimumLinkThickness <<
        "'>"
        "<tr>"
        "<td class=left>Additional thickness per read"
        "<td class=centered><input type=text name=additionalLinkThicknessPerRead size=8 style='text-align:center'"
        " value='" << additionalLinkThicknessPerRead <<
        "'>"
        "</table>"

        "</table>"


        ;

}



// Write the local assembly graph in gfa format.
void LocalAssemblyGraph::writeGfa(const string& fileName) const
{
    ofstream gfa(fileName);
    writeGfa(gfa);
}
void LocalAssemblyGraph::writeGfa(ostream& gfa) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;

    // Write the header.
    gfa << "H\tVN:Z:1.0\n";

    // Write the segments.
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const uint64_t segmentId = localAssemblyGraph[v].segmentId;
        const auto path = assemblyGraph.paths[segmentId];
        gfa <<
            "S\t" << segmentId << "\t" <<
            "*\tLN:i:" << path.size() << "\n";
    }


    // Write the links.
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const uint64_t linkId = localAssemblyGraph[e].linkId;
        const Link& link = assemblyGraph.links[linkId];
        gfa << "L\t" <<
            link.segmentId0 << "\t+\t" <<
            link.segmentId1 << "\t+\t0M\n";
    }

}



#endif
