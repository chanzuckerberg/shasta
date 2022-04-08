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
#include <boost/geometry/algorithms/make.hpp>
#include <boost/geometry/algorithms/length.hpp>
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


    // If necessary, compute a map containing a SegmentPairInformation object
    // containing pair information between the reference segment
    // and each segment in the local assembly graph.
    const bool doSegmentPairComputations = true;
    std::map<vertex_descriptor, mode3::AssemblyGraph::SegmentPairInformation> segmentPairInformationTable;
    mode3::AssemblyGraph::SegmentOrientedReadInformation referenceSegmentInfo;
    if(doSegmentPairComputations) {

        // Find oriented reads in the reference segment.
        assemblyGraph.getOrientedReadsOnSegment(options.referenceSegmentId, referenceSegmentInfo);

        // Loop over segments in the localAssemblyGraph.
        BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph){
            mode3::AssemblyGraph::SegmentOrientedReadInformation segmentInfo;
            assemblyGraph.getOrientedReadsOnSegment(
                localAssemblyGraph[v].segmentId, segmentInfo);

            mode3::AssemblyGraph::SegmentPairInformation segmentPairInformation;
            assemblyGraph.analyzeSegmentPair(
                options.referenceSegmentId, localAssemblyGraph[v].segmentId,
                referenceSegmentInfo, segmentInfo,
                assemblyGraph.markers, segmentPairInformation);

            segmentPairInformationTable.insert(make_pair(v, segmentPairInformation));
        }
    }



    using boost::geometry::add_point;
    using boost::geometry::expand;
    using boost::geometry::make_inverse;
    using boost::geometry::multiply_value;
    using boost::geometry::subtract_point;
    using Box = boost::geometry::model::box<Point>;

    // Compute the view box.
    Box box = make_inverse<Box>();
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex = localAssemblyGraph[v];
        SHASTA_ASSERT(vertex.position.size() >= 2);
        const Point& p1 = vertex.position.front();
        const Point& p2 = vertex.position.back();

        expand(box, p1);
        expand(box, p2);
    }
    Point minCorner = box.min_corner();
    Point maxCorner = box.max_corner();

    // Add a bit of extra space.
    Point delta = maxCorner;
    subtract_point(delta, minCorner);
    multiply_value(delta, 0.05);
    subtract_point(minCorner, delta);
    add_point(maxCorner, delta);



    // Figure out the required size of the viewbox.
    Point diagonal = maxCorner;
    subtract_point(diagonal, minCorner);

    // Begin the svg.
    const string svgId = "LocalAssemblyGraph";
    svg << "\n<svg id='" << svgId <<
        "' width='" <<  options.sizePixels <<
        "' height='" << options.sizePixels <<
        "' viewbox='" << minCorner.x() << " " << minCorner.y() << " " <<
        diagonal.x() << " " <<
        diagonal.y() << "'"
        " style='border-style:solid;border-color:Black;'"
        ">\n";



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
        const Point& p1 = vertex1.position.back();
        const Point& p2 = vertex2.position.front();
        const double length = boost::geometry::distance(p1, p2);

        // Get the tangents and compute the control points.
        const double controlPointDistance = 0.25 * length;
        const Point& t1 = vertex1.t2;
        const Point& t2 = vertex2.t1;
        Point q1 = t1;
        multiply_value(q1, controlPointDistance);
        add_point(q1, p1);
        Point q2 = t2;
        multiply_value(q2, controlPointDistance);
        add_point(q2, p2);

        const double linkThickness =
            options.minimumLinkThickness +
            options.additionalLinkThicknessPerRead * double(assemblyGraph.linkCoverage(linkId) - 1);

        const string dash =
            areConsecutivePaths ? "" :
            " stroke-dasharray='0 " + to_string(1.5 * linkThickness) + "'";


        svg <<
            "<g>"
            // "<a href='exploreMode3AssemblyGraphLink?linkId=" << linkId << "'>"
            "<title>"
            "Link " << linkId <<
            " from segment " << segmentId1 <<
            " to segment " << segmentId2 <<
            ", coverage " << assemblyGraph.linkCoverage(linkId) <<
            "</title>"
            "<path d="
            "'M " << p1.x() << " " << p1.y() <<
            " C " << q1.x() << " " << q1.y() << ", "
                  << q2.x() << " " << q2.y() << ","
                  << p2.x() << " " << p2.y() << "'"
            " stroke='" << options.linkColor << "'" <<
            dash <<
            " stroke-width='" << linkThickness << "'"
            " stroke-linecap='round'"
            " fill='transparent'"
            // " vector-effect='non-scaling-stroke'"
            " onclick='if(event.ctrlKey) {location.href=\"exploreMode3AssemblyGraphLink?linkId=" << linkId << "\";}'"
            "/>"
            // "</a>"
            "</g>\n";

    }
    svg << "</g>\n";



    // Write the segments.
    svg << "<g id='" << svgId << "-segments'>\n";
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        const LocalAssemblyGraphVertex& vertex = localAssemblyGraph[v];
        const uint64_t distance = localAssemblyGraph[v].distance;

        // Get the positions of the ends of this segment.
        SHASTA_ASSERT(vertex.position.size() >= 2);
        const Point& p1 = vertex.position.front();
        const Point& p2 = vertex.position.back();
        const double length = boost::geometry::distance(p1, p2);

        // Get the tangents and compute the control points.
        const double controlPointDistance = 0.25 * length;
        const Point& t1 = vertex.t1;
        const Point& t2 = vertex.t2;
        Point q1 = t1;
        multiply_value(q1, -controlPointDistance);
        add_point(q1, p1);
        Point q2 = t2;
        multiply_value(q2, -controlPointDistance);
        add_point(q2, p2);

        const uint64_t segmentId = localAssemblyGraph[v].segmentId;



        // Decide the color for this segment.
        string color;
        if(distance == maxDistance) {
            color = options.segmentAtMaxDistanceColor;
        } else {
            if(options.segmentColoring == "random") {
                color = randomSegmentColor(segmentId);
            } else if(options.segmentColoring == "uniform") {
                color = options.segmentColor;
            } else if(options.segmentColoring == "byCommonReads") {
                const uint64_t commonReadCount = segmentPairInformationTable[v].commonOrientedReadCount;
                double fraction;
                if(options.greenThreshold) {
                    fraction = min(1., double(commonReadCount) / double(options.greenThreshold));
                } else {
                    fraction = double(commonReadCount) / double(referenceSegmentInfo.infos.size());
                }
                const uint64_t hue = uint64_t(std::round(fraction * 120.));
                color = "hsl(" + to_string(hue) + ",100%, 50%)";
            } else if(options.segmentColoring == "byMissingFractionOnDisplayedSegment") {
                const double fraction = 1. - segmentPairInformationTable[v].missingFraction(1);
                const uint64_t hue = uint64_t(std::round(fraction * 120.));
                color = "hsl(" + to_string(hue) + ",100%, 50%)";
            } else if(options.segmentColoring == "byMissingFractionOnReferenceSegment") {
                const double fraction = 1. - segmentPairInformationTable[v].missingFraction(0);
                const uint64_t hue = uint64_t(std::round(fraction * 120.));
                color = "hsl(" + to_string(hue) + ",100%, 50%)";
            } else {
                color = "Black";
            }
        }



        // Get the oriented reads and average edge coverage.
        vector<OrientedReadId> orientedReadIds;
        const double averageEdgeCoverage = assemblyGraph.findOrientedReadsOnSegment(segmentId, orientedReadIds);

        // Create a marker to show the arrow for this segment.
        const string arrowMarkerName = "arrow" + to_string(segmentId);
        svg <<
            "<defs>\n"
            "<marker id='" << arrowMarkerName <<
            "' viewBox='0 0 0.6 1'\n"
            "refX='0.1' refY='0.5'\n"
            "markerUnits='strokeWidth'\n"
            "markerWidth='0.6' markerHeight='1'\n"
            "orient='auto'>\n"
            "<path d='M 0 0 L 0.1 0 L 0.6 0.5 L 0.1 1 L 0 1 z' "
            "fill='" << color << "' "
            "/>\n"
            "</marker>\n"
            "</defs>\n";

        // Add this segment to the svg.
        const auto& segmentPairInfo = segmentPairInformationTable[v];
        const auto oldPrecision = svg.precision(1);
        const auto oldFlags = svg.setf(std::ios_base::fixed, std::ios_base::floatfield);
        /*
        svg <<
            "<g>"
            // "<a href='exploreMode3AssemblyGraphSegment?segmentId=" << segmentId << "'>"
            "<title>"
            "Segment " << segmentId <<
            ", distance from start segment " << distance <<
            ", path length " << assemblyGraph.paths.size(segmentId) <<
            ", average marker graph edge coverage " << averageEdgeCoverage <<
            ", number of distinct oriented reads " << orientedReadIds.size();
        if(doSegmentPairComputations) {
            svg << ", number of common oriented reads " << segmentPairInfo.commonOrientedReadCount <<
                " of " << referenceSegmentInfo.infos.size();
        }
        */
        svg <<
            // "</title>"
            "<path id='Segment-" << segmentId << "'"
            " onmouseenter='onMouseEnterSegment(" <<
            segmentId << "," <<
            distance << "," <<
            assemblyGraph.paths.size(segmentId) << "," <<
            averageEdgeCoverage << "," <<
            orientedReadIds.size() << "," <<
            segmentPairInfo.commonOrientedReadCount << "," <<
            segmentPairInfo.tooShortCount << "," <<
            segmentPairInfo.missingOrientedReadCount[0] << "," <<
            segmentPairInfo.missingOrientedReadCount[1] << ")'" <<
            " onmouseleave='onMouseExitSegment()'" <<
#if 0
            common, tooShort, missingFromReference, missingFromDisplayed


            "' d='M " <<
            p1.x() << " " << p1.y() << " C " <<
            q1.x() << " " << q1.y() << ", " <<
            q2.x() << " " << q2.y() << ", " <<
            p2.x() << " " << p2.y() << "'" <<
#endif
            " d='M " <<
            p1.x() << " " << p1.y() << " L " <<
            p2.x() << " " << p2.y() << "'" <<
            " stroke='" << color << "'"
            " stroke-width='" <<
            options.minimumSegmentThickness + averageEdgeCoverage * options.additionalSegmentThicknessPerUnitCoverage << "'"
            " fill='none'"
            " marker-end='url(#" <<
            arrowMarkerName <<
            ")'"
            " onclick='if(event.ctrlKey) {location.href=\"exploreMode3AssemblyGraphSegment?segmentId=" << segmentId << "\";}'"
            "/>"
            // "</a>"
            // "</g>"
            "\n";
        svg.precision(oldPrecision);
        svg.flags(oldFlags);
    }
    svg << "</g>\n";



    // End the svg.
    svg << "</svg>\n";
}



void mode3::LocalAssemblyGraph::computeLayout(
    const SvgOptions& options,
    double timeout)
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
    ComputeLayoutReturnCode returnCode = ComputeLayoutReturnCode::Success;
    if(options.layoutMethod == "neato") {
        returnCode = shasta::computeLayoutGraphviz(g, "neato", timeout, positionMap, "", &edgeLengthMap);
    } else if(options.layoutMethod == "custom") {
        returnCode = shasta::computeLayoutCustom(g, edgeLengthMap, positionMap, timeout);
    } else {
        throw runtime_error("Invalid layout method specified: " + options.layoutMethod);
    }
    if(returnCode == ComputeLayoutReturnCode::Timeout) {
        throw runtime_error("Graph layout took too long. "
            "Increase the timeout or decrease the maximum distance.");
    }
    if(returnCode != ComputeLayoutReturnCode::Success) {
        throw runtime_error("Graph layout failed.");
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
            const array<double, 2>& p = jt->second;
            vertex.position.push_back(Point(p[0], p[1]));
        }
    }
}



void LocalAssemblyGraph::computeSegmentTangents()
{
    LocalAssemblyGraph& localAssemblyGraph = *this;
    BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
        computeSegmentTangents(v);
    }
}




void LocalAssemblyGraph::computeSegmentTangents(vertex_descriptor v0)
{
    LocalAssemblyGraph& localAssemblyGraph = *this;
    LocalAssemblyGraphVertex& vertex0 = localAssemblyGraph[v0];
    SHASTA_ASSERT(vertex0.position.size() >= 2);
    const Point& vertex0Start = vertex0.position.front();
    const Point& vertex0End = vertex0.position.back();

    Point t = vertex0End;
    boost::geometry::subtract_point(t, vertex0Start);
    const double length = sqrt(t.x() * t.x() + t.y() * t.y());
    boost::geometry::multiply_value(t, 1. / length);
    vertex0.t2 = t;
    boost::geometry::multiply_value(t, -1.);
    vertex0.t1 = t;


#if 0
    // This is used if we display segments as Bezier cubics.


    // To compute t1, average the unit vectors of the backward links.
    array<double, 2> direction = {0., 0.};
    uint64_t n = 0;
    BGL_FORALL_INEDGES(v0, e, localAssemblyGraph, LocalAssemblyGraph) {
        const vertex_descriptor v1 = source(e, localAssemblyGraph);
        LocalAssemblyGraphVertex& vertex1 = localAssemblyGraph[v1];
        SHASTA_ASSERT(vertex1.position.size() >= 2);
        const Point& vertex1Start = vertex1.position.front();

        const double dx = vertex1Start.x() - vertex0End.x();
        const double dy = vertex1Start.y() - vertex0End.y();
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
    if(dLength == 0.) {
        direction[0] = vertex0Start.x() - vertex0End.x();
        direction[1] = vertex0Start.y() - vertex0End.y();
        dLength = sqrt(direction[0] * direction[0] + direction[1] * direction[1]);
    }
    direction[0] /= dLength;
    direction[1] /= dLength;

    vertex0.t1.x(direction[0]);
    vertex0.t1.y(direction[1]);



    // To compute the second control point, q2,
    // average the unit vectors of the forward links.
    direction = {0., 0.};
    n = 0;
    BGL_FORALL_OUTEDGES(v0, e, localAssemblyGraph, LocalAssemblyGraph) {
        const vertex_descriptor v1 = target(e, localAssemblyGraph);
        LocalAssemblyGraphVertex& vertex1 = localAssemblyGraph[v1];
        SHASTA_ASSERT(vertex1.position.size() >= 2);
        const Point& vertex1Start = vertex1.position.front();

        const double dx = vertex1Start.x() - vertex0End.x();
        const double dy = vertex1Start.y() - vertex0End.y();
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
    if(dLength == 0.) {
        direction[0] = vertex0End.x() - vertex0Start.x();
        direction[1] = vertex0End.y() - vertex0Start.y();
        dLength = sqrt(direction[0] * direction[0] + direction[1] * direction[1]);
    }
    direction[0] /= dLength;
    direction[1] /= dLength;

    vertex0.t2.x(direction[0]);
    vertex0.t2.y(direction[1]);
#endif
}



// Return the svg color for a segment.
string LocalAssemblyGraph::randomSegmentColor(uint64_t segmentId)
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
    HttpServer::getParameterValue(request, "sizePixels", sizePixels);
    HttpServer::getParameterValue(request, "layoutMethod", layoutMethod);

    // Segment length and thickness.
    HttpServer::getParameterValue(request, "minimumSegmentLength", minimumSegmentLength);
    HttpServer::getParameterValue(request, "additionalSegmentLengthPerMarker", additionalSegmentLengthPerMarker);
    HttpServer::getParameterValue(request, "minimumSegmentThickness", minimumSegmentThickness);
    HttpServer::getParameterValue(request, "additionalSegmentThicknessPerUnitCoverage", additionalSegmentThicknessPerUnitCoverage);

    // Segment coloring
    HttpServer::getParameterValue(request, "segmentColoring", segmentColoring);
    HttpServer::getParameterValue(request, "segmentColor", segmentColor);
    HttpServer::getParameterValue(request, "greenThreshold", greenThreshold);
    HttpServer::getParameterValue(request, "referenceSegmentId", referenceSegmentId);

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
        "<td>Graphics size in pixels"
        "<td class=centered><input type=text name=sizePixels size=8 style='text-align:center'"
        " value='" << sizePixels <<
        "'>"

        "<tr>"
        "<td>Graph layout method"
        "<td class=left>"
        "<input type=radio name=layoutMethod value=neato"
        << (layoutMethod=="neato" ? " checked=checked" : "") <<
        ">Graphviz neato (slow for large graphs)<br>"
        "<input type=radio name=layoutMethod value=custom"
        << (layoutMethod=="custom" ? " checked=checked" : "") <<
        ">Custom (user-provided command <code>customLayout</code>)<br>"

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
        "<td class=left>Minimum thickness"
        "<td class=centered><input type=text name=minimumSegmentThickness size=8 style='text-align:center'"
        " value='" << minimumSegmentThickness <<
        "'>"
        "<tr>"
        "<td class=left>Additional thickness per unit coverage"
        "<td class=centered><input type=text name=additionalSegmentThicknessPerUnitCoverage size=8 style='text-align:center'"
        " value='" << additionalSegmentThicknessPerUnitCoverage <<
        "'>"



        // Segment coloring.
        "<tr>"
        "<td class = left>Color"
        "<td class=left>"

        // Random segment coloring.
        "<input type=radio name=segmentColoring value=random"
        << (segmentColoring=="random" ? " checked=checked" : "") <<
        ">Random<br>"

        // Uniform segment coloring.
        "<input type=radio name=segmentColoring value=uniform"
        << (segmentColoring=="uniform" ? " checked=checked" : "") <<
        ">"
        "<input type=text name=segmentColor size=8 style='text-align:center'"
                " value='" << segmentColor << "'>"
        "<br>"

        // Segment coloring by number of common reads with the reference segment.
        "<input type=radio name=segmentColoring value=byCommonReads"
        << (segmentColoring=="byCommonReads" ? " checked=checked" : "") <<
        ">By number of common supporting oriented reads with reference segment"
        "<div style='text-indent:3em'>"
        "Green if at least "
        "<input type=text name=greenThreshold size=4 style='text-align:center'"
        " value='" << greenThreshold <<
        "'>"        " common reads (0 = automatic)"
        "</div>"

        // Segment coloring by missing fraction on the displayed segment versus the reference segment.
        "<input type=radio name=segmentColoring value=byMissingFractionOnDisplayedSegment"
        << (segmentColoring=="byMissingFractionOnDisplayedSegment" ? " checked=checked" : "") <<
        ">By missing fraction on the displayed segment versus the reference segment"
        "<br>"

        // Segment coloring by missing fraction on the reference segment versus the displayed segment.
        "<input type=radio name=segmentColoring value=byMissingFractionOnReferenceSegment"
        << (segmentColoring=="byMissingFractionOnReferenceSegment" ? " checked=checked" : "") <<
        ">By missing fraction on the reference segment versus the displayed segment"
        "<br>"

        "Reference segment&nbsp;<input type=text name=referenceSegmentId size=8 style='text-align:center'"
                " value='" << referenceSegmentId << "'>"

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
