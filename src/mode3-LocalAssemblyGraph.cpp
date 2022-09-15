// Shasta.
#include "mode3-LocalAssemblyGraph.hpp"
#include "mode3-AssemblyPath.hpp"
#include "mode3-SegmentPairInformation.hpp"
#include "computeLayout.hpp"
#include "html.hpp"
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
            const mode3::AssemblyGraph::Link& link = assemblyGraph.links[linkId];
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
            const mode3::AssemblyGraph::Link& link = assemblyGraph.links[linkId];
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
            const mode3::AssemblyGraph::Link& link = assemblyGraph.links[linkId];
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



void mode3::LocalAssemblyGraph::writeHtml(ostream& html, const SvgOptions& options) const
{
    // Write the svg object.
    html << "<div style='display: inline-block; vertical-align:top'>";
    vector<mode3::AssemblyGraph::AnalyzeSubgraphClasses::Cluster> clusters;
    writeSvg(html, options, clusters);
    html << "</div>";
    addSvgDragAndZoom(html);

    // Side panel.
    html << "<div style='display: inline-block'>";



    // Highlight a segment.
    html << R"stringDelimiter(
        <script>
        function highlightSegment()
        {
            // Get the segment id from the input field.
            inputField = document.getElementById("highlightInputField");
            segmentId = inputField.value;
            inputField.value = "";

            // Make it dashed and wider.
            var element = document.getElementById("Segment-" + segmentId);
            var thickness = element.getAttribute("stroke-width");
            element.style.strokeDasharray = 0.2 * thickness;
            element.setAttribute("stroke-width", 2. * thickness);
        }
        </script>
        Highlight segment
        <input id=highlightInputField type=text onchange="highlightSegment()" size=10>
        )stringDelimiter";



    // Zoom to a segment.
    html << R"stringDelimiter(
        <script>
        function zoomToSegment()
        {
            // Get the segment id from the input field.
            inputField = document.getElementById("zoomInputField");
            segmentId = inputField.value;
            inputField.value = "";

            zoomToGivenSegment(segmentId);
        }

        function zoomToGivenSegment(segmentId)
        {

            // Find the bounding box and its center.
            var element = document.getElementById("Segment-" + segmentId);
            var box = element.getBBox();
            var xCenter = box.x + 0.5 * box.width;
            var yCenter = box.y + 0.5 * box.height;

            // Change the viewbox of the svg to be a bit larger than a square
            // containing the bounding box.
            var enlargeFactor = 5.;
            var size = enlargeFactor * Math.max(box.width, box.height);
            width = size;
            height = size;
            x = xCenter - 0.5 * size;
            y = yCenter - 0.5 * size;
            var svg = document.querySelector('svg');
            svg.setAttribute('viewBox', `${x} ${y} ${size} ${size}`);
            ratio = size / svg.getBoundingClientRect().width;

        }
        </script>
        <p>Zoom to segment
        <input id=zoomInputField type=text onchange="zoomToSegment()" size=10>
        )stringDelimiter";



    // Initial zoom to segment of interest.
    if(options.segmentColoring == "path") {
        html << "\n<script>zoomToGivenSegment(" << options.pathStart << ");</script>\n";
    }
    if(
        options.segmentColoring == "byCommonReads" or
        options.segmentColoring == "byJaccard" or
        options.segmentColoring == "byUnexplainedFractionOnReferenceSegment" or
        options.segmentColoring == "byUnexplainedFractionOnDisplayedSegment"
        ) {
        html << "\n<script>zoomToGivenSegment(" << options.referenceSegmentId << ");</script>\n";
    }



    // Tables that will be automatically updated when the mouse is on a segment.
    html << R"zzz(  
<p>
Hover on a segment to populate the tables below.
<p>
<table style='font-size:9'>
<tr><th class='left'>Segment id<td id='segmentIdCell' class=centered style='width:8em'>
<tr><th class='left'>Distance from start segment<td id='distanceCell' class=centered style='width:8em'>
<tr><th class='left'>Path length<td id='pathLengthCell' class=centered style='width:8em'>
<tr><th class='left'>Average edge coverage<td id='coverageCell' class=centered style='width:8em'>
<tr><th class='left'>Cluster id<td id='clusterIdCell' class=centered style='width:8em'>
</table>
<p>
Comparison of read compositions
<p>
<table>

<tr>
<td>
<th>Reference<br>segment
<th>Displayed<br>segment

<tr>
<th class='left'>Total
<th id='totalReferenceCell'>
<th id='totalDisplayedCell'>

<tr>
<th class='left'>Common
<th id='commonReferenceCell'>
<th id='commonDisplayedCell'>

<tr>
<th class='left'>Short
<th id='shortReferenceCell'>
<th id='shortDisplayedCell'>

<tr>
<th class='left'>Jaccard
<th id='jaccardReferenceCell'>
<th id='jaccardDisplayedCell'>

<tr>
<th class='left'>Unexplained
<th id='unexplainedReferenceCell'>
<th id='unexplainedDisplayedCell'>

<tr>
<th class='left'>Unexplained fraction
<th id='unexplainedFractionReferenceCell'>
<th id='unexplainedFractionDisplayedCell'>

</table>

<script>
function onMouseEnterSegment(id, distance, pathLength, coverage, clusterId,
    totalReference, totalDisplayed,
    shortReference, shortDisplayed,
    common, 
    unexplainedReference, unexplainedDisplayed)
{
    document.getElementById('segmentIdCell').innerHTML = id;
    document.getElementById('distanceCell').innerHTML = distance;
    document.getElementById('pathLengthCell').innerHTML = pathLength;
    document.getElementById('coverageCell').innerHTML = coverage;
    if(clusterId != 18446744073709551615) {
        document.getElementById('clusterIdCell').innerHTML = clusterId;
    }

    document.getElementById('totalReferenceCell').innerHTML = totalReference;
    document.getElementById('totalDisplayedCell').innerHTML = totalDisplayed;
    document.getElementById('commonReferenceCell').innerHTML = common;
    document.getElementById('commonDisplayedCell').innerHTML = common;

    if(common > 0) {
        document.getElementById('shortReferenceCell').innerHTML = shortReference;
        document.getElementById('shortDisplayedCell').innerHTML = shortDisplayed;
        jaccard = (common / (common + unexplainedReference + unexplainedDisplayed)).toFixed(2)
        document.getElementById('jaccardReferenceCell').innerHTML = jaccard;
        document.getElementById('jaccardDisplayedCell').innerHTML = jaccard;
        document.getElementById('unexplainedReferenceCell').innerHTML = unexplainedReference;
        document.getElementById('unexplainedDisplayedCell').innerHTML = unexplainedDisplayed;
        document.getElementById('unexplainedFractionReferenceCell').innerHTML = 
            (unexplainedReference / (common + unexplainedReference)).toFixed(2);
        document.getElementById('unexplainedFractionDisplayedCell').innerHTML = 
            (unexplainedDisplayed / (common + unexplainedDisplayed)).toFixed(2);
    }   
}
function onMouseExitSegment()
{
    document.getElementById('segmentIdCell').innerHTML = '';
    document.getElementById('distanceCell').innerHTML = '';
    document.getElementById('pathLengthCell').innerHTML = '';
    document.getElementById('coverageCell').innerHTML = '';
    document.getElementById('clusterIdCell').innerHTML = '';

    document.getElementById('totalReferenceCell').innerHTML = '';
    document.getElementById('totalDisplayedCell').innerHTML = '';
    document.getElementById('shortReferenceCell').innerHTML = '';
    document.getElementById('shortDisplayedCell').innerHTML = '';
    document.getElementById('commonReferenceCell').innerHTML = '';
    document.getElementById('commonDisplayedCell').innerHTML = '';
    document.getElementById('jaccardReferenceCell').innerHTML = '';
    document.getElementById('jaccardDisplayedCell').innerHTML = '';
    document.getElementById('unexplainedReferenceCell').innerHTML = '';
    document.getElementById('unexplainedDisplayedCell').innerHTML = '';
    document.getElementById('unexplainedFractionReferenceCell').innerHTML = '';
    document.getElementById('unexplainedFractionDisplayedCell').innerHTML = '';
}
</script>
    )zzz";



    // Change segment thickness
    html << R"stringDelimiter(
    <p><table>
    <tr><th class=left>Segment thickness<td>
    <button type='button' onClick='segmentThickness(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='segmentThickness(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='segmentThickness(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='segmentThickness(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='segmentThickness(2.)' style='width:3em'>++</button>
    <button type='button' onClick='segmentThickness(10.)' style='width:3em'>+++</button>
        <script>
        function segmentThickness(factor)
        {
            const group = document.getElementById('LocalAssemblyGraph-segments');
            descendants = group.querySelectorAll("path");
            for (let i=0; i<descendants.length; i++) {
                path = descendants[i];
                path.setAttribute('stroke-width', factor * path.getAttribute('stroke-width'));
            }
        }
        </script>
        )stringDelimiter";



    // Change link thickness
    html << R"stringDelimiter(
    <tr><th class=left>Link thickness<td>
    <button type='button' onClick='linkThickness(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='linkThickness(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='linkThickness(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='linkThickness(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='linkThickness(2.)' style='width:3em'>++</button>
    <button type='button' onClick='linkThickness(10.)' style='width:3em'>+++</button>
        <script>
        function linkThickness(factor)
        {
            const group1 = document.getElementById('LocalAssemblyGraph-links');
            for (let i=0; i<group1.children.length; i++) {
                group2 = group1.children[i];
                if(group2.tagName == 'g') {
                    for (let j=0; j<group2.children.length; j++) {
                        path = group2.children[j];
                        if(path.tagName == 'path') {
                            path.setAttribute('stroke-width', factor * path.getAttribute('stroke-width'));
                        }
                    }
                }
            }
        }
        </script>
        )stringDelimiter";



    // Zoom buttons.
    html << R"stringDelimiter(
    <tr title='Or use the mouse wheel.'><th class=left>Zoom<td>
    <button type='button' onClick='zoomSvg(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='zoomSvg(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='zoomSvg(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='zoomSvg(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='zoomSvg(2.)' style='width:3em'>++</button>
    <button type='button' onClick='zoomSvg(10.)' style='width:3em'>+++</button>
     </table>
        )stringDelimiter";


    // Code to display one local cluster at a time, with a button
    // to cycle through them.
    if(options.segmentColoring == "byLocalCluster") {
        html <<
            "<br>Found " << clusters.size() << " clusters. "
            "Displaying cluster <span id='currentCluster'></span>"
            "<br><button onClick='previousCluster()'>Previous<br>cluster</button>"
            "<button onClick='nextCluster()'>Next<br>cluster</button>"
            "<script>\n"
            "var clusters = [";
        for(uint64_t i=0; i<clusters.size(); i++) {
            html << "[";
            const auto & cluster = clusters[i];
            for(uint64_t j=0; j<cluster.segments.size(); j++) {
                html << cluster.segments[j].first;
                if(j != cluster.segments.size() - 1) {
                    html << ",";
                }
            }
            html << "]";
            if(i != clusters.size() -1) {
                html << ",";
            }
        }
        html << "];\n";

        html << R"stringDelimiter(

        function clusterColor(clusterId)
        {
            var ratio = clusterId / clusters.length;
            return 'hsl(' + Math.round(360*ratio) + ', 85%, 70%)';
        }

        function highlightCluster(clusterId, color)
        {
            for(i=0; i<clusters[clusterId].length; i++) {
                segmentId = clusters[clusterId][i];
                document.getElementById("Segment-" + segmentId).style.stroke = color;
                document.getElementById("marker" + segmentId).style.fill = color;
            }
        }
        var currentCluster = 0;
        highlightCluster(currentCluster, clusterColor(currentCluster));
        document.getElementById("currentCluster").innerHTML = currentCluster;
        function nextCluster()
        {
            highlightCluster(currentCluster, "Black");
            currentCluster = currentCluster + 1;
            if(currentCluster == clusters.length) {
                currentCluster = 0;
            }
            highlightCluster(currentCluster, clusterColor(currentCluster));
            document.getElementById("currentCluster").innerHTML = currentCluster;
        }
        function previousCluster()
        {
            highlightCluster(currentCluster, "Black");
            if(currentCluster == 0) {
                currentCluster = clusters.length;
            }
            currentCluster = currentCluster - 1;
            highlightCluster(currentCluster, clusterColor(currentCluster));
            document.getElementById("currentCluster").innerHTML = currentCluster;
        }
        </script>

        )stringDelimiter";
    }

    // End of side panel.
    html << "</div>";

}



void mode3::LocalAssemblyGraph::writeSvg(
    const string& fileName,
    const SvgOptions& options,
    vector<mode3::AssemblyGraph::AnalyzeSubgraphClasses::Cluster>& clusters) const
{
    ofstream svg(fileName);
    writeSvg(svg, options, clusters);
}
void mode3::LocalAssemblyGraph::writeSvg(
    ostream& svg,
    const SvgOptions& options,
    vector<mode3::AssemblyGraph::AnalyzeSubgraphClasses::Cluster>& clusters
    ) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;


    // If necessary, compute a map containing a SegmentPairInformation object
    // containing pair information between the reference segment
    // and each segment in the local assembly graph.
    const bool doSegmentPairComputations = true;
    std::map<vertex_descriptor, SegmentPairInformation> segmentPairInformationTable;
    mode3::AssemblyGraph::SegmentOrientedReadInformation referenceSegmentInfo;
    if(doSegmentPairComputations) {

        // Find oriented reads in the reference segment.
        assemblyGraph.getOrientedReadsOnSegment(options.referenceSegmentId, referenceSegmentInfo);

        // Loop over segments in the localAssemblyGraph.
        BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph){
            mode3::AssemblyGraph::SegmentOrientedReadInformation segmentInfo;
            assemblyGraph.getOrientedReadsOnSegment(
                localAssemblyGraph[v].segmentId, segmentInfo);

            SegmentPairInformation segmentPairInformation;
            assemblyGraph.analyzeSegmentPair(
                options.referenceSegmentId, localAssemblyGraph[v].segmentId,
                referenceSegmentInfo, segmentInfo,
                assemblyGraph.markers, segmentPairInformation);

            segmentPairInformationTable.insert(make_pair(v, segmentPairInformation));
        }
    }


    std::map<uint64_t, vector<pair<uint64_t, bool> > > pathSegments; // map(segmentId, (positionsInPath, is referenceSegment)).
    AssemblyPath path;
    if(options.segmentColoring == "path") {
        if(options.pathDirection=="forward" or options.pathDirection=="backward") {
            // Forward or backward.
            assemblyGraph.createAssemblyPath(options.pathStart,
                (options.pathDirection == "forward") ? 0 : 1, path);
            if(options.pathDirection == "backward") {
                reverse(path.segments.begin(), path.segments.end());
            }
        } else {
            // Bidirectional.
            AssemblyPath forwardPath;
            AssemblyPath backwardPath;
            assemblyGraph.createAssemblyPath(options.pathStart, 0, forwardPath);
            assemblyGraph.createAssemblyPath(options.pathStart, 1, backwardPath);
            // Stitch them together, making sure not to repeat the starting segment.
            path.segments.clear();
            copy(backwardPath.segments.rbegin(), backwardPath.segments.rend(), back_inserter(path.segments));
            copy(forwardPath.segments.begin() + 1, forwardPath.segments.end(), back_inserter(path.segments));
        }
        for(uint64_t position=0; position<path.segments.size(); position++) {
            const AssemblyPathSegment& segment = path.segments[position];
            const uint64_t segmentId = segment.id;
            pathSegments[segmentId].push_back(make_pair(position, segment.isPrimary));
        }
        svg << "\nPath of length " << path.segments.size() << " starting at segment " << path.segments.front().id <<
            " and ending at segment " << path.segments.back().id << "<br>";

        ofstream csv("Path.csv");
        csv << "Position,SegmentId,Reference\n";
        for(uint64_t position=0; position<path.segments.size(); position++) {
            const AssemblyPathSegment& segment = path.segments[position];
            csv << position << "," << segment.id << "," << int(segment.isPrimary) << "\n";
        }

        // If requested, assemble path sequence.
        if(options.assemblePathSequence) {
            path.assemble(assemblyGraph);
        }
    }



    // If coloring by local cluster, call mode3::AssemblyGraph::analyzeSubgraph,
    // passing as input all the segments in the LocalAssemblyGraph
    // except those at maximum distance.
    if(options.segmentColoring == "byLocalCluster") {
        vector<uint64_t> segmentIds;
        BGL_FORALL_VERTICES(v, localAssemblyGraph, LocalAssemblyGraph) {
            const LocalAssemblyGraphVertex& vertex = localAssemblyGraph[v];
            if(vertex.distance != maxDistance) {
                segmentIds.push_back(vertex.segmentId);
            }
        }
        assemblyGraph.analyzeSubgraph(segmentIds, clusters, true);

    }



    // If coloring by cluster id only some clusters, create a color map
    // for the clusters to be colored.
    std::map<uint64_t, string> clusterColorMap;
    if(options.segmentColoring == "byCluster") {
        const uint64_t clusterCount = options.clustersToBeColored.size();
        if(clusterCount > 0) {
            for(uint64_t i=0; i<clusterCount; i++) {
                const uint64_t hue = uint64_t(std::round(double(i) * 360. / double(clusterCount)));
                const string color = "hsl(" + to_string(hue) + ",100%, 50%)";
                const uint64_t clusterId = options.clustersToBeColored[i];
                clusterColorMap.insert(make_pair(clusterId, color));
            }
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
        const AssemblyGraph::Link& link = assemblyGraph.links[linkId];

        // Access the LocalAssemblyGraph vertices corresponding to
        // the two segments of this Link and extract some information
        // from them.
        const vertex_descriptor v1 = source(e, localAssemblyGraph);
        const vertex_descriptor v2 = target(e, localAssemblyGraph);
        const LocalAssemblyGraphVertex& vertex1 = localAssemblyGraph[v1];
        const LocalAssemblyGraphVertex& vertex2 = localAssemblyGraph[v2];
        const uint64_t segmentId1 = vertex1.segmentId;
        const uint64_t segmentId2 = vertex2.segmentId;

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
            link.segmentsAreAdjacent ? "" :
            " stroke-dasharray='0 " + to_string(1.5 * linkThickness) + "'";

        // If the link participates in a path, color it consistently with the
        // segments is joins.
        string linkColor = options.linkColor;
        if(options.segmentColoring == "path") {
            const auto it1 = pathSegments.find(segmentId1);
            if(it1 != pathSegments.end()) {
                const auto positions1 = it1->second;
                SHASTA_ASSERT(not positions1.empty());
                const auto it2 = pathSegments.find(segmentId2);
                if(it2 != pathSegments.end()) {
                    const auto positions2 = it2->second;
                    SHASTA_ASSERT(not positions2.empty());
                    if(positions1.size()==1 and positions2.size()==1) {
                        const uint64_t position1 = positions1.front().first;
                        const uint64_t position2 = positions2.front().first;
                        if(position2 == position1 + 1) {
                            const uint32_t hue = uint32_t(
                                std::round(120. * double(position1 + position2) / double(path.segments.size())));
                            linkColor = "hsl(" + to_string(hue) + ",100%, 20%)";
                        }
                    } else {
                        linkColor = "Fuchsia";
                    }
                }
            }
        }

        svg <<
            "<g>"
            // "<a href='exploreMode3AssemblyGraphLink?linkId=" << linkId << "'>"
            "<title>"
            "Link " << linkId <<
            " from segment " << segmentId1 <<
            " to segment " << segmentId2 <<
            ", coverage " << assemblyGraph.linkCoverage(linkId) <<
            ", separation " << link.separation <<
            "</title>"
            "<path d="
            "'M " << p1.x() << " " << p1.y() <<
            " C " << q1.x() << " " << q1.y() << ", "
                  << q2.x() << " " << q2.y() << ","
                  << p2.x() << " " << p2.y() << "'"
            " stroke='" << linkColor << "'" <<
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
                const uint64_t commonCount = segmentPairInformationTable[v].commonCount;
                double fraction;
                if(options.greenThreshold) {
                    fraction = min(1., double(commonCount) / double(options.greenThreshold));
                } else {
                    fraction = double(commonCount) / double(referenceSegmentInfo.infos.size());
                }
                const uint64_t hue = uint64_t(std::round(fraction * 120.));
                color = "hsl(" + to_string(hue) + ",100%, 50%)";
            } else if(options.segmentColoring == "byJaccard") {
                const auto& pairInfo = segmentPairInformationTable[v];
                if(pairInfo.commonCount > 0) {
                    const double jaccard = pairInfo.jaccard();
                    const uint64_t hue = uint64_t(std::round(jaccard * 120.));
                    color = "hsl(" + to_string(hue) + ",100%, 50%)";
                } else {
                    color = "blue";
                }
            } else if(options.segmentColoring == "byUnexplainedFractionOnReferenceSegment") {
                const auto& pairInfo = segmentPairInformationTable[v];
                if(pairInfo.commonCount > 0) {
                    const double fraction = 1. - pairInfo.unexplainedFraction(0);
                    const uint64_t hue = uint64_t(std::round(fraction * 120.));
                    color = "hsl(" + to_string(hue) + ",100%, 50%)";
                } else {
                    color = "blue";
                }
            } else if(options.segmentColoring == "byUnexplainedFractionOnDisplayedSegment") {
                const auto& pairInfo = segmentPairInformationTable[v];
                if(pairInfo.commonCount > 0) {
                    const double fraction = 1. - pairInfo.unexplainedFraction(1);
                    const uint64_t hue = uint64_t(std::round(fraction * 120.));
                    color = "hsl(" + to_string(hue) + ",100%, 50%)";
                } else {
                    color = "blue";
                }
            } else if(options.segmentColoring == "byCluster") {
                const uint64_t clusterId = assemblyGraph.clusterIds[segmentId];
                if(clusterId == std::numeric_limits<uint64_t>::max()) {
                    color = "Gray";
                } else {
                    if(options.clustersToBeColored.empty()) {
                        // We are coloring all cluster. Use a hash function to decide the color.
                        const uint32_t hashValue = MurmurHash2(&clusterId, sizeof(clusterId), uint32_t(options.hashSeed));
                        const uint32_t hue = hashValue % 360;
                        color = "hsl(" + to_string(hue) + ",100%, 50%)";
                    } else {
                        // We are only coloring some segments.
                        auto it = clusterColorMap.find(clusterId);
                        if(it == clusterColorMap.end()) {
                            color = "Black";
                        } else {
                            color = it->second;
                        }
                    }
                }
            } else if(options.segmentColoring == "path") {
                auto it = pathSegments.find(segmentId);
                if(it == pathSegments.end()) {
                    color = "Black";
                } else {
                    const auto positions = it->second;
                    SHASTA_ASSERT(not positions.empty());
                    if(positions.size() == 1) {
                        const auto& p = positions.front();
                        const uint64_t positionInPath = p.first;
                        const bool isReferenceSegment = p.second;
                        const uint32_t hue = uint32_t(
                            std::round(240. * double(positionInPath) / double(path.segments.size())));
                        color = "hsl(" + to_string(hue) + ",100%, " + (isReferenceSegment ? "40%" : "70%") + ")";
                    } else {
                        // This segment appears more than once on the path.
                        color = "Fuchsia";
                    }
                }
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
            "<path id='marker" << segmentId << "' d='M 0 0 L 0.1 0 L 0.6 0.5 L 0.1 1 L 0 1 z' "
            "fill='" << color << "' "
            "/>\n"
            "</marker>\n"
            "</defs>\n";

        // Add this segment to the svg.
        const auto& segmentPairInfo = segmentPairInformationTable[v];
        const auto oldPrecision = svg.precision(1);
        const auto oldFlags = svg.setf(std::ios_base::fixed, std::ios_base::floatfield);

        if(options.segmentColoring == "path") {
            svg << "<g>";
            auto it = pathSegments.find(segmentId);
            if(it != pathSegments.end()) {
                const auto positions = it->second;
                SHASTA_ASSERT(not positions.empty());
                svg << "<title>";
                for(const auto& p: positions) {
                    svg << p.first << " ";
                }
                svg << "</title>";
            }
        }

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
            assemblyGraph.markerGraphPaths.size(segmentId) << "," <<
            averageEdgeCoverage << "," <<
            assemblyGraph.clusterIds[segmentId] << "," <<
            segmentPairInfo.totalCount[0] << "," <<
            segmentPairInfo.totalCount[1] << "," <<
            segmentPairInfo.shortCount[0] << "," <<
            segmentPairInfo.shortCount[1] << "," <<
            segmentPairInfo.commonCount << "," <<
            segmentPairInfo.unexplainedCount[0] << "," <<
            segmentPairInfo.unexplainedCount[1] << ")'" <<
            " onmouseleave='onMouseExitSegment()'" <<

#if 0
            // Old code that displays the segment as a cubic spline.
            // This can create artifacts when the segment is very thick.
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
        if(options.segmentColoring == "path") {
            svg << "</g>";
        }
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

        const uint64_t pathLength = assemblyGraph.markerGraphPaths.size(segmentId);
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
        const LocalAssemblyGraphEdge& edge = localAssemblyGraph[e];
        const uint64_t linkId = edge.linkId;
        const AssemblyGraph::Link& link = assemblyGraph.links[linkId];

        double edgeLength;
        if(link.segmentsAreAdjacent) {
            edgeLength = options.minimumLinkLength;
        } else {
            const int32_t linkSeparation = max(link.separation, 0);
            edgeLength = 3. * options.minimumLinkLength + double(linkSeparation) * options.additionalLinkLengthPerMarker;
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

    const auto path0 = assemblyGraph.markerGraphPaths[segmentId0];
    const auto path1 = assemblyGraph.markerGraphPaths[segmentId1];

    const MarkerGraphEdgeId edgeId0 = path0.back();
    const MarkerGraphEdgeId edgeId1 = path1.front();

    const MarkerGraph::Edge& edge0 = markerGraph.edges[edgeId0];
    const MarkerGraph::Edge& edge1 = markerGraph.edges[edgeId1];

    return edge0.target == edge1.source;
}



// Return the average link separation for the Link
// described by an edge.
int32_t LocalAssemblyGraph::linkSeparation(edge_descriptor e) const
{
    const LocalAssemblyGraph& localAssemblyGraph = *this;
    const uint64_t linkId = localAssemblyGraph[e].linkId;
    return assemblyGraph.links[linkId].separation;
}



// Construct the svg options from an html request.
LocalAssemblyGraph::SvgOptions::SvgOptions(const vector<string>& request)
{
    // The initial layout method if set to "custom" if
    // command "customLayout" is available, "neato" otherwise.
    static bool firstTime = true;
    static string layoutDefaultMethod = "neato";
    if(firstTime) {
        firstTime = false;
        const string command = "which customLayout";
        const int returnCode = system(command.c_str());
        if(returnCode == 0) {
            layoutDefaultMethod = "custom";
        }
    }
    layoutMethod = layoutDefaultMethod;

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
    HttpServer::getParameterValue(request, "hashSeed", hashSeed);
    HttpServer::getParameterValue(request, "pathStart", pathStart);
    HttpServer::getParameterValue(request, "pathDirection", pathDirection);

    string clustersToBeColoredString;
    HttpServer::getParameterValue(request, "clustersToBeColored", clustersToBeColoredString);
    clustersToBeColored.clear();
    if(not clustersToBeColoredString.empty()) {
        vector<string> tokens;
        boost::algorithm::split(tokens, clustersToBeColoredString, boost::algorithm::is_any_of(" "));
        for(const string& token: tokens) {
            try {
                const uint64_t clusterId =std::stoi(token);
                clustersToBeColored.push_back(clusterId);
            } catch(const std::exception&) {
                // Neglect it.
            }
        }
    }

    // Flag to turn on sequence assembly when coloring a path.
    string assemblePathSequenceString;
    assemblePathSequence = HttpServer::getParameterValue(request, "assemblePathSequence", assemblePathSequenceString);

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
        ">Random<hr>"

        // Uniform segment coloring.
        "<input type=radio name=segmentColoring value=uniform"
        << (segmentColoring=="uniform" ? " checked=checked" : "") <<
        ">"
        "<input type=text name=segmentColor size=8 style='text-align:center'"
                " value='" << segmentColor << "'>"
        "<hr>"

        // Segment coloring by Jaccard similarity with the reference segment.
        "<input type=radio name=segmentColoring value=byJaccard"
        << (segmentColoring=="byJaccard" ? " checked=checked" : "") <<
        ">By Jaccard similarity with reference segment, without counting short reads"
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

        // Segment coloring by unexplained fraction on the reference segment.
        "<input type=radio name=segmentColoring value=byUnexplainedFractionOnReferenceSegment"
        << (segmentColoring=="byUnexplainedFractionOnReferenceSegment" ? " checked=checked" : "") <<
        ">By unexplained fraction on the reference segment"
        "<br>"

        // Segment coloring by unexplained fraction on the displayed segment.
        "<input type=radio name=segmentColoring value=byUnexplainedFractionOnDisplayedSegment"
        << (segmentColoring=="byUnexplainedFractionOnDisplayedSegment" ? " checked=checked" : "") <<
        ">By unexplained fraction on the displayed segment"
        "<br>"

        "Reference segment&nbsp;<input type=text name=referenceSegmentId size=8 style='text-align:center'"
                " value='" << referenceSegmentId << "'><hr>"

        // Segment coloring by cluster id.
        "<input type=radio name=segmentColoring value=byCluster"
        << (segmentColoring=="byCluster" ? " checked=checked" : "") <<
        ">By cluster"
        "<br>"
        "Hash seed&nbsp;<input type=text name=hashSeed size=8 style='text-align:center'"
                " value='" << hashSeed << "'><br>"
        "Only color clusters&nbsp;<input type=text name=clustersToBeColored size=8 style='text-align:center'"
                " value='";
     for(const uint64_t clusterId: clustersToBeColored) {
         html << clusterId << " ";
     }
     html <<  "'><hr>"

         // Segment coloring by local cluster
         // (computed by analyzeSubgraph using as input only the segments at
         // distance less than maxDistance).
         "<input type=radio name=segmentColoring value=byLocalCluster"
         << (segmentColoring=="byLocalCluster" ? " checked=checked" : "") <<
         ">By local cluster"
         "<br>";

         // Segment coloring using a path.
         html <<
             "<hr>"
             "<input type=radio name=segmentColoring value=path"
             << (segmentColoring=="path" ? " checked=checked" : "") <<
             ">Color an assembly path"
             "<br>"
             "Start the path at segment &nbsp;<input type=text name=pathStart size=8 style='text-align:center'"
                     " value='" << pathStart << "'>"
             "<br><input type=radio name=pathDirection value=forward" <<
             (pathDirection=="forward" ? " checked=checked" : "") << "> Forward"
             "<br><input type=radio name=pathDirection value=backward" <<
             (pathDirection=="backward" ? " checked=checked" : "") << "> Backward"
             "<br><input type=radio name=pathDirection value=bidirectional" <<
             (pathDirection=="bidirectional" ? " checked=checked" : "") << "> Both directions" <<
             "<br><input type=checkbox name=assemblePathSequence" <<
             (assemblePathSequence ? " checked=checked" : "") <<
             "> Assemble path sequence.";


        html << "</table>"



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



// Return true if there were no changes in the options
// that affect graph layout changed, compared to another
// SvgOptions object.
bool LocalAssemblyGraph::SvgOptions::hasSameLayoutOptions(const SvgOptions& that) const
{
    return
        (layoutMethod == that.layoutMethod) and
        (minimumSegmentLength == that.minimumSegmentLength) and
        (additionalSegmentLengthPerMarker == that.additionalSegmentLengthPerMarker) and
        (minimumLinkLength == that.minimumLinkLength) and
        (additionalLinkLengthPerMarker == that.additionalLinkLengthPerMarker)
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
        const auto path = assemblyGraph.markerGraphPaths[segmentId];
        gfa <<
            "S\t" << segmentId << "\t" <<
            "*\tLN:i:" << path.size() << "\n";
    }


    // Write the links.
    BGL_FORALL_EDGES(e, localAssemblyGraph, LocalAssemblyGraph) {
        const uint64_t linkId = localAssemblyGraph[e].linkId;
        const mode3::AssemblyGraph::Link& link = assemblyGraph.links[linkId];
        gfa << "L\t" <<
            link.segmentId0 << "\t+\t" <<
            link.segmentId1 << "\t+\t0M\n";
    }

}

