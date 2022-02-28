#ifndef SHASTA_MODE3_LOCAL_ASSEMBLY_GRAPH_HPP
#define SHASTA_MODE3_LOCAL_ASSEMBLY_GRAPH_HPP

#include "mode3.hpp"

#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/arithmetic/arithmetic.hpp>


namespace shasta {
    namespace mode3 {

        class LocalAssemblyGraph;
        class LocalAssemblyGraphEdge;
        class LocalAssemblyGraphVertex;

        using Point = boost::geometry::model::d2::point_xy<double>;
    }

}


// Classes used to display in the http server a local portion of the AssemblyGraph.
class shasta::mode3::LocalAssemblyGraphVertex {
public:
    uint64_t segmentId;
    uint64_t distance;  // From the start vertex.
    LocalAssemblyGraphVertex(
        uint64_t segmentId,
        uint64_t distance);
    LocalAssemblyGraphVertex();

    // The positions of the auxiliary graph vertices corresponding
    // to this segment.
    vector<Point> position;

    // Unit vectors for the outward pointing tangents at the two ends of the segment.
    // The are computed as averages of the directions of the
    // incoming/outgoing links.
    // They are used to display the segment as a cubic spline.
    Point t1;
    Point t2;
};



class shasta::mode3::LocalAssemblyGraphEdge {
public:
    uint64_t linkId;
    LocalAssemblyGraphEdge(uint64_t linkId=0) :
        linkId(linkId)
        {}
};



class shasta::mode3::LocalAssemblyGraph :
    public boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS,
    LocalAssemblyGraphVertex, LocalAssemblyGraphEdge> {
public:

    LocalAssemblyGraph(
        const MarkerGraph&,
        const AssemblyGraph&,
        uint64_t startSegmentId,
        uint64_t maxDistance);

    const MarkerGraph& markerGraph;
    const AssemblyGraph& assemblyGraph;
    uint64_t maxDistance;

    vertex_descriptor addVertex(
        uint64_t segmentId,
        uint64_t distance);



    class SvgOptions {
    public:

        double pixelsPerUnitLength = 20.;
        string layoutMethod = "custom";



        // Segment length and thickness.

        // The display length of a segment is computed as
        // minimumSegmentLength + (n-1) * additionalSegmentLengthPerMarker
        // where n is the path length of the segment, in markers.
        double minimumSegmentLength = 1.;
        double additionalSegmentLengthPerMarker = 0.2;

        double segmentThickness = 0.3;



        // Link length and thickness.

        // The display length of a link is computed as follows:
        // - For a link between segments that are consecutive in the marker graph:
        //   linkLength = minimumLinkLength
        // - For a link between segments that are not consecutive in the marker graph:
        //   linkLength = 3 * minimumLinkLength + linkSeparation * additionalLinkLengthPerMarker
        //   (with the linkSeperation replaced with zero if it is negative).
        double minimumLinkLength = 1;
        double additionalLinkLengthPerMarker = 0.2;

        // The display thickness of a link is computed as
        // minimumLinkThickness + (n-1) * additionalSegmentLengthPerMarker
        // where n is the path length of the segment, in markers.
        double minimumLinkThickness = 0.05;
        double additionalLinkThicknessPerRead = 0.005;



        // Colors.
        string segmentAtMaxDistanceColor = "LightGray";
        string linkColor = "Black";

        // Construct the options from an html request.
        SvgOptions(const vector<string>& request);

        // Add rows to the html request form.
        void addFormRows(ostream& html);
    };
    void writeSvg(const string& fileName, const SvgOptions&) const;
    void writeSvg(ostream&, const SvgOptions&) const;
    void computeLayout(const SvgOptions&, double timeout);
    void computeSegmentTangents();
    void computeSegmentTangents(vertex_descriptor);

    // Return the svg color for a segment.
    static string segmentColor(uint64_t segmentId);



    bool haveConsecutivePaths(
        vertex_descriptor v1,
        vertex_descriptor v2) const;

    // Return the average link separation for the Link
    // described by an edge.
    double linkSeparation(edge_descriptor) const;

    // Write the local assembly graph in gfa format.
    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
};
#endif

