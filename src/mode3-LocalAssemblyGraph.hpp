#ifndef SHASTA_MODE3_LOCAL_ASSEMBLY_GRAPH_HPP
#define SHASTA_MODE3_LOCAL_ASSEMBLY_GRAPH_HPP

#include "mode3.hpp"

namespace shasta {
    namespace mode3 {

        class LocalAssemblyGraph;
        class LocalAssemblyGraphEdge;
        class LocalAssemblyGraphVertex;
    }

}


// Classes used to display in the http server a local portion of the AssemblyGraph.
class shasta::mode3::LocalAssemblyGraphVertex {
public:
    uint64_t segmentId;
    uint64_t distance;  // From the start vertex.
    array<double, 2> position;
    LocalAssemblyGraphVertex(
        uint64_t segmentId,
        uint64_t distance);
    LocalAssemblyGraphVertex();
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
        string segmentColor = "Green";
        string segmentAtMaxDistanceColor = "LightGray";
        string linkColor = "Black";

        // Construct the options from an html request.
        SvgOptions(const vector<string>& request);

        // Add rows to the html request form.
        void addFormRows(ostream& html);
    };
    void writeSvg(const string& fileName, const SvgOptions&) const;
    void writeSvg(ostream&, const SvgOptions&) const;



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

