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



    // Svg output.
    // There are two flavors:
    // - The writeSvg1 functions use sfdp, which does not honor the specified edge length.
    //   So, we add additional vertices to achieve approximate display lengths of
    //   segment and links.
    // - The writeSvg2 functions use neato, which honors the specified edge lengths.
    //   So we only need to generate 3 vertices for each segments.
    // The writeSvg1 functions generally achieve better layouts.
    class SvgOptions {
    public:

        uint64_t sizePixels = 800;
        double segmentLengthScalingFactor = 2.;
        double segmentThickness = 6.;
        string segmentColor = "Green";
        string segmentAtZeroDistanceColor = "LightGreen";
        string segmentAtMaxDistanceColor = "Cyan";
        double nonConsecutiveLinkLengthScalingFactor = 2.;
        double minimumLinkThickness = 1.;
        double linkThicknessScalingFactor = 0.1;
        string linkColor = "black";

        // Construct the options from an html request.
        SvgOptions(const vector<string>& request);

        // Add rows to the html request form.
        void addFormRows(ostream& html);
    };
    void writeSvg1(const string& fileName, const SvgOptions&) const;
    void writeSvg1(ostream&, const SvgOptions&) const;
    void writeSvg2(const string& fileName, const SvgOptions&) const;
    void writeSvg2(ostream&, const SvgOptions&) const;

    // Layout computation for writeSvg2.
    using AuxiliaryGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
    static void computeLayout(
        const AuxiliaryGraph&,
        const std::map<AuxiliaryGraph::edge_descriptor, double>& edgeLength,
        std::map<AuxiliaryGraph::vertex_descriptor, array<double, 2> >& positionMap);



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

