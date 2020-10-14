// Shasta.
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "MarkerGraph2.hpp"
#include "orderPairs.hpp"
#include "seqan.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/disjoint_sets.hpp>

// Seqan.
#include <seqan/align.h>

// Standard library.
#include "array.hpp"
#include "fstream.hpp"
#include <map>
#include <set>



// An experimental De Bruijn graph used below to align pseudo-paths.
// Each vertex stores a sequence of k segment ids that occur
// consecutively in one or more pseudo-paths.
// An directed edge is v0->v1 is created between vertices v0 and v1
// if the last k-1 segments of v0 are the same as the first
// k-1 segment of v1.
namespace shasta {
    namespace pseudoPaths {
        template<int k> class DeBruijnGraph;
        template<int k> class DeBruijnGraphVertex;
        template<int k> class DeBruijnGraphEdge;
        template<int k> using DeBruijnGraphBaseClass =
            boost::adjacency_list<
            boost::listS,
            boost::listS,
            boost::bidirectionalS,
            DeBruijnGraphVertex<k>,
            DeBruijnGraphEdge<k> >;
    }
}
using namespace shasta::pseudoPaths;

template<int k> class shasta::pseudoPaths::DeBruijnGraphVertex {
public:

    // The k segment ids associated with this vertex.
    array<uint64_t, k> segmentIds;

    // The oriented reads that have this sequence of segment ids
    // on their pseudo-path, each stored with the starting position
    // in the pseudo-path.
    vector< pair<OrientedReadId, uint64_t > > orientedReadIds;

    void writeGraphviz(ostream& s) const
    {
        s << "\"";
        for(uint64_t i=0; i<k; i++) {
            s << segmentIds[i];
            if(i != k-1) {
                s << "\\n";
            }
        }
        s << "\"";
    }
    void write(ostream& s) const
    {
        for(uint64_t i=0; i<k; i++) {
            s << segmentIds[i];
            if(i != k-1) {
                s << ",";
            }
        }
    }
};



template<int k> class shasta::pseudoPaths::DeBruijnGraphEdge {
public:

    // The k-1 segment ids associated with this edge.
    array<uint64_t, k-1> segmentIds;

    // The oriented reads that have this sequence of segment ids
    // on their pseudo-path, each stored with the starting position
    // in the pseudo-path.
    vector< pair<OrientedReadId, uint64_t > > orientedReadIds;

    DeBruijnGraphEdge(const array<uint64_t, k-1>& segmentIds) :
        segmentIds(segmentIds) {}
};



template<int k> class shasta::pseudoPaths::DeBruijnGraph :
    public DeBruijnGraphBaseClass<k> {
public:
    using Graph = DeBruijnGraph;
    using SegmentId = AssemblyGraph::EdgeId;
    using BaseClass = DeBruijnGraphBaseClass<k>;
    using vertex_descriptor = typename BaseClass::vertex_descriptor;
    using edge_descriptor = typename BaseClass::edge_descriptor;

    // The vertices, keyed by the k segment ids.
    std::map< array<uint64_t, k>, vertex_descriptor> vertexMap;



    // Add vertices for an oriented read.
    void addOrientedReadVertices(
        OrientedReadId orientedReadId,
        const Assembler::PseudoPath& pseudoPath)
    {
        Graph& graph = *this;

        // Loop over possible starting positions such that the k segments
        // starting there are on the pseudo-path.
        for(uint64_t startPosition=0; startPosition+k<=pseudoPath.size(); startPosition++) {

            // Extract the k segments from the pseudo-path
            array<uint64_t, k> segmentIds;
            for(uint64_t i=0; i<k; i++) {
                segmentIds[i] = pseudoPath[startPosition+i].segmentId;
            }

            // Get the vertex corresponding to these k segments, creating it
            // if necessary.
            vertex_descriptor v;
            auto it = vertexMap.find(segmentIds);
            if(it == vertexMap.end()) {
                v = add_vertex(graph);
                graph[v].segmentIds = segmentIds;
                vertexMap.insert(make_pair(segmentIds, v));
            } else {
                v = it->second;
            }
            graph[v].orientedReadIds.push_back(make_pair(orientedReadId, startPosition));
        }
    }


    void createEdges()
    {
        Graph& graph = *this;

        // Index the vertices by their first k-1 segment ids.
        std::map< array<uint64_t, k-1>, vector<vertex_descriptor> > vertexIndex;
        BGL_FORALL_VERTICES_T(v, graph, Graph) {
            const array<uint64_t, k>& segmentIds = graph[v].segmentIds;
            array<uint64_t, k-1> firstSegmentIds;
            const auto begin = segmentIds.begin();
            const auto end = begin + (k-1);
            copy(begin, end, firstSegmentIds.begin());
            vertexIndex[firstSegmentIds].push_back(v);
        }

        // Use the index to create the edges.
        BGL_FORALL_VERTICES_T(v0, graph, Graph) {
            const array<uint64_t, k>& segmentIds0 = graph[v0].segmentIds;
            array<uint64_t, k-1> lastSegmentIds0;
            const auto begin = segmentIds0.begin() + 1;
            const auto end = segmentIds0.end();
            copy(begin, end, lastSegmentIds0.begin());

            for(const vertex_descriptor v1: vertexIndex[lastSegmentIds0]) {
                edge_descriptor e;
                tie(e, ignore) = add_edge(v0, v1, DeBruijnGraphEdge<k>(lastSegmentIds0), graph);
                findOrientedReadIds(e);
            }
        }
    }



    void findOrientedReadIds(edge_descriptor e)
    {
        Graph& graph = *this;
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const auto& orientedReadIds0 = graph[v0].orientedReadIds;
        const auto& orientedReadIds1 = graph[v1].orientedReadIds;

        // This could be done better.
        for(const auto& p0: orientedReadIds0) {
            auto p1 = p0;
            ++p1.second;
            if(binary_search(orientedReadIds1.begin(), orientedReadIds1.end(), p1)) {
                graph[e].orientedReadIds.push_back(p0);
            }
        }
    }


    void removeLowCoverageEdges(uint64_t minCoverage)
    {
        Graph& graph = *this;
        vector<edge_descriptor> edgesToBeRemoved;
        BGL_FORALL_EDGES_T(e, graph, Graph) {
            if(graph[e].orientedReadIds.size() < minCoverage) {
                edgesToBeRemoved.push_back(e);
            }
        }

        for(const edge_descriptor e: edgesToBeRemoved) {
            remove_edge(e, graph);
        }
    }



    void writeGraphviz() const
    {
        const Graph& graph = *this;

        ofstream out("DeBruijnGraph.dot");
        out << "digraph DeBruijnGraph {\n"
            "node[shape=rectangle];\n";
        BGL_FORALL_EDGES_T(e, graph, Graph) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            const auto coverage = graph[e].orientedReadIds.size();
            graph[v0].writeGraphviz(out);
            out << "->";
            graph[v1].writeGraphviz(out);
            out << "[";
            out << "penwidth=\"" << sqrt(double(coverage + 1)) << "\"";
            out << " tooltip=\"(";
            graph[v0].write(out);
            out << ")->(";
            graph[v1].write(out);
            out << ") ";
            out << coverage << "\"";
            if(coverage == 0) {
                out << " style=dotted";
            }
            out << "]";
            out << ";\n";
        }
        out << "}\n";

    }
};



void Assembler::computePseudoPath(
    OrientedReadId orientedReadId,

    // The marker graph path computed using computeOrientedReadMarkerGraphPath.
    // This is computed by this function - it does not need to be filled in
    // in advance.
    vector<MarkerGraph::EdgeId>& path,
    vector< pair<uint32_t, uint32_t> >& pathOrdinals,

    // The pseudo-path computed by this function.
    PseudoPath& pseudoPath) const
{
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    using SegmentId = AssemblyGraph::EdgeId;

    // Compute the marker graph path.
    const uint64_t markerCount = markers.size(orientedReadId.getValue());
    if(markerCount < 2) {
        pathOrdinals.clear();
        path.clear();
    } else {
        computeOrientedReadMarkerGraphPath(
            orientedReadId,
            0, uint32_t(markerCount - 1),
            path, pathOrdinals);
        SHASTA_ASSERT(path.size() == pathOrdinals.size());
    }



    // Now compute the pseudo-path.
    pseudoPath.clear();
    pseudoPath.clear();
    PseudoPathEntry pseudoPathEntry;
    pseudoPathEntry.segmentId = std::numeric_limits<SegmentId>::max();
    for(uint64_t i=0; i<path.size(); i++) {
        const MarkerGraph::EdgeId markerGraphEdgeId = path[i];
        const pair<uint32_t, uint32_t>& ordinals = pathOrdinals[i];

        // Get the corresponding assembly graph segments.
        const span<const pair<SegmentId, uint32_t> > v =
            assemblyGraph.markerToAssemblyTable[markerGraphEdgeId];

        // If no segments, skip.
        if(v.size() == 0) {
            continue;
        }

        // If detangling was used, there can be more than one,
        // and we don't want this here.
        SHASTA_ASSERT(v.size() == 1);

        // There is only one segment.
        const SegmentId segmentId = v.front().first;
        const uint32_t positionInSegment = v.front().second;

        // Same as the previous.
        if(segmentId == pseudoPathEntry.segmentId) {
            pseudoPathEntry.lastOrdinal = ordinals.second;
            pseudoPathEntry.lastPosition = positionInSegment;
            ++pseudoPathEntry.markerGraphEdgeCount;
            continue;
        }

        // This is the next segment edge encountered
        // by this oriented read along its marker graph path.
        if(pseudoPathEntry.segmentId != std::numeric_limits<SegmentId>::max()) {
            pseudoPath.push_back(pseudoPathEntry);
        }
        pseudoPathEntry.segmentId = segmentId;
        pseudoPathEntry.firstOrdinal = ordinals.first;
        pseudoPathEntry.lastOrdinal = ordinals.second;
        pseudoPathEntry.firstPosition = positionInSegment;
        pseudoPathEntry.lastPosition = positionInSegment;
        pseudoPathEntry.markerGraphEdgeCount = 1;
    }

    // Add the last entry.
    if(pseudoPathEntry.segmentId != std::numeric_limits<SegmentId>::max()) {
        pseudoPath.push_back(pseudoPathEntry);
    }
}



// Write the pseudo-path of an oriented read to a csv file.
// The pseudo-path is the sequence of assembly graph edges
// (not necsssarily all adjacent, so not necessatily a path)
// encountered by the oriented read.
void Assembler::writePseudoPath(ReadId readId, Strand strand) const
{
    // Compute the pseudo path.
    const OrientedReadId orientedReadId(readId, strand);
    vector<MarkerGraph::EdgeId> path;
    vector< pair<uint32_t, uint32_t> > pathOrdinals;
    PseudoPath pseudoPath;
    computePseudoPath(orientedReadId, path, pathOrdinals, pseudoPath);

    // Write it out.
    ofstream csv("PseudoPath.csv");
    csv << "Segment id,First ordinal,Last ordinal,"
        "First position in segment,Last position in segment, Marker graph edge count\n";
    for(const PseudoPathEntry& pseudoPathEntry: pseudoPath) {
        csv << pseudoPathEntry.segmentId << ",";
        csv << pseudoPathEntry.firstOrdinal << ",";
        csv << pseudoPathEntry.lastOrdinal << ",";
        csv << pseudoPathEntry.firstPosition << ",";
        csv << pseudoPathEntry.lastPosition << ",";
        csv << pseudoPathEntry.markerGraphEdgeCount << "\n";
    }
}



// Get the vector of segments corresponding to a PseudoPath.
void Assembler::getPseudoPathSegments(
    const PseudoPath& pseudoPath,
    vector<AssemblyGraph::EdgeId>& segmentIds)
{
    segmentIds.clear();
    for(const PseudoPathEntry& pseudoPathEntry: pseudoPath) {
        segmentIds.push_back(pseudoPathEntry.segmentId);
    }
}



void Assembler::alignPseudoPaths(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1)
{
    using SegmentId = AssemblyGraph::EdgeId;
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;

    // Parameters that control the process below. EXPOSE WHEN CODE STABILIZES. *********
    const int matchScore = 1;
    const int mismatchScore = -1;
    const int gapScore = -1;

    // Gather the oriented read ids.
    const array<OrientedReadId, 2> orientedReadIds =
        {OrientedReadId(readId0, strand0), OrientedReadId(readId1, strand1)};
    cout << "Aligning pseudo-paths of " << orientedReadIds[0] <<
        " and " << orientedReadIds[1] << endl;


    // Compute the two pseudo-paths.
    vector<MarkerGraph::EdgeId> path;
    vector< pair<uint32_t, uint32_t> > pathOrdinals;
    PseudoPath pseudoPath;
    array<vector<SegmentId>, 2> pseudoPathSegments;
    for(uint64_t i=0; i<2; i++) {
            computePseudoPath(orientedReadIds[i], path, pathOrdinals,
                pseudoPath);
            getPseudoPathSegments(pseudoPath, pseudoPathSegments[i]);
        cout << "The pseudo-path of " << orientedReadIds[i] <<
            " has " << pseudoPathSegments[i].size() << " segments." << endl;
    }

    // Align them.
    vector< pair<bool, bool> > alignment;
    const uint64_t alignmentScore = shasta::seqanAlign(
        pseudoPathSegments[0].begin(), pseudoPathSegments[0].end(),
        pseudoPathSegments[1].begin(), pseudoPathSegments[1].end(),
        matchScore,
        mismatchScore,
        gapScore,
        true, true,
        alignment);
    cout << "Alignment score " << alignmentScore << endl;
    cout << "Alignment length " << alignment.size() << endl;



    // Write out the alignment.
    uint64_t position0 = 0;
    uint64_t position1 = 0;
    uint64_t weakMatchCount =0;
    uint64_t strongMatchCount =0;
    uint64_t mismatchCount =0;
    uint64_t gapCount =0;
    uint64_t leftUnalignedCount =0;
    uint64_t rightUnalignedCount =0;
    ofstream csv("PseudoPathsAlignment.csv");
    for(const auto& p: alignment) {
        if(p.first) {
            const SegmentId segment0 = pseudoPathSegments[0][position0];
            csv << segment0;
            }
        csv << ",";
        if(p.second) {
            const SegmentId segment1 = pseudoPathSegments[1][position1];
            csv << segment1;
            }
        csv << ",";

        // Write an annotation column.
        if(p.first and p.second) {
            if(pseudoPathSegments[0][position0] != pseudoPathSegments[1][position1]) {
                csv << "Mismatch";
                ++mismatchCount;
            } else {
                // Match.
                // Decide if it is a strong or weak match.
                const SegmentId segmentId = pseudoPathSegments[0][position0];
                const AssemblyGraph::Edge& edge = assemblyGraph.edges[segmentId];
                const AssemblyGraph::VertexId v0 = edge.source;
                const AssemblyGraph::VertexId v1 = edge.target;
                const auto out0 = assemblyGraph.outDegree(v0);
                const auto in1 = assemblyGraph.inDegree(v1);
                if(out0==1 and in1==1) {
                    csv << "Weak match";
                    ++weakMatchCount;
                } else {
                    csv << "Strong match";
                    ++strongMatchCount;
                }
            }
        } else if(position0 == 0 or position1==0) {
            csv << "Left unaligned portion";
            ++leftUnalignedCount;
        } else if(
            position0 == pseudoPathSegments[0].size() or
            position1 == pseudoPathSegments[1].size()) {
            csv << "Right unaligned portion";
            ++rightUnalignedCount;
        } else if(not (p.first and p.second)) {
            csv << "Gap";
            ++gapCount;
        }
        csv << "\n";

        if(p.first) {
            ++position0;
        }
        if(p.second) {
            ++position1;
        }
    }
    SHASTA_ASSERT(position0 == pseudoPathSegments[0].size());
    SHASTA_ASSERT(position1 == pseudoPathSegments[1].size());

    const uint64_t matchCount = weakMatchCount + strongMatchCount;
    cout << "Total match "<< matchCount << endl;
    cout << "Strong match "<< strongMatchCount << endl;
    cout << "Weak match "<< weakMatchCount << endl;
    cout << "Mismatch "<< mismatchCount << endl;
    cout << "Gap "<< gapCount << endl;
    cout << "Left unaligned "<< leftUnalignedCount << endl;
    cout << "Right unaligned "<< rightUnalignedCount << endl;
    cout << "Mismatch/match ratio " << double(mismatchCount)/double(matchCount) << endl;
}
