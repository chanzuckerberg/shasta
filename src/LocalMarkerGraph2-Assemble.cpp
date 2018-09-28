// Shasta
#include "LocalMarkerGraph2.hpp"
#include "ConsensusCaller.hpp"
#include "Histogram.hpp"
#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/icl/interval_map.hpp>

// Standard libraries.
#include "iterator.hpp"



// Compute the clipped optimal spanning tree path
// and use it to assemble its dominant sequence.
void LocalMarkerGraph2::assembleDominantSequence(
    int maxDistance,
    vector< pair<Base, int> >& sequence)
{
    CZI_ASSERT(!optimalSpanningTreeBestPath.empty());
    computeLocalAssemblyPath(maxDistance);
    assembleDominantSequence(localAssemblyPath, sequence);

}


// Assemble the dominant sequence for a given path.
void LocalMarkerGraph2::assembleDominantSequence(
    const vector<edge_descriptor>& path,
    vector< pair<Base, int> >& sequence)
{
    CZI_ASSERT(!path.empty());
    const LocalMarkerGraph2& graph = *this;

    // This is only used with the run-length representation of reads.
    CZI_ASSERT(useRunLengthReads);

    // Compute SeqAn alignments for all edges on the local assembly path.
    computeSeqanAlignments();

    // Start with empty sequence.
    sequence.clear();

    // Add the sequence of the first vertex, with its coverage
    const vertex_descriptor v0 = source(path.front(), graph);
    const KmerId kmerId0 = getKmerId(v0);
    const Kmer kmer0(kmerId0, k);
    const LocalMarkerGraph2Vertex& vertex0 = graph[v0];
    const auto coverage0 = vertex0.markerInfos.size();
    for(size_t i=0; i<k; i++) {
        sequence.push_back(make_pair(kmer0[i], coverage0));
    }

    // Loop over all edges in the path.
    for(const edge_descriptor e: path) {

        // Add the edge sequence, if any.
        const LocalMarkerGraph2Edge& edge = graph[e];
        CZI_ASSERT(!edge.infos.empty());
        const auto& p = edge.infos.front();
        const auto coverage = p.second.size();
        const auto& edgeSequence = p.first.sequence;
        for(const Base base: edgeSequence) {
            sequence.push_back(make_pair(base, coverage));
        }

        // Add the sequence of the target vertex of this edge.
        const vertex_descriptor v1 = target(e, graph);
        const KmerId kmerId1 = getKmerId(v1);
        const Kmer kmer1(kmerId1, k);
        const LocalMarkerGraph2Vertex& vertex1 = graph[v1];
        const auto coverage1 = vertex1.markerInfos.size();

        // First, handle any overlapping bases.
        for(size_t i=0; i<size_t(p.first.overlappingBaseCount); i++) {
            const size_t positionInAssembledSequence = sequence.size() - p.first.overlappingBaseCount + i;
            CZI_ASSERT(sequence[positionInAssembledSequence].first == kmer1[i]);
            sequence[positionInAssembledSequence].second =
                max(sequence[positionInAssembledSequence].second, int(coverage1));
        }

        // Handle the non-overlapping bases.
        for(size_t i=size_t(p.first.overlappingBaseCount); i<k; i++) {
            sequence.push_back(make_pair(kmer1[i], coverage1));
        }
#if 0
        cout << "After processing edge ";
        cout << graph[source(e, graph)].vertexId << " ";
        cout << graph[target(e, graph)].vertexId << " ";
        cout << sequence.size() << endl;
        for(const auto& p: sequence) {
            const auto coverage = p.second;
            if(coverage < 10) {
                cout << char(std::tolower(p.first.character()));
            } else {
                cout << p.first;
            }
        }
        cout << endl;
        for(const auto& p: sequence) {
            const auto coverage = p.second;
            if(coverage<10) {
                cout << coverage;
            } else {
                cout << " ";
            }
        }
        cout << endl;
#endif
    }

}



// Version of assembleDominantSequence for the case where
// we use a run-length representation of the reads.
// This  assumes that clippedOptimalSpanningTreeBestPath was
// already computed and only creates html output.
void LocalMarkerGraph2::assembleDominantSequence(ostream& html) const
{
    // Shorthands.
    const LocalMarkerGraph2& graph = *this;
    const vector<edge_descriptor>& path = localAssemblyPath;
    using MarkerInfo = LocalMarkerGraph2Vertex::MarkerInfo;
    using InfoWithRepeatCounts = LocalMarkerGraph2Edge::InfoWithRepeatCounts;
    using Sequence = LocalMarkerGraph2Edge::Sequence;

    // If the path is empty, do nothing.
    if(path.empty()) {
        html << "<p>The assembly path is empty.";
        return;
    }



    // Create two vectors with the vertices and edges
    // contributing to assembled sequence and their positions.
    // Edges with overlapping sequence are not included.
    vector< pair<vertex_descriptor, uint32_t> > pathVertices;
    vector< pair<edge_descriptor, uint32_t> > pathEdges;
    uint32_t assembledPosition = 0;
    for(const edge_descriptor e: path) {
        const vertex_descriptor v0 = source(e, graph);
        pathVertices.push_back(make_pair(v0, assembledPosition));
        assembledPosition += k;
        const LocalMarkerGraph2Edge& edge = graph[e];
        CZI_ASSERT(!edge.infos.empty());
        const pair<Sequence, vector<InfoWithRepeatCounts> > p = edge.infos.front();
        const Sequence& sequence = p.first;
        if(sequence.sequence.size() == 0) {
            // This edge does not contribute
            assembledPosition -= sequence.overlappingBaseCount;
        } else {
            pathEdges.push_back(make_pair(e, assembledPosition));
            assembledPosition += uint32_t(sequence.sequence.size());
        }
    }
    const vertex_descriptor vLast = target(path.back(), graph);
    pathVertices.push_back(make_pair(vLast, assembledPosition));
    const uint32_t assembledLength = assembledPosition + k;


    // Write out the contributing vertices and edges.
    if(false) {
        cout << "Vertices contributing to assembled sequence:" << endl;
        for(const auto& p: pathVertices) {
            const vertex_descriptor v = p.first;
            const uint32_t position = p.second;
            const LocalMarkerGraph2Vertex& vertex = graph[v];
            cout << "Vertex " << vertex.vertexId;
            cout << " rank " << vertex.rank;
            cout << " position " << position << endl;
        }
        cout << "Edges contributing to assembled sequence:" << endl;
        for(const auto& p: pathEdges) {
            const edge_descriptor e = p.first;
            const uint32_t position = p.second;
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            const LocalMarkerGraph2Vertex& vertex0 = graph[v0];
            const LocalMarkerGraph2Vertex& vertex1 = graph[v1];
            cout << "Edge " << vertex0.vertexId << "->";
            cout << vertex1.vertexId << "->";
            cout << " ranks " << vertex0.rank << " " << vertex1.rank;
            cout << " position " << position << endl;
        }
        cout << "Assembled length " << assembledLength << endl;
    }


    // Assemble sequence.
    vector<Base> sequence(assembledLength);
    vector<bool> sequenceFlag(assembledLength, false);
    for(const auto& p: pathVertices) {
        const vertex_descriptor v = p.first;
        const uint32_t position = p.second;
        const KmerId kmerId = graph.getKmerId(v);
        const Kmer kmer(kmerId, k);
        for(size_t i=0; i<k; i++) {
            const Base base = kmer[i];
            if(sequenceFlag[position+i]) {
                CZI_ASSERT(sequence[position+i] == base);
            } else {
                sequence[position+i] = base;
                sequenceFlag[position+i] = true;
            }
        }
    }
    for(const auto& p: pathEdges) {
        const edge_descriptor e = p.first;
        const uint32_t position = p.second;
        const LocalMarkerGraph2Edge& edge = graph[e];
        CZI_ASSERT(!edge.infos.empty());
        const pair<Sequence, vector<InfoWithRepeatCounts> > q = edge.infos.front();
        const vector<Base>& edgeSequence = q.first.sequence;
        for(size_t i=0; i<edgeSequence.size(); i++) {
            const Base base = edgeSequence[i];
            CZI_ASSERT(!sequenceFlag[position+i]) ;
            sequence[position+i] = base;
            sequenceFlag[position+i] = true;
        }
    }
    CZI_ASSERT(find(sequenceFlag.begin(), sequenceFlag.end(), false) == sequenceFlag.end());
    if(false) {
        cout << "Assembled sequence:" << endl;
        copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(cout));
        cout << endl;
    }


    // For each position in assembled sequence, store unique
    // combinations OrientedReadId/position/repeatCount.
    // Indexed by position in assembled sequence.
    vector< std::set< tuple<OrientedReadId, uint32_t, uint8_t> > > repeatCountTable(assembledLength);
    for(const auto& p: pathVertices) {
        const vertex_descriptor v = p.first;
        const uint32_t position = p.second;
        const LocalMarkerGraph2Vertex& vertex = graph[v];
        for(const MarkerInfo& markerInfo: vertex.markerInfos) {
            const MarkerId markerId = markerInfo.markerId;
            const CompressedMarker& marker = markers.begin()[markerId];
            const vector<uint8_t> counts = getRepeatCounts(markerInfo);
            for(uint32_t i=0; i<k; i++) {
                repeatCountTable[position+i].insert(
                    make_tuple(markerInfo.orientedReadId, marker.position+i, counts[i]));
            }
        }
    }
    for(const auto& p: pathEdges) {
        const edge_descriptor e = p.first;
        const uint32_t position = p.second;
        const LocalMarkerGraph2Edge& edge = graph[e];
        CZI_ASSERT(!edge.infos.empty());
        const pair<Sequence, vector<InfoWithRepeatCounts> > q = edge.infos.front();
        const Sequence& sequence = q.first;
        const vector<InfoWithRepeatCounts>& infos = q.second;
        CZI_ASSERT(sequence.sequence.size() > 0);
        for(const InfoWithRepeatCounts& info: infos) {
            CZI_ASSERT(info.repeatCounts.size() == sequence.sequence.size());
            const CompressedMarker& marker = markers[info.orientedReadId.getValue()][info.ordinals[0]];
            for(uint32_t i=0; i<info.repeatCounts.size(); i++) {
                repeatCountTable[position+i].insert(
                    make_tuple(info.orientedReadId, marker.position+k+i, info.repeatCounts[i]));

            }

        }
    }
    if(false) {
        for(uint32_t position=0; position<assembledLength; position++) {
            const auto & s = repeatCountTable[position];
            for(const auto& x: s) {
                cout << "Assembled position " << position;
                cout << ": " << get<0>(x) << " " << get<1>(x) << " " << int(get<2>(x)) << endl;
            }
        }
    }



    // For each assembled position, create a histogram of repeat counts.
    vector<Histogram> histogram(assembledLength);
    for(uint32_t position=0; position<assembledLength; position++) {
        std::set< tuple<OrientedReadId, uint32_t, uint8_t> >& s = repeatCountTable[position];
        Histogram& h = histogram[position];
        for(const tuple<OrientedReadId, uint32_t, uint8_t>& t: s) {
            h.increment(get<2>(t));
        }
    }
    if(false) {
        for(uint32_t position=0; position<assembledLength; position++) {
            const Histogram& h = histogram[position];
            for(size_t repeatCount=0; repeatCount<h.size(); repeatCount++) {
                const size_t frequency = h[repeatCount];
                if(frequency) {
                    cout << position << " " << repeatCount << " " << frequency << endl;
                }
            }
        }
    }



    // Begin the big table describing assembled sequence.
    html << "<table style='font-family:monospace;font-size:10px'>";



    // Assembled sequence in run-length representation.
    html << "<tr title='Position on assembled sequence in run-length representation'>"
        "<th class=left style='white-space: nowrap'>Position (run-length)";
    for(uint32_t position=0; position<assembledLength; position++) {
        html << "<td class=centered>" << position;
    }
    html << "<tr title='Assembled sequence in run-length representation'>"
        "<th class=left style='white-space:nowrap'>Assembled sequence (run-length)";
    for(uint32_t position=0; position<assembledLength; position++) {
        const Base base = sequence[position];
        const size_t bestCoverage = histogram[position].bestFrequency();
        char c = base.character();
        if(bestCoverage < 10) {
            c = char(tolower(c));
        }
        html << "<td class=centered>" << c;
    }



    // Best repeat count.
    html << "<tr title='Most frequent repeat count'>"
        "<th class=left style='white-space:nowrap'>Most frequent repeat count";
    for(uint32_t position=0; position<assembledLength; position++) {
        html << "<td  class=centered>" << histogram[position].bestValue();
    }

    // Coverage for best repeat count.
    html << "<tr title='Coverage for most frequent repeat count'>"
        "<th class=left style='white-space:nowrap'>Coverage for most frequent repeat count";
    for(uint32_t position=0; position<assembledLength; position++) {
        html << "<td  class=centered>" << histogram[position].bestFrequency();
    }


    // Total coverage.
    html << "<tr title='Total coverage for all repeat counts'>"
        "<th class=left style='white-space:nowrap'>Total coverage for all repeat counts";
    for(uint32_t position=0; position<assembledLength; position++) {
        const size_t totalCoverage = histogram[position].sum();
        html << "<td  class=centered>" << totalCoverage;
    }


    // Coverage for each repeat count.
    std::set<size_t> activeRepeatCounts;
    for(uint32_t position=0; position<assembledLength; position++) {
        const Histogram& h = histogram[position];
        for(size_t repeatCount=0; repeatCount<h.size(); repeatCount++) {
            if(h[repeatCount]) {
                activeRepeatCounts.insert(repeatCount);
            }
        }
    }
    for(const size_t repeatCount: activeRepeatCounts) {
        html << "<tr title='Coverage for repeat count " << repeatCount << "'>"
            "<th class=left style='white-space:nowrap'>Coverage for repeat count " << repeatCount;
        for(uint32_t position=0; position<assembledLength; position++) {
            const Histogram& h = histogram[position];
            html << "<td  class=centered>";
            if(repeatCount < h.size()) {
                const size_t coverage = h[repeatCount];
                if(coverage) {
                    html << coverage;
                }
            }
        }
    }



    html << "<tr title='Assembled sequence'>"
        "<th class=left style='white-space:nowrap'>Assembled sequence";
    string rawAssembledSequence;
    for(uint32_t position=0; position<assembledLength; position++) {
        const Base base = sequence[position];
        const size_t bestCoverage = histogram[position].bestFrequency();
        const size_t bestRepeatCount = histogram[position].bestValue();
        char c = base.character();
        if(bestCoverage < 10) {
            c = char(tolower(c));
        }
        html << "<td class=centered>";
        for(size_t i=0; i<bestRepeatCount; i++) {
            html << c;
            rawAssembledSequence.push_back(c);
        }
    }



    html << "<tr title='Position in assembled sequence'>"
        "<th class=left style='white-space:nowrap'>Position in assembled sequence";
    size_t totalBaseCount = 0;
    for(uint32_t position=0; position<assembledLength; position++) {
        const size_t bestRepeatCount = histogram[position].bestValue();
        html << "<td class=centered>" << totalBaseCount;
        totalBaseCount += bestRepeatCount;
    }

    // End the table.
    html << "</table>";
    html << "Assembled " << assembledLength << " bases in run-length representation "
        "corresponding to a total " << totalBaseCount << " bases:";




    // Write assembled sequence outside the table.
    html << "<pre style='margin:0px'>";
    for(size_t position=0; position<rawAssembledSequence.size(); position+=10) {
        const string label = to_string(position);
        html << label;
        for(size_t i=0; i<10-label.size(); i++) {
            html << " ";
        }
    }
    html << "</pre>";
    html << "<pre style='margin:0px'>";
    for(size_t position=0; position<rawAssembledSequence.size(); position++) {
        if((position%10)==0) {
            html << "|";
        } else if((position%5)==0) {
            html << "+";
        } else {
            html << ".";
        }
    }
    html << "</pre>";
    html << "<pre style='margin:0px'>" << rawAssembledSequence << "</pre>";

    html <<
        "<a id=fastaDownload>Download in FASTA format</a><br>"
        "<script>"
        "var element = document.getElementById('fastaDownload');"
        "element.setAttribute('href', 'data:text/plain;charset=utf-8,' +"
        "encodeURIComponent('>AssembledSequence " << rawAssembledSequence.size() <<
        "\\n" <<  rawAssembledSequence << "\\n'));"
        "element.setAttribute('download', 'AssembledSequence.fa');"
        "</script>";
}



// Version of assembleDominantSequence for the case where
// we use a run-length representation of the reads,
// and we use SeqAn alignments stored in edges to assemble sequence.
// This  assumes that clippedOptimalSpanningTreeBestPath was
// already computed and only creates html output.
void LocalMarkerGraph2::assembleDominantSequenceUsingSeqan(ostream& html) const
{
    // Shorthands.
    const LocalMarkerGraph2& graph = *this;
    const vector<edge_descriptor>& path = localAssemblyPath;
    using InfoWithRepeatCounts = LocalMarkerGraph2Edge::InfoWithRepeatCounts;
    using Sequence = LocalMarkerGraph2Edge::Sequence;

    // If the path is empty, do nothing.
    if(path.empty()) {
        html << "<p>The assembly path is empty.";
        return;
    }



    // Create two vectors with the vertices and edges
    // contributing to assembled sequence and their positions.
    // Edges with overlapping sequence are not included.
    vector< pair<vertex_descriptor, uint32_t> > pathVertices;
    vector< pair<edge_descriptor, uint32_t> > pathEdges;
    uint32_t assembledPosition = 0;
    for(const edge_descriptor e: path) {
        const vertex_descriptor v0 = source(e, graph);
        pathVertices.push_back(make_pair(v0, assembledPosition));
        assembledPosition += k;
        const LocalMarkerGraph2Edge& edge = graph[e];

        if(edge.seqanAlignmentWasComputed) {
            uint32_t baseCount = 0;
            for(const Coverage& consensusInfo: edge.coverages) {
                if(!consensusCaller(consensusInfo).base.isGap()) {
                    ++baseCount;
                }
            }
            if(baseCount) {
                pathEdges.push_back(make_pair(e, assembledPosition));
                assembledPosition += baseCount;
            }

        } else {
            CZI_ASSERT(!edge.infos.empty());
            const pair<Sequence, vector<InfoWithRepeatCounts> > p = edge.infos.front();
            const Sequence& sequence = p.first;
            if(sequence.sequence.size() == 0) {
                // This edge does not contribute
                assembledPosition -= sequence.overlappingBaseCount;
            } else {
                pathEdges.push_back(make_pair(e, assembledPosition));
                assembledPosition += uint32_t(sequence.sequence.size());
            }
        }
    }
    const vertex_descriptor vLast = target(path.back(), graph);
    pathVertices.push_back(make_pair(vLast, assembledPosition));
    const uint32_t assembledLength = assembledPosition + k;



    // Write out the contributing vertices and edges.
    if(false) {
        cout << "Vertices contributing to assembled sequence:" << endl;
        for(const auto& p: pathVertices) {
            const vertex_descriptor v = p.first;
            const uint32_t position = p.second;
            const LocalMarkerGraph2Vertex& vertex = graph[v];
            cout << "Vertex " << vertex.vertexId;
            cout << " position " << position << endl;
        }
        cout << "Edges contributing to assembled sequence:" << endl;
        for(const auto& p: pathEdges) {
            const edge_descriptor e = p.first;
            const uint32_t position = p.second;
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            const LocalMarkerGraph2Vertex& vertex0 = graph[v0];
            const LocalMarkerGraph2Vertex& vertex1 = graph[v1];
            cout << "Edge " << vertex0.vertexId << "->";
            cout << vertex1.vertexId;
            cout << " position " << position << endl;
        }
        cout << "Assembled length " << assembledLength << endl;
    }



    // Some assembled positions have more than one contributing vertex.
    // For each position, pick the vertex with the maximum coverage.
    // Vector verticesByPosition stores pairs (vertex_descriptor, positionInVertex)
    // for each position.
    vector< pair<vertex_descriptor, size_t> > verticesByPosition(
        assembledLength, make_pair(null_vertex(), 0));
    for(const auto& p: pathVertices) {
        const vertex_descriptor v = p.first;
        const size_t startPosition = p.second;
        const size_t coverage = graph[v].markerInfos.size();
        for(size_t i=0; i<k; i++) {
            const size_t position = startPosition + i;
            auto& p = verticesByPosition[position];
            if(p.first==null_vertex() || coverage > graph[p.first].markerInfos.size()) {
                p.first = v;
                p.second = i;
            }
        }
    }

    // Also construct a similar table of edges by position.
    // Here, the second element in the pair is an index into edge.seqanConsensus
    // (that is, it includes gap positions in the alignment).
    vector< pair<edge_descriptor, size_t> > edgesByPosition(
        assembledLength);
    for(const auto& p: pathEdges) {
        const edge_descriptor e = p.first;
        const uint32_t startPosition = p.second;
        const LocalMarkerGraph2Edge& edge = graph[e];

        if(edge.seqanAlignmentWasComputed) {
            size_t position = startPosition;
            for(size_t offset=0; offset<edge.coverages.size(); offset++) {
                const Coverage& coverage = edge.coverages[offset];
                if(consensusCaller(coverage).base.isGap()) {
                    continue;
                }
                CZI_ASSERT(verticesByPosition[position].first == null_vertex());
                edgesByPosition[position].first = e;
                edgesByPosition[position].second = offset;
                ++position;
            }

        } else {
            CZI_ASSERT(!edge.infos.empty());
            const pair<Sequence, vector<InfoWithRepeatCounts> > p = edge.infos.front();
            const Sequence& sequence = p.first;
            for(size_t offset=0; offset<sequence.sequence.size(); offset++) {
                const size_t position = startPosition + offset; // In this case there are no gaps.
                CZI_ASSERT(verticesByPosition[position].first == null_vertex());
                edgesByPosition[position].first = e;
                edgesByPosition[position].second = offset;
            }
        }
    }

    if(false) {
        cout << "Vertex or edge used for each assembly position:" << endl;
        for(size_t position=0; position<assembledLength; position++) {
            const auto& p = verticesByPosition[position];
            cout << position;
            if(p.first != null_vertex()) {
                cout << " vertex " << graph[p.first].vertexId << " " << p.second;
            } else {
                const auto& q = edgesByPosition[position];
                const edge_descriptor e = q.first;
                const vertex_descriptor v0 = source(e, graph);
                const vertex_descriptor v1 = target(e, graph);
                cout << " edge " << graph[v0].vertexId << "->" << graph[v1].vertexId << " " << q.second;
            }
            cout << endl;
        }
    }



    // Create an icl::interval_map containing the vertex or edge
    // corresponding to each interval of assembled positions.
    // The vertex_descriptor is null_vertex() for intervals
    // covered by an edge.
    boost::icl::interval_map<size_t, pair<vertex_descriptor, edge_descriptor> > intervalMap;
    edge_descriptor unitializedEdgeDescriptor;
    for(size_t position=0; position<assembledLength; position++) {
        const auto& p = verticesByPosition[position];
        if(p.first != null_vertex()) {
            intervalMap.insert(make_pair(
                boost::icl::interval<size_t>::right_open(position, position+1),
                make_pair(p.first, unitializedEdgeDescriptor)));
        } else {
            const auto& q = edgesByPosition[position];
            const edge_descriptor e = q.first;
            intervalMap.insert(make_pair(
                boost::icl::interval<size_t>::right_open(position, position+1),
                make_pair(null_vertex(), e)));
        }
    }
    if(false) {
        cout << "Vertex or edge corresponding to each assembly interval:" << endl;
        for(const auto& p: intervalMap) {
            const auto& interval = p.first;
            cout << interval << " ";
            const auto& q = p.second;
            const vertex_descriptor v = q.first;
            if(v == null_vertex()) {
                const edge_descriptor e = q.second;
                const vertex_descriptor v0 = source(e, graph);
                const vertex_descriptor v1 = target(e, graph);
                cout << graph[v0].vertexId << "->" << graph[v1].vertexId;
            } else {
                cout << graph[v].vertexId;
            }
            cout << endl;
        }
    }



    // Create a ConsensusInfo for each assembled positions.
    vector<Coverage> coverages(assembledLength);
    for(size_t position=0; position<assembledLength; position++) {
        const auto& p = verticesByPosition[position];
        const vertex_descriptor v = p.first;
        if(v != null_vertex()) {

            // Get the Coverage for this position from the vertex.
            coverages[position] = graph[v].coverages[p.second];
        } else {

            // Get the Coverage for this position from the edge.
            const auto& q = edgesByPosition[position];
            const edge_descriptor e = q.first;
            const size_t positionInEdge = q.second;
            const LocalMarkerGraph2Edge& edge = graph[e];

            if(edge.seqanAlignmentWasComputed) {
                coverages[position] = edge.coverages[positionInEdge];
            } else {

                // This is an end case where the SeqAn alignment for the edge
                // was not computed. Use the most frequent sequence instead.

                // Start by gathering the information we need.
                const LocalMarkerGraph2Edge& edge = graph[e];
                CZI_ASSERT(!edge.infos.empty());
                const auto& p = edge.infos.front();
                const LocalMarkerGraph2Edge::Sequence& s = p.first;
                CZI_ASSERT(s.overlappingBaseCount == 0);
                const vector<Base>& sequence = s.sequence;
                const vector<LocalMarkerGraph2Edge::InfoWithRepeatCounts>& infos = p.second;

                // Loop over the infos.
                Coverage& coverage = coverages[position];
                for(const auto& info: infos) {
                    coverage.addRead(
                        AlignedBase(sequence[positionInEdge]),
                        info.orientedReadId.getStrand(),
                        info.repeatCounts[positionInEdge]);
                }
            }
        }
    }


    // Create a vector of consensus bases and repeat counts
    // computed using the consensus caller.
    vector<Consensus> consensus(coverages.size());
    for(size_t position=0; position<consensus.size(); position++) {
        consensus[position] = consensusCaller(coverages[position]);
    }



    // Gather raw assembled sequence.
    vector<Base> rawAssembledSequence;
    for(auto& c: consensus) {
        const AlignedBase base = c.base;
        const size_t repeatCount = c.repeatCount;
        for(size_t k=0; k<repeatCount; k++) {
            rawAssembledSequence.push_back(Base(base));
        }
    }



    // Write assembled sequence.
    html << "<p>Assembed " << rawAssembledSequence.size();
    html << " bases ( " << coverages.size() << " bases in run-length representation):<br>";
    html << "<pre style='margin:0px'>";
    for(size_t position=0; position<rawAssembledSequence.size(); position+=10) {
        const string label = to_string(position);
        html << label;
        for(size_t i=0; i<10-label.size(); i++) {
            html << " ";
        }
    }
    html << "</pre>";
    html << "<pre style='margin:0px'>";
    for(size_t position=0; position<rawAssembledSequence.size(); position++) {
        if((position%10)==0) {
            html << "|";
        } else if((position%5)==0) {
            html << "+";
        } else {
            html << ".";
        }
    }
    html << "</pre>";
    html << "<pre style='margin:0px'>";
    copy(rawAssembledSequence.begin(), rawAssembledSequence.end(), ostream_iterator<Base>(html));
    html << "</pre>";

    html <<
        "<a id=fastaDownload>Download in FASTA format</a><br>"
        "<script>"
        "var element = document.getElementById('fastaDownload');"
        "element.setAttribute('href', 'data:text/plain;charset=utf-8,' +"
        "encodeURIComponent('>AssembledSequence " << rawAssembledSequence.size() <<
        "\\n";
    copy(rawAssembledSequence.begin(), rawAssembledSequence.end(), ostream_iterator<Base>(html));
    html << "\\n'));"
        "element.setAttribute('download', 'AssembledSequence.fa');"
        "</script>";



    // Write out to a big table the information contained in the consensusInfos vector.
    // Make the table hidden for now and display it later, to avoid flickering.
    html << "<p><table id=assemblyTable style='text-align:center;font-size:12px;visibility:hidden'>";



    // Table row containing positions in the run-length assembled sequence,
    // which are also indices in the consensusInfo vector.
    html <<
        "<tr style='background-color:#ffe6e6' title='Position (run-length)'>"
        "<th style='min-width:200px;text-align:left'>Position (run-length)";
    for(size_t position=0; position<coverages.size(); position++) {
        html << "<td>" << position;
    }



    // Table row containing consensus bases in run-length sequence.
    html <<
        "<tr style='background-color:#ffcccc' title='Consensus base (run-length)'>"
        "<th style='min-width:200px;text-align:left'>Consensus base (run-length)";
    for(const auto& c: consensus) {
        html << "<td>" << c.base;
    }



    // Table row containing consensus repeat counts.
    html <<
        "<tr style='background-color:#ffffe6' title='Repeat count consensus'>"
        "<th style='min-width:200px;text-align:left'>Repeat count consensus";
    for(const auto& c: consensus) {
        html << "<td>" << c.repeatCount;
    }



    // Table row containing position in raw sequence.
    size_t rawPosition = 0;
    html <<
        "<tr style='background-color:#e6ffe6' title='Position (raw)'>"
        "<th style='min-width:200px;text-align:left'>Position (raw)";
    for(const auto& c: consensus) {
        html << "<td>" << rawPosition;
        rawPosition += c.repeatCount;
    }



    // Table row containing consensus bases in raw sequence.
    html <<
        "<tr style='background-color:#ccffcc' title='Consensus bases (raw)'>"
        "<th style='min-width:200px;text-align:left'>Consensus bases (raw)";
    for(const auto& c: consensus) {
        const AlignedBase base = c.base;
        const size_t repeatCount = c.repeatCount;
        html << "<td>";
        for(size_t i=0; i<repeatCount; i++) {
            html << base;
        }
    }



    // Four rows with coverage for ACGT.
    for(size_t b=0; b<4; b++) {
        const Base base = Base::fromInteger(b);
        const string rowTitle = string("Coverage for ") + base.character();
        html <<
            "<tr style='background-color:#e6f2ff' title='" << rowTitle << "'>"
            "<th style='min-width:200px;text-align:left'>" << rowTitle;
        for(const Coverage& coverage: coverages) {
            const size_t baseCoverage = coverage.coverage(AlignedBase(base));
            html << "<td>";
            if(baseCoverage) {
                html << baseCoverage;
            }
        }
    }



    // Table row containing coverage for the consensus base.
    html <<
        "<tr style='background-color:#cce6ff' title='Coverage for consensus base'>"
        "<th style='min-width:200px;text-align:left'>"
        "Coverage for consensus base";
    for(size_t position=0; position<coverages.size(); position++) {
        const Coverage& coverage = coverages[position];
        const Consensus& c = consensus[position];
        const AlignedBase base = c.base;
        html << "<td>" << coverage.coverage(base);
    }



    // One row with coverage for each represented repeat count.
    const std::set<size_t> repeatCounts = consensusCaller.findRepeatCounts(coverages);
    for(const size_t repeatCount: repeatCounts) {
        const string rowTitle = "Coverage for repeat count " + to_string(repeatCount);
        html <<
            "<tr style='background-color:#fff2e6' title='" << rowTitle << "'>"
            "<th style='min-width:200px;text-align:left'>" << rowTitle;
        for(size_t position=0; position<consensus.size(); position++) {
            const Coverage& coverage = coverages[position];
            const Consensus& c = consensus[position];
            const AlignedBase base = c.base;
            const size_t repeatCountCoverage = coverage.coverage(base, repeatCount);
            html << "<td>";
            if(repeatCountCoverage) {
                html << repeatCountCoverage;
            }
        }
    }



    // One row with coverage for the consensus repeat count.
    string rowTitle = "Coverage for repeat count consensus";
    html <<
        "<tr style='background-color:#ffe6cc' title='" << rowTitle << "'>"
        "<th style='min-width:200px;text-align:left'>" << rowTitle;
    for(size_t position=0; position<consensus.size(); position++) {
        const Coverage& coverage = coverages[position];
        const Consensus& c = consensus[position];
        const AlignedBase base = c.base;
        const size_t repeatCount = c.repeatCount;
        const size_t repeatCountCoverage = coverage.coverage(base, repeatCount);
        html << "<td>";
        if(repeatCountCoverage) {
            html << repeatCountCoverage;
        }
    }



    // Table row with links to vertices or edges.
    rowTitle = "Vertex or edge";
    html <<
        "<tr style='background-color:#e6e6e6' title='" << rowTitle << "'>"
        "<th style='min-width:200px;text-align:left'>" << rowTitle;
    for(const auto& p: intervalMap) {
        const auto& interval = p.first;
        html << "<td colspan=" << boost::icl::length(interval) << ">";
        const auto& q = p.second;
        const vertex_descriptor v = q.first;
        if(v == null_vertex()) {
            const edge_descriptor e = q.second;
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            html <<
                "<span style='cursor:pointer' onClick='positionAtVertex(" << graph[v0].vertexId << ")'>" <<
                graph[v0].vertexId << "</span>" <<
                "&rarr;"
                "<span style='cursor:pointer' onClick='positionAtVertex(" << graph[v1].vertexId << ")'>" <<
                graph[v1].vertexId << "</span>";
        } else {
            html <<
                "<span style='cursor:pointer' onClick='positionAtVertex(" << graph[v].vertexId << ")'>" <<
                graph[v].vertexId << "</span>";
        }
    }



    // Finish the table. Make it visible now. This avoids flickering.
    html << "</table>"
        "<script>"
        "document.getElementById('assemblyTable').style.visibility = 'visible';"
        "</script>";
}
