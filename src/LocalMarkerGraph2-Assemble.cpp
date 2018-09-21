// Shasta
#include "LocalMarkerGraph2.hpp"
#include "Histogram.hpp"
#include "Marker.hpp"
#include "MemoryMappedVectorOfVectors.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

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
