// Shasta.
#include "AssembledSegment.hpp"
#include "MarkerGraph.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/icl/interval.hpp>
#include <boost/icl/right_open_interval.hpp>

// Standard library.
#include "iterator.hpp"
#include "fstream.hpp"



void AssembledSegment::clear()
{
    k = 0;
    vertexCount = 0;
    edgeCount = 0;

    vertexIds.clear();
    edgeIds.clear();
    vertexCoverage.clear();
    edgeCoverage.clear();

    vertexSequences.clear();
    vertexRepeatCounts.clear();

    edgeSequences.clear();
    edgeRepeatCounts.clear();
    edgeOverlappingBaseCounts.clear();

    vertexOffsets.clear();

    runLengthSequence.clear();
    repeatCounts.clear();
    rawSequence.clear();

    vertexRunLengthRange.clear();
    vertexRawRange.clear();
    edgeRunLengthRange.clear();
    edgeRawRange.clear();

    vertexCoverageData.clear();
    edgeCoverageData.clear();
    assembledCoverageData.clear();
}



void AssembledSegment::computeVertexOffsets()
{
    vertexOffsets.resize(vertexCount);
    vertexOffsets[0] = 0;

    for(size_t i=0; i<edgeCount; i++) {
        const uint8_t overlap = edgeOverlappingBaseCounts[i];
        if(overlap > 0) {
            SHASTA_ASSERT(edgeSequences[i].empty());
            SHASTA_ASSERT(edgeRepeatCounts[i].empty());
            vertexOffsets[i+1] = uint32_t(vertexOffsets[i] + k - overlap);
        } else {
            vertexOffsets[i+1] = uint32_t(vertexOffsets[i] + k + edgeSequences[i].size());
        }
    }
}



void AssembledSegment::computeVertexAssembledPortion()
{
    // Compute, for each vertex, the portion of vertex sequence that contributes
    // to the assembly. This is the portion that does not overlap a vertex with greater coverage.
    // (Break ties using vertex ids).
    // An edge with overlapping markers does not contribute to the assembly.
    // An edge with at least one intervening base contributes all of its bases
    // to the assembly.

    vertexAssembledPortion.resize(vertexCount);

    for(int i=0; i<int(vertexCount); i++) {

        // Check previous vertices.
        vertexAssembledPortion[i].first = 0;
        for(int j=i-1; j>=0; j--) {
            if(vertexOffsets[j]+k < vertexOffsets[i]) {
                break;
            }
            if(vertexCoverage[j]>vertexCoverage[i] ||
                (vertexCoverage[j]==vertexCoverage[i] && vertexIds[j]<vertexIds[i])) {
                vertexAssembledPortion[i].first =
                    vertexOffsets[j] + uint32_t(k) - vertexOffsets[i];
                break;
            }
        }

        // Check following vertices.
        vertexAssembledPortion[i].second = uint32_t(k);
        for(int j=i+1; j<int(vertexCount); j++) {
            if(vertexOffsets[i]+k < vertexOffsets[j]) {
                break;
            }
            if(vertexCoverage[j]>vertexCoverage[i] ||
                (vertexCoverage[j]==vertexCoverage[i] && vertexIds[j]<vertexIds[i])) {
                vertexAssembledPortion[i].second = vertexOffsets[j] - vertexOffsets[i];
                break;
            }
        }

        // Handle the case of a vertex that contributes nothing.
        if(vertexAssembledPortion[i].second <= vertexAssembledPortion[i].first) {
            vertexAssembledPortion[i].first = 0;
            vertexAssembledPortion[i].second = 0;
        }
        SHASTA_ASSERT(vertexAssembledPortion[i].second <= k);
    }
}



void AssembledSegment::assemble()
{
    // Figure out if we should store coverage data for assembled sequence.
    const bool storeCoverageData =
        vertexCoverageData.size() == vertexCount &&
        edgeCoverageData.size()   == edgeCount;

    vertexRunLengthRange.resize(vertexCount);
    vertexRawRange.resize(vertexCount);
    edgeRunLengthRange.resize(edgeCount);
    edgeRawRange.resize(edgeCount);

    for(size_t i=0; ; i++) {

        // Vertex.
        vertexRunLengthRange[i].first = uint32_t(runLengthSequence.size());
        vertexRawRange[i].first = uint32_t(rawSequence.size());
        for(uint32_t j=vertexAssembledPortion[i].first; j!=vertexAssembledPortion[i].second; j++) {
            const Base base = vertexSequences[i][j];
            const uint32_t repeatCount = vertexRepeatCounts[i][j];
            SHASTA_ASSERT(repeatCount > 0);
            runLengthSequence.push_back(base);
            repeatCounts.push_back(repeatCount);
            if(storeCoverageData) {
                assembledCoverageData.push_back(vertexCoverageData[i][j]);
            }
            for(uint32_t k=0; k!=repeatCount; k++) {
                rawSequence.push_back(base);
            }
        }
        vertexRunLengthRange[i].second = uint32_t(runLengthSequence.size());
        vertexRawRange[i].second = uint32_t(rawSequence.size());

        // This was the last vertex.
        if(i == edgeCount) {
            break;
        }

        // Edge.
        edgeRunLengthRange[i].first = uint32_t(runLengthSequence.size());
        edgeRawRange[i].first = uint32_t(rawSequence.size());
        if(edgeSequences[i].size() > 0) {
            for(uint32_t j=0; j!=uint32_t(edgeSequences[i].size()); j++) {
                const Base base = edgeSequences[i][j];
                const uint32_t repeatCount = edgeRepeatCounts[i][j];
                SHASTA_ASSERT(repeatCount > 0);
                runLengthSequence.push_back(base);
                repeatCounts.push_back(repeatCount);
                if(storeCoverageData) {
                    assembledCoverageData.push_back(edgeCoverageData[i][j]);
                }
                for(uint32_t k=0; k!=repeatCount; k++) {
                    rawSequence.push_back(base);
                }
            }
        }
        edgeRunLengthRange[i].second = uint32_t(runLengthSequence.size());
        edgeRawRange[i].second = uint32_t(rawSequence.size());
    }
}



// Write out details in html.
void AssembledSegment::writeHtml(
    ostream& html,
    bool showSequence,
    bool showDetails,
    uint32_t begin,
    uint32_t end) const
{
    // Write a title.
    html << "<h1>Assembled segment</h1>";

    // Validate begin, end.
    if(showSequence || showDetails) {
        if(begin>rawSequence.size() || end>rawSequence.size()+1 || end<=begin) {
            html <<  "<p>Invalid begin/end positions in raw sequence. "
                "Raw sequence length is " << rawSequence.size() << " bases.";
            return;
        }
    }

    if(showSequence) {
        writeRawSequenceHtml(html, begin, end);
    }
    if(showDetails) {
        writeDetailHtml(html, begin, end);
    }
}



void AssembledSegment::writeRawSequenceHtml(
    ostream& html,
    uint32_t begin,
    uint32_t end) const
{

    if(begin==0 and end==rawSequence.size()) {
        html << "<p>" << rawSequence.size() <<
            " bases of assembled raw sequence:<br><pre style='font-family:courier'>\n";
    } else {
        html << "<p>Bases " << begin << " to " << end << " (" << end-begin <<
            " bases) of " << rawSequence.size() <<
            " bases of assembled raw sequence :<br><pre style='font-family:courier'>\n";
    }

    // Write the position labels.
    for(size_t i=begin; i<end; ) {
        if((i % 10) == 0) {
            const string label = to_string(i);
            html << label;
            for(size_t j=0; j<10-label.size(); j++) {
                html << " ";
            }
            i += 10;
        } else {
            html << " ";
            i++;
        }
    }
    html <<"\n";

    // Write the position scale.
    for(size_t i=begin; i!=end; i++) {
        if((i % 10) == 0) {
            html << "|";
        } else if((i % 5) == 0) {
            html << "+";
        } else {
            html << ".";
        }
    }
    html <<"\n";

    // Write the sequence.
    for(size_t i=begin; i!=end; i++) {
        html << rawSequence[i];
    }
    html << "</pre>";



    // Button to Blat this portion of this assembled segment.
    // We cannot use a simple <a> because we need to do a POST
    // (the GET request fails when the read is too long).
    html <<
        "<p><form action='https://genome.ucsc.edu/cgi-bin/hgBlat' method=post target='_blank'>"
        "<input type=submit value='Blat "
        "this sequence in the UCSC browser'>"
        "<input type=text hidden name=type value=DNA>"
        // Don't specify the genome.
        // UCSC browser will Blat against last used genome (stored in cookies).
        // "<input type=text hidden name=type value=DNA>"
        // "<input type=text hidden name=name value=Human>"
        // "<input type=text hidden name=db value=hg38>"
        "<input type=text hidden name=userSeq value=";
        for(size_t i=begin; i!=end; i++) {
            html << rawSequence[i];
        }
    html << "></form>";
}



// Write a table with a row for each marker graph vertex or edge
// in the marker graph chain.
void AssembledSegment::writeDetailHtml(
    ostream& html,
    uint32_t begin,
    uint32_t end) const
{

    html <<
        "<p>This vertex of the assembly graph corresponds to a chain of " <<
        edgeIds.size() << " edges in the marker graph. "
        "The table below shows consensus sequences "
        "for the vertices and edges of this chain of the marker graph. "
        "All vertex and edge ids in the table refer to the marker graph."
        "<p><table><tr>"
        "<th rowspan=2>Vertex<br>or<br>edge"
        "<th rowspan=2>Vertex<br>index<br>in<br>chain"
        "<th rowspan=2>Edge<br>index<br>in<br>chain"
        "<th rowspan=2>Global<br>id"
        "<th rowspan=2>Coverage"
        "<th colspan=4>Run-length"
        "<th colspan=3>Raw"
        "<tr>"
        "<th>Offset"
        "<th>Begin"
        "<th>End"
        "<th>Sequence"
        "<th>Begin"
        "<th>End"
        "<th>Sequence";
    for(size_t i=0; ; i++) {

        // Vertex.
        using Roi = boost::icl::right_open_interval<uint32_t>;
        if(boost::icl::intersects(
            Roi(begin, end),
            Roi(vertexRawRange[i].first, vertexRawRange[i].second))) {
            const MarkerGraph::VertexId vertexId = vertexIds[i];
            const vector<Base>& vertexSequence = vertexSequences[i];
            const vector<uint32_t>& vertexRepeatCount = vertexRepeatCounts[i];
            // const uint32_t maxVertexRepeatCount =
            //     *std::max_element(vertexRepeatCount.begin(), vertexRepeatCount.end());
            html <<
                "<tr><td>Vertex" <<
                "<td class=centered title='Vertex index in chain'>" << i << "<td>"
                "<td class=centered><a href='exploreMarkerGraphVertex?vertexId=" <<
                vertexId << "'>" << vertexId << "</a>"
                "<td class=centered title='Vertex coverage'>" << vertexCoverage[i] <<
                "<td class=centered title='Vertex offset in assembled RLE sequence'>" <<
                vertexOffsets[i] <<
                "<td class=centered "
                "title='Begin offset of RLE sequence contributed by this vertex'>" <<
                vertexRunLengthRange[i].first <<
                "<td class=centered "
                "title='End offset (one past) of RLE sequence contributed by this vertex'>" <<
                vertexRunLengthRange[i].second <<
                "<td style='font-family:courier' title="
                "'Vertex consensus RLE sequence. Portion contributed to assembly is highlighted.'>";

            // Vertex RLE sequence.
            for(size_t j=0; j<vertexSequence.size(); j++) {
                if(j==vertexAssembledPortion[i].first &&
                    vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "<span style='background-color:LightGreen'>";
                }
                html << vertexSequence[j];
                if(j==vertexAssembledPortion[i].second-1  &&
                    vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "</span>";
                }
            }
            html << "<br>";

            // Vertex repeat counts.
            for(size_t j=0; j<vertexSequence.size(); j++) {
                const uint32_t repeatCount = vertexRepeatCount[j];
                if(j==vertexAssembledPortion[i].first &&
                    vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "<span style='background-color:LightGreen'>";
                }
                if(repeatCount < 10) {
                    html << repeatCount % 10;
                } else {
                    html << "*";
                }
                if(j==vertexAssembledPortion[i].second-1 &&
                    vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "</span>";
                }
            }
            html <<
                "<td class=centered "
                "title='Begin offset of raw sequence contributed by this vertex'>" <<
                vertexRawRange[i].first <<
                "<td class=centered "
                "title='End offset (one past) of raw sequence contributed by this vertex'>" <<
                vertexRawRange[i].second <<
                "<td style='font-family:courier' "
                "title='Vertex consensus raw sequence. Portion contributed to assembly is highlighted.'>";

            // Vertex raw sequence.
            for(size_t j=0; j<vertexSequence.size(); j++) {
                if(j==vertexAssembledPortion[i].first &&
                    vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "<span style='background-color:LightGreen'>";
                }
                const Base b = vertexSequence[j];
                const uint32_t repeatCount = vertexRepeatCount[j];
                for(uint32_t k=0; k<repeatCount; k++) {
                    html << b;
                }
                if(j==vertexAssembledPortion[i].second-1 &&
                    vertexAssembledPortion[i].first!=vertexAssembledPortion[i].second) {
                    html << "</span>";
                }
            }
        }

        // This was the last vertex.
        if(i == edgeCount) {
            break;
        }



        // Edge.
        if(boost::icl::intersects(
            Roi(begin, end),
            Roi(edgeRawRange[i].first, edgeRawRange[i].second))) {
            const MarkerGraph::EdgeId edgeId = edgeIds[i];
            // const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
            const vector<Base>& edgeSequence = edgeSequences[i];
            const vector<uint32_t>& edgeRepeatCount = edgeRepeatCounts[i];
            const size_t edgeSequenceLength = edgeSequence.size();
            SHASTA_ASSERT(edgeRepeatCount.size() == edgeSequenceLength);
            // const uint32_t maxEdgeRepeatCount =
            //    *std::max_element(edgeRepeatCount.begin(), edgeRepeatCount.end());
            // SHASTA_ASSERT(maxEdgeRepeatCount < 10);  // For now. Add additional code when this fails.
            html <<
                "<tr><td>Edge"
                "<td><td class=centered title='Edge index in chain'>" << i <<
                "<td class=centered><a href='exploreMarkerGraphEdge?edgeId=" << edgeId << "'>" << edgeId << "</a>"
                "<td class=centered title='Edge coverage'>" << edgeCoverage[i] <<
                "<td class=centered>" <<
                "<td class=centered "
                "title='Begin offset of RLE sequence contributed by this edge'>";
            if(edgeRunLengthRange[i].first != edgeRunLengthRange[i].second) {
                html << edgeRunLengthRange[i].first;
            }
            html << "<td class=centered "
                "title='End offset (one past) of RLE sequence contributed by this edge'>";
            if(edgeRunLengthRange[i].first != edgeRunLengthRange[i].second) {
                html << edgeRunLengthRange[i].second;
            }
            html << "<td style='font-family:courier' "
                "title='Edge consensus RLE sequence. Portion contributed to assembly is highlighted.'>";

            // Edge RLE sequence.
            if(edgeSequenceLength > 0) {
                html << "<span style='background-color:pink'>";
                for(size_t j=0; j<edgeSequenceLength; j++) {
                    html << edgeSequence[j];
                }
                html << "</span>";
            }
            html << "<br>";

            // Edge repeat counts.
            if(edgeSequenceLength > 0) {
                html << "<span style='background-color:pink'>";
                for(size_t j=0; j<edgeSequenceLength; j++) {
                    const uint32_t repeatCount = edgeRepeatCount[j];
                    if(repeatCount < 10) {
                        html << repeatCount % 10;
                    } else {
                        html << "*";
                    }
                }
                html << "</span>";
            }
            html << "<td class=centered "
                "title='Begin offset of raw sequence contributed by this edge'>";
            if(edgeRawRange[i].first != edgeRawRange[i].second) {
                html << edgeRawRange[i].first;
            }
            html << "<td class=centered "
                "title='End offset (one past) of raw sequence contributed by this edge'>";
            if(edgeRawRange[i].first != edgeRawRange[i].second) {
                html << edgeRawRange[i].second;
            }
            html << "<td style='font-family:courier' "
                "title='Edge consensus raw sequence. Portion contributed to assembly is highlighted.'>";

            // Edge raw sequence.
            if(edgeSequenceLength > 0) {
                html << "<span style='background-color:pink'>";
                for(size_t j=0; j<edgeSequenceLength; j++) {
                    const Base b = edgeSequence[j];
                    const uint32_t repeatCount = edgeRepeatCount[j];
                    for(uint32_t k=0; k<repeatCount; k++) {
                        html << b;
                    }
                }
                html << "</span>";
            }
         }
    }
}



void AssembledSegment::writeCoverageDataCsv() const
{
    ofstream csv("Coverage.csv");

    for(uint32_t position=0; position<size(); position++) {
        csv << position << ",";
        csv << getBase(position) << ",";
        csv << getRepeatCount(position) << ",";
        for(const auto& coverageData: getCoverageData(position)) {
            csv <<
                coverageData.getBase() <<
                coverageData.getRepeatCount() <<
                coverageData.getStrand() << " " <<
                coverageData.getFrequency() << ",";
        }
        csv << "\n";
    }
}
