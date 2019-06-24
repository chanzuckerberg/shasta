#ifndef SHASTA_STATIC_EXECUTABLE

// Shasta.
#include "Assembler.hpp"
#include "ConsensusCaller.hpp"
#include "LocalMarkerGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Spoa.
#include "spoa/spoa.hpp"

// Standard library.
#include "chrono.hpp"



void Assembler::exploreMarkerGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the request parameters.
    LocalMarkerGraphRequestParameters requestParameters;
    getLocalMarkerGraphRequestParameters(request, requestParameters);

    // Write the form.
    requestParameters.writeForm(html, markerGraph.vertices.size());

    // If any required values are missing, stop here.
    if(requestParameters.hasMissingRequiredParameters()) {
        return;
    }



    // Validity checks.
    if(requestParameters.vertexId > markerGraph.vertices.size()) {
        html << "<p>Invalid vertex id " << requestParameters.vertexId;
        html << ". Must be between 0 and " << markerGraph.vertices.size()-1 << " inclusive.";
        return;
    }



    // Create the local marker graph.
    LocalMarkerGraph graph(
        uint32_t(assemblerInfo->k),
        reads,
        readRepeatCounts,
        markers,
        markerGraph.vertexTable,
        *consensusCaller);
    const auto createStartTime = steady_clock::now();
    if(!extractLocalMarkerGraphUsingStoredConnectivity(
        requestParameters.vertexId,
        requestParameters.maxDistance,
        requestParameters.timeout,
        requestParameters.useWeakEdges,
        requestParameters.usePrunedEdges,
        requestParameters.useSuperBubbleEdges,
        graph)) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    if(num_vertices(graph) == 0) {
        html << "<p>The local marker graph is empty.";
        return;
    }
    vector< pair<shasta::Base, int> > sequence;
    const auto createFinishTime = steady_clock::now();
    if(requestParameters.timeout>0 && seconds(createFinishTime - createStartTime) > requestParameters.timeout) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }

    html <<
        "<script>\n"
        "function positionAtVertex(vertexId) {\n";
    if(requestParameters.detailed) {
        html <<
            "var element = document.getElementById('a_vertexDistance' + vertexId);\n";
    } else {
        html <<
            "var element = document.getElementById('vertex' + vertexId);\n";
    }
    html <<
        "var r = element.getBoundingClientRect();\n"
        "window.scrollBy((r.left + r.right - window.innerWidth) / 2, (r.top + r.bottom - window.innerHeight) / 2);\n"
        "}\n"
        "</script>\n";



    // Write it out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = "/dev/shm/" + uuid + ".dot";
    graph.write(
        dotFileName,
        requestParameters.maxDistance,
        requestParameters.detailed);

    // Compute layout in svg format.
    const string command =
        "timeout " + to_string(requestParameters.timeout - int(seconds(createFinishTime - createStartTime))) +
        " dot -O -T svg " + dotFileName +
        " -Gsize=" + to_string(requestParameters.sizePixels/72.);
    const int commandStatus = ::system(command.c_str());
    if(WIFEXITED(commandStatus)) {
        const int exitStatus = WEXITSTATUS(commandStatus);
        if(exitStatus == 124) {
            html << "<p>Timeout for graph layout exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
            filesystem::remove(dotFileName);
            return;
        }
        else if(exitStatus!=0 && exitStatus!=1) {    // sfdp returns 1 all the time just because of the message about missing triangulation.
            filesystem::remove(dotFileName);
            throw runtime_error("Error " + to_string(exitStatus) + " running graph layout command: " + command);
        }
    } else if(WIFSIGNALED(commandStatus)) {
        const int signalNumber = WTERMSIG(commandStatus);
        throw runtime_error("Signal " + to_string(signalNumber) + " while running graph layout command: " + command);
    } else {
        throw runtime_error("Abnormal status " + to_string(commandStatus) + " while running graph layout command: " + command);

    }
    // Remove the .dot file.
    filesystem::remove(dotFileName);



    // Write the graph.

    // Write the legend.
    const string legendName =
        requestParameters.detailed ?
        "MarkerGraphLegend-Detailed.html" :
        "MarkerGraphLegend-Compact.html";
    html <<
        "<h2>Marker graph near marker graph vertex " << requestParameters.vertexId <<
        "</h2>";

    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << svgFile.rdbuf();
    svgFile.close();

    // Remove the .svg file.
    filesystem::remove(svgFileName);

    // Make the vertices clickable to recompute the graph with the
    // same parameters, but starting at the clicked vertex.
    // For a detailed graph, only the "Distance" label of each vertex
    // is made clickable.
    html << "<script>\n";
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph) {
        const LocalMarkerGraphVertex& vertex = graph[v];
        CZI_ASSERT(!vertex.markerInfos.empty());
        const auto& markerInfo = vertex.markerInfos.front();
        const string url =
            "exploreMarkerGraph?vertexId=" + to_string(vertex.vertexId) +
            "&maxDistance=" + to_string(requestParameters.maxDistance) +
            "&sizePixels=" + to_string(requestParameters.sizePixels) +
            "&timeout=" + to_string(requestParameters.timeout) +
            (requestParameters.detailed ? "&detailed=on" : "") +
            (requestParameters.useWeakEdges ? "&useWeakEdges=on" : "") +
            (requestParameters.usePrunedEdges ? "&usePrunedEdges=on" : "") +
            (requestParameters.useSuperBubbleEdges ? "&useSuperBubbleEdges=on" : "");
        if(requestParameters.detailed) {
            html <<
                "document.getElementById('a_vertexDistance" << vertex.vertexId <<
                "').onclick = function() {location.href='" << url << "';};\n";
        } else {
            html <<
                "document.getElementById('vertex" << vertex.vertexId <<
                "').onclick = function() {location.href='" << url << "';};\n";

            // We are displaying the graph in compact mode.
            // Add a right click to recenter and show detailed.
            const string detailUrl =
                "exploreMarkerGraph?readId=" + to_string(markerInfo.orientedReadId.getReadId()) +
                "&strand=" + to_string(markerInfo.orientedReadId.getStrand()) +
                "&ordinal="  + to_string(markerInfo.ordinal) +
                "&maxDistance=1" +
                "&sizePixels=" + to_string(requestParameters.sizePixels) +
                "&timeout=" + to_string(requestParameters.timeout) +
                "&detailed=on";
            html <<
                "document.getElementById('vertex" << vertex.vertexId <<
                "').oncontextmenu = function() {location.href='" << detailUrl << "';"
                "return false;};\n";
        }
    }
    html << "</script>\n";



    // Position the start vertex at the center of the window.
    html <<
        "<script>\n"
        "positionAtVertex(" << requestParameters.vertexId << ");\n"
        "</script>\n";
}



// Extract  from the request the parameters for the display
// of the local marker graph.
void Assembler::getLocalMarkerGraphRequestParameters(
    const vector<string>& request,
    LocalMarkerGraphRequestParameters& parameters) const
{
    parameters.vertexId = 0;
    parameters.vertexIdIsPresent = getParameterValue(
        request, "vertexId", parameters.vertexId);

    parameters.maxDistance = 0;
    parameters.maxDistanceIsPresent = getParameterValue(
        request, "maxDistance", parameters.maxDistance);

    string detailedString;
    parameters.detailed = getParameterValue(
        request, "detailed", detailedString);

    string useWeakEdgesString;
    parameters.useWeakEdges = getParameterValue(
        request, "useWeakEdges", useWeakEdgesString);

    string usePrunedEdgesString;
    parameters.usePrunedEdges = getParameterValue(
        request, "usePrunedEdges", usePrunedEdgesString);

    string useSuperBubbleEdgesString;
    parameters.useSuperBubbleEdges = getParameterValue(
        request, "useSuperBubbleEdges", useSuperBubbleEdgesString);

    parameters.sizePixels = 800;
    parameters.sizePixelsIsPresent = getParameterValue(
        request, "sizePixels", parameters.sizePixels);

    parameters.timeout = 30;
    parameters.timeoutIsPresent = getParameterValue(
        request, "timeout", parameters.timeout);

}



void Assembler::LocalMarkerGraphRequestParameters::writeForm(
    ostream& html,
    MarkerGraph::VertexId vertexCount) const
{
    html <<
        "<h3>Display a local subgraph of the global marker graph</h3>"
        "<form>"

        "<table>"

        "<tr title='Start vertex id between 0 and " << vertexCount << "'>"
        "<td>Start vertex id"
        "<td><input type=text required name=vertexId size=8 style='text-align:center'"
        << (vertexIdIsPresent ? ("value='"+to_string(vertexId)+"'") : "") <<
        ">"

        "<tr title='Maximum distance from start vertex (number of edges)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        << (maxDistanceIsPresent ? ("value='" + to_string(maxDistance)+"'") : " value='6'") <<
        ">"

        "<tr title='Check for detailed graph with labels'>"
        "<td>Detailed"
        "<td class=centered><input type=checkbox name=detailed"
        << (detailed ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to include in the local marker graph "
        "edges that were removed during transitive reduction'>"
        "<td>Edges removed during transitive reduction"
        "<td class=centered>"
        "<input type=checkbox name=useWeakEdges" <<
        (useWeakEdges ? " checked=checked" : "") << ">"

        "<tr title='Check to include in the local marker graph "
        "edges that were removed during pruning'>"
        "<td>Edges removed during pruning"
        "<td class=centered>"
        "<input type=checkbox name=usePrunedEdges" <<
        (usePrunedEdges ? " checked=checked" : "") << ">"

        "<tr title='Check to include in the local marker graph "
        "edges that were removed during bubble/superbubble detection'>"
        "<td>Edges removed during bubble/superbubble detection"
        "<td class=centered>"
        "<input type=checkbox name=useSuperBubbleEdges" <<
        (useSuperBubbleEdges ? " checked=checked" : "") << ">"

        "<tr title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "<td>Graphics size in pixels"
        "<td><input type=text required name=sizePixels size=8 style='text-align:center'"
        << (sizePixelsIsPresent ? (" value='" + to_string(sizePixels)+"'") : " value='800'") <<
        ">"

        "<tr title='Maximum time allowed (seconds) for graph creation and layout, or 0 if unlimited'>"
        "<td>Timeout (seconds) for graph creation and layout"
        "<td><input type=text required name=timeout size=8 style='text-align:center'"
        << (timeoutIsPresent ? (" value='" + to_string(timeout)+"'") : " value='30'") <<
        ">"
        "</table>"



        "<br><input type=submit value='Display'>"
        "</form>";
}



bool Assembler::LocalMarkerGraphRequestParameters::hasMissingRequiredParameters() const
{
    return
        !vertexIdIsPresent ||
        !maxDistanceIsPresent ||
        !timeoutIsPresent;
}



void Assembler::exploreMarkerGraphVertex(const vector<string>& request, ostream& html)
{
    // Get the vertex id.
    MarkerGraph::VertexId vertexId = 0;
    const bool vertexIdIsPresent = getParameterValue(request, "vertexId", vertexId);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Show details for marker graph vertex'> "
        "<input type=text name=vertexId required" <<
        (vertexIdIsPresent ? (" value=" + to_string(vertexId)) : "") <<
        " size=8 title='Enter a vertex id between 0 and " << markerGraph.vertices.size()-1 << "'>";
    html << "</form>";

    // If the vertex id missing or invalid, stop here.
    if(!vertexIdIsPresent || !vertexIdIsPresent) {
        return;
    }
    if(vertexId >= markerGraph.vertices.size()) {
        html << "<p>Invalid vertex id. Must be less than " << markerGraph.vertices.size() << ".";
        return;
    }

    // Access the markers of this vertex.
    MemoryAsContainer<MarkerId> markerIds = markerGraph.vertices[vertexId];
    const size_t markerCount = markerIds.size();
    CZI_ASSERT(markerCount > 0);

    // Get the marker sequence.
    const KmerId kmerId = markers.begin()[markerIds[0]].kmerId;
    const size_t k = assemblerInfo->k;
    const Kmer kmer(kmerId, k);



    // Extract the information we need.
    vector<OrientedReadId> orientedReadIds(markerCount);
    vector<uint32_t> ordinals(markerCount);
    vector< vector<uint8_t> > repeatCounts(markerCount, vector<uint8_t>(k));
    for(size_t j=0; j<markerCount; j++) {
        const MarkerId markerId = markerIds[j];
        const CompressedMarker& marker = markers.begin()[markerId];
        tie(orientedReadIds[j], ordinals[j]) = findMarkerId(markerId);

        // Get the repeat count for this marker at each of the k positions.
        for(size_t i=0; i<k; i++) {
            Base base;
            tie(base, repeatCounts[j][i]) =
                getOrientedReadBaseAndRepeatCount(orientedReadIds[j], uint32_t(marker.position+i));
            CZI_ASSERT(base == kmer[i]);
        }
    }



    // Find all the repeat counts represented.
    std::set<size_t> repeatCountsSet;
    for(const auto& v: repeatCounts) {
        for(const auto r: v) {
            repeatCountsSet.insert(r);
        }
    }



    // Compute consensus repeat counts at each of the k positions.
    vector<size_t> consensusRepeatCounts(k);
    const auto& storedConsensusRepeatCounts =
        markerGraph.vertexRepeatCounts.begin() + k * vertexId;
    for(size_t i=0; i<k; i++) {

        Coverage coverage;
        for(size_t j=0; j<markerCount; j++) {
            coverage.addRead(
                AlignedBase(kmer[i]),
                orientedReadIds[j].getStrand(),
                repeatCounts[j][i]);
        }

        const Consensus consensus = (*consensusCaller)(coverage);
        CZI_ASSERT(Base(consensus.base) == kmer[i]);
        consensusRepeatCounts[i] = consensus.repeatCount;

        // Check that this repeat count agrees with what was
        // computed during the assembly.
        CZI_ASSERT(consensusRepeatCounts[i] == storedConsensusRepeatCounts[i]);
    }




    // Compute concordant and discordant coverage at each position.
    vector<size_t> concordantCoverage(k, 0);
    vector<size_t> discordantCoverage(k, 0);
    for(size_t i=0; i<k; i++) {
        for(size_t j=0; j<markerCount; j++) {
            if(repeatCounts[j][i] == consensusRepeatCounts[i]) {
                ++concordantCoverage[i];
            } else {
                ++discordantCoverage[i];
            }
        }
    }


    // Page title.
    html << "<h1>Marker graph vertex "<< vertexId << "</h1>";
    html << "<p>Vertex coverage (number of markers) is " << markerCount << ".";



    // Write a table with one row for each marker.
    html <<
        "<table>"
        "<tr>"
        "<th>Oriented<br>read"
        "<th title='Ordinal of the marker in the oriented read'>Ordinal"
        "<th>Repeat<br>counts";
    for(size_t j=0; j<markerIds.size(); j++) {
        const OrientedReadId orientedReadId = orientedReadIds[j];
        const ReadId readId = orientedReadId.getReadId();
        const Strand strand = orientedReadId.getStrand();
        const uint32_t ordinal = ordinals[j];

        // Oriented read id.
        html <<
            "<tr>"
            "<td class=centered>"
            "<a href='exploreRead"
            "?readId=" << readId <<
            "&strand=" << strand << "'>" <<
            orientedReadId << "</a>";

        // Marker ordinal.
        html <<
            "<td class=centered>" <<
            "<a href='exploreRead"
            "?readId=" << readId <<
            "&strand=" << strand <<
            "&highlightMarker=" << ordinal <<
            "'>" <<
            ordinal << "</a>";

        // Repeat counts.
        html << "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<k; i++) {
            const uint8_t repeatCount = repeatCounts[j][i];
            if(repeatCount < 10) {
                html << int(repeatCount);
            } else {
                html << "*";
            }
        }
    }



    // Write a row with consensus repeat counts.
    html <<
        "<tr><th colspan=2 class=left>Consensus repeat counts"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<k; i++) {
        const size_t repeatCount = consensusRepeatCounts[i];
        if(repeatCount < 10) {
            html << repeatCount;
        } else {
            html << "*";
        }
    }



    // Write a row with the marker sequence (run-length).
    html <<
        "<tr><th colspan=2 class=left>Run-length sequence"
        "<td class=centered style='font-family:monospace'>";
    kmer.write(html, assemblerInfo->k);



    // Write rows with coverage information for each represented repeat value.
    for(size_t repeatCount: repeatCountsSet) {
        html <<
            "<tr><th colspan=2 class=left>Coverage for repeat count " << repeatCount <<
            "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<k; i++) {

            // Compute coverage for this repeat count, at this position.
            size_t coverage = 0;
            for(size_t j=0; j<markerCount; j++) {
                if(repeatCounts[j][i] == repeatCount) {
                    coverage++;
                }
            }

            // Write it out.
            if(coverage == 0) {
                html << ".";
            } else if(coverage < 10) {
                html << coverage;
            } else {
                html << "*";
            }
        }
    }



    // Write a row with concordant coverage.
    html <<
        "<tr><th colspan=2 class=left>Concordant repeat count coverage"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<assemblerInfo->k; i++) {
        const size_t coverage = concordantCoverage[i];
        if(coverage == 0) {
            html << ".";
        } else if(coverage < 10) {
            html << coverage;
        } else {
            html << "*";
        }
    }



    // Write a row with discordant coverage.
    html <<
        "<tr><th colspan=2 class=left>Discordant repeat count coverage"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<assemblerInfo->k; i++) {
        const size_t coverage = discordantCoverage[i];
        if(coverage == 0) {
            html << ".";
        } else if(coverage < 10) {
            html << coverage;
        } else {
            html << "*";
        }
    }



    // Write a row with the consensus raw sequence.
    html <<
        "<tr><th colspan=2 class=left>Consensus raw sequence"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<k; i++) {
        const Base base = kmer[i];
        const size_t repeatCount = consensusRepeatCounts[i];
        for(size_t r=0; r<repeatCount; r++) {
            html << base;
        }
    }

    html << "</table>";
}



void Assembler::exploreMarkerGraphEdge(const vector<string>& request, ostream& html)
{
    // Get the edge id.
    MarkerGraph::EdgeId edgeId = 0;
    const bool edgeIdIsPresent = getParameterValue(request, "edgeId", edgeId);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Show details for marker graph edge'> "
        "<input type=text name=edgeId required" <<
        (edgeIdIsPresent ? (" value=" + to_string(edgeId)) : "") <<
        " size=8 title='Enter an edge id between 0 and " << markerGraph.edges.size()-1 << "'>";
    html << "</form>";

    // If the edge id missing or invalid, stop here.
    if(!edgeIdIsPresent || !edgeIdIsPresent) {
        return;
    }
    if(edgeId >= markerGraph.edges.size()) {
        html << "<p>Invalid edge id. Must be less than " << markerGraph.edges.size() << ".";
        return;
    }

    // Access the edge.
    const MarkerGraph::Edge& edge = markerGraph.edges[edgeId];
    array<MarkerGraph::VertexId, 2> vertexIds = {edge.source, edge.target};
    const size_t markerCount = edge.coverage;

    // The marker intervals of this edge.
    const MemoryAsContainer<MarkerInterval> markerIntervals = markerGraph.edgeMarkerIntervals[edgeId];
    CZI_ASSERT(markerIntervals.size() == markerCount);

    // The length of each marker sequence.
    const size_t k = assemblerInfo->k;



    // Extract the sequences and repeat counts for each marker interval.
    // This includes sequences and repeat counts for the flanking markers.
    vector< vector<Base> > sequences(markerCount);
    vector< vector<uint8_t> > repeatCounts(markerCount);
    for(size_t j=0; j!=markerCount; j++) {
        const MarkerInterval& markerInterval = markerIntervals[j];
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const auto orientedReadMarkers = markers[orientedReadId.getValue()];

        // Get the two markers.
        const CompressedMarker& marker0 = orientedReadMarkers[markerInterval.ordinals[0]];
        const CompressedMarker& marker1 = orientedReadMarkers[markerInterval.ordinals[1]];

        // Get the position range, including the flanking markers.
        const uint32_t positionBegin = marker0.position;
        const uint32_t positionEnd = marker1.position + uint32_t(k);

        // Store the bases and repeat counts in this interval.
        for(uint32_t position=positionBegin; position!=positionEnd; ++position) {
            Base base;
            uint8_t repeatCount;
            tie(base, repeatCount) = getOrientedReadBaseAndRepeatCount(orientedReadId, position);
            sequences[j].push_back(base);
            repeatCounts[j].push_back(repeatCount);
        }
    }



    // Initialize a spoa multiple sequence alignment.
    const spoa::AlignmentType alignmentType = spoa::AlignmentType::kNW;
    const int8_t match = 1;
    const int8_t mismatch = -1;
    const int8_t gap = -1;
    auto alignmentEngine = spoa::createAlignmentEngine(alignmentType, match, mismatch, gap);
    auto alignmentGraph = spoa::createGraph();


    // Add all the sequences to the alignment,
    // including the flanking markers.
    string sequenceString;
    vector<string> msa;
    for(size_t j=0; j!=markerCount; j++) {

        // Get the sequence.
        sequenceString.clear();
        for(const Base base: sequences[j]) {
            sequenceString.push_back(base.character());
        }

        // Add it to the alignment.
        auto alignment = alignmentEngine->align(sequenceString, alignmentGraph);
        alignmentGraph->add_alignment(alignment, sequenceString);
    }

    // Use spoa to compute the multiple sequence alignment.
    msa.clear();
    alignmentGraph->generate_multiple_sequence_alignment(msa);

    // The length of the alignment, including gaps.
    const size_t alignmentLength = msa.front().size();



    // Compute coverage for each base (including "-")
    // at each position of the alignment.
    vector< array<size_t, 5> > baseCoverage(alignmentLength);
    for(size_t i=0; i<alignmentLength; i++) {
        fill(baseCoverage[i].begin(), baseCoverage[i].end(), 0);
        for(size_t j=0; j<markerCount; j++) {
            const AlignedBase base = AlignedBase::fromCharacter(msa[j][i]);
            ++(baseCoverage[i][base.value]);
        }
    }



    // Use the consensus caller to compute the consensus base
    // and repeat count at each position of the alignment.
    vector<AlignedBase> consensusBase(alignmentLength);
    vector<size_t> consensusRepeatCount(alignmentLength);
    vector<size_t> positions(markerCount, 0);
    for(size_t i=0; i<alignmentLength; i++) {
        Coverage coverage;
        for(size_t j=0; j<markerCount; j++) {
            const Strand strand = markerIntervals[j].orientedReadId.getStrand();
            const AlignedBase base = AlignedBase::fromCharacter(msa[j][i]);
            if(base.isGap()) {
                coverage.addRead(base, strand, 0);  // Gap always gets repeat count of 0.
            } else {
                const AlignedBase base = AlignedBase(sequences[j][positions[j]]);
                const size_t repeatCount = repeatCounts[j][positions[j]];
                ++(positions[j]);
                coverage.addRead(base, strand, repeatCount);
            }

        }
        const Consensus consensus = (*consensusCaller)(coverage);
        consensusBase[i] = consensus.base;
        consensusRepeatCount[i] = consensus.repeatCount;
    }



    // Count the number of gaps in the consensus bases.
    const size_t consensusGapCount =
        count(consensusBase.begin(), consensusBase.end(), AlignedBase::gap());

    // Compute concordant and discordant base coverage at each position.
    vector<size_t> concordantBaseCoverage(alignmentLength, 0);
    vector<size_t> discordantBaseCoverage(alignmentLength, 0);
    for(size_t i=0; i<alignmentLength; i++) {
        for(size_t b=0; b<5; b++) {
            const size_t coverage = baseCoverage[i][b];
            if(consensusBase[i] == AlignedBase::fromInteger(b)) {
                concordantBaseCoverage[i] += coverage;
            } else {
                discordantBaseCoverage[i] += coverage;
            }
        }
    }



    // Page title.
    html << "<h1>Marker graph edge "<< edgeId << "</h1>";

    // Basic information about this edge.
    html << "<p>Source vertex " << vertexIds[0];
    html << ", target vertex " << vertexIds[1] << ".";
    html << "<p>Edge coverage (number of markers) is " << markerCount << ".";

    // Begin the table.
    html <<
        "<table>"
        "<tr>"
        "<th>Oriented<br>read"
        "<th>Ordinal0"
        "<th>Ordinal1"
        "<th>Sequence<br>(includes flanking markers)"
        "<th>Repeat<br>counts";




    // Write one row for each oriented read.
    for(size_t j=0; j<markerCount; j++) {
        const MarkerInterval& markerInterval = markerIntervals[j];
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;
        const ReadId readId = orientedReadId.getReadId();
        const Strand strand = orientedReadId.getStrand();

        // Oriented read id.
        html <<
            "<tr>"
            "<td class=centered>"
            "<a href='exploreRead"
            "?readId=" << readId <<
            "&strand=" << strand << "'>" <<
            orientedReadId << "</a>";

        // Ordinals.
        for(size_t m=0; m<2; m++) {
            html <<
                "<td class=centered>" <<
                "<a href='exploreRead"
                "?readId=" << readId <<
                "&strand=" << strand <<
                "&highlightMarker=" << markerInterval.ordinals[m] <<
                "'>" <<
                markerInterval.ordinals[m] << "</a>";
        }

        // Sequence, written aligned.
        html << "<td class=centered style='font-family:monospace'>";
        size_t position = 0;
        for(size_t i=0; i<alignmentLength; i++) {
            const char alignmentCharacter = msa[j][i];
            if(alignmentCharacter == '-') {
                html << "-";
            } else {
                if(sequences[j].size() > 2*k && position==k) {
                    html << "<span style='background-color:LightGreen'>";
                }
                const Base base = sequences[j][position];
                ++position;
                CZI_ASSERT(base == Base::fromCharacter(alignmentCharacter));
                html << base;
                if(sequences[j].size() > 2*k && position==sequences[j].size()-k) {
                    html << "</span>";
                }
            }
        }
        CZI_ASSERT(position == sequences[j].size());
        CZI_ASSERT(position == repeatCounts[j].size());



        // Repeat counts, also written aligned.
        html << "<td class=centered style='font-family:monospace'>";
        position = 0;
        for(size_t i=0; i<alignmentLength; i++) {
            const char alignmentCharacter = msa[j][i];
            if(alignmentCharacter == '-') {
                html << "-";
            } else {
                if(sequences[j].size() > 2*k && position==k) {
                    html << "<span style='background-color:LightGreen'>";
                }
                const Base base = sequences[j][position];
                const uint8_t repeatCount = repeatCounts[j][position];
                ++position;
                CZI_ASSERT(base == Base::fromCharacter(alignmentCharacter));
                if(repeatCount < 10) {
                    html << int(repeatCount);
                } else {
                    html << "*";
                }
                if(sequences[j].size() > 2*k && position==sequences[j].size()-k) {
                    html << "</span>";
                }
            }
        }
        CZI_ASSERT(position == sequences[j].size());
        CZI_ASSERT(position == repeatCounts[j].size());
    }



    // Consensus bases and repeat counts.
    html <<
        "<tr>"
        "<th colspan=3 class=left>Consensus (run-length)"
        "<td class=centered style='font-family:monospace'>";
    size_t position = 0;
    for(size_t i=0; i<alignmentLength; i++) {
        if(consensusBase[i].isGap()) {
            html << "-";
        } else {
            if(alignmentLength-consensusGapCount > 2*k && position==k) {
                html << "<span style='background-color:LightGreen'>";
            }
            html << consensusBase[i];
            ++position;
            if(alignmentLength-consensusGapCount > 2*k && position==alignmentLength-consensusGapCount-k) {
                html << "</span>";
            }
        }
    }
    CZI_ASSERT(position == alignmentLength-consensusGapCount);
    html << "<td class=centered style='font-family:monospace'>";
    position = 0;
    for(size_t i=0; i<alignmentLength; i++) {
        if(consensusBase[i].isGap()) {
            html << "-";
        } else {
            if(alignmentLength-consensusGapCount > 2*k && position==k) {
                html << "<span style='background-color:LightGreen'>";
            }
            const size_t repeatCount = consensusRepeatCount[i];
            if(repeatCount < 10) {
                html << repeatCount;
            } else {
                html << "*";
            }
            ++position;
            if(alignmentLength-consensusGapCount > 2*k && position==alignmentLength-consensusGapCount-k) {
                html << "</span>";
            }
        }
    }
    CZI_ASSERT(position == alignmentLength-consensusGapCount);



    // Concordant base coverage.
    html <<
        "<tr>"
        "<th colspan=3 class=left>Concordant base coverage"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<alignmentLength; i++) {
        const size_t coverage = concordantBaseCoverage[i];
        if(coverage == 0) {
            html << ".";
        } else if(coverage < 10) {
            html << coverage;
        } else {
            html << "*";
        }
    }
    html << "<td>";



    // Discordant base coverage.
    html <<
        "<tr>"
        "<th colspan=3 class=left>Discordant base coverage"
        "<td class=centered style='font-family:monospace'>";
    for(size_t i=0; i<alignmentLength; i++) {
        const size_t coverage = discordantBaseCoverage[i];
        if(coverage == 0) {
            html << ".";
        } else if(coverage < 10) {
            html << coverage;
        } else {
            html << "*";
        }
    }
    html << "<td>";



    // Base coverage.
    for(size_t b=0; b<5; b++) {
        html <<
            "<tr>"
            "<th colspan=3 class=left>Coverage for " << AlignedBase::fromInteger(b) <<
            "<td class=centered style='font-family:monospace'>";
        for(size_t i=0; i<alignmentLength; i++) {
            const size_t coverage = baseCoverage[i][b];
            if(coverage == 0) {
                html << ".";
            } else if(coverage < 10) {
                html << coverage;
            } else {
                html << "*";
            }
        }
        html << "<td>";
    }



    // Raw consensus sequence.
    html <<
        "<tr>"
        "<th colspan=3 class=left>Consensus (raw)"
        "<td class=centered style='font-family:monospace'>";
    position = 0;
    for(size_t i=0; i<alignmentLength; i++) {
        const AlignedBase base = consensusBase[i];

        if(base.isGap()) {
            continue;
        }

        if(alignmentLength-consensusGapCount > 2*k && position==k) {
            html << "<span style='background-color:LightGreen'>";
        }

        const size_t repeatCount = consensusRepeatCount[i];
        for(size_t k=0; k<repeatCount; k++) {
            html << base;
        }
        ++position;

        if(alignmentLength-consensusGapCount > 2*k && position==alignmentLength-consensusGapCount-k) {
            html << "</span>";
        }
    }
    html << "<td>";



    // End the table.
    html << "</table>";
}

#endif

