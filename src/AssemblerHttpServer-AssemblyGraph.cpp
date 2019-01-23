// Shasta.
#include "Assembler.hpp"
#include "LocalAssemblyGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "chrono.hpp"



void Assembler::exploreAssemblyGraph(
    const vector<string>& request,
    ostream& html)
{


    // Get the request parameters.
    LocalAssemblyGraphRequestParameters requestParameters;
    getLocalAssemblyGraphRequestParameters(request, requestParameters);

    // Write the form.
    requestParameters.writeForm(html, assemblyGraph.edges.size());

    // If any required values are missing, stop here.
    if(requestParameters.hasMissingRequiredParameters()) {
        return;
    }

    // Validity check.
    if(requestParameters.edgeId > assemblyGraph.edges.size()) {
        html << "<p>Invalid edge id " << requestParameters.edgeId;
        html << ". Must be between 0 and " << assemblyGraph.edges.size()-1 << " inclusive.";
        return;
    }



    // Create the local assembly graph.
    html << "<h1>Local assembly graph</h1>";
    LocalAssemblyGraph graph(assemblyGraph);
    const auto createStartTime = steady_clock::now();
    if(!extractLocalAssemblyGraph(
        requestParameters.edgeId,
        requestParameters.maxDistance,
        requestParameters.timeout,
        graph)) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    html << "<p>The local assembly graph has " << num_vertices(graph);
    html << " vertices and " << num_edges(graph) << " edges.";

    const auto createFinishTime = steady_clock::now();
    if(seconds(createFinishTime - createStartTime) > requestParameters.timeout) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }

    // Write it out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = "/dev/shm/" + uuid + ".dot";
    graph.write(
        dotFileName,
        requestParameters.maxDistance,
        requestParameters.useDotLayout,
        requestParameters.showVertexLabels,
        requestParameters.showEdgeLabels);



    // Compute graph layout in svg format.
    const string command =
        "timeout " + to_string(requestParameters.timeout - seconds(createFinishTime - createStartTime)) +
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
    // filesystem::remove(dotFileName);

    // Copy the svg to html.
    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << "<br>" << svgFile.rdbuf();
    svgFile.close();

    // Remove the .svg file.
    filesystem::remove(svgFileName);

}



// Extract  from the request the parameters for the display
// of the local assembly graph.
void Assembler::getLocalAssemblyGraphRequestParameters(
    const vector<string>& request,
    LocalAssemblyGraphRequestParameters& parameters) const
{
    parameters.edgeId = 0;
    parameters.edgeIdIsPresent = getParameterValue(
        request, "edgeId", parameters.edgeId);

    parameters.maxDistance = 0;
    parameters.maxDistanceIsPresent = getParameterValue(
        request, "maxDistance", parameters.maxDistance);

    string useDotLayoutString;
    parameters.useDotLayout = getParameterValue(
        request, "useDotLayout", useDotLayoutString);

    string showVertexLabelsString;
    parameters.showVertexLabels = getParameterValue(
        request, "showVertexLabels", showVertexLabelsString);

    string showEdgeLabelsString;
    parameters.showEdgeLabels = getParameterValue(
        request, "showEdgeLabels", showEdgeLabelsString);

    parameters.sizePixels = 600;
    parameters.sizePixelsIsPresent = getParameterValue(
        request, "sizePixels", parameters.sizePixels);

    parameters.timeout = 30;
    parameters.timeoutIsPresent = getParameterValue(
        request, "timeout", parameters.timeout);

}



void Assembler::LocalAssemblyGraphRequestParameters::writeForm(
    ostream& html,
    AssemblyGraph::EdgeId edgeCount) const
{
    html <<
        "<h3>Display a local subgraph of the global assembly graph</h3>"
        "<form>"

        "<table>"

        "<tr title='Edge id between 0 and " << edgeCount << "'>"
        "<td>Edge id"
        "<td><input type=text required name=edgeId size=8 style='text-align:center'"
        << (edgeIdIsPresent ? ("value='"+to_string(edgeId)+"'") : "") <<
        ">"

        "<tr title='Maximum distance from start edge (number of edges)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        << (maxDistanceIsPresent ? ("value='" + to_string(maxDistance)+"'") : " value='6'") <<
        ">"

        "<tr title='Check to use graphviz dot layout. Uncheck to use graphviz sfdp layout.'>"
        "<td>Use dot layout"
        "<td class=centered><input type=checkbox name=useDotLayout"
        << (useDotLayout ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to show vertex labels'>"
        "<td>Label vertices"
        "<td class=centered><input type=checkbox name=showVertexLabels"
        << (showVertexLabels ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to show edge labels'>"
        "<td>Label edges"
        "<td class=centered><input type=checkbox name=showEdgeLabels"
        << (showEdgeLabels ? " checked=checked" : "") <<
        ">"

        "<tr title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "<td>Graphics size in pixels"
        "<td><input type=text required name=sizePixels size=8 style='text-align:center'"
        << (sizePixelsIsPresent ? (" value='" + to_string(sizePixels)+"'") : " value='800'") <<
        ">"

        "<tr>"
        "<td>Timeout (seconds) for graph creation and layout"
        "<td><input type=text required name=timeout size=8 style='text-align:center'"
        << (timeoutIsPresent ? (" value='" + to_string(timeout)+"'") : " value='30'") <<
        ">"
        "</table>"

        "<br><input type=submit value='Display'>"
        "</form>";
}



bool Assembler::LocalAssemblyGraphRequestParameters::hasMissingRequiredParameters() const
{
    return
        !edgeIdIsPresent ||
        !maxDistanceIsPresent ||
        !timeoutIsPresent;
}



void Assembler::exploreAssemblyGraphEdge(const vector<string>& request, ostream& html)
{
    html << "<h2>Show details about an edge of the assembly graph</h2>";

    // Get the request parameters.
    AssemblyGraph::EdgeId edgeId = 0;
    const bool edgeIdIsPresent = getParameterValue(
        request, "edgeId", edgeId);
    string showDetailsString;
    getParameterValue(request, "showDetails", showDetailsString);
    const bool showDetails = (showDetailsString == "on");
    cout << "showDetailsString " << showDetailsString << " " << int(showDetails) << endl;

    // Write the form to get the edge id.
    html <<
        "<form>"
        "<br>Assembly graph edge id: <input type=text name=edgeId" <<
        (edgeIdIsPresent ? (" value='" + to_string(edgeId)) + "'" : "") <<
        " title='Enter an assembly graph edge id between 0 and " << assemblyGraph.edges.size()-1 << " inclusive'"
        "><br>Show assembly details <input type=checkbox name=showDetails" <<
        (showDetails ? " checked=checked" : "") <<
        ">(this can be slow)<br><input type=submit value='Go'>"
        "</form>";

    // If the edge id is missing or invalid, don't do anything.
    if(!edgeIdIsPresent) {
        return;
    }
    if(edgeId >= assemblyGraph.edges.size()) {
        html <<
            "<p>Invalid edge id " << edgeId <<
            ". Enter an assembly graph edge id between 0 and " <<
            assemblyGraph.edges.size()-1 << " inclusive." << endl;
        return;
    }



    // Assemble the sequence and output detailed information to html.
    // Note that this always uses spoa, not marginPhase, regardless of
    // what was done during assembly.
    if(showDetails) {
        vector<Base> sequence;
        vector<uint32_t> repeatCounts;
        const uint32_t markerGraphEdgeLengthThresholdForConsensus = 10000;
        assembleAssemblyGraphEdge(edgeId, markerGraphEdgeLengthThresholdForConsensus, false, sequence, repeatCounts, &html);
    } else {

        // Assembly details were not requested but a global assembly
        // is not available.
        if( !assemblyGraph.sequences.isOpen() ||
            !assemblyGraph.repeatCounts.isOpen()) {
            html << "<p>A global assembly is not available. You can check "
                "the showDetails checkbox to have the sequence computed on the fly. "
                "This could be slow for long sequence.";
            return;
        } else {

            // Assembly details were not requested and a global assembly
            // is available. Get the sequence from the global assembly.

            // Write a title.
            html <<
                "<h1>Assembly graph edge <a href="
                "'exploreAssemblyGraph?edgeId=" << edgeId <<
                "&maxDistance=6&detailed=on&sizePixels=1600&timeout=30'>" <<
                edgeId << "</a></h1>";



            // Assembled run-length sequence.
            const LongBaseSequenceView runLengthSequence = assemblyGraph.sequences[edgeId];
            const MemoryAsContainer<uint8_t> repeatCounts = assemblyGraph.repeatCounts[edgeId];
            CZI_ASSERT(repeatCounts.size() == runLengthSequence.baseCount);
            html << "<p>Assembled run-length sequence (" << runLengthSequence.baseCount <<
                " bases):<br><span style='font-family:courier'>";
            html << runLengthSequence;



            // Write repeat counts in 3 lines each containing
            // a decimal digit of the repeat count.
            uint32_t maxRepeatCount = 0;
            for(size_t j=0; j<repeatCounts.size(); j++) {
                maxRepeatCount = max(maxRepeatCount, uint32_t(repeatCounts[j]));
            }
            html << "<br>";
            for(size_t j=0; j<repeatCounts.size(); j++) {
                html << (repeatCounts[j] % 10);
            }
            if(maxRepeatCount >= 10) {
                html << "<br>";
                for(size_t j=0; j<repeatCounts.size(); j++) {
                    const uint32_t digit = uint32_t(repeatCounts[j] / 10) % 10;
                    if(digit == 0) {
                        html << "&nbsp;";
                    } else {
                        html << digit;
                    }
                }
            }
            if(maxRepeatCount >= 100) {
                html << "<br>";
                for(size_t j=0; j<repeatCounts.size(); j++) {
                    const uint32_t digit = uint32_t(repeatCounts[j] / 100) % 10;
                    if(digit == 0) {
                        html << "&nbsp;";
                    } else {
                        html << digit;
                    }
                }
            }
            html << "</span>";



            // Assembled raw sequence.
            size_t rawSequenceSize = 0;
            for(size_t j=0; j<repeatCounts.size(); j++) {
                const uint32_t repeatCount = repeatCounts[j];
                rawSequenceSize += repeatCount;
            }
            html << "<p>Assembled raw sequence (" << rawSequenceSize <<
                " bases):<br><span style='font-family:courier'>";
            for(size_t i=0; i<runLengthSequence.baseCount; i++) {
                const Base base = runLengthSequence[i];
                const uint8_t repeatCount = repeatCounts[i];
                for(uint8_t k=0; k<repeatCount; k++) {
                    html << base;
                }
            }
            html << "</span>";
        }
    }
}
