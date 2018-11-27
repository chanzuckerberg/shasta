// Shasta.
#include "Assembler.hpp"
#include "LocalAssemblyGraph.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

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
    requestParameters.writeForm(html, assemblyGraph.vertices.size());

    // If any required values are missing, stop here.
    if(requestParameters.hasMissingRequiredParameters()) {
        return;
    }

    // Validity check.
    if(requestParameters.vertexId > assemblyGraph.vertices.size()) {
        html << "<p>Invalid vertex id " << requestParameters.vertexId;
        html << ". Must be between 0 and " << assemblyGraph.vertices.size()-1 << " inclusive.";
        return;
    }



    // Create the local assembly graph.
    LocalAssemblyGraph graph(assemblyGraph);
    const auto createStartTime = steady_clock::now();


    html << "<h1>Local assembly graph</h1>";
}



// Extract  from the request the parameters for the display
// of the local assembly graph.
void Assembler::getLocalAssemblyGraphRequestParameters(
    const vector<string>& request,
    LocalAssemblyGraphRequestParameters& parameters) const
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

    parameters.sizePixels = 800;
    parameters.sizePixelsIsPresent = getParameterValue(
        request, "sizePixels", parameters.sizePixels);

    parameters.timeout = 30;
    parameters.timeoutIsPresent = getParameterValue(
        request, "timeout", parameters.timeout);

}



void Assembler::LocalAssemblyGraphRequestParameters::writeForm(
    ostream& html,
    AssemblyGraph::VertexId vertexCount) const
{
    html <<
        "<h3>Display a local subgraph of the global assembly graph</h3>"
        "<form>"

        "<table>"

        "<tr title='Vertex id between 0 and " << vertexCount << "'>"
        "<td>Vertex id"
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


        "<tr title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "<td>Graphics size in pixels"
        "<td><input type=text required name=sizePixels size=8 style='text-align:center'"
        << (sizePixelsIsPresent ? (" value='" + to_string(sizePixels)+"'") : " value='1600'") <<
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
        !vertexIdIsPresent ||
        !maxDistanceIsPresent ||
        !timeoutIsPresent;
}

