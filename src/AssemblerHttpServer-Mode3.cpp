#include "Assembler.hpp"
#include "mode3.hpp"
#include "mode3-LocalAssemblyGraph.hpp"
using namespace shasta;
using namespace mode3;


void Assembler::exploreMode3AssemblyGraph(
    const vector<string>& request,
    ostream& html)
{
    SHASTA_ASSERT(assemblyGraph3Pointer);

    // Get the parameters for the request.
    uint64_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    uint64_t startSegmentId;
    const bool startSegmentIdIsPresent = getParameterValue(request, "startSegmentId", startSegmentId);

    uint32_t sizePixels = 600;
    getParameterValue(request, "sizePixels", sizePixels);



    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<td>Start segment"
        "<td><input type=text required name=startSegmentId size=8 style='text-align:center'"
        " value='" << (startSegmentIdIsPresent ? to_string(startSegmentId) : "") <<
        "'>"

        "<tr>"
        "<td>Maximum distance"
        "<td><input type=text name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance <<
        "'>"

        "<tr>"
        "<td>Graphics size in pixels"
        "<td><input type=text name=sizePixels size=8 style='text-align:center'"
        " value='" << sizePixels <<
        "'>"

        "</table>"
        "<br><input type=submit value='Display'>"
        "</form>";




    if(not startSegmentIdIsPresent) {
        return;
    }


    html << "<h1>Local assembly graph near segment " << startSegmentId << "</h1></p>";


    // Create the local assembly graph and write it to html in svg format.
    mode3::LocalAssemblyGraph localAssemblyGraph(
        *assemblyGraph3Pointer, startSegmentId, maxDistance);
    // localAssemblyGraph.writeGraphviz("LocalAssemblyGraph.dot");
    localAssemblyGraph.writeSvg1(html, sizePixels);

}
