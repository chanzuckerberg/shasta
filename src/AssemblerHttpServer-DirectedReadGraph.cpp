#include "Assembler.hpp"
using namespace shasta;

void Assembler::exploreDirectedReadGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);

    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);

    uint32_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    uint32_t sizePixels = 600;
    getParameterValue(request, "sizePixels", sizePixels);

    string format = "png";
    getParameterValue(request, "format", format);

    double timeout= 30;
    getParameterValue(request, "timeout", timeout);

    string addBlastAnnotationsString;
    const bool addBlastAnnotations = getParameterValue(request, "addBlastAnnotations", addBlastAnnotationsString);

    string saveDotFileString;
    const bool saveDotFile = getParameterValue(request, "saveDotFile", saveDotFileString);



    // Write the form.
    html <<
        "<h3>Display a local subgraph of the global read graph</a></h3>"
        "<form>"

        "<table>"

        "<tr title='Read id between 0 and " << reads.size()-1 << "'>"
        "<td>Read id"
        "<td><input type=text required name=readId size=8 style='text-align:center'"
        << (readIdIsPresent ? ("value='"+to_string(readId)+"'") : "") <<
        ">";

    html << "<tr><td>Strand<td class=centered>";
        writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);



    html <<

        "<tr title='Maximum distance from start vertex (number of edges)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance <<
        "'>"

        "<tr title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "<td>Graphics size in pixels"
        "<td><input type=text required name=sizePixels size=8 style='text-align:center'" <<
        " value='" << sizePixels <<
        "'>"

        "<tr>"
        "<td>Graphics format"
        "<td class=centered>"
        "svg <input type=radio required name=format value='svg'" <<
        (format == "svg" ? " checked=on" : "") <<
        ">"
        "<td class=centered>png <input type=radio required name=format value='png'" <<
        (format == "png" ? " checked=on" : "") <<
        ">"

        "<tr title='Maximum time (in seconds) allowed for graph creation and layout'>"
        "<td>Timeout (seconds) for graph creation and layout"
        "<td><input type=text required name=timeout size=8 style='text-align:center'" <<
        " value='" << timeout <<
        "'>"

        "<tr title='Add to each vertex tooltip summary information on the best alignment "
        "to the reference stored in file reference.fa'>"
        "<td>Add Blast annotations"
        "<td class=centered><input type=checkbox name=addBlastAnnotations" <<
        (addBlastAnnotations ? " checked" : "") <<
        ">"

        "<tr title='Save the Graphviz dot file representing this local read graph'>"
        "<td>Save the Graphviz dot file"
        "<td class=centered><input type=checkbox name=saveDotFile" <<
        (saveDotFile ? " checked" : "") <<
        ">"
        "</table>"

        "<br><input type=submit value='Display'>"
        "</form>";



    // If any necessary values are missing, stop here.
    if(! (readIdIsPresent && strandIsPresent)) {
        return;
    }



    // Validity checks.
    if(readId > reads.size()) {
        html << "<p>Invalid read id " << readId;
        html << ". Must be between 0 and " << reads.size()-1 << ".";
        return;
    }
}
