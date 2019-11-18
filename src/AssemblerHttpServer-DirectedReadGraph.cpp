// Shasta.
#include "Assembler.hpp"
#include "LocalDirectedReadGraph.hpp"
#include "platformDependent.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "chrono.hpp"



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

    double vertexScalingFactor = 1.;
    getParameterValue(request, "vertexScalingFactor", vertexScalingFactor);

    double edgeThicknessScalingFactor = 1.;
    getParameterValue(request, "edgeThicknessScalingFactor", edgeThicknessScalingFactor);

    double edgeArrowScalingFactor = 1.;
    getParameterValue(request, "edgeArrowScalingFactor", edgeArrowScalingFactor);

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
        "<td>Vertex scaling factor"
        "<td><input type=text required name=vertexScalingFactor size=8 style='text-align:center'" <<
        " value='" << vertexScalingFactor <<
        "'>"

        "<tr>"
        "<td>Edge thickness scaling factor"
        "<td><input type=text required name=edgeThicknessScalingFactor size=8 style='text-align:center'" <<
        " value='" << edgeThicknessScalingFactor <<
        "'>"

        "<tr>"
        "<td>Edge arrow scaling factor"
        "<td><input type=text required name=edgeArrowScalingFactor size=8 style='text-align:center'" <<
        " value='" << edgeArrowScalingFactor <<
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

    // Create the local subgraph.
    const OrientedReadId orientedReadId(readId, strand);
    LocalDirectedReadGraph graph;
    const auto createStartTime = steady_clock::now();
    if(not directedReadGraph.extractLocalSubgraph(orientedReadId, maxDistance, timeout, graph)) {
        html << "<p>Timeout for graph creation exceeded. "
            "Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    const auto createFinishTime = steady_clock::now();
    html << "<p>The local read graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges.";

    // Write it out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    graph.write(dotFileName, maxDistance, vertexScalingFactor, edgeThicknessScalingFactor, edgeArrowScalingFactor);



    // Display the graph in svg format.
    if(format == "svg") {
        const string command =
            timeoutCommand() + " " + to_string(timeout - seconds(createFinishTime - createStartTime)) +
            " sfdp -O -T svg " + dotFileName +
            " -Gsize=" + to_string(sizePixels/72.);
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
        if(not saveDotFile) {
            filesystem::remove(dotFileName);
        }



        // Write a title and display the graph.
        html <<
            "<h1 style='line-height:10px'>Read graph near oriented read " << orientedReadId << "</h1>"
            "Color legend: "
            "<span style='background-color:green'>start vertex</span> "
            "<span style='background-color:cyan'>vertices at maximum distance (" << maxDistance <<
            ") from the start vertex</span> "
            ".<br>";
        if(saveDotFile) {
            html << "<p>Graphviz dot file saved as " << dotFileName << "<br>";
        }


        // Display the graph.
        const string svgFileName = dotFileName + ".svg";
        ifstream svgFile(svgFileName);
        html << svgFile.rdbuf();
        svgFile.close();



        // Add to each vertex a cursor that shows you can click on it.
        html <<
            "<script>"
            "var vertices = document.getElementsByClassName('node');"
            "for (var i=0;i<vertices.length; i++) {"
            "    vertices[i].style.cursor = 'pointer';"
            "}"
            "</script>";



        // Remove the .svg file.
        filesystem::remove(svgFileName);
    }



    // Display the graph in png format.
    else if(format == "png") {

        // Run graphviz to create the png file and the cmapx file.
        // The cmapx file is used to create links on the image.
        // See here for more information:
        // https://www.graphviz.org/doc/info/output.html#d:imap
        const string command =
            timeoutCommand() + " " + to_string(timeout - seconds(createFinishTime - createStartTime)) +
            " sfdp -O -T png -T cmapx " + dotFileName +
            " -Gsize=" + to_string(sizePixels/72.);
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
            filesystem::remove(dotFileName);
            throw runtime_error("Signal " + to_string(signalNumber) + " while running graph layout command: " + command);
        } else {
            filesystem::remove(dotFileName);
            throw runtime_error("Abnormal status " + to_string(commandStatus) + " while running graph layout command: " + command);

        }
        // Remove the .dot file.
        if(!saveDotFile) {
            filesystem::remove(dotFileName);
        } else {
            html << "<p>Graphviz dot file saved as " << dotFileName;
        }

        // Get the names of the files we created.
        const string pngFileName = dotFileName + ".png";
        const string cmapxFileName = dotFileName + ".cmapx";

        // Create a base64 version of the png file.
        const string base64FileName = pngFileName + ".base64";
        const string base64Command = "base64 " + pngFileName + " > " +
            base64FileName;
        ::system(base64Command.c_str());


        // Write a title.
        html <<
            "<h1 style='line-height:10px'>Read graph near oriented read " << orientedReadId << "</h1>"
            "Color legend: "
            "<span style='background-color:green'>start vertex</span> "
            "<span style='background-color:cyan'>vertices at maximum distance (" << maxDistance <<
            ") from the start vertex</span> "
            ".<br>";

        // Write out the png image.
        html << "<p><img usemap='#G' src=\"data:image/png;base64,";
        ifstream png(base64FileName);
        html << png.rdbuf();
        html << "\"/>";
        ifstream cmapx(cmapxFileName);
        html << cmapx.rdbuf();

        // Remove the files we created.
        filesystem::remove(pngFileName);
        filesystem::remove(cmapxFileName);
        filesystem::remove(base64FileName);
    }



    // If got here, the format string is not one of the ones
    // we support.
    else {
        html << "Invalid format " << format << " specified";
        filesystem::remove(dotFileName);
    }
}
