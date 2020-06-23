// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "LocalReadGraph.hpp"
#include "platformDependent.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "chrono.hpp"
#include "iterator.hpp"


void Assembler::exploreReadGraph(
    const vector<string>& request,
    ostream& html)
{
    if(readGraph.edges.isOpen && readGraph.connectivity.isOpen()) {
        exploreUndirectedReadGraph(request, html);
    } else if(directedReadGraph.isOpen()) {
        exploreDirectedReadGraph(request, html);
    } else {
        html << "The read graph is not available." << endl;
    }

}


bool parseCommaSeparatedReadIDs(string& commaSeparatedReadIDs, vector<OrientedReadId>& readIds, ostream& html){
    readIds.clear();
    string token;

    for (auto& c: commaSeparatedReadIDs){
        if (c == ','){
            try {
                OrientedReadId readID(token);
                readIds.emplace_back(readID);
                token.clear();
            }
            catch(exception& e){
                html << "<p>Invalid read id or read strand: '" << token << "'</p>";
                html << "<p>" << e.what() << "</p>";
                return false;
            }
        }
        else{
            token += c;
        }
    }

    // Place final token which may not have a comma after it
    if (not token.empty()) {
        try {
            OrientedReadId readID(token);
            readIds.emplace_back(readID);
        }
        catch(exception& e){
            html << "<p>Invalid read id or read strand: '" << token << "'</p>";
            html << "<p>" << e.what() << "</p>";
            return false;
        }
    }

    return true;
}


void Assembler::exploreUndirectedReadGraph(
    const vector<string>& request,
    ostream& html) {
    // Get the parameters.
    vector<OrientedReadId> readIds;
    string readIDsString;
    const bool readIdsArePresent = getParameterValue(request, "readId", readIDsString);
    const bool readStringsAreValid = parseCommaSeparatedReadIDs(readIDsString, readIds, html);

    uint32_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    uint32_t maxTrim = httpServerData.assemblerOptions->alignOptions.maxTrim;
    getParameterValue(request, "maxTrim", maxTrim);

    string allowChimericReadsString;
    const bool allowChimericReads = getParameterValue(request, "allowChimericReads", allowChimericReadsString);

    string allowCrossStrandEdgesString;
    const bool allowCrossStrandEdges = getParameterValue(request, "allowCrossStrandEdges", allowCrossStrandEdgesString);

    string layoutMethod = "sfdp";
    getParameterValue(request, "layoutMethod", layoutMethod);

    uint32_t sizePixels = 1200;
    getParameterValue(request, "sizePixels", sizePixels);

    double vertexScalingFactor = 2.;
    getParameterValue(request, "vertexScalingFactor", vertexScalingFactor);

    double edgeThicknessScalingFactor = 6.;
    getParameterValue(request, "edgeThicknessScalingFactor", edgeThicknessScalingFactor);

    double edgeArrowScalingFactor = 1.;
    getParameterValue(request, "edgeArrowScalingFactor", edgeArrowScalingFactor);

    string format = "svg";
    getParameterValue(request, "format", format);

    double timeout = 30;
    getParameterValue(request, "timeout", timeout);

    string addBlastAnnotationsString;
    const bool addBlastAnnotations = getParameterValue(request, "addBlastAnnotations", addBlastAnnotationsString);

    string saveDotFileString;
    const bool saveDotFile = getParameterValue(request, "saveDotFile", saveDotFileString);


    // Write the form.
    html <<
         "<h3>Display a local subgraph of the <a href='docs/ReadGraph.html'>global read graph</a></h3>"
         "<form>"
         "<div style='clear:both; display:table;'>"
         "<div style='float:left;margin:10px;'>"
         "<table>"

         "<tr title='Read id between 0 and " << reads.size() - 1 << "'>"
         "<td style=\"white-space:pre-wrap; word-wrap:break-word\">"
         "Start vertex reads\n"
         "The oriented read should be in the form <code>readId-strand</code>\n"
         "where strand is 0 or 1. For example, <code>\"1345871-1</code>\".\n"
         "To add multiple start points, use a comma separator."
         "<td><input type=text required name=readId size=8 style='text-align:center'"
         << (readIdsArePresent ? ("value='" + readIDsString + "'") : "") <<
         ">";


    html <<
         "<tr title='Maximum distance from start vertex (number of edges)'>"
         "<td>Maximum distance"
         "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
         " value='" << maxDistance <<
         "'>"

         "<tr title='Maximum trim (markers) used to define containment'>"
         "<td>Maximum trim"
         "<td><input type=text required name=maxTrim size=8 style='text-align:center'"
         " value='" << maxTrim <<
         "'>"

         "<tr title='Allow reads marked as chimeric to be included in the local read graph.'>"
         "<td>Allow chimeric reads"
         "<td class=centered><input type=checkbox name=allowChimericReads" <<
         (allowChimericReads ? " checked" : "") <<
         ">"

         "<tr title='Allow edges that skip across strands.'>"
         "<td>Allow cross-strand edges"
         "<td class=centered><input type=checkbox name=allowCrossStrandEdges" <<
         (allowCrossStrandEdges ? " checked" : "") <<
         ">"

         "<tr>"
         "<td>Layout method"
         "<td class=centered>"
         "<input type=radio required name=layoutMethod value='sfdp'" <<
         (layoutMethod == "sfdp" ? " checked=on" : "") <<
         ">sfdp"
         "<br><input type=radio required name=layoutMethod value='fdp'" <<
         (layoutMethod == "fdp" ? " checked=on" : "") <<
         ">fdp"
         "<br><input type=radio required name=layoutMethod value='neato'" <<
         (layoutMethod == "neato" ? " checked=on" : "") <<
         ">neato"

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
         "<input type=radio required name=format value='svg'" <<
         (format == "svg" ? " checked=on" : "") <<
         ">svg"
         "<br><input type=radio required name=format value='png'" <<
         (format == "png" ? " checked=on" : "") <<
         ">png"

         "<tr title='Maximum time (in seconds) allowed for graph creation and layout'>"
         "<td>Timeout (seconds) for graph layout"
         "<td><input type=text required name=timeout size=8 style='text-align:center'" <<
         " value='" << timeout <<
         "'>"

         "<tr title='Add to each vertex tooltip summary information on the best alignment to the reference'>"
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
         "</div>"
         "</div>"
         "<br><input type=submit value='Display'>"
         "</form>";


    // If any necessary values are missing, stop here.
    if (not readIdsArePresent) {
        return;
    }

    // If there was a failure to parse comma separated readId-strand tokens, stop here
    if (not readStringsAreValid) {
        return;
    }

    // Validity checks.
    for (auto &readId: readIds){
        if (readId.getReadId() > reads.size()) {
            html << "<p>Invalid read id " << readId;
            html << ". Must be between 0 and " << reads.size() - 1 << ".";
            return;
        }
    }



    // Create the local read graph.
    LocalReadGraph graph;
    const auto createStartTime = steady_clock::now();
    if(!createLocalReadGraph(readIds,
        maxDistance, allowChimericReads, allowCrossStrandEdges, maxTrim, timeout, graph)) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    html << "<p>The local read graph has " << num_vertices(graph);
    html << " vertices and " << num_edges(graph) << " edges.";
    const auto createFinishTime = steady_clock::now();



    // Add Blast annotations, if requested.
    if(addBlastAnnotations) {

        // Create a fasta file containing the sequences of all the oriented reads
        // in the local read graph.
        const string uuid = to_string(boost::uuids::random_generator()());
        const string fastaFileName = tmpDirectory() + uuid + ".fa";
        ofstream fastaFile(fastaFileName);
        BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
            const LocalReadGraphVertex& vertex = graph[v];
            const vector<Base> sequence = getOrientedReadRawSequence(vertex.orientedReadId);
            const auto readName = readNames[vertex.orientedReadId.getReadId()];
            fastaFile << ">" << vertex.orientedReadId << " ";
            copy(readName.begin(), readName.end(), ostream_iterator<char>(fastaFile));
            fastaFile << "\n";
            copy(sequence.begin(), sequence.end(), ostream_iterator<Base>(fastaFile));
            fastaFile << "\n";
        }

        // Create the blast command and run it.
        const string blastOptions =
            "-outfmt '10 qseqid sseqid sstart send' "
            "-evalue 1e-200 -max_hsps 1 -max_target_seqs 1 "
            "-num_threads " + to_string(std::thread::hardware_concurrency());
        const string blastOutputFileName = tmpDirectory() + uuid + ".txt";
        const string blastErrFileName = tmpDirectory() + uuid + ".errtxt";
        const string command = "blastn -task megablast -subject " + httpServerData.referenceFastaFileName +
            " -query " + fastaFileName + " 1>" + blastOutputFileName + " 2>" + blastErrFileName +
            " " + blastOptions;
        cout << timestamp << "Running Blast command." << endl;
        ::system(command.c_str());
        cout << timestamp << "Blast command completed." << endl;

        // Copy any error output to html.
        if(filesystem::fileSize(blastErrFileName)) {
            ifstream blastErrFile(blastErrFileName);
            html << "<pre style='font-size:10px'>";
            html << blastErrFile.rdbuf();
            html << "</pre>";
            blastErrFile.close();
        }

        // Store alignments, keyed by OrientedReadId.
        // For each OrientedReadId we store all the alignments we found,
        // already tokenized.
        using Separator = boost::char_separator<char>;
        using Tokenizer = boost::tokenizer<Separator>;
        const Separator separator(",");
        std::map<OrientedReadId, vector< vector<string> > > alignments;
        ifstream blastOutputFile(blastOutputFileName);
        string line;
        vector<string> tokens;
        while(true) {

            // Get a line.
            string line;
            std::getline(blastOutputFile, line);
            if(!blastOutputFile) {
                break;
            }

            // Tokenize it.
            Tokenizer tokenizer(line, separator);
            tokens.clear();
            tokens.insert(tokens.begin(), tokenizer.begin(), tokenizer.end());

            // Extract the OrientedReadId.
            SHASTA_ASSERT(!tokens.empty());
            const OrientedReadId orientedReadId = OrientedReadId(tokens.front());

            // Store it.
            alignments[orientedReadId].push_back(tokens);
        }

        // Remove the files we created.
        filesystem::remove(fastaFileName);
        filesystem::remove(blastOutputFileName);
        filesystem::remove(blastErrFileName);

        // Now store the alignments as additional text in the vertices tooltips.
        BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
            LocalReadGraphVertex& vertex = graph[v];
            const auto& vertexAlignments = alignments[vertex.orientedReadId];
            for(const auto& alignment: vertexAlignments) {
                SHASTA_ASSERT(alignment.size() == 4);
                vertex.additionalToolTipText += " " + alignment[1] + ":" + alignment[2] + "-" + alignment[3];
            }
        }
    }


    // Write it out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";

    graph.write(dotFileName,
                layoutMethod,
                maxDistance,
                vertexScalingFactor,
                edgeThicknessScalingFactor,
                edgeArrowScalingFactor,
                httpServerData.assemblerOptions->alignOptions.maxTrim);



    // Display the graph in svg format.
    if(format == "svg") {
        const string command =
            timeoutCommand() + " " + to_string(timeout - seconds(createFinishTime - createStartTime)) +
            " dot -O -T svg " + dotFileName +
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
        if(!saveDotFile) {
            filesystem::remove(dotFileName);
        }



        // Write a title.
        html <<
             "<h1 style='line-height:10px'>Read graph near oriented read(s) " << readIDsString << "</h1>"
             "Color legend: "
             "<span style='background-color:green'>start vertex</span> "
             "<span style='background-color:cyan'>vertices at maximum distance (" << maxDistance <<
             ") from the start vertex</span> "
             ".<br>";

        if(saveDotFile) {
            html << "<p>Graphviz dot file saved as " << dotFileName << "<br>";
        }

        // Allow manually highlighting selected vertices.
        html << R"stringDelimiter(
            <script>
            function highlight_vertex()
            {
                vertex = document.getElementById("highlight").value;
                document.getElementById("highlight").value = "";
                element = document.getElementById("Vertex-" + vertex);
                ellipse = element.children[1].children[0].children[0];
                ellipse.setAttribute("fill", "#ff00ff");
                ellipse.setAttribute("stroke", "#ff00ff");
            }
            </script>
            <p>
            <input id=highlight type=text onchange="highlight_vertex()" size=10>
            Enter an oriented read to highlight, then press Enter. The oriented read should be
            in the form <code>readId-strand</code> where strand is 0 or 1 (for example, <code>"1345871-1</code>").
            To highlight multiple oriented reads, enter them one at a time in the same way.
            <p>
            )stringDelimiter";

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
            " dot -O -T png -T cmapx " + dotFileName +
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
