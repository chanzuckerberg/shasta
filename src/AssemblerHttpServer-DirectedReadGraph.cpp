// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "LocalDirectedReadGraph.hpp"
#include "platformDependent.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "chrono.hpp"
#include "iterator.hpp"



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

    uint32_t minAlignedMarkerCount = httpServerData.assemblerOptions->alignOptions.minAlignedMarkerCount;
    getParameterValue(request, "minAlignedMarkerCount", minAlignedMarkerCount);

    uint32_t maxOffsetAtCenter = 1000000;
    getParameterValue(request, "maxOffsetAtCenter", maxOffsetAtCenter);

    double minAlignedFraction = httpServerData.assemblerOptions->alignOptions.minAlignedFraction;
    getParameterValue(request, "minAlignedFraction", minAlignedFraction);

    float minTransitiveCoverage = 0.;
    getParameterValue(request, "minTransitiveCoverage", minTransitiveCoverage);

    string allowEdgesNotKeptString;
    bool allowEdgesNotKept = getParameterValue(request,
        "allowEdgesNotKept", allowEdgesNotKeptString);

    string excludeConflictEdgesString;
    bool excludeConflictEdges = getParameterValue(request,
        "excludeConflictEdges", excludeConflictEdgesString);

    uint32_t sizePixels = 600;
    getParameterValue(request, "sizePixels", sizePixels);

    double vertexScalingFactor = 1.;
    getParameterValue(request, "vertexScalingFactor", vertexScalingFactor);

    double edgeThicknessScalingFactor = 1.;
    getParameterValue(request, "edgeThicknessScalingFactor", edgeThicknessScalingFactor);

    double edgeArrowScalingFactor = 0.;
    getParameterValue(request, "edgeArrowScalingFactor", edgeArrowScalingFactor);

    string colorEdgeArrowsString;
    const bool colorEdgeArrows = getParameterValue(request, "colorEdgeArrows", colorEdgeArrowsString);

    string dashedContainmentEdgesString;
    const bool dashedContainmentEdges = getParameterValue(request, "dashedContainmentEdges", dashedContainmentEdgesString);

    string vertexColoringMethodString = "ByDistance";
    getParameterValue(request, "vertexColoringMethod", vertexColoringMethodString);
    LocalDirectedReadGraph::VertexColoringMethod vertexColoringMethod;
    if(vertexColoringMethodString == "ByConflictCount") {
        vertexColoringMethod = LocalDirectedReadGraph::VertexColoringMethod::ByConflictCount;
    } else if(vertexColoringMethodString == "ByCluster") {
        vertexColoringMethod = LocalDirectedReadGraph::VertexColoringMethod::ByCluster;
    } else {
        vertexColoringMethod = LocalDirectedReadGraph::VertexColoringMethod::None;
    }

    string highlightConflicting;
    getParameterValue(request, "highlightConflicting", highlightConflicting);

    string format = "svg";
    getParameterValue(request, "format", format);

    double timeout= 30;
    getParameterValue(request, "timeout", timeout);

    double blastAnnotationsSampling = 0.;
    getParameterValue(request, "blastAnnotationsSampling", blastAnnotationsSampling);

    string saveDotFileString;
    const bool saveDotFile = getParameterValue(request, "saveDotFile", saveDotFileString);



    // Write the form.
    html <<
        "<h3>Display a local subgraph of the global read graph</a></h3>"
        "<form>"
        "<div style='clear:both; display:table;'>"
        "<div style='float:left;margin:10px;'>"
        "<table>"

        "<tr title='Read id between 0 and " << reads.size()-1 << "'>"
        "<td>Start vertex read id"
        "<td><input type=text required name=readId size=8 style='text-align:center'"
        << (readIdIsPresent ? ("value='"+to_string(readId)+"'") : "") <<
        ">";

    html << "<tr><td>Start vertex strand<td class=centered>";
        writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);



    html <<

        "<tr title='Maximum distance from start vertex (number of edges)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance <<
        "'>"

        "<tr>"
        "<td>Minimum number of aligned markers"
        "<td><input type=text required name=minAlignedMarkerCount size=8 style='text-align:center'"
        " value='" << minAlignedMarkerCount <<
        "'>"

        "<tr>"
        "<td>Maximum offset at centers (markers)"
        "<td><input type=text required name=maxOffsetAtCenter size=8 style='text-align:center'"
        " value='" << maxOffsetAtCenter <<
        "'>"

        "<tr>"
        "<td>Minimum aligned fraction"
        "<td><input type=text required name=minAlignedFraction size=8 style='text-align:center'"
        " value='" << minAlignedFraction <<
        "'>"

        "<tr>"
        "<td>Include edges not kept for marker graph creation"
        "<td class=centered><input type=checkbox name=allowEdgesNotKept" <<
        (allowEdgesNotKept ? " checked" : "") <<
        ">"

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
        "<td>Color edge arrows by direction"
        "<td class=centered><input type=checkbox name=colorEdgeArrows" <<
        (colorEdgeArrows ? " checked" : "") <<
        ">"

        "<tr>"
        "<td>Draw containment edges dashed"
        "<td class=centered><input type=checkbox name=dashedContainmentEdges" <<
        (dashedContainmentEdges ? " checked" : "") <<
        ">"

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

        "<tr title='Fraction of vertices that should be mapped to the reference using Blast. "
        "Requires Blast to be installed and reference.fa to be present in the assembly directory.'>"
        "<td>Blast annotations sampling factor"
        "<td><input type=text required name=blastAnnotationsSampling size=8 style='text-align:center'" <<
        " value='" << blastAnnotationsSampling <<
        "'>"

        "<tr title='Save the Graphviz dot file representing this local read graph'>"
        "<td>Save the Graphviz dot file"
        "<td class=centered><input type=checkbox name=saveDotFile" <<
        (saveDotFile ? " checked" : "") <<
        ">"
        "</table>"
        "</div>";

    // If the conflict read graph is open, also show the second table,
    // containing parameters to control the display of conflict information.
    if(conflictReadGraph.isOpen()) {
        html <<
            "<div style='float:left;margin:10px;'>"
            "<table>"

            "<tr>"
            "<td>Exclude conflict edges"
            "<td class=centered><input type=checkbox name=excludeConflictEdges" <<
            (excludeConflictEdges ? " checked" : "") <<
            ">"

            "<tr title='Enter oriented reads in the format readId-strand'>"
            "<td>Highlight vertices conflicting with"
            "<td class=centered><input type=text name=highlightConflicting size=20 " <<
            " value=\"" << highlightConflicting << "\""
            ">"

            "<tr>"
            "<td>Vertex coloring method"
            "<td>"
            "<input type=radio required name=vertexColoringMethod value='None'" <<
            (vertexColoringMethod == LocalDirectedReadGraph::VertexColoringMethod::None ? " checked=on" : "") <<
            ">None"
            "<br><input type=radio required name=vertexColoringMethod value='ByConflictCount'" <<
            (vertexColoringMethod == LocalDirectedReadGraph::VertexColoringMethod::ByConflictCount ? " checked=on" : "") <<
            ">By number of conflicts"
            "<br><input type=radio required name=vertexColoringMethod value='ByCluster'" <<
            (vertexColoringMethod == LocalDirectedReadGraph::VertexColoringMethod::ByCluster ? " checked=on" : "") <<
            ">By cluster"

            "</table>"
            "</div>";
    }

    html <<
        "</div>"
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
    if(not directedReadGraph.extractLocalSubgraph(
        orientedReadId, maxDistance,
        minAlignedMarkerCount, maxOffsetAtCenter, minAlignedFraction,
        allowEdgesNotKept,
        excludeConflictEdges,
        timeout, graph)) {
        html << "<p>Timeout for graph creation exceeded. "
            "Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    const auto createFinishTime = steady_clock::now();
    html << "<p>The local read graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges.";



    // If the conflict read graph is available, add
    // conflict information to the vertices.
    if(conflictReadGraph.isOpen()) {

        BGL_FORALL_VERTICES(v, graph, LocalDirectedReadGraph) {
            LocalDirectedReadGraphVertex& vertex = graph[v];
            const OrientedReadId orientedReadId = vertex.orientedReadId;
            const ConflictReadGraph::VertexId cVertexId =
                ConflictReadGraph::getVertexId(orientedReadId);
            const ConflictReadGraphVertex& cVertex =
                conflictReadGraph.getVertex(cVertexId);
            vertex.componentId = cVertex.componentId;
            vertex.color = cVertex.color;
            vertex.wasRemoved = cVertex.wasRemoved;
            vertex.conflictCount = 0;
            for(ConflictReadGraph::EdgeId edgeId: conflictReadGraph.incidentEdges(cVertexId)) {
                if(not conflictReadGraph.getEdge(edgeId).wasRemoved) {
                    ++vertex.conflictCount;
                }
            }
        }



        // Highlight in green vertices specified in highlightConflicting
        // and in red the vertices that have a conflict with them.
        vector<string> tokens;
        boost::algorithm::split(tokens, highlightConflicting,
            boost::algorithm::is_any_of(" "),
            boost::algorithm::token_compress_on);
        for(const string& token: tokens) {

            // Parse the token as an OrientedReadId (readId-Strand).
            vector<string> tokens2;
            boost::algorithm::split(tokens2, token,
                boost::algorithm::is_any_of("-"));
            if(tokens2.size() != 2) {
                continue;
            }
            ReadId readId0;
            Strand strand0;
            try {
                readId0 = boost::lexical_cast<ReadId>(tokens2[0]);
                strand0 = boost::lexical_cast<ReadId>(tokens2[1]);
            } catch(const boost::bad_lexical_cast&) {
                continue;
            }
            if(strand0 > 1) {
                continue;
            }
            const OrientedReadId orientedReadId0(readId0, strand0);

            // Find the corresponding vertex in the local directed read graph.
            const auto v0 = graph.getVertex(orientedReadId0);
            if(v0 != graph.null_vertex()) {
                graph[v0].isConflictingGreen = true;
            }

            // Now, use the ConflictReadGraph to find conflicting oriented reads.
            const ConflictReadGraph::VertexId cv0 = ConflictReadGraph::getVertexId(orientedReadId0);
            if(cv0 > conflictReadGraph.vertices.size()) {
                continue;
            }
            const auto incidentEdges = conflictReadGraph.incidentEdges(cv0);
            for(const auto cEdgeId: incidentEdges) {
                const auto cv1 = conflictReadGraph.otherVertex(cEdgeId, cv0);
                const OrientedReadId orientedReadId1 = ConflictReadGraph::getOrientedReadId(cv1);
                const auto v1 = graph.getVertex(orientedReadId1);
                if(v1 != graph.null_vertex() and not graph[v1].wasRemoved) {
                    graph[v1].isConflictingRed = true;
                }
                            }
        }

    }



    // Add Blast annotations, if requested.
    if(blastAnnotationsSampling > 0.) {
        const uint32_t hashThreshold = uint32_t(blastAnnotationsSampling * double(std::numeric_limits<uint32_t>::max()));

        // Create a fasta file containing the sequences of vertices in the graph.
        // Each vertex is included with probability blastAnnotationsSampling.
        const string uuid = to_string(boost::uuids::random_generator()());
        const string fastaFileName = tmpDirectory() + uuid + ".fa";
        ofstream fastaFile(fastaFileName);
        BGL_FORALL_VERTICES(v, graph, LocalDirectedReadGraph) {
            const LocalDirectedReadGraphVertex& vertex = graph[v];
            const OrientedReadId orientedReadId = vertex.orientedReadId;
            // Only include it with probability with probability blastAnnotationsSampling.
            if(MurmurHash2(&orientedReadId, sizeof(orientedReadId), 117) > hashThreshold) {
                continue;
            }
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
        if(::system(command.c_str()) != 0) {
            // Copy any error output to html.
            if(filesystem::fileSize(blastErrFileName)) {
                ifstream blastErrFile(blastErrFileName);
                html << "<pre style='font-size:10px'>";
                html << blastErrFile.rdbuf();
                html << "</pre>";
                blastErrFile.close();
            }
            html << "Blast command failed. "
                "Make sure Blast is installed and reference.fa is present in the assembly directory.";
            return;
        }
        cout << timestamp << "Blast command completed." << endl;


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
        BGL_FORALL_VERTICES(v, graph, LocalDirectedReadGraph) {
            LocalDirectedReadGraphVertex& vertex = graph[v];
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
        maxDistance,
        vertexScalingFactor,
        edgeThicknessScalingFactor,
        edgeArrowScalingFactor,
        colorEdgeArrows,
        dashedContainmentEdges,
        httpServerData.assemblerOptions->alignOptions.maxTrim,
        vertexColoringMethod);



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



        // Write a title.
        html <<
            "<h1 style='line-height:10px'>Read graph near oriented read " << orientedReadId << "</h1>"
            "Color legend: "
            "<span style='background-color:green'>start vertex</span> "
            "<span style='background-color:cyan'>vertices at maximum distance (" << maxDistance <<
            ") from the start vertex</span> "
            ".<br>";
        if(blastAnnotationsSampling > 0.) {
            html << "Diamond shaped vertex indicates the presence of Blast annotations.<br>";
        }
        if(saveDotFile) {
            html << "<p>Graphviz dot file saved as " << dotFileName << "<br>";
        }



        // Allow manually highlighting selected vertices.
        html << R"stringDelimiter(
<p>
<input id=highlight type=text onchange="highlight()" size=10>
Enter an oriented read to highlight, then press Enter. The oriented read should be
in the form <code>readId-strand</code> where strand is 0 or 1 (for example, <code>"1345871-1</code>").
To highlight multiple oriented reads, enter them one at a time in the same way.
<script>
function highlight()
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



        // Write a title.
        html <<
            "<h1 style='line-height:10px'>Read graph near oriented read " << orientedReadId << "</h1>"
            "Color legend: "
            "<span style='background-color:green'>start vertex</span> "
            "<span style='background-color:cyan'>vertices at maximum distance (" << maxDistance <<
            ") from the start vertex</span> "
            ".<br>";

        // Write out the png image.
        writePngToHtml(html, pngFileName, "#G");
        ifstream cmapx(cmapxFileName);
        html << cmapx.rdbuf();

        // Remove the files we created.
        filesystem::remove(pngFileName);
        filesystem::remove(cmapxFileName);
    }



    // If got here, the format string is not one of the ones
    // we support.
    else {
        html << "Invalid format " << format << " specified";
        filesystem::remove(dotFileName);
    }
}
