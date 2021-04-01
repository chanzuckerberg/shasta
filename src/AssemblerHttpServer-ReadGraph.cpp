#ifdef SHASTA_HTTP_SERVER

// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "LocalReadGraph.hpp"
#include "orderPairs.hpp"
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
    } else {
        html << "The read graph is not available." << endl;
    }

}


bool Assembler::parseCommaSeparatedReadIDs(string& commaSeparatedReadIds, vector<OrientedReadId>& readIds, ostream& html){
    readIds.clear();
    string token;

    for (auto& c: commaSeparatedReadIds){
        if (c == ','){
            try {
                OrientedReadId readId(token);
                readIds.push_back(readId);
                token.clear();
            }
            catch(exception& e){
                html << "<p>Invalid oriented read id: '" << token << "'</p>";
                html << "<p>Specify one or more comma separated oriented read ids of the form ReadId-strand "
                    "where strand is 0 or 1. For example: 757-1,1048-0";
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
            OrientedReadId readId(token);
            readIds.push_back(readId);
        }
        catch(exception& e){
            html << "<p>Invalid oriented read id: '" << token << "'</p>";
            html << "<p>Specify one or more comma separated oriented read ids of the form ReadId-strand "
                "where strand is 0 or 1. For example: 757-1,1048-0";
            return false;
        }
    }

    return true;
}



// Write to html buttons to resize the svg locally (in the browser).
// This assumes that the page contains a single svg object.
void Assembler::addScaleSvgButtons(ostream& html)
{
    html << R"stringDelimiter(
        <script>
        function svgLarger()
        {
            var element = document.getElementsByTagName("svg")[0];
            element.setAttribute("width", 1.25*element.getAttribute("width"));
            element.setAttribute("height", 1.25*element.getAttribute("height"));
            document.getElementById("largerButton").focus();
        }
        function svgSmaller()
        {
            var element = document.getElementsByTagName("svg")[0];
            element.setAttribute("width", 0.8*element.getAttribute("width"));
            element.setAttribute("height", 0.8*element.getAttribute("height"));
            document.getElementById("smallerButton").focus();
        }
        </script>
        <button type="button" id=largerButton onclick='svgLarger()'>Larger</button>
        &nbsp;
        <button type="button" id=smallerButton onclick='svgSmaller()'>Smaller</button>
        <br>
    )stringDelimiter";    
}



void Assembler::exploreUndirectedReadGraph(
    const vector<string>& request,
    ostream& html) {
    // Get the parameters.
    vector<OrientedReadId> readIds;
    string readIdsString;
    const bool readIdsArePresent = getParameterValue(request, "readId", readIdsString);
    const bool readStringsAreValid = parseCommaSeparatedReadIDs(readIdsString, readIds, html);

    uint32_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    string allowChimericReadsString;
    const bool allowChimericReads = getParameterValue(request, "allowChimericReads", allowChimericReadsString);

    string allowCrossStrandEdgesString;
    const bool allowCrossStrandEdges = getParameterValue(request, "allowCrossStrandEdges", allowCrossStrandEdgesString);

    string layoutMethod = "sfdp";
    getParameterValue(request, "layoutMethod", layoutMethod);

    uint32_t sizePixels = 600;
    getParameterValue(request, "sizePixels", sizePixels);

    double vertexScalingFactor = 1.;
    getParameterValue(request, "vertexScalingFactor", vertexScalingFactor);

    double edgeThicknessScalingFactor = 1.;
    getParameterValue(request, "edgeThicknessScalingFactor", edgeThicknessScalingFactor);

    double timeout = 30;
    getParameterValue(request, "timeout", timeout);

    string addBlastAnnotationsString;
    const bool addBlastAnnotations = getParameterValue(request, "addBlastAnnotations", addBlastAnnotationsString);

    string analyzeString;
    const bool analyze = getParameterValue(request, "analyze", analyzeString);


    // Write the form.
    string readGraphHeading;
    if (httpServerData.docsDirectory.empty()) {
        readGraphHeading = "<h3>Display a local subgraph of the global alignment graph</h3>";
    } else {
        readGraphHeading =
            "<h3>Display a local subgraph of the <a href='docs/ComputationalMethods.html#ReadGraph'>read graph</a></h3>";
    }
    html << readGraphHeading << 
         "<form>"
         "<div style='clear:both; display:table;'>"
         "<div style='float:left;margin:10px;'>"
         "<table>"

         "<tr title='Read id between 0 and " << reads->readCount() - 1 << "'>"
         "<td style=\"white-space:pre-wrap; word-wrap:break-word\">"
         "Start vertex reads\n"
         "The oriented read should be in the form <code>readId-strand</code>\n"
         "where strand is 0 or 1. For example, <code>\"1345871-1</code>\".\n"
         "To add multiple start points, use a comma separator."
         "<td><input type=text required name=readId size=8 style='text-align:center'"
         << (readIdsArePresent ? ("value='" + readIdsString + "'") : "") <<
         ">";


    html <<
         "<tr title='Maximum distance from start vertex (number of edges)'>"
         "<td>Maximum distance"
         "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
         " value='" << maxDistance <<
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

         "<tr>"
         "<td>Perform least square analysis"
         "<td class=centered><input type=checkbox name=analyze" <<
         (analyze ? " checked" : "") <<
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
        if (readId.getReadId() > reads->readCount()) {
            html << "<p>Invalid read id " << readId;
            html << ". Must be between 0 and " << reads->readCount() - 1 << ".";
            return;
        }
    }



    // Create the local read graph.
    LocalReadGraph graph;
    if(!createLocalReadGraph(readIds,
        maxDistance, allowChimericReads, allowCrossStrandEdges, timeout, graph)) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    html << "<p>The local read graph has " << num_vertices(graph);
    html << " vertices and " << num_edges(graph) << " edges.";



    // Add Blast annotations, if requested.
    if(addBlastAnnotations) {

        // Create a fasta file containing the sequences of all the oriented reads
        // in the local read graph.
        const string uuid = to_string(boost::uuids::random_generator()());
        const string fastaFileName = tmpDirectory() + uuid + ".fa";
        ofstream fastaFile(fastaFileName);
        BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
            const LocalReadGraphVertex& vertex = graph[v];
            const vector<Base> sequence = reads->getOrientedReadRawSequence(vertex.orientedReadId);
            const auto readName = reads->getReadName(vertex.orientedReadId.getReadId());
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



    // Analyze the local read graph, if requested.
    if(analyze) {
        analyzeLocalReadGraph(graph);
    }


    // Write a title.
    html <<
         "<h1 style='line-height:10px'>Read graph near oriented read(s) " << readIdsString << "</h1>";


    // Allow manually highlighting selected vertices.
    html << R"stringDelimiter(
        <script>
        function highlight_vertex()
        {
            vertex = document.getElementById("highlight").value;
            document.getElementById("highlight").value = "";
            element = document.getElementById("Vertex-" + vertex);
            element.setAttribute("fill", "#ff00ff");
        }
        </script>
        <p>
        <input id=highlight type=text onchange="highlight_vertex()" size=10>
        Enter an oriented read to highlight, then press Enter. The oriented read should be
        in the form <code>readId-strand</code> where strand is 0 or 1 (for example, <code>"1345871-1</code>").
        To highlight multiple oriented reads, enter them one at a time in the same way.
        <p>
        )stringDelimiter";

    // Buttons to resize the svg locally.
    addScaleSvgButtons(html);

    // Write the graph to svg directly, without using Graphviz rendering.
    ComputeLayoutReturnCode returnCode = graph.computeLayout(layoutMethod, timeout);
    if(returnCode == ComputeLayoutReturnCode::Timeout){
        html << "<p>Timeout exceeded for computing graph layout. Try longer timeout or different parameters.</p>";
    }
    else if (returnCode != ComputeLayoutReturnCode::Success){
        html << "<p>ERROR: graph layout failed </p>";
    }
    else{
        graph.writeSvg("svg",
                       sizePixels,
                       sizePixels,
                       vertexScalingFactor,
                       edgeThicknessScalingFactor,
                       maxDistance,
                       html);
    }



    if(analyze) {
        using vertex_descriptor = LocalReadGraph::vertex_descriptor;
        using edge_descriptor = LocalReadGraph::edge_descriptor;

        // Sort vertices by OrientedReadId.
        vector< pair<OrientedReadId, vertex_descriptor> > sortedVertices;
        BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
            sortedVertices.push_back(make_pair(graph[v].orientedReadId, v));
        }
        sort(sortedVertices.begin(), sortedVertices.end(),
            OrderPairsByFirstOnly<OrientedReadId, vertex_descriptor>());

        // Write least square positions of the vertices.
        html << "<h2>Least square analysis</h2>"
            "<h3>Vertices</h3>"
            "<table><tr><th>Oriented<br>Read Id<th>Least<br>square<br>position";
        for(const auto& p: sortedVertices) {
            const vertex_descriptor v = p.second;
            const LocalReadGraphVertex& vertex = graph[v];
            html << "<tr><td class=centered>" << vertex.orientedReadId <<
                "<td class=centered>" << vertex.leastSquarePosition;
        }
        html << "</table>";


        // Sort the edges by absolute value of residual.
        vector< pair<double, edge_descriptor> > sortedEdges;
        BGL_FORALL_EDGES(e, graph, LocalReadGraph) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            const double x0 = graph[v0].leastSquarePosition;
            const double x1 = graph[v1].leastSquarePosition;
            const double residual = (x1 - x0) - graph[e].averageAlignmentOffset;
            sortedEdges.push_back(make_pair(abs(residual), e));
        }
        sort(sortedEdges.begin(), sortedEdges.end(),
            OrderPairsByFirstOnlyGreater<double, edge_descriptor>());


        // Write edge information.
        const auto oldPrecision = html.precision(1);
        const auto oldFlags = html.setf(std::ios_base::fixed, std::ios_base::floatfield);

        html <<
            "<h3>Edges</h3>"
            "<table><tr>"
            "<th>Id0"
            "<th>Id1"
            "<th>Least<br>square<br>position0"
            "<th>Least<br>square<br>position1"
            "<th>Alignment<br>offset"
            "<th>Least<br>square<br>offset"
            "<th>Least<br>square<br>residual";
        for(const auto& p: sortedEdges) {
            const edge_descriptor e = p.second;
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            const double x0 = graph[v0].leastSquarePosition;
            const double x1 = graph[v1].leastSquarePosition;
            const double residual = (x1 - x0) - graph[e].averageAlignmentOffset;
            html << "<tr>"
                "<td class=centered>" << graph[v0].orientedReadId <<
                "<td class=centered>" << graph[v1].orientedReadId <<
                "<td class=centered>" << x0 <<
                "<td class=centered>" << x1 <<
                "<td class=centered>" << graph[e].averageAlignmentOffset <<
                "<td class=centered>" << x1 - x0 <<
                "<td class=centered>" << residual;
        }
        html << "</table>";
        html.precision(oldPrecision);
        html.setf(oldFlags);

    }

}

#endif
