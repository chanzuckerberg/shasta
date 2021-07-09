#ifdef SHASTA_HTTP_SERVER

// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "LocalReadGraph.hpp"
#include "orderPairs.hpp"
#include "platformDependent.hpp"
#include "Reads.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "chrono.hpp"
#include "iterator.hpp"
#include <filesystem>

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
void Assembler::addScaleSvgButtons(ostream& html, uint64_t sizePixels)
{
    html << "<script>var sizePixels = " << sizePixels << ";</script>\n";

    html << R"stringDelimiter(
        <script>
        function svgLarger()
        {
            var element = document.getElementsByTagName("svg")[0];
            width = element.getAttribute("width");
            height = element.getAttribute("height");
            element.setAttribute("width", 1.25*width);
            element.setAttribute("height", 1.25*height);
            document.getElementById("largerButton").focus();
            sizePixels = sizePixels * 1.25;
        }
        function svgSmaller()
        {
            var element = document.getElementsByTagName("svg")[0];
            width = element.getAttribute("width");
            height = element.getAttribute("height");
            element.setAttribute("width", 0.8*width);
            element.setAttribute("height", 0.8*height);
            document.getElementById("smallerButton").focus();
            sizePixels = sizePixels * 0.8
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

    using vertex_descriptor = LocalReadGraph::vertex_descriptor;
    using edge_descriptor = LocalReadGraph::edge_descriptor;

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

    string allowInconsistentAlignmentEdgesString;
    const bool allowInconsistentAlignmentEdges = getParameterValue(request,
        "allowInconsistentAlignmentEdges", allowInconsistentAlignmentEdgesString);

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

    string alignmentAnalysis = "none";
    getParameterValue(request, "alignmentAnalysis", alignmentAnalysis);

    double residualThreshold = 100;
    getParameterValue(request, "residualThreshold", residualThreshold);



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
         "<td>Allow edges with inconsistent alignments"
         "<td class=centered><input type=checkbox name=allowInconsistentAlignmentEdges" <<
         (allowInconsistentAlignmentEdges ? " checked" : "") <<
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
         "<td>Alignment analysis"
         "<td>"
         "<input type=radio name=alignmentAnalysis value=none" <<
         (alignmentAnalysis=="none" ? " checked=on" : "") <<
         ">None"
         "<br><input type=radio name=alignmentAnalysis value=triangles" <<
         (alignmentAnalysis=="triangles" ? " checked=on" : "") <<
         ">Triangles"
         "<br><input type=radio name=alignmentAnalysis value=leastSquare" <<
         (alignmentAnalysis=="leastSquare" ? " checked=on" : "") <<
         ">Least square"

         "<tr title='Edges with least square residual greater than this value will be colored red'>"
         "<td>Residual threshold for alignment analysis"
         "<td><input type=text required name=residualThreshold size=8 style='text-align:center'" <<
         " value='" << residualThreshold <<
         "'>"

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
        maxDistance,
        allowChimericReads, allowCrossStrandEdges, allowInconsistentAlignmentEdges,
        timeout, graph)) {
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
        if(std::filesystem::file_size(blastErrFileName)) {
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



    // Triangle analysis of the local read graph, if requested.
    vector< pair<array<LocalReadGraph::edge_descriptor, 3>, int32_t> > triangles;
    if(alignmentAnalysis == "triangles") {
        triangleAnalysis(graph, triangles);

        // Color all edges that are part of at least one triangle
        // with error greater than residualThreshold.
        BGL_FORALL_EDGES(e, graph, LocalReadGraph) {
            graph[e].color = "hsl(120, 60%, 50%)";
        }
        for(const auto& p: triangles) {
            const auto& edges = p.first;
            const int32_t offsetError = p.second;
            if(abs(offsetError) < residualThreshold) {
                continue;
            }

            for(uint64_t i=0; i<3; i++) {
                graph[edges[i]].color = "hsl(0, 60%, 50%)";
            }
        }
    }



    // Least square analysis of the local read graph, if requested.
    double maxLeastSquareResidual = 0.;
    double rmsLeastSquareResidual = 0.;
    vector<double> singularValues;
    if(alignmentAnalysis == "leastSquare") {
        leastSquareAnalysis(graph, singularValues);

        // Set edge colors bases on residual.
        // The hue is set to green for zero residual,
        // and red for residual greater than residualThreshold.
        uint64_t edgeCount = 0;
        BGL_FORALL_EDGES(e, graph, LocalReadGraph) {
            ++edgeCount;
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            const double x0 = graph[v0].leastSquarePosition;
            const double x1 = graph[v1].leastSquarePosition;
            const double residual = abs((x1 - x0) - graph[e].averageAlignmentOffset);
            maxLeastSquareResidual = max(maxLeastSquareResidual, residual);
            rmsLeastSquareResidual += residual * residual;
            const double ratio = min(1., residual/residualThreshold);
            const int hue = int(120. * (1. - ratio));
            graph[e].color = "hsl(" + to_string(hue) + ", 60%, 50%)";
        }
        rmsLeastSquareResidual = sqrt(rmsLeastSquareResidual / double(edgeCount));
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
    addScaleSvgButtons(html, sizePixels);

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




    if(alignmentAnalysis == "triangles" and not triangles.empty()) {
        html << "<h2>Triangle analysis</h2>"
            "The worst offset error is " << triangles.front().second <<
            " markers."
            "<p><table>"
            "<tr>"
            "<th>0"
            "<th>1"
            "<th>2"
            "<th>Offset01"
            "<th>Offset12"
            "<th>Offset20"
            "<th>OffsetError";

        // Write a table row for each triangle.
        for(const auto& p: triangles) {
            const auto& edges = p.first;
            const int32_t offsetError = p.second;

            const edge_descriptor e01 = edges[0];
            const edge_descriptor e12 = edges[1];
            const edge_descriptor e20 = edges[2];

            const vertex_descriptor v0 = source(e01, graph);
            const vertex_descriptor v1 = source(e12, graph);
            const vertex_descriptor v2 = source(e20, graph);
            SHASTA_ASSERT(target(e20, graph) == v0);

            const OrientedReadId orientedReadId0 = graph[v0].orientedReadId;
            const OrientedReadId orientedReadId1 = graph[v1].orientedReadId;
            const OrientedReadId orientedReadId2 = graph[v2].orientedReadId;

            const ReadGraphEdge& globalEdge01 = readGraph.edges[graph[e01].globalEdgeId];
            const ReadGraphEdge& globalEdge12 = readGraph.edges[graph[e12].globalEdgeId];
            const ReadGraphEdge& globalEdge20 = readGraph.edges[graph[e20].globalEdgeId];

            const uint64_t alignmentId01 = globalEdge01.alignmentId;
            const uint64_t alignmentId12 = globalEdge12.alignmentId;
            const uint64_t alignmentId20 = globalEdge20.alignmentId;

            const AlignmentInfo alignmentInfo01 = alignmentData[alignmentId01].orient(orientedReadId0, orientedReadId1);
            const AlignmentInfo alignmentInfo12 = alignmentData[alignmentId12].orient(orientedReadId1, orientedReadId2);
            const AlignmentInfo alignmentInfo20 = alignmentData[alignmentId20].orient(orientedReadId2, orientedReadId0);

            html << "<tr>"
                "<td class=centered>" << orientedReadId0 <<
                "<td class=centered>" << orientedReadId1 <<
                "<td class=centered>" << orientedReadId2 <<
                "<td class=centered>" << - alignmentInfo01.averageOrdinalOffset <<
                "<td class=centered>" << - alignmentInfo12.averageOrdinalOffset <<
                "<td class=centered>" << - alignmentInfo20.averageOrdinalOffset <<
                "<td class=centered>" << offsetError;
        }
        html << "</table>";
    }



    if(alignmentAnalysis == "leastSquare") {

        // Sort vertices by OrientedReadId.
        vector< pair<OrientedReadId, vertex_descriptor> > sortedVertices;
        BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
            sortedVertices.push_back(make_pair(graph[v].orientedReadId, v));
        }
        sort(sortedVertices.begin(), sortedVertices.end(),
            OrderPairsByFirstOnly<OrientedReadId, vertex_descriptor>());

        // Write least square positions of the vertices.
        const auto oldPrecision = html.precision(1);
        const auto oldFlags = html.setf(std::ios_base::fixed, std::ios_base::floatfield);
        html << "<h2>Least square analysis</h2>"
            "Maximum absolute value of least square residual is " <<
            maxLeastSquareResidual << " markers."
            "<br>Root mean square value of least square residual is " <<
            rmsLeastSquareResidual << " markers."
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
        html.flags(oldFlags);

        html <<
            "<h3>Singular values</h3>"
            "<table><tr>";
        for(const double s: singularValues) {
            html << "<tr><td class=centered>" << s;
        }
        html << "</table>";


    }

}

#endif
