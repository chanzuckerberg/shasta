
// Shasta.
#include "Assembler.hpp"
#include "buildId.hpp"
#include "platformDependent.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#ifdef SHASTA_HTTP_SERVER

#define SHASTA_ADD_TO_FUNCTION_TABLE(name) httpServerData.functionTable[string("/") + #name ] = &Assembler::name



// Associate http keywords with member functions.
void Assembler::fillServerFunctionTable()
{
    httpServerData.functionTable[""]        = &Assembler::exploreSummary;
    httpServerData.functionTable["/"]       = &Assembler::exploreSummary;
    httpServerData.functionTable["/index"]  = &Assembler::exploreSummary;

    SHASTA_ADD_TO_FUNCTION_TABLE(exploreSummary);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreRead);
    SHASTA_ADD_TO_FUNCTION_TABLE(blastRead);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAlignments);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAlignment);
    SHASTA_ADD_TO_FUNCTION_TABLE(computeAllAlignments);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAlignmentGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(displayAlignmentMatrix);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreReadGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerGraphVertex);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerGraphEdge);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerCoverage);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerGraphInducedAlignment);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAssemblyGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAssemblyGraphEdge);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAssemblyGraphEdgesSupport);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreCompressedAssemblyGraph);

}
#undef SHASTA_ADD_TO_FUNCTION_TABLE



void Assembler::processRequest(
    const vector<string>& request,
    ostream& html,
    const BrowserInformation&)
{
    // Process a documentation request.
    const string& keyword = request.front();
    if(keyword.size()>6 && keyword.substr(0, 6)=="/docs/") {

        // Extract the file name.
        const string name = keyword.substr(6);

        // If it contains "/", reject it.
        if(name.find('/') != string::npos) {
            writeHtmlBegin(html);
            html << "Unknown documentation file " << name;
            writeHtmlEnd(html);
            return;
        }

        // Construct the full file name and open it.
        const string fileName = httpServerData.docsDirectory + "/" + name;
        ifstream file(fileName);
        if(!file) {
            writeHtmlBegin(html);
            html << "Could not open " << fileName;
            writeHtmlEnd(html);
        }

        // Send it to html.
        html << "\r\n" << file.rdbuf();
        return;
    }



    // Look up the keyword to find the function that will process this request.
    // Note that the keyword includes the initial "/".
    const auto it = httpServerData.functionTable.find(keyword);
    if(it == httpServerData.functionTable.end()) {
        writeHtmlBegin(html);
        html << "Unsupported keyword " << keyword;
        writeHtmlEnd(html);
        return;
    }


    // We found the keyword. Call the function that processes this keyword.
    // The processing function is only responsible for writing the html body.
    writeHtmlBegin(html);
    try {
        const auto function = it->second;
        (this->*function)(request, html);
    } catch(const std::exception& e) {
        html << "<br><br><span style='color:purple'>" << e.what() << "</span>";
    }
    writeHtmlEnd(html);
}
#endif



void Assembler::writeStyle(ostream& html)
{
    html << R"%(
<style>
    body {
        font-family: Arial;
    }
    pre {
        font-family: courier;
    }
    p, input {
        font-size: 16px;
    }
    h1, h2, h3 {
        color: DarkSlateBlue;
    }
    table {
        border-collapse: collapse;
    }
    th, td {
        border: 1px solid #b8b5c7d9;
        padding: 2px;
    }
    th {
        font-weight: bold;
        text-align: center;
    }
    th.left {
        text-align: left;
    }
    td.centered {
        text-align: center;
    }
    td.right {
        text-align: right;
    }
    td.smaller {
        font-size: smaller;
    }
    a {
        color: DarkSlateBlue;
    }
    ul.navigationMenu {
        list-style-type: none;
        margin: 0px 0px 12px 0px;
        padding: 0;
        overflow: hidden;
        background-color: #404040;
    }
    
    div.navigationButton {
        display: inline-block;
        color: white;
        text-align: center;
        padding: 14px 16px;
        text-decoration: none;
        // min-width: 120px;
    }
    
    .navigationMenuEntry:hover .navigationButton {
        background-color: black;
    }
    
    li.navigationMenuEntry {
        display: inline-block;
    }
    
    .navigationItems {
        display: none;
        position: absolute;
        background-color: DodgerBlue;
        // min-width: 120px;
        box-shadow: 0px 8px 16px 0px rgba(0,0,0,0.2);
        z-index: 1;
    }
    
    a.navigationItem {
        color: black;
        padding: 12px 16px;
        text-decoration: none;
        display: block;
        text-align: left;
    }
    
    .navigationItems a:hover {background-color: SteelBlue}
    
    .navigationMenuEntry:hover .navigationItems {
        display: block;
    }

    input[type=submit] {
        background-color: #89bef2;
        padding: 4px;
    }

    input[type=button] {
        padding: 4px;
    }

    input[type=text], input[type=radio] {
        background-color: #ecf1f0;
        border-width: thin;
    }
</style>
    )%";
}



void Assembler::writeHtmlBegin(ostream& html, bool navigation) const
{
    html <<
        "\r\n"
        "<!DOCTYPE html>"
        "<html>"
        "<head>"
        "<link rel=icon href=docs/CZI-new-logo.png />"
        "<meta charset='UTF-8'>"
        "<title>Shasta assembler</title>";
    writeStyle(html);
    // writeMakeAllTablesSelectable(html);
    html <<
        "</head>"
        ;// "<body onload='makeAllTablesSelectableByDoubleClick()'>";

    if(navigation) {
        writeNavigation(html);
    }
}



void Assembler::writeHtmlEnd(ostream& html) const
{
    html << "</body>";
    html << "</html>";
}



#ifdef SHASTA_HTTP_SERVER
void Assembler::writeMakeAllTablesSelectable(ostream& html) const
{
    html << R"###(
<script>

// Make all tables selectable by double click.
// This must be called after all tables have
// already been created, so it can be called during onload.

// This function is called when the user double clicks on a table.
function selectElement(table)
{
    var selection = window.getSelection();
    selection.removeAllRanges();
    var range = document.createRange();
    range.selectNode(table);
    selection.addRange(range);
}

// Attach the above function to the double click event
// for all tables in the document.
// Also add to each table a title that displays a tooltip 
// explaining that the table can be selected via double click.
function makeAllTablesSelectableByDoubleClick()
{
    var allTables = document.getElementsByTagName("table");
    for (var i=0; i<allTables.length; i++) {
        var table = allTables[i];
        table.ondblclick = function() {selectElement(this);};
        table.setAttribute("title", 
        "Double click to select the entire table. You can then paste it into a spreadsheet.");
    }
}
</script>
    )###";
}
#endif



void Assembler::writeNavigation(ostream& html) const
{
    html << "<ul class=navigationMenu>";

    writeNavigation(html, "Assembly information", {
        {"Summary", "exploreSummary"},
        });
    writeNavigation(html, "Reads", {
        {"Reads", "exploreRead"},
        });
    writeNavigation(html, "Alignments", {
        {"Stored alignments", "exploreAlignments"},
        {"Align two reads", "exploreAlignment"},
        {"Align one read with all", "computeAllAlignments"},
        {"Alignment graph", "exploreAlignmentGraph"},
        {"Alignment matrix", "displayAlignmentMatrix"},
        });
    writeNavigation(html, "Read graph", {
        {"Read graph", "exploreReadGraph"},
        });
    writeNavigation(html, "Marker graph", {
        {"Local marker graph", "exploreMarkerGraph?useBubbleReplacementEdges=on"},
        {"Marker graph vertices", "exploreMarkerGraphVertex"},
        {"Marker graph edges", "exploreMarkerGraphEdge"},
        {"Marker coverage", "exploreMarkerCoverage"},
        {"Induced alignments", "exploreMarkerGraphInducedAlignment"},
        });
    writeNavigation(html, "Assembly graph", {
        {"Local assembly graph", "exploreAssemblyGraph"},
        {"Assembly graph edges", "exploreAssemblyGraphEdge"},
        {"Assembly graph edges support", "exploreAssemblyGraphEdgesSupport"},
        {"Compressed assembly graph", "exploreCompressedAssemblyGraph"},
        });
    writeNavigation(html, "Help", {
        {"Documentation", "docs/index.html"},
        });

    html << "</ul>";
}



void Assembler::writeNavigation(
    ostream& html,
    const string& title,
    const vector<pair <string, string> >& items) const
{
    html <<
        "<li class=navigationMenuEntry>"
        "<div class=navigationButton>" << title << "</div>"
        "<div class=navigationItems>";

    for(const auto& item: items) {
        html << "<a class=navigationItem href=" << item.second << ">" << item.first << "</a>";
    }

    html << "</div></li>";

}



// Write to html an img tag displaying a png file.
void Assembler::writePngToHtml(
    ostream& html,
    const string& pngFileName,
    const string useMap)
{
    // Convert the png file to base64.
    const string base64FileName = tmpDirectory() + to_string(boost::uuids::random_generator()());
    const string base64Command = "base64 " + pngFileName + " > " +
        base64FileName;
    const int errorCode = ::system(base64Command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " +
            to_string(errorCode) + " " + strerror(errorCode) +
            "\nrunning command: " + base64Command);
    }

    // Write the base64 file to html in an img tag.
    html << "<p><img ";
    if(not useMap.empty()) {
        html << "usemap='" << useMap << "'";
    }
    html << " src=\"data:image/png;base64,";
    ifstream png(base64FileName);
    SHASTA_ASSERT(png);
    html << png.rdbuf();
    html << "\"/>";

    // Remove the base64 file.
    filesystem::remove(base64FileName);

}



void Assembler::writeGnuPlotPngToHtml(
    ostream& html,
    int width,
    int height,
    const string& gnuplotCommands)
{

    // Create a file to contain gnuplot commands.
    const string gnuplotFileName = tmpDirectory() + to_string(boost::uuids::random_generator()());
    const string pngFileName = tmpDirectory() + to_string(boost::uuids::random_generator()());
    {
        ofstream gnuplotFile(gnuplotFileName);
        gnuplotFile <<
            "set terminal pngcairo size " << width << "," << height <<
            " font 'Noto Serif'\n"
            "set output '" << pngFileName << "'\n" <<
            gnuplotCommands;
    }

    // Invoke gnuplot.
    const string command = "gnuplot " + gnuplotFileName;
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " +
            to_string(errorCode) + " " + strerror(errorCode) +
            "\nrunning command: " + command);
    }

    // Write the png file to html.
    writePngToHtml(html, pngFileName);
}


#ifdef SHASTA_HTTP_SERVER

// Access all available assembly data, without throwing exceptions
void Assembler::accessAllSoft()
{

    bool allDataAreAvailable = true;

    try {
        accessReadFlags(false);
    } catch(const exception& e) {
        cout << "Read flags are not accessible." << endl;
        allDataAreAvailable = false;
    }


    try {
        accessKmers();
    } catch(const exception& e) {
        cout << "K-mers are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessMarkers();
    } catch(const exception& e) {
        cout << "Markers are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessAlignmentCandidates();
    } catch(const exception& e) {
        cout << "Alignment candidates are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessReadLowHashStatistics();
    } catch(const exception& e) {
        cout << "Read alignment statistics are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessAlignmentData();
    } catch(const exception& e) {
        cout << "Alignments are not accessible." << endl;
        allDataAreAvailable = false;
    }


    // Read graph.
    // Try accessing the undirected one first.
    // if that is not there, try the undirected one.
    try {
        accessReadGraph();
    } catch(const exception& e) {

        // We don't have the undirected read graph. Try the directed one.
        try {
            accessDirectedReadGraphReadOnly();
        } catch(const exception& e) {
            cout << "The read graph is not accessible." << endl;
            allDataAreAvailable = false;
        }

        try {
            accessConflictReadGraph();
        } catch(const exception& e) {
        }
    }



    try {
        accessMarkerGraphVertices();
    } catch(const exception& e) {
        cout << "Marker graph vertices are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessMarkerGraphEdges(false);
    } catch(const exception& e) {
        cout << "Marker graph edges are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessMarkerGraphConsensus();
    } catch(const exception& e) {
        cout << "MarkerGraph graph consensus is accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessAssemblyGraphVertices();
    } catch(const exception& e) {
        cout << "Assembly graph vertices are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessAssemblyGraphEdges();
    } catch(const exception& e) {
        cout << "Assembly graph edges are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessAssemblyGraphEdgeLists();
    } catch(const exception& e) {
        cout << "Assembly graph edge lists are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessAssemblyGraphSequences();
    } catch(const exception& e) {
        cout << "Assembly graph sequences are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessPhasingData();
    } catch(const exception& e) {
    	// Don't threat it as missing because this does not get created in all cases.
        // cout << "Assembly graph sequences are not accessible." << endl;
        // allDataAreAvailable = false;
    }

    if(!allDataAreAvailable) {
        cout << "Not all assembly data are accessible." << endl;
        cout << "Some functionality is not available." << endl;
    }
}



void Assembler::exploreSummary(
    const vector<string>& request,
    ostream& html)
{
    writeAssemblySummaryBody(html);
}

#endif




void Assembler::writeAssemblySummary(ostream& html)
{
    writeHtmlBegin(html, false);
    writeAssemblySummaryBody(html);
    writeHtmlEnd(html);
}

void Assembler::writeAssemblySummaryBody(ostream& html)
{
    using std::setprecision;
    AssemblyGraph& assemblyGraph = *assemblyGraphPointer;


    // Compute the number of run-length k-mers used as markers.
    uint64_t totalRleKmerCount = 0;
    uint64_t markerRleKmerCount = 0;
    for(const auto& tableEntry: kmerTable) {
        if(tableEntry.isRleKmer) {
            ++totalRleKmerCount;
            if(tableEntry.isMarker) {
                ++markerRleKmerCount;
            }
        }
    }

    const uint64_t totalDiscardedReadCount =
        assemblerInfo->discardedInvalidBaseReadCount +
        assemblerInfo->discardedShortReadReadCount +
        assemblerInfo->discardedBadRepeatCountReadCount;
    const uint64_t totalDiscardedBaseCount =
        assemblerInfo->discardedInvalidBaseBaseCount +
        assemblerInfo->discardedShortReadBaseCount +
        assemblerInfo->discardedBadRepeatCountBaseCount;


    html <<
        "<h1>Shasta assembly summary</h1>"
        "<h3>Shasta version</h3>" <<
        buildId() <<



        "<h3>Reads used in this assembly</h3>"
        "<table>"
        "<tr><td>Number of reads"
        "<td class=right>" << assemblerInfo->readCount <<
        "<tr><td>Number of raw sequence bases"
        "<td class=right>" << assemblerInfo->baseCount <<
        "<tr><td>Average read length (for raw read sequence)"
        "<td class=right>" << assemblerInfo->baseCount / assemblerInfo->readCount <<
        "<tr><td>Read N50 (for raw read sequence)"
        "<td class=right>" << assemblerInfo->readN50 <<
        "<tr><td>Number of run-length encoded bases"
        "<td class=right>" << readRepeatCounts.totalSize() <<
        "<tr><td>Average length ratio of run-length encoded sequence over raw sequence"
        "<td class=right>" << setprecision(4) << double(readRepeatCounts.totalSize()) / double(assemblerInfo->baseCount) <<
        "<tr><td>Number of reads flagged as palindromic"
        "<td class=right>" << assemblerInfo->palindromicReadCount <<
        "<tr><td>Number of reads flagged as chimeric"
        "<td class=right>" << assemblerInfo->chimericReadCount <<
        "</table>"
        "<ul>"
        "<li>Here and elsewhere, &quot;raw&quot; refers to the original read sequence, "
        "as opposed to run-length encoded sequence."
        "<li>Reads discarded on input are not included in the above table (see "
        "<a href='#discarded'>below</a>)."
        "<li>See ReadLengthHistogram.csv and Binned-ReadLengthHistogram.csv "
        "for details of the read length distribution of reads used in this assembly.</ul>"


        "<h3 id=discarded>Reads discarded on input</h3>"
        "<table>"
        "<tr><th><th>Reads<th>Bases"
        "<tr><td>Reads discarded on input because they contained invalid bases"
        "<td class=right>" << assemblerInfo->discardedInvalidBaseReadCount <<
        "<td class=right>" << assemblerInfo->discardedInvalidBaseBaseCount <<
        "<tr><td>Reads discarded on input because they were too short"
        "<td class=right>" << assemblerInfo->discardedShortReadReadCount <<
        "<td class=right>" << assemblerInfo->discardedShortReadBaseCount <<
        "<tr><td>Reads discarded on input because they contained repeat counts greater than 255"
        "<td class=right>" << assemblerInfo->discardedBadRepeatCountReadCount <<
        "<td class=right>" << assemblerInfo->discardedBadRepeatCountBaseCount <<
        "<tr><td>Reads discarded on input, total"
        "<td class=right>" <<totalDiscardedReadCount <<
        "<td class=right>" <<totalDiscardedBaseCount <<
        "<tr><td>Fraction of reads discarded on input over total present in input files"
        "<td class=right>" <<
        double(totalDiscardedReadCount) /
        double(totalDiscardedReadCount+assemblerInfo->readCount)
        <<
        "<td class=right>" <<
        double(totalDiscardedBaseCount) /
        double(totalDiscardedBaseCount+assemblerInfo->baseCount)
        <<
        "</table>"
        "<ul><li>Base counts in the above table are raw sequence bases."
        "<li>Here and elsewhere, &quot;raw&quot; refers to the original read sequence, "
        "as opposed to run-length encoded sequence.</ul>"


        "<h3>Marker <i>k</i>-mers</h3>"
        "<table>"
        "<tr><td>Length <i>k</i> of <i>k</i>-mers used as markers"
        "<td class=right>" << assemblerInfo->k <<
        "<tr><td>Total number of <i>k</i>-mers"
        "<td class=right>" << totalRleKmerCount <<
        "<tr><td>Number of <i>k</i>-mers used as markers"
        "<td class=right>" << markerRleKmerCount <<
        "<tr><td>Fraction of <i>k</i>-mers used as markers"
        "<td class=right>" << setprecision(3) << double(markerRleKmerCount) / double(totalRleKmerCount) <<
        "</table>"
        "<ul><li>In the above table, all <i>k</i>-mer counts only include run-length encoded <i>k</i>-mers, "
        "that is, <i>k</i>-mers without repeated bases.</ul>"



        "<h3>Markers</h3>"
        "<table>"
        "<tr><td>Total number of markers on all reads, one strand"
        "<td class=right>" << markers.totalSize()/2 <<
        "<tr><td>Total number of markers on all reads, both strands"
        "<td class=right>" << markers.totalSize() <<
        "<tr><td>Average number of markers per raw base"
        "<td class=right>" << setprecision(4) << double(markers.totalSize()/2)/double(assemblerInfo->baseCount) <<
        "<tr><td>Average number of markers per run-length encoded base"
        "<td class=right>" << setprecision(4) << double(markers.totalSize()/2)/double(readRepeatCounts.totalSize()) <<
        "<tr><td>Average base offset between markers in raw sequence"
        "<td class=right>" << setprecision(4) << double(assemblerInfo->baseCount)/double(markers.totalSize()/2) <<
        "<tr><td>Average base offset between markers in run-length encoded sequence"
        "<td class=right>" << setprecision(4) << double(readRepeatCounts.totalSize())/double(markers.totalSize()/2) <<
        "<tr><td>Average base gap between markers in run-length encoded sequence"
        "<td class=right>" << setprecision(4) <<
        double(readRepeatCounts.totalSize())/double(markers.totalSize()/2) - double(assemblerInfo->k) <<
        "</table>"
        "<ul><li>Here and elsewhere, &quot;raw&quot; refers to the original read sequence, "
        "as opposed to run-length encoded sequence.</ul>"



        "<h3>Alignments</h3>"
        "<table>"
        "<tr><td>Number of alignment candidates found by the LowHash algorithm"
        "<td class=right>" << alignmentCandidates.candidates.size() <<
        "<tr><td>Number of good alignments"
        "<td class=right>" << alignmentData.size() <<
        "<tr><td>Number of good alignments kept in the read graph"
        "<td class=right>" << (directedReadGraph.isOpen() ?
            directedReadGraph.edges.size()/2 :
            readGraph.edges.size()/2) <<
        "</table>"



        "<h3>Read graph</h3>"
        "<table>"
        "<tr><td>Number of vertices"
        "<td class=right>" << (directedReadGraph.isOpen() ?
            directedReadGraph.vertices.size() :
            readGraph.connectivity.size()) <<
        "<tr><td>Number of edges"
        "<td class=right>" << (directedReadGraph.isOpen() ?
            directedReadGraph.edges.size() :
            readGraph.edges.size()) <<
        "</table>"
        "<ul>"
        "<li>The read graph contains both strands. Each read generates two vertices."
        "<li>Isolated reads in the read graph don't contribute to the assembly. "
        "See the table below for a summary of isolated reads in the read graph. "
        "Each isolated read corresponds to two isolated vertices in the read graph, one for each strand."
        "</ul>"
        "<table>"
        "<tr><th><th>Reads<th>Bases"
        "<tr><td>Isolated reads"
        "<td class=centered>" << assemblerInfo->isolatedReadCount <<
        "<td class=centered>" << assemblerInfo->isolatedReadBaseCount <<
        "<tr><td>Non-isolated reads"
        "<td class=centered>" << assemblerInfo->readCount - assemblerInfo->isolatedReadCount <<
        "<td class=centered>" << assemblerInfo->baseCount - assemblerInfo->isolatedReadBaseCount <<
        "<tr><td>Isolated reads fraction"
        "<td class=centered>" << double(assemblerInfo->isolatedReadCount)/double(assemblerInfo->readCount) <<
        "<td class=centered>" << double(assemblerInfo->isolatedReadBaseCount)/double(assemblerInfo->baseCount) <<
        "<tr><td>Non-isolated reads fraction"
        "<td class=centered>" << double(assemblerInfo->readCount-assemblerInfo->isolatedReadCount)/double(assemblerInfo->readCount) <<
        "<td class=centered>" << double(assemblerInfo->baseCount-assemblerInfo->isolatedReadBaseCount)/double(assemblerInfo->baseCount) <<
        "</table>"



        "<h3>Marker graph</h3>"
        "<table>"
        "<tr><td>Total number of vertices"
        "<td class=right>" << markerGraph.vertexCount() <<
        "<tr><td>Total number of edges"
        "<td class=right>" << markerGraph.edges.size() <<
        "<tr><td>Number of vertices that are not isolated after edge removal"
        "<td class=right>" << assemblerInfo->markerGraphVerticesNotIsolatedCount <<
        "<tr><td>Number of edges that were not removed"
        "<td class=right>" << assemblerInfo->markerGraphEdgesNotRemovedCount <<
        "</table>"
        "<ul><li>The marker graph contains both strands.</ul>"



        "<h3>Assembly graph</h3>"
        "<table>"
        "<tr><td>Number of vertices"
        "<td class=right>" << assemblyGraph.vertices.size() <<
        "<tr><td>Number of edges"
        "<td class=right>" << assemblyGraph.edges.size() <<
        "<tr><td>Number of edges assembled"
        "<td class=right>" << assemblerInfo->assemblyGraphAssembledEdgeCount <<
        "</table>"
        "<ul><li>The assembly graph contains both strands.</ul>"



        "<h3>Assembled segments (&quot;contigs&quot;)</h3>"
        "<table>"
        "<tr><td>Number of segments assembled"
        "<td class=right>" << assemblerInfo->assemblyGraphAssembledEdgeCount <<
        "<tr><td>Total assembled segment length"
        "<td class=right>" << assemblerInfo->totalAssembledSegmentLength <<
        "<tr><td>Longest assembled segment length"
        "<td class=right>" << assemblerInfo->longestAssembledSegmentLength <<
        "<tr><td>Assembled segments N<sub>50</sub>"
        "<td class=right>" << assemblerInfo->assembledSegmentN50 <<
        "</table>"
        "<ul><li>Shasta uses GFA terminology "
        "(<i>segment</i> instead of the most common <i>contig</i>). "
        "A contiguous section of assembled sequence can consist of multiple segments, "
        "for example in the presence of heterozygous bubbles."
        "<li>See AssemblySummary.csv for lengths of assembled segments."
        "</ul>"



        "<h3>Performance</h3>"
        "<table>"
        "<tr><td>Elapsed time (seconds)"
        "<td class=right>" << assemblerInfo->assemblyElapsedTimeSeconds <<
        "<tr><td>Elapsed time (minutes)"
        "<td class=right>" << assemblerInfo->assemblyElapsedTimeSeconds/60. <<
        "<tr><td>Elapsed time (hours)"
        "<td class=right>" << assemblerInfo->assemblyElapsedTimeSeconds/3600. <<
        "<tr><td>Average CPU utilization"
        "<td class=right>" << assemblerInfo->averageCpuUtilization <<
        "<tr><td>Peak Memory utilization (bytes)"
        "<td class=right>" <<
	assemblerInfo->peakMemoryUsageForSummaryStats() <<
        "</table>"
        ;
}



void Assembler::writeAssemblySummaryJson(ostream& json)
{
    AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    using std::setprecision;



    // Compute the number of run-length k-mers used as markers.
    uint64_t totalRleKmerCount = 0;
    uint64_t markerRleKmerCount = 0;
    for(const auto& tableEntry: kmerTable) {
        if(tableEntry.isRleKmer) {
            ++totalRleKmerCount;
            if(tableEntry.isMarker) {
                ++markerRleKmerCount;
            }
        }
    }

    const uint64_t totalDiscardedReadCount =
        assemblerInfo->discardedInvalidBaseReadCount +
        assemblerInfo->discardedShortReadReadCount +
        assemblerInfo->discardedBadRepeatCountReadCount;
    const uint64_t totalDiscardedBaseCount =
        assemblerInfo->discardedInvalidBaseBaseCount +
        assemblerInfo->discardedShortReadBaseCount +
        assemblerInfo->discardedBadRepeatCountBaseCount;


    json <<
        "{\n"
        "  \"Comment\": \"See AssemblySummary.html for a human-readable version of this file\",\n"



        "  \"Shasta version\": \"" << buildId() << "\",\n"



        "  \"Reads used in this assembly\":\n"
        "  {\n"
        "    \"Number of reads\": " << assemblerInfo->readCount << ",\n"
        "    \"Number of raw sequence bases\": " << assemblerInfo->baseCount << ",\n"
        "    \"Average read length (for raw read sequence)\": " <<
        assemblerInfo->baseCount / assemblerInfo->readCount << ",\n"
        "    \"Read N50 (for raw read sequence)\": " << assemblerInfo->readN50 << ",\n"
        "    \"Number of run-length encoded bases\": " << readRepeatCounts.totalSize() << ",\n"
        "    \"Average length ratio of run-length encoded sequence over raw sequence\": " <<
        setprecision(4) << double(readRepeatCounts.totalSize()) / double(assemblerInfo->baseCount) << ",\n"
        "    \"Number of reads flagged as palindromic\": " << assemblerInfo->palindromicReadCount << ",\n"
        "    \"Number of reads flagged as chimeric\": " << assemblerInfo->chimericReadCount << "\n"
        "  },\n"



        "  \"Reads discarded on input\":\n"
        "  {\n"
        "    \"Reads discarded on input because they contained invalid bases\":\n"
        "    {\n"
        "      \"Reads\": " << assemblerInfo->discardedInvalidBaseReadCount << ",\n"
        "      \"Bases\": " << assemblerInfo->discardedInvalidBaseBaseCount << "\n"
        "    },\n"
        "    \"Reads discarded on input because they were too short\":\n"
        "    {\n"
        "      \"Reads\": " << assemblerInfo->discardedShortReadReadCount << ",\n"
        "      \"Bases\": " << assemblerInfo->discardedShortReadBaseCount << "\n"
        "    },\n"
        "    \"Reads discarded on input because they contained repeat counts greater than 255\":\n"
        "    {\n"
        "      \"Reads\": " << assemblerInfo->discardedBadRepeatCountReadCount << ",\n"
        "      \"Bases\": " << assemblerInfo->discardedBadRepeatCountBaseCount << "\n"
        "    },\n"
        "    \"Reads discarded on input, total\":\n"
        "    {\n"
        "      \"Reads\": " << totalDiscardedReadCount << ",\n"
        "      \"Bases\": " << totalDiscardedBaseCount << "\n"
        "    },\n"
        "    \"Fraction of reads discarded on input over total present in input files\":\n"
        "    {\n"
        "      \"Reads\": " <<
        double(totalDiscardedReadCount) /
        double(totalDiscardedReadCount + assemblerInfo->readCount)
        << ",\n"
        "      \"Bases\": " <<
        double(totalDiscardedBaseCount) /
        double(totalDiscardedBaseCount + assemblerInfo->baseCount)
        << "\n"
        "    }\n"
        "  },\n"



        "  \"Marker k-mers\":\n"
        "  {\n"
        "    \"Length k of k-mers used as markers\": " << assemblerInfo->k << ",\n"
        "    \"Total number of k-mers\": " << totalRleKmerCount << ",\n"
        "    \"Number of k-mers used as markers\": " << markerRleKmerCount << ",\n"
        "    \"Fraction of k<-mers used as markers\": " <<
        setprecision(3) << double(markerRleKmerCount) / double(totalRleKmerCount) <<
        "  },\n"



        "  \"Markers\":\n"
        "  {\n"
        "    \"Total number of markers on all reads, one strand\": "
        << markers.totalSize()/2 << ",\n"
        "    \"Total number of markers on all reads, both strands\": "
        << markers.totalSize() << ",\n"
        "    \"Average number of markers per raw base\": "
        << setprecision(4) << double(markers.totalSize()/2)/double(assemblerInfo->baseCount) << ",\n"
        "    \"Average number of markers per run-length encoded base\": "
        << setprecision(4) << double(markers.totalSize()/2)/double(readRepeatCounts.totalSize()) << ",\n"
        "    \"Average base offset between markers in raw sequence\": "
        << setprecision(4) << double(assemblerInfo->baseCount)/double(markers.totalSize()/2) << ",\n"
        "    \"Average base offset between markers in run-length encoded sequence\": "
        << setprecision(4) << double(readRepeatCounts.totalSize())/double(markers.totalSize()/2) << ",\n"
        "    \"Average base gap between markers in run-length encoded sequence\": "
        << setprecision(4) <<
        double(readRepeatCounts.totalSize())/double(markers.totalSize()/2) - double(assemblerInfo->k) << "\n"
        "  },\n"


        "  \"Alignments\":\n"
        "  {\n"
        "    \"Number of alignment candidates found by the LowHash algorithm\": " <<
        alignmentCandidates.candidates.size() << ",\n"
        "    \"Number of good alignments\": " << alignmentData.size() << ",\n"
        "    \"Number of good alignments kept in the read graph\": " << (directedReadGraph.isOpen() ?
            directedReadGraph.edges.size()/2 :
            readGraph.edges.size()/2) << "\n"
        "  },\n"



        "  \"Read graph\":\n"
        "  {\n"
        "    \"Number of vertices\": " << (directedReadGraph.isOpen() ?
            directedReadGraph.vertices.size() :
            readGraph.connectivity.size()) << ",\n"
        "    \"Number of edges\": " << (directedReadGraph.isOpen() ?
            directedReadGraph.edges.size() :
            readGraph.edges.size()) << ",\n"
        "    \"Isolated reads\":\n"
        "    {\n"
        "      \"Reads\": " << assemblerInfo->isolatedReadCount << ",\n"
        "      \"Bases\": " << assemblerInfo->isolatedReadBaseCount << "\n"
        "    },\n"
        "    \"Non-isolated reads\":\n"
        "    {\n"
        "      \"Reads\": " << assemblerInfo->readCount - assemblerInfo->isolatedReadCount << ",\n"
        "      \"Bases\": " << assemblerInfo->baseCount - assemblerInfo->isolatedReadBaseCount << "\n"
        "    },\n"
        "    \"Isolated reads fraction\":\n"
        "    {\n"
        "      \"Reads\": " << double(assemblerInfo->isolatedReadCount)/double(assemblerInfo->readCount) << ",\n"
        "      \"Bases\": " << double(assemblerInfo->isolatedReadBaseCount)/double(assemblerInfo->baseCount) << "\n"
        "    },\n"
        "    \"Non-isolated reads fraction\":\n"
        "    {\n"
        "      \"Reads\": " << double(assemblerInfo->readCount-assemblerInfo->isolatedReadCount)/double(assemblerInfo->readCount) << ",\n"
        "      \"Bases\": " << double(assemblerInfo->baseCount-assemblerInfo->isolatedReadBaseCount)/double(assemblerInfo->baseCount) << "\n"
        "    }\n"
        "  },\n"



        "  \"Marker graph\":\n"
        "  {\n"
        "    \"Total number of vertices\": " << markerGraph.vertexCount() << ",\n"
        "    \"Total number of edges\": " << markerGraph.edges.size() << ",\n"
        "    \"Number of vertices that are not isolated after edge removal\": " <<
        assemblerInfo->markerGraphVerticesNotIsolatedCount << ",\n"
        "    \"Number of edges that were not removed\": " << assemblerInfo->markerGraphEdgesNotRemovedCount << "\n"
        "  },\n"



        "  \"Assembly graph\":\n"
        "  {\n"
        "    \"Number of vertices\": " << assemblyGraph.vertices.size() << ",\n"
        "    \"Number of edges\": " << assemblyGraph.edges.size() << ",\n"
        "    \"Number of edges assembled\": " << assemblerInfo->assemblyGraphAssembledEdgeCount << "\n"
        "  },\n"



        "  \"Assembled segments\":\n"
        "  {\n"
        "    \"Number of segments assembled\": " << assemblerInfo->assemblyGraphAssembledEdgeCount << ",\n"
        "    \"Total assembled segment length\": " << assemblerInfo->totalAssembledSegmentLength << ",\n"
        "    \"Longest assembled segment length\": " << assemblerInfo->longestAssembledSegmentLength << ",\n"
        "    \"Assembled segments N50\": " << assemblerInfo->assembledSegmentN50 << "\n"
        "  },\n"



        "  \"Performance\":\n"
        "  {\n"
        "    \"Elapsed time (seconds)\": " << assemblerInfo->assemblyElapsedTimeSeconds << ",\n"
        "    \"Elapsed time (minutes)\": " << assemblerInfo->assemblyElapsedTimeSeconds/60. << ",\n"
        "    \"Elapsed time (hours)\": " << assemblerInfo->assemblyElapsedTimeSeconds/3600. << ",\n"
        "    \"Average CPU utilization\": " << assemblerInfo->averageCpuUtilization << ",\n"
        "    \"Peak Memory utilization (bytes)\": " <<
	assemblerInfo->peakMemoryUsageForSummaryStats() <<
	"\n"
        "  }\n"



        "}";
}



void Assembler::writeAssemblyIndex(ostream& html) const
{
    writeHtmlBegin(html, false);

    const string s = R"ABCDE(
    <body>
    <h1>Assembly output files</h1>
    <table>

    <tr>
    <td><a href='AssemblySummary.html'>AssemblySummary.html</a>
    <td>Assembly summary information.

    <tr>
    <td><a href='Assembly.fasta'>Assembly.fasta</a>
    <td>Assembly in Fasta format (one strand only).

    <tr>
    <td><a href='Assembly.gfa'>Assembly.gfa</a>
    <td>Assembly in gfa format (one strand only).

    <tr>
    <td><a href='Assembly-BothStrands.gfa'>Assembly-BothStrands.gfa</a>
    <td>Assembly in gfa format (both strands).

    <tr>
    <td><a href='AssemblySummary.csv'>AssemblySummary.csv</a>
    <td>List of assembled segments in order of decreasing length.

    <tr>
    <td><a href='Binned-ReadLengthHistogram.csv'>Binned-ReadLengthHistogram.csv</a>
    <td>Read length distribution in 1 kb bins.

    <tr>
    <td><a href='ReadLengthHistogram.csv'>ReadLengthHistogram.csv</a>
    <td>Detailed read length distribution.

    <tr>
    <td><a href='ReadSummary.csv'>ReadSummary.csv</a>
    <td>Summary file containing one line of information for each read.

    <tr>
    <td><a href='LowHashBucketHistogram.csv'>LowHashBucketHistogram.csv</a>
    <td>MinHash bucket population histogram.

    <tr>
    <td><a href='DisjointSetsHistogram.csv'>DisjointSetsHistogram.csv</a>
    <td>Coverage histogram for all disjoint sets.

    <tr>
    <td><a href='MarkerGraphVertexCoverageHistogram.csv'>MarkerGraphVertexCoverageHistogram.csv</a>
    <td>Coverage histogram for marker graph vertices.

    <tr>
    <td><a href='MarkerGraphEdgeCoverageHistogram.csv'>MarkerGraphEdgeCoverageHistogram.csv</a>
    <td>Coverage histogram for marker graph edges.

    <tr>
    <td><a href='ReadLengthHistogram.csv'>ReadLengthHistogram.csv</a>
    <td>Detailed read length distribution.

    <tr>
    <td><a href='SuppressedAlignmentCandidates.csv'>SuppressedAlignmentCandidates.csv</a>
    <td>Details of suppressed alignment candidates.

    </table>
    </body>
)ABCDE";

    html << s;
    writeHtmlEnd(html);
}



#ifdef SHASTA_HTTP_SERVER


void Assembler::blastRead(
    const vector<string>& request,
    ostream& html)
{

    if(!filesystem::isRegularFile(httpServerData.referenceFastaFileName)) {
        html << "<p>The fasta sequence " << httpServerData.referenceFastaFileName <<
            " to be used as the reference (Blast subject) does not exist.";
        return;
    }



    // Get the ReadId and Strand from the request.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);
    if(!(readIdIsPresent && strandIsPresent)) {
        return;
    }

    // Get the begin and end position.
    uint32_t beginPosition = 0;
    getParameterValue(request, "beginPosition", beginPosition);
    uint32_t endPosition = 0;
    const bool endPositionIsPresent = getParameterValue(request, "endPosition", endPosition);

    // Get blast options.
    string blastOptions;
    getParameterValue(request, "blastOptions", blastOptions);
    string summaryString;
    const bool isSummary = getParameterValue(request, "summary", summaryString);
    if(isSummary) {
        blastOptions =
            "-outfmt '10 bitscore qstart qend sseqid sstart send length pident' "
            "-evalue 1e-200";
            // The following is used to avoid breaking up alignments too much.
            // But it also slows down the search a lot.
            // "-reward 3 -penalty -2 -gapopen 5 -gapextend 5";
    }



    // Access the read.
    if(readId >= reads.size()) {
        html << "<p>Invalid read id.";
        return;
    }
    if(strand!=0 && strand!=1) {
        html << "<p>Invalid strand.";
        return;
    }
    const OrientedReadId orientedReadId(readId, strand);
    const vector<Base> rawOrientedReadSequence = getOrientedReadRawSequence(orientedReadId);
    if(!endPositionIsPresent) {
        endPosition = uint32_t(rawOrientedReadSequence.size());
    }
    if(endPosition <= beginPosition) {
        html << "<p>Invalid choice of begin and end position.";
        return;
    }



    // Write a title.
    html << "<h1>Blast results for oriented read " << orientedReadId;
    html << ", position range " << beginPosition << " " << endPosition;
    html << " (" << endPosition-beginPosition << " bases)</h1>";



    // Create a fasta file with this sequence.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string fastaFileName = tmpDirectory() + uuid + ".fa";
    ofstream fastaFile(fastaFileName);
    fastaFile << ">" << OrientedReadId(readId, strand);
    fastaFile << "-" << beginPosition << "-" << endPosition<< "\n";
    copy(rawOrientedReadSequence.begin() + beginPosition,
        rawOrientedReadSequence.begin() + endPosition,
        ostream_iterator<Base>(fastaFile));
    fastaFile << "\n";
    fastaFile.close();



    // Create the blast command and run it.
    const string blastOutputFileName = tmpDirectory() + uuid + ".txt";
    const string blastErrFileName = tmpDirectory() + uuid + ".errtxt";
    const string command = "blastn -task megablast -subject " + httpServerData.referenceFastaFileName +
        " -query " + fastaFileName + " 1>" + blastOutputFileName + " 2>" + blastErrFileName +
        " " + blastOptions;
    ::system(command.c_str());



    // Copy any error output to html.
    if(filesystem::fileSize(blastErrFileName)) {
        ifstream blastErrFile(blastErrFileName);
        html << "<pre style='font-size:10px'>";
        html << blastErrFile.rdbuf();
        html << "</pre>";
        blastErrFile.close();
    }



    // Output to html.
    if(isSummary) {

        html << "<br>Blast options used: " << blastOptions << "<p>";

        // Tokenize and gather the output, each line with its score.
        using Separator = boost::char_separator<char>;
        using Tokenizer = boost::tokenizer<Separator>;
        const Separator separator(",");
        vector< pair<double, vector<string> > > alignments;
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

            // Extract the score.
            SHASTA_ASSERT(!tokens.empty());
            const double score = std::stod(tokens.front());;

            // Store it.
            alignments.push_back(make_pair(score, tokens));
        }

        // Sort by score.
        sort(alignments.begin(), alignments.end(),
            std::greater< pair<double, vector<string> > >());

        // Write it out.
        html <<
            "<table><tr>"
            "<th rowspan=2>Bit<br>score"
            "<th colspan=3>In " << orientedReadId <<
            "<th colspan=5>In " << httpServerData.referenceFastaFileName <<
            "<th rowspan=2>Alignment<br>length"
            "<th rowspan=2>Identity<br>(%)"
            "<tr>"
            "<th>Begin"
            "<th>End"
            "<th>Length"
            "<th>Strand"
            "<th>Name"
            "<th>Begin"
            "<th>End"
            "<th>Length";
        for(const auto& p: alignments) {
            const auto& tokens = p.second;
            // bitscore qstart qend sseqid sstart send length pident
            SHASTA_ASSERT(tokens.size() == 8);
            const string& bitscore = tokens[0];
            const size_t qstart = std::stoi(tokens[1]) + beginPosition;
            const size_t qend = std::stoi(tokens[2]) + beginPosition;
            const string& sseqid = tokens[3];
            size_t sstart = std::stoi(tokens[4]);
            size_t send = std::stoi(tokens[5]);
            const string& length = tokens[6];
            const string& pident = tokens[7];
            Strand strand = 0;
            if(send < sstart) {
                swap(sstart, send);
                strand = 1;
            }
            html <<
                "<tr style='text-align:center'>"
                "<td>" << bitscore <<
                "<td>" << qstart <<
                "<td>" << qend <<
                "<td>" << qend-qstart <<
                "<td>" << (strand==0 ? "+" : "-") << " (" << strand << ")"
                "<td>" << sseqid <<
                "<td>" << sstart <<
                "<td>" << send <<
                "<td>" << send-sstart <<
                "<td>" << length <<
                "<td>" << pident;
        }
        html << "</table>";

    } else {

        // This is not summary output.
        // Just copy Blast output to html.
        ifstream blastOutputFile(blastOutputFileName);
        html << "<pre style='font-size:10px'>";
        html << blastOutputFile.rdbuf();
        html << "</pre>";
    }



    // Remove the files we created.
    filesystem::remove(fastaFileName);
    filesystem::remove(blastOutputFileName);
    filesystem::remove(blastErrFileName);
}



void shasta::writeStrandSelection(
    ostream& html,          // The html stream to write the form to.
    const string& name,     // The selection name.
    bool select0,           // Whether strand 0 is selected.
    bool select1)           // Whether strand 1 is selected.
{
    html <<
        "<select name=" << name <<
        " title='Choose 0 (+) for the input read or 1 (-) for its reverse complement'>"

        "<option value=0"
        << (select0 ? " selected" : "") <<
        ">0 (+)</option>"

        "<option value=1"
        << (select1 ? " selected" : "") <<
        ">1 (-)</option>"

        "</select>";

}

#endif
