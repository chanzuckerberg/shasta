
// Shasta.
// PngImage.hpp must be included first because of png issues on Ubuntu 16.04.
#include "PngImage.hpp"
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "AlignmentGraph.hpp"
#include "buildId.hpp"
#include "deduplicate.hpp"
#include "LocalAlignmentGraph.hpp"
#include "platformDependent.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/tokenizer.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Seqan
#ifdef SHASTA_HTTP_SERVER
#include <seqan/align.h>
#endif



// Standard library.
#include "chrono.hpp"
#include <iomanip>
#include "iterator.hpp"

#ifdef SHASTA_HTTP_SERVER


#define SHASTA_ADD_TO_FUNCTION_TABLE(name) httpServerData.functionTable[string("/") + #name ] = &Assembler::name



// Associate http keywords wth member functions.
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



void Assembler::exploreAlignments(
    const vector<string>& request,
    ostream& html)
{
    // Get the ReadId and Strand from the request.
    ReadId readId0 = 0;
    const bool readId0IsPresent = getParameterValue(request, "readId", readId0);
    Strand strand0 = 0;
    const bool strand0IsPresent = getParameterValue(request, "strand", strand0);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Show alignments involving read'> "
        "<input type=text name=readId required" <<
        (readId0IsPresent ? (" value=" + to_string(readId0)) : "") <<
        " size=8 title='Enter a read id between 0 and " << reads.size()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand", strand0IsPresent && strand0==0, strand0IsPresent && strand0==1);
    html << "</form>";

    // If the readId or strand are missing, stop here.
    if(!readId0IsPresent || !strand0IsPresent) {
        return;
    }

    // Page title.
    const OrientedReadId orientedReadId0(readId0, strand0);
    html <<
        "<h1>Alignments involving oriented read "
        "<a href='exploreRead?readId=" << readId0  << "&strand=" << strand0 << "'>"
        << OrientedReadId(readId0, strand0) << "</a>"
        << " (" << markers[orientedReadId0.getValue()].size() << " markers)"
        "</h1>";


#if 0
    // Begin the table.
    html <<
        "<table><tr>"
        "<th rowspan=2>Other<br>oriented<br>read"
        "<th rowspan=2 title='The number of aligned markers. Click on a cell in this column to see more alignment details.'>Aligned<br>markers"
        "<th colspan=5>Markers on oriented read " << OrientedReadId(readId0, strand0) <<
        "<th colspan=5>Markers on other oriented read"
        "<tr>";
    for(int i=0; i<2; i++) {
        html <<
            "<th title='Number of aligned markers on the left of the alignment'>Left<br>unaligned"
            "<th title='Number of markers in the aligned range'>Alignment<br>range"
            "<th title='Number of aligned markers on the right of the alignment'>Right<br>unaligned"
            "<th title='Total number of markers on the oriented read'>Total"
            "<th title='Fraction of aligned markers in the alignment range'>Aligned<br>fraction";
    }
#endif


    // Loop over the alignments that this oriented read is involved in, with the proper orientation.
    const vector< pair<OrientedReadId, AlignmentInfo> > alignments =
        findOrientedAlignments(orientedReadId0);
    if(alignments.empty()) {
        html << "<p>No alignments found.";
    } else {
        html << "<p>Found " << alignments.size() << " alignments.";
        displayAlignments(orientedReadId0, alignments, html);
    }

}


void Assembler::displayAlignment(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const AlignmentInfo& alignment,
    ostream& html) const
{
    vector< pair<OrientedReadId, AlignmentInfo> > alignments;
    alignments.push_back(make_pair(orientedReadId1, alignment));
    displayAlignments(orientedReadId0, alignments, html);
}



// Display alignments in an html table.
void Assembler::displayAlignments(
    OrientedReadId orientedReadId0,
    const vector< pair<OrientedReadId, AlignmentInfo> >& alignments,
    ostream& html) const
{
    const ReadId readId0 = orientedReadId0.getReadId();
    const Strand strand0 = orientedReadId0.getStrand();
    const uint32_t markerCount0 = uint32_t(markers[orientedReadId0.getValue()].size());


    // Compute the maximum number of markers that orientedReadId1
    // hangs out of orientedReadId0 on the left and right.
    uint32_t maxLeftHang = 0;
    uint32_t maxRightHang = 0;
    for(size_t i=0; i<alignments.size(); i++) {
        const auto& p = alignments[i];

        // Access information for this alignment.
        const AlignmentInfo& alignmentInfo = p.second;
        const uint32_t leftTrim0  = alignmentInfo.data[0].leftTrim ();
        const uint32_t leftTrim1  = alignmentInfo.data[1].leftTrim ();
        const uint32_t rightTrim0 = alignmentInfo.data[0].rightTrim();
        const uint32_t rightTrim1 = alignmentInfo.data[1].rightTrim();

        // Update the maximum left hang.
        if(leftTrim1 > leftTrim0) {
            maxLeftHang = max(maxLeftHang, leftTrim1 - leftTrim0);
        }

        // Update the maximum left hang.
        if(rightTrim1 > rightTrim0) {
            maxRightHang = max(maxRightHang, rightTrim1 - rightTrim0);
        }
    }



    // Buttons to scale the alignment sketches.
    html <<
        "<script>"
        "function scale(factor)"
        "{"
        "    var elements = document.getElementsByClassName('sketch');"
        "    for (i=0; i<elements.length; i++) {"
        "        elements[i].style.width = factor * parseFloat(elements[i].style.width) + 'px'"
        "    }"
        "}"
        "function larger() {scale(1.5);}"
        "function smaller() {scale(1./1.5);}"
        "</script>";
    if(alignments.size() > 1) {
        html <<
        "&nbsp;<button onclick='larger()'>Make alignment sketches larger</button>"
        "&nbsp;<button onclick='smaller()'>Make alignment sketches smaller</button>"
        ;
    } else {
        html <<
        "&nbsp;<button onclick='larger()'>Make alignment sketch larger</button>"
        "&nbsp;<button onclick='smaller()'>Make alignment sketch smaller</button>"
        ;
    }

    // Begin the table.
    const double markersPerPixel = 50.; // Controls the scaling of the alignment sketch.
    html <<
        "<p><table>"
        "<tr>"
        "<th rowspan=2>Index"
        "<th rowspan=2>Other<br>oriented<br>read"
        "<th rowspan=2 title='The number of aligned markers. Click on a cell in this column to see more alignment details.'>Aligned<br>markers"
        "<th colspan=3>Ordinal offset"
        "<th rowspan=2 title='The marker offset of the centers of the two oriented reads.'>Center<br>offset"
        "<th colspan=5>Markers on oriented read " << orientedReadId0;
    if(alignments.size() > 1) {
        html << "<th colspan=5>Markers on other oriented read";
    } else {
        html << "<th colspan=5>Markers on oriented read " <<
            alignments.front().first;
    }
    html <<
        "<th rowspan=2>Alignment sketch"
        "<tr>"
        "<th>Min"
        "<th>Ave"
        "<th>Max";
    for(int i=0; i<2; i++) {
        html <<
            "<th title='Number of aligned markers on the left of the alignment'>Left<br>unaligned"
            "<th title='Number of markers in the aligned range'>Alignment<br>range"
            "<th title='Number of aligned markers on the right of the alignment'>Right<br>unaligned"
            "<th title='Total number of markers on the oriented read'>Total"
            "<th title='Fraction of aligned markers in the alignment range'>Aligned<br>fraction";
    }



    // Loop over the alignments.
    for(size_t i=0; i<alignments.size(); i++) {
        const auto& p = alignments[i];

        // Access information for this alignment.
        const OrientedReadId orientedReadId1 = p.first;
        const AlignmentInfo& alignmentInfo = p.second;
        const ReadId readId1 = orientedReadId1.getReadId();
        const ReadId strand1 = orientedReadId1.getStrand();
        const uint32_t markerCount1 = uint32_t(markers[orientedReadId1.getValue()].size());

        const uint32_t leftTrim0 = alignmentInfo.data[0].leftTrim();
        const uint32_t leftTrim1 = alignmentInfo.data[1].leftTrim();
        const uint32_t rightTrim0 = alignmentInfo.data[0].rightTrim();
        const uint32_t rightTrim1 = alignmentInfo.data[1].rightTrim();

        // Write a row in the table for this alignment.
        html <<
            "<tr>"
            "<td class=centered>" << i <<
            "<td class=centered><a href='exploreRead?readId=" << readId1  << "&strand=" << strand1 <<
            "' title='Click to see this read'>" << orientedReadId1 << "</a>"
            "<td class=centered>"
            "<a href='exploreAlignment"
            "?readId0=" << readId0 << "&strand0=" << strand0 <<
            "&readId1=" << readId1 << "&strand1=" << strand1 <<
            "' title='Click to see the alignment'>" << alignmentInfo.markerCount << "</a>"
            "<td>" << alignmentInfo.minOrdinalOffset <<
            "<td>" << alignmentInfo.maxOrdinalOffset <<
            "<td>" << alignmentInfo.averageOrdinalOffset <<
            "<td class=centered>" << std::setprecision(6) << alignmentInfo.offsetAtCenter() <<
            "<td class=centered>" << alignmentInfo.leftTrim(0) <<
            "<td class=centered>" << alignmentInfo.range(0) <<
            "<td class=centered>" << alignmentInfo.rightTrim(0) <<
            "<td class=centered>" << markerCount0 <<
            "<td class=centered>" << std::setprecision(2) <<
            alignmentInfo.alignedFraction(0) <<
            "<td class=centered>" << alignmentInfo.leftTrim(1) <<
            "<td class=centered>" << alignmentInfo.range(1) <<
            "<td class=centered>" << alignmentInfo.rightTrim(1) <<
            "<td class=centered>" << markerCount1 <<
            "<td class=centered>" << std::setprecision(2) <<
            alignmentInfo.alignedFraction(1);



        // Write the alignment sketch.
        html <<
            "<td class=centered style='line-height:8px;white-space:nowrap'>"

            // Oriented read 0.
            "<div class=sketch style='display:inline-block;margin:0px;padding:0px;"
            "background-color:white;height:6px;width:" << double(maxLeftHang)/markersPerPixel <<
            "px;'></div>"
            "<div class=sketch title='Oriented read " << orientedReadId0 <<
            "' style='display:inline-block;margin:0px;padding:0px;"
            "background-color:blue;height:6px;width:" << double(markerCount0)/markersPerPixel <<
            "px;'></div>"
            "<div class=sketch style='display:inline-block;margin:0px;padding:0px;"
            "background-color:white;height:6px;width:" << double(maxRightHang)/markersPerPixel <<
            "px;'></div>"

            // Aligned portion.
            "<br>"
            "<div class=sketch style='display:inline-block;margin:0px;padding:0px;"
            "background-color:white;height:6px;width:" << double(maxLeftHang+leftTrim0)/markersPerPixel <<
            "px;'></div>"
            "<div class=sketch title='Aligned portion'"
            " style='display:inline-block;margin:0px;padding:0px;"
            "background-color:red;height:6px;width:" << double(markerCount0-leftTrim0-rightTrim0)/markersPerPixel <<
            "px;'></div>"
            "<div class=sketch style='display:inline-block;margin:0px;padding:0px;"
            "background-color:white;height:6px;width:" << double(maxRightHang+rightTrim0)/markersPerPixel <<
            "px;'></div>"

            // Oriented read 1.
            "<br>"
            "<div class=sketch style='display:inline-block;margin:0px;padding:0px;"
            "background-color:white;height:6px;width:" << double(maxLeftHang+leftTrim0-leftTrim1)/markersPerPixel <<
            "px;'></div>"
            "<div class=sketch title='Oriented read " << orientedReadId1 <<
            "' style='display:inline-block;margin:0px;padding:0px;"
            "background-color:green;height:6px;width:" << double(markerCount1)/markersPerPixel <<
            "px;'></div>"
            "<div class=sketch style='display:inline-block;margin:0px;padding:0px;"
            "background-color:white;height:6px;width:" << double(maxRightHang+rightTrim0-rightTrim1)/markersPerPixel <<
            "px;'></div>"
             ;
    }

    html << "</table>";
}



void Assembler::exploreAlignment(
    const vector<string>& request,
    ostream& html)
{
    // Get the read ids and strands from the request.
    ReadId readId0 = 0;
    const bool readId0IsPresent = getParameterValue(request, "readId0", readId0);
    Strand strand0 = 0;
    const bool strand0IsPresent = getParameterValue(request, "strand0", strand0);
    ReadId readId1 = 0;
    const bool readId1IsPresent = getParameterValue(request, "readId1", readId1);
    Strand strand1 = 0;
    const bool strand1IsPresent = getParameterValue(request, "strand1", strand1);

    // Get alignment parameters.
    int method = httpServerData.assemblerOptions->alignOptions.alignMethod;
    getParameterValue(request, "method", method);
    size_t maxSkip = httpServerData.assemblerOptions->alignOptions.maxSkip;
    getParameterValue(request, "maxSkip", maxSkip);
    size_t maxDrift = httpServerData.assemblerOptions->alignOptions.maxDrift;
    getParameterValue(request, "maxDrift", maxDrift);
    uint32_t maxMarkerFrequency = httpServerData.assemblerOptions->alignOptions.maxMarkerFrequency;
    getParameterValue(request, "maxMarkerFrequency", maxMarkerFrequency);
    int matchScore = httpServerData.assemblerOptions->alignOptions.matchScore;
    getParameterValue(request, "matchScore", matchScore);
    
    uint32_t minAlignedMarkerCount = httpServerData.assemblerOptions->alignOptions.minAlignedMarkerCount;
    getParameterValue(request, "minAlignedMarkerCount", minAlignedMarkerCount);
    double minAlignedFraction = httpServerData.assemblerOptions->alignOptions.minAlignedFraction;
    getParameterValue(request, "minAlignedFraction", minAlignedFraction);
    uint32_t maxTrim = httpServerData.assemblerOptions->alignOptions.maxTrim;
    getParameterValue(request, "maxTrim", maxTrim);
    
    int mismatchScore = httpServerData.assemblerOptions->alignOptions.mismatchScore;
    getParameterValue(request, "mismatchScore", mismatchScore);
    int gapScore = httpServerData.assemblerOptions->alignOptions.gapScore;
    getParameterValue(request, "gapScore", gapScore);
    double downsamplingFactor = httpServerData.assemblerOptions->alignOptions.downsamplingFactor;
    getParameterValue(request, "downsamplingFactor", downsamplingFactor);
    int bandExtend = httpServerData.assemblerOptions->alignOptions.bandExtend;
    getParameterValue(request, "bandExtend", bandExtend);



    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Compute marker alignment'>"
        "&nbsp of read &nbsp"
        "<input type=text name=readId0 required size=8 " <<
        (readId0IsPresent ? "value="+to_string(readId0) : "") <<
        " title='Enter a read id between 0 and " << reads.size()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand0", strand0IsPresent && strand0==0, strand0IsPresent && strand0==1);
    html <<
        "&nbsp and read <input type=text name=readId1 required size=8 " <<
        (readId1IsPresent ? "value="+to_string(readId1) : "") <<
        " title='Enter a read id between 0 and " << reads.size()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand1", strand1IsPresent && strand1==0, strand1IsPresent && strand1==1);

    renderEditableAlignmentConfig(
        method,
        maxSkip,
        maxDrift,
        maxMarkerFrequency,
        minAlignedMarkerCount,
        minAlignedFraction,
        maxTrim,
        matchScore,
        mismatchScore,
        gapScore,
        downsamplingFactor,
        bandExtend,
        html
    );

    html << "</form>";


    // If the readId's or strand's are missing, stop here.
    if(!readId0IsPresent || !strand0IsPresent || !readId1IsPresent || !strand1IsPresent) {
        return;
    }



    // Page title.
    const OrientedReadId orientedReadId0(readId0, strand0);
    const OrientedReadId orientedReadId1(readId1, strand1);
    html <<
        "<h1>Marker alignment of oriented reads " <<
        "<a href='exploreRead?readId=" << readId0 << "&strand=" << strand0 << "'>" << orientedReadId0 << "</a>" <<
        " and " <<
        "<a href='exploreRead?readId=" << readId1 << "&strand=" << strand1 << "'>" << orientedReadId1 << "</a>" <<
        "</h1>";



    // Compute the alignment.
    // This creates file Alignment.png.
    Alignment alignment;
    AlignmentInfo alignmentInfo;
    if(method == 0) {
        array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
        getMarkersSortedByKmerId(orientedReadId0, markersSortedByKmerId[0]);
        getMarkersSortedByKmerId(orientedReadId1, markersSortedByKmerId[1]);
        AlignmentGraph graph;
        const bool debug = true;
        alignOrientedReads(
            markersSortedByKmerId,
            maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);
            
        if(alignment.ordinals.empty()) {
            html << "<p>The alignment is empty (it has no markers).";
            return;
        }
    } else if(method == 1) {
        alignOrientedReads1(
            orientedReadId0, orientedReadId1,
            matchScore, mismatchScore, gapScore, alignment, alignmentInfo);
    } else if(method == 3) {
        alignOrientedReads3(
            orientedReadId0, orientedReadId1,
            matchScore, mismatchScore, gapScore,
            downsamplingFactor, bandExtend,
            alignment, alignmentInfo);
    } else {
        SHASTA_ASSERT(0);
    }


    // Make sure we have Alignment.png to display.
    if(method != 0) {
        vector<MarkerWithOrdinal> sortedMarkers0;
        vector<MarkerWithOrdinal> sortedMarkers1;
        getMarkersSortedByKmerId(orientedReadId0, sortedMarkers0);
        getMarkersSortedByKmerId(orientedReadId1, sortedMarkers1);
        AlignmentGraph::writeImage(
            sortedMarkers0,
            sortedMarkers1,
            alignment,
            "Alignment.png");
    }

    if (alignment.ordinals.empty()) {
        html << "<p>The computed alignment is empty.";
        return;
    }
    if (alignment.ordinals.size() < minAlignedMarkerCount) {
        html << "<p>Alignment has fewer than " << minAlignedMarkerCount << " markers.";
        return;
    }
    if (alignmentInfo.minAlignedFraction() < minAlignedFraction) {
        html << "<p>Min aligned fraction is smaller than " << minAlignedFraction << ".";
        return;
    }

    // If the alignment has too much trim, skip it.
    uint32_t leftTrim;
    uint32_t rightTrim;
    tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();
    if(leftTrim>maxTrim || rightTrim>maxTrim) {
        html << "<p>Alignment has too much trim. Left trim = " << leftTrim
            << " Right trim = " << rightTrim;
        return;
    }

    
    // Write summary information for this alignment.
    html << "<h3>Alignment summary</h3>";
    displayAlignment(
        orientedReadId0,
        orientedReadId1,
        alignmentInfo,
        html);
    html << "<br>See below for alignment details.";


    // Create a base64 version of the png file.
    const string command = "base64 Alignment.png > Alignment.png.base64";
    ::system(command.c_str());


    // Write out the picture with the alignment.
    html <<
        "<h3>Alignment matrix</h3>"
        "<p>In the picture, horizontal positions correspond to marker ordinals on " <<
        orientedReadId0 << " (marker 0 is on left) "
        "and vertical positions correspond to marker ordinals on " <<
        orientedReadId1 << " (marker 0 is on top). "
        "Each faint line corresponds to 10 markers."
        "<p><img id=\"alignmentMatrix\" onmousemove=\"updateTitle(event)\" "
        "src=\"data:image/png;base64,";
    ifstream png("Alignment.png.base64");
    html << png.rdbuf();
    html << "\"/>"
        "<script>"
        "function updateTitle(e)"
        "{"
        "    var element = document.getElementById(\"alignmentMatrix\");"
        "    var rectangle = element.getBoundingClientRect();"
        "    var x = e.clientX - Math.round(rectangle.left);"
        "    var y = e.clientY - Math.round(rectangle.top);"
        "    element.title = " <<
        "\"" << orientedReadId0 << " marker \" + x + \", \" + "
        "\"" << orientedReadId1 << " marker \" + y;"
        "}"
        "</script>";



    // Write out details of the alignment.
    html <<
        "<h3>Alignment details</h3>"
        "<table>"

        "<tr>"
        "<th rowspan=2>K-mer"
        "<th colspan=3>Ordinals"
        "<th colspan=2>Positions<br>(RLE)"

        "<tr>"
        "<th>" << orientedReadId0 <<
        "<th>" << orientedReadId1 <<
        "<th>Offset"
        "<th>" << orientedReadId0 <<
        "<th>" << orientedReadId1;

    const auto markers0 = markers[orientedReadId0.getValue()];
    const auto markers1 = markers[orientedReadId1.getValue()];
    for(const auto& ordinals: alignment.ordinals) {
        const auto ordinal0 = ordinals[0];
        const auto ordinal1 = ordinals[1];
        const auto& marker0 = markers0[ordinal0];
        const auto& marker1 = markers1[ordinal1];
        const auto kmerId = marker0.kmerId;
        SHASTA_ASSERT(marker1.kmerId == kmerId);
        const Kmer kmer(kmerId, assemblerInfo->k);

        html << "<tr><td style='font-family:monospace'>";
        kmer.write(html, assemblerInfo->k);
        html <<

            "<td class=centered>"
            "<a href=\"exploreRead?readId=" << orientedReadId0.getReadId() <<
            "&amp;strand=" << orientedReadId0.getStrand() <<
            "&amp;highlightMarker=" << ordinal0 <<
            "#" << ordinal0 << "\">" << ordinal0 << "</a>"

            "<td class=centered>"
            "<a href=\"exploreRead?readId=" << orientedReadId1.getReadId() <<
            "&amp;strand=" << orientedReadId1.getStrand() <<
            "&amp;highlightMarker=" << ordinal1 <<
            "#" << ordinal1 << "\">" << ordinal1 << "</a>"

            "<td class=centered>" << int32_t(ordinal0) - int32_t(ordinal1) <<
            "<td class=centered>" << marker0.position <<
            "<td class=centered>" << marker1.position;

    }

    html << "</table>";
}



// Display a base-by-base alignment matrix between two given sequences.
void Assembler::displayAlignmentMatrix(
    const vector<string>& request,
    ostream& html)
{
#ifndef __linux__
    html << "<p>This functionality is only available on Linux.";
    return;
#else

    html << "<h1>Base-by-base alignment of two sequences</h1>"
        "<p>This page does not use run-length representation of sequences. "
        "It also does not use markers. "
        "Alignments computed and displayed here are standard "
        "base-by-base alignments.";

    // Get the request parameters.
    string sequenceString0;
    getParameterValue(request, "sequence0", sequenceString0);
    string sequenceString1;
    getParameterValue(request, "sequence1", sequenceString1);
    int zoom = 1;
    getParameterValue(request, "zoom", zoom);
    string clip0String;
    getParameterValue(request, "clip0", clip0String);
    const bool clip0 = (clip0String == "on");
    string clip1String;
    getParameterValue(request, "clip1", clip1String);
    const bool clip1 = (clip1String == "on");
    string showAlignmentString;
    getParameterValue(request, "showAlignment", showAlignmentString);
    const bool showAlignment = (showAlignmentString == "on");
    string showGridString;
    getParameterValue(request, "showGrid", showGridString);
    const bool showGrid = (showGridString == "on");


    // Get the zoom factor.

    // Write the form.
    html <<
        "<p>Display a base-by-base alignment of these two sequences:"
        "<form>"
        "<input style='font-family:monospace' type=text name=sequence0 required size=64 value='" << sequenceString0 << "'>"
        "<br><input style='font-family:monospace' type=text name=sequence1 required size=64 value='" << sequenceString1 << "'>"
        "<br><input type=checkbox name=clip0" << (clip0 ? " checked" : "") << "> Allow clipping on both ends of first sequence."
        "<br><input type=checkbox name=clip1" << (clip1 ? " checked" : "") << "> Allow clipping on both ends of second sequence."
        "<br><input type=checkbox name=showAlignment" << (showAlignment ? " checked" : "") << "> Show the alignment and highlight it in the alignment matrix."
        "<br><input type=checkbox name=showGrid" << (showGrid ? " checked" : "") << "> Show a grid on the alignment matrix."
        "<br>Zoom factor: <input type=text name=zoom required value=" << zoom << ">"
        "<br><input type=submit value='Display'>"
        "</form>";

    // If either sequence is empty, do nothing.
    if(sequenceString0.empty() || sequenceString1.empty()) {
        return;
    }

    // Convert to base sequences, discarding all characters that
    // don't represent a base.
    vector<Base> sequence0;
    for(const char c: sequenceString0) {
        try {
            const Base b = Base::fromCharacter(c);
            sequence0.push_back(b);
        } catch (const std::exception&) {
            // Just discard the character.
        }
    }
    vector<Base> sequence1;
    for(const char c: sequenceString1) {
        try {
            const Base b = Base::fromCharacter(c);
            sequence1.push_back(b);
        } catch (const std::exception&) {
            // Just discard the character.
        }
    }

    // If either sequence is empty, do nothing.
    if(sequenceString0.empty() || sequenceString1.empty()) {
        return;
    }
    if(sequence0.empty() || sequence1.empty()) {
        return;
    }



    // If getting here, we have two non-empty sequences and
    // we can display they alignment matrix.


    // Create the image, which gets initialized to black.
    const int n0 = int(sequence0.size());
    const int n1 = int(sequence1.size());
    PngImage image(n0*zoom, n1*zoom);




    // Display a position grid.
    if(showGrid) {

        // Every 10.
        for(int i0=0; i0<n0; i0+=10) {
            for(int i1=0; i1<n1; i1++) {
                const int begin0 = i0 *zoom;
                const int end0 = begin0 + zoom;
                const int begin1 = i1 *zoom;
                const int end1 = begin1 + zoom;
                for(int j0=begin0; j0!=end0; j0++) {
                    for(int j1=begin1; j1!=end1; j1++) {
                        image.setPixel(j0, j1, 128, 128, 128);
                    }
                }
            }
        }
        for(int i1=0; i1<n1; i1+=10) {
            for(int i0=0; i0<n0; i0++) {
                const int begin0 = i0 *zoom;
                const int end0 = begin0 + zoom;
                const int begin1 = i1 *zoom;
                const int end1 = begin1 + zoom;
                for(int j0=begin0; j0!=end0; j0++) {
                    for(int j1=begin1; j1!=end1; j1++) {
                    	image.setPixel(j0, j1, 128, 128, 128);
                    }
                }
            }
        }

        // Every 100.
        for(int i0=0; i0<n0; i0+=100) {
            for(int i1=0; i1<n1; i1++) {
                const int begin0 = i0 *zoom;
                const int end0 = begin0 + zoom;
                const int begin1 = i1 *zoom;
                const int end1 = begin1 + zoom;
                for(int j0=begin0; j0!=end0; j0++) {
                    for(int j1=begin1; j1!=end1; j1++) {
                    	image.setPixel(j0, j1, 192, 192, 192);
                    }
                }
            }
        }
        for(int i1=0; i1<n1; i1+=100) {
            for(int i0=0; i0<n0; i0++) {
                const int begin0 = i0 *zoom;
                const int end0 = begin0 + zoom;
                const int begin1 = i1 *zoom;
                const int end1 = begin1 + zoom;
                for(int j0=begin0; j0!=end0; j0++) {
                    for(int j1=begin1; j1!=end1; j1++) {
                    	image.setPixel(j0, j1, 192, 192, 192);
                    }
                }
            }
        }
    }



    // Fill in pixel values.
    for(int i0=0; i0<n0; i0++) {
        const Base base0 = sequence0[i0];
        const int begin0 = i0 *zoom;
        const int end0 = begin0 + zoom;
        for(int i1=0; i1<n1; i1++) {
            const Base base1 = sequence1[i1];
            if(!(base1 == base0)) {
                continue;
            }
            const int begin1 = i1 *zoom;
            const int end1 = begin1 + zoom;
            for(int j0=begin0; j0!=end0; j0++) {
                for(int j1=begin1; j1!=end1; j1++) {
                    image.setPixel(j0, j1,0, 255, 0);
                }
            }
        }
    }



    // Use SeqAn to compute an alignment free at both ends
    // and highlight it in the image.
    // https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html
    if(showAlignment){
        using namespace seqan;
        using seqan::Alignment; // Hide shasta::Alignment.

        typedef String<char> TSequence;
        typedef StringSet<TSequence> TStringSet;
        typedef StringSet<TSequence, Dependent<> > TDepStringSet;
        typedef Graph<Alignment<TDepStringSet> > TAlignGraph;
        // typedef Align<TSequence, ArrayGaps> TAlign;

        TSequence seq0;
        for(const Base b: sequence0) {
            appendValue(seq0, b.character());
        }
        TSequence seq1;
        for(const Base b: sequence1) {
            appendValue(seq1, b.character());
        }

        TStringSet sequences;
        appendValue(sequences, seq0);
        appendValue(sequences, seq1);


        // Call the globalAlignment with AlignConfig arguments
        // determined by clip0 and clip1.
        TAlignGraph graph(sequences);
        int score;
        if(clip0) {
            if(clip1) {
                score = globalAlignment(
                    graph,
                    Score<int, Simple>(1, -1, -1),
                    AlignConfig<true, true, true, true>(),
                    LinearGaps());
            } else {
                score = globalAlignment(
                    graph,
                    Score<int, Simple>(1, -1, -1),
                    AlignConfig<true, false, true, false>(),
                    LinearGaps());
            }
        } else {
            if(clip1) {
                score = globalAlignment(
                    graph,
                    Score<int, Simple>(1, -1, -1),
                    AlignConfig<false, true, false, true>(),
                    LinearGaps());
            } else {
                score = globalAlignment(
                    graph,
                    Score<int, Simple>(1, -1, -1),
                    AlignConfig<false, false, false, false>(),
                    LinearGaps());
            }

        }



        // Extract the alignment from the graph.
        // This creates a single sequence consisting of the two rows
        // of the alignment, concatenated.
        TSequence align;
        convertAlignment(graph, align);
        const int totalAlignmentLength = int(seqan::length(align));
        SHASTA_ASSERT((totalAlignmentLength % 2) == 0);    // Because we are aligning two sequences.
        const int alignmentLength = totalAlignmentLength / 2;

        // Extract the two rows of the alignment.
        array<vector<AlignedBase>, 2> alignment;
        alignment[0].resize(alignmentLength);
        alignment[1].resize(alignmentLength);
        for(int i=0; i<alignmentLength; i++) {
            alignment[0][i] = AlignedBase::fromCharacter(align[i]);
            alignment[1][i] = AlignedBase::fromCharacter(align[i + alignmentLength]);
        }
        html << "<br>Sequence lengths: " << n0 << " " << n1 <<
            "<br>Optimal alignment has length " << alignmentLength <<
            ", score " << score <<
            ":<div style='font-family:monospace'>";
        for(size_t i=0; i<2; i++) {
            html << "<br>";
            for(int j=0; j<alignmentLength; j++) {
                html << alignment[i][j];
            }
        }
        html << "</div>";


        int i0 = 0;
        int i1 = 0;
        for(int position=0; position<alignmentLength; position++) {
            const AlignedBase b0 = alignment[0][position];
            const AlignedBase b1 = alignment[1][position];

            if(!(b0.isGap() || b1.isGap())) {

                // This pixel is part of the optimal alignment
                const int begin0 = i0 *zoom;
                const int end0 = begin0 + zoom;
                const int begin1 = i1 *zoom;
                const int end1 = begin1 + zoom;
                for(int j0=begin0; j0!=end0; j0++) {
                    for(int j1=begin1; j1!=end1; j1++) {
                        if(b0 == b1) {
                        	image.setPixel(j0, j1, 255, 0, 0);
                        } else {
                        	image.setPixel(j0, j1, 255, 255, 0);
                        }
                    }
                }
            }

            if(!b0.isGap()) {
                ++i0;
            }
            if(!b1.isGap()) {
                ++i1;
            }
        }

    }



    // Write it out.
    image.write("AlignmentMatrix.png");

    // Create a base64 version of the png file.
    const string command = "base64 AlignmentMatrix.png > AlignmentMatrix.png.base64";
    ::system(command.c_str());


    // Write out the png file.
    html << "<p><img src=\"data:image/png;base64,";
    ifstream png("AlignmentMatrix.png.base64");
    html << png.rdbuf();
    html << "\"/>";

#endif
}

void Assembler::renderEditableAlignmentConfig(
    const int method,
    const uint64_t maxSkip,
    const uint64_t maxDrift,
    const uint32_t maxMarkerFrequency,
    const uint64_t minAlignedMarkerCount,
    const double minAlignedFraction,
    const uint64_t maxTrim,
    const int matchScore,
    const int mismatchScore,
    const int gapScore,
    const double downsamplingFactor,
    const uint32_t bandExtend,
    ostream& html
) {
    const auto& descriptions = httpServerData.assemblerOptions->allOptionsDescription;

    html << "<p><table >";

    html << "<tr><th class=left>[Align]<th class=center>Value<th class=left>Description";
        
    html << "<tr><th class=left>alignMethod<td>"
        "<input type=radio name=method value=0" <<
        (method==0 ? " checked=checked" : "") << "> 0 (Shasta)<br>"
        "<input type=radio name=method value=1" <<
        (method==1 ? " checked=checked" : "") << "> 1 (SeqAn)<br>"
        "<input type=radio name=method value=3" <<
        (method==3 ? " checked=checked" : "") << "> 3 (SeqAn, banded)"
        "<td class=smaller>" << descriptions.find("Align.alignMethod", false).description();

    html << "<tr><th class=left>maxSkip"
        "<td class=centered>"
            "<input type=text style='text-align:center;border:none' name=maxSkip size=16 value=" << maxSkip << ">"
        "<td class=smaller>" << descriptions.find("Align.maxSkip", false).description();

    html << "<tr>"
        "<th class=left>maxDrift"
        "<td class=centered>"
            "<input type=text style='text-align:center;border:none' name=maxDrift size=16 value=" << maxDrift << ">"
        "<td class=smaller>" << descriptions.find("Align.maxDrift", false).description();

    html << "<tr>"
        "<th class=left>maxMarkerFrequency"
        "<td class=centered>"
            "<input type=text style='text-align:center;border:none' name=maxMarkerFrequency size=16 value=" << maxMarkerFrequency << ">"
        "<td class=smaller>" << descriptions.find("Align.maxMarkerFrequency", false).description();

    html << "<tr>"
        "<th class=left>matchScore"
        "<td class=centered>"
            "<input type=text style='text-align:center;border:none' name=matchScore size=16 value=" << matchScore << ">"
        "<td class=smaller>" << descriptions.find("Align.matchScore", false).description();

    html << "<tr>"
        "<th class=left>mismatchScore "
        "<td class=centered>"
            "<input type=text style='text-align:center;border:none' name=mismatchScore size=16 value=" << mismatchScore << ">"
        "<td class=smaller>" << descriptions.find("Align.mismatchScore", false).description();

    html << "<tr>"
        "<th class=left>gapScore"
        "<td class=centered>"
            "<input type=text style='text-align:center;border:none' name=gapScore size=16 value=" << gapScore << ">"
        "<td class=smaller>" << descriptions.find("Align.gapScore", false).description();

    html << "<tr>"
        "<th class=left>downsamplingFactor"
        "<td class=centered>"
            "<input type=text style='text-align:center;border:none' name=downsamplingFactor size=16 value=" << downsamplingFactor << ">"
        "<td class=smaller>" << descriptions.find("Align.downsamplingFactor", false).description();

    html << "<tr>"
        "<th class=left>bandExtend"
        "<td class=centered>"
            "<input type=text style='text-align:center;border:none' name=bandExtend size=16 value=" << bandExtend << ">"
        "<td class=smaller>" << descriptions.find("Align.bandExtend", false).description();

    html << "<tr>"
        "<th class=left>minAlignedMarkers"
        "<td class=centered>"
            "<input type=text style='text-align:center;border:none' name=minAlignedMarkerCount size=16 value=" << minAlignedMarkerCount << ">"
        "<td class=smaller>" << descriptions.find("Align.minAlignedMarkerCount", false).description();

    html << "<tr>"
        "<th class=left>minAlignedFraction"
        "<td class=centered>"
            "<input type=text style='text-align:center;border:none' name=minAlignedFraction size=16 value=" << minAlignedFraction << ">"
        "<td class=smaller>" << descriptions.find("Align.minAlignedFraction", false).description();

    html << "<tr>"
        "<th class=left>maxTrim"
        "<td class=centered>"
            "<input type=text style='text-align:center;border:none' name=maxTrim size=16 value=" << maxTrim << ">"
        "<td class=smaller>" << descriptions.find("Align.maxTrim", false).description();

    html << "</table>";
}

// Compute alignments on an oriented read against
// all other oriented reads.
void Assembler::computeAllAlignments(
    const vector<string>& request,
    ostream& html)
{
    // Get the read id and strand from the request.
    ReadId readId0 = 0;
    const bool readId0IsPresent = getParameterValue(request, "readId0", readId0);
    Strand strand0 = 0;
    const bool strand0IsPresent = getParameterValue(request, "strand0", strand0);

    // Get alignment parameters.
    computeAllAlignmentsData.method = httpServerData.assemblerOptions->alignOptions.alignMethod;
    getParameterValue(request, "method", computeAllAlignmentsData.method);
    computeAllAlignmentsData.minMarkerCount = 0;
    getParameterValue(request, "minMarkerCount", computeAllAlignmentsData.minMarkerCount);
    computeAllAlignmentsData.maxSkip = httpServerData.assemblerOptions->alignOptions.maxSkip;
    getParameterValue(request, "maxSkip", computeAllAlignmentsData.maxSkip);
    computeAllAlignmentsData.maxDrift = httpServerData.assemblerOptions->alignOptions.maxDrift;
    getParameterValue(request, "maxDrift", computeAllAlignmentsData.maxDrift);
    computeAllAlignmentsData.maxMarkerFrequency = httpServerData.assemblerOptions->alignOptions.maxMarkerFrequency;
    getParameterValue(request, "maxMarkerFrequency", computeAllAlignmentsData.maxMarkerFrequency);
    computeAllAlignmentsData.minAlignedMarkerCount = httpServerData.assemblerOptions->alignOptions.minAlignedMarkerCount;
    getParameterValue(request, "minAlignedMarkerCount", computeAllAlignmentsData.minAlignedMarkerCount);
    computeAllAlignmentsData.minAlignedFraction = httpServerData.assemblerOptions->alignOptions.minAlignedFraction;
    getParameterValue(request, "minAlignedFraction", computeAllAlignmentsData.minAlignedFraction);
    computeAllAlignmentsData.maxTrim = httpServerData.assemblerOptions->alignOptions.maxTrim;
    getParameterValue(request, "maxTrim", computeAllAlignmentsData.maxTrim);
    computeAllAlignmentsData.matchScore = httpServerData.assemblerOptions->alignOptions.matchScore;
    getParameterValue(request, "matchScore", computeAllAlignmentsData.matchScore);
    computeAllAlignmentsData.mismatchScore = httpServerData.assemblerOptions->alignOptions.mismatchScore;
    getParameterValue(request, "mismatchScore", computeAllAlignmentsData.mismatchScore);
    computeAllAlignmentsData.gapScore = httpServerData.assemblerOptions->alignOptions.gapScore;
    getParameterValue(request, "gapScore", computeAllAlignmentsData.gapScore);
    computeAllAlignmentsData.downsamplingFactor = httpServerData.assemblerOptions->alignOptions.downsamplingFactor;
    getParameterValue(request, "downsamplingFactor", computeAllAlignmentsData.downsamplingFactor);
    computeAllAlignmentsData.bandExtend = httpServerData.assemblerOptions->alignOptions.bandExtend;
    getParameterValue(request, "bandExtend", computeAllAlignmentsData.bandExtend);


    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Compute marker alignments'>"
        "&nbsp of oriented read &nbsp"
        "<input type=text name=readId0 required size=8 " <<
        (readId0IsPresent ? "value="+to_string(readId0) : "") <<
        " title='Enter a read id between 0 and " << reads.size()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand0", strand0IsPresent && strand0==0, strand0IsPresent && strand0==1);

    renderEditableAlignmentConfig(
        computeAllAlignmentsData.method,
        computeAllAlignmentsData.maxSkip,
        computeAllAlignmentsData.maxDrift,
        computeAllAlignmentsData.maxMarkerFrequency,
        computeAllAlignmentsData.minAlignedMarkerCount,
        computeAllAlignmentsData.minAlignedFraction,
        computeAllAlignmentsData.maxTrim,
        computeAllAlignmentsData.matchScore,
        computeAllAlignmentsData.mismatchScore,
        computeAllAlignmentsData.gapScore,
        computeAllAlignmentsData.downsamplingFactor,
        computeAllAlignmentsData.bandExtend,
        html
    );
    
    html << "</form>";

    
    // If the readId or strand are missing, stop here.
    if(!readId0IsPresent || !strand0IsPresent) {
        return;
    }
 
    const OrientedReadId orientedReadId0(readId0, strand0);

    // Vectors to contain markers sorted by kmerId.
    vector<MarkerWithOrdinal> markers0SortedByKmerId;
    vector<MarkerWithOrdinal> markers1SortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markers0SortedByKmerId);


    // Compute the alignments in parallel.
    computeAllAlignmentsData.orientedReadId0 = orientedReadId0;
    const size_t threadCount =std::thread::hardware_concurrency();
    computeAllAlignmentsData.threadAlignments.resize(threadCount);
    const size_t batchSize = 1000;
    setupLoadBalancing(reads.size(), batchSize);
    const auto t0 = std::chrono::steady_clock::now();
    runThreads(&Assembler::computeAllAlignmentsThreadFunction, threadCount);
    const auto t1 = std::chrono::steady_clock::now();
    html << "<p>Alignment computation using " << threadCount << " threads took " <<
        1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count()) << "s.";

    // Gather the alignments found by each thread.
    vector< pair<OrientedReadId, AlignmentInfo> > alignments;
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        const vector< pair<OrientedReadId, AlignmentInfo> >& threadAlignments =
            computeAllAlignmentsData.threadAlignments[threadId];
        copy(threadAlignments.begin(), threadAlignments.end(), back_inserter(alignments));
    }
    computeAllAlignmentsData.threadAlignments.clear();
    sort(alignments.begin(), alignments.end(),
        OrderPairsByFirstOnly<OrientedReadId, AlignmentInfo>());

#if 0
    // Reusable data structures for alignOrientedReads.
    AlignmentGraph graph;
    Alignment alignment;


    // Loop over oriented reads.
    // Eventually this should be multithreaded if we want to use
    // it for large assemblies.
    size_t computedAlignmentCount = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for(ReadId readId1=0; readId1<reads.size(); readId1++) {
        if((readId1 % 10000) == 0) {
            cout << timestamp << readId1 << "/" << reads.size() << " " << alignments.size() << endl;
        }
        for(Strand strand1=0; strand1<2; strand1++) {

            // Skip alignments with self in the same orientation.
            // Allow alignment with self reverse complemented.
            if(readId0==readId1 && strand0==strand1) {
                continue;
            }

            // If this read has less than the required number of markers, skip.
            const OrientedReadId orientedReadId1(readId1, strand1);
            if(markers[orientedReadId1.getValue()].size() < minMarkerCount) {
                continue;
            }

            // Get markers sorted by kmer id.
            getMarkersSortedByKmerId(orientedReadId1, markers1SortedByKmerId);

            // Compute the alignment.
            ++computedAlignmentCount;
            const bool debug = false;
            alignOrientedReads(
                markers0SortedByKmerId,
                markers1SortedByKmerId,
                maxSkip, maxVertexCountPerKmer, debug, graph, alignment);

            // If the alignment has too few markers skip it.
            if(alignment.ordinals.size() < minAlignedMarkerCount) {
                continue;
            }

            // Compute the AlignmentInfo.
            AlignmentInfo alignmentInfo;
            alignmentInfo.create(alignment);

            // If the alignment has too much trim, skip it.
            uint32_t leftTrim;
            uint32_t rightTrim;
            tie(leftTrim, rightTrim) = computeTrim(
                orientedReadId0,
                orientedReadId1,
                alignmentInfo);
            if(leftTrim>maxTrim || rightTrim>maxTrim) {
                continue;
            }
            alignments.push_back(make_pair(orientedReadId1, alignmentInfo));
        }
    }
    const auto t1 = std::chrono::steady_clock::now();
    html << "<p>Computed " << computedAlignmentCount << " alignments.";
    html << "<p>Found " << alignments.size() <<
        " alignments satisfying the given criteria.";
    html << "<p>Alignment computation took " <<
        1.e-9 * double((std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count()) << "s.";
#endif

    html << "<p>Found " << alignments.size() <<
        " alignments satisfying the given criteria.";
    if(alignments.empty()) {
        html << "<p>No alignments found.";
    } else {
        displayAlignments(orientedReadId0, alignments, html);
    }
}



void Assembler::computeAllAlignmentsThreadFunction(size_t threadId)
{
    // Get the first oriented read.
    const OrientedReadId orientedReadId0 = computeAllAlignmentsData.orientedReadId0;
    const ReadId readId0 = orientedReadId0.getReadId();
    const Strand strand0 = orientedReadId0.getStrand();

    // Get parameters for alignment computation.
    const size_t minMarkerCount = computeAllAlignmentsData.minMarkerCount;
    const size_t maxSkip = computeAllAlignmentsData.maxSkip;
    const size_t maxDrift = computeAllAlignmentsData.maxDrift;
    const uint32_t maxMarkerFrequency = computeAllAlignmentsData.maxMarkerFrequency;
    const size_t minAlignedMarkerCount = computeAllAlignmentsData.minAlignedMarkerCount;
    const double minAlignedFraction = computeAllAlignmentsData.minAlignedFraction;
    const size_t maxTrim = computeAllAlignmentsData.maxTrim;
    const int method = computeAllAlignmentsData.method;
    const int matchScore = computeAllAlignmentsData.matchScore;
    const int mismatchScore = computeAllAlignmentsData.mismatchScore;
    const int gapScore = computeAllAlignmentsData.gapScore;
    const double downsamplingFactor = computeAllAlignmentsData.downsamplingFactor;
    const uint32_t bandExtend = computeAllAlignmentsData.bandExtend;

    // Vector where this thread will store the alignments it finds.
    vector< pair<OrientedReadId, AlignmentInfo> >& alignments =
        computeAllAlignmentsData.threadAlignments[threadId];

    // Reusable data structures for alignOrientedReads.
    AlignmentGraph graph;
    Alignment alignment;
    AlignmentInfo alignmentInfo;

    // Vectors to contain markers sorted by kmerId.
    array<vector<MarkerWithOrdinal>, 2> markersSortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markersSortedByKmerId[0]);

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over reads assigned to this batch.
        for(ReadId readId1=ReadId(begin); readId1!=ReadId(end); ++readId1) {

            // Loop over strands.
            for(Strand strand1=0; strand1<2; strand1++) {

                // Skip alignments with self in the same orientation.
                // Allow alignment with self reverse complemented.
                if(readId0==readId1 && strand0==strand1) {
                    continue;
                }

                // If this read has less than the required number of markers, skip.
                const OrientedReadId orientedReadId1(readId1, strand1);
                if(markers[orientedReadId1.getValue()].size() < minMarkerCount) {
                    continue;
                }

                // Get markers sorted by kmer id.
                getMarkersSortedByKmerId(orientedReadId1, markersSortedByKmerId[1]);

                // Compute the alignment.
                if (method == 0) {
                    const bool debug = false;
                    alignOrientedReads(
                        markersSortedByKmerId,
                        maxSkip, maxDrift, maxMarkerFrequency, debug, graph, alignment, alignmentInfo);
                } else if (method == 1) {
                    alignOrientedReads1(
                        orientedReadId0, orientedReadId1,
                        matchScore, mismatchScore, gapScore, alignment, alignmentInfo);
                } else if (method == 3) {
                    alignOrientedReads3(
                        orientedReadId0, orientedReadId1,
                        matchScore, mismatchScore, gapScore,
                        downsamplingFactor, bandExtend,
                        alignment, alignmentInfo);
                } else {
                    SHASTA_ASSERT(0);
                }

                // If the alignment is poor, skip it.
                if ((alignment.ordinals.size() < minAlignedMarkerCount) ||
                    (alignmentInfo.minAlignedFraction() < minAlignedFraction)) {
                    continue;
                }

                // If the alignment has too much trim, skip it.
                uint32_t leftTrim;
                uint32_t rightTrim;
                tie(leftTrim, rightTrim) = alignmentInfo.computeTrim();
                if(leftTrim>maxTrim || rightTrim>maxTrim) {
                    continue;
                }
                alignments.push_back(make_pair(orientedReadId1, alignmentInfo));
            }
        }
    }
}



void Assembler::exploreAlignmentGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);

    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);

    size_t minAlignedMarkerCount = httpServerData.assemblerOptions->alignOptions.minAlignedMarkerCount;
    getParameterValue(request, "minAlignedMarkerCount", minAlignedMarkerCount);

    size_t maxTrim = httpServerData.assemblerOptions->alignOptions.maxTrim;
    getParameterValue(request, "maxTrim", maxTrim);

    uint32_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    uint32_t sizePixels = 1200;
    getParameterValue(request, "sizePixels", sizePixels);

    double timeout= 30;
    getParameterValue(request, "timeout", timeout);



    // Write the form.
    html <<
        "<h3>Display a local subgraph of the <a href='docs/ReadGraph.html'>global alignment graph</a></h3>"
        "<form>"

        "<table>"

        "<tr title='Read id between 0 and " << reads.size()-1 << "'>"
        "<td>Read id"
        "<td><input type=text required name=readId size=8 style='text-align:center'"
        << (readIdIsPresent ? ("value='"+to_string(readId)+"'") : "") <<
        ">"

        "<tr title='Choose 0 (+) for the input read or 1 (-) for its reverse complement'>"
        "<td>Strand"
        "<td class=centered>";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);

    html <<

        "<tr title='Maximum distance from start vertex (number of edges)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        " value='" << maxDistance <<
        "'>"

        "<tr title='The minimum number of aligned markers "
        "in order for an edge to be generated'>"
        "<td>Minimum number of aligned markers"
        "<td><input type=text required name=minAlignedMarkerCount size=8 style='text-align:center'"
        " value='" << minAlignedMarkerCount <<
        "'>"

        "<tr title='The maximum number of trimmed bases on either side "
        "in order for an edge to be generated'>"
        "<td>Minimum alignment trim"
        "<td><input type=text required name=maxTrim size=8 style='text-align:center'"
        " value='" << maxTrim <<
        "'>"

        "<tr title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "<td>Graphics size in pixels"
        "<td><input type=text required name=sizePixels size=8 style='text-align:center'" <<
        " value='" << sizePixels <<
        "'>"

        "<tr title='Maximum time (in seconds) allowed for graph creation and layout'>"
        "<td>Timeout (seconds) for graph creation and layout"
        "<td><input type=text required name=timeout size=8 style='text-align:center'" <<
        " value='" << timeout <<
        "'>"

        "</table>"

        "<input type=submit value='Display'>"
        "</form>";



    // If any necessary values are missing, stop here.
    if(!readIdIsPresent || !strandIsPresent) {
        return;
    }



    // Validity checks.
    if(readId > reads.size()) {
        html << "<p>Invalid read id " << readId;
        html << ". Must be between 0 and " << reads.size()-1 << ".";
        return;
    }
    if(strand>1) {
        html << "<p>Invalid strand " << strand;
        html << ". Must be 0 or 1.";
        return;
    }
    const OrientedReadId orientedReadId(readId, strand);



    // Create the local alignment graph.
    LocalAlignmentGraph graph;
    const auto createStartTime = steady_clock::now();
    if(!createLocalAlignmentGraph(orientedReadId,
        minAlignedMarkerCount, maxTrim, maxDistance, timeout, graph)) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    const auto createFinishTime = steady_clock::now();

    // Write it out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    graph.write(dotFileName, maxDistance);

    // Compute layout in svg format.
    const string command =
        timeoutCommand() + " " + to_string(timeout - seconds(createFinishTime - createStartTime)) +
        " sfdp -O -T svg " + dotFileName +
        " -Gsize=" + to_string(sizePixels/72.);
    const auto layoutStartTime = steady_clock::now();
    const int commandStatus = ::system(command.c_str());
    const auto layoutFinishTime = steady_clock::now();
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
    filesystem::remove(dotFileName);



    // Write a title and display the graph.
    html <<
        "<h1 style='line-height:10px'>Alignment graph near oriented read " << orientedReadId << "</h1>"
        "Color legend: "
        "<span style='background-color:LightGreen'>start vertex</span> "
        "<span style='background-color:cyan'>vertices at maximum distance (" << maxDistance <<
        ") from the start vertex</span>.";


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

    // Write additional graph information.
    html <<
        "<br>This portion of the alignment graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." <<
        "<br>Graph creation took " <<
        std::setprecision(2) << seconds(createFinishTime-createStartTime) <<
        " s.<br>Graph layout took " <<
        std::setprecision(2) << seconds(layoutFinishTime-layoutStartTime) << " s.";

    // Write a histogram of the number of vertices by distance.
    vector<int> histogram(maxDistance+1, 0);
    BGL_FORALL_VERTICES(v, graph, LocalAlignmentGraph) {
        ++histogram[graph[v].distance];
    }
    html <<
        "<h4>Vertex count by distance from start vertex</h4>"
        "<table>"
        "<tr><th>Distance<th>Count";
    for(uint32_t distance=0; distance<=maxDistance; distance++) {
        html << "<tr><td class=centered>" << distance << "<td class=centered>" << histogram[distance];
    }
    html << "</table>";

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
