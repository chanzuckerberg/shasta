
// Shasta.
#include "Assembler.hpp"
#include "AssemblyGraph.hpp"
#include "Coverage.hpp"
#include "buildId.hpp"
#include "filesystem.hpp"
#include "platformDependent.hpp"
#include "Reads.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <filesystem>


// A map containing descriptions of output files.
namespace shasta {
    std::map<string, string> outputDescriptionTable = {
        {
            "index.html",
            "Html file containing a list of assembly output files and their descriptions."
        },
        {
            "stdout.log",
            "Assembly log output."
        },
        {
            "performance.log",
            "Performance log output. Contains messages that are mnostly useful for "
            "performance analysis."
        },
        {
            "Data",
            "Directory containing Shasta binary data."
        },
        {
            "AssemblySummary.json",
            "Assembly summary information in json format."
        },
        {
            "AssemblySummary.html",
            "Assembly summary information in html format."
        },
        {
            "Assembly.fasta",
            "Assembly in FASTA format (one strand only)."
        },
        {
            "Assembly.gfa",
            "Assembly in GFA format (one strand only)."
        },
        {
            "Assembly-BothStrands.gfa",
            "Assembly in GFA format (both strands)."
        },
        {
            "Assembly-BothStrands-NoSequence.gfa",
            "Assembly in GFA format (compact output without sequence, both strands)."
        },
        {
            "AssemblySummary.csv",
            "List of assembled segments in order of decreasing length."
        },
        {
            "Binned-ReadLengthHistogram.csv",
            "Read length distribution in 1 kb bins."
        },
        {
            "ReadLengthHistogram.csv",
            "Detailed read length distribution."
        },
        {
            "ReadSummary.csv",
            "Summary file containing one line of information for each read."
        },
        {
            "LowHashBucketHistogram.csv",
            "MinHash bucket population histogram."
        },
        {
            "DisjointSetsHistogram.csv",
            "Coverage histogram for all disjoint sets."
        },
        {
            "MarkerGraphVertexCoverageHistogram.csv",
            "Coverage histogram for marker graph vertices."
        },
        {
            "MarkerGraphEdgeCoverageHistogram.csv",
            "Coverage histogram for marker graph edges."
        },
        {
            "SuppressedAlignmentCandidates.csv",
            "Details of suppressed alignment candidates."
        },
        {
            "ReadGraphComponents.csv",
            "Information about connected components of the read graph."
        },
        {
            "shasta.conf",
            "Configuration file containing options in effect for this assembly."
        },
        {
            "ReadLowHashStatistics.csv",
            "MinHash/LowHash statistics for each read."
        },
        {
            "AlignedFractionHistogram.csv",
            "Histogram of aligned fraction for all computed alignments."
        },
        {
            "AlignmentDriftHistogram.csv",
            "Histogram of marker drift for all computed alignments."
        },
        {
            "AlignmentMarkerCountHistogram.csv",
            "Histogram of number of aligned markers for all computed alignments."
        },
        {
            "AlignmentSkipHistogram.csv",
            "Histogram of marker skip for all computed alignments."
        },
        {
            "AlignmentTrimHistogram.csv",
            "Histogram of marker trim for all computed alignments."
        },
        {
            "BubbleChains.csv",
            "Information about bubble chains in the assembly."
        },
        {
            "PhasingRegions.csv",
            "Information about Phasing regions in the assembly."
        },
        {
            "Assembly-Detailed-NoSequence.gfa",
            "Detailed assembly representation with small bubbles: "
            "compact GFA file without sequence."
        },
        {
            "Assembly-Detailed.csv",
            "Detailed assembly representation with small bubbles: "
            "csv companion for the GFA files."
        },
        {
            "Assembly-Detailed.fasta",
            "Detailed assembly representation with small bubbles: "
            "FASTA file."
        },
        {
            "Assembly-Detailed.gfa",
            "Detailed assembly representation with small bubbles: "
            "complete GFA file."
        },
        {
            "Assembly-Haploid-NoSequence.gfa",
            "Haploid assembly representation: "
            "compact GFA file without sequence."
        },
        {
            "Assembly-Haploid.csv",
            "Haploid assembly representation: "
            "csv companion for the GFA files."
        },
        {
            "Assembly-Haploid.fasta",
            "Haploid assembly representation: "
            "FASTA file."
        },
        {
            "Assembly-Haploid.gfa",
            "Haploid assembly representation: "
            "complete GFA file."
        },
        {
            "Assembly-Phased-NoSequence.gfa",
            "Phased assembly representation with large bubbles: "
            "compact GFA file without sequence."
        },
        {
            "Assembly-Phased.csv",
            "Phased assembly representation with large bubbles: "
            "csv companion for the GFA files."
        },
        {
            "Assembly-Phased.fasta",
            "Phased assembly representation with large bubbles: "
            "FASTA file."
        },
        {
            "Assembly-Phased.gfa",
            "Phased assembly representation with large bubbles: "
            "complete GFA file."
        }
    };
}



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
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAlignmentCoverage);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAlignment);
    SHASTA_ADD_TO_FUNCTION_TABLE(computeAllAlignments);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAlignmentCandidateGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAlignmentGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(alignSequencesInBaseRepresentation);
    SHASTA_ADD_TO_FUNCTION_TABLE(alignSequencesInMarkerRepresentation);
    SHASTA_ADD_TO_FUNCTION_TABLE(assessAlignments);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreReadGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerGraphVertex);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerGraphEdge);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerCoverage);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerGraphInducedAlignment);
    SHASTA_ADD_TO_FUNCTION_TABLE(followReadInMarkerGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMarkerConnectivity);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAssemblyGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAssemblyGraphEdge);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreAssemblyGraphEdgesSupport);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreCompressedAssemblyGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMode3AssemblyGraph);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMode3AssemblyGraphSegment);
    SHASTA_ADD_TO_FUNCTION_TABLE(exploreMode3AssemblyGraphLink);

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
            writeNavigation(html);
            html << "Unknown documentation file " << name;
            writeHtmlEnd(html);
            return;
        }

        // Construct the full file name and open it.
        const string fileName = httpServerData.docsDirectory + "/" + name;
        ifstream file(fileName);
        if(!file) {
            writeHtmlBegin(html);
            writeNavigation(html);
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
        writeNavigation(html);
        html << "Unsupported keyword " << keyword;
        writeHtmlEnd(html);
        return;
    }


    // We found the keyword. Call the function that processes this keyword.
    // The processing function is only responsible for writing the html body.
    writeHtmlBegin(html);
    writeNavigation(html);
    try {
        const auto function = it->second;
        (this->*function)(request, html);
    } catch(const std::exception& e) {
        html << "<br><br><span style='color:purple'>" << e.what() << "</span>";
    }
    writeHtmlEnd(html);
}
#endif



void Assembler::writeMakeAllTablesCopyable(ostream& html) const
{
    html << R"###(
    <script>

    // Copy to the clipboard the table that generated the event.
    function copyToClipboard(event)
    {
        // If the CTRL key is not pressed, don't do anything.
        if(!event.ctrlKey) {
            return;
        }

        // Prevent default behavior.
        // event.preventDefault();
        // event.stopPropagation();
        // event.returnValue = false;
        
        // Get the table element.
        var element = event.currentTarget;
         
        // Remove any previous selection.
        var selection = window.getSelection();
        selection.removeAllRanges();
        
        // Select the table.
        var range = document.createRange();
        range.selectNodeContents(element);
        selection.addRange(range);
        
        // Copy it to the clipboard.
        document.execCommand("copy");

        // Unselect it.
        selection.removeAllRanges();

        window.alert("The table was copied to the clipboard");
    }

    // Make a table copyable by Ctrl-click.
    function makeCopyable(element)
    {
        element.addEventListener('click', copyToClipboard);
        element.title = 'Ctrl-click anywhere on the table to copy the entire table to the clipboard';
    }

    // Make all tables copyable by Ctrl-click.
    function makeAllTablesCopyable()
    {
        var tables = document.getElementsByTagName('table');
        var i;
        for(i=0; i<tables.length; i++) {
            makeCopyable(tables[i]);
        }
    }
    </script>
    )###";


#if 0
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
#endif
}



#ifdef SHASTA_HTTP_SERVER


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
        {"Candidate graph", "exploreAlignmentCandidateGraph"},
        {"Stored alignments", "exploreAlignments"},
        {"Alignment coverage", "exploreAlignmentCoverage"},
        {"Align two reads", "exploreAlignment"},
        {"Align one read with all", "computeAllAlignments"},
        {"Alignment graph", "exploreAlignmentGraph"},
        {"Assess alignments", "assessAlignments"},
        {"Align sequences in base representation", "alignSequencesInBaseRepresentation"},
        // {"Align sequences in marker representation", "alignSequencesInMarkerRepresentation"},
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
        {"Follow a read in the marker graph", "followReadInMarkerGraph"},
        {"Marker connectivity", "exploreMarkerConnectivity"},
        });
    if(assemblerInfo->assemblyMode == 0) {
        writeNavigation(html, "Assembly graph", {
            {"Local assembly graph", "exploreAssemblyGraph"},
            {"Assembly graph edges", "exploreAssemblyGraphEdge"},
            {"Assembly graph edges support", "exploreAssemblyGraphEdgesSupport"},
            {"Compressed assembly graph", "exploreCompressedAssemblyGraph"},
            });
    }
    if(assemblerInfo->assemblyMode == 3) {
        writeNavigation(html, "Assembly graph", {
            {"Local assembly graph", "exploreMode3AssemblyGraph"},
            {"Local assembly graph segments", "exploreMode3AssemblyGraphSegment"},
            {"Local assembly graph links", "exploreMode3AssemblyGraphLink"},
            });
    }
    
    if (!httpServerData.docsDirectory.empty()) {
        writeNavigation(html, "Help", {
            {"Documentation", "docs/index.html"},
            });
    }

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



void Assembler::exploreSummary(
    const vector<string>& request,
    ostream& html)
{
    writeAssemblySummaryBody(html);
}

#endif


// Access all available assembly data, without throwing exceptions
void Assembler::accessAllSoft()
{

    bool allDataAreAvailable = true;

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
        accessAlignmentCandidateTable();
    } catch(const exception& e) {
        cout << "Alignment candidate table is not accessible." << endl;
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
        cout << "The read graph is not accessible." << endl;
        allDataAreAvailable = false;
    }



    try {
        accessMarkerGraphVertices();
    } catch(const exception& e) {
        cout << "Marker graph vertices are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessMarkerGraphReverseComplementVertex();
    } catch(const exception& e) {
        cout << "Marker graph reverse complement vertices are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessMarkerGraphEdges(false);
    } catch(const exception& e) {
        cout << "Marker graph edges are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessMarkerGraphReverseComplementEdge();
    } catch(const exception& e) {
        cout << "Marker graph reverse complement edges are not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessMarkerGraphConsensus();
    } catch(const exception& e) {
        cout << "MarkerGraph graph consensus is not accessible." << endl;
        allDataAreAvailable = false;
    }

    try {
        accessCompressedAlignments();
    } catch(const exception& e) {
        cout << "Alignments are not accessible." << endl;
        allDataAreAvailable = false;
    }



    // Data specific to assembly mode 0.
    if(assemblerInfo->assemblyMode == 0) {
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

    }



    // Data specific to assembly mode 3.
    if(assemblerInfo->assemblyMode == 3) {
        try {
            accessMode3AssemblyGraph();
        } catch(const exception& e) {
            cout << "The mode 3 assembly graph is not accessible." << endl;
            allDataAreAvailable = false;
        }
    }



    if(!allDataAreAvailable) {
        cout << "Not all assembly data are accessible." << endl;
        cout << "Some functionality is not available." << endl;
    }
}




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
    td.left {
        text-align: left;
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
        margin: 2px;
        border-radius: 8px;
    }

    input[type=button] {
        padding: 4px;
    }

    input[type=text], input[type=radio] {
        background-color: #ecf1f0;
        border-width: thin;
    }

    button {
        background-color: #89bef2;
        padding: 4px;
        margin: 2px;
        border-radius: 8px;
    }

</style>
    )%";
}


void Assembler::writeHtmlBegin(ostream& html) const
{
    html <<
        "\r\n"
        "<!DOCTYPE html>"
        "<html>"
        "<head>"
        "<meta charset='UTF-8'>"
        "<title>Shasta assembler</title>";
    writeStyle(html);
    writeMakeAllTablesCopyable(html);
    html <<
        "</head>"
        "<body onload='makeAllTablesCopyable()'>";
}



void Assembler::writeHtmlEnd(ostream& html) const
{
    html << "</body>";
    html << "</html>";
}





void Assembler::writeAssemblySummary(ostream& html)
{
    writeHtmlBegin(html);
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
        "<tr><td>Read representation"
        "<td class=right>" << (assemblerInfo->readRepresentation==1 ? "1 (RLE)" : "0 (Raw - no RLE)") <<
        "<tr><td>Minimum read length"
        "<td class=right>" << assemblerInfo->minReadLength <<
        "<tr><td>Number of reads"
        "<td class=right>" << assemblerInfo->readCount <<
        "<tr><td>Number of read bases"
        "<td class=right>" << assemblerInfo->baseCount <<
        "<tr><td>Average read length"
        "<td class=right>" << assemblerInfo->baseCount / assemblerInfo->readCount <<
        "<tr><td>Read N50"
        "<td class=right>" << assemblerInfo->readN50;

    if(assemblerInfo->readRepresentation == 1) {
        html <<
        "<tr><td>Number of run-length encoded bases"
        "<td class=right>" << reads->getRepeatCountsTotalSize() <<
        "<tr><td>Average length ratio of run-length encoded sequence over raw sequence"
        "<td class=right>" << setprecision(4) << double(reads->getRepeatCountsTotalSize()) / double(assemblerInfo->baseCount);
    }

    html <<
        "<tr><td>Number of reads flagged as palindromic by self alignment"
        "<td class=right>" << assemblerInfo->palindromicReadCount <<
        "<tr><td>Number of reads flagged as chimeric"
        "<td class=right>" << assemblerInfo->chimericReadCount <<
        "</table>"
        "<ul>"
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
        "as opposed to run-length encoded sequence.</ul>";



    html <<
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
        "<td class=right>" << markers.totalSize();

    if(assemblerInfo->readRepresentation == 1) {
        html <<
            "<tr><td>Average number of markers per raw base"
            "<td class=right>" << setprecision(4) << double(markers.totalSize()/2)/double(assemblerInfo->baseCount) <<
            "<tr><td>Average number of markers per run-length encoded base"
            "<td class=right>" << setprecision(4) << double(markers.totalSize()/2)/double(reads->getRepeatCountsTotalSize()) <<
            "<tr><td>Average base offset between markers in raw sequence"
            "<td class=right>" << setprecision(4) << double(assemblerInfo->baseCount)/double(markers.totalSize()/2) <<
            "<tr><td>Average base offset between markers in run-length encoded sequence"
            "<td class=right>" << setprecision(4) << double(reads->getRepeatCountsTotalSize())/double(markers.totalSize()/2) <<
            "<tr><td>Average base gap between markers in run-length encoded sequence"
            "<td class=right>" << setprecision(4) <<
            double(reads->getRepeatCountsTotalSize())/double(markers.totalSize()/2) - double(assemblerInfo->k);
    } else {
        html <<
            "<tr><td>Average number of markers per base"
            "<td class=right>" << setprecision(4) << double(markers.totalSize()/2)/double(assemblerInfo->baseCount) <<
            "<tr><td>Average base offset between markers "
            "<td class=right>" << setprecision(4) << double(assemblerInfo->baseCount)/double(markers.totalSize()/2) <<
            "<tr><td>Average base gap between markers"
            "<td class=right>" << setprecision(4) <<
            double(assemblerInfo->baseCount)/double(markers.totalSize()/2) - double(assemblerInfo->k);
    }
    html << "</table>";



    html <<
        "<h3>Alignments</h3>"
        "<table>"
        "<tr><td>Number of alignment candidates found by the LowHash algorithm"
        "<td class=right>" << alignmentCandidates.candidates.size() <<
        "<tr><td>Number of good alignments"
        "<td class=right>" << alignmentData.size() <<
        "<tr><td>Number of good alignments kept in the read graph"
        "<td class=right>" << readGraph.edges.size()/2 <<
        "</table>"



        "<h3>Alignment criteria actually used for creation of the read graph</h3>"
        "<table>"
        "<tr><td>minAlignedMarkerCount<td class=right>" << assemblerInfo->actualMinAlignedMarkerCount <<
        "<tr><td>minAlignedFraction<td class=right>" << assemblerInfo->actualMinAlignedFraction <<
        "<tr><td>maxSkip<td class=right>" << assemblerInfo->actualMaxSkip <<
        "<tr><td>maxDrift<td class=right>" << assemblerInfo->actualMaxDrift <<
        "<tr><td>maxTrim<td class=right>" << assemblerInfo->actualMaxTrim <<
        "</table>"



        "<h3>Read graph</h3>"
        "<table>"
        "<tr><td>Number of vertices"
        "<td class=right>" << readGraph.connectivity.size() <<
        "<tr><td>Number of edges"
        "<td class=right>" << readGraph.edges.size() <<
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
        "<ul><li>The marker graph contains both strands.</ul>";


        if(assemblyGraphPointer) {
            html <<
                "<h3>Assembly graph</h3>"
                "<table>"
                "<tr><td>Number of vertices"
                "<td class=right>" << assemblyGraph.vertices.size() <<
                "<tr><td>Number of edges"
                "<td class=right>" << assemblyGraph.edges.size() <<
                "<tr><td>Number of edges assembled"
                "<td class=right>" << assemblerInfo->assemblyGraphAssembledEdgeCount <<
                "</table>"
                "<ul><li>The assembly graph contains both strands.</ul>";
        }



        if(assemblerInfo->assemblyMode == 0) {
            html <<
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
            "</ul>";
        }



        if(assemblerInfo->assemblyMode == 2) {
            const AssemblyGraph2Statistics& statistics = assemblerInfo->assemblyGraph2Statistics;
            html <<
                "<h3>Phased assembly statistics</h3>"
                "<table>"
                "<tr><th><th>Length<th>N<sub>50</sub>"
                "<tr><td>Bubble chains"
                "<td class=right>" << statistics.totalBubbleChainLength <<
                "<td class=right>" << statistics.bubbleChainN50 <<
                "<tr><td>Diploid sequence (per haplotype)"
                "<td class=right>" << statistics.totalDiploidLengthBothHaplotypes / 2 <<
                "<td class=right>" << statistics.diploidN50 <<
                "<tr><td>Haploid sequence"
                "<td class=right>" << statistics.totalHaploidLength <<
                "<td class=right>" << statistics.haploidN50 <<
                "<tr><td>Total sequence assembled in bubble chains, per haplotype"
                "<td class=right>" <<
                statistics.totalDiploidLengthBothHaplotypes / 2 +
                statistics.totalHaploidLength <<
                "<td class=right>"
                "<tr><td>Sequence outside bubble chains"
                "<td class=right>" << statistics.outsideBubbleChainsLength <<
                "<td class=right>"
                "</table>"

                "<p>"
                "<table>"
                "<tr><td>Number of bubbles that describe a single SNP (transition)"
                "<td class=right>" << statistics.simpleSnpBubbleTransitionCount <<
                "<tr><td>Number of bubbles that describe a single SNP (transversion)"
                "<td class=right>" << statistics.simpleSnpBubbleTransversionCount <<
                "<tr><td>Number of bubbles that describe a single SNP (total)"
                "<td class=right>" <<
                statistics.simpleSnpBubbleTransitionCount +
                statistics.simpleSnpBubbleTransversionCount <<
                "<tr><td>Transition/transversion ratio for bubbles that describe a single SNP"
                "<td class=right>" <<
                double(statistics.simpleSnpBubbleTransitionCount) /
                double(statistics.simpleSnpBubbleTransversionCount) <<
                "<tr><td>Number of bubbles that describe indels or more than one SNP"
                "<td class=right>" << statistics.nonSimpleSnpBubbleCount <<
                "</table>";
        }


    html <<
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
        "    \"Read representation\": \"" << (assemblerInfo->readRepresentation==1 ? "1 (RLE)" : "0 (Raw - no RLE)") << "\",\n"
        "    \"Minimum read length\": " << assemblerInfo->minReadLength << ",\n"
        "    \"Number of reads\": " << assemblerInfo->readCount << ",\n"
        "    \"Number of read bases\": " << assemblerInfo->baseCount << ",\n"
        "    \"Average read length\": " <<
        assemblerInfo->baseCount / assemblerInfo->readCount << ",\n"
        "    \"Read N50\": " << assemblerInfo->readN50 << ",\n";

    if(assemblerInfo->readRepresentation == 1) {
        json <<
            "    \"Number of run-length encoded bases\": " << reads->getRepeatCountsTotalSize() << ",\n"
            "    \"Average length ratio of run-length encoded sequence over raw sequence\": " <<
            setprecision(4) <<
            double(reads->getRepeatCountsTotalSize()) / double(assemblerInfo->baseCount) << ",\n";
    }

    json <<
        "    \"Number of reads flagged as palindromic by self alignment\": " << assemblerInfo->palindromicReadCount << ",\n"
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
        "  },\n";


    json <<
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
        << markers.totalSize() << ",\n";

    if(assemblerInfo->readRepresentation == 1) {
        json <<
            "    \"Average number of markers per raw base\": "
            << setprecision(4) << double(markers.totalSize()/2)/double(assemblerInfo->baseCount) << ",\n"
            "    \"Average number of markers per run-length encoded base\": "
            << setprecision(4) << double(markers.totalSize()/2)/double(reads->getRepeatCountsTotalSize()) << ",\n"
            "    \"Average base offset between markers in raw sequence\": "
            << setprecision(4) << double(assemblerInfo->baseCount)/double(markers.totalSize()/2) << ",\n"
            "    \"Average base offset between markers in run-length encoded sequence\": "
            << setprecision(4) << double(reads->getRepeatCountsTotalSize())/double(markers.totalSize()/2) << ",\n"
            "    \"Average base gap between markers in run-length encoded sequence\": "
            << setprecision(4) <<
            double(reads->getRepeatCountsTotalSize())/double(markers.totalSize()/2) - double(assemblerInfo->k) << "\n";
    } else {
        json <<
            "    \"Average number of markers per base\": "
            << setprecision(4) << double(markers.totalSize()/2)/double(assemblerInfo->baseCount) << ",\n"
            "    \"Average base offset between markers\": "
            << setprecision(4) << double(assemblerInfo->baseCount)/double(markers.totalSize()/2) << ",\n"
            "    \"Average base gap between markers\": "
            << setprecision(4) <<
            double(assemblerInfo->baseCount)/double(markers.totalSize()/2) - double(assemblerInfo->k) << "\n";

    }



    json <<
        "  },\n"


        "  \"Alignments\":\n"
        "  {\n"
        "    \"Number of alignment candidates found by the LowHash algorithm\": " <<
        alignmentCandidates.candidates.size() << ",\n"
        "    \"Number of good alignments\": " << alignmentData.size() << ",\n"
        "    \"Number of good alignments kept in the read graph\": " << readGraph.edges.size()/2 << "\n"
        "  },\n"



        "  \"Alignment criteria actually used for creation of the read graph\":\n"
        "  {\n"
        "    \"minAlignedMarkerCount\": " << assemblerInfo->actualMinAlignedMarkerCount << ",\n"
        "    \"minAlignedFraction\": " << assemblerInfo->actualMinAlignedFraction << ",\n"
        "    \"maxSkip\": " << assemblerInfo->actualMaxSkip << ",\n"
        "    \"maxDrift\": " << assemblerInfo->actualMaxDrift << ",\n"
        "    \"maxTrim\": " << assemblerInfo->actualMaxTrim << "\n"
        "  },\n"



        "  \"Read graph\":\n"
        "  {\n"
        "    \"Number of vertices\": " << readGraph.connectivity.size() << ",\n"
        "    \"Number of edges\": " << readGraph.edges.size() << ",\n"
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
        "  },\n";


    if(assemblyGraphPointer) {
        json <<
            "  \"Assembly graph\":\n"
            "  {\n"
            "    \"Number of vertices\": " << assemblyGraph.vertices.size() << ",\n"
            "    \"Number of edges\": " << assemblyGraph.edges.size() << ",\n"
            "    \"Number of edges assembled\": " << assemblerInfo->assemblyGraphAssembledEdgeCount << "\n"
            "  },\n";
    }


    if(assemblerInfo->assemblyMode == 0) {
        json <<
            "  \"Assembled segments\":\n"
            "  {\n"
            "    \"Number of segments assembled\": " << assemblerInfo->assemblyGraphAssembledEdgeCount << ",\n"
            "    \"Total assembled segment length\": " << assemblerInfo->totalAssembledSegmentLength << ",\n"
            "    \"Longest assembled segment length\": " << assemblerInfo->longestAssembledSegmentLength << ",\n"
            "    \"Assembled segments N50\": " << assemblerInfo->assembledSegmentN50 << "\n"
            "  },\n";
    }



    if(assemblerInfo->assemblyMode == 2) {
        const AssemblyGraph2Statistics& statistics = assemblerInfo->assemblyGraph2Statistics;

        json <<
            "  \"Phased assembly statistics\":\n"
            "  {\n"
            "    \"Bubble chains\":\n"
            "    {\n"
            "      \"Length\": " << statistics.totalBubbleChainLength << ",\n"
            "      \"N50\": " << statistics.bubbleChainN50 << "\n"
            "    },\n"
            "    \"Diploid sequence (per haplotype)\":\n"
            "    {\n"
            "      \"Length\": " << statistics.totalDiploidLengthBothHaplotypes / 2 << ",\n"
            "      \"N50\": " << statistics.diploidN50 << "\n"
            "    },\n"
            "    \"Haploid sequence\":\n"
            "    {\n"
            "      \"Length\": " << statistics.totalHaploidLength << ",\n"
            "      \"N50\": " << statistics.haploidN50 << "\n"
            "    },\n"
            "    \"Total sequence assembled in bubble chains, per haplotype\": " <<
            statistics.totalDiploidLengthBothHaplotypes / 2 +
            statistics.totalHaploidLength << ",\n"
            "    \"Sequence outside bubble chains\": " <<
            statistics.outsideBubbleChainsLength << ",\n"
            "    \"Number of bubbles that describe a single SNP (transition)\": " <<
            statistics.simpleSnpBubbleTransitionCount << ",\n"
            "    \"Number of bubbles that describe a single SNP (transversion)\": " <<
            statistics.simpleSnpBubbleTransversionCount << ",\n"
            "    \"Number of bubbles that describe a single SNP (total)\": " <<
            statistics.simpleSnpBubbleTransitionCount + statistics.simpleSnpBubbleTransversionCount << ",\n"
            "    \"Transition/transversion ratio for bubbles that describe a single SNP\": " <<
            double(statistics.simpleSnpBubbleTransitionCount) /
            double(statistics.simpleSnpBubbleTransversionCount) << ",\n"
            "    \"Number of bubbles that describe indels or more than one SNP\": " <<
            statistics.nonSimpleSnpBubbleCount << "\n"
            "  },\n";
    }



    json <<
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
    writeHtmlBegin(html);
    html << "<body><h1>Shasta assembly output files</h1><table>";

    // Loop over files in the assembly directory,
    // in case-insensitive alphabetic order.
    vector<string> assemblyFiles = shasta::filesystem::directoryContents(".");
    sort(assemblyFiles.begin(), assemblyFiles.end(),
        [](const string& s1, const string& s2) {
        return lexicographical_compare(
        s1.begin(), s1.end(),
        s2.begin(), s2.end(),
        [](const char& c1, const char& c2) {
        return tolower(c1) < tolower(c2);
        });});
    for(string file: assemblyFiles) {

        // Take out "./"
        file = file.substr(2);

        // Get the description.
        string description;
        auto it = outputDescriptionTable.find(file);
        if(it != outputDescriptionTable.end()) {
            description = it->second;
        }

        // Write a row in the table for this file.
        html <<
            "<tr><td><a href='" << file << "'>" << file <<
            "</a><td>" << description << endl;
    }

    writeHtmlEnd(html);
}



#ifdef SHASTA_HTTP_SERVER


void Assembler::blastRead(
    const vector<string>& request,
    ostream& html)
{

    if(!std::filesystem::is_regular_file(httpServerData.referenceFastaFileName)) {
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
    if(readId >= reads->readCount()) {
        html << "<p>Invalid read id.";
        return;
    }
    if(strand!=0 && strand!=1) {
        html << "<p>Invalid strand.";
        return;
    }
    const OrientedReadId orientedReadId(readId, strand);
    const vector<Base> rawOrientedReadSequence = reads->getOrientedReadRawSequence(orientedReadId);
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
    if(std::filesystem::file_size(blastErrFileName)) {
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
