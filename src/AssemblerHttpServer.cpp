// Nanopore2.
#include "Assembler.hpp"
#include "LocalMarkerGraph2.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include <iomanip>
#include "iterator.hpp"


#define CZI_ADD_TO_FUNCTION_TABLE(name) httpServerData.functionTable[string("/") + #name ] = &Assembler::name



// Associate http keywords wth member functions.
void Assembler::fillServerFunctionTable()
{
    httpServerData.functionTable[""]        = &Assembler::exploreSummary;
    httpServerData.functionTable["/"]       = &Assembler::exploreSummary;
    httpServerData.functionTable["/index"]  = &Assembler::exploreSummary;

    CZI_ADD_TO_FUNCTION_TABLE(exploreSummary);
    CZI_ADD_TO_FUNCTION_TABLE(exploreRead);
    CZI_ADD_TO_FUNCTION_TABLE(exploreAlignment);
    CZI_ADD_TO_FUNCTION_TABLE(exploreMarkerGraph);

}
#undef CZI_ADD_TO_FUNCTION_TABLE



void Assembler::processRequest(
    const vector<string>& request,
    ostream& html,
    const BrowserInformation&)
{
    // Look up the keyword to find the function that will process this request.
    // Note that the keyword includes the initial "/".
    const string& keyword = request.front();
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
    } catch(std::exception& e) {
        html << e.what();
    }
    writeHtmlEnd(html);
}



void Assembler::writeHtmlBegin(ostream& html) const
{
    html <<
        "\r\n"
        "<!DOCTYPE html>"
        "<html>"
        "<head>"
        "<link rel=icon href=\"https://s0.wp.com/wp-content/themes/vip/czi/images/build/favicon.ico\" />"
        "<meta charset='UTF-8'>";
    writeStyle(html);
    // writeMakeAllTablesSelectable(html);
    html <<
        "</head>"
        ;// "<body onload='makeAllTablesSelectableByDoubleClick()'>";
    writeNavigation(html);
}



void Assembler::writeHtmlEnd(ostream& html) const
{
    html << "</body>";
    html << "</html>";
}




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



void Assembler::writeNavigation(ostream& html) const
{
    html << "<ul class=navigationMenu>";

    writeNavigation(html, "Run information", {
        {"Summary", "exploreSummary"},
        });
    writeNavigation(html, "Reads", {
        {"Reads", "exploreRead"},
        {"Align two reads", "exploreAlignment"},
        });
    writeNavigation(html, "Marker graph", {
        {"Marker graph", "exploreMarkerGraph"},
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



void Assembler::exploreSummary(
    const vector<string>& request,
    ostream& html)
{
    using std::setprecision;

    // Compute the total number of bases.
    uint64_t totalBaseCount = 0;
    for(ReadId readId=0; readId<reads.size(); readId++) {
        totalBaseCount += reads[readId].baseCount;
    }

    // Compute the number of k-mers used as markers.
    uint64_t markerKmerCount = 0;
    for(const auto& tableEntry: kmerTable) {
        if(tableEntry.isMarker) {
            ++ markerKmerCount;
        }
    }


    html <<
        "<h1>Run summary</h1>"
        "<table>"

        "<tr><td title='Total number of input reads'>Reads"
        "<td class=right>" << reads.size() <<

        "<tr><td title='Total number of reads on both strands"
        " (equal to twice the number of reads)'>Oriented reads"
        "<td class=right>" << 2*reads.size() <<

        "<tr><td title='Total number of input bases'>Bases"
        "<td class=right>" << totalBaseCount <<

        "<tr><td title='Average number of bases in a read'>Average read length"
        "<td class=right>" << int(0.5 + double(totalBaseCount) / double(reads.size())) <<

        "<tr><td title='The length of k-mers used as markers'>Marker length k"
        "<td class=right>" << assemblerInfo->k <<

        "<tr><td title='The total number of k-mers of length k'>Total k-mers"
        "<td class=right>" << kmerTable.size() <<

        "<tr><td title='The number of k-mers of length k used as markers'>Marker k-mers"
        "<td class=right>" << markerKmerCount <<

        "<tr><td title='The fraction of k-mers of length k used as markers'>Marker fraction"
        "<td class=right>" << setprecision(3) << double(markerKmerCount) / double(kmerTable.size()) <<

        "<tr><td title='Total number of markers on both strands'>Oriented markers"
        "<td class=right>" << markers.totalSize() <<

        "<tr><td title='The average number of markers per base'>Marker density"
        "<td class=right>" << setprecision(4) << double(markers.totalSize()) / (2.*double(totalBaseCount)) <<

        "<tr><td title='The average shift between consecutive markers in a read'>Marker average shift"
        "<td class=right>" << setprecision(4) << (2.*double(totalBaseCount)) / double(markers.totalSize())  <<

        "<tr><td title='The average gap between consecutive markers in a read'>Marker average gap"
        "<td class=right>" << setprecision(4) <<
        (2.*double(totalBaseCount)) / double(markers.totalSize()) - double(assemblerInfo->k) <<

        "<tr><td title='Number of candidate overlaps found by the MinHash algorithm'>Overlaps"
        "<td class=right>" << overlaps.size() <<

        "<tr><td title='Number of vertices in the global marker graph'>Marker graph vertices"
        "<td class=right>" << globalMarkerGraphVertices.size() <<

        "</table>";
}



void Assembler::exploreRead(
    const vector<string>& request,
    ostream& html)
{
    // Get the ReadId and Strand from the request.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);

    // Get the set of ordinal for markers that should be highlighted.
    vector<string> highlightedMarkerStrings;
    getParameterValues(request, "highlightMarker", highlightedMarkerStrings);
    std::set<uint32_t> highlightedMarkers;
    for(const string& s: highlightedMarkerStrings) {
        try {
            highlightedMarkers.insert(boost::lexical_cast<uint32_t>(s));
        } catch(std::exception&) {
            // Ignore.
        }
    }

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Show read'> "
        "<input type=text name=readId required" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " size=8 title='Enter a read id between 0 and " << reads.size()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);
    html << "</form>";

    // If the readId or strand are missing, stop here.
    if(!readIdIsPresent || !strandIsPresent) {
        return;
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
    const auto readSequence = reads[readId];
    const auto readName = readNames[readId];
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];



    // Page title.
    html << "<h1 title='Read " << readId << " on strand " << strand;
    if(strand == 0) {
        html << " (input read without reverse complementing)";
    } else {
        html << " (reverse complement of input read)";
    }
    html << "'>Oriented read " << orientedReadId << "</h1>";

    // Read name.
    html << "<p>Read name on input: ";
    copy(readName.begin(), readName.end(), ostream_iterator<char>(html));

    // Read length.
    html << "<p>This read is " << readSequence.baseCount;
    html << " bases long and has " << orientedReadMarkers.size() << " markers.";

    // Hyperlink to Blat this read.
    html << "<p><a href='https://genome.ucsc.edu/cgi-bin/hgBlat?userSeq=";
    readSequence.write(html, strand==1);
    html << "&type=DNA&name=Human&db=hg38'>Blat this read against human reference hg38</a>.";

    // Read sequence.
    html << "<p><div style='font-family:monospace'>";
    html << "<br>";
    for(uint32_t i=0; i<readSequence.baseCount; i+=10) {
        const string label = to_string(i);
        html << label;
        for(size_t j=0; j<10-label.size(); j++) {
            html << "&nbsp;";
        }
    }
    html << "<br>";
    for(uint32_t i=0; i<readSequence.baseCount; i++) {
        if((i%10)==0) {
            html << "|";
        } else if((i%5)==0) {
            html << "+";
        } else {
            html << ".";
        }
    }
    html << "<br>";
    readSequence.write(html, strand==1);



    // Write the markers on k rows.
    const size_t k = assemblerInfo->k;
    for(size_t markerRow=0; markerRow<k; markerRow++) {
        html << "<br>";
        size_t position = 0;
        for(uint32_t ordinal=uint32_t(markerRow);
            ordinal<uint32_t(orientedReadMarkers.size()); ordinal+=uint32_t(k)) {
            const CompressedMarker& marker = orientedReadMarkers[ordinal];
            const string url = "exploreMarkerGraph?readId=" + to_string(readId) +
                "&strand=" + to_string(strand) +
                "&ordinal=" + to_string(ordinal) +
                "&maxDistance=2&detailed=on&minCoverage=3&minConsensus=3&sizePixels=3200&timeout=30";
            const Kmer kmer(marker.kmerId, k);
            while(position < marker.position) {
                html << "&nbsp;";
                ++position;
            }
            html << "<a id=" << ordinal << " style='font-family:monospace;";
            if(highlightedMarkers.find(ordinal) != highlightedMarkers.end()) {
                html << "background-color:LightPink;";
            }
            html << "' title='Marker ordinal " << ordinal << "'";
            html << "href='" << url << "'>";
            kmer.write(html, k);
            html << "</a>";
            position += assemblerInfo->k;

        }
    }


    html << "</div>";
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
    size_t maxSkip = 30;
    getParameterValue(request, "maxSkip", maxSkip);
    size_t maxVertexCountPerKmer = 100;
    getParameterValue(request, "maxVertexCountPerKmer", maxVertexCountPerKmer);

    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Compute a marker alignment'> of these two reads:"
        "<br><input type=text name=readId0 required size=8 " <<
        (readId0IsPresent ? "value="+to_string(readId0) : "") <<
        " title='Enter a read id between 0 and " << reads.size()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand0", strand0IsPresent && strand0==0, strand0IsPresent && strand0==1);
    html <<
         "<br><input type=text name=readId1 required size=8 " <<
         (readId1IsPresent ? "value="+to_string(readId1) : "") <<
         " title='Enter a read id between 0 and " << reads.size()-1 << "'>"
        " on strand ";
    writeStrandSelection(html, "strand1", strand1IsPresent && strand1==0, strand1IsPresent && strand1==1);
    html <<
        "<br>Maximum ordinal skip allowed: " <<
        "<input type=text name=maxSkip required size=8 value=" << maxSkip << ">";
    html <<
        "<br>Maximum number of alignment graph vertices per k-mer: " <<
        "<input type=text name=maxVertexCountPerKmer required size=8 value=" << maxVertexCountPerKmer << ">";
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
        "</h1>"
        "<p>This alignment was computed allowing a skip of up to " << maxSkip << " markers "
        "and up to " << maxVertexCountPerKmer << " alignment graph vertices for each k-mer."
        "<p>In the picture, horizontal positions correspond to marker ordinals on " <<
        orientedReadId0 << " (marker 0 is on left) "
        "and vertical positions correspond to marker ordinals on " <<
        orientedReadId0 << " (marker 0 is on top). "
        "Each faint line corresponds to 10 markers.";


    // Compute the alignment.
    // This creates file Alignment.png.
    alignOrientedReads(orientedReadId0, orientedReadId1, maxSkip, maxVertexCountPerKmer);

    // Create a base64 version of the png file.
    const string command = "base64 Alignment.png > Alignment.png.base64";
    ::system(command.c_str());


    // Write out the picture with the alignment.
    html << "<p><img src=\"data:image/png;base64,";
    ifstream png("Alignment.png.base64");
    html << png.rdbuf();
    html << "\"/>";
}



void Assembler::exploreMarkerGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);
    uint32_t ordinal = 0;
    const bool ordinalIsPresent = getParameterValue(request, "ordinal", ordinal);
    uint32_t maxDistance = 0;
    const bool maxDistanceIsPresent = getParameterValue(request, "maxDistance", maxDistance);
    string detailedString;
    const bool detailed = getParameterValue(request, "detailed", detailedString);
    uint32_t minCoverage = 0;
    const bool minCoverageIsPresent = getParameterValue(request, "minCoverage", minCoverage);
    uint32_t minConsensus = 0;
    const bool minConsensusIsPresent = getParameterValue(request, "minConsensus", minConsensus);
    uint32_t sizePixels;
    const bool sizePixelsIsPresent = getParameterValue(request, "sizePixels", sizePixels);
    uint32_t timeout;
    const bool timeoutIsPresent = getParameterValue(request, "timeout", timeout);



    // Write the form.
    html <<
        "<h3>Display a local subgraph of the global marker graph</h3>"
        "<form action='#startVertex'>"

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
        "<tr title='Ordinal for the desired marker in the specified oriented read.'>"
        "<td>Marker ordinal"
        "<td><input type=text required name=ordinal size=8 style='text-align:center'"
        << (ordinalIsPresent ? ("value='"+to_string(ordinal)+"'") : "") <<
        ">"

        "<tr title='Maximum distance from start vertex (number of edges)'>"
        "<td>Maximum distance"
        "<td><input type=text required name=maxDistance size=8 style='text-align:center'"
        << (maxDistanceIsPresent ? ("value='" + to_string(maxDistance)+"'") : " value='6'") <<
        ">"

        "<tr title='Check for detailed graph with labels'>"
        "<td>Detailed"
        "<td class=centered><input type=checkbox name=detailed"
        << (detailed ? " checked=checked" : "") <<
        ">"

        "<tr title='Minimum coverage (number of markers) for a vertex to be considered strong. Used to color vertices.'>"
        "<td>Vertex coverage threshold"
        "<td><input type=text required name=minCoverage size=8 style='text-align:center'"
        << (minCoverageIsPresent ? (" value='" + to_string(minCoverage)+"'") : " value='3'") <<
        ">"

        "<tr title='Minimum consensus (number or reads that agree on sequence) "
        "for an edge to be considered strong. Used to color edges.'>"
        "<td>Edge consensus threshold"
        "<td><input type=text required name=minConsensus size=8 style='text-align:center'"
        << (minConsensusIsPresent ? (" value='" + to_string(minConsensus)+"'") : " value='3'") <<
        ">"

        "<tr title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "<td>Graphics size in pixels"
        "<td><input type=text required name=sizePixels size=8 style='text-align:center'"
        << (sizePixelsIsPresent ? (" value='" + to_string(sizePixels)+"'") : " value='1600'") <<
        ">"

        "<tr title='Maximum time (in seconds) allowed for graph layout'>"
        "<td>Graph layout timeout"
        "<td><input type=text required name=timeout size=8 style='text-align:center'"
        << (timeoutIsPresent ? (" value='" + to_string(timeout)+"'") : " value='30'") <<
        ">"

        "</table>"

        "<input type=submit value='Display'>"

        " <span style='background-color:#e0e0e0' title='"
        "Fill this form to display a local subgraph of the global marker graph starting at the "
        "vertex containing the specified marker and extending out to a given distance "
        "(number of edges) from the start vertex."
        "The marker is specified by its oriented read id (read id, strand) and marker ordinal. "
        "The marker ordinal is the sequential index of the marker in the specified oriented read "
        "(the first marker is marker 0, the second marker is marker 1, and so on).'>"
        "Mouse here to explain form</span>"
        "</form>";



    // If any values are missing, stop here.
    if(!readIdIsPresent || !strandIsPresent || !ordinalIsPresent
        || !maxDistanceIsPresent || !minCoverageIsPresent || !minConsensusIsPresent
        || !timeoutIsPresent) {
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
    const auto orientedReadMarkerCount = markers.size(orientedReadId.getValue());
    if(ordinal >= orientedReadMarkerCount) {
        html <<
            "<p>Invalid marker ordinal. "
            "Oriented read " << orientedReadId <<
            " has "  << orientedReadMarkerCount <<
            " markers, and therefore the ordinal must be"
            " between 0 and " << orientedReadMarkerCount-1 << ".";
        return;
    }


    // Create the local marker graph.
    LocalMarkerGraph2 graph(uint32_t(assemblerInfo->k), reads, markers);
    extractLocalMarkerGraph(orientedReadId, ordinal, maxDistance, graph);

    // Write it out in graphviz format.
    // const boost::uuids::uuid uuid = boost::uuids::random_generator()();
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = "/dev/shm/" + uuid + ".dot";
    graph.write(dotFileName, minCoverage, minConsensus, maxDistance, detailed);

    // Compute layout in svg format.
    const string command =
        "timeout " + to_string(timeout) +
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
    // filesystem::remove(dotFileName);

    // Finally, we can display it.
    html << "<h1>Marker graph near marker " << ordinal << " of oriented read " << orientedReadId << "</h1>";
    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << svgFile.rdbuf();
    svgFile.close();

    // Remove the .svg file.
    // filesystem::remove(svgFileName);

    // Make the vertices clickable to recompute the graph with the
    // same parameters, but starting at the clicked vertex.
    // For a detailed graph, only the "Distance" label of each vertex
    // is made clickable.
    html << "<script>\n";
    BGL_FORALL_VERTICES(v, graph, LocalMarkerGraph2) {
        const LocalMarkerGraph2Vertex& vertex = graph[v];
        CZI_ASSERT(!vertex.markerInfos.empty());
        const auto& markerInfo = vertex.markerInfos.front();
        const string url =
            "exploreMarkerGraph?readId=" + to_string(markerInfo.orientedReadId.getReadId()) +
            "&strand=" + to_string(markerInfo.orientedReadId.getStrand()) +
            "&ordinal="  + to_string(markerInfo.ordinal) +
            "&maxDistance=" + to_string(maxDistance) +
            "&minCoverage=" + to_string(minCoverage) +
            "&minConsensus=" + to_string(minConsensus) +
            "&sizePixels=" + to_string(sizePixels) +
            "&timeout=" + to_string(timeout) +
            (detailed ? "&detailed=on" : "");
        if(detailed) {
            html <<
                "document.getElementById('a_vertexDistance" << vertex.vertexId <<
                "').onclick = function() {location.href='" << url << "';};\n";
        } else {
            html <<
                "document.getElementById('vertex" << vertex.vertexId <<
                "').onclick = function() {location.href='" << url << "';};\n";
        }
    }
    html << "</script>\n";



    // Position the start vertex at the center of the window.
    const GlobalMarkerGraphVertexId startVertexId =
        getGlobalMarkerGraphVertex(orientedReadId, ordinal);
    html << "<script>\n";
    if(detailed) {
        html <<
            "var element = document.getElementById('a_vertexDistance" << startVertexId << "');\n";
    } else {
        html <<
            "var element = document.getElementById('vertex" << startVertexId << "');\n";
    }
    html <<
        "var r = element.getBoundingClientRect();\n"
        "window.scrollBy((r.left + r.right - window.innerWidth) / 2, (r.top + r.bottom - window.innerHeight) / 2);\n"
        "</script>\n";
}



void Assembler::writeStrandSelection(
    ostream& html,
    const string& name,
    bool select0,
    bool select1) const
{
    html <<
        "<select name=" << name << " title='Choose 0 (+) for the input read or 1 (-) for its reverse complement'>"
        "<option value=0"
        << (select0 ? " selected" : "") <<
        ">0 (+)</option>"
        "<option value=1"
        << (select1 ? " selected" : "") <<
        ">1 (-)</option>"
        "</select>";

}

