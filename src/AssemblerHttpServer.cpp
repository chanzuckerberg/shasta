// shasta.
#include "Assembler.hpp"
#include "AlignmentGraph.hpp"
#include "LocalReadGraph.hpp"
#include "LocalMarkerGraph2.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "chrono.hpp"
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
    CZI_ADD_TO_FUNCTION_TABLE(exploreOverlappingReads);
    CZI_ADD_TO_FUNCTION_TABLE(exploreAlignment);
    CZI_ADD_TO_FUNCTION_TABLE(exploreReadGraph);
    CZI_ADD_TO_FUNCTION_TABLE(exploreMarkerGraph);

}
#undef CZI_ADD_TO_FUNCTION_TABLE

void Assembler::setDocsDirectory(const string& docsDirectoryArgument)
{
    httpServerData.docsDirectory = docsDirectoryArgument;
}



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
        {"Overlapping reads", "exploreOverlappingReads"},
        {"Align two reads", "exploreAlignment"},
        });
    writeNavigation(html, "Read graph", {
        {"Read graph", "exploreReadGraph"},
        });
    writeNavigation(html, "Marker graph", {
        {"Marker graph", "exploreMarkerGraph"},
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

        "<tr><td title='Number of overlaps with good alignments'>Alignments"
        "<td class=right>" << alignmentData.size() <<

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

    // Get the begin and end position.
    uint32_t beginPosition = 0;
    const bool beginPositionIsPresent = getParameterValue(request, "beginPosition", beginPosition);
    uint32_t endPosition = 0;
    const bool endPositionIsPresent = getParameterValue(request, "endPosition", endPosition);

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
    html << "<br><input type=text name=beginPosition size=8";
    if(beginPositionIsPresent) {
        html << " value=" << beginPosition;
    }
    html <<
        ">Begin display at this base position (leave blank to begin at beginning of read)."
        "<br><input type=text name=endPosition size=8";
    if(endPositionIsPresent) {
        html << " value=" << endPosition;
    }
    html <<
        ">End display at this base position (leave blank to end at end of read)."
        "</form>";

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
    if(!beginPositionIsPresent) {
        beginPosition = 0;
    }
    if(!endPositionIsPresent) {
        endPosition = uint32_t(readSequence.baseCount);
    }



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
    if(beginPositionIsPresent || endPositionIsPresent) {
        html <<
            " Displaying only " << endPosition-beginPosition <<
            " bases begining at base position " << beginPosition <<
            " and ending at base position " << endPosition <<
            " . Markers partially contained in this interval will not be displayed.";
    }



    // Button to Blat this read.
    // We cannot use a simple <a> because we need to do a POST
    // (the GET request fails when the read is too long).
    html <<
        "<p><form action='https://genome.ucsc.edu/cgi-bin/hgBlat' method=post>"
        "<input type=submit value='Blat ";
    if(beginPositionIsPresent || endPositionIsPresent) {
        html << " this portion of ";
    }
    html <<
        "this read against human reference hg38'>"
        "<input type=text hidden name=type value=DNA>"
        "<input type=text hidden name=type value=DNA>"
        "<input type=text hidden name=name value=Human>"
        "<input type=text hidden name=db value=hg38>"
        "<input type=text hidden name=userSeq value=";
    readSequence.write(html, strand==1, beginPosition, endPosition);
    html << "></form>";



    // Link to align this read against another read.
    html <<
        "<p><a href='exploreAlignment?readId0=" << readId << "&strand0=" << strand <<
        "'>Compute a marker alignment of this read with another read.</a>";

    // Link to show overlapping reads.
    html <<
        "<p><a href='exploreOverlappingReads?readId=" << readId << "&strand=" << strand <<
        "'>Find other reads that overlap this read.</a>";



    // Decide on which row each marker gets displayed
    // (first row of markers is row 0).
    const size_t k = assemblerInfo->k;
    vector<int> markerRow(orientedReadMarkers.size(), -1);
    vector<uint32_t> nextAvailableCharacterPosition(k, 0);
    for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
        const uint32_t position = orientedReadMarkers[ordinal].position;
        for(int row=0; row<int(k); row++) {
            if(position >= nextAvailableCharacterPosition[row]) {
                markerRow[ordinal] = row;
                // Require one character space to next marker.
                nextAvailableCharacterPosition[row] = position + uint32_t(k) + 1;
                // cout << "Marker " << ordinal << " at position " << position << " placed on row " << row << endl;
                break;
            }
        }
        /*
        if(markerRow[ordinal] == -1) {
            cout << "Marker " << ordinal << " at position " << position << " could not be placed." << endl;
            cout << "nextAvailableCharacterPosition: ";
            for(const auto rowNext: nextAvailableCharacterPosition) {
                cout << " " << rowNext;
            }
            cout << endl;
        }
        */
        CZI_ASSERT(markerRow[ordinal] != -1);
    }
    const int markerRowCount = *std::max_element(markerRow.begin(), markerRow.end());



    // Use an svg object to display the read sequence and the markers.
    // To ensure correct positioning and alignment, we use
    // a textLength attribute on every <text> element.
    // (The older code, ifdef'ed out, uses a separate <text>
    // element for each character).
    const int monospaceFontSize = 12;
    const int horizontalSpacing = 7;
    const int verticalSpacing = 13;
    const int charactersPerLine = endPosition - beginPosition + 10; // Add space for labels
    const int svgLineCount = int(3 + markerRowCount); // Labels, scale, sequence, markers.
    const int svgWidth = horizontalSpacing * charactersPerLine;
    const int svgHeight = verticalSpacing * svgLineCount;
    const int highlightedMarkerVerticalOffset = 2;
    html <<
        "<p><svg width=" << svgWidth << " height=" << svgHeight << ">"
        "<style>"
        ".mono{font-family:monospace; font-size:" << monospaceFontSize << "px;}"
        ".blueMono{font-family:monospace; font-size:" << monospaceFontSize << "px; fill:blue;}"
        "</style>";



    // Labels for position scale.
    const int firstLabelPosition =
        (beginPosition % 10 == 0) ?
        beginPosition :
        beginPosition + (10 - beginPosition % 10);
    for(int labelPosition=firstLabelPosition; labelPosition<int(endPosition); labelPosition+=10) {
        const string label = to_string(labelPosition);
#if 1
        // This code uses a single <text> element,
        // with a textLength attribute for exact alignment.
        html <<
            "<text class='mono'" <<
            " x='" << (labelPosition-beginPosition)*horizontalSpacing << "'" <<
            " y='" << verticalSpacing << "'" <<
            " textLength='" << label.size()*horizontalSpacing << "px'>" <<
            label << "</text>";
#else
        // This code uses one <text> element per character.
        for(size_t j=0; j<label.size(); j++) {
            html <<
                "<text class='mono'" <<
                " x='" << (labelPosition-beginPosition+j)*horizontalSpacing << "'" <<
                " y='" << verticalSpacing << "'>" <<
                label[j] << "</text>";
        }
#endif
    }



    // Position scale.
#if 0
    // This code uses a single <text> element,
    // with a textLength attribute for exact alignment.
    html <<
        "<text class='mono'" <<
        " x='0'" <<
        " y='" << 2*verticalSpacing << "'" <<
        " textLength ='" << (endPosition-beginPosition) * horizontalSpacing << "px'>";
    for(size_t position=beginPosition; position!=endPosition; position++) {
        if((position%10)==0) {
            html << "|";
        } else if((position%5)==0) {
            html << "+";
        } else {
            html << ".";
        }
    }
    html << "</text>";
#endif
#if 0
    // This code uses one <text> element per character.
    for(size_t position=beginPosition; position!=endPosition; position++) {
        html <<
            "<text class='mono'" <<
            " x='" << (position-beginPosition)*horizontalSpacing << "'" <<
            " y='" << 2*verticalSpacing << "'>";
        if((position%10)==0) {
            html << "|";
        } else if((position%5)==0) {
            html << "+";
        } else {
            html << ".";
        }
        html << "</text>";
    }
#endif
#if 1
    // This code uses one <text> element for every blockSize characters.
    // This way you can select sequence text without getting a
    // new line after each character, while still achieving good
    // alignment.
    // It is a reasonable compromise between the two extreme choices above.
    const size_t blockSize = 100;
    for(size_t blockBegin=beginPosition; blockBegin<endPosition; blockBegin+=blockSize) {
        const size_t blockEnd = min(blockBegin+blockSize, size_t(endPosition));
        html <<
            "<text class='mono'" <<
            " x='" << (blockBegin-beginPosition)*horizontalSpacing << "'" <<
            " y='" << 2*verticalSpacing << "'"
            " textLength='" << (blockEnd-blockBegin) * horizontalSpacing<< "'>";
        for(size_t position=blockBegin; position!=blockEnd; position++) {
            if((position%10)==0) {
                html << "|";
            } else if((position%5)==0) {
                html << "+";
            } else {
                html << ".";
            }
        }
        html << "</text>";
    }
#endif



    // Read sequence.
#if 0
    // This code uses a single <text> element,
    // with a textLength attribute for exact alignment.
    // Unfortunately this does not achieve exact alignment on Chrome.
    html <<
        "<text class='mono'" <<
        " x='0'" <<
        " y='" << 3*verticalSpacing << "'" <<
        " textLength='" << (endPosition-beginPosition) * horizontalSpacing << "px'>" <<
        readSequence <<
        "</text>";
#endif
#if 0
    // This code uses one <text> element per character.
    for(size_t position=beginPosition; position!=endPosition; position++) {
        html <<
            "<text class='mono'" <<
            " x='" << (position-beginPosition)*horizontalSpacing << "'" <<
            " y='" << 3*verticalSpacing << "'>";
        html << readSequence[position];
        html << "</text>";
    }
#endif
#if 1
    // This code uses one <text> element for every blockSize characters.
    // This way you can select sequence text without getting a
    // new line after each character, while still achieving good
    // alignment.
    // It is a reasonable compromise between the two extreme choices above.
    for(size_t blockBegin=beginPosition; blockBegin<endPosition; blockBegin+=blockSize) {
        const size_t blockEnd = min(blockBegin+blockSize, size_t(endPosition));
        html <<
            "<text class='mono'" <<
            " x='" << (blockBegin-beginPosition)*horizontalSpacing << "'" <<
            " y='" << 3*verticalSpacing << "'"
            " textLength='" << (blockEnd-blockBegin) * horizontalSpacing<< "'>";
        for(size_t position=blockBegin; position!=blockEnd; position++) {
            html << readSequence[position];
        }
        html << "</text>";
    }
#endif



    // Draw a rectangle for each highlighted marker.
    for(const uint32_t ordinal: highlightedMarkers) {
        const CompressedMarker& marker = orientedReadMarkers[ordinal];
        html <<
            "<rect" <<
            " x='" << (marker.position-beginPosition)*horizontalSpacing << "'"
            " y='" << (3 + markerRow[ordinal])*verticalSpacing + highlightedMarkerVerticalOffset << "'"
            " height='" << verticalSpacing << "'"
            " width='" << k * horizontalSpacing << "'"
            " style='fill:pink; stroke:none;'"
            "/>";
    }



    // Markers.
    for(uint32_t ordinal=0; ordinal<uint32_t(orientedReadMarkers.size()); ordinal++) {
        const CompressedMarker& marker = orientedReadMarkers[ordinal];

        // Only show the marker if it is completely contained in the
        // base interval we are displaying.
        if(marker.position < beginPosition) {
            continue;
        }
        if(marker.position + k > endPosition) {
            break;
        }

        // See if this marker is contained in a vertex of the marker graph.
        const GlobalMarkerGraphVertexId vertexId =
            getGlobalMarkerGraphVertex(orientedReadId, ordinal);
        const bool hasMarkerGraphVertex =
            (vertexId != invalidCompressedGlobalMarkerGraphVertexId);



        // Write the k-mer of this marker.
        const Kmer kmer(marker.kmerId, k);
        html << "<a xlink:title='Marker " << ordinal << ", position " << marker.position;
        if(hasMarkerGraphVertex) {
            html << ", coverage " << globalMarkerGraphVertices.size(vertexId);
        }
        html << "' id='" << ordinal << "'";
        if(hasMarkerGraphVertex) {
            // Add a hyperlink to the marker graph vertex
            // that contains this marker.
            const string url = "exploreMarkerGraph?readId=" + to_string(readId) +
                "&strand=" + to_string(strand) +
                "&ordinal=" + to_string(ordinal) +
                "&maxDistance=2&detailed=on&minCoverage=3&minConsensus=3&sizePixels=3200&timeout=30";
            html << " xlink:href='" << url << "' style='cursor:pointer'";
        }
        html << ">";
#if 0
        // This code uses a single <text> element,
        // with a textLength attribute for exact alignment.
        html << "<text class='";
        if(hasMarkerGraphVertex) {
            html << "blueMono";
        } else {
            html << "mono";
        }
        html << "'" <<
            " x='" << (marker.position-beginPosition)*horizontalSpacing << "px'" <<
            " y='" << (4+markerRow[ordinal])*verticalSpacing << "'"
            " textLength='" << k*horizontalSpacing << "px'>";
        kmer.write(html, k);
        html << "</text>";
#else
        // This code uses one <text> element per character.
        for(size_t positionInMarker=0; positionInMarker<k; positionInMarker++) {
            html << "<text class='";
            if(hasMarkerGraphVertex) {
                html << "blueMono";
            } else {
                html << "mono";
            }
            html << "'" <<
                " x='" << (marker.position+positionInMarker-beginPosition)*horizontalSpacing << "'" <<
                " y='" << (4+markerRow[ordinal])*verticalSpacing << "'>";
            html << kmer[positionInMarker];
            html << "</text>";
        }
#endif
        html << "</a>";

    }


    // Finish the svg object.
    html << "</svg>";

    // Scroll to the first highlighted marker.
    if(!highlightedMarkers.empty()) {
        const uint32_t ordinal = *highlightedMarkers.begin();
        html <<
            "<script>"
            "var element = document.getElementById('" << ordinal << "');"
            "var rectangle = element.getBoundingClientRect();"
            "window.scroll(rectangle.left-100, rectangle.top-100);"
            "</script>";
    }



#if 0
    // THE IFDEF CONTAINS THE OLD CODE THAT USES MONOSPACE CHARACTERS
    // TO DISPLAY THE READ SEQUENCE AND THE MARKERS.
    // THIS DOES NOT WORK WELL ON CHROME: FOR LONG READS,
    // THE MARKERS ARE NOT CORRECTLY ALIGNED WITH THE READ SEQUENCE.

    // Labels for position scale.
    html << "<p><div style='font-family:monospace'>";
    html << "<br>";
    if(beginPositionIsPresent || endPositionIsPresent) {
        // Portion of the read.
        const uint32_t skip = (10 - beginPosition%10) % 10;
        for(uint32_t i=0; i<skip; i++) {
            html << "&nbsp;";
        }
        for(uint32_t i=beginPosition+skip; i<endPosition; i+=10) {
            const string label = to_string(i);
            html << label;
            for(size_t j=0; j<10-label.size(); j++) {
                html << "&nbsp;";
            }
        }
    } else {
        // Entire read.
        for(uint32_t i=0; i<readSequence.baseCount; i+=10) {
            const string label = to_string(i);
            html << label;
            for(size_t j=0; j<10-label.size(); j++) {
                html << "&nbsp;";
            }
        }
    }
    html << "<br>";



    // Position scale.
    for(uint32_t i=beginPosition; i<endPosition; i++) {
        if((i%10)==0) {
            html << "|";
        } else if((i%5)==0) {
            html << "+";
        } else {
            html << ".";
        }
    }



    // Read sequence.
    html << "<br>";
    readSequence.write(html, strand==1, beginPosition, endPosition);



    // Write the markers on k rows.
    const size_t k = assemblerInfo->k;
    for(size_t markerRow=0; markerRow<k; markerRow++) {
        html << "<br>";
        size_t position = beginPosition;
        for(uint32_t ordinal=uint32_t(markerRow);
            ordinal<uint32_t(orientedReadMarkers.size()); ordinal+=uint32_t(k)) {
            const CompressedMarker& marker = orientedReadMarkers[ordinal];

            // Only show the marker if it is completely contained in the
            // base interval we are displaying.
            if(marker.position < position) {
                continue;
            }
            if(marker.position + k > endPosition) {
                break;
            }

            const Kmer kmer(marker.kmerId, k);
            while(position < marker.position) {
                html << "&nbsp;";
                ++position;
            }

            // See if this marker is contained in a vertex of the marker graph.
            const GlobalMarkerGraphVertexId vertexId =
                getGlobalMarkerGraphVertex(orientedReadId, ordinal);
            const bool hasMarkerGraphVertex =
                (vertexId != invalidCompressedGlobalMarkerGraphVertexId);

            // If it has a vertex, write it as link. Otherwise, write it as text.
            if(hasMarkerGraphVertex) {
                const string url = "exploreMarkerGraph?readId=" + to_string(readId) +
                    "&strand=" + to_string(strand) +
                    "&ordinal=" + to_string(ordinal) +
                    "&maxDistance=2&detailed=on&minCoverage=3&minConsensus=3&sizePixels=3200&timeout=30";
                html << "<a id=" << ordinal << " style='font-family:monospace;";
                if(highlightedMarkers.find(ordinal) != highlightedMarkers.end()) {
                    html << "background-color:LightPink;";
                }
                html << "' title='Marker ordinal " << ordinal << ", coverage "  << globalMarkerGraphVertices.size(vertexId) << "'";
                html << "href='" << url << "'>";
                kmer.write(html, k);
                html << "</a>";
            } else {
                html << "<span id =" << ordinal;
                html << " title='Marker ordinal " << ordinal << "'";
                html << " style='font-family:monospace;";
                if(highlightedMarkers.find(ordinal) != highlightedMarkers.end()) {
                    html << "background-color:LightPink;'";
                }
                html << "'>";
                kmer.write(html, k);
                html << "</span>";
            }

            position += assemblerInfo->k;

        }
    }


    html << "</div>";
#endif


    html <<
        "<p>You can click on a blue marker above "
        "to see the global marker graph around that marker. "
        "Black markers correspond to a vertex of the marker graph "
        "that was removed because of low coverage.";
}



void Assembler::exploreOverlappingReads(
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
        "<input type=submit value='Show reads that overlap read'> "
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
        "<h1>Reads that overlap oriented read "
        "<a href='exploreRead?readId=" << readId0  << "&strand=" << strand0 << "'>"
        << orientedReadId0 << "</a></h1>";



    // Loop over Alignment objects that this oriented read in involved in.
    const auto overlapIndexes = alignmentTable[orientedReadId0.getValue()];
    html <<
        "<table><tr>"
        "<th>Oriented<br>read"
        "<th title='The number of aligned markers. Click on a cell in this column to see more alignment details.'>Aligned<br>markers";
    for(const auto i: overlapIndexes) {
        const AlignmentData& ad = alignmentData[i];

        ReadId readId1;
        if(ad.readIds[0] == readId0) {
            readId1 = ad.readIds[1];
        } else {
            CZI_ASSERT(ad.readIds[1] == readId0);
            readId1 = ad.readIds[0];
        }
        const Strand strand1 = ad.isSameStrand ? strand0 : 1-strand0;
        const OrientedReadId orientedReadId1(readId1, strand1);

        html <<
            "<tr>"
            "<td class=centered><a href='exploreRead?readId=" << readId1  << "&strand=" << strand1 <<
            "' title='Click to see this read'>" << orientedReadId1 << "</a>"
            "<td class=centered>"
            "<a href='exploreAlignment"
            "?readId0=" << readId0 << "&strand0=" << strand0 <<
            "&readId1=" << readId1 << "&strand1=" << strand1 <<
            "' title='Click to see the alignment'>" << ad.info.markerCount << "</a>";
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
        "and up to " << maxVertexCountPerKmer << " alignment graph vertices for each k-mer.";



    // Compute the alignment.
    // This creates file Alignment.png.
    vector<MarkerWithOrdinal> markers0SortedByKmerId;
    vector<MarkerWithOrdinal> markers1SortedByKmerId;
    getMarkersSortedByKmerId(orientedReadId0, markers0SortedByKmerId);
    getMarkersSortedByKmerId(orientedReadId1, markers1SortedByKmerId);
    AlignmentGraph graph;
    Alignment alignment;
    const bool debug = true;
    alignOrientedReads(
        markers0SortedByKmerId,
        markers1SortedByKmerId,
        maxSkip, maxVertexCountPerKmer, debug, graph, alignment);
    if(alignment.ordinals.empty()) {
        html << "<p>The alignment is empty (it has no markers).";
        return;
    }
    const AlignmentInfo alignmentInfo(alignment);
    uint32_t leftTrim;
    uint32_t rightTrim;
    tie(leftTrim, rightTrim) = computeTrim(
        orientedReadId0,
        orientedReadId1,
        alignmentInfo);



    // Write out a table with some information on the alignment.
    const auto markers0 = markers[orientedReadId0.getValue()];
    const auto markers1 = markers[orientedReadId1.getValue()];
    const auto markerCount0 = markers0.size();
    const auto markerCount1 = markers1.size();
    const auto baseCount0 = reads[orientedReadId0.getReadId()].baseCount;
    const auto baseCount1 = reads[orientedReadId1.getReadId()].baseCount;
    const auto firstOrdinal0 = alignment.ordinals.front().first;
    const auto firstOrdinal1 = alignment.ordinals.front().second;
    const auto lastOrdinal0 = alignment.ordinals.back().first;
    const auto lastOrdinal1 = alignment.ordinals.back().second;
    const auto& firstMarker0 = markers0[firstOrdinal0];
    const auto& firstMarker1 = markers1[firstOrdinal1];
    const auto& lastMarker0 = markers0[lastOrdinal0];
    const auto& lastMarker1 = markers1[lastOrdinal1];
    html <<
        "<h3>Alignment summary</h3>"
        "<table>"
        "<tr>"
        "<th rowspan=2>"
        "<th colspan=2>Markers"
        "<th colspan=2>Bases"

        "<tr>"
        "<th>" << orientedReadId0 <<
        "<th>" << orientedReadId1 <<
        "<th>" << orientedReadId0 <<
        "<th>" << orientedReadId1 <<

        "<tr>"
        "<td title='Total number of markers or bases in this read'>Total"
        "<td class=centered>" << markerCount0 <<
        "<td class=centered>" << markerCount1 <<
        "<td class=centered>" << baseCount0 <<
        "<td class=centered>" << baseCount1 <<

        "<tr>"
        "<td title='Number of unaligned markers or bases to the left of the aligned portion'>Unaligned on left"
        "<td class=centered>" << firstOrdinal0 <<
        "<td class=centered>" << firstOrdinal1 <<
        "<td class=centered>" << firstMarker0.position <<
        "<td class=centered>" << firstMarker1.position <<

        "<tr>"
        "<td title='Number of unaligned markers or bases to the right of the aligned portion'>Unaligned on right"
        "<td class=centered>" << markerCount0 - 1 - lastOrdinal0 <<
        "<td class=centered>" << markerCount1 - 1 - lastOrdinal1 <<
        "<td class=centered>" << baseCount0 - 1 - lastMarker0.position <<
        "<td class=centered>" << baseCount1 - 1 - lastMarker1.position <<

        "<tr>"
        "<td title='Number of aligned markers or bases in the aligned portion'>Aligned range"
        "<td class=centered>" << lastOrdinal0 + 1 - firstOrdinal0 <<
        "<td class=centered>" << lastOrdinal1 + 1 - firstOrdinal1 <<
        "<td class=centered>" << lastMarker0.position + 1 - firstMarker0.position <<
        "<td class=centered>" << lastMarker1.position + 1 - firstMarker1.position <<

        "<tr>"
        "<td title='Number of aligned markers'>Aligned"
        "<td class=centered>" << alignment.ordinals.size() <<
        "<td class=centered>" << alignment.ordinals.size() <<
        "<td>"
        "<td>"

        "<tr>"
        "<td title='Fraction of aligned markers in the aligned portion'>Aligned fraction"
        "<td class=centered>" << std::setprecision(2) << double(alignment.ordinals.size()) / double(lastOrdinal0 + 1 - firstOrdinal0) <<
        "<td class=centered>" << std::setprecision(2) << double(alignment.ordinals.size()) / double(lastOrdinal1 + 1 - firstOrdinal1) <<
        "<td>"
        "<td>"

        "</table>"
        "<p>See bottom of this page for alignment details.";



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
        "Each faint line corresponds to 10 markers.";
    html << "<p><img src=\"data:image/png;base64,";
    ifstream png("Alignment.png.base64");
    html << png.rdbuf();
    html << "\"/>";



    // Write out details of the alignment.
    html <<
        "<h3>Alignment details</h3>"
        "<table>"

        "<tr>"
        "<th rowspan=2>K-mer"
        "<th colspan=2>Ordinals"
        "<th colspan=2>Positions"

        "<tr>"
        "<th>" << orientedReadId0 <<
        "<th>" << orientedReadId1 <<
        "<th>" << orientedReadId0 <<
        "<th>" << orientedReadId1;

    for(const auto& ordinals: alignment.ordinals) {
        const auto ordinal0 = ordinals.first;
        const auto ordinal1 = ordinals.second;
        const auto& marker0 = markers0[ordinal0];
        const auto& marker1 = markers1[ordinal1];
        const auto kmerId = marker0.kmerId;
        CZI_ASSERT(marker1.kmerId == kmerId);
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

            "<td class=centered>" << marker0.position <<
            "<td class=centered>" << marker1.position;

    }

    html << "</table>";
}



void Assembler::exploreReadGraph(
    const vector<string>& request,
    ostream& html)
{
    // Get the parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);

    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);

    size_t minAlignedMarkerCount = 100;
    getParameterValue(request, "minAlignedMarkerCount", minAlignedMarkerCount);

    size_t maxTrim = 200;
    getParameterValue(request, "maxTrim", maxTrim);

    uint32_t maxDistance = 2;
    getParameterValue(request, "maxDistance", maxDistance);

    uint32_t sizePixels = 1200;
    getParameterValue(request, "sizePixels", sizePixels);

    double timeout= 30;
    getParameterValue(request, "timeout", timeout);



    // Write the form.
    html <<
        "<h3>Display a local subgraph of the global read graph</h3>"
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



    // Create the LocalReadGraph.
    LocalReadGraph graph;
    const auto createStartTime = steady_clock::now();
    if(!createLocalReadGraph(orientedReadId,
        minAlignedMarkerCount, maxTrim, maxDistance, timeout, graph)) {
        html << "<p>Timeout for graph creation exceeded. Increase the timeout or reduce the maximum distance from the start vertex.";
        return;
    }
    const auto createFinishTime = steady_clock::now();

    // Write it out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = "/dev/shm/" + uuid + ".dot";
    graph.write(dotFileName, maxDistance);

    // Compute layout in svg format.
    const string command =
        "timeout " + to_string(timeout - seconds(createFinishTime - createStartTime)) +
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
        "<h1 style='line-height:10px'>Read graph near oriented read " << orientedReadId << "</h1>"
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
        "<br>This portion of the read graph has " << num_vertices(graph) <<
        " vertices and " << num_edges(graph) << " edges." <<
        "<br>Graph creation took " <<
        std::setprecision(2) << seconds(createFinishTime-createStartTime) <<
        " s.<br>Graph layout took " <<
        std::setprecision(2) << seconds(layoutFinishTime-layoutStartTime) << " s.";

    // Write a histogram of the number of vertices by distance.
    vector<int> histogram(maxDistance+1, 0);
    BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
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
    string showVertexIdString;
    const bool showVertexId = getParameterValue(request, "showVertexId", showVertexIdString);
    string useStoredConnectivityString;
    const bool useStoredConnectivity = getParameterValue(request, "useStoredConnectivity", useStoredConnectivityString);
    uint32_t minCoverage = 0;
    const bool minCoverageIsPresent = getParameterValue(request, "minCoverage", minCoverage);
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

        "<tr title='Check to show vertex ids (only useful for debugging)'>"
        "<td>Show vertex ids"
        "<td class=centered><input type=checkbox name=showVertexId"
        << (showVertexId ? " checked=checked" : "") <<
        ">"

        "<tr title='Check to use stored connectivity of the marker graph "
        "(for testing only - if everything works, this should not affect the display'>"
        "<td>Use stored connectivity"
        "<td class=centered><input type=checkbox name=useStoredConnectivity"
        << (useStoredConnectivity ? " checked=checked" : "") <<
        ">"


        "<tr title='Minimum coverage (number of markers) for a vertex or edge to be considered strong. "
        "Affects the display of vertices and edges.'>"
        "<td>Coverage threshold"
        "<td><input type=text required name=minCoverage size=8 style='text-align:center'"
        << (minCoverageIsPresent ? (" value='" + to_string(minCoverage)+"'") : " value='3'") <<
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
        || !maxDistanceIsPresent || !minCoverageIsPresent
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
    LocalMarkerGraph2 graph(uint32_t(assemblerInfo->k), reads, markers, globalMarkerGraphVertex);
    extractLocalMarkerGraph(orientedReadId, ordinal, maxDistance, useStoredConnectivity, graph);
    if(num_vertices(graph) == 0) {
        html << "<p>The specified marker does not correspond to a vertex of the marker graph.";
        return;
    }

    // Write it out in graphviz format.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = "/dev/shm/" + uuid + ".dot";
    graph.write(dotFileName, minCoverage, maxDistance, detailed, showVertexId);

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
    filesystem::remove(dotFileName);

    // Finally, we can display it.
    html << "<h1>Marker graph near marker " << ordinal << " of oriented read " << orientedReadId << "</h1>";
    const string svgFileName = dotFileName + ".svg";
    ifstream svgFile(svgFileName);
    html << svgFile.rdbuf();
    svgFile.close();

    // Remove the .svg file.
    filesystem::remove(svgFileName);

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
            "&sizePixels=" + to_string(sizePixels) +
            "&timeout=" + to_string(timeout) +
            (detailed ? "&detailed=on" : "") +
            (showVertexId ? "&showVertexId=on" : "") +
            (useStoredConnectivity ? "&useStoredConnectivity=on" : "");
        if(detailed) {
            html <<
                "document.getElementById('a_vertexDistance" << vertex.vertexId <<
                "').onclick = function() {location.href='" << url << "';};\n";
        } else {
            html <<
                "document.getElementById('vertex" << vertex.vertexId <<
                "').onclick = function() {location.href='" << url << "';};\n";

            // We are displaying the graph in compact mode.
            // Add a right click to recenter and show detailed.
            const string detailUrl =
                "exploreMarkerGraph?readId=" + to_string(markerInfo.orientedReadId.getReadId()) +
                "&strand=" + to_string(markerInfo.orientedReadId.getStrand()) +
                "&ordinal="  + to_string(markerInfo.ordinal) +
                "&maxDistance=1" +
                "&minCoverage=" + to_string(minCoverage) +
                "&sizePixels=" + to_string(sizePixels) +
                "&timeout=" + to_string(timeout) +
                "&detailed=on" +
                (showVertexId ? "&showVertexId=on" : "") +
                (useStoredConnectivity ? "&useStoredConnectivity=on" : "");
            html <<
                "document.getElementById('vertex" << vertex.vertexId <<
                "').oncontextmenu = function() {location.href='" << detailUrl << "';"
                "return false;};\n";
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

