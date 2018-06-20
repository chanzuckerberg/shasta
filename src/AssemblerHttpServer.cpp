#include "Assembler.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

#include <iomanip>



void Assembler::fillServerFunctionTable()
{
    // Summary.
    httpServerData.functionTable[""]        = &Assembler::exploreSummary;
    httpServerData.functionTable["/"]       = &Assembler::exploreSummary;
    httpServerData.functionTable["/index"]  = &Assembler::exploreSummary;

}



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
    writeMakeAllTablesSelectable(html);
    html <<
        "</head>"
        "<body onload='makeAllTablesSelectableByDoubleClick()'>";
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

}



void Assembler::exploreSummary (
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

        "<tr><td title='The number of k-mers of length k used as marker'>Marker k-mers"
        "<td class=right>" << markerKmerCount <<

        "<tr><td title='The fraction of k-mers of length k used as marker'>Marker fraction"
        "<td class=right>" << setprecision(4) << double(markerKmerCount) / double(kmerTable.size()) <<

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


