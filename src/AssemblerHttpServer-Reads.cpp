#ifdef SHASTA_HTTP_SERVER

// Shasta.
#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
#include "Reads.hpp"
using namespace shasta;



void Assembler::exploreRead(
    const vector<string>& request,
    ostream& html)
{
    if(assemblerInfo->readRepresentation == 0) {
        exploreReadRaw(request, html);
    } else {
        exploreReadRle(request, html);
    }
}



void Assembler::exploreReadRle(
    const vector<string>& request,
    ostream& html)
{
    SHASTA_ASSERT(assemblerInfo->readRepresentation == 1);

    // Get the ReadId, read name, and Strand from the request.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    string requestReadName;
    getParameterValue(request, "readName", requestReadName);
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

    // Get the arguments to select what to display.
    string showMarkersString;
    const bool showMarkers = getParameterValue( request, "showMarkers", showMarkersString);
    string showMarkerTableString;
    const bool showMarkerTable = getParameterValue( request, "showMarkerTable", showMarkerTableString);
    string showPositionTableString;
    const bool showPositionTable = getParameterValue( request, "showPositionTable", showPositionTableString);
    string showMarkerFrequencyTableString;
    const bool showMarkerFrequencyTable = getParameterValue( request, "showMarkerFrequencyTable", showMarkerFrequencyTableString);



    // Write the form.
    html <<
        "<form>"
        "<input type=submit value='Show read'> " <<
        "numeric read id " <<
        "<input type=text name=readId" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " size=8 title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"
        " or read name <input type=text name=readName" <<
        (requestReadName.empty() ? "" : " value='" + requestReadName + "'") <<
        "> on strand ";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);
    
    html << "<font color=grey style='font-size:smaller'>";
    html << "&nbsp&nbsp starting at&nbsp";
    html << "<input type=text name=beginPosition size=8 "
        "title='Begin display of raw sequence at this base position (leave blank to begin at beginning of read).'";
    if(beginPositionIsPresent) {
        html << " value=" << beginPosition;
    }
    html << "> &nbsp ending at&nbsp";
    
    html << "<input type=text name=endPosition size=8 "
        "title='End display of raw sequence at this base position (leave blank to end at end of read).'";
    if(endPositionIsPresent) {
        html << " value=" << endPosition;
    }
    html << "> &nbsp with marker&nbsp";
    html << "<input type=text name=highlightMarker size=5 "
        "title='Highlight specified marker string (form accepts only one. URL can be modified to highlight multiple)'";
    if(highlightedMarkerStrings.size() > 0) {
        html << " value=" << highlightedMarkerStrings[0];
    }
    html << "> highlighted.";
    html << "</font>";
    


    // Checkboxes to choose what we want to display.
    html <<
        "<p>Display the following:"
        "<br>&nbsp;<input type=checkbox name=showMarkers" <<
        (showMarkers ? " checked=checked" : "") <<
        "> RLE sequence and markers."
        "<br>&nbsp;<input type=checkbox name=showMarkerTable" <<
        (showMarkerTable ? " checked=checked" : "") <<
        "> Marker table."
        "<br>&nbsp;<input type=checkbox name=showPositionTable" <<
        (showPositionTable ? " checked=checked" : "") <<
        "> Position table."
        "<br>&nbsp;<input type=checkbox name=showMarkerFrequencyTable" <<
        (showMarkerFrequencyTable ? " checked=checked" : "") <<
        "> Marker frequency table."
        "</form>";

    // Check that one and only one of readId and readName was entered.
    if(readIdIsPresent and not requestReadName.empty()) {
        html << "Specify either a numeric read id or a read name, but not both.";
        return;
    }
    if(not readIdIsPresent and requestReadName.empty()) {
        html << "Specify a numeric read id or a read name.";
        return;
    }

    // If a read name was specified, get the read id.
    if(not requestReadName.empty()) {
        readId = getReads().getReadId(requestReadName);
        if(readId == invalidReadId) {
            html << "A read with that name was not found. See ReadSummary.csv.";
            return;
        }
    }



    // If the strand is missing, stop here.
    if(not strandIsPresent) {
        return;
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
    const auto readStoredSequence = reads->getRead(readId);
    const auto readName = reads->getReadName(readId);
    const auto metaData = reads->getReadMetaData(readId);
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];
    if(!beginPositionIsPresent) {
        beginPosition = 0;
    }
    if(!endPositionIsPresent) {
        endPosition = uint32_t(rawOrientedReadSequence.size());
    }else{
        endPosition++; // To include the base at `endPosition`.
    }
    if(endPosition <= beginPosition) {
        html << "<p>Invalid choice of begin and end position.";
        return;
    }



    // Page title.
    html << "<h2 title='Read " << readId << " on strand " << strand;
    if(strand == 0) {
        html << " (input read without reverse complementing)";
    } else {
        html << " (reverse complement of input read)";
    }
    html << "'>Oriented read " << orientedReadId << "</h2>";

    html << "<div style='display:flex;margin-top:10px'>";
    html << "<div style='flex:50%;margin-right:10px'>"; // start column 1
        html << "<div style='background-color:lightgrey;padding:10px'>";
        html << "<font style='font-size:small'>";
        // Read name.
        html << "Read name: ";
        copy(readName.begin(), readName.end(), ostream_iterator<char>(html));

        // Read meta data.
        html << "<br/>Read meta data: ";
        copy(metaData.begin(), metaData.end(), ostream_iterator<char>(html));
        html << "</font>";

        // Read length.
        html << "<br/><br/> This read is <b>" << rawOrientedReadSequence.size() << "</b> bases long";
        html << " (" << readStoredSequence.baseCount << " bases in run-length representation)";
        html << " and has <b>" << orientedReadMarkers.size() << "</b> markers.";

        html << "</div>";
    html << "</div>"; // end column 1
    html << "<div style='flex:50%'>"; // start column 2
        // Button to Blat this read or portion of a read.
        // We cannot use a simple <a> because we need to do a POST
        // (the GET request fails when the read is too long).
        html <<
            "<form style='padding-bottom:8px' action='https://genome.ucsc.edu/cgi-bin/hgBlat' method=post>"
            "<input style='font-size:12px' type=submit value='Blat ";
        if(beginPositionIsPresent || endPositionIsPresent) {
            html << "this portion of ";
        }
        html <<
            "this read in the UCSC browser'>"
            "<input type=text hidden name=type value=DNA>"
            // Don't specify the genome.
            // UCSC browser will Blat again last used genome (stored in cookies).
            // "<input type=text hidden name=type value=DNA>"
            // "<input type=text hidden name=name value=Human>"
            // "<input type=text hidden name=db value=hg38>"
            "<input type=text hidden name=userSeq value=";
        copy(
            rawOrientedReadSequence.begin() + beginPosition,
            rawOrientedReadSequence.begin() + endPosition,
            ostream_iterator<Base>(html));
        html << "></form>";

        // Button to Blast this read or portion of a read.
        html <<
            "<form style='padding-bottom:8px' action='blastRead'>"
            "<input style='font-size:12px' type=submit value='Blast ";
        if(beginPositionIsPresent || endPositionIsPresent) {
            html << "this portion of ";
        }
        html <<
            "this read against " << httpServerData.referenceFastaFileName << " using Blast options: '>"
            "<input type=text hidden name=readId value=" << readId << ">" <<
            "<input type=text hidden name=strand value=" << strand << ">" <<
            "<input type=text hidden name=beginPosition value=" << beginPosition << ">" <<
            "<input type=text hidden name=endPosition value=" << endPosition << ">"
            "<input type=text size=80 name=blastOptions>"
            "</form>";

        // Button to Blast this read or portion of a read (summary output).
        html <<
            "<form action='blastRead'>"
            "<input style='font-size:12px' type=submit value='Blast ";
        if(beginPositionIsPresent || endPositionIsPresent) {
            html << "this portion of ";
        }
        html <<
            "this read against " << httpServerData.referenceFastaFileName << " (summary output)'>"
            "<input type=text hidden name=readId value=" << readId << ">" <<
            "<input type=text hidden name=strand value=" << strand << ">" <<
            "<input type=text hidden name=beginPosition value=" << beginPosition << ">" <<
            "<input type=text hidden name=endPosition value=" << endPosition << ">"
            "<input type=checkbox checked hidden name=summary>"
            "</form>";

        // Link to align this read against another read.
        html <<
            "<button style='background-color:lightsteelblue;font-size:12px;margin-top:16px' " 
            "onclick=\"window.location.href = 'exploreAlignment?readId0=" << readId << "&strand0=" << strand <<
            "';\">Compute a marker alignment of this read with another read</button>";

        html << "<br/>";

        // Link to show overlapping reads.
        html <<
            "<button style='background-color:lightsteelblue;font-size:12px;margin-top:4px' " 
            "onclick=\"window.location.href = 'exploreAlignments?readId=" << readId << "&strand=" << strand <<
            "';\">Find other reads that overlap this read</button>";

    html << "</div>"; // end column 2
    html << "</div>"; // end row

    // Begin/end position (in raw sequence).
    html << "<br/><br/>";
    if(beginPositionIsPresent || endPositionIsPresent) {
        html <<
            " Displaying only " << endPosition-beginPosition << " bases";
        html << " of raw read sequences";
        html << " beginning at base position " << beginPosition <<
            " and ending at base position " << endPosition-1 <<
            " .";
        html << " For sequence in run-length representation see below.";
    }


    // Display the selected portion of raw sequence.
    const bool partialSequenceRequested =  beginPositionIsPresent || endPositionIsPresent;
    if(true) {
        html << "<h3>";
        if(partialSequenceRequested) {
            html << "Selected portion of raw sequence of this oriented read";
        } else {
            html << "Raw sequence of this oriented read";
        }
        html << "</h3>";

        // Here we don't have to worry about using an svg object like we do below,
        // because we are just writing text without html, and so there will
        // be no alignment problems.

        // Labels for position scale.
        html << "<pre style='font-family:monospace;margin:0'";
        html << " title='Position in raw read sequence'";
        html<< ">";
        for(size_t position=beginPosition; position<endPosition; ) {
            if((position%10)==0) {
                const string label = to_string(position);
                html << label;
                for(size_t i=0; i<10-label.size(); i++) {
                    html << " ";
                }
                position += 10;
            } else {
                html << " ";
                ++position;
            }
        }
        html<< "\n";

        // Position scale.
        for(size_t position=beginPosition; position<endPosition; position++) {
            if((position%10)==0) {
                html << "|";
            } else if((position%5)==0) {
                html << "+";
            } else {
                html << ".";
            }
        }
        html << "</pre>";



        // Sequence.
        html << "<pre style='font-family:monospace;margin:0'>";
        for(uint32_t position=beginPosition; position!=endPosition; position++) {
            html << rawOrientedReadSequence[position];
        }
        html << "</pre>";



        // Also write a position scale for positions in the run-length representation.
        if(true) {
            html << "<pre style='font-family:monospace;margin:0'";
            html << " title='Position in run-length read sequence'>";

            const vector<uint32_t> rawPositions = reads->getRawPositions(orientedReadId);

            // Scale.
            bool firstTime = true;
            for(int runLengthPosition=0; runLengthPosition<int(rawPositions.size()); runLengthPosition++) {
                const int rawPosition = rawPositions[runLengthPosition];
                // cout << runLengthPosition << " " << rawPosition << endl;
                if(rawPosition >= int(endPosition)) {
                    break;
                }
                if(rawPosition < int(beginPosition)) {
                    continue;
                }
                uint32_t skip;
                if(firstTime) {
                    skip = rawPosition - beginPosition;
                } else {
                    skip = rawPosition - rawPositions[runLengthPosition-1] - 1;
                }
                for(uint32_t i=0; i<skip; i++) {
                    html << " ";
                }
                firstTime = false;
                if((runLengthPosition % 10) == 0) {
                    html << "|";
                    //cout << "|";
                } else if((runLengthPosition % 5) == 0) {
                    html << "+";
                    //cout << "+";
                } else {
                    html << ".";
                    //cout << ".";
                }
            }
            html << "\n";

            // Labels.
            firstTime = true;
            for(int runLengthPosition=0; runLengthPosition<int(rawPositions.size()); runLengthPosition+=10) {
                const int rawPosition = rawPositions[runLengthPosition];
                if(rawPosition >= int(endPosition)) {
                    break;
                }
                if(rawPosition < int(beginPosition)) {
                    continue;
                }

                uint32_t skip;
                if(firstTime) {
                    skip = rawPosition - beginPosition;
                } else {
                    skip = rawPosition - rawPositions[runLengthPosition-10] - 10;
                }
                for(uint32_t i=0; i<skip; i++) {
                    html << " ";
                }
                firstTime = false;

                const string label = to_string(runLengthPosition);
                html << label;
                for(size_t i=0; i<10-label.size(); i++) {
                    html << " ";
                }

            }


            html << "</pre>";
            //cout << endl;
        }


        // Button to download the sequence to a fasta file
        html <<
            "<a id=fastaDownload>Download in FASTA format</a><br>"
            "<script>"
            "var element = document.getElementById('fastaDownload');"
            "element.setAttribute('href', 'data:text/plain;charset=utf-8,' +"
            "encodeURIComponent('>" << orientedReadId <<
            "-" << beginPosition << "-" << endPosition << " " << endPosition-beginPosition <<
            " ";
        copy(readName.begin(), readName.end(), ostream_iterator<char>(html));
        html << "\\n";
        for(uint32_t position=beginPosition; position!=endPosition; position++) {
            html << rawOrientedReadSequence[position];
        }
        html << "\\n'));"
            "element.setAttribute('download', '" << orientedReadId << "-" <<
            beginPosition << "-" << endPosition <<
            ".fa');"
            "</script>";
    }



    // If there are no markers, stop here.
    if(orientedReadMarkers.empty()) {
        html << "<p>This read has no markers.";
        return;
    }



    // Compute begin/end in RLE coordinates.
    size_t beginRlePosition, endRlePosition;
    const vector<uint32_t> rawPositions = reads->getRawPositions(orientedReadId);

    if (beginPositionIsPresent) {
        beginRlePosition = std::lower_bound(rawPositions.begin(), rawPositions.end(), beginPosition) - rawPositions.begin();
        while(beginRlePosition > 0 && beginRlePosition % 10 != 0) {
            // Start with a multiple of 10 so that labeling positions is easier.
            beginRlePosition--;
        }
    } else {
        beginRlePosition = beginPosition;
    }

    if (endPositionIsPresent) {
        endRlePosition = std::lower_bound(rawPositions.begin(), rawPositions.end(), endPosition) - rawPositions.begin();
    } else{
        endRlePosition = size_t(readStoredSequence.baseCount);
    }



    // Display RLE sequence and markers.
    const size_t k = assemblerInfo->k;
    if(showMarkers) {

        // Decide on which row each marker gets displayed
        // (first row of markers is row 0).
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
            SHASTA_ASSERT(markerRow[ordinal] != -1);
        }
        const int markerRowCount = *std::max_element(markerRow.begin(), markerRow.end());


        // Title for the next portion of the display, which shows the markers.
        if (beginPositionIsPresent || endPositionIsPresent) {
            html << "<h3>Selected portion of run-length representation of oriented read sequence and its markers</h3>";
        } else {
            html << "<h3>Run-length representation of oriented read sequence and its markers</h3>";
        }

        // Use an svg object to display the read sequence as stored and the markers.
        // To ensure correct positioning and alignment, we use
        // a textLength attribute on every <text> element.
        // (The older code, ifdef'ed out, uses a separate <text>
        // element for each character).

        const int monospaceFontSize = 12;
        const int horizontalSpacing = 7;
        const int verticalSpacing = 13;
        const size_t charactersPerLine = endRlePosition - beginRlePosition + 10; // Add space for labels
        int svgLineCount = int(3 + markerRowCount); // Labels, scale, sequence, markers.
        svgLineCount++;     // Add a line with the repeat counts.
        const size_t svgWidth = horizontalSpacing * charactersPerLine;
        const size_t svgHeight = verticalSpacing * svgLineCount;
        const int highlightedMarkerVerticalOffset = 2;
        html <<
            "<p><svg width=" << svgWidth << " height=" << svgHeight << ">"
            "<style>"
            ".mono{font-family:monospace; font-size:" << monospaceFontSize << "px;}"
            ".blueMono{font-family:monospace; font-size:" << monospaceFontSize << "px; fill:blue;}"
            "</style>";


        // Labels for position scale.
        for(size_t position=beginRlePosition; position<endRlePosition; position+=10) {
            const string label = to_string(position);

            // Use a single <text> element with a textLength attribute for exact alignment.
            html <<
                "<text class='mono'" <<
                " x='" << (position-beginRlePosition) * horizontalSpacing << "'" <<
                " y='" << verticalSpacing << "'" <<
                " textLength='" << label.size() * horizontalSpacing << "px'>" <<
                label << "</text>";
        }


        // Position scale.
        // This code uses one <text> element for every blockSize characters.
        // This way you can select sequence text without getting a
        // new line after each character, while still achieving good
        // alignment.
        const uint64_t blockSize = 100;
        for(uint64_t blockBegin=0; blockBegin<(endRlePosition-beginRlePosition); blockBegin+=blockSize) {
            const uint64_t blockEnd = min(
                uint64_t(blockBegin+blockSize),
                uint64_t(endRlePosition-beginRlePosition)
            );
            html <<
                "<text class='mono'" <<
                " x='" << blockBegin*horizontalSpacing << "'" <<
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



        // Repeat counts.
        for(size_t position=beginRlePosition; position!=endRlePosition; position++) {
            Base base;
            uint8_t repeatCount;
            tie(base, repeatCount) = reads->getOrientedReadBaseAndRepeatCount(orientedReadId, uint32_t(position));

            html <<
                "<text class='mono'" <<
                " x='" << (position-beginRlePosition)*horizontalSpacing << "'" <<
                " y='" << 3*verticalSpacing << "'"
                " textLength='" << horizontalSpacing<< "'>";
            if(repeatCount < 10) {
                html << int(repeatCount);
            } else {
                html << "*";
            }
            html << "<title>" << base << " at run-length position " << position <<
                " is repeated " << int(repeatCount) << " times</title>";
            html << "</text>";
        }

        
        // Read sequence in run length encoding.
        // This code uses one <text> element for every blockSize characters.
        // This way you can select sequence text without getting a
        // new line after each character, while still achieving good
        // alignment.
        const uint32_t readSequenceLine = 4;
        for(uint64_t blockBegin=0; blockBegin<(endRlePosition-beginRlePosition); blockBegin+=blockSize) {
            const uint64_t blockEnd = min(
                uint64_t(blockBegin+blockSize),
                uint64_t(endRlePosition-beginRlePosition)
            );
            html <<
                "<text class='mono'" <<
                " x='" << blockBegin*horizontalSpacing << "'" <<
                " y='" << readSequenceLine*verticalSpacing << "'"
                " textLength='" << (blockEnd-blockBegin) * horizontalSpacing<< "'>";
            for(size_t position=beginRlePosition+blockBegin; position<beginRlePosition+blockEnd; position++) {
                html << reads->getOrientedReadBase(orientedReadId, uint32_t(position));
            }
            html << "</text>";
        }
    


        // Draw a rectangle for each highlighted marker.
        for(const uint32_t ordinal: highlightedMarkers) {
            const CompressedMarker& marker = orientedReadMarkers[ordinal];
            html <<
                "<rect" <<
                " x='" << (marker.position-beginRlePosition)*horizontalSpacing << "'"
                " y='" << (readSequenceLine + markerRow[ordinal])*verticalSpacing + highlightedMarkerVerticalOffset << "'"
                " height='" << verticalSpacing << "'"
                " width='" << k * horizontalSpacing << "'"
                " style='fill:pink; stroke:none;'"
                "/>";
        }



        // Markers.
        if(markers.isOpen() and markerGraph.vertices().isOpen()) {
            for(uint32_t ordinal=0; ordinal<uint32_t(orientedReadMarkers.size()); ordinal++) {
                const CompressedMarker& marker = orientedReadMarkers[ordinal];
                if (marker.position < beginRlePosition || marker.position > endRlePosition-k) {
                    continue;
                }

                // See if this marker is contained in a vertex of the marker graph.
                const MarkerGraph::VertexId vertexId =
                    getGlobalMarkerGraphVertex(orientedReadId, ordinal);
                const bool hasMarkerGraphVertex =
                    (vertexId != MarkerGraph::invalidCompressedVertexId);



                // Write the k-mer of this marker.
                const Kmer kmer(marker.kmerId, k);
                html << "<a xlink:title='Marker " << ordinal <<
                    ", position " << marker.position <<
                    ", k-mer id " << marker.kmerId;
                if(hasMarkerGraphVertex) {
                    html << ", coverage " << markerGraph.vertexCoverage(vertexId);
                }
                html << "' id='" << ordinal << "'";
                if(hasMarkerGraphVertex) {
                    // Add a hyperlink to the marker graph vertex
                    // that contains this marker.
                    const string url = "exploreMarkerGraph?vertexId=" + to_string(vertexId) +
                        "&maxDistance=2&detailed=on&minCoverage=3&minConsensus=3&sizePixels=320&timeout=30";
                    html << " xlink:href='" << url << "' style='cursor:pointer'";
                }
                html << ">";

                // This code uses one <text> element per character.
                for(size_t positionInMarker=0; positionInMarker<k; positionInMarker++) {
                    html << "<text class='";
                    if(hasMarkerGraphVertex) {
                        html << "blueMono";
                    } else {
                        html << "mono";
                    }
                    html << "'" <<
                        " x='" << (marker.position-beginRlePosition+positionInMarker)*horizontalSpacing << "'" <<
                        " y='" << (readSequenceLine+1+markerRow[ordinal])*verticalSpacing << "'>";
                    html << kmer[positionInMarker];
                    html << "</text>";
                }
                html << "</a>";

            }
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



        html <<
            "<p>You can click on a blue marker above "
            "to see the global marker graph around that marker. "
            "Black markers correspond to a vertex of the marker graph "
            "that was removed because of low coverage.";
    }



    // Position table.
    // It can be used to convert RLE positions to raw positions.
    if(showPositionTable) {
        html << "<h3>Position table</h3>"
            "<table><tr><th>RLE<br>position<th>Raw<br>position\n";
        for(uint32_t rlePosition=0; rlePosition<rawPositions.size(); rlePosition++) {
            const uint32_t rawPosition = rawPositions[rlePosition];
            if(beginPositionIsPresent and rawPosition<beginPosition) {
                continue;
            }
            if(endPositionIsPresent and rawPosition>endPosition) {
                continue;
            }
            html << "<tr><td class=centered>" << rlePosition <<
                "<td class=centered>" << rawPositions[rlePosition] << "\n";
        }
        html << "</table>\n";
    }



    // Details of markers.
    if(showMarkerTable) {
        if (beginPositionIsPresent || endPositionIsPresent) {
            html << "<h3>Detail of markers on this portion of this oriented read</h3>";
        } else {
            html << "<h3>Detail of markers on this oriented read</h3>";
        }
        html << "<table><tr><th>Marker<br>ordinal<th>k-mer"
            "<th>RLE<br>position<th>Raw<br>position"
            "<th>Marker<br>graph<br>vertex"
            "<th>Marker<br>graph<br>vertex<br>coverage";

        // Loop over all markers on this oriented read.
        for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
            const CompressedMarker& marker = orientedReadMarkers[ordinal];
            const Kmer kmer(marker.kmerId, k);
            const uint32_t rlePosition = marker.position;
            const uint32_t rawPosition = rawPositions[rlePosition];

            if(beginPositionIsPresent and rawPosition<beginPosition) {
                continue;
            }
            if(endPositionIsPresent and rawPosition>endPosition) {
                continue;
            }

            const MarkerGraph::VertexId vertexId =
                getGlobalMarkerGraphVertex(orientedReadId, ordinal);
            const bool hasMarkerGraphVertex =
                (vertexId != MarkerGraph::invalidCompressedVertexId);

            html << "<tr><td class=centered>"<< ordinal << "<td class=centered><code>";
            kmer.write(html, k);
            html << "</code><td class=centered>" << rlePosition << "<td class=centered>" << rawPosition;

            if(hasMarkerGraphVertex) {
                const string url = "exploreMarkerGraph?vertexId=" + to_string(vertexId) +
                    "&maxDistance=2&detailed=on&minCoverage=3&minConsensus=3&sizePixels=320&timeout=30";
                html << "<td class=centered><a href='" << url << "'>" << vertexId << "</a>"
                    "<td class=centered>" << markerGraph.vertexCoverage(vertexId);
            } else {
                html << "<td><td>";
            }


            html << "\n";
        }

        html << "</table>";
    }



    // Frequency of markers on this oriented read.
    if(showMarkerFrequencyTable) {
        vector<KmerId> kmers;
        for(uint32_t ordinal=0; ordinal<uint32_t(orientedReadMarkers.size()); ordinal++) {
            const CompressedMarker& marker = orientedReadMarkers[ordinal];
            if (marker.position >= beginRlePosition && marker.position <= endRlePosition - k) {
                kmers.push_back(marker.kmerId);
            }
        }
        vector<uint32_t> kmerFrequency;
        deduplicateAndCount(kmers, kmerFrequency);
        vector< pair<KmerId, uint32_t> > markerFrequencyTable;
        for(uint32_t i=0; i<kmers.size(); i++) {
            markerFrequencyTable.push_back(make_pair(kmers[i], kmerFrequency[i]));
        }
        sort(markerFrequencyTable.begin(), markerFrequencyTable.end(),
            OrderPairsBySecondOnlyGreater<KmerId, uint32_t>());

        if(beginPositionIsPresent || endPositionIsPresent) {
            html << "<h2>Frequency of markers in the selected portion in this oriented read</h2>";
        } else {
            html << "<h2>Frequency of markers in this oriented read</h2>";
        }
        html << "<table><tr><th>Marker<th>Frequency";
        for(const auto& p: markerFrequencyTable) {
            const KmerId kmerId = p.first;
            const uint32_t frequency = p.second;
            const Kmer kmer(kmerId, k);
            html << "<tr><td>";
            kmer.write(html, k);
            html << "<td>" << frequency;
        }
        html << "<table>";
    }


}



void Assembler::exploreReadRaw(
    const vector<string>& request,
    ostream& html)
{
    SHASTA_ASSERT(assemblerInfo->readRepresentation == 0);

    // Get the request parameters.
    ReadId readId = 0;
    const bool readIdIsPresent = getParameterValue(request, "readId", readId);
    string requestReadName;
    getParameterValue(request, "readName", requestReadName);
    Strand strand = 0;
    const bool strandIsPresent = getParameterValue(request, "strand", strand);
    uint32_t beginPosition = 0;
    const bool beginPositionIsPresent = getParameterValue(request, "beginPosition", beginPosition);
    uint32_t endPosition = 0;
    const bool endPositionIsPresent = getParameterValue(request, "endPosition", endPosition);


    // Write the form.
    html <<
        "<form>"
        "<table>"

        "<tr>"
        "<th class=left>Numeric read id"
        "<td><input type=text name=readId" <<
        (readIdIsPresent ? (" value=" + to_string(readId)) : "") <<
        " title='Enter a read id between 0 and " << reads->readCount()-1 << "'>"

        "<tr>"
        "<th class=left>Read name"
        "<td><input type=text name=readName" <<
        (requestReadName.empty() ? "" : " value='" + requestReadName + "'") << ">"

        "<tr>"
        "<th class=left>Strand"
        "<td>";
    writeStrandSelection(html, "strand", strandIsPresent && strand==0, strandIsPresent && strand==1);

    html <<
        "<tr>"
        "<th class=left>Begin position"
        "<td><input type=text name=beginPosition"
        " title='Leave blank to begin display at beginning of read.'";
    if(beginPositionIsPresent) {
        html << " value=" << beginPosition;
    }
    html << ">";

    html <<
        "<tr>"
        "<th class=left>End position"
        "<td><input type=text name=endPosition"
        " title='Leave blank to end display at end of read.'";
    if(endPositionIsPresent) {
        html << " value=" << endPosition;
    }
    html << ">";


    html <<
        "</table>"
        "<input type=submit value='Display'>"
        "</form>";



    // Check that one and only one of readId and readName was entered.
    if(readIdIsPresent and not requestReadName.empty()) {
        html << "Specify either a numeric read id or a read name, but not both.";
        return;
    }
    if(not readIdIsPresent and requestReadName.empty()) {
        html << "Specify a numeric read id or a read name.";
        return;
    }

    // If a read name was specified, get the read id.
    if(not requestReadName.empty()) {
        readId = getReads().getReadId(requestReadName);
        if(readId == invalidReadId) {
            html << "A read with that name was not found. See ReadSummary.csv.";
            return;
        }
    }

    // If the strand is missing, stop here.
    if(not strandIsPresent) {
        return;
    }

    // Sanity checks.
    if(readId >= reads->readCount()) {
        html << "<p>Invalid read id.";
        return;
    }
    if(strand!=0 && strand!=1) {
        html << "<p>Invalid strand.";
        return;
    }


    // Access the read information we need.
    const OrientedReadId orientedReadId(readId, strand);
    const auto sequence = reads->getRead(readId);
    const auto readName = reads->getReadName(readId);
    const auto metaData = reads->getReadMetaData(readId);

    // Adjust the position range, if necessary.
    if(!beginPositionIsPresent) {
        beginPosition = 0;
    }
    if(!endPositionIsPresent) {
        endPosition = uint32_t(sequence.baseCount);
    } else {
        endPosition++; // To include the base at `endPosition`.
    }
    if(endPosition <= beginPosition) {
        html << "<p>Invalid choice of begin and end position.";
        return;
    }



    // Page title.
    html << "<h2 title='Read " << readId << " on strand " << strand;
    if(strand == 0) {
        html << " (input read without reverse complementing)";
    } else {
        html << " (reverse complement of input read)";
    }
    html << "'>Oriented read " << orientedReadId << "</h2>";



    // Read information.
    html << "<table>";

    html << "<tr><th class=left>Read id<td>" << readId;

    html << "<tr><th class=left>Read name<td>";
    copy(readName.begin(), readName.end(), ostream_iterator<char>(html));

    html << "<tr><th class=left>Read meta data<td>";
    copy(metaData.begin(), metaData.end(), ostream_iterator<char>(html));

    html << "<tr><th class=left>Length<td>" << sequence.baseCount;

    html << "<tr><th class=left>Length displayed<td>" << endPosition - beginPosition;

    html << "</table>";



    // Sequence.
    html << "<p><pre style='font-family:monospace;margin:0'>";
    for(uint64_t position=beginPosition; position!=endPosition; position++) {
        html << sequence[position];
    }
    html<< "\n";

    // Position scale
    for(size_t position=beginPosition; position<endPosition; position++) {
        if((position%10)==0) {
            html << "|";
        } else if((position%5)==0) {
            html << "+";
        } else {
            html << ".";
        }
    }
    html<< "\n";

    // Position scale labels.
    for(size_t position=beginPosition; position<endPosition; ) {
        if((position%10)==0) {
            const string label = to_string(position);
            html << label;
            for(size_t i=0; i<10-label.size(); i++) {
                html << " ";
            }
            position += 10;
        } else {
            html << " ";
            ++position;
        }
    }
    html<< "\n";
    html << "</pre>";

}


#endif
