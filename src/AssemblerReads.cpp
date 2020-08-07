// Shasta.
#include "Assembler.hpp"
#include "ReadLoader.hpp"
using namespace shasta;

// Standard libraries.
#include "algorithm.hpp"
#include "iterator.hpp"


// Add reads.
// The reads are added to those already previously present.
void Assembler::addReads(
    const string& fileName,
    uint64_t minReadLength,
    uint64_t desiredCoverage,
    bool noCache,
    const size_t threadCount)
{
    reads.checkReadsAreOpen();
    reads.checkReadNamesAreOpen();

    ReadLoader readLoader(
        fileName,
        minReadLength,
        desiredCoverage,
        noCache,
        threadCount,
        largeDataFileNamePrefix,
        largeDataPageSize,
        reads);

    reads.checkSanity();

    cout << "Discarded read statistics for file " << fileName << ":" << endl;
    cout << "    Discarded " << readLoader.discardedInvalidBaseReadCount <<
        " reads containing invalid bases for a total " <<
        readLoader.discardedInvalidBaseBaseCount << " valid bases." << endl;
    cout << "    Discarded " << readLoader.discardedShortReadReadCount <<
        " reads shorter than " << minReadLength <<
        " bases for a total " << readLoader.discardedShortReadBaseCount << " bases." << endl;
    cout << "    Discarded " << readLoader.discardedBadRepeatCountReadCount <<
        " reads containing repeat counts 256 or more" <<
        " for a total " << readLoader.discardedBadRepeatCountBaseCount << " bases." << endl;

    // Increment the discarded reads statistics.
    assemblerInfo->discardedInvalidBaseReadCount += readLoader.discardedInvalidBaseReadCount;
    assemblerInfo->discardedInvalidBaseBaseCount += readLoader.discardedInvalidBaseBaseCount;
    assemblerInfo->discardedShortReadReadCount += readLoader.discardedShortReadReadCount;
    assemblerInfo->discardedShortReadBaseCount += readLoader.discardedShortReadBaseCount;
    assemblerInfo->discardedBadRepeatCountReadCount += readLoader.discardedBadRepeatCountReadCount;
    assemblerInfo->discardedBadRepeatCountBaseCount += readLoader.discardedBadRepeatCountBaseCount;
}


// Create a histogram of read lengths.
// All lengths here are raw sequence lengths
// (length of the original read), not lengths
// in run-length representation.
void Assembler::histogramReadLength(const string& fileName)
{
    reads.writeReadLengthHistogram(fileName);

    cout << "Discarded read statistics for all input files:" << endl;;
    cout << "    Discarded " << assemblerInfo->discardedInvalidBaseReadCount <<
        " reads containing invalid bases for a total " <<
        assemblerInfo->discardedInvalidBaseBaseCount << " valid bases." << endl;
    cout << "    Discarded " << assemblerInfo->discardedShortReadReadCount <<
        " short reads for a total " <<
        assemblerInfo->discardedShortReadBaseCount << " bases." << endl;
    cout << "    Discarded " << assemblerInfo->discardedBadRepeatCountReadCount <<
        " reads containing repeat counts 256 or more" <<
        " for a total " << assemblerInfo->discardedBadRepeatCountBaseCount << " bases." << endl;

    cout << "Read statistics for reads that will be used in this assembly:" << endl;
    cout << "    Total number of reads is " << reads.readCount() << "." << endl;
    cout << "    Total number of raw bases is " << reads.getTotalBaseCount() << "." << endl;
    cout << "    Average read length is " << double(reads.getTotalBaseCount()) / double(reads.readCount());
    cout << " bases." << endl;
    cout << "    N50 for read length is " << reads.getN50() << " bases." << endl;
    cout << "    The above statistics only include reads that will be used in this assembly." << endl;
    cout << "    Read discarded because they contained invalid bases, were too short or contained repeat counts 256"
        " or more are not counted." << endl;

    // Store read statistics in AssemblerInfo.
    assemblerInfo->readCount = reads.readCount();
    assemblerInfo->baseCount = reads.getTotalBaseCount();
    assemblerInfo->readN50 = reads.getN50();
}


// Write a csv file with summary information for each read.
void Assembler::writeReadsSummary()
{
    reads.checkReadsAreOpen();
    reads.checkReadNamesAreOpen();
    SHASTA_ASSERT(markers.isOpen());

    // Count the number of alignment candidates for each read.
    vector<uint64_t> alignmentCandidatesCount(reads.readCount(), 0);
    for(const OrientedReadPair& p: alignmentCandidates.candidates) {
        ++alignmentCandidatesCount[p.readIds[0]];
        ++alignmentCandidatesCount[p.readIds[1]];
    }


    ofstream csv("ReadSummary.csv");
    csv << "Id,Name,RawLength,RleLength,RawOverRleLengthRatio,"
        "MarkerCount,MarkerDensity,MaximumMarkerOffset,"
        "Palindromic,Chimeric,"
        "AlignmentCandidates,ReadGraphNeighbors,"
        "VertexCount,VertexDensity,runid,sampleid,read,ch,start_time,\n";
    for(ReadId readId=0; readId!=reads.readCount(); readId++) {
        const OrientedReadId orientedReadId(readId, 0);

        // Read id.
        csv << readId << ",";

        // Read name.
        const auto readName = reads.getReadName(readId);
        copy(readName.begin(), readName.end(), ostream_iterator<char>(csv));
        csv << ",";

        // Number of raw bases.
        const auto repeatCounts = reads.getReadRepeatCounts(readId);
        uint64_t rawBaseCount = 0;
        for(const auto repeatCount: repeatCounts) {
            rawBaseCount += repeatCount;
        }
        csv << rawBaseCount << ",";

        // Number of RLE bases.
        const uint64_t rleBaseCount = reads.getRead(readId).baseCount;
        csv << rleBaseCount << ",";

        // Ratio of raw over RLE base count.
        csv << double(rawBaseCount)/double(rleBaseCount) << ",";

        // Number of markers.
        const uint64_t markerCount = markers.size(orientedReadId.getValue());
        csv << markerCount << ",";

        // Marker density.
        const double markerDensity = double(markerCount) / double(rleBaseCount);
        csv << markerDensity << ",";

        // Maximum marker offset (offset between consecutive markers).
        const auto readMarkers = markers[orientedReadId.getValue()];
        uint64_t position = 0;
        uint64_t maximumMarkerOffset = 0;
        for(const auto marker: readMarkers) {
            const uint64_t offset = marker.position - position;
            maximumMarkerOffset = max(maximumMarkerOffset, offset);
            position = marker.position;
        }
        maximumMarkerOffset = max(maximumMarkerOffset, rleBaseCount-position);
        csv << maximumMarkerOffset << ",";

        // Palindromic flag.
        csv << (reads.getFlags(readId).isPalindromic ? "Yes" : "No") << ",";

        // Chimeric flag.
        csv << (reads.getFlags(readId).isChimeric ? "Yes" : "No") << ",";

        // Alignment candidates
        csv << alignmentCandidatesCount[readId] << ",";

        // Number of read graph neighbors.
        if(directedReadGraph.isOpen()) {
            csv << directedReadGraph.totalDegree(orientedReadId.getValue()) << ",";
        } else {
            csv << readGraph.connectivity.size(orientedReadId.getValue()) << ",";
        }

        // Number and fraction of markers associated with a marker graph vertex.
        uint64_t vertexCount = 0;
        for(uint32_t ordinal=0; ordinal<readMarkers.size(); ordinal++) {
            const MarkerId markerId = getMarkerId(orientedReadId, ordinal);
            const MarkerGraph::VertexId vertexId = markerGraph.vertexTable[markerId];
            if(vertexId != MarkerGraph::invalidCompressedVertexId) {
                ++vertexCount;
            }
        }
        csv << vertexCount << ",";
        csv << double(vertexCount) / double(markerCount) << ",";

        // Oxford Nanopore standard metadata.
        csv << reads.getMetaData(readId, "runid") << ",";
        csv << reads.getMetaData(readId, "sampleid") << ",";
        csv << reads.getMetaData(readId, "read") << ",";
        csv << reads.getMetaData(readId, "ch") << ",";
        csv << reads.getMetaData(readId, "start_time") << ",";

        // End of the line for this read.
        csv << "\n";
     }
}

