// Shasta.
#include "Assembler.hpp"
#include "ReadLoader.hpp"
using namespace shasta;

// Standard libraries.
#include "iterator.hpp"



void Assembler::checkReadsAreOpen() const
{
    if(!reads.isOpen()) {
        throw runtime_error("Reads are not accessible.");
    }
    if(!readRepeatCounts.isOpen()) {
        throw runtime_error("Read repeat counts are not accessible.");
    }
}
void Assembler::checkReadNamesAreOpen() const
{
    if(!readNames.isOpen()) {
        throw runtime_error("Read names are not accessible.");
    }
}
void Assembler::checkReadMetaDataAreOpen() const
{
    if(!readMetaData.isOpen()) {
        throw runtime_error("Read metadata are not accessible.");
    }
}
void Assembler::checkReadId(ReadId readId) const
{
    if(readId >= reads.size()) {
        throw runtime_error("Read id " + to_string(readId) +
            " is not valid. Must be between 0 and " + to_string(reads.size()) +
            " inclusive.");
    }
}



// Add reads.
// The reads are added to those already previously present.
void Assembler::addReads(
    const string& fileName,
    size_t minReadLength,
    const size_t threadCount)
{
    checkReadsAreOpen();
    checkReadNamesAreOpen();

    ReadLoader readLoader(
        fileName,
        minReadLength,
        threadCount,
        largeDataFileNamePrefix,
        largeDataPageSize,
        reads,
        readNames,
        readMetaData,
        readRepeatCounts);

    // Sanity checks.
    SHASTA_ASSERT(readNames.size() == reads.size());
    SHASTA_ASSERT(readMetaData.size() == reads.size());

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
    // Check that we have what we need.
    checkReadsAreOpen();

    // Create the histogram.
    // It contains the number of reads of each length.
    // Indexed by the length.
    vector<size_t> histogram;
    const ReadId totalReadCount = readCount();
    size_t totalBaseCount = 0;
    for(ReadId readId=0; readId<totalReadCount; readId++) {
        const size_t length = getReadRawSequenceLength(readId);
        totalBaseCount += length;
        if(histogram.size() <= length) {
            histogram.resize(length+1, 0);
        }
        ++(histogram[length]);
    }

    // Write it out.
    size_t n50 = 0;
    {
        ofstream csv(fileName);
        csv << "Length,Reads,Bases,CumulativeReads,CumulativeBases,"
            "FractionalCumulativeReads,FractionalCumulativeBases,\n";
        size_t cumulativeReadCount = totalReadCount;
        size_t cumulativeBaseCount = totalBaseCount;
        for(size_t length=0; length<histogram.size(); length++) {
            const size_t frequency = histogram[length];
            if(frequency) {
                const  size_t baseCount = frequency * length;
                const double cumulativeReadFraction =
                    double(cumulativeReadCount)/double(totalReadCount);
                const double comulativeBaseFraction =
                    double(cumulativeBaseCount)/double(totalBaseCount);
                csv << length << "," << frequency << "," << baseCount << ",";
                csv << cumulativeReadCount << "," << cumulativeBaseCount << ",";
                csv << cumulativeReadFraction << ",";
                csv << comulativeBaseFraction << "\n";
                cumulativeReadCount -= frequency;
                cumulativeBaseCount -= baseCount;
                if(comulativeBaseFraction > 0.5) {
                    n50 = length;
                }
            }
        }
        SHASTA_ASSERT(cumulativeReadCount == 0);
        SHASTA_ASSERT(cumulativeBaseCount == 0);

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
        cout << "    Total number of reads is " << totalReadCount << "." << endl;
        cout << "    Total number of raw bases is " << totalBaseCount << "." << endl;
        cout << "    Average read length is " << double(totalBaseCount) / double(totalReadCount);
        cout << " bases." << endl;
        cout << "    N50 for read length is " << n50 << " bases." << endl;
        cout << "    The above statistics only include reads that will be used in this assembly." << endl;
        cout << "    Read discarded because they contained invalid bases, were too short or contained repeat counts 256"
            " or more are not counted." << endl;

        // Store read statistics in AssemblerInfo.
        assemblerInfo->readCount = totalReadCount;
        assemblerInfo->baseCount = totalBaseCount;
        assemblerInfo->readN50 = n50;

    }



    // Also write out a histogram of number of reads in 1 Kb bins.
    {
        // Each entry of the binned histogram contains pair(read count, bases)
        // for that bin.
        const size_t binWidth = 1000;
        vector<pair <size_t, size_t> > binnedHistogram;

        for(size_t length=0; length<histogram.size(); length++) {
            const size_t readCount = histogram[length];
            if(readCount) {
                const size_t bin = length / binWidth;
                if(binnedHistogram.size() <= bin) {
                    binnedHistogram.resize(bin+1, make_pair(0, 0));
                }
                binnedHistogram[bin].first += readCount;
                binnedHistogram[bin].second += readCount * length;
            }
        }

        ofstream csv("Binned-" + fileName);
        csv << "LengthBegin,LengthEnd,Reads,Bases,CumulativeReads,CumulativeBases,"
            "FractionalCumulativeReads,FractionalCumulativeBases,\n";
        size_t cumulativeReadCount = totalReadCount;
        size_t cumulativeBaseCount = totalBaseCount;
        for(size_t bin=0; bin<binnedHistogram.size(); bin++) {
            const auto& histogramBin = binnedHistogram[bin];
            const size_t readCount = histogramBin.first;
            const size_t baseCount = histogramBin.second;
            const double cumulativeReadFraction =
                double(cumulativeReadCount)/double(totalReadCount);
            const double comulativeBaseFraction =
                double(cumulativeBaseCount)/double(totalBaseCount);
            csv << bin*binWidth << ",";
            csv << (bin+1)*binWidth << ",";
            csv << readCount << "," << baseCount << ",";
            csv << cumulativeReadCount << "," << cumulativeBaseCount << ",";
            csv << cumulativeReadFraction << ",";
            csv << comulativeBaseFraction << "\n";
            cumulativeReadCount -= readCount;
            cumulativeBaseCount -= baseCount;
        }
        SHASTA_ASSERT(cumulativeReadCount == 0);
        SHASTA_ASSERT(cumulativeBaseCount == 0);
    }

    cout << "See " << fileName << " and Binned-" << fileName <<
        " for details of the read length distribution." << endl;

}



// Function to write one or all reads in Fasta format.
void Assembler::writeReads(const string& fileName)
{
    ofstream file(fileName);
    for(ReadId readId=0; readId<readCount(); readId++) {
        writeRead(readId, file);
    }
}



void Assembler::writeRead(ReadId readId, const string& fileName)
{
    ofstream file(fileName);
    writeRead(readId, file);
}


void Assembler::writeRead(ReadId readId, ostream& file)
{
    checkReadsAreOpen();
    checkReadNamesAreOpen();
    checkReadId(readId);

    const vector<Base> rawSequence = getOrientedReadRawSequence(OrientedReadId(readId, 0));
    const auto readName = readNames[readId];
    const auto metaData = readMetaData[readId];

    file << ">";
    copy(readName.begin(), readName.end(), ostream_iterator<char>(file));
    file << " " << readId;
    file << " " << rawSequence.size();
    if(metaData.size() > 0) {
        file << " ";
        copy(metaData.begin(), metaData.end(), ostream_iterator<char>(file));
    }
    file << "\n";
    copy(rawSequence.begin(), rawSequence.end(), ostream_iterator<Base>(file));
    file << "\n";

}



void Assembler::writeOrientedRead(ReadId readId, Strand strand, const string& fileName)
{
    writeOrientedRead(OrientedReadId(readId, strand), fileName);
}



void Assembler::writeOrientedRead(OrientedReadId orientedReadId, const string& fileName)
{
    ofstream file(fileName);
    writeOrientedRead(orientedReadId, file);
}



void Assembler::writeOrientedRead(OrientedReadId orientedReadId, ostream& file)
{
    checkReadsAreOpen();
    checkReadNamesAreOpen();

    const vector<Base> rawSequence = getOrientedReadRawSequence(orientedReadId);
    const auto readName = readNames[orientedReadId.getReadId()];

    file << ">" << orientedReadId;
    file << " " << rawSequence.size() << " ";
    copy(readName.begin(), readName.end(), ostream_iterator<char>(file));
    file << "\n";
    copy(rawSequence.begin(), rawSequence.end(), ostream_iterator<Base>(file));
    file << "\n";

}



// Return a vector containing the raw sequence of an oriented read.
vector<Base> Assembler::getOrientedReadRawSequence(OrientedReadId orientedReadId)
{
    // The sequence we will return;
    vector<Base> sequence;

    // The number of bases stored, in run-length representation.
    const uint32_t storedBaseCount = uint32_t(reads[orientedReadId.getReadId()].baseCount);


    // We are storing a run-length representation of the read.
    // Expand it base by base to create the raw representation.
    for(uint32_t position=0; position<storedBaseCount; position++) {
        Base base;
        uint8_t count;
        tie(base, count) = getOrientedReadBaseAndRepeatCount(orientedReadId, position);
        for(uint32_t i=0; i<uint32_t(count); i++) {
            sequence.push_back(base);
        }
    }

    return sequence;
}



// Return the length of the raw sequence of a read.
// If using the run-length representation of reads, this counts each
// base a number of times equal to its repeat count.
size_t Assembler::getReadRawSequenceLength(ReadId readId)
{

        // We are using the run-length representation.
        // The number of raw bases equals the sum of all
        // the repeat counts.
        // Don't use std::accumulate to compute the sum,
        // otherwise the sum is computed using uint8_t!
        const auto& counts = readRepeatCounts[readId];
        size_t sum = 0;;
        for(uint8_t count: counts) {
            sum += count;
        }
        return sum;

}



// Get a vector of the raw read positions
// corresponding to each position in the run-length
// representation of an oriented read.
vector<uint32_t> Assembler::getRawPositions(OrientedReadId orientedReadId) const
{
    const ReadId readId = orientedReadId.getReadId();
    const ReadId strand = orientedReadId.getStrand();
    const auto repeatCounts = readRepeatCounts[readId];
    const size_t n = repeatCounts.size();

    vector<uint32_t> v;

    uint32_t position = 0;
    for(size_t i=0; i<n; i++) {
        v.push_back(position);
        uint8_t count;
        if(strand == 0) {
            count = repeatCounts[i];
        } else {
            count = repeatCounts[n-1-i];
        }
        position += count;
    }

    return v;
}



void Assembler::initializeReadFlags()
{
    readFlags.createNew(largeDataName("ReadFlags"), largeDataPageSize);
    readFlags.resize(reads.size());
}
void Assembler::accessReadFlags(bool readWriteAccess)
{
    readFlags.accessExisting(largeDataName("ReadFlags"), readWriteAccess);
}



// Write a csv file with summary information for each read.
void Assembler::writeReadsSummary()
{
    SHASTA_ASSERT(reads.isOpen());
    SHASTA_ASSERT(readNames.isOpen());
    SHASTA_ASSERT(markers.isOpen());

    // Count the number of alignment candidates for each read.
    vector<uint64_t> alignmentCandidatesCount(reads.size(), 0);
    for(const OrientedReadPair& p: alignmentCandidates.candidates) {
        ++alignmentCandidatesCount[p.readIds[0]];
        ++alignmentCandidatesCount[p.readIds[1]];
    }


    ofstream csv("ReadSummary.csv");
    csv << "Id,Name,RawLength,RleLength,RawOverRleLengthRatio,"
        "MarkerCount,MarkerDensity,MaximumMarkerOffset,"
        "Palindromic,Chimeric,"
        "AlignmentCandidates,ReadGraphNeighbors,"
        "VertexCount,VertexDensity,\n";
    for(ReadId readId=0; readId!=reads.size(); readId++) {
        const OrientedReadId orientedReadId(readId, 0);

        // Read id.
        csv << readId << ",";

        // Read name.
        const auto readName = readNames[readId];
        copy(readName.begin(), readName.end(), ostream_iterator<char>(csv));
        csv << ",";

        // Number of raw bases.
        const auto repeatCounts = readRepeatCounts[readId];
        uint64_t rawBaseCount = 0;
        for(const auto repeatCount: repeatCounts) {
            rawBaseCount += repeatCount;
        }
        csv << rawBaseCount << ",";

        // Number of RLE bases.
        const uint64_t rleBaseCount = reads[readId].baseCount;
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
        csv << (readFlags[readId].isPalindromic ? "Yes" : "No") << ",";

        // Chimeric flag.
        csv << (readFlags[readId].isChimeric ? "Yes" : "No") << ",";

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

        // End of the line for this read.
        csv << "\n";
     }
}
