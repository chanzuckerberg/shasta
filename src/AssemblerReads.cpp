// shasta.
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
void Assembler::checkReadId(ReadId readId) const
{
    if(readId >= reads.size()) {
        throw runtime_error("Read id " + to_string(readId) +
            " is not valid. Must be between 0 and " + to_string(reads.size()) +
            " inclusive.");
    }
}



// Add reads from a fasta file.
// The reads are added to those already previously present.
void Assembler::addReadsFromFasta(
    const string& fileName,
    size_t minReadLength,
    size_t blockSize,
    const size_t threadCountForReading,
    const size_t threadCountForProcessing)
{
    checkReadsAreOpen();
    checkReadNamesAreOpen();

    ReadLoader readLoader(
        fileName,
        minReadLength,
        blockSize,
        threadCountForReading,
        threadCountForProcessing,
        largeDataFileNamePrefix,
        largeDataPageSize,
        reads,
        readNames,
        readRepeatCounts);

    cout << "Discarded read statistics for file " << fileName << ":" << endl;;
    cout << "    Discarded " << readLoader.discardedShortReadReadCount <<
        " reads shorter than " << minReadLength <<
        " bases for a total " << readLoader.discardedShortReadBaseCount << " bases." << endl;
    cout << "    Discarded " << readLoader.discardedBadRepeatCountReadCount <<
        " reads containing repeat counts 256 or more" <<
        " for a total " << readLoader.discardedBadRepeatCountBaseCount << " bases." << endl;

    // Increment the discarded reads statistics.
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
        cout << "    Read discarded because they were too short or contained repeat counts 256"
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
    // As written, this only works when using raw read representation.
    SHASTA_ASSERT(0);

    ofstream file(fileName);
    for(ReadId readId=0; readId<readCount(); readId++) {
        writeRead(readId, file);
    }

}



void Assembler::writeRead(ReadId readId, const string& fileName)
{
    // As written, this only works when using raw read representation.
    SHASTA_ASSERT(0);

    ofstream file(fileName);
    writeRead(readId, file);
}


void Assembler::writeRead(ReadId readId, ostream& file)
{
    // As written, this only works when using raw read representation.
    SHASTA_ASSERT(0);

    checkReadsAreOpen();
    checkReadNamesAreOpen();
    checkReadId(readId);

    const auto readSequence = reads[readId];
    const auto readName = readNames[readId];

    file << ">" << readId;
    file << " " << readSequence.baseCount << " ";
    copy(readName.begin(), readName.end(), ostream_iterator<char>(file));
    file << "\n";
    file << readSequence << "\n";

}



void Assembler::writeOrientedRead(ReadId readId, Strand strand, const string& fileName)
{
    // As written, this only works when using raw read representation.
    SHASTA_ASSERT(!0);

    writeOrientedRead(OrientedReadId(readId, strand), fileName);
}



void Assembler::writeOrientedRead(OrientedReadId orientedReadId, const string& fileName)
{
    // As written, this only works when using raw read representation.
    SHASTA_ASSERT(0);

    ofstream file(fileName);
    writeOrientedRead(orientedReadId, file);
}



void Assembler::writeOrientedRead(OrientedReadId orientedReadId, ostream& file)
{
    // As written, this only works when using raw read representation.
    SHASTA_ASSERT(0);

    checkReadsAreOpen();
    checkReadNamesAreOpen();

    const ReadId readId = orientedReadId.getReadId();
    checkReadId(readId);
    const Strand strand = orientedReadId.getStrand();
    const auto readSequence = reads[readId];
    const auto readName = readNames[readId];

    file << ">" << readId << "-" << strand;
    file << " " << readSequence.baseCount << " ";
    copy(readName.begin(), readName.end(), ostream_iterator<char>(file));
    file << "\n";
    readSequence.write(file, strand==1);
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
