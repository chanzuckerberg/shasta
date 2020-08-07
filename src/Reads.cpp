// Shasta
#include "Reads.hpp"

// Standard Library
#include "fstream.hpp"

using namespace shasta;


void Reads::createNew(
    const string& readsDataName,
    const string& readNamesDataName,
    const string& readMetaDataDataName,
    const string& readRepeatCountsDataName,
    const string& readFlagsDataName,
    uint64_t largeDataPageSize)
{
    reads = make_unique<LongBaseSequences>();
    readNames = make_unique<ReadNamesType>();
    readMetaData = make_unique<ReadMetaDataType>();
    readRepeatCounts = make_unique<ReadRepeatCountsType>();
    readFlags = make_unique<ReadFlagsType>();

    reads->createNew(readsDataName, largeDataPageSize);
    readNames->createNew(readNamesDataName, largeDataPageSize);
    readMetaData->createNew(readMetaDataDataName, largeDataPageSize);
    readRepeatCounts->createNew(readRepeatCountsDataName, largeDataPageSize);
    readFlags->createNew(readFlagsDataName, largeDataPageSize);

    this->largeDataPageSize = largeDataPageSize;
}

void Reads::access(
    const string& readsDataName,
    const string& readNamesDataName,
    const string& readMetaDataDataName,
    const string& readRepeatCountsDataName,
    const string& readFlagsDataName)
{
    reads = make_unique<LongBaseSequences>();
    readNames = make_unique<ReadNamesType>();
    readMetaData = make_unique<ReadMetaDataType>();
    readRepeatCounts = make_unique<ReadRepeatCountsType>();
    readFlags = make_unique<ReadFlagsType>();

    reads->accessExistingReadWrite(readsDataName);
    readNames->accessExistingReadWrite(readNamesDataName);
    readMetaData->accessExistingReadWrite(readMetaDataDataName);
    readRepeatCounts->accessExistingReadWrite(readRepeatCountsDataName);
    readFlags->accessExistingReadWrite(readFlagsDataName);
}


void Reads::checkIfAChimericIsAlsoInSmallComponent() const {
    for (const ReadFlags& flags: *readFlags) {
        if (flags.isChimeric) {
            SHASTA_ASSERT(flags.isInSmallComponent);
        }
    }
}


// Return a base of an oriented read.
Base Reads::getOrientedReadBase(
    OrientedReadId orientedReadId,
    uint32_t position) const
{
    const auto& read = (*reads)[orientedReadId.getReadId()];
    if(orientedReadId.getStrand() == 0) {
        return read[position];
    } else {
        return read[read.baseCount-1-position].complement();
    }
}

// Same as above, but also returns the repeat count.
pair<Base, uint8_t> Reads::getOrientedReadBaseAndRepeatCount(
    OrientedReadId orientedReadId,
    uint32_t position) const
{

    // Extract the read id and strand.
    const ReadId readId = orientedReadId.getReadId();
    const Strand strand = orientedReadId.getStrand();

    // Access the bases and repeat counts for this read.
    const auto& read = (*reads)[readId];
    const auto& counts = (*readRepeatCounts)[readId];

    // Compute the position as stored, depending on strand.
    uint32_t orientedPosition = position;
    if(strand == 1) {
        orientedPosition = uint32_t(read.baseCount) - 1 - orientedPosition;
    }

    // Extract the base and repeat count at this position.
    pair<Base, uint8_t> p = make_pair(read[orientedPosition], counts[orientedPosition]);

    // Complement the base, if necessary.
    if(strand == 1) {
        p.first = p.first.complement();
    }

    return p;
}

// Return a vector containing the raw sequence of an oriented read.
vector<Base> Reads::getOrientedReadRawSequence(OrientedReadId orientedReadId) const
{
    // The sequence we will return;
    vector<Base> sequence;

    // The number of bases stored, in run-length representation.
    const uint32_t storedBaseCount = uint32_t((*reads)[orientedReadId.getReadId()].baseCount);


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
uint64_t Reads::getReadRawSequenceLength(ReadId readId) const
{
    // We are using the run-length representation.
    // The number of raw bases equals the sum of all
    // the repeat counts.
    // Don't use std::accumulate to compute the sum,
    // otherwise the sum is computed using uint8_t!
    const auto& counts = (*readRepeatCounts)[readId];
    uint64_t sum = 0;;
    for(uint8_t count: counts) {
        sum += count;
    }
    return sum;
}



// Get a vector of the raw read positions
// corresponding to each position in the run-length
// representation of an oriented read.
vector<uint32_t> Reads::getRawPositions(OrientedReadId orientedReadId) const
{
    const ReadId readId = orientedReadId.getReadId();
    const ReadId strand = orientedReadId.getStrand();
    const auto repeatCounts = (*readRepeatCounts)[readId];
    const uint64_t n = repeatCounts.size();

    vector<uint32_t> v;

    uint32_t position = 0;
    for(uint64_t i=0; i<n; i++) {
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


// Return a meta data field for a read, or an empty string
// if that field is missing. This treats the meta data
// as a space separated sequence of Key=Value,
// without embedded spaces in each Key=Value pair.
span<const char> Reads::getMetaData(ReadId readId, const string& key) const
{
    SHASTA_ASSERT(readId < readMetaData->size());
    const uint64_t keySize = key.size();
    char* keyBegin = const_cast<char*>(&key[0]);
    char* keyEnd = keyBegin + keySize;
    const char* begin = readMetaData->begin(readId);
    const char* end = readMetaData->end(readId);


    const char* p = begin;
    while(p != end) {

        // Look for the next space or line end.
        const char* q = p;
        while(q != end and not isspace(*q)) {
            ++q;
        }

        // When getting here, the interval [p, q) contains
        // a possible (Key,Value) pair.

        // Check if we have the key we are looking for,
        // immediately followed by an equal sign.
        // If so, return what follows.
        if(q > p + keySize + 1) {
            if(std::equal(keyBegin, keyEnd, p)) {
                if(p[keySize] == '=') {
                    const char* valueBegin = p + keySize + 1;
                    const char* valueEnd = q;
                    return span<const char>(valueBegin, valueEnd);
                }
            }
        }

        // If we reached the end of our meta data, stop here.
        if(q == end) {
            break;
        }

        // Look for the next non-space.
        p = q;
        while(p != end and isspace(*p)) {
            ++p;
        }
    }

    // If getting here, we didn't find this keyword.
    // Return an empty string.
    return span<const char>();
}



// Function to write one or all reads in Fasta format.
void Reads::writeReads(const string& fileName)
{
    ofstream file(fileName);
    for(ReadId readId=0; readId<reads->size(); readId++) {
        writeRead(readId, file);
    }
}


void Reads::writeRead(ReadId readId, const string& fileName)
{
    ofstream file(fileName);
    writeRead(readId, file);
}


void Reads::writeRead(ReadId readId, ostream& file)
{
    checkReadsAreOpen();
    checkReadNamesAreOpen();
    checkReadId(readId);

    const vector<Base> rawSequence = getOrientedReadRawSequence(OrientedReadId(readId, 0));
    const auto readName = (*readNames)[readId];
    const auto metaData = (*readMetaData)[readId];

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


void Reads::writeOrientedRead(ReadId readId, Strand strand, const string& fileName)
{
    writeOrientedRead(OrientedReadId(readId, strand), fileName);
}


void Reads::writeOrientedRead(OrientedReadId orientedReadId, const string& fileName)
{
    ofstream file(fileName);
    writeOrientedRead(orientedReadId, file);
}


void Reads::writeOrientedRead(OrientedReadId orientedReadId, ostream& file)
{
    checkReadsAreOpen();
    checkReadNamesAreOpen();

    const vector<Base> rawSequence = getOrientedReadRawSequence(orientedReadId);
    const auto readName = (*readNames)[orientedReadId.getReadId()];

    file << ">" << orientedReadId;
    file << " " << rawSequence.size() << " ";
    copy(readName.begin(), readName.end(), ostream_iterator<char>(file));
    file << "\n";
    copy(rawSequence.begin(), rawSequence.end(), ostream_iterator<Base>(file));
    file << "\n";
}

// Create a histogram of read lengths.
// All lengths here are raw sequence lengths
// (length of the original read), not lengths
// in run-length representation.
void Reads::computeReadLengthHistogram() {
    checkReadsAreOpen();
    histogram.clear();
    // Create the histogram.
    // It contains the number of reads of each length.
    // Indexed by the length.
    const ReadId totalReadCount = readCount();
    totalBaseCount = 0;
    
    for(ReadId readId=0; readId<totalReadCount; readId++) {
        const uint64_t length = getReadRawSequenceLength(readId);
        totalBaseCount += length;
        if(histogram.size() <= length) {
            histogram.resize(length+1, 0);
        }
        ++(histogram[length]);
    }

    // Create binned histogram.
    binnedHistogram.clear();
    const uint64_t binWidth = 1000;
    for(uint64_t length=0; length<histogram.size(); length++) {
        const uint64_t readCount = histogram[length];
        if(readCount) {
            const uint64_t bin = length / binWidth;
            if(binnedHistogram.size() <= bin) {
                binnedHistogram.resize(bin+1, make_pair(0, 0));
            }
            binnedHistogram[bin].first += readCount;
            binnedHistogram[bin].second += readCount * length;
        }
    }

}

void Reads::writeReadLengthHistogram(const string& fileName) {
    checkReadsAreOpen();
    uint64_t totalReadCount = readCount();

    n50 = 0;
    {
        ofstream csv(fileName);
        csv << "Length,Reads,Bases,CumulativeReads,CumulativeBases,"
            "FractionalCumulativeReads,FractionalCumulativeBases,\n";
        uint64_t cumulativeReadCount = totalReadCount;
        uint64_t cumulativeBaseCount = totalBaseCount;
        for(uint64_t length=0; length<histogram.size(); length++) {
            const uint64_t frequency = histogram[length];
            if(frequency) {
                const  uint64_t baseCount = frequency * length;
                const double cumulativeReadFraction =
                    double(cumulativeReadCount)/double(totalReadCount);
                const double comulativeBaseFraction =
                    double(cumulativeBaseCount)/double(baseCount);
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
    }

    // Binned Histogram.
    {
        const uint64_t binWidth = 1000;
        ofstream csv("Binned-" + fileName);
        csv << "LengthBegin,LengthEnd,Reads,Bases,CumulativeReads,CumulativeBases,"
            "FractionalCumulativeReads,FractionalCumulativeBases,\n";
        uint64_t cumulativeReadCount = totalReadCount;
        uint64_t cumulativeBaseCount = totalBaseCount;
        for(uint64_t bin=0; bin<binnedHistogram.size(); bin++) {
            const auto& histogramBin = binnedHistogram[bin];
            const uint64_t readCount = histogramBin.first;
            const uint64_t baseCount = histogramBin.second;
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
