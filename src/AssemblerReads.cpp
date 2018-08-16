// shasta.
#include "Assembler.hpp"
#include "ReadLoader.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Standard libraries.
#include "iterator.hpp"



void Assembler::accessReadsReadOnly()
{
    reads.accessExistingReadOnly(largeDataName("Reads"));
    if(assemblerInfo->useRunLengthReads) {
        readRepeatCounts.accessExistingReadOnly(largeDataName("ReadRepeatCounts"));
    }
}
void Assembler::accessReadsReadWrite()
{
    reads.accessExistingReadWrite(largeDataName("Reads"));
    if(assemblerInfo->useRunLengthReads) {
        readRepeatCounts.accessExistingReadWrite(largeDataName("ReadRepeatCounts"));
    }
}
void Assembler::accessReadNamesReadOnly()
{
    readNames.accessExistingReadOnly(largeDataName("ReadNames"));
}
void Assembler::accessReadNamesReadWrite()
{
    readNames.accessExistingReadWrite(largeDataName("ReadNames"));
}
void Assembler::checkReadsAreOpen() const
{
    if(!reads.isOpen()) {
        throw runtime_error("Reads are not accessible.");
    }
    if(assemblerInfo->useRunLengthReads && !readRepeatCounts.isOpen()) {
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

    ReadLoader(
        fileName,
        minReadLength,
        assemblerInfo->useRunLengthReads,
        blockSize,
        threadCountForReading,
        threadCountForProcessing,
        largeDataFileNamePrefix,
        largeDataPageSize,
        reads,
        readNames,
        readRepeatCounts);

}


// Create a histogram of read lengths.
void Assembler::histogramReadLength(const string& fileName)
{
    // As written, this only works when using raw read representation.
    CZI_ASSERT(!assemblerInfo->useRunLengthReads);

    checkReadsAreOpen();

    // Create the histogram.
    vector<size_t> histogram;
    for(ReadId readId=0; readId<readCount(); readId++) {
        const size_t length = reads[readId].baseCount;
        if(histogram.size() <= length) {
            histogram.resize(length+1, 0);
        }
        ++(histogram[length]);
    }

    // Write it out.
    ofstream csv(fileName);
    csv << "Length,Frequency,Bases,CumulativeFrequency,CumulativeBases\n";
    size_t cumulativeFrequency = 0;
    size_t cumulativeBaseCount = 0;
    for(size_t length=0; length<histogram.size(); length++) {
        const size_t frequency = histogram[length];
        if(frequency) {
            const  size_t baseCount = frequency * length;
            cumulativeFrequency += frequency;
            cumulativeBaseCount += baseCount;
            csv << length << "," << frequency << "," << baseCount << ",";
            csv<< cumulativeFrequency << "," << cumulativeBaseCount << "\n";
        }
    }

}



// Function to write one or all reads in Fasta format.
void Assembler::writeReads(const string& fileName)
{
    // As written, this only works when using raw read representation.
    CZI_ASSERT(!assemblerInfo->useRunLengthReads);

    ofstream file(fileName);
    for(ReadId readId=0; readId<readCount(); readId++) {
        writeRead(readId, file);
    }

}



void Assembler::writeRead(ReadId readId, const string& fileName)
{
    // As written, this only works when using raw read representation.
    CZI_ASSERT(!assemblerInfo->useRunLengthReads);

    ofstream file(fileName);
    writeRead(readId, file);
}


void Assembler::writeRead(ReadId readId, ostream& file)
{
    // As written, this only works when using raw read representation.
    CZI_ASSERT(!assemblerInfo->useRunLengthReads);

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
    CZI_ASSERT(!assemblerInfo->useRunLengthReads);

    writeOrientedRead(OrientedReadId(readId, strand), fileName);
}



void Assembler::writeOrientedRead(OrientedReadId orientedReadId, const string& fileName)
{
    // As written, this only works when using raw read representation.
    CZI_ASSERT(!assemblerInfo->useRunLengthReads);

    ofstream file(fileName);
    writeOrientedRead(orientedReadId, file);
}



void Assembler::writeOrientedRead(OrientedReadId orientedReadId, ostream& file)
{
    // As written, this only works when using raw read representation.
    CZI_ASSERT(!assemblerInfo->useRunLengthReads);

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

    /// The number of bases stored.
    const uint32_t storedBaseCount = uint32_t(reads[orientedReadId.getReadId()].baseCount);

    if(assemblerInfo->useRunLengthReads) {

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

    } else {

        // We are storing raw read sequence, so we can just copy the stored sequence.
        sequence.resize(storedBaseCount);
        for(size_t position=0; position<storedBaseCount; position++) {
            sequence[position] = getOrientedReadBase(orientedReadId, uint32_t(position));
        }

    }

    return sequence;
}



// Get a vector of the raw read positions
// corresponding to each position in the run-length
// representation of an oriented read.
vector<uint32_t> Assembler::getRawPositions(OrientedReadId orientedReadId) const
{
    CZI_ASSERT(assemblerInfo->useRunLengthReads);
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
