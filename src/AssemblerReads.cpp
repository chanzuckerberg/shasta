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
}
void Assembler::accessReadsReadWrite()
{
    reads.accessExistingReadWrite(largeDataName("Reads"));
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
        blockSize,
        threadCountForReading,
        threadCountForProcessing,
        largeDataFileNamePrefix,
        largeDataPageSize,
        reads,
        readNames);

}


// Create a histogram of read lengths.
void Assembler::histogramReadLength(const string& fileName)
{
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


