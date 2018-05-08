// Nanopore2.
#include "Assembler.hpp"
#include "ReadLoader.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Standard libraries.
#include "stdexcept.hpp"



Assembler::Assembler(
    const string& smallDataFileNamePrefix,
    const string& largeDataFileNamePrefix,
    size_t smallDataPageSize,
    size_t largeDataPageSize) :
    MultithreadedObject(*this),
    smallDataFileNamePrefix(smallDataFileNamePrefix),
    largeDataFileNamePrefix(largeDataFileNamePrefix),
    smallDataPageSize(smallDataPageSize),
    largeDataPageSize(largeDataPageSize)
{
    try {
        assemblerInfo.accessExistingReadWrite(smallDataName("Info"));
    } catch(runtime_error) {
        assemblerInfo.createNew(smallDataName("Info"));
    }
}



// Add reads from a fasta file.
// The reads are added to those already previously present.
void Assembler::addReadsFromFasta(
    const string& fileName,
    size_t blockSize,
    const size_t threadCountForReading,
    const size_t threadCountForProcessing)
{
    reads.accessExistingReadWriteOrCreateNew(largeDataName("Reads"), largeDataPageSize);
    readNames.accessExistingReadWriteOrCreateNew(largeDataName("ReadNames"), largeDataPageSize);

    ReadLoader(
        fileName,
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
    // Access the reads.
    reads.accessExistingReadWriteOrCreateNew(largeDataName("Reads"), largeDataPageSize);

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

}



void Assembler::writeRead(ReadId, const string& fileName)
{

}


