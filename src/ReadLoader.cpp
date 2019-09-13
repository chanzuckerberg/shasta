// Shasta.
#include "ReadLoader.hpp"
#include "CompressedRunnieReader.hpp"
using namespace shasta;



// Load reads from a fastq or fasta file.
ReadLoader::ReadLoader(
    const string& fileName,
    size_t minReadLength,
    size_t threadCountForReadingArgument,
    size_t threadCountForProcessingArgument,
    const string& dataNamePrefix,
    size_t pageSize,
    LongBaseSequences& reads,
    MemoryMapped::VectorOfVectors<char, uint64_t>& readNames,
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& readRepeatCounts) :

    MultithreadedObject(*this),
    fileName(fileName),
    minReadLength(minReadLength),
    threadCountForReading(threadCountForReadingArgument),
    threadCountForProcessing(threadCountForProcessingArgument),
    dataNamePrefix(dataNamePrefix),
    pageSize(pageSize),
    reads(reads),
    readNames(readNames),
    readRepeatCounts(readRepeatCounts)
{
    adjustThreadCounts();

    // Get the file extension.
    string extension;
    try {
        extension = filesystem::extension(fileName);
    } catch (...) {
        throw runtime_error("Input file " + fileName +
            " must have an extension consistent with its format.");
    }

    // Runnie compressed file.
    if(extension=="rq" || extension=="RQ") {
        processCompressedRunnieFile();
        return;
    }

    // If getting here, the file extension is not supported.
    throw runtime_error("File extension " + extension + " is not supported. "
        "Supported file extensions are .fasta, .fa, .FASTA, .FA.");
}

void ReadLoader::adjustThreadCounts()
{
    if(threadCountForReading == 0) {
        threadCountForReading = 1;
    }
    if(threadCountForProcessing == 0) {
        threadCountForProcessing = std::thread::hardware_concurrency();
    }
}

void ReadLoader::processCompressedRunnieFile()
{
    SHASTA_ASSERT(0);
    /*
    CompressedRunnieReader reader(fileName);
    const uint64_t readCountInFile = reader.countReads();
    cout << "File " << fileName << " contains " << readCountInFile << " reads." << endl;

    CompressedRunnieSequence sequence;
    for(uint64_t i=0; i<readCountInFile; i++) {
        const CompressedRunnieIndex& sequenceInformation = reader.indexes[i];
        reader.readSequence(sequence, i);
    }
    */
}

