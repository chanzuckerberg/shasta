// Shasta.
#include "ReadLoader.hpp"
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
    minReadLength(minReadLength),
    threadCountForReading(threadCountForReadingArgument),
    threadCountForProcessing(threadCountForProcessingArgument)
{
    // Get the file extension.
    string extension;
    try {
        extension = filesystem::extension(fileName);
    } catch (...) {
        throw runtime_error("Input file " + fileName +
            " must have an extension consistent with its format.");
    }


    // If getting here, the file extension is not supported.
    throw runtime_error("File extension " + extension + " is not supported. "
        "Supported file extensions are .fasta, .fa, .FASTA, .FA.");
}
