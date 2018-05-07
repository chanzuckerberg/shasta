// Nanopore2.
#include "ReadLoader.hpp"
#include "splitRange.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Standard library.
#include "tuple.hpp"



// Load reads from a fastq or fasta file.
ReadLoader::ReadLoader(
    const string& fileName,
    size_t blockSize,
    size_t threadCountForReading,
    size_t threadCountForProcessing,
    LongBaseSequences& reads,
    MemoryMapped::VectorOfVectors<char, uint64_t>& readNames) :
    MultithreadedObject(*this),
    blockSize(blockSize),
    threadCountForProcessing(threadCountForProcessing)
{
    cout << timestamp << "Loading reads from " << fileName << "." << endl;
    cout << "Input file block size: " << blockSize << " bytes." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCountForReading == 0) {
        threadCountForReading = 1;
    }
    if(threadCountForProcessing == 0) {
        threadCountForProcessing = std::thread::hardware_concurrency();
    }
    cout << "Using " << threadCountForReading << " threads for reading and ";
    cout << threadCountForProcessing << " threads for processing." << endl;

    // Allocate space to keep a block of the file.
    // Reserve extra space to allow for reads spanning blocks.
    const size_t extraReservedSpace = 6*1024*1024;
    buffer.reserve(blockSize + extraReservedSpace);

    // Open the input file.
    fileDescriptor = ::open(fileName.c_str(), O_RDONLY);
    if(fileDescriptor == -1) {
        throw runtime_error("Error opening " + fileName + " for read.");
    }

    // Find the size of the input file.
    getFileSize();
    cout << "Input file size is " << fileSize << " bytes." << endl;



    // Main loop over blocks in the input file.
    for(blockBegin=0; blockBegin<fileSize; ) {

        // Read this block.
        blockEnd = min(blockBegin+blockSize, fileSize);
        cout << "Working on block " << blockBegin << " " << blockEnd << endl;
        readBlock(threadCountForReading);
        cout << "Block updated to " << blockBegin << " " << blockEnd << endl;

        // Process this block in parallel.
        runThreads(&ReadLoader::threadFunction, threadCountForProcessing);

        // Prepare to process the next block.
        blockBegin = blockEnd;
    }


    // Close the input file.
    ::close(fileDescriptor);

    cout << timestamp << "Done loading reads." << endl;
}



void ReadLoader::getFileSize()
{
    CZI_ASSERT(fileDescriptor != -1);

    struct stat buffer;
    if(::fstat(fileDescriptor, &buffer)) {
        throw runtime_error("Error from fstat.");
    }
    fileSize = buffer.st_size;
}



// Read one block into the above buffer.
// This reads into the buffer the portion of the input file
// at offset in [blockBegin, blockEnd).
// If this is not the last block in the file,
// it then keeps only the portion until the last '>'
// preceded by a '\n'. The discarded characters are stored
// as leftover and will be used when processing the next block.
void ReadLoader::readBlock(size_t threadCount)
{
    if(threadCount <= 1) {
        readBlockSequential();
    } else {
        readBlockParallel(threadCount);
    }

    if(buffer.front() != '>') {
        throw runtime_error("Expected '>' at beginning of a block.");
    }

    if(blockEnd != fileSize) {
        // Go back to the last '>' preceded by '\n'.
        size_t bufferIndex=buffer.size()-1;
        for(; bufferIndex>0; bufferIndex--) {
            if(buffer[bufferIndex] == '>' && buffer[bufferIndex-1]=='\n') {
                break;
            }
        }
        CZI_ASSERT(buffer[bufferIndex]=='>' && (bufferIndex==0 || buffer[bufferIndex-1]=='\n'));

        // Throw away what follows, store it as leftover.
        leftOver.clear();
        leftOver.resize(buffer.size() - bufferIndex);
        copy(buffer.begin()+bufferIndex, buffer.end(), leftOver.begin());
        buffer.resize(bufferIndex);
    }
}



void ReadLoader::readBlockSequential()
{

    size_t bytesToRead = blockEnd - blockBegin;
    buffer.resize(leftOver.size() + bytesToRead);
    copy(leftOver.begin(), leftOver.end(), buffer.begin());
    char* bufferPointer = buffer.data() + leftOver.size();
    while(bytesToRead) {
        const ssize_t byteCount = ::read(fileDescriptor, bufferPointer, bytesToRead);
        if(byteCount <= 0) {
            throw runtime_error("Error during read.");
        }
        bytesToRead -= byteCount;
        bufferPointer += byteCount;
    }

}



void ReadLoader::readBlockParallel(size_t threadCount)
{
    CZI_ASSERT(0);
}



void ReadLoader::threadFunction(size_t threadId)
{

    // Get the slice of the buffer assigned to this thread.
    size_t sliceBegin, sliceEnd;
    tie(sliceBegin, sliceEnd) = splitRange(0, buffer.size(), threadCountForProcessing, threadId);

    cout << "Thread function called " << threadId << " " << sliceBegin << " " << sliceEnd << endl;


    // Locate the first read that begins in the slice assigned to this thread.
    // We only process reads that begin in this slice.
    size_t bufferIndex = sliceBegin;
    while(bufferIndex<buffer.size() && !readBeginsHere(bufferIndex)) {
        cout << bufferIndex << endl;
        ++bufferIndex;
    }
    if(threadId == 0) {
        CZI_ASSERT(bufferIndex == 0);
    }


    // Main loop over the buffer slice assigned to this thread.
    cout << "Thread main loop " << threadId << " " << bufferIndex << endl;
    string readName;
    while(bufferIndex < buffer.size()) {

        // Skip the '>' that introduces the new read.
        CZI_ASSERT(buffer[bufferIndex++] == '>');

        // Extract the read name and discard the rest of the line.
        bool blankFound = false;
        readName.clear();
        while(bufferIndex < buffer.size()) {
            const char c = buffer[bufferIndex++];
            if(c == '\n') {
                break;
            }
            if(c==' ') {
                blankFound = true;
            }
            if(!blankFound) {
                readName.push_back(c);
            }
        }
        cout << readName << endl;


        // Read the base characters.
        while(bufferIndex < buffer.size()) {
            const char c = buffer[bufferIndex++];
            if(c == '\n') {
                break;
            }
            // cout << c;
        }
        // cout << endl;

    }
}



// Return true if a read begins at this position in the buffer.
bool ReadLoader::readBeginsHere(size_t bufferIndex) const
{
    const char c = buffer[bufferIndex];
    if(bufferIndex == 0) {
        CZI_ASSERT(c == '>');
        return true;
    } else {
        return c=='>' && buffer[bufferIndex-1]=='\n';
    }
}


