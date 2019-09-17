// Shasta.
#include "ReadLoader.hpp"
#include "CompressedRunnieReader.hpp"
#include "computeRunLengthRepresentation.hpp"
#include "splitRange.hpp"
using namespace shasta;

// Standard library.
#include "chrono.hpp"
#include "iterator.hpp"



// Load reads from a fastq or fasta file.
ReadLoader::ReadLoader(
    const string& fileName,
    size_t minReadLength,
    size_t threadCount,
    const string& dataNamePrefix,
    size_t pageSize,
    LongBaseSequences& reads,
    MemoryMapped::VectorOfVectors<char, uint64_t>& readNames,
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& readRepeatCounts) :

    MultithreadedObject(*this),
    fileName(fileName),
    minReadLength(minReadLength),
    threadCount(threadCount),
    dataNamePrefix(dataNamePrefix),
    pageSize(pageSize),
    reads(reads),
    readNames(readNames),
    readRepeatCounts(readRepeatCounts)
{
    cout << timestamp << "Getting reads from " << fileName << endl;

    adjustThreadCount();

    // Get the file extension.
    string extension;
    try {
        extension = filesystem::extension(fileName);
    } catch (...) {
        throw runtime_error("Input file " + fileName +
            " must have an extension consistent with its format.");
    }

    // Fasta file. ReadLoader is more forgiving than OldFastaReadLoader.
    // This code is not active because in main.cpp we use OldFastaReadLoader
    // (via Assembler::addReadsFromFasta) to read fasta files.
    if(extension=="fasta" || extension=="fa" || extension=="FASTA" || extension=="FA") {
        processFastaFile();
        return;
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

void ReadLoader::adjustThreadCount()
{
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
}



void ReadLoader::processFastaFile()
{
    const auto t0 = std::chrono::steady_clock::now();

    // Read the entire fasta file.
    readFile();

    // Each thread stores reads in its own data structures.
    allocatePerThreadDataStructures();
    runThreads(&ReadLoader::processFastaFileThreadFunction, threadCount);

    // Store the reads computed by each thread and free
    // the per-thread data structures.
    storeReads();
    const auto t1 = std::chrono::steady_clock::now();

    cout << timestamp << "Input file read and processed in " <<
        seconds(t1-t0) << " s." << endl;

}



// Each thread processes the reads whose initial ">" character
// is in the file block assigned to the read.
void ReadLoader::processFastaFileThreadFunction(size_t threadId)
{
    // Allocate and access the data structures where this thread will store the
    // reads it finds.
    allocatePerThreadDataStructures(threadId);
    MemoryMapped::VectorOfVectors<char, uint64_t>& thisThreadReadNames = *threadReadNames[threadId];
    LongBaseSequences& thisThreadReads = *threadReads[threadId];
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& thisThreadReadRepeatCounts =
        *threadReadRepeatCounts[threadId];


    // Compute the file block assigned to this thread.
    uint64_t begin, end;
    tie(begin, end) = splitRange(0, buffer.size(), threadCount, threadId);

    // Locate the first read that begins in this block.
    uint64_t offset = begin;
    if(offset == 0) {
        if(!fastaReadBeginsHere(offset)) {
            throw runtime_error("Fasta file " + fileName + " does not begin with a \">\".");
        }
    } else {
        while(!fastaReadBeginsHere(offset)) {
            ++offset;
            if(offset == end) {
                // We reached the end of the block assigned to this thread
                // without finding any reads.
                return;
            }
        }
    }

    // Ok, we are at the ">" for the first read assigned to this thread.

    // Main loop over the reads in the file block allocated to this thread.
    string readName;
    vector<Base> read;
    vector<Base> runLengthRead;
    vector<uint8_t> readRepeatCount;
    while(offset < end) {
        SHASTA_ASSERT(fastaReadBeginsHere(offset));

        // Consume the ">".
        ++offset;

        // Read the name.
        readName.clear();
        while(true) {
            if(offset == buffer.size()) {
                throw runtime_error("Reached end of file while processing a read name.");
            }
            const char c = buffer[offset++];
            if(isspace(c)) {
                break;
            } else {
                readName.push_back(c);
            }
        }

        // Read the rest of the line.
        while(true) {
            if(offset == buffer.size()) {
                throw runtime_error("Reached end of file while processing header line for read "
                    + readName);
            }
            if(buffer[offset] == '\n') {
                break;
            }
            ++offset;
        }

        // Consume the line end at the end of the header line.
        ++offset;



        // Read the bases.
        // Note that here we can go past the file block assigned to this thread.
        read.clear();
        uint64_t invalidBaseCount = 0;
        while(offset != buffer.size()) {
            const char c = buffer[offset];

            // Skip white space.
            if(c==' ' || c=='\t' || c=='\n' || c=='\r') {
                ++offset;
                continue;
            }

            // If we reached the beginning of another read, stop here.
            if(c=='>' && fastaReadBeginsHere(offset))  {
                break;
            }

            // Get a base from this character.
            const Base base = Base::fromCharacterNoException(c);
            if(base.isValid()) {
                read.push_back(base);
            } else {
                ++invalidBaseCount;
            }
            ++offset;
        }


        // If we found invalid bases, skip this read.
        if(invalidBaseCount) {
            __sync_fetch_and_add(&discardedInvalidBaseReadCount, 1);
            __sync_fetch_and_add(&discardedInvalidBaseReadCount, read.size());
            continue;
        }

        // If the read is too short, skip it.
        if(read.size() < minReadLength) {
            __sync_fetch_and_add(&discardedShortReadReadCount, 1);
            __sync_fetch_and_add(&discardedShortReadBaseCount, read.size());
            continue;
        }

        // Store the read bases.
        if(computeRunLengthRepresentation(read, runLengthRead, readRepeatCount)) {
            thisThreadReadNames.appendVector(readName.begin(), readName.end());
            thisThreadReads.append(runLengthRead);
            thisThreadReadRepeatCounts.appendVector(readRepeatCount);
        } else {
            __sync_fetch_and_add(&discardedBadRepeatCountReadCount, 1);
            __sync_fetch_and_add(&discardedBadRepeatCountBaseCount, read.size());
        }
    }

}



// Function that returns true if a read begins
// at this position in Fasta format.
bool ReadLoader::fastaReadBeginsHere(uint64_t offset) const
{
    if(buffer[offset] == '>') {
        if(offset == 0) {
            return true;
        } else {
            return buffer[offset-1] == '\n';
        }
    } else {
        return false;
    }
}



// Read an entire file into a buffer,
// using threadCountForReading threads.
void ReadLoader::readFile()
{
    // Resize our buffer to contain the entire file.
    int64_t bytesToRead = filesystem::fileSize(fileName);
    buffer.resize(bytesToRead);

    // Open the input file.
    const int fileDescriptor = ::open(fileName.c_str(), O_RDONLY);
    if(fileDescriptor == -1) {
        throw runtime_error("Error opening " + fileName + " for read.");
    }

    // Read it in.
    const auto t0 = std::chrono::steady_clock::now();
    char* bufferPointer = buffer.data();
    while(bytesToRead) {
        const int64_t bytesRead = ::read(fileDescriptor, bufferPointer, bytesToRead);
        if(bytesRead == -1) {
            ::close(fileDescriptor);
            throw runtime_error("Error reading from " + fileName + " near offset " +
                to_string(buffer.size()-bytesToRead));
        }
        bufferPointer += bytesRead;
        bytesToRead -= bytesRead;
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = seconds(t1 - t0);
    cout <<  "File size is " << buffer.size() << " bytes." << endl;
    cout << "Read in " << t01 <<
        " s at " << double(buffer.size()) / t01 << " bytes/s." << endl;


    ::close(fileDescriptor);
}



// Create the name to be used for a MemoryMapped
// object to be used by a thread.
string ReadLoader::threadDataName(
    size_t threadId,
    const string& dataName) const
{
    if(dataNamePrefix.empty()) {
        return "";
    } else {
        return dataNamePrefix + "tmp-OldFastaReadLoader-" + dataName + "-" + to_string(threadId);
    }

}

// Allocate space for the data structures where
// each thread stores the reads it found and their names.
void ReadLoader::allocatePerThreadDataStructures()
{
    threadReadNames.resize(threadCount);
    threadReads.resize(threadCount);
    threadReadRepeatCounts.resize(threadCount);
}
void ReadLoader::allocatePerThreadDataStructures(size_t threadId)
{
    threadReadNames[threadId] = make_shared< MemoryMapped::VectorOfVectors<char, uint64_t> >();
    threadReadNames[threadId]->createNew(
        threadDataName(threadId, "ReadNames"), pageSize);
    threadReads[threadId] = make_shared<LongBaseSequences>();
    threadReads[threadId]->createNew(
        threadDataName(threadId, "Reads"), pageSize);
    threadReadRepeatCounts[threadId] = make_shared< MemoryMapped::VectorOfVectors<uint8_t, uint64_t> >();
    threadReadRepeatCounts[threadId]->createNew(
        threadDataName(threadId, "ReadRepeatCounts"), pageSize);
}



// Compressed runnie file, via class CompressedRunnieReader.
// For now this is sequential but it could be made faster
// using CompressedRunnieReader functionality.
void ReadLoader::processCompressedRunnieFile()
{
    // Create the CompressedRunnieReader.
    CompressedRunnieReader reader(fileName);
    const uint64_t readCountInFile = reader.getReadCount();
    cout << "File " << fileName << " contains " << readCountInFile << " reads." << endl;

    // Store the reads one by one.
    NamedCompressedRunnieSequence read;
    for(uint64_t i=0; i!=readCountInFile; i++) {
        reader.getSequenceData(read, i);

        readNames.appendVector(read.name.begin(), read.name.begin());
        reads.append(read.sequence.size());
        LongBaseSequenceView storedSequence = reads[reads.size()-1];
        for(uint64_t i=0; i<read.sequence.size(); i++) {
            storedSequence.set(i, Base::fromCharacter(read.sequence[i]));
        }
        readRepeatCounts.appendVector(read.encoding);
    }
}



// Store the reads computed by each thread and free
// the per-thread data structures.
void ReadLoader::storeReads()
{
    // Loop over all threads.
    for(size_t threadId=0; threadId<threadCount; threadId++) {

        // Access the names.
        MemoryMapped::VectorOfVectors<char, uint64_t>& thisThreadReadNames =
            *(threadReadNames[threadId]);

        // Access the reads.
        LongBaseSequences& thisThreadReads = *(threadReads[threadId]);
        const size_t n = thisThreadReadNames.size();
        SHASTA_ASSERT(thisThreadReads.size() == n);

        // Access the repeat counts.
        MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& thisThreadReadRepeatCounts =
            *threadReadRepeatCounts[threadId];
        SHASTA_ASSERT(thisThreadReadRepeatCounts.size() == n);

        // Store the reads.
        for(size_t i=0; i<n; i++) {
            readNames.appendVector(thisThreadReadNames.begin(i), thisThreadReadNames.end(i));
            reads.append(thisThreadReads[i]);
            const size_t j = readRepeatCounts.size();
            readRepeatCounts.appendVector(thisThreadReadRepeatCounts.size(i));
            copy(
                thisThreadReadRepeatCounts.begin(i),
                thisThreadReadRepeatCounts.end(i),
                readRepeatCounts.begin(j));
        }

        // Remove the data structures used by this thread.
        thisThreadReadNames.remove();
        thisThreadReads.remove();
        thisThreadReadRepeatCounts.remove();
    }

    // Clear the per-thread data structures.
    threadReadNames.clear();
    threadReads.clear();
    threadReadRepeatCounts.clear();

}

