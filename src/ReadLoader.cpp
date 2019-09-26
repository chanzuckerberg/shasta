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
    cout << timestamp << "Loading reads from " << fileName << endl;

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
    if(extension=="fasta" || extension=="fa" || extension=="FASTA" || extension=="FA") {
        processFastaFile();
        return;
    }

#if 0
    // Fastq file.
    if(extension=="fastq" || extension=="fq" || extension=="FASTQ" || extension=="FQ") {
        processFastqFile();
        return;
    }
#endif

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

    // Read the entire fasta file.
    const auto t0 = std::chrono::steady_clock::now();
    readFile();

    // Each thread stores reads in its own data structures.
    const auto t1 = std::chrono::steady_clock::now();
    allocatePerThreadDataStructures();
    runThreads(&ReadLoader::processFastaFileThreadFunction, threadCount);
    const auto t2 = std::chrono::steady_clock::now();

    // Store the reads computed by each thread and free
    // the per-thread data structures.
    storeReads();
    buffer.remove();
    const auto t3 = std::chrono::steady_clock::now();


    cout << "Time to process this file:\n" <<
        "Allocate buffer + read: " << seconds(t1-t0) << " s.\n" <<
        "Parse: " << seconds(t2-t1) << " s.\n"
        "Store: " << seconds(t3-t2) << " s.\n"
        "Total: " << seconds(t3-t0) << " s." << endl;


}



// Each thread processes the reads whose initial ">" character
// is in the file block assigned to the read.
void ReadLoader::processFastaFileThreadFunction(size_t threadId)
{
    const char* bufferPointer = &buffer[0];
    const uint64_t bufferSize = buffer.size();

    // Allocate and access the data structures where this thread will store the
    // reads it finds.
    allocatePerThreadDataStructures(threadId);
    MemoryMapped::VectorOfVectors<char, uint64_t>& thisThreadReadNames = *threadReadNames[threadId];
    LongBaseSequences& thisThreadReads = *threadReads[threadId];
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& thisThreadReadRepeatCounts =
        *threadReadRepeatCounts[threadId];


    // Compute the file block assigned to this thread.
    uint64_t begin, end;
    tie(begin, end) = splitRange(0, bufferSize, threadCount, threadId);
    if(begin == end) {
        return;
    }

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
            if(offset == bufferSize) {
                throw runtime_error("Reached end of file while processing a read name.");
            }
            const char c = bufferPointer[offset];
            if(isspace(c)) {
                break;
            } else {
                readName.push_back(c);
            }
            offset++;
        }

        // Read the rest of the line.
        while(true) {
            if(offset == bufferSize) {
                throw runtime_error("Reached end of file while processing header line for read "
                    + readName);
            }
            if(bufferPointer[offset] == '\n') {
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
        while(offset != bufferSize) {
            const char c = bufferPointer[offset];

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



// Process a fastq file under the following assumptions:
// - 4 lines per read.
// - For each read:
//   * The first line must begin with '@'.
//   * The third line must consist of just a '+' with nothing else.
//   * The second and fourth line must have the same number of bytes,
//     also equal to the number of bases in the read. The fourth
//     line is otherwise ignored.
//   * The second line can only contain ACGT. Invalid bases are not allowed.
// - No Windows line ends.
void ReadLoader::processFastqFile()
{

    // Read the entire fastq file.
    const auto t0 = std::chrono::steady_clock::now();
    readFile();

    // Find all line ends in the file.
    const auto t1 = std::chrono::steady_clock::now();
    findLineEnds();
    cout << "Found " << lineEnds.size() << " lines in this file." << endl;

    // Check that the number of lines is a multiple of 4
    // (there must be exactly 4 mlines per read per the above assumptions).
    if((lineEnds.size() %4) != 0) {
        throw runtime_error("File has " + to_string(lineEnds.size()) +
            " lines. Expected a multiple of 4. "
            "Only fastq files with each read on exactly 4 lines are supported.");
    }

    // Each thread stores reads in its own data structures.
    const auto t2 = std::chrono::steady_clock::now();
    allocatePerThreadDataStructures();
    runThreads(&ReadLoader::processFastqFileThreadFunction, threadCount);


    // Store the reads computed by each thread and free
    // the per-thread data structures.
    const auto t3 = std::chrono::steady_clock::now();
    storeReads();
    buffer.remove();
    const auto t4 = std::chrono::steady_clock::now();


    cout << "Time to process this file:\n" <<
        "Allocate buffer + read: " << seconds(t1-t0) << " s.\n" <<
        "Locate: " << seconds(t2-t1) << " s.\n"
        "Parse: " << seconds(t3-t2) << " s.\n"
        "Store: " << seconds(t4-t3) << " s.\n"
        "Total: " << seconds(t4-t0) << " s." << endl;
}



void ReadLoader::processFastqFileThreadFunction(size_t threadId)
{
    SHASTA_ASSERT(0);
#if 0
    const char* bufferPointer = &buffer[0];
    const uint64_t bufferSize = buffer.size();

    // Allocate and access the data structures where this thread will store the
    // reads it finds.
    allocatePerThreadDataStructures(threadId);
    MemoryMapped::VectorOfVectors<char, uint64_t>& thisThreadReadNames = *threadReadNames[threadId];
    LongBaseSequences& thisThreadReads = *threadReads[threadId];
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& thisThreadReadRepeatCounts =
        *threadReadRepeatCounts[threadId];

    // Find the total number of reads in the file.
    SHASTA_ASSERT((lineEnds.size() % 4) == 0); // We already checked for that.
    const ReadId readCountInFile = ReadId(lineEnds.size() / 4);

    // Compute the range of reads in the file assigned to this thread.
    uint64_t begin, end;
    tie(begin, end) = splitRange(0, readCountInFile, threadCount, threadId);
    if(begin == end) {
        return;
    }



    // Loop over this range of reads.
    for(uint64_t i=begin; i!=end; i++) {

        // Locate the line ends corresponding to this read.
        /*
        const auto thisReadLineEnds = lineEnds.begin() + i * 4;
        const auto headerOffset = thisReadLineEnds[0];
        const auto basesOffset = thisReadLineEnds[1];
        const auto plusOffset = thisReadLineEnds[2];
        const auto scoresOffset = thisReadLineEnds[3];
        */

    }
#if 0
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
            if(offset == bufferSize) {
                throw runtime_error("Reached end of file while processing a read name.");
            }
            const char c = bufferPointer[offset];
            if(isspace(c)) {
                break;
            } else {
                readName.push_back(c);
            }
            offset++;
        }

        // Read the rest of the line.
        while(true) {
            if(offset == bufferSize) {
                throw runtime_error("Reached end of file while processing header line for read "
                    + readName);
            }
            if(bufferPointer[offset] == '\n') {
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
        while(offset != bufferSize) {
            const char c = bufferPointer[offset];

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
#endif
#endif
}



// Read an entire file into a buffer,
// using threadCountForReading threads.
void ReadLoader::readFile()
{
    // Resize our buffer to contain the entire file.
    const auto t0 = std::chrono::steady_clock::now();
    int64_t bytesToRead = filesystem::fileSize(fileName);
    buffer.createNew(dataName("tmp-FastaBuffer"), pageSize);
    buffer.resize(bytesToRead);

    // Open the input file.
    const int fileDescriptor = ::open(fileName.c_str(), O_RDONLY);
    if(fileDescriptor == -1) {
        throw runtime_error("Error opening " + fileName + " for read.");
    }

    // Read it in.
    const auto t1 = std::chrono::steady_clock::now();
    char* bufferPointer = &buffer[0];
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
    ::close(fileDescriptor);
    const auto t2 = std::chrono::steady_clock::now();
    const double t01 = seconds(t1 - t0);
    const double t12 = seconds(t2 - t1);

    cout <<  "File size: " << buffer.size() << " bytes." << endl;
    cout << "Allocate buffer time: " << t01 << " s." << endl;
    cout << "Read time: " << t12 << " s." << endl;
    cout << "Read rate: " << double(buffer.size()) / t12 << " bytes/s." << endl;


}



// Find all of the line ends ('\n') in the buffer.
void ReadLoader::findLineEnds()
{
    // Each thread finds line ends in a block of the file
    // and stores them in its own vector.
    threadLineEnds.resize(threadCount);
    runThreads(&ReadLoader::findLineEndsThreadFunction, threadCount);

    // Combine the line ends founds by all threads.
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        const vector<uint64_t>& thisThreadLineEnds = threadLineEnds[threadId];
        lineEnds.insert(lineEnds.end(),
            thisThreadLineEnds.begin(), thisThreadLineEnds.end());
    }
    threadLineEnds.clear();
}
void ReadLoader::findLineEndsThreadFunction(size_t threadId)
{
    // Access the vector where this thread will store the line ends it finds.
    vector<uint64_t>& thisThreadLineEnds = threadLineEnds[threadId];

    // Compute the file block assigned to this thread.
    uint64_t begin, end;
    tie(begin, end) = splitRange(0, buffer.size(), threadCount, threadId);
    if(begin == end) {
        return;
    }

    // Look for line ends in this block.
    for(uint64_t offset=begin; offset!=end; offset++) {
        if(buffer[offset] == '\n') {
            thisThreadLineEnds.push_back(offset);
        }
    }
}



// Create the name to be used for a MemoryMapped object.
string ReadLoader::dataName(
    const string& dataName) const
{
    if(dataNamePrefix.empty()) {
        return "";
    } else {
        return dataNamePrefix + dataName;
    }

}
string ReadLoader::threadDataName(
    size_t threadId,
    const string& dataName) const
{
    if(dataNamePrefix.empty()) {
        return "";
    } else {
        return dataNamePrefix + "tmp-ReadLoader-" + dataName + "-" + to_string(threadId);
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
    const auto t0 = std::chrono::steady_clock::now();

    // Create the CompressedRunnieReader.
    compressedRunnieReader = make_shared<CompressedRunnieReader>(fileName);
    CompressedRunnieReader& reader = *compressedRunnieReader;
    const uint64_t readCountInFile = reader.getReadCount();
    cout << "Input file contains " << readCountInFile << " reads." << endl;

    // Use single-threaded code to create the space.
    readIdTable.resize(readCountInFile);
    ReadId readId = ReadId(reads.size());
    for(uint64_t i=0; i!=readCountInFile; i++) {
        const uint64_t baseCount = reader.getLength(i);
        if(baseCount >= minReadLength) {
            readNames.appendVector(reader.getReadName(i).size());
            reads.append(baseCount);
            readRepeatCounts.appendVector(baseCount);
            readIdTable[i] = readId++;
        } else {
            discardedShortReadReadCount++;
            discardedShortReadBaseCount += baseCount;
            readIdTable[i] = invalidReadId;
        }
    }

    // Use multithreaded code to store the reads.
    setupLoadBalancing(readIdTable.size(), 1000);
    runThreads(&ReadLoader::processCompressedRunnieFileThreadFunction, threadCount);
    compressedRunnieReader = 0;

    const auto t1 = std::chrono::steady_clock::now();
    cout << "Input file read and processed in " <<
        seconds(t1-t0) << " s." << endl;
}



void ReadLoader::processCompressedRunnieFileThreadFunction(size_t threadId)
{
    CompressedRunnieReader& reader = *compressedRunnieReader;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    NamedCompressedRunnieSequence read;
    while(getNextBatch(begin, end)) {

        // Loop over all reads in this batch.
        for(uint64_t i=begin; i!=end; i++) {
            const ReadId readId = readIdTable[i];
            if(readId == invalidReadId) {
                continue;
            }
            reader.getSequenceData(read, i);
            copy(read.name.begin(), read.name.end(), readNames.begin(readId));
            LongBaseSequenceView storedSequence = reads[readId];
            for(uint64_t j=0; j<read.sequence.size(); j++) {
                storedSequence.set(j, Base::fromCharacter(read.sequence[j]));
            }
            copy(read.encoding.begin(), read.encoding.end(), readRepeatCounts.begin(readId));
        }
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

