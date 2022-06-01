#include "Assembler.hpp"
#include "deduplicate.hpp"
#include "MurmurHash2.hpp"
#include "Reads.hpp"
using namespace shasta;

#include <random>



void Assembler::accessKmers()
{
    kmerTable.accessExistingReadOnly(largeDataName("Kmers"));
    if(kmerTable.size() != (1ULL<< (2*assemblerInfo->k))) {
        throw runtime_error("Size of k-mer vector is inconsistent with stored value of k.");
    }
}

void Assembler::checkKmersAreOpen()const
{
    if(!kmerTable.isOpen) {
        throw runtime_error("Kmers are not accessible.");
    }
}



// Randomly select the k-mers to be used as markers.
void Assembler::randomlySelectKmers(
    size_t k,           // k-mer length.
    double probability, // The probability that a k-mer is selected as a marker.
    int seed            // For random number generator.
)
{
    // Sanity check on the value of k, then store it.
    if(k > Kmer::capacity) {
        throw runtime_error("K-mer capacity exceeded.");
    }
    assemblerInfo->k = k;

    // The total number of k-mers of this length.
    // This includes both RLE and non-RLE k-mers.
    const size_t kmerCount = 1ULL << (2ULL*k);

    // Sanity check on the requested fraction.
    // It can be 1 at most. If it is 1, all k-mers
    // are guaranteed to be selected (because we use <=
    // comparison in the main loop.
    if(probability<0. || probability>1.) {
        throw runtime_error("Invalid k-mer probability " +
            to_string(probability) + " requested.");
    }



    // Fill in the fields of the k-mer table
    // that depends only on k.
    initializeKmerTable();



    // Compute the probability p with which we select
    // each k-mer and its reverse complement
    // in order to achieve the required k-mer fraction.
    // For k-mers that are not self-complementary,
    // the probability of not being selected
    // is (1-p)^2 (there are two chances,
    // with and without reverse complement).
    // So, probability = 1 - (1-p)^2, and therefore p=1-sqrt(1-probability).
    // For simplicity, we use the same p for k-mers that are
    // self-complementary. They are a small minority, and because
    // of this they are chose with lower probability.
    // Probably a good thing anyway.
    const double p = 1. - sqrt(1. - probability);
    if(probability == 1.) {
        SHASTA_ASSERT(p == 1.);
    }



    // Prepare to generate uniformly distributed numbers between 0 and 1.
    std::mt19937 randomSource(seed);
    std::uniform_real_distribution<> uniformDistribution;

    // Pick each k-mer and its reverse complement with probability p.
    // Use <= comparison, so if probability=1, p=1, all k-mers are kept.
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        const double x = uniformDistribution(randomSource);
        if(x <= p) {
            kmerTable[kmerId].isMarker = true;
            const uint64_t reverseComplementedKmerId = kmerTable[kmerId].reverseComplementedKmerId;
            kmerTable[reverseComplementedKmerId].isMarker = true;
        }
    }


    // Do some counting.
    uint64_t rleKmerCount = 0;
    uint64_t markerKmerCount = 0;
    uint64_t markerRleKmerCount = 0;
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        const KmerInfo& kmerInfo = kmerTable[kmerId];
        if(kmerInfo.isRleKmer) {
            ++rleKmerCount;
        }
        if(kmerInfo.isMarker) {
            ++markerKmerCount;
        }
        if(kmerInfo.isRleKmer and kmerInfo.isMarker) {
            ++markerRleKmerCount;
        }
    }



    // Summary messages.
    if(assemblerInfo->readRepresentation == 0) {

        // We are using the raw representation of the reads.
        cout << "Total number of k-mers of length " << k << " is " << kmerCount << endl;
        cout << "Of those, " << markerKmerCount << " will be used as markers." << endl;
        cout << "Fraction of k-mers used as markers: requested " <<
            probability << ", actual " <<
            double(markerKmerCount)/double(kmerCount) << "." << endl;


    } else {

        // We are using the RLE representation of the reads.
        cout << "Total number of k-mers of length " << k << " is " << kmerCount << endl;
        cout << "Number of RLE k-mers of length " << k << " is " << rleKmerCount << endl;
        cout << "Of those, " << markerRleKmerCount << " will be used as markers." << endl;
        cout << "Fraction of k-mers used as markers: requested " <<
            probability << ", actual " <<
            double(markerRleKmerCount)/double(rleKmerCount) << "." << endl;

    }

}



void Assembler::initializeKmerTable()
{
    // Create the kmer table with the necessary size.
    kmerTable.createNew(largeDataName("Kmers"), largeDataPageSize);
    const size_t k = assemblerInfo->k;
    const size_t kmerCount = 1ULL << (2ULL*k);
    kmerTable.resize(kmerCount);

    // Store the reverse complement of each k-mer.
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        const Kmer kmer(kmerId, k);
        const Kmer reverseComplementedKmer = kmer.reverseComplement(k);
        kmerTable[kmerId].reverseComplementedKmerId = KmerId(reverseComplementedKmer.id(k));
    }
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        const uint64_t reverseComplementedKmerId = kmerTable[kmerId].reverseComplementedKmerId;
        SHASTA_ASSERT(kmerTable[reverseComplementedKmerId].reverseComplementedKmerId == kmerId);
    }



    // Figure out which k-mers are run-length-encoded sequence.
    // They are the ones without repeated bases.
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        kmerTable[kmerId].isRleKmer = true;
        const Kmer kmer(kmerId, k);
        for(size_t i=1; i<k; i++) {
            if(kmer[i-1] == kmer[i]) {
                kmerTable[kmerId].isRleKmer = false;
                break;
            }
        }
    }


    // Fill in hash values used for downsampling.
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        const uint64_t n = kmerId + kmerTable[kmerId].reverseComplementedKmerId;
        kmerTable[kmerId].hash = MurmurHash2(&n, sizeof(n), 13477);
    }

}



void Assembler::writeKmers(const string& fileName) const
{
    checkKmersAreOpen();

    // Get the k-mer length.
    const size_t k = assemblerInfo->k;
    const size_t kmerCount = 1ULL << (2ULL*k);
    SHASTA_ASSERT(kmerTable.size() == kmerCount);

    // Open the output file and write the header line.
    ofstream file(fileName);
    file << "KmerId,Kmer,IsMarker,ReverseComplementedKmerId,ReverseComplementedKmer\n";

    // Write a line for each k-mer.
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        file << kmerId << ",";
        file << Kmer(kmerId, k) << ",";
        file << int(kmerTable[kmerId].isMarker) << ",";
        file << kmerTable[kmerId].reverseComplementedKmerId << ",";
        file << Kmer(kmerTable[kmerId].reverseComplementedKmerId, k) << "\n";
    }
}



// Select marker k-mers randomly, but excluding
// the ones that have high frequency in the reads.
void Assembler::selectKmersBasedOnFrequency(

    // k-mer length.
    size_t k,

    // The desired marker density
    double markerDensity,

    // Seed for random number generator.
    int seed,

    // Exclude k-mers enriched by more than this amount.
    // Enrichment is the ratio of k-mer frequency in reads
    // over what a random distribution would give.
    double enrichmentThreshold,

    size_t threadCount
)
{

    // Sanity check on the value of k, then store it.
    if(k > Kmer::capacity) {
        throw runtime_error("K-mer capacity exceeded.");
    }
    assemblerInfo->k = k;

    // Sanity check.
    if(markerDensity<0. || markerDensity>1.) {
        throw runtime_error("Invalid marker density " +
            to_string(markerDensity) + " requested.");
    }

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Fill in the fields of the k-mer table
    // that depends only on k.
    initializeKmerTable();

    // Compute the frequency of all k-mers in oriented reads.
    setupLoadBalancing(reads->readCount(), 1000);
    runThreads(&Assembler::computeKmerFrequency, threadCount);

    // Compute the total number of k-mer occurrences in reads
    // and the number of k-mers that can possibly occur.
    // This is done by counting all k-mers
    // when using the raw read representation and
    // only RLE k-mers when using the RLE read representation.
    uint64_t totalKmerOccurrences = 0;
    uint64_t possibleKmerCount = 0;
    for(uint64_t kmerId=0; kmerId!=kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        totalKmerOccurrences += info.frequency;
        if(assemblerInfo->readRepresentation == 0) {
            ++possibleKmerCount;
        } else {
            if(info.isRleKmer) {
                ++possibleKmerCount;
            }
        }
    }
    const double averageOccurrenceCount =
        double(totalKmerOccurrences) / double(possibleKmerCount);



    if(assemblerInfo->readRepresentation == 0) {

        // We are using raw read representation.
        cout <<
            "K-mer length k " << k << "\n"
            "Total number of k-mers " << kmerTable.size() << "\n"
            "Total number of k-mer occurrences in all oriented reads " << totalKmerOccurrences << "\n"
            "Average number of occurrences per k-mer " <<
            averageOccurrenceCount << endl;

    } else {

        // We are using RLE read representation.
        cout <<
            "K-mer length k " << k << "\n"
            "Total number of k-mers " << kmerTable.size() << "\n"
            "Total number of RLE k-mers " << possibleKmerCount << "\n"
            "Total number of k-mer occurrences in all oriented reads " << totalKmerOccurrences << "\n"
            "Average number of occurrences per RLE k-mer " <<
            averageOccurrenceCount << endl;
    }



    // Convert the enrichment threshold to a frequency.
    const uint64_t frequencyThreshold =
        uint64_t(enrichmentThreshold * averageOccurrenceCount);



    // Write out what we found.
    ofstream csv("KmerFrequencies.csv");
    csv << "KmerId,Kmer,ReverseComplementedKmerId,ReverseComplementedKmer,Frequency,Enrichment,Overenriched?\n";
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        const uint64_t frequency = info.frequency;

        const Kmer kmer(kmerId, k);
        const Kmer reverseComplementedKmer(info.reverseComplementedKmerId, k);
        csv << kmerId << ",";
        kmer.write(csv, k);
        csv << ",";
        reverseComplementedKmer.write(csv, k);
        csv << ",";
        csv << info.reverseComplementedKmerId << ",";
        csv << frequency << ",";
        csv << double(frequency) / averageOccurrenceCount;
        csv << ",";
        if(frequency > frequencyThreshold) {
            csv << "Yes";
        } else {
            csv << "No";
        }
        csv << "\n";
    }


    // Gather k-mers that are not overenriched.
    vector<KmerId> candidateKmers;
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        if((assemblerInfo->readRepresentation==1) and  (not info.isRleKmer)) {
            continue;
        }
        const uint64_t frequency = info.frequency;
        if(frequency > frequencyThreshold) {
            continue;
        }
        candidateKmers.push_back(KmerId(kmerId));
    }
    cout << possibleKmerCount - candidateKmers.size() << " k-mers were found to be "
        "enriched by more than a factor of " << enrichmentThreshold <<
        " and will not be used as markers." << endl;
    cout << "Markers will be chosen randomly among the remaining pool of " <<
        candidateKmers.size() << " k-mers." << endl;



    // Prepare to generate a random index into this vector of candidate k-mers.
    std::mt19937 randomSource(seed);
    std::uniform_int_distribution<uint64_t> uniformDistribution(0, candidateKmers.size()-1);

    // Flag all k-mers as not markers.
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        kmerTable[kmerId].isMarker = false;
    }


    // Now randomly generate markers from this pool of k-mers
    // until we have enough.
    uint64_t kmerOccurrencesCount = 0;
    uint64_t kmerCount = 0;
    const uint64_t desiredKmerOccurrencesCount =
        uint64_t(markerDensity * double(totalKmerOccurrences));
    while(kmerOccurrencesCount < desiredKmerOccurrencesCount) {

        // Generate a random index into the candidateKmers vector.
        const uint64_t index = uniformDistribution(randomSource);

        // Check that this k-mer is not already selected as a marker.
        const KmerId kmerId = candidateKmers[index];
        KmerInfo& info = kmerTable[kmerId];
        if(info.isMarker) {
            continue;
        }

        // This k-mer is not already selected as a marker.
        // Let's add it.
        info.isMarker = true;
        kmerOccurrencesCount += info.frequency;
        ++kmerCount;

        // If this k-mer is palindromic, we are done.
        if(info.reverseComplementedKmerId == kmerId) {
            continue;
        }

        // This k-mer is not palindromic, so we also add its reverse complement.
        KmerInfo& reverseComplementedInfo = kmerTable[info.reverseComplementedKmerId];
        SHASTA_ASSERT(!reverseComplementedInfo.isMarker);
        SHASTA_ASSERT(reverseComplementedInfo.frequency == info.frequency);
        reverseComplementedInfo.isMarker = true;
        kmerOccurrencesCount += reverseComplementedInfo.frequency;
        ++kmerCount;
    }
    cout << "Selected " << kmerCount << " k-mers as markers." << endl;


}


void Assembler::computeKmerFrequency(size_t threadId)
{
    // Create a frequency vector for this thread.
    MemoryMapped::Vector<uint64_t> frequency;
    frequency.createNew(
        largeDataName("tmp-KmerFrequency-" + to_string(threadId)),
        largeDataPageSize);
    frequency.resize(kmerTable.size());
    fill(frequency.begin(), frequency.end(), 0);



    // Loop over all batches assigned to this thread.
    const size_t k = assemblerInfo->k;
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            // Access the sequence of this read.
            const LongBaseSequenceView read = reads->getRead(readId);

            // If the read is pathologically short, it has no k-mers.
            if(read.baseCount < k) {
                continue;
            }

            // Loop over k-mers of this read.
            Kmer kmer;
            for(size_t position=0; position<k; position++) {
                kmer.set(position, read[position]);
            }
            for(uint32_t position=0; /*The check is done later */; position++) {

                // Get the k-mer id.
                const KmerId kmerId = KmerId(kmer.id(k));

                // Increment its frequency.
                ++frequency[kmerId];

                // Also increment the frequency of the reverse complemented k-mer.
                ++frequency[kmerTable[kmerId].reverseComplementedKmerId];

                // Check if we reached the end of the read.
                if(position+k == read.baseCount) {
                    break;
                }

                // Update the k-mer.
                kmer.shiftLeft();
                kmer.set(k-1, read[position+k]);
            }
        }
    }


    // Update the frequency in the k-mer table.
    {
        std::lock_guard<std::mutex> lock(mutex);
        for(uint64_t kmerId=0; kmerId!=frequency.size(); kmerId++) {
            kmerTable[kmerId].frequency += frequency[kmerId];
        }
    }



    // Remove the frequency vector for this thread.
    frequency.remove();
}



// Read the k-mers from file.
void Assembler::readKmersFromFile(uint64_t k, const string& fileName)
{
    // Sanity check on the value of k, then store it.
    if(k > Kmer::capacity) {
        throw runtime_error("K-mer capacity exceeded.");
    }
    assemblerInfo->k = k;

    // Fill in the fields of the k-mer table
    // that depends only on k.
    initializeKmerTable();

    // Open the file.
    ifstream file(fileName);
    if(not file) {
        throw runtime_error("Error opening " + fileName);
    }



    // Read one k-mer per line.
    uint64_t lineCount = 0;
    while(true) {

        // Read a line.
        string line;
        std::getline(file, line);
        if(not file) {
            break;
        }

        // Check the length.
        if(line.size() != k) {
            throw runtime_error("Unexpected line length in " + fileName + ":\n" + line + "\n" +
                "Expected " + to_string(k) + " characters, found " + to_string(line.size()));
        }

        // Read the k-mer.
        Kmer kmer;
        for(uint64_t i=0; i<k; i++) {
            const char c = line[i];
            const Base base = Base::fromCharacterNoException(c);
            if(not base.isValid()) {
                throw runtime_error("Unexpected base character in " + fileName + ":\n" + line);
            }
            kmer.set(i, base);
        }

        // Sanity checks.
        const KmerId kmerId = KmerId(kmer.id(k));
        SHASTA_ASSERT(kmerId < kmerTable.size());
        KmerInfo& kmerInfo = kmerTable[kmerId];
        if((assemblerInfo->readRepresentation==1) and (not kmerInfo.isRleKmer)) {
            throw runtime_error("Non-RLE k-mer (duplicate consecutive bases) in " +
                fileName + ":\n" + line);
        }

        // Flag it as a marker, together with its reverse complement.
        kmerInfo.isMarker = 1;
        kmerTable[kmerInfo.reverseComplementedKmerId].isMarker = 1;
        ++lineCount;

    }
    cout << "Processed " << lineCount << " lines of " << fileName << endl;

    // Count the number of k-mers flagged as markers.
    uint64_t usedKmerCount = 0;
    uint64_t possibleKmerCount = 0;
    for(const KmerInfo& kmerInfo: kmerTable) {
        if(kmerInfo.isMarker) {
            ++usedKmerCount;
        }
        if(assemblerInfo->readRepresentation == 0) {
            ++possibleKmerCount;
        } else {
            if(kmerInfo.isRleKmer) {
                ++possibleKmerCount;
            }
        }
    }
    cout << "Flagged as markers " << usedKmerCount << " out of " << possibleKmerCount <<
        " possible k-mers of length " << k << endl;
}



// In this version, marker k-mers are selected randomly, but excluding
// any k-mer that is over-enriched even in a single oriented read.
void Assembler::selectKmers2(

    // k-mer length.
    size_t k,

    // The desired marker density
    double markerDensity,

    // Seed for random number generator.
    int seed,

    // Exclude k-mers enriched by more than this amount,
    // even in a single oriented read.
    // Enrichment is the ratio of k-mer frequency in reads
    // over what a random distribution would give.
    double enrichmentThreshold,

    size_t threadCount
)
{
    // Sanity check on the value of k, then store it.
    if(k > Kmer::capacity) {
        throw runtime_error("K-mer capacity exceeded.");
    }
    assemblerInfo->k = k;

    // Sanity check.
    if(markerDensity<0. || markerDensity>1.) {
        throw runtime_error("Invalid marker density " +
            to_string(markerDensity) + " requested.");
    }

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Fill in the fields of the k-mer table
    // that depends only on k.
    initializeKmerTable();

    // Store the enrichmentThreshold so all threads can see it.
    selectKmers2Data.enrichmentThreshold = enrichmentThreshold;

    // For each KmerId that is an RLE k-mer, compute the
    // global frequency (total number of occurrences in all
    // oriented reads) and the number of reads in
    // which the k-mer is over-enriched.
    selectKmers2Data.globalFrequency.createNew(
        largeDataName("tmp-SelectKmers2-GlobalFrequency"),  largeDataPageSize);
    selectKmers2Data.overenrichedReadCount.createNew(
        largeDataName("tmp-SelectKmers2-OverenrichedReadCount"),  largeDataPageSize);
    selectKmers2Data.globalFrequency.resize(kmerTable.size());
    selectKmers2Data.overenrichedReadCount.resize(kmerTable.size());
    fill(
        selectKmers2Data.globalFrequency.begin(),
        selectKmers2Data.globalFrequency.end(), 0);
    fill(
        selectKmers2Data.overenrichedReadCount.begin(),
        selectKmers2Data.overenrichedReadCount.end(), 0);
    setupLoadBalancing(reads->readCount(), 100);
    runThreads(&Assembler::selectKmers2ThreadFunction, threadCount);



    // Compute the total number of k-mer occurrences
    // and the number of possible kmers.
    uint64_t totalKmerOccurrences = 0;
    uint64_t possibleKmerCount = 0;
    for(uint64_t kmerId=0; kmerId!=kmerTable.size(); kmerId++) {
        totalKmerOccurrences += selectKmers2Data.globalFrequency[kmerId];
        if(assemblerInfo->readRepresentation == 0) {
            ++ possibleKmerCount;
        } else {
            if(kmerTable[kmerId].isRleKmer) {
                ++possibleKmerCount;
            }
        }
    }
    const double averageOccurrenceCount =
        double(totalKmerOccurrences) / double(possibleKmerCount);



    // Write out what we found.
    ofstream csv("KmerFrequencies.csv");
    csv << "KmerId,Kmer,ReverseComplementedKmerId,ReverseComplementedKmer,"
        "GlobalFrequency,GlobalEnrichment,NumberOfReadsOverenriched\n";
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        const uint64_t frequency = selectKmers2Data.globalFrequency[kmerId];

        const Kmer kmer(kmerId, k);
        const Kmer reverseComplementedKmer(info.reverseComplementedKmerId, k);
        csv << kmerId << ",";
        kmer.write(csv, k);
        csv << ",";
        reverseComplementedKmer.write(csv, k);
        csv << ",";
        csv << info.reverseComplementedKmerId << ",";
        csv << frequency << ",";
        csv << double(frequency) / averageOccurrenceCount;
        csv << ",";
        csv << selectKmers2Data.overenrichedReadCount[kmerId];

        csv << "\n";
    }


    // Gather k-mers that are not overenriched in any read and therefore
    // can be used as markers.
    vector<KmerId> candidateKmers;
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        if(kmerTable[kmerId].isRleKmer and selectKmers2Data.overenrichedReadCount[kmerId] == 0) {
            candidateKmers.push_back(KmerId(kmerId));
        }
    }
    cout << "Out of a total " << possibleKmerCount << " possible k-mers, " <<
        possibleKmerCount - candidateKmers.size() <<
        " were found to be over-enriched by more than a factor of " <<
        enrichmentThreshold <<
        " in at least one read and will not be used as markers." << endl;
    cout << "Markers will be chosen randomly from the remaining pool of " <<
        candidateKmers.size() << " k-mers." << endl;
    cout << "The enrichment threshold of " << enrichmentThreshold <<
        " corresponds to one occurrence every " <<
        double(possibleKmerCount) / enrichmentThreshold <<
        " bases";
    if(assemblerInfo->readRepresentation == 1) {
        cout << " (in RLE representation)";
    }
    cout << "." << endl;


    // Prepare to generate a random index into this vector of candidate k-mers.
    std::mt19937 randomSource(seed);
    std::uniform_int_distribution<uint64_t> uniformDistribution(0, candidateKmers.size()-1);

    // Flag all k-mers as not markers.
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        kmerTable[kmerId].isMarker = false;
    }



    // Now randomly generate markers from this pool of k-mers
    // until we have enough.
    uint64_t kmerOccurrencesCount = 0;
    uint64_t kmerCount = 0;
    const uint64_t desiredKmerOccurrencesCount =
        uint64_t(markerDensity * double(totalKmerOccurrences));
    while(kmerOccurrencesCount < desiredKmerOccurrencesCount) {

        // Generate a random index into the candidateKmers vector.
        const uint64_t index = uniformDistribution(randomSource);

        // Check that this k-mer is not already selected as a marker.
        const KmerId kmerId = candidateKmers[index];
        KmerInfo& info = kmerTable[kmerId];
        if(info.isMarker) {
            continue;
        }

        // This k-mer is not already selected as a marker.
        // Let's add it.
        info.isMarker = true;
        kmerOccurrencesCount += selectKmers2Data.globalFrequency[kmerId];
        ++kmerCount;

        // If this k-mer is palindromic, we are done.
        if(info.reverseComplementedKmerId == kmerId) {
            continue;
        }

        // This k-mer is not palindromic, so we also add its reverse complement.
        KmerInfo& reverseComplementedInfo = kmerTable[info.reverseComplementedKmerId];
        SHASTA_ASSERT(!reverseComplementedInfo.isMarker);
        SHASTA_ASSERT(reverseComplementedInfo.frequency == info.frequency);
        reverseComplementedInfo.isMarker = true;
        kmerOccurrencesCount += selectKmers2Data.globalFrequency[info.reverseComplementedKmerId];
        ++kmerCount;
    }
    cout << "Selected " << kmerCount << " k-mers as markers." << endl;
    cout << "These k-mers have a total " << kmerOccurrencesCount <<
        " occurrences out of a total " << totalKmerOccurrences <<
        " in all oriented reads." << endl;

}



void Assembler::selectKmers2ThreadFunction(size_t threadId)
{
    // Initialize globalFrequency for this thread.
    MemoryMapped::Vector<uint64_t> globalFrequency;
    globalFrequency.createNew(
        largeDataName("tmp-SelectKmers2-GlobalFrequency-" + to_string(threadId)),
        largeDataPageSize);
    globalFrequency.resize(kmerTable.size());
    fill(globalFrequency.begin(), globalFrequency.end(), 0);

    // Initialize overenrichedReadCount for this thread.
    MemoryMapped::Vector<ReadId> overenrichedReadCount;
    overenrichedReadCount.createNew(
        largeDataName("tmp-SelectKmers2-OverenrichedReadCount-" + to_string(threadId)),
        largeDataPageSize);
    overenrichedReadCount.resize(kmerTable.size());
    fill(overenrichedReadCount.begin(), overenrichedReadCount.end(), 0);

    // Vectors to hold KmerIds and their frequencies for a single read.
    vector<KmerId> readKmerIds;
    vector<uint32_t> readKmerIdFrequencies;

    // Access the enrichmentThreshold.
    const double enrichmentThreshold = selectKmers2Data.enrichmentThreshold;

    // Compute the total number of possible k-mers.
    // It is needed below for overenrichment computations.
    uint64_t possibleKmerCount = 0;
    for(const KmerInfo& kmerInfo: kmerTable) {
        if(assemblerInfo->readRepresentation == 0) {
            ++possibleKmerCount;
        } else {
            if(kmerInfo.isRleKmer) {
                ++possibleKmerCount;
            }
        }
    }


    // Loop over all batches assigned to this thread.
    const size_t k = assemblerInfo->k;
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            // Access the sequence of this read.
            const LongBaseSequenceView read = reads->getRead(readId);

            // If the read is pathologically short, it has no k-mers.
            if(read.baseCount < k) {
                continue;
            }

            // Loop over k-mers of this read.
            Kmer kmer;
            for(size_t position=0; position<k; position++) {
                kmer.set(position, read[position]);
            }
            readKmerIds.clear();
            for(uint32_t position=0; /*The check is done later */; position++) {

                // Get the k-mer id.
                const KmerId kmerId = KmerId(kmer.id(k));
                readKmerIds.push_back(kmerId);

                // Increment its global frequency.
                ++globalFrequency[kmerId];

                // Also increment the frequency of the reverse complemented k-mer.
                ++globalFrequency[kmerTable[kmerId].reverseComplementedKmerId];

                // Check if we reached the end of the read.
                if(position+k == read.baseCount) {
                    break;
                }

                // Update the k-mer.
                kmer.shiftLeft();
                kmer.set(k-1, read[position+k]);
            }

            // Compute k-mer frequencies for this read.
            deduplicateAndCount(readKmerIds, readKmerIdFrequencies);

            // Compute the k-mer frequency threshold for an over-enriched k-mer
            // in this read.
            const uint32_t readKmerCount = uint32_t(read.baseCount + 1 - k);
            const uint32_t frequencyThreshold =
                uint32_t(enrichmentThreshold * double(readKmerCount) / double(possibleKmerCount));

            // See if any k-mers are over-enriched in this read.
            SHASTA_ASSERT(readKmerIds.size() == readKmerIdFrequencies.size());
            for(uint64_t i=0; i<readKmerIds.size(); i++) {
                const KmerId kmerId = readKmerIds[i];
                const uint32_t frequency = readKmerIdFrequencies[i];
                if(frequency > frequencyThreshold) {
                    ++overenrichedReadCount[kmerId];
                    ++overenrichedReadCount[kmerTable[kmerId].reverseComplementedKmerId];
                }
            }
        }
    }


    // Add our globalFrequency and overenrichedReadCount
    // to the values computer by the other threads.
    {
        std::lock_guard<std::mutex> lock(mutex);
        for(uint64_t kmerId=0; kmerId!=globalFrequency.size(); kmerId++) {
            selectKmers2Data.globalFrequency[kmerId] += globalFrequency[kmerId];
            selectKmers2Data.overenrichedReadCount[kmerId] += overenrichedReadCount[kmerId];
        }
    }
    globalFrequency.remove();
    overenrichedReadCount.remove();
}




// In this version, marker k-mers are selected randomly, but excluding
// k-mers that appear repeated at short distances in any oriented read.
// More precisely, for each k-mer we compute the minimum distance
// (in RLE bases) at which any two copies of that k-mer appear in any oriented read.
// K-mers for which this minimum distance is less than distanceThreshold
// are not used as markers. Marker k-mers are selected randomly among the
// remaining k-mers, until the desired marker density is achieved.
void Assembler::selectKmers4(

    // k-mer length.
    uint64_t k,

    // The desired marker density
    double markerDensity,

    // Seed for random number generator.
    uint64_t seed,

    // Exclude k-mers that appear in any read in two copies,
    // with the two copies closer than this distance (in RLE bases).
    uint64_t distanceThreshold,

    size_t threadCount
)
{
    const bool debug = false;
    cout << timestamp << "Begin selectKmers4." << endl;

    // Sanity check on the value of k, then store it.
    if(k > Kmer::capacity) {
        throw runtime_error("K-mer capacity exceeded.");
    }
    assemblerInfo->k = k;

    // Sanity check.
    if(markerDensity<0. || markerDensity>1.) {
        throw runtime_error("Invalid marker density " +
            to_string(markerDensity) + " requested.");
    }

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Fill in the fields of the k-mer table
    // that depends only on k.
    initializeKmerTable();

    // Initialize the global frequency of all k-mers.
    selectKmers4Data.globalFrequency.createNew(
        largeDataName("tmp-SelectKmers4-GlobalFrequency"),  largeDataPageSize);
    selectKmers4Data.globalFrequency.resize(kmerTable.size());
    fill(
        selectKmers4Data.globalFrequency.begin(),
        selectKmers4Data.globalFrequency.end(), 0);


    // Initialize the minimumDistance vector, which stores
    // the minimum RLE distance between any two copies of each k-mer
    // in any oriented read.
    selectKmers4Data.minimumDistance.createNew(
        largeDataName("tmp-selectKmers4-minimumDistance"), largeDataPageSize);
    const uint64_t kmerCount = kmerTable.size();
    selectKmers4Data.minimumDistance.resize(kmerCount);
    for(uint64_t i=0; i<kmerCount; i++) {
        selectKmers4Data.minimumDistance[i].second = std::numeric_limits<uint32_t>::max();
    }

    // Compute the minimumDistance vector.
    setupLoadBalancing(reads->readCount(), 100);
    runThreads(&Assembler::selectKmers4ThreadFunction, threadCount);



    // Write out what we found.
    if(debug) {
        const uint64_t totalFrequency = std::accumulate(
            selectKmers4Data.globalFrequency.begin(),
            selectKmers4Data.globalFrequency.end(), 0ULL);
        cout << "Total number of k-mer occurrences in all oriented reads is " << totalFrequency << endl;
        ofstream csv("KmerInfo.csv");
        csv << "KmerId,Kmer,KmerIdRc,KmerRc,Frequency,FrequencyRc,TotalFrequency,"
            "MinDist,MinDistRc,MinMinDist\n";
        for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
            const KmerInfo& info = kmerTable[kmerId];
            if(!info.isRleKmer) {
                continue;
            }

            const uint64_t frequency = selectKmers4Data.globalFrequency[kmerId];
            const uint64_t frequencyReverseComplement = selectKmers4Data.globalFrequency[info.reverseComplementedKmerId];
            const uint64_t totalFrequency = frequency + frequencyReverseComplement;

            const uint32_t minimumDistance = selectKmers4Data.minimumDistance[kmerId].second;
            const uint32_t minimumDistanceReverseComplement =
                selectKmers4Data.minimumDistance[info.reverseComplementedKmerId].second;

            const Kmer kmer(kmerId, k);
            const Kmer reverseComplementedKmer(info.reverseComplementedKmerId, k);
            csv << kmerId << ",";
            kmer.write(csv, k);
            csv << ",";
            csv << info.reverseComplementedKmerId << ",";
            reverseComplementedKmer.write(csv, k);
            csv << ",";
            csv << frequency << ",";
            csv << frequencyReverseComplement << ",";
            csv << totalFrequency << ",";
            csv << minimumDistance << ",";
            csv << minimumDistanceReverseComplement << ",";
            csv << min(minimumDistance, minimumDistanceReverseComplement) << "\n";
        }
    }



    // Compute the total number of k-mer occurrences
    // and the number of RLE kmers.
    uint64_t totalKmerOccurrences = 0;
    uint64_t rleKmerCount = 0;
    for(uint64_t kmerId=0; kmerId!=kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        if(not info.isRleKmer) {
            SHASTA_ASSERT(selectKmers4Data.globalFrequency[kmerId] == 0);
            continue;
        }
        totalKmerOccurrences += selectKmers4Data.globalFrequency[kmerId];
        if(kmerTable[kmerId].isRleKmer) {
            ++rleKmerCount;
        }
    }
    cout << "K-mer length k " << k << endl;
    cout << "Distance threshold " << distanceThreshold << " RLE bases." << endl;
    cout << "Total number of distinct RLE k-mers " << rleKmerCount << endl;
    cout << "Total number of RLE k-mers in all oriented reads " << totalKmerOccurrences << endl;
    cout << "Requested marker density " << markerDensity << endl;
    const uint64_t requiredMarkerOccurrences = uint64_t(markerDensity * double(totalKmerOccurrences));
    cout << "Required number of marker occurrences in all oriented reads " << requiredMarkerOccurrences << endl;



    // Gather k-mers for which the minimum distance between two copies
    // equals at least distanceThreshold. Exclude palindromic k-mers.
    vector<KmerId> candidateKmers;
    uint64_t candidateFrequency = 0;
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        const KmerId kmerIdRc = info.reverseComplementedKmerId;
        if(not info.isRleKmer) {
            continue;
        }
        if(kmerIdRc == kmerId) {
            // Palindromic. Exclude.
            continue;
        }
        // Only store the lower KmerId in the pair.
        if(kmerId > kmerIdRc) {
            continue;
        }
        if(selectKmers4Data.minimumDistance[kmerId].second < distanceThreshold) {
            // Too close. skip.
            continue;
        }
        if(selectKmers4Data.minimumDistance[kmerIdRc].second < distanceThreshold) {
            // Too close. Skip.
            continue;
        }

        candidateKmers.push_back(KmerId(kmerId));
        candidateFrequency += selectKmers4Data.globalFrequency[kmerId];
        candidateFrequency += selectKmers4Data.globalFrequency[kmerIdRc];
    }
    cout << "Markers will be chosen randomly from the a pool of " <<
        2*candidateKmers.size() << " RLE k-mers." << endl;
    cout << "RLE k-mers in this pool occur " <<
        candidateFrequency << " times in all oriented reads." << endl;
    cout << "This is sufficient to achieve marker density up to " <<
        double(candidateFrequency) / double(totalKmerOccurrences) << endl;

    // If these candidates don't have sufficient frequency, we
    // can't achieve the required marker density.
    if(candidateFrequency < requiredMarkerOccurrences) {
        throw runtime_error("Cannot achieve required marker density. "
            "Increase k, or decrease marker density, or decrease distance threshold.");
    }

    // Flag all k-mers as not markers.
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        kmerTable[kmerId].isMarker = false;
    }



    // Randomly pick markers in this vector of candidate k-mers.
    std::mt19937_64 randomSource(seed);
    std::uniform_real_distribution<> uniformDistribution;
    uint64_t markerOccurrencesCount = 0;
    uint64_t markerCount = 0;
    while(true) {

        // Pick a random index in the candidateKmers vector.
        const double x = uniformDistribution(randomSource);  // In [0,1)
        const uint64_t i = uint64_t(x * double(candidateKmers.size()));
        SHASTA_ASSERT(i < candidateKmers.size());

        // This KmerId and its reverse complement  will be used as markers.
        const KmerId kmerId = candidateKmers[i];
        const KmerId kmerIdRc = kmerTable[kmerId].reverseComplementedKmerId;

        kmerTable[kmerId].isMarker = true;
        kmerTable[kmerIdRc].isMarker = true;

        // Increment counters.
        markerCount += 2;
        markerOccurrencesCount += selectKmers4Data.globalFrequency[kmerId];
        markerOccurrencesCount += selectKmers4Data.globalFrequency[kmerIdRc];

        // Remove kmerId from the vector of candidates.
        if(i != candidateKmers.size()-1) {
            candidateKmers[i] = candidateKmers.back();
        }
        candidateKmers.resize(candidateKmers.size() - 1);

        if(markerOccurrencesCount >= requiredMarkerOccurrences) {
            break;
        }
    }
    cout << "Selected " << markerCount << " k-mers as markers." << endl;
    cout << "Actual marker density " << double(markerOccurrencesCount) / double(totalKmerOccurrences) << endl;



    // Clean up.
    selectKmers4Data.minimumDistance.remove();
    selectKmers4Data.globalFrequency.remove();

    // Done.
    cout << timestamp << "End selectKmers4." << endl;
}



void Assembler::selectKmers4ThreadFunction(size_t threadId)
{
    // K-mer length.
    const size_t k = assemblerInfo->k;

    // Initialize globalFrequency for this thread.
    // Having all threads accumulate atomically on the global frequency vector is too slow.
    MemoryMapped::Vector<uint64_t> globalFrequency;
    globalFrequency.createNew(
        largeDataName("tmp-SelectKmers4-GlobalFrequency-" + to_string(threadId)),
        largeDataPageSize);
    globalFrequency.resize(kmerTable.size());
    fill(globalFrequency.begin(), globalFrequency.end(), 0);

    // Vector to hold pairs(KmerId, RLE position) for one read.
    vector< pair<KmerId, uint32_t> > readKmers;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all reads assigned to this batch.
        for(ReadId readId=ReadId(begin); readId!=ReadId(end); readId++) {

            // Access the sequence of this read.
            const LongBaseSequenceView read = reads->getRead(readId);

            // If the read is pathologically short, it has no k-mers.
            if(read.baseCount < k) {
                continue;
            }

            // Loop over k-mers of this read.
            Kmer kmer;
            for(size_t position=0; position<k; position++) {
                kmer.set(position, read[position]);
            }
            readKmers.clear();
            for(uint32_t position=0; /*The check is done later */; position++) {

                // Get the k-mer id.
                const KmerId kmerId = KmerId(kmer.id(k));
                readKmers.push_back(make_pair(kmerId, position));

                // Update the frequency of this k-mer.
                ++globalFrequency[kmerId];
                ++globalFrequency[kmerTable[kmerId].reverseComplementedKmerId];

                // Check if we reached the end of the read.
                if(position+k == read.baseCount) {
                    break;
                }

                // Update the k-mer.
                kmer.shiftLeft();
                kmer.set(k-1, read[position+k]);
            }

            // Sort by k-mer, then by position.
            sort(readKmers.begin(), readKmers.end());

            // Update minDistance for each pair of repeated k-mers.
            for(uint64_t i=1; i<readKmers.size(); i++) {
                const auto& p0 = readKmers[i-1];
                const auto& p1 = readKmers[i];
                const KmerId kmerId0 = p0.first;
                const KmerId kmerId1 = p1.first;
                if(kmerId0 != kmerId1) {
                    continue;
                }
                const uint32_t distance = p1.second - p0.second;

                pair<std::mutex, uint32_t>& p = selectKmers4Data.minimumDistance[kmerId0];
                std::lock_guard<std::mutex> lock(p.first);;
                p.second = min(p.second, distance);
            }
        }
    }

    // Add our globalFrequency to the values computed by the other threads.
    {
        std::lock_guard<std::mutex> lock(mutex);
        for(uint64_t kmerId=0; kmerId!=globalFrequency.size(); kmerId++) {
            selectKmers4Data.globalFrequency[kmerId] += globalFrequency[kmerId];
        }
    }
    globalFrequency.remove();
}

