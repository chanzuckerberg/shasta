#include "Assembler.hpp"
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
    //of this they are chose with lower probability.
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
    size_t usedKmerCount = 0;
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        if(kmerTable[kmerId].isMarker) {
            ++usedKmerCount;
        }
    }
    cout << "Selected " << usedKmerCount << " " << k << "-mers as markers out of ";
    cout << kmerCount << " total." << endl;
    cout << "Requested inclusion probability: " << probability << "." << endl;
    cout << "Actual fraction of marker k-mers: ";
    cout << double(usedKmerCount)/double(kmerCount) << "." << endl;
    cout << "The above statistics include all k-mers, not just those present in "
        "run-length encoded sequence." << endl;

    if(probability == 1.) {
        SHASTA_ASSERT(usedKmerCount == kmerCount);
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
    setupLoadBalancing(reads.size(), 1000);
    runThreads(&Assembler::computeKmerFrequency, threadCount);

    // Compute the total number of k-mer occurrences
    // and the number of RLE kmers.
    uint64_t totalKmerOccurrences = 0;
    uint64_t rleKmerCount = 0;
    for(uint64_t kmerId=0; kmerId!=kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        totalKmerOccurrences += info.frequency;
        if(info.isRleKmer) {
            ++rleKmerCount;
        }
    }
    const double averageOccurrenceCount =
        double(totalKmerOccurrences) / double(rleKmerCount);
    cout <<
        "K-mer length k " << k << "\n"
        "Total number of k-mers " << kmerTable.size() << "\n"
        "Total number of RLE k-mers " << rleKmerCount << "\n"
        "Total number of k-mer occurrences " << totalKmerOccurrences << "\n"
        "Average number of occurrences per RLE k-mer " <<
        averageOccurrenceCount << endl;

    // Convert the enrichment threshold to a frequency.
    const uint64_t frequencyThreshold =
        uint64_t(enrichmentThreshold * averageOccurrenceCount);



    // Write out what we found.
    ofstream csv("KmerFrequencies.csv");
    csv << "KmerId,Kmer,ReverseComplementedKmerId,ReverseComplementedKmer,Frequency,Enrichment,Overenriched?\n";
    for(uint64_t kmerId=0; kmerId<kmerTable.size(); kmerId++) {
        const KmerInfo& info = kmerTable[kmerId];
        const uint64_t frequency = info.frequency;
        if(!info.isRleKmer) {
            SHASTA_ASSERT(frequency == 0);
            continue;
        }

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
        if(!info.isRleKmer) {
            continue;
        }
        const uint64_t frequency = info.frequency;
        if(frequency > frequencyThreshold) {
            continue;
        }
        candidateKmers.push_back(KmerId(kmerId));
    }
    cout << rleKmerCount - candidateKmers.size() << " k-mers were found to be "
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


    // Now randomly generate markers fromthis pool of k-mers
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
            const LongBaseSequenceView read = reads[readId];

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


