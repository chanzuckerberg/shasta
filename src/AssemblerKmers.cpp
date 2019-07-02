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

    // Sanity check on the requested fraction.
    // It can be 1 at most. If it is 1, all k-mers
    // are guaranteed to be selected (because we use <=
    // comparison in the main loop.
    if(probability<0. || probability>1.) {
        throw runtime_error("Invalid k-mer probability " +
            to_string(probability) + " requested.");
    }



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



    // Create the kmer table with the necessary size.
    kmerTable.createNew(largeDataName("Kmers"), largeDataPageSize);
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


