#include "Coverage.hpp"
using namespace ChanZuckerberg;
using namespace shasta;



// CoverageData constructor.
// If the base is '-', repeatCount must be zero.
// Otherwise, it must not be zero.
CoverageData::CoverageData(AlignedBase base, Strand strand, size_t repeatCount) :
    base(base), strand(strand), repeatCount(repeatCount)
{
    if(base.isGap()) {
        CZI_ASSERT(repeatCount == 0);
    } else {
        CZI_ASSERT(repeatCount > 0);
    }
}



// Coverage default constructor.
Coverage::Coverage()
{
    // Zero out the baseCoverage.
    for(size_t base=0; base<5; base++) {
        auto& v = baseCoverage[base];
        v[0] = 0;
        v[1] = 0;
    }

    // The detailedCoverage structure does not need to be zero,
    // as each entry is an empty vector.
}



// Add information about a supporting read.
// If the AlignedBase is '-',repeatCount must be zero.
// Otherwise, it must not be zero.
// This is the only public non-const function.
void Coverage::addRead(AlignedBase base, Strand strand, size_t repeatCount)
{
    // Sanity check on the base.
    const size_t baseValue = base.value;
    CZI_ASSERT(baseValue < 5);

    // Sanity check on the repeat count.
    if(base.isGap()) {
        CZI_ASSERT(repeatCount == 0);
    } else {
        CZI_ASSERT(repeatCount > 0);
    }

    // Store a CoverageData for this read.
    readCoverageData.push_back(CoverageData(base, strand, repeatCount));

    // Increment detailed coverage.
    auto& c = detailedCoverage[baseValue][strand];
    if(c.size() <= repeatCount) {
        c.resize(repeatCount + 1);
    }
    ++c[repeatCount];

    // Increment base coverage.
    ++baseCoverage[baseValue][strand];

}



// Return the base with the most coverage.
// This can return ACGT or '-'.
AlignedBase Coverage::mostFrequentBase() const
{
    uint8_t bestBaseValue = 4;
    size_t bestBaseCoverage = 0;

    for(uint8_t baseValue=0; baseValue<5; baseValue++) {
        const size_t baseCoverage = coverage(AlignedBase::fromInteger(baseValue));
        if(baseCoverage > bestBaseCoverage) {
            bestBaseValue = baseValue;
            bestBaseCoverage = baseCoverage;
        }
    }
    return AlignedBase::fromInteger(bestBaseValue);
}



// Get the repeat count with the most coverage for a given base.
// Note that, if the base is '-', this will always return 0.
size_t Coverage::mostFrequentRepeatCount(AlignedBase base) const
{
    size_t countEnd =repeatCountEnd(base);
    size_t bestCount = 0;
    size_t bestCountCoverage = 0;

    for(size_t repeatCount=0; repeatCount<countEnd; repeatCount++) {
        const size_t coverageForRepeatCount = coverage(base, repeatCount);

        // Changed ">" to ">=" here to break ties in favor
        // of the longer alternative.
        // This improves sequence identity a bit.
        // but it is still biased in favor of deletions.
        if(coverageForRepeatCount >= bestCountCoverage) {
            bestCount = repeatCount;
            bestCountCoverage = coverageForRepeatCount;
        }
    }

    return bestCount;
}



// Get the repeat count with the most coverage for the base
// with the most coverage.
// The should only be called if the base with the best coverage
// is not '-'.
size_t Coverage::mostFrequentBaseMostFrequentRepeatCount() const
{
    return mostFrequentRepeatCount(mostFrequentBase());
}



// Represent a coverage value with a single character.
char Coverage::coverageCharacter(size_t coverage)
{
    if(coverage == 0) {
        return '.';
    } else if(coverage < 10) {
        const string coverageString = to_string(coverage);
        CZI_ASSERT(coverageString.size() == 1);
        return coverageString[0];
    } else {
        return '*';
    }
}



// Get coverage for a given base, for all repeat counts,
// summing over both strands.
size_t Coverage::coverage(AlignedBase base) const
{
    // Extract the base value and check it.
    const uint8_t baseValue = base.value;
    CZI_ASSERT(baseValue < 5);

    // Return total coverage for this base,
    // for both strands and all repeat counts.
    return baseCoverage[baseValue][0] + baseCoverage[baseValue][1];

}



// Same as above, but return a single character representing coverage.
char Coverage::coverageCharacter(AlignedBase base) const
{
    return coverageCharacter(coverage(base));
}



// Get coverage for a given base and repeat count,
// summing over both strands.
size_t Coverage::coverage(AlignedBase base, size_t repeatCount) const
{
    // Extract the base value and check it.
    const uint8_t baseValue = base.value;
    CZI_ASSERT(baseValue < 5);

    // Access the coverage vector for this base.
    const auto& baseDetailedCoverage = detailedCoverage[baseValue];

    // Return coverage for the given repeat count, summing over both strands.
    size_t c = 0;
    for(Strand strand=0; strand<2; strand++) {
        const auto& baseAndStrandDetailedCoverage = baseDetailedCoverage[strand];
        if(repeatCount < baseAndStrandDetailedCoverage.size()) {
            c += baseAndStrandDetailedCoverage[repeatCount];
        }
    }
    return c;
}



// Same as above, but return a single character representing coverage.
char Coverage::coverageCharacter(AlignedBase base, size_t repeatCount) const
{
    return coverageCharacter(coverage(base, repeatCount));
}



// Get base coverage for the best base.
size_t Coverage::mostFrequentBaseCoverage() const
{
    return coverage(mostFrequentBase());
}



// Same as above, but return a single character representing coverage.
char Coverage::mostFrequentBaseCoverageCharacter() const
{
    return coverageCharacter(mostFrequentBaseCoverage());
}



// Get, for a given base, the first repeat count for which
// coverage becomes permanently zero.
// This can be used to loop over repeat coutns for that base.
// Note that, if the base is '-', this will always return 0.
size_t Coverage::repeatCountEnd(AlignedBase base) const
{
    const size_t baseValue = base.value;
    CZI_ASSERT(baseValue < 5);

    const auto& c = detailedCoverage[baseValue];
    return max(c[0].size(), c[1].size());
}



// Get a vector of CompressedCoverageData.
void Coverage::count(vector<CompressedCoverageData>& compressedCoverageData) const
{
    compressedCoverageData.clear();
    for(uint8_t base=0; base<5; base++) {
        const array<vector<size_t>, 2>& x = detailedCoverage[base];
        for(uint8_t strand=0; strand<2; strand++) {
            const vector<size_t>& y = x[strand];
            for(size_t repeatCount=0; repeatCount<y.size(); repeatCount++) {
                const size_t frequency = y[repeatCount];
                if(frequency > 0) {
                    CompressedCoverageData c;
                    c.base = base & 7;
                    c.strand = strand & 1;
                    c.repeatCount = uint8_t(min(size_t(255), repeatCount));
                    c.frequency = uint8_t(min(size_t(255), frequency));
                    compressedCoverageData.push_back(c);
                }
            }
        }
    }
}



