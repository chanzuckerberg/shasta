#ifndef SHASTA_COVERAGE_HPP
#define SHASTA_COVERAGE_HPP

// Shasta.
#include "Base.hpp"
#include "SHASTA_ASSERT.hpp"
#include "ReadId.hpp"

// Standard library.
#include "algorithm.hpp"
#include "array.hpp"
#include "vector.hpp"



/*******************************************************************************

Class CoverageData stores coverage information for a single
read at a single position of a multiple sequence alignment.

Class Coverage stores coverage information for all reads at a single
position of a multiple sequence alignment.

*******************************************************************************/



namespace ChanZuckerberg {
    namespace shasta {
        class Coverage;
        class CoverageData;
        class CompressedCoverageData;
    }
}



// Class CoverageData stores coverage information for a single
// read at a single position of a multiple sequence alignment.
class ChanZuckerberg::shasta::CoverageData {
public:
    AlignedBase base;   // ACGT or "-" for a gap
    Strand strand;      // 0 for + strand or 1 for - strand.
    size_t repeatCount; // The repeat count found in this read.

    // Constructor.
    // If the base is '-', repeatCount must be zero.
    // Otherwise, it must not be zero.
    CoverageData(AlignedBase base, Strand strand, size_t repeatCount);
};



class ChanZuckerberg::shasta::CompressedCoverageData {
public:
    uint8_t base: 4;     // Used to code the AlignedBase.
    uint8_t strand: 4;
    uint8_t repeatCount; // Clipped at 255.
    uint8_t frequency;   // How many times present. Clipped at 255.

    // Python-callable accessors.
    char getBase() const
    {
        const AlignedBase alignedBase = AlignedBase::fromInteger(base);
        if(alignedBase.isGap()) {
            return '_';
        } else {
            return alignedBase.character();
        }
    }
    char getStrand() const
    {
        if(strand == 0) {
            return '+';
        } else if(strand == 1) {
            return '-';
        } else {
            SHASTA_ASSERT(0);
        }
    }
    int getRepeatCount() const
    {
        return repeatCount;
    }
    int getFrequency() const
    {
        return frequency;
    }
};
static_assert(sizeof(ChanZuckerberg::shasta::CompressedCoverageData) == 3*sizeof(uint8_t),
    "Unexpected size of CompressedCoverageData");



// Class Coverage stores coverage information for all reads at a single
// position of a multiple sequence alignment.
class ChanZuckerberg::shasta::Coverage {
public:

    // Default constructor.
    Coverage();

    // Add information about a supporting read.
    // If the AlignedBase is '-',repeatCount must be zero.
    // Otherwise, it must not be zero.
    // This is the only public non-const function.
    void addRead(AlignedBase, Strand, size_t repeatCount);

    // Return the list detailing coverage from each read.
    const vector<CoverageData>& getReadCoverageData() const
    {
        return readCoverageData;
    }


    // Return the base with the most coverage.
    // This can return ACGT or '-'.
    AlignedBase mostFrequentBase() const;

    // Get the repeat count with the most coverage for a given base.
    size_t mostFrequentRepeatCount(AlignedBase) const;

    // Get the repeat count with the most coverage for the base
    // with the most coverage.
    size_t mostFrequentBaseMostFrequentRepeatCount() const;



    // Represent a coverage value with a single character.
    static char coverageCharacter(size_t);

    // Get coverage for a given base,
    // summing over all repeats and both strands.
    size_t coverage(AlignedBase) const;
    char coverageCharacter(AlignedBase) const;

    // Get coverage for a given base and repeat count,
    // summing over both strands.
    size_t coverage(AlignedBase, size_t repeatCount) const;
    char coverageCharacter(AlignedBase, size_t repeatCount) const;

    // Get base coverage for the best base.
    size_t mostFrequentBaseCoverage() const;
    char mostFrequentBaseCoverageCharacter() const;

    // Get, for a given base, the first repeat count for which
    // coverage becomes permanently zero.
    // This can be used to loop over repeat counts for that base.
    // Note that, if the base is '-', this will always return 0.
    size_t repeatCountEnd(AlignedBase) const;

    // Get a vector of CompressedCoverageData.
    void count(vector<CompressedCoverageData>&) const;

private:

    // An entry for each read in the alignment.
    vector<CoverageData> readCoverageData;

    // Coverage for each base (ACGT or '-'), strand, and repeat count.
    // Indexed by [AlignedBase::value][strand][repeatCount].
    array< array<vector<size_t>, 2>, 5> detailedCoverage;

    // Coverage for each base (ACGT or '-') and strand.
    // Indexed by [AlignedBase::value][strand].
    // This contains the same information in repeatCountCoverage,
    // summed over all repeat counts.
    array< array<size_t, 2>, 5> baseCoverage;

};



#endif

