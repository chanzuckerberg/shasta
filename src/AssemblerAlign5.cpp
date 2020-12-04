#include "Assembler.hpp"
#include "Align5.hpp"
#include "html.hpp"
using namespace shasta;


// Python-callable version.
void Assembler::alignOrientedReads5(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1,
    uint64_t m,
    uint64_t deltaX,
    uint64_t deltaY,
    int64_t matchScore,
    int64_t mismatchScore,
    int64_t gapScore) const
{
    Align5::Options options;
    options.m = m;
    options.deltaX = deltaX;
    options.deltaY = deltaY;
    options.matchScore = matchScore;
    options.mismatchScore = mismatchScore;
    options.gapScore = gapScore;

    // Align5 work area.
    MemoryMapped::VectorOfVectors<Align5::MatrixEntry, uint64_t> matrix;
    matrix.createNew(largeDataName("tmp-Align5Matrix"), largeDataPageSize);

    Alignment alignment;
    AlignmentInfo alignmentInfo;

    const bool debug = true;

    // Compute the alignment.
    for(uint64_t i=0; i<3; i++) {
        cout << "Start computing alignment, attempt " << i << endl;
        alignOrientedReads5(
            OrientedReadId(readId0, strand0),
            OrientedReadId(readId1, strand1),
            options, matrix, alignment, alignmentInfo, debug);
    }

    matrix.remove();
}



// Align two reads using alignment method 5.
void Assembler::alignOrientedReads5(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    const Align5::Options& options,
    MemoryMapped::VectorOfVectors<Align5::MatrixEntry, uint64_t>& matrix,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug) const
{
    const auto markers0 = markers[orientedReadId0.getValue()];
    const auto markers1 = markers[orientedReadId1.getValue()];

    align5(markers0, markers1,
        options, matrix, alignment, alignmentInfo, debug);
}





