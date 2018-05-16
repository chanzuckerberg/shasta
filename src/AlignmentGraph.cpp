#include "AlignmentGraph.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;



// Ccompute an alignment of the markers of two oriented reads.
void ChanZuckerberg::Nanopore2::align(

    // Markers of the two oriented reads to be aligned, sorted by KmerId.
    const vector<Marker>& markers0,
    const vector<Marker>& markers1,

    // The maximum ordinal skip to be tolerated between successive markers
    // in the alignment.
    int maxSkip,

    // The AlignmentGraph can be reused.
    // For performance, it should be reused when doing many alignments.
    AlignmentGraph& graph,

    // Flag to control various types of debug output.
    bool debug,

    // The computed alignment.
    // This should also be reused when performance is important.
    Alignment& alignment
    )
{
    graph.create(markers0, markers1, maxSkip, debug, alignment);
}
