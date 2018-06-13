// Nanopore2.
#include "Assembler.hpp"
#include "MarkerFinder.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;



void Assembler::findMarkers(size_t threadCount)
{
    checkReadsAreOpen();
    checkKmersAreOpen();

    markers0.createNew(largeDataName("Markers0"), largeDataPageSize);
    markers.createNew(largeDataName("Markers"), largeDataPageSize);
    MarkerFinder markerFinder(
        assemblerInfo->k,
        kmerTable,
        reads,
        markers0,
        markers,
        threadCount);

}



void Assembler::accessMarkers()
{
    markers0.accessExistingReadOnly(largeDataName("Markers0"));
    markers.accessExistingReadOnly(largeDataName("Markers"));
}

void Assembler::checkMarkersAreOpen() const
{
    if(!markers0.isOpen()) {
        throw runtime_error("Markers are not accessible.");
    }
    if(!markers.isOpen()) {
        throw runtime_error("Markers are not accessible.");
    }
}


void Assembler::writeMarkers(ReadId readId, Strand strand, const string& fileName)
{
    // Check that we have what we need.
    checkKmersAreOpen();
    checkReadsAreOpen();
    checkMarkersAreOpen();
    checkReadId(readId);

    // Get the markers.
    const OrientedReadId orientedReadId(readId, strand);
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];

    // Write them out.
    ofstream csv(fileName);
    csv << "GlobalMarkerId,GlobalOrientedMarkerId,Ordinal,KmerId,Kmer,Position\n";
    for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
        const CompressedMarker& marker = orientedReadMarkers[ordinal];
        const OrientedMarkerId orientedMarkerId =
            getGlobalOrientedMarkerId(orientedReadId, ordinal);
        CZI_ASSERT(orientedMarkerId.getStrand() == strand);
        csv << orientedMarkerId.getMarkerId() << ",";
        csv << orientedMarkerId.getValue() << ",";
        csv << ordinal << ",";
        csv << marker.kmerId << ",";
        csv << Kmer(marker.kmerId, assemblerInfo->k) << ",";
        csv << marker.position << "\n";
    }
}



// Get markers sorted by KmerId for a given OrientedReadId.
void Assembler::getMarkersSortedByKmerId(
    OrientedReadId orientedReadId,
    vector<MarkerWithOrdinal>& markersSortedByKmerId) const
{
    const auto compressedMarkers = markers[orientedReadId.getValue()];
    markersSortedByKmerId.clear();
    markersSortedByKmerId.resize(compressedMarkers.size());

    for(uint32_t ordinal=0; ordinal<compressedMarkers.size(); ordinal++) {
        const CompressedMarker& compressedMarker = compressedMarkers[ordinal];
        markersSortedByKmerId[ordinal] = MarkerWithOrdinal(compressedMarker, ordinal);
    }

    // Sort by kmerId.
    sort(markersSortedByKmerId.begin(), markersSortedByKmerId.end());
}



// Given a marker by its OrientedReadId and ordinal,
// return the corresponding global marker id.
OrientedMarkerId Assembler::getGlobalOrientedMarkerId(
    OrientedReadId orientedReadId, uint32_t ordinal) const
{

    if(orientedReadId.getStrand() == 0) {
        // Count forward from the first marker of the read.
        return OrientedMarkerId(
            (markers0.begin(orientedReadId.getReadId()) - markers0.begin()) + ordinal,
            0);

    } else {
        // Count backwards from the last marker of the read
        // (which is the marker preceding the first marker of the next read).
        return OrientedMarkerId(
            (markers0.begin(orientedReadId.getReadId()+1ULL) - markers0.begin()) - 1ULL - ordinal,
            1);
    }

}


// Inverse of the above: given a global marker id,
// return its OrientedReadId and ordinal.
// This requires a binary search in the markers toc.
// This could be avoided, at the cost of storing
// an additional 4 bytes per marker.
pair<OrientedReadId, uint32_t>
    Assembler::findGlobalOrientedMarkerId(OrientedMarkerId orientedMarkerId) const
{
    const MarkerId markerId = orientedMarkerId.getMarkerId();
    const Strand strand = Strand(orientedMarkerId.getStrand());
    ReadId readId;
    uint32_t ordinal;
    tie(readId, ordinal) = markers0.find(markerId);

    if(strand == 1) {
        ordinal = uint32_t(markers0.size(readId)) - 1 - ordinal;
    }

    return make_pair(OrientedReadId(readId, strand), ordinal);
}

