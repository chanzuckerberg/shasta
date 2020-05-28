// shasta.
#include "Assembler.hpp"
#include "findMarkerId.hpp"
#include "MarkerFinder.hpp"
using namespace shasta;



void Assembler::findMarkers(size_t threadCount)
{
    checkReadsAreOpen();
    checkKmersAreOpen();

    markers.createNew(largeDataName("Markers"), largeDataPageSize);
    MarkerFinder markerFinder(
        assemblerInfo->k,
        kmerTable,
        reads,
        markers,
        threadCount);

}



void Assembler::accessMarkers()
{
    markers.accessExistingReadOnly(largeDataName("Markers"));
}

void Assembler::checkMarkersAreOpen() const
{
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
    csv << "MarkerId,Ordinal,KmerId,Kmer,Position\n";
    for(uint32_t ordinal=0; ordinal<orientedReadMarkers.size(); ordinal++) {
        const CompressedMarker& marker = orientedReadMarkers[ordinal];
        const MarkerId markerId = getMarkerId(orientedReadId, ordinal);
        csv << markerId << ",";
        csv << ordinal << ",";
        csv << marker.kmerId << ",";
        csv << Kmer(marker.kmerId, assemblerInfo->k) << ",";
        csv << marker.position << "\n";
    }
}



vector<KmerId> Assembler::getMarkers(ReadId readId, Strand strand)
{
    const OrientedReadId orientedReadId(readId, strand);
    const auto orientedReadMarkers = markers[orientedReadId.getValue()];

    vector<KmerId> v;
    for(const CompressedMarker& marker: orientedReadMarkers) {
        v.push_back(marker.kmerId);
    }
    return v;
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
MarkerId Assembler::getMarkerId(
    OrientedReadId orientedReadId, uint32_t ordinal) const
{
    return
        (markers.begin(orientedReadId.getValue()) - markers.begin())
        + ordinal;
}



// Inverse of the above: given a global marker id,
// return its OrientedReadId and ordinal.
// This requires a binary search in the markers toc.
// This could be avoided, at the cost of storing
// an additional 4 bytes per marker.
pair<OrientedReadId, uint32_t>
    Assembler::findMarkerId(MarkerId markerId) const
{
    return shasta::findMarkerId(markerId, markers);
}



// Given a MarkerId, compute the MarkerId of the
// reverse complemented marker.
MarkerId Assembler::findReverseComplement(MarkerId markerId) const
{
	// Find the oriented read id and marker ordinal.
	OrientedReadId orientedReadId;
	uint32_t ordinal;
	tie(orientedReadId, ordinal) = findMarkerId(markerId);

	// Reverse complement.
	ordinal = uint32_t(markers.size(orientedReadId.getValue()) - 1 - ordinal);
	orientedReadId.flipStrand();

	// Return the corresponding Markerid.
	return getMarkerId(orientedReadId, ordinal);
}



// Write the frequency of markers in oriented reads.
void Assembler::writeMarkerFrequency()
{
    const uint64_t k = assemblerInfo->k;
    const uint64_t kmerCount = 1ULL << (2ULL*k);
    SHASTA_ASSERT(markers.isOpen());
    vector<uint64_t> frequency(kmerCount, 0);

    const CompressedMarker* compressedMarker = markers.begin();
    const CompressedMarker* end = markers.end();
    for(; compressedMarker!=end; ++compressedMarker) {
        ++frequency[compressedMarker->kmerId];
    }

    ofstream csv("MarkerFrequency.csv");
    for(uint64_t kmerId=0; kmerId<kmerCount; kmerId++) {
        const uint64_t n = frequency[kmerId];
        if(n== 0) {
            continue;
        }
        const Kmer kmer(kmerId, k);
        kmer.write(csv, k);
        csv << "," << n << "\n";
    }
}
