#include "Align4.hpp"
#include "html.hpp"
#include "PngImage.hpp"
using namespace shasta;

#include "algorithm.hpp"
#include "fstream.hpp"
#include <queue>



// Align two arbitrary sequences  using alignment method 4.
// If debug is true, detailed output to html is produced.
// Otherwise, html is not used.
void shasta::align4(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const Align4Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug)
{
    switch(options.m) {
    case 1:
        align4<1>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    case 2:
        align4<2>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    case 3:
        align4<3>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    case 4:
        align4<4>(markers0, markers1, options, alignment, alignmentInfo, debug);
        return;
    default:
        SHASTA_ASSERT(0);
    }
}



// Version templated on m, the number of markers that define
// a "feature" used in the alignment.
template<uint64_t m> void shasta::align4(
    const span<const CompressedMarker>& markers0,
    const span<const CompressedMarker>& markers1,
    const Align4Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug)
{
    Align4<m> graph(markers0, markers1,
        options, alignment, alignmentInfo,
        debug);
}



template<uint64_t m> shasta::Align4<m>::Align4(
    const Sequence& sequence0,
    const Sequence& sequence1,
    const Align4Options& options,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo,
    bool debug)
{
    // Check that we are in the templated version consistent with
    /// the options.
    SHASTA_ASSERT(options.m == m);

    if(debug) {
        cout << "Computing a marker alignment of two sequences with " <<
            sequence0.size() << " and " << sequence1.size() << " markers." << endl;
    }

    // Create the FeatureMap for sequence0.
    const uint64_t inverseLoadFactor = 2;
    FeatureMap featureMap0(inverseLoadFactor * sequence0.size());
    fillFeatureMap(sequence0, featureMap0);

    // Create the alignment matrix.
    const uint32_t nx = uint32_t(sequence0.size()-(m-1));
    const uint32_t ny = uint32_t(sequence1.size()-(m-1));
    fillAlignmentMatrix(featureMap0, sequence1,
        nx, ny,
        uint32_t(options.deltaX), uint32_t(options.deltaY));
    computeReachability(int32_t(options.deltaX), int32_t(options.deltaY));
    if(debug) {
        cout << "Before removing unreachable entries, the alignment matrix has " <<
            alignmentMatrix.size() << " entries." << endl;
    }
    removeUnreachable();
    if(debug) {
        cout << "After removing unreachable entries, the alignment matrix has " <<
            alignmentMatrix.size() << " entries." << endl;
    }

    if(debug) {
        writeMatrixCsv("Align4-Matrix.csv");
        writeMatrixPng(nx, ny, "Align4-Matrix.png");
    }

    // Find alignment candidates.
    findCandidateAlignments(int32_t(options.deltaX), int32_t(options.deltaY), debug);
    if(debug) {
        writeCandidateAlignmentsCsv("Align4-CandidateAlignments.csv");
    }
}



template<uint64_t m> void shasta::Align4<m>::fillFeatureMap(
    const Sequence& sequence,
    FeatureMap& featureMap)
{
    SHASTA_ASSERT(featureMap.empty());

    // If the sequence is too short, it has no features.
    if(sequence.size() < m) {
        return;
    }

    // Start with the first m kmerIds.
    Feature feature;
    uint64_t nextOrdinal = 0;
    for(uint64_t i=0; i<m; i++) {
        feature[i] = sequence[nextOrdinal++].kmerId;
    }

    // Add the features.
    for(uint32_t ordinal=0; ; ordinal++) {
        featureMap.insert(make_pair(feature, ordinal));

        // Check if done.
        if(nextOrdinal == sequence.size()) {
            break;
        }

        // Shift by one.
        for(uint64_t i=1; i<m; i++) {
            feature[i-1] = feature[i];
        }
        feature.back() = sequence[nextOrdinal++].kmerId;
    }

    if(sequence.size() >= m) {
        SHASTA_ASSERT(featureMap.size() == sequence.size() + 1 - m);
    } else {
        SHASTA_ASSERT(featureMap.empty());
    }
}



template<uint64_t m> void shasta::Align4<m>::fillAlignmentMatrix(
    const FeatureMap& featureMap0,
    const Sequence& sequence1,
    int32_t nx,
    int32_t ny,
    int32_t deltaX,
    int32_t deltaY)
{
    // Start with the first m kmerIds.
    Feature feature;
    uint64_t nextOrdinal1 = 0;
    for(uint64_t i=0; i<m; i++) {
        feature[i] = sequence1[nextOrdinal1++].kmerId;
    }

    // Add the features.
    for(int32_t y=0; ; y++) {

        // Look up this feature in the feature map for sequence 0.
        typename FeatureMap::const_iterator begin, end;
        tie(begin, end) = featureMap0.equal_range(feature);

        // Loop over all the hits. Each hit generates an
        // entry in the alignment matrix.
        for(auto it=begin; it!=end; ++it) {
            const int32_t x = it->second;
            const int32_t X = x + y;
            const int32_t Y = y + nx - 1 - x;
            const int32_t iX = X / deltaX;
            const int32_t iY = Y / deltaY;
            AlignmentMatrixEntry alignmentMatrixEntry;
            alignmentMatrixEntry.xy = make_pair(x, y);
            alignmentMatrixEntry.XY = make_pair(X, Y);
            alignmentMatrixEntry.isNearLeft = (x < deltaX) ? 1 : 0;
            alignmentMatrixEntry.isNearTop  = (y < deltaX) ? 1 : 0;
            alignmentMatrixEntry.isNearRight = (nx-1-x < deltaX) ? 1 : 0;
            alignmentMatrixEntry.isNearBottom = (ny-1-y < deltaX) ? 1 : 0;
            alignmentMatrixEntry.isForwardReachableFromTopOrLeft = 0;
            alignmentMatrixEntry.isBackwardReachableFromBottomOrRight = 0;
            alignmentMatrix.insert(make_pair(Coordinates(iX, iY), alignmentMatrixEntry));
        }

        // Check if done.
        if(nextOrdinal1 == sequence1.size()) {
            break;
        }

        // Shift by one.
        for(uint64_t i=1; i<m; i++) {
            feature[i-1] = feature[i];
        }
        feature.back() = sequence1[nextOrdinal1++].kmerId;
    }
}


template<uint64_t m> void shasta::Align4<m>::writeMatrixCsv(
    const string& fileName)
{
    ofstream csv(fileName);
    csv << "x,y,X,Y,iX,iY,Near left,Near top,Near right,Near bottom\n";

    for(const auto& p: alignmentMatrix) {
        const Coordinates& iXY = p.first;
        const AlignmentMatrixEntry& alignmentMatrixEntry = p.second;
        csv <<
            alignmentMatrixEntry.xy.first << "," <<
            alignmentMatrixEntry.xy.second << "," <<
            alignmentMatrixEntry.XY.first << "," <<
            alignmentMatrixEntry.XY.second << "," <<
            iXY.first << "," <<
            iXY.second << "," <<
            int(alignmentMatrixEntry.isNearLeft) << "," <<
            int(alignmentMatrixEntry.isNearTop) << "," <<
            int(alignmentMatrixEntry.isNearRight) << "," <<
            int(alignmentMatrixEntry.isNearBottom) << "\n";
    }
}



template<uint64_t m> void shasta::Align4<m>::writeMatrixPng(
    uint32_t nx, uint32_t ny,
    const string& fileName)
{
    // Create the image, which gets initialized to black.
    PngImage image(nx, ny);

    // Write a grid.
    vector<int> gridSpacing;
    vector< array<int, 3> > gridRgb;
    gridSpacing.push_back(   10); gridRgb.push_back({ 15,  15,  15});  // Grey
    gridSpacing.push_back(   50); gridRgb.push_back({ 30,  30,  30});  // Grey
    gridSpacing.push_back(  100); gridRgb.push_back({ 90,  90,  90});  // Grey
    gridSpacing.push_back(  500); gridRgb.push_back({160, 160, 160});  // Grey
    gridSpacing.push_back( 1000); gridRgb.push_back({255, 255, 255});  // White
    gridSpacing.push_back( 5000); gridRgb.push_back({255, 120, 255});  // Purple
    gridSpacing.push_back(10000); gridRgb.push_back({255, 255,  60});  // Yellow
    gridSpacing.push_back(50000); gridRgb.push_back({255, 255, 120});  // Yellow
    for(size_t i=0; i<gridSpacing.size(); i++) {
        const int spacing = gridSpacing[i];

        const array<int, 3>& rgb = gridRgb[i];
        for(uint32_t x=0; x<nx; x+=spacing) {
            for(uint32_t y=0; y<ny; y++) {
                image.setPixel(x, y, rgb[0], rgb[1], rgb[2]);
            }
        }
        for(uint32_t y=0; y<ny; y+=spacing) {
            for(uint32_t x=0; x<nx; x++) {
                image.setPixel(x, y, rgb[0], rgb[1], rgb[2]);
            }
        }
    }

    // Write the alignment matrix.
    for(const auto& p: alignmentMatrix) {
        const AlignmentMatrixEntry& entry = p.second;
        const Coordinates& coordinates = entry.xy;
        if(entry.isForwardReachableFromTopOrLeft and entry.isBackwardReachableFromBottomOrRight) {
            image.setPixel(coordinates.first, coordinates.second, 0, 255, 0);
        } else {
            image.setPixel(coordinates.first, coordinates.second, 255, 0, 0);
        }
    }

    // Write it out.
    image.write(fileName);
}



template<uint64_t m> void shasta::Align4<m>::clearDiscoveredFlags()
{
    for(auto& p: alignmentMatrix) {
        p.second.wasDiscovered = 0;
    }
}



// Compute the reachability flags in the alignment matrix
template<uint64_t m> void shasta::Align4<m>::computeReachability(
    int32_t deltaX,
    int32_t deltaY)
{
    using iterator = typename AlignmentMatrix::iterator;
    std::queue<iterator> q;
    vector<iterator> neighbors;



    // Do a forward BFS starting from all matrix entries near the
    // top or left of the alignment matrix.
    clearDiscoveredFlags();
    for(iterator it=alignmentMatrix.begin(); it!=alignmentMatrix.end(); ++it) {
        AlignmentMatrixEntry& entry = it->second;
        if(entry.isNearTop or entry.isNearLeft) {
            q.push(it);
            entry.wasDiscovered = 1;
            entry.isForwardReachableFromTopOrLeft = 1;
        }
    }
    while(not q.empty()) {
        auto it0 = q.front();
        q.pop();
        findAndFlagUndiscoveredChildren(it0, deltaX, deltaY, neighbors);
        for(const auto& it1: neighbors) {
            q.push(it1);
            AlignmentMatrixEntry& entry1 = it1->second;
            entry1.isForwardReachableFromTopOrLeft = 1;
        }
    }



    // Do a backward BFS starting from all matrix entries near the
    // top or left of the alignment matrix.
    SHASTA_ASSERT(q.empty());
    clearDiscoveredFlags();
    for(iterator it=alignmentMatrix.begin(); it!=alignmentMatrix.end(); ++it) {
        AlignmentMatrixEntry& entry = it->second;
        if(entry.isNearBottom or entry.isNearRight) {
            q.push(it);
            entry.wasDiscovered = 1;
            entry.isBackwardReachableFromBottomOrRight = 1;
        }
    }
    while(not q.empty()) {
        auto it0 = q.front();
        q.pop();
        findAndFlagUndiscoveredParents(it0, deltaX, deltaY, neighbors);
        for(const auto& it1: neighbors) {
            q.push(it1);
            AlignmentMatrixEntry& entry1 = it1->second;
            entry1.isBackwardReachableFromBottomOrRight = 1;
        }
    }

}



// Remove alignment matrix entries that are not reachable in both directions.
template<uint64_t m> void shasta::Align4<m>:: removeUnreachable()
{
    // The loop must be done with care to avoid using invalid iterators
    // pointing to removed elements.
    for(auto it=alignmentMatrix.begin(); it!=alignmentMatrix.end();) {

        const AlignmentMatrixEntry& entry = it->second;
        if(not(entry.isForwardReachableFromTopOrLeft and entry.isBackwardReachableFromBottomOrRight)) {
            auto jt = it;
            ++it;
            alignmentMatrix.erase(jt);
        } else {
            ++it;
        }
    }
}



// Do a BFS to find candidate alignments.
// See comments at the top of Align4.hpp for details.
template<uint64_t m> void shasta::Align4<m>::findCandidateAlignments(
    int32_t deltaX,
    int32_t deltaY,
    bool debug)
{
    using iterator = typename AlignmentMatrix::iterator;
    std::queue<iterator> q;
    vector<iterator> neighbors;
    vector<iterator> candidateAlignment;
    clearDiscoveredFlags();

    // Loop over possible starting vertices for the BFS.
    for(auto itStart=alignmentMatrix.begin(); itStart!=alignmentMatrix.end(); itStart++) {
        const AlignmentMatrixEntry& start = itStart->second;

        // If already discovered,skip.
        if(start.wasDiscovered) {
            continue;
        }

        // If not near the top of left, skip.
        if(not(start.isNearLeft or start.isNearTop)) {
            continue;
        }

        // The candidate alignment we find starting here is only
        // usable if it reaches at least one vertex near the
        // right or bottom of the alignment matrix.
        bool isUsable = false;
        candidateAlignment.clear();

        // Do a BFS starting here.
        SHASTA_ASSERT(q.empty());
        q.push(itStart);
        itStart->second.wasDiscovered = true;
        candidateAlignment.push_back(itStart);
        while(not q.empty()) {

            // Dequeue a vertex.
            auto it0 = q.front();
            q.pop();

            // Find its neighbors.
            findAndFlagUndiscoveredNeighbors(it0, deltaX, deltaY, neighbors);

            // Queue them.
            for(const auto& it1: neighbors) {
                q.push(it1);
                candidateAlignment.push_back(it1);
                const AlignmentMatrixEntry& entry1 = it1->second;
                if(entry1.isNearRight or entry1.isNearBottom) {
                    isUsable = true;
                }
            }
        }

        // If not usable, skip it.
        if(not isUsable) {
            continue;
        }


        // Store it.
        sort(candidateAlignment.begin(), candidateAlignment.end(), OrderByIncreasingXThenY());
        const iterator itFirst = candidateAlignment.front();
        const AlignmentMatrixEntry& firstEntry = itFirst->second;
        if(debug) {
        cout << "Candidate alignment " << candidateAlignments.size() << " start at " <<
            firstEntry.xy.first << " " <<
            firstEntry.xy.second << ", " <<
            candidateAlignment.size() << " alignment matrix entries." << endl;
        }

        candidateAlignments.resize(candidateAlignments.size() + 1);
        candidateAlignment.swap(candidateAlignments.back());
    }

    if(debug) {
        cout << "Found " << candidateAlignments.size() << " candidate alignments." << endl;
    }
}



template<uint64_t m> void shasta::Align4<m>::findAndFlagUndiscoveredNeighbors(
    typename AlignmentMatrix::iterator it0,
    int32_t deltaX,
    int32_t deltaY,
    vector<typename AlignmentMatrix::iterator>& neighbors)
{
    using iterator = typename AlignmentMatrix::iterator;
    neighbors.clear();

    const AlignmentMatrixEntry& entry0 = it0->second;
    const Coordinates& iXY0 = it0->first;
    const int32_t iX0 = iXY0.first;
    const int32_t iY0 = iXY0.second;

    // Loop over the 9 cells centered here.
    for(int32_t diX=-1; diX<=1; diX++) {
        const int iX1 = iX0 + diX;
        for(int32_t diY=-1; diY<=1; diY++) {
            const int iY1 = iY0 + diY;

            // Loop over all alignment matrix entries in this cell.
            iterator begin, end;
            tie(begin, end) = alignmentMatrix.equal_range(Coordinates(iX1, iY1));
            for(iterator it1=begin; it1!=end; ++it1) {
                if(it1 == it0) {
                    continue;
                }
                AlignmentMatrixEntry& entry1 = it1->second;
                if(entry1.wasDiscovered) {
                    continue;
                }
                if(not entry0.isNeighbor(entry1, deltaX, deltaY)) {
                    continue;
                }
                entry1.wasDiscovered = true;
                neighbors.push_back(it1);
            }
        }
    }
}



template<uint64_t m> void shasta::Align4<m>::findAndFlagUndiscoveredChildren(
    typename AlignmentMatrix::iterator it0,
    int32_t deltaX,
    int32_t deltaY,
    vector<typename AlignmentMatrix::iterator>& neighbors)
{
    using iterator = typename AlignmentMatrix::iterator;
    neighbors.clear();

    const AlignmentMatrixEntry& entry0 = it0->second;
    const Coordinates& iXY0 = it0->first;
    const int32_t iX0 = iXY0.first;
    const int32_t iY0 = iXY0.second;

    // Loop over the 6 cells that could contain children.
    for(int32_t diX=0; diX<=1; diX++) {
        const int iX1 = iX0 + diX;
        for(int32_t diY=-1; diY<=1; diY++) {
            const int iY1 = iY0 + diY;

            // Loop over all alignment matrix entries in this cell.
            iterator begin, end;
            tie(begin, end) = alignmentMatrix.equal_range(Coordinates(iX1, iY1));
            for(iterator it1=begin; it1!=end; ++it1) {
                if(it1 == it0) {
                    continue;
                }
                AlignmentMatrixEntry& entry1 = it1->second;
                if(entry1.wasDiscovered) {
                    continue;
                }
                if(not entry1.isChild(entry0, deltaX, deltaY)) {
                    continue;
                }
                entry1.wasDiscovered = true;
                neighbors.push_back(it1);
            }
        }
    }
}



template<uint64_t m> void shasta::Align4<m>::findAndFlagUndiscoveredParents(
    typename AlignmentMatrix::iterator it0,
    int32_t deltaX,
    int32_t deltaY,
    vector<typename AlignmentMatrix::iterator>& neighbors)
{
    using iterator = typename AlignmentMatrix::iterator;
    neighbors.clear();

    const AlignmentMatrixEntry& entry0 = it0->second;
    const Coordinates& iXY0 = it0->first;
    const int32_t iX0 = iXY0.first;
    const int32_t iY0 = iXY0.second;

    // Loop over the 6 cells that could contain parents.
    for(int32_t diX=-1; diX<=0; diX++) {
        const int iX1 = iX0 + diX;
        for(int32_t diY=-1; diY<=1; diY++) {
            const int iY1 = iY0 + diY;

            // Loop over all alignment matrix entries in this cell.
            iterator begin, end;
            tie(begin, end) = alignmentMatrix.equal_range(Coordinates(iX1, iY1));
            for(iterator it1=begin; it1!=end; ++it1) {
                if(it1 == it0) {
                    continue;
                }
                AlignmentMatrixEntry& entry1 = it1->second;
                if(entry1.wasDiscovered) {
                    continue;
                }
                if(not entry1.isParent(entry0, deltaX, deltaY)) {
                    continue;
                }
                entry1.wasDiscovered = true;
                neighbors.push_back(it1);
            }
        }
    }
}



template<uint64_t m> void shasta::Align4<m>::writeCandidateAlignmentsCsv(
    const string& fileName)
{
    ofstream csv(fileName);
    csv << "Candidate,x,y,X,Y\n";
    for(uint64_t i=0; i<candidateAlignments.size(); i++) {
        const auto& candidateAlignment = candidateAlignments[i];
        for(const typename AlignmentMatrix::const_iterator it: candidateAlignment) {
            const AlignmentMatrixEntry& entry = it->second;
            csv <<
                i << "," <<
                entry.xy.first << "," <<
                entry.xy.second << "," <<
                entry.XY.first << "," <<
                entry.XY.second << "\n";
        }
    }
}



template<uint64_t m> bool shasta::Align4<m>::AlignmentMatrixEntry::isNeighbor(
    const AlignmentMatrixEntry& entry1,
    int32_t deltaX, int32_t deltaY) const
{
    const AlignmentMatrixEntry& entry0 = *this;
    const Coordinates& XY0 = entry0.XY;
    const Coordinates& XY1 = entry1.XY;
    const int32_t X0 = XY0.first;
    const int32_t X1 = XY1.first;
    if(abs(X0-X1) > deltaX) {
        return false;
    }
    const int32_t Y0 = XY0.second;
    const int32_t Y1 = XY1.second;
    if(abs(Y0-Y1) > deltaY) {
        return false;
    }
    return true;
}



// Return true if (*this) is a child of that.
template<uint64_t m> bool shasta::Align4<m>::AlignmentMatrixEntry::isChild(
    const AlignmentMatrixEntry& that,
    int32_t deltaX, int32_t deltaY) const
{
    const Coordinates& thisXY = XY;
    const Coordinates& thatXY = that.XY;

    const int32_t thisX = thisXY.first;
    const int32_t thatX = thatXY.first;
    if(thisX <= thatX) {
        return false;
    }
    if(thisX >= thatX + deltaX) {
        return false;
    }

    const int32_t thisY = thisXY.second;
    const int32_t thatY = thatXY.second;
    if(thisY <= thatY - deltaY) {
        return false;
    }
    if(thisY >= thatY + deltaY) {
        return false;
    }

    return true;
}



// Return true if (*this) is a parent of that.
template<uint64_t m> bool shasta::Align4<m>::AlignmentMatrixEntry::isParent(
    const AlignmentMatrixEntry& that,
    int32_t deltaX, int32_t deltaY) const
{
    const Coordinates& thisXY = XY;
    const Coordinates& thatXY = that.XY;

    const int32_t thisX = thisXY.first;
    const int32_t thatX = thatXY.first;
    if(thisX > thatX) {
        return false;
    }
    if(thisX <= thatX - deltaX) {
        return false;
    }

    const int32_t thisY = thisXY.second;
    const int32_t thatY = thatXY.second;
    if(thisY <= thatY - deltaY) {
        return false;
    }
    if(thisY >= thatY + deltaY) {
        return false;
    }

    return true;
}

