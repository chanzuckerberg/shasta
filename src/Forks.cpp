#include "Forks.hpp"
#include "html.hpp"
#include "orderPairs.hpp"
using namespace shasta;

#include <set>


Forks::Forks(
    const MemoryMapped::VectorOfVectors<CompressedMarker, uint64_t>& markers,
    const MarkerGraph& markerGraph) :
    markers(markers),
    markerGraph(markerGraph)
{

    // Sanity checks.
    SHASTA_ASSERT(markerGraph.edges.isOpen);
    SHASTA_ASSERT(markerGraph.edgesBySource.isOpen());
    SHASTA_ASSERT(markerGraph.edgesByTarget.isOpen());
    const VertexId vertexCount = markerGraph.edgesBySource.size();
    SHASTA_ASSERT(markerGraph.edgesByTarget.size() == vertexCount);


    // Create the forks.
    for(VertexId vertexId=0; vertexId<vertexCount; vertexId++) {
        if(markerGraph.edgesBySource.size(vertexId) > 1) {
            forks.push_back(Fork(vertexId, ForkDirection::Forward, markerGraph));
        }
        if(markerGraph.edgesByTarget.size(vertexId) > 1) {
            forks.push_back(Fork(vertexId, ForkDirection::Backward, markerGraph));
        }
    }

    // Create the data structures that, given a Fork, allow us to locate
    // other nearby Forks.
    constructEdgeInfos();
    constructOrientedReadEdgeInfos();

}



Forks::Fork::Fork(
    VertexId vertexId,
    ForkDirection direction,
    const MarkerGraph& markerGraph) :
    vertexId(vertexId),
    direction(direction)
{
    // Access the edges of this Fork.
    // These are the out-edges or in-edges of vertexId,
    // depending on the Fork direction.
    const span<const Uint40> edgeIds =
        (direction == ForkDirection::Forward) ?
            markerGraph.edgesBySource[vertexId] :
            markerGraph.edgesByTarget[vertexId];

    // Each edge generates a branch.
    for(const Uint40 edgeId: edgeIds) {
        branches.push_back(Branch(EdgeId(edgeId), markerGraph));
    }

    // Find the branches that each OrientedReadId is on.
    // Normally each OrientedReadId is on only one branch.
    // Note we can assume that each OrientedReadId appears
    // no more than once on each Branch because the Branch
    // constructor has already taken care of that.
    std::map<OrientedReadId, vector<uint64_t> > orientedReadBranches;
    for(uint64_t i=0; i<branches.size(); i++) {
        const Branch& branch = branches[i];
        for(const MarkerInterval& markerInterval: branch.markerIntervals) {
            orientedReadBranches[markerInterval.orientedReadId].push_back(i);
        }
    }



    // Remove marker intervals for OrientedReadIds that appear
    // in more than one branch.
    for(const auto& p: orientedReadBranches) {
        const OrientedReadId orientedReadId = p.first;
        const vector<uint64_t>& orientedReadBranches = p.second;

        // If this OrientedReadId appears in only one branch, do nothing.
        if(orientedReadBranches.size() == 1) {
            continue;
        }

        // This OrientedReadId appears in more than one branch.
        // Remove it from all branches.
        // This could be written better but it is unlikely to become
        // a performance problem because this is an unusual case.
        for(const uint64_t i: orientedReadBranches) {
            Branch& branch = branches[i];
            for(MarkerInterval& markerInterval: branch.markerIntervals) {
                if(markerInterval.orientedReadId == orientedReadId) {
                    branch.markerIntervals.erase(branch.markerIntervals.begin() + i);
                    break;
                }
            }
        }
    }
}



Forks::Branch::Branch(
    EdgeId edgeId,
    const MarkerGraph& markerGraph) :
    edgeId(edgeId)
{
    // Get the marker intervals from the edge,
    // but don't store the ones for OrientedReadId's
    // that appear more than once.
    const span<const MarkerInterval> edgeMarkerIntervals = markerGraph.edgeMarkerIntervals[edgeId];

    for(uint64_t i=0; i<edgeMarkerIntervals.size(); i++) {
        const MarkerInterval& markerInterval = edgeMarkerIntervals[i];
        const OrientedReadId orientedReadId = markerInterval.orientedReadId;

        // Skip if same OrientedReadId as previous marker.
        if(i!=0 and edgeMarkerIntervals[i-1].orientedReadId == orientedReadId) {
            continue;
        }

        // Skip if same OrientedReadId as next marker.
        if(i!=edgeMarkerIntervals.size()-1 and edgeMarkerIntervals[i+1].orientedReadId == orientedReadId) {
            continue;
        }

        // Store this marker interval.
        markerIntervals.push_back(markerInterval);
    }
}



void Forks::constructEdgeInfos()
{
    for(const Fork& fork: forks) {
        for(uint64_t i=0; i<fork.branches.size(); i++) {
            const Branch& branch = fork.branches[i];
            edgeInfos.push_back(EdgeInfo(branch.edgeId, &fork, i));
        }
    }

    sort(edgeInfos.begin(), edgeInfos.end());
}



void Forks::constructOrientedReadEdgeInfos()
{
    const bool debug = true;

    const uint64_t orientedReadCount = markers.size();
    SHASTA_ASSERT((orientedReadCount % 2) == 0);

    orientedReadEdgeInfos.resize(orientedReadCount);
    for(const EdgeInfo& edgeInfo: edgeInfos) {
        const Branch& branch = edgeInfo.fork->branches[edgeInfo.branch];
        for(const MarkerInterval& markerInterval: branch.markerIntervals) {
            orientedReadEdgeInfos[markerInterval.orientedReadId.getValue()].push_back(
                OrientedReadEdgeInfo(&edgeInfo, markerInterval.ordinals));
        }
    }

    // Sort them.
    for(vector<OrientedReadEdgeInfo>& v: orientedReadEdgeInfos) {
        sort(v.begin(), v.end());
    }



    // Write them out.
    if(debug) {
        ofstream csv("OrientedReadEdgeInfos.csv");
        csv << "OrientedReadId,Ordinal0,Ordinal1,EdgeId,VertexId0,VertexId1,ForkVertexId,Direction,BranchIndex\n";
        const ReadId readCount = ReadId(orientedReadCount / 2);
        for(ReadId readId=0; readId<readCount; readId++) {
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);
                const vector<OrientedReadEdgeInfo>& infos = orientedReadEdgeInfos[orientedReadId.getValue()];
                for(const OrientedReadEdgeInfo& info: infos) {
                    const EdgeInfo& edgeInfo = *info.edgeInfo;
                    const MarkerGraph::Edge& edge = markerGraph.edges[edgeInfo.edgeId];
                    const Branch& branch = edgeInfo.fork->branches[edgeInfo.branch];
                    SHASTA_ASSERT(branch.edgeId == edgeInfo.edgeId);
                    csv << orientedReadId << ",";
                    csv << info.ordinals[0] << ",";
                    csv << info.ordinals[1] << ",";
                    csv << edgeInfo.edgeId << ",";
                    csv << edge.source << ",";
                    csv << edge.target << ",";
                    csv << (edgeInfo.fork->direction==ForkDirection::Forward ? edge.source : edge.target) << ",";
                    csv << directionString(edgeInfo.fork->direction) << ",";
                    csv << edgeInfo.branch << "\n";
                }
            }

        }
    }
}




void Forks::analyze(
    VertexId vertexId,
    ForkDirection direction,
    uint32_t maxDistance) const
{
    const bool debug = true;
    ofstream html;
    if(debug) {
        html.open("Fork.html");
        writeHtmlBegin(html,"Marker graph fork analysis");
        html << "<body>";
    }

    // Locate the Fork.
    const Fork* fork0Pointer = findFork(vertexId, direction);
    if(fork0Pointer == 0) {
        throw runtime_error("Fork not found.");
    }
    const Fork& fork0 = *fork0Pointer;

    // Debug output.
    if(debug) {
        fork0.writeHtml(html);
    }

    // Find nearby Forks.
    // The position of a Fork in  this vector will be used as the local Fork index.
    vector< pair<const Fork*, int32_t> > nearbyForks;
    findNearbyForks(fork0, maxDistance, nearbyForks);

    // Debug output.
    if(debug) {
        html << "<h2>Nearby forks</h2>"
            "<p>There are " << nearbyForks.size() << " forks within a maximum distance of " << maxDistance <<
            " markers from this fork."
            "<table><tr><th>Index<th>VertexId<th>Direction<th>Marker<br>offset";
        for(uint64_t id=0; id<nearbyForks.size(); id++) {
            const auto& p = nearbyForks[id];
            const Fork* fork1 = p.first;
            const int32_t offset = p.second;
            html << "<tr>"
                "<td class=centered>" << id <<
                "<td class=centered>" << fork1->vertexId <<
                "<td class=centered>" << directionString(fork1->direction) <<
                "<td class=centered>" << offset << endl;
        }
        html << "</table>";
    }



    // Gather the oriented reads in the nearby forks.
    std::set<OrientedReadId> nearbyForksOrientedReadIdsSet;
    for(const auto& p: nearbyForks) {
        const Fork* fork1 = p.first;
        for(const Branch& branch: fork1->branches) {
            for(const MarkerInterval& markerInterval: branch.markerIntervals) {
                nearbyForksOrientedReadIdsSet.insert(markerInterval.orientedReadId);
            }
        }
    }

    // Convert to a sorted vector.
    // The position of an oriented read in this vector will be used as the local
    // oriented read id index.
    vector<OrientedReadId> nearbyForksOrientedReadIds;
    copy(
        nearbyForksOrientedReadIdsSet.begin(),
        nearbyForksOrientedReadIdsSet.end(),
        back_inserter(nearbyForksOrientedReadIds));

    // Debug output.
    if(debug) {
        html << "<h2>Oriented reads in nearby forks</h2>" << endl;
        html << "<p>There are " << nearbyForksOrientedReadIds.size() <<
            " oriented reads in nearby forks." <<
            "<table><tr><th>Index<th>Oriented<br>Read Id";
        for(uint64_t id=0; id<nearbyForksOrientedReadIds.size(); id++) {
            const OrientedReadId orientedReadId = nearbyForksOrientedReadIds[id];
            html << "<tr><td class=centered>" << id <<
                "<td class=centered>" << orientedReadId;
        }
        html << "</table>";
    }



    // Create a table that tells us, for each oriented read, which branch
    // of which forks the oriented read appears in.
    // Indexed by [local oriented read id index][local Fork index]
    // (see above for definitionof those).
    // Contains the branch index (position in fork.branches or noBranch).
    const uint64_t noBranch = std::numeric_limits<uint64_t>::max();
    vector< vector<uint64_t> > branchTable(
        nearbyForksOrientedReadIds.size(),
        vector<uint64_t>(nearbyForks.size(), noBranch));
    for(uint64_t forkIndex=0; forkIndex<nearbyForks.size(); forkIndex++) {
        const Fork& fork = *nearbyForks[forkIndex].first;
        for(uint64_t branchIndex=0; branchIndex<fork.branches.size(); branchIndex++) {
            const Branch& branch = fork.branches[branchIndex];
            for(const MarkerInterval& markerInterval: branch.markerIntervals) {
                const OrientedReadId orientedReadId = markerInterval.orientedReadId;
                auto it = std::lower_bound(
                    nearbyForksOrientedReadIds.begin(),
                    nearbyForksOrientedReadIds.end(),
                    orientedReadId);
                SHASTA_ASSERT(it != nearbyForksOrientedReadIds.end());
                if(*it != orientedReadId) {
                    cout << "Assertion at " << forkIndex << " " << branchIndex << " " << orientedReadId << endl;
                }
                SHASTA_ASSERT(*it == orientedReadId);
                const uint64_t orientedReadIndex = it - nearbyForksOrientedReadIds.begin();
                SHASTA_ASSERT(branchTable[orientedReadIndex][forkIndex] == noBranch);
                branchTable[orientedReadIndex][forkIndex] = branchIndex;
            }
        }
    }



    // Debug output.
    if(debug) {
        const string symbol = "&#11044;";
        html << "<h2>Branch table</h2>"
            "<table>";

        // Header.
        html << "<tr><td><td>";
        for(uint64_t forkIndex=0; forkIndex<nearbyForks.size(); forkIndex++) {
            const Fork* fork1 = nearbyForks[forkIndex].first;
            html << "<td class=centered style='transform:rotate(270deg);";
            if(fork1 == &fork0) {
                html << "background-color:Beige'";
            }
            html << "'>" << forkIndex;
        }

        // Table body.
        for(uint64_t orientedReadIndex=0; orientedReadIndex<branchTable.size(); orientedReadIndex++) {
            const OrientedReadId orientedReadId = nearbyForksOrientedReadIds[orientedReadIndex];
            const auto& v = branchTable[orientedReadIndex];
            html << "<tr><td class=centered>" << orientedReadIndex <<
                "<td class=centered>" <<orientedReadId;
            for(uint64_t forkIndex=0; forkIndex<v.size(); forkIndex++) {
                const Fork* fork1 = nearbyForks[forkIndex].first;
                html << "<td class=centered title='" << orientedReadId << " " <<
                    fork1->vertexId << " " << directionString(fork1->direction) << "'";
                if(nearbyForks[forkIndex].first == &fork0) {
                    html << " style='background-color:Beige'";
                }
                html << ">";
                const uint64_t branchIndex = branchTable[orientedReadIndex][forkIndex];
                if(branchIndex != noBranch) {
                    const string color = "hsl(" + to_string(branchIndex*90) + ",100%, 50%)" ;
                    html << "<span style='color:" << color << "'>" << symbol;
                }

            }
        }
        html << "</table>";

    }



    // Debug output.
    if(debug) {
        html << "<body>";
        writeHtmlEnd(html);
    }
}



const Forks::Fork* Forks::findFork(VertexId vertexId, ForkDirection direction) const
{
    for(const Fork& fork: forks) {
        if(fork.vertexId == vertexId and fork.direction == direction) {
            return &fork;
        }
    }
    return 0;
}



void Forks::Fork::write(ostream& s) const
{
    s << ((direction == ForkDirection::Forward) ? "Forward" : "Backward") <<
        " fork at vertex " << vertexId << "\n";
    for(const Branch& branch: branches) {
        branch.write(s);
    }


    s << flush;
}


void Forks::Fork::writeHtml(ostream& html) const
{
    html <<
        "<h1>" << directionString(direction) <<
        " fork at vertex " << vertexId << "</h1>"
        "<p>This fork has " << branches.size() << " branches.";

    for(const Branch& branch: branches) {
        branch.writeHtml(html);
    }
}


void Forks::Branch::write(ostream& s) const
{
    s << "Branch at edge " << edgeId << "\n";
    for(const MarkerInterval& markerInterval: markerIntervals) {
        s << markerInterval.orientedReadId << " " <<
            markerInterval.ordinals[0] << " " <<
            markerInterval.ordinals[1] << "\n";
    }
}



void Forks::Branch::writeHtml(ostream& html) const
{
    html << "<h2>Branch at edge " << edgeId << "</h2>"
        "<p>This branch has " << markerIntervals.size() << " oriented reads."
        "<table><tr><th>Oriented<br>Read id<th>Ordinal0<th>Ordinal1";
    for(const MarkerInterval& markerInterval: markerIntervals) {
        html << "<tr><td class=centered>" << markerInterval.orientedReadId <<
            "<td class=centered>" << markerInterval.ordinals[0] <<
            "<td class=centered>" << markerInterval.ordinals[1];
    }
    html << "</table>";
}



// Find Forks that are close to a given Fork,
// and return them with their approximate marker offsets.
void Forks::findNearbyForks(
    const Fork& fork0,                              // Our starting Fork
    uint32_t maxDistance,                           // The maximum distance in markers
    vector< pair<const Fork*, int32_t> >& nearbyForks     // pairs(Fork, marker offset)
    ) const
{
    const int32_t doubleMaxDistance = 2 * maxDistance;

    std::map<const Fork*, vector<int32_t> > doubleOffsetMap;

    // Follow all the reads in this fork.
    for(const Branch& branch0: fork0.branches) {
        for(const MarkerInterval& markerInterval0: branch0.markerIntervals) {
            const OrientedReadId orientedReadId = markerInterval0.orientedReadId;
            // cout << "Following " << orientedReadId << endl;
            const int32_t doubleOrdinal0 = markerInterval0.ordinals[0] + markerInterval0.ordinals[1];

            // Loop over the OrientedReadEdgeInfos for this oriented read.
            const vector<OrientedReadEdgeInfo>& infos = orientedReadEdgeInfos[orientedReadId.getValue()];
            for(const OrientedReadEdgeInfo& info: infos) {
                const int32_t doubleOrdinal1 = info.ordinals[0] + info.ordinals[1];
                const int32_t doubleOffset = doubleOrdinal1 - doubleOrdinal0;
                // cout << info.ordinals[0] << " " << info.ordinals[1] << endl;

                // Check the distance.
                if(abs(doubleOffset) > doubleMaxDistance) {
                    // cout << "Too far." << endl;
                    continue;
                }

                const Fork* fork1 = info.edgeInfo->fork;
                doubleOffsetMap[fork1].push_back(doubleOffset);
            }
        }
    }



    // For each of the Forks we found, compute an average offset.
    nearbyForks.clear();
    for(const auto& p: doubleOffsetMap) {
        const Fork* fork1 = p.first;

        // Compute average offset.
        const auto& v = p.second;
        int64_t sum = 0;
        for(const auto& x: v) {
            sum += x;
        }
        const uint64_t n = v.size();
        const int32_t offset = int32_t(std::round(double(sum)/(double(2*n))));
        nearbyForks.push_back(make_pair(fork1, offset));
    }
    sort(nearbyForks.begin(), nearbyForks.end(),
        OrderPairsBySecondOnly<const Fork*, int32_t>());
}
