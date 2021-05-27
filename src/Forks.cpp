#include "Forks.hpp"
using namespace shasta;


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
        csv << "OrientedReadId,Ordinal0,Ordinal1,EdgeId,VertexId0,VertexId1,Direction\n";
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
                    csv << (edgeInfo.fork->direction == ForkDirection::Forward ? "Forward" : "Backward") << "\n";
                }
            }

        }
    }
}

