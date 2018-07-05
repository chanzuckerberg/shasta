// Nanopore2.
#include "Assembler.hpp"
#include "LocalReadGraph.hpp"
#include "orderPairs.hpp"
#include "timestamp.hpp"
using namespace ChanZuckerberg;
using namespace Nanopore2;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard libraries.
#include <queue>
#include "tuple.hpp"



// Compute connected components of the global read graph.
// We only keep components that are large enough.
void Assembler::computeReadGraphComponents(
    size_t minComponentSize,        // Minimum size for a connected component to be kept.
    size_t minAlignedMarkerCount,
    size_t maxTrim
    )
{
    checkReadsAreOpen();
    checkKmersAreOpen();
    checkMarkersAreOpen();
    checkOverlapsAreOpen();
    checkAlignmentDataAreOpen();
    const OrientedReadId::Int orientedReadCount = OrientedReadId::Int(2*reads.size());

    // Initialize the disjoint set data structures.
    vector<ReadId> rank(orientedReadCount);
    vector<ReadId> parent(orientedReadCount);
    boost::disjoint_sets<ReadId*, ReadId*> disjointSets(&rank[0], &parent[0]);
    for(OrientedReadId::Int orientedReadId=0; orientedReadId<orientedReadCount; orientedReadId++) {
        disjointSets.make_set(orientedReadId);
    }



    // Loop over alignments which satisfy our criteria.
    for(size_t i=0; i<alignmentData.size(); i++) {
        if((i%1000000) == 0) {
            cout << timestamp << i << " of " << alignmentData.size() << endl;
        }
        const AlignmentData& ad = alignmentData[i];

        // If the number of markers in the alignment is too small, skip.
        if(ad.info.markerCount < minAlignedMarkerCount) {
            continue;
        }

        // Compute the left and right trim (bases excluded from
        // the alignment).
        OrientedReadId orientedReadId0(ad.readIds[0], 0);
        OrientedReadId orientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
        uint32_t leftTrim;
        uint32_t rightTrim;
        tie(leftTrim, rightTrim) = computeTrim(
            orientedReadId0,
            orientedReadId1,
            ad.info);

        // If the trim is too high, skip.
        if(leftTrim>maxTrim || rightTrim>maxTrim) {
            continue;
        }

        // Join the connected components that these two
        // oriented reads belong to.
        disjointSets.union_set(
            orientedReadId0.getValue(),
            orientedReadId1.getValue());
        orientedReadId0.flipStrand();
        orientedReadId1.flipStrand();
        disjointSets.union_set(
            orientedReadId0.getValue(),
            orientedReadId1.getValue());
    }



    // Store the component that each oriented read belongs to.
    overlapGraphComponent.createNew(largeDataName("OverlapGraphComponent"), largeDataPageSize);
    overlapGraphComponent.resize(orientedReadCount);
    for(OrientedReadId::Int orientedReadId=0; orientedReadId<orientedReadCount; orientedReadId++) {
        overlapGraphComponent[orientedReadId] = disjointSets.find_set(orientedReadId);
    }

    // Gather the vertices (oriented read ids) of each component.
    // The vertices of each component are sorted by construction.
    std::map<OrientedReadId::Int, vector<OrientedReadId::Int> > componentMap;
    for(OrientedReadId::Int orientedReadId=0; orientedReadId<orientedReadCount; orientedReadId++) {
        componentMap[overlapGraphComponent[orientedReadId]].push_back(orientedReadId);
    }



    // Sort the components by decreasing size.
    // We only keep components that are large enough.
    vector< pair<OrientedReadId::Int, OrientedReadId::Int> > componentTable;
    for(const auto& p: componentMap) {
        const auto componentId = p.first;
        const auto& component = p.second;
        const auto componentSize = component.size();
        CZI_ASSERT(componentSize > 0);
        if(componentSize >= minComponentSize) {
            // Keep this component.
            componentTable.push_back(make_pair(componentId, componentSize));

        }
    }
    sort(componentTable.begin(), componentTable.end(),
        OrderPairsBySecondOnlyGreater<ReadId, ReadId>());



    // Renumber the components.
    std::map<ReadId, ReadId> componentNumberMap;
    for(ReadId newComponentId=0; newComponentId<componentTable.size(); newComponentId++) {
        const auto& p = componentTable[newComponentId];
        const auto oldComponentId = p.first;
        componentNumberMap.insert(make_pair(oldComponentId, newComponentId));
    }
    for(ReadId i=0; i<orientedReadCount; i++) {
        ReadId& readComponentId = overlapGraphComponent[i];
        const auto it = componentNumberMap.find(readComponentId);
        if(it == componentNumberMap.end()) {
            readComponentId = std::numeric_limits<ReadId>::max();
        } else {
            readComponentId = it->second;
        }
    }

    // Permanently store the renumbered components.
    overlapGraphComponents.createNew(largeDataName("OverlapGraphComponents"), largeDataPageSize);
    for(ReadId newComponentId=0; newComponentId<componentTable.size(); newComponentId++) {
        const auto& p = componentTable[newComponentId];
        const auto oldComponentId = p.first;
        const auto& component = componentMap[oldComponentId];
        CZI_ASSERT(component.size() == p.second);
        overlapGraphComponents.appendVector();
        for(const auto& orientedRead: component) {
            overlapGraphComponents.append(OrientedReadId(orientedRead));
        }
    }



    // Write out the size of each component we kept.
    for(ReadId newComponentId=0; newComponentId<overlapGraphComponents.size(); newComponentId++) {
        cout << "Component " << newComponentId << " has size ";
        cout << overlapGraphComponents.size(newComponentId);
        cout << "." << endl;

    }


}



// Create a local read graph starting from a given oriented read
// and walking out a given distance on the global read graph.
void Assembler::createLocalReadGraph(
    ReadId readId,
    Strand strand,
    size_t minFrequency,            // Minimum number of minHash hits to generate an edge.
    size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
    size_t maxTrim,                 // Maximum left/right trim to generate an edge.
    uint32_t distance                 // How far to go from starting oriented read.
)
{
    createLocalReadGraph(OrientedReadId(readId, strand),
        minFrequency, minAlignedMarkerCount, maxTrim,
        distance);
}
void Assembler::createLocalReadGraph(
    OrientedReadId orientedReadId,
    size_t minFrequency,            // Minimum number of minHash hits to generate an edge.
    size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
    size_t maxTrim,                 // Maximum left/right trim to generate an edge.
    uint32_t maxDistance              // How far to go from starting oriented read.
)
{
    // Check that we have what we need.
    checkReadsAreOpen();
    checkReadNamesAreOpen();
    checkKmersAreOpen();
    checkMarkersAreOpen();
    checkAlignmentDataAreOpen();

    // Create the LocalReadGraph.
    LocalReadGraph graph;
    createLocalReadGraph(orientedReadId,
        minFrequency, minAlignedMarkerCount, maxTrim, maxDistance,
        graph);

    cout << "The local read graph has " << num_vertices(graph);
    cout << " vertices and " << num_edges(graph) << " edges." << endl;
    graph.write("LocalReadGraph.dot", maxDistance);
    writeLocalReadGraphToFasta(graph, "LocalReadGraph.fasta");

}



// Create a local read graph starting from a given oriented read
// and walking out a given distance on the global read graph.
void Assembler::createLocalReadGraph(
    OrientedReadId orientedReadIdStart,
    size_t minFrequency,            // Minimum number of minHash hits to generate an edge.
    size_t minAlignedMarkerCount,   // Minimum number of alignment markers to generate an edge.
    size_t maxTrim,                 // Maximum left/right trim to generate an edge.
    uint32_t maxDistance,           // How far to go from starting oriented read.
    LocalReadGraph& graph)
{
    // Add the starting vertex.
    graph.addVertex(orientedReadIdStart,
        uint32_t(reads[orientedReadIdStart.getReadId()].baseCount), 0);

    // Initialize a BFS starting at the start vertex.
    std::queue<OrientedReadId> q;
    q.push(orientedReadIdStart);



    // Do the BFS.
    while(!q.empty()) {

        // Dequeue a vertex.
        const OrientedReadId orientedReadId0 = q.front();
        // cout << "Dequeued " << orientedReadId0;
        // cout << " with " << overlapTable.size(orientedReadId0.getValue()) << " overlaps." << endl;
        q.pop();
        const uint32_t distance0 = graph.getDistance(orientedReadId0);
        const uint32_t distance1 = distance0 + 1;

        // Loop over overlaps/alignments involving this vertex.
        for(const uint64_t i: alignmentTable[orientedReadId0.getValue()]) {
            CZI_ASSERT(i < alignmentData.size());
            const AlignmentData& ad = alignmentData[i];

            // If the alignment involves too few markers, skip.
            if(ad.info.markerCount < minAlignedMarkerCount) {
                continue;
            }

            // To compute the trim, keep into account the fact
            // that the stored AlignmentInfo was computed for
            // the ReadId's stored in the Overlap, with the first one on strand 0.
            const OrientedReadId overlapOrientedReadId0(ad.readIds[0], 0);
            const OrientedReadId overlapOrientedReadId1(ad.readIds[1], ad.isSameStrand ? 0 : 1);
            uint32_t leftTrim;
            uint32_t rightTrim;
            tie(leftTrim, rightTrim) = computeTrim(
                overlapOrientedReadId0,
                overlapOrientedReadId1,
                ad.info);
            if(leftTrim>maxTrim || rightTrim>maxTrim) {
                continue;
            }

            // The overlap and the alignment satisfy our criteria.
            // Get the other oriented read involved in this overlap.
            const OrientedReadId orientedReadId1 = ad.getOther(orientedReadId0);

            // Add the vertex for orientedReadId1, if necessary.
            // Also add orientedReadId1 to the queue, unless
            // we already reached the maximum distance.
            if(!graph.vertexExists(orientedReadId1)) {
                graph.addVertex(orientedReadId1,
                    uint32_t(reads[orientedReadId1.getReadId()].baseCount), distance1);
                if(distance1 < maxDistance) {
                    q.push(orientedReadId1);
                    // cout << "Enqueued " << orientedReadId1 << endl;
                }
            }

            // Add the edge.
            graph.addEdge(orientedReadId0, orientedReadId1,
                ad.info);
        }

    }

}



// Write in fasta format the sequences of the vertices of a local read graph.
void Assembler::writeLocalReadGraphToFasta(
    const LocalReadGraph& graph,
    const string& fileName)
{
    ofstream fasta(fileName);
    BGL_FORALL_VERTICES(v, graph, LocalReadGraph) {
        const OrientedReadId orientedReadId(graph[v].orientedReadId);
        writeOrientedRead(orientedReadId, fasta);
    }
}
