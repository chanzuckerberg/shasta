// Shasta.
#include "ReadGraph.hpp"
#include "deduplicate.hpp"
#include "orderPairs.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/pending/disjoint_sets.hpp>

// Standard library.
#include "fstream.hpp"
#include <map>
#include <queue>
#include <random>

const uint32_t ReadGraph::infiniteDistance = std::numeric_limits<uint32_t>::max();



void ReadGraph::unreserve() {
    if (edges.isOpenWithWriteAccess) edges.unreserve();
    if (connectivity.isOpenWithWriteAccess()) connectivity.unreserve();
}

// Compute a shortest path, disregarding edges flagged as cross-strand edges.
void ReadGraph::computeShortPath(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1,
    size_t maxDistance,

    // Edge ids of the shortest path starting at orientedReadId0 and
    // ending at orientedReadId1.
    vector<uint32_t>& path,

    // Work areas.
    vector<uint32_t>& distance, // One per vertex, equals infiniteDistance before and after.
    vector<OrientedReadId>& reachedVertices,   // For which distance is not infiniteDistance.
    vector<uint32_t>& parentEdges  // One per vertex
    )
{
    const bool debug = false;
    if(debug) {
        cout << "Looking for a shortest path between " << orientedReadId0 <<
            " " << orientedReadId1 << endl;
    }

    path.clear();

    // Initialize the BFS.
    std::queue<OrientedReadId> queuedVertices;
    queuedVertices.push(orientedReadId0);
    distance[orientedReadId0.getValue()] = 0;
    reachedVertices.clear();
    reachedVertices.push_back(orientedReadId0);


    // Do the BFS.
    while(!queuedVertices.empty()) {

        // Dequeue a vertex.
        const OrientedReadId vertex0 = queuedVertices.front();
        queuedVertices.pop();
        const uint32_t distance0 = distance[vertex0.getValue()];
        const uint32_t distance1 = distance0 + 1;
        if(debug) {
            cout << "Dequeued " << vertex0 << " at distance " << distance0 << endl;
        }

        // Loop over adjacent vertices.
        bool pathFound = false;
        for(const uint32_t edgeId: connectivity[vertex0.getValue()]) {
            const ReadGraphEdge& edge = edges[edgeId];
            if(edge.crossesStrands) {
                continue;
            }
            const OrientedReadId vertex1 = edge.getOther(vertex0);

            // If we did not encounter this vertex before, process it.
            if(distance[vertex1.getValue()] == infiniteDistance) {
                distance[vertex1.getValue()] = distance1;
                reachedVertices.push_back(vertex1);
                parentEdges[vertex1.getValue()] = edgeId;
                if(distance1 < maxDistance) {
                    if(debug) {
                        cout << "Enqueued " << vertex1 << endl;
                    }
                    queuedVertices.push(vertex1);
                }
            }

            // If we reached orientedReadId1, construct the path.
            if(vertex1 == orientedReadId1) {
                pathFound = true;
                // We have already cleared the path above.
                OrientedReadId vertex = vertex1;
                while(vertex != orientedReadId0) {
                    const uint32_t edgeId = parentEdges[vertex.getValue()];
                    path.push_back(edgeId);
                    vertex = edges[edgeId].getOther(vertex);
                }
                std::reverse(path.begin(), path.end());
                break;
            }

        }
        if(pathFound) {
            break;
        }

    }



    // Clean up.
    for(const OrientedReadId orientedReadId: reachedVertices) {
        distance[orientedReadId.getValue()] = infiniteDistance;
    }
    reachedVertices.clear();
}



// Find neighbors to distance 1.
// The neighbors are returned sorted and do not include orientedReadId0.
void ReadGraph::findNeighbors(
    OrientedReadId orientedReadId0,
    vector<OrientedReadId>& neighbors) const
{
    neighbors.clear();
    for(const uint32_t edgeId: connectivity[orientedReadId0.getValue()]) {
        const ReadGraphEdge& edge = edges[edgeId];
        const OrientedReadId orientedReadId1 = edge.getOther(orientedReadId0);
        neighbors.push_back(orientedReadId1);
    }
    sort(neighbors.begin(), neighbors.end());
}



// Find neighbors to specified maximum distance.
// The neighbors are returned sorted and do not include orientedReadId0.
void ReadGraph::findNeighbors(
    OrientedReadId orientedReadId,
    uint64_t maxDistance,
    vector<OrientedReadId>& neighbors) const
{
    // Initialize the BFS.
    std::queue<OrientedReadId> q;
    q.push(orientedReadId);
    std::map<OrientedReadId, uint64_t> distanceMap;
    distanceMap.insert(make_pair(orientedReadId, 0));



    // Do the BFS to the specified maximum distance.
    neighbors.clear();
    while(not q.empty()) {

        // Dequeue a vertex.
        const OrientedReadId orientedReadId0 = q.front();
        q.pop();
        const auto it0 = distanceMap.find(orientedReadId0);
        SHASTA_ASSERT(it0 != distanceMap.end());
        const uint64_t distance0 = it0->second;
        const uint64_t distance1 = distance0 + 1;
        SHASTA_ASSERT(distance1 <= maxDistance);

        // Loop over its neighbors.
        const span<const uint32_t> adjacentEdges = connectivity[orientedReadId0.getValue()];
        for(const uint32_t edgeId: adjacentEdges) {
            const ReadGraphEdge& edge = edges[edgeId];
            const OrientedReadId orientedReadId1 = edge.getOther(orientedReadId0);
            if(distanceMap.find(orientedReadId1) != distanceMap.end()) {
                // We already found orientedReadId1.
                continue;
            }
            neighbors.push_back(orientedReadId1);
            distanceMap.insert(make_pair(orientedReadId1, distance1));
            if(distance1 < maxDistance) {
                q.push(orientedReadId1);
            }
        }
    }



    // Sort the neighbors.
    sort(neighbors.begin(), neighbors.end());
}



// Find "bridges" from the read graph.
// Takes as input a vector<bool> that says, for each alignmentId,
// whether that alignment is used in the read graph.
// Updates that vector to set to false the entries corresponding
// to read graph "bridges".
void ReadGraph::findBridges(vector<bool>& keepAlignment)
{
    const uint64_t maxDistance = 2;

    // Vector to flag edges that will be removed.
    vector<bool> isToBeRemoved(edges.size(), false);

    // Work vectors to contain neighbors and
    // the corresponding edge ids.
    vector<OrientedReadId> neighbors;
    vector<uint32_t> neighborEdges;

    vector<uint64_t> rank;
    vector<uint64_t> parent;

    // Loop over reads. We only consider vertices corresponding
    // reads on strand 0, then for each edge to be removed
    // also flag its reverse complement.
    const uint32_t orientedReadCount = uint32_t(connectivity.size());
    SHASTA_ASSERT((orientedReadCount % 2) == 0);
    const uint32_t readCount = orientedReadCount / 2;
    for(ReadId readId0=0; readId0<readCount; readId0++) {
        const OrientedReadId orientedReadId0(readId0, 0);
        // cout << "Working on " << orientedReadId0 << endl;

        // Find neighbors within the specified distance.
        findNeighbors(orientedReadId0, maxDistance, neighbors);
        const uint64_t n = neighbors.size();
        if(n == 0) {
            continue;
        }

        /*
        cout << "Neighbors:" << endl;
        for(uint64_t i1=0; i1<n; i1++) {
            cout << i1 << " " << neighbors[i1] << endl;
        }
        */

        // Initialize the disjoint set data structure.
        rank.resize(n);
        parent.resize(n);
        boost::disjoint_sets<uint64_t*, uint64_t*> disjointSets(&rank[0], &parent[0]);
        for(uint64_t i=0; i<n; i++) {
            disjointSets.make_set(i);
        }

        // Compute connected components of the subgraph formed by the neighbors.
        for(uint64_t i1=0; i1<n; i1++) {
            const OrientedReadId orientedReadId1 = neighbors[i1];
            for(const uint32_t edgeId12: connectivity[orientedReadId1.getValue()]) {
                const ReadGraphEdge& edge12 = edges[edgeId12];
                const OrientedReadId orientedReadId2 = edge12.getOther(orientedReadId1);
                SHASTA_ASSERT(orientedReadId1 != orientedReadId2);

                // Only consider each edge once.
                if(orientedReadId2 < orientedReadId1) {
                    continue;
                }

                // Look up the other oriented read in our vector of neighbors.
                const auto it2 = lower_bound(neighbors.begin(), neighbors.end(), orientedReadId2);
                if((it2 == neighbors.end()) or (*it2 != orientedReadId2)) {
                    // Not a neighbor. Ignore.
                    continue;
                }
                SHASTA_ASSERT(*it2 == orientedReadId2);
                const uint64_t i2 = it2 - neighbors.begin();

                // Update our disjoint set data structure.
                // cout << "Updating: " << i1 << " " << i2 << endl;
                disjointSets.union_set(i1, i2);
            }
        }

        /*
        for(uint64_t i1=0; i1<n; i1++) {
            cout << i1 << " " << neighbors[i1] << " " << disjointSets.find_set(i1) << endl;
        }
        */

        // Gather connected components.
        vector< vector<uint64_t> > components(n);
        for(uint64_t i1=0; i1<n; i1++) {
            const uint64_t componentId = disjointSets.find_set(i1);
            components[componentId].push_back(i1);
        }

        // Sort them by size.
        vector<pair<uint64_t, uint64_t> > componentTable; // pair(componentId, size).
        for(uint64_t componentId=0; componentId<n;componentId++) {
            const auto& component = components[componentId];
            const uint64_t componentSize = component.size();
            if(componentSize > 0) {
                componentTable.push_back(make_pair(componentId, componentSize));
            }
        }
        sort(componentTable.begin(), componentTable.end(),
            OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());
        const uint64_t largestComponentId = componentTable.front().first;
        const vector<uint64_t>& largestComponent = components[largestComponentId];
        SHASTA_ASSERT(largestComponent.size() == componentTable.front().second);

        /*
        cout << "Largest component:" << endl;
        for(const uint64_t i1: largestComponent) {
            cout << i1 << " " << neighbors[i1] << endl;
        }
        */

        vector<bool> keep(n, false);
        for(const uint64_t i1: largestComponent) {
            keep[i1] = true;
        }

        /*
        for (uint64_t i1=0; i1<n; i1++) {
            if(not keep[i1]) {
                cout << "Remove " << orientedReadId0 << " " << neighbors[i1] << endl;
            }
        }
        */

        for (uint64_t i1=0; i1<n; i1++) {
            if(not keep[i1]) {
                const uint32_t edgeId = connectivity[orientedReadId0.getValue()][i1];
                const ReadGraphEdge& edge = edges[edgeId];
                keepAlignment[edge.alignmentId] = false;
            }
        }
    }

}



// Quick and dirty label propagation, without attempting to keep the clustering
// invariant under reverse complementing.
void ReadGraph::clustering() const
{
    // Initialize each vertex to its own cluster.
    const ReadId vertexCount = ReadId(connectivity.size());
    vector<ReadId> cluster(vertexCount);
    for(ReadId i=0; i<vertexCount; i++) {
        cluster[i] = i;
    }

    // Random number generator.
    const uint32_t seed = 231;
    std::mt19937 randomSource(seed);
    std::uniform_int_distribution<ReadId> uniformDistribution(0, vertexCount-1);

    // Iterate.
    const uint64_t sweepCount = 100;
    const uint64_t iterationCount = sweepCount * vertexCount;
    vector<ReadId> neighborLabels;
    vector<ReadId> labelFrequencies;
    for(uint64_t iteration=0; iteration<iterationCount; iteration++) {

        // Randomly pick a vertex to update.
        const ReadId vertexId0 = uniformDistribution(randomSource);
        const OrientedReadId orientedReadId0 = OrientedReadId(vertexId0);

        // Get the clusters of its neighbors.
        neighborLabels.clear();
        for(const uint32_t edgeId: connectivity[vertexId0]) {
            const ReadGraphEdge& edge = edges[edgeId];
            const OrientedReadId orientedReadId1 = edge.getOther(orientedReadId0);
            const ReadId vertexId1 = orientedReadId1.getValue();
            const ReadId label1 = cluster[vertexId1];
            neighborLabels.push_back(label1);
        }

        if(neighborLabels.empty()) {
            continue;
        }

        // Count the occurrences of each label.
        deduplicateAndCount(neighborLabels, labelFrequencies);

        // Find the most frequent label.
        ReadId bestLabel = neighborLabels.front();
        ReadId bestLabelFrequency = labelFrequencies.front();
        for(uint64_t i=1; i<neighborLabels.size(); i++) {
            if(labelFrequencies[i] > bestLabelFrequency) {
                bestLabel = neighborLabels[i];
                bestLabelFrequency = labelFrequencies[i];
            }
        }

        // Assign to this vertex the most frequent cluster in its neighbors.
        cluster[vertexId0] = bestLabel;

    }



    // Summarize the clusters.
    std::map<ReadId, ReadId> clusterMap;    // (cluster, number of vertices)
    for(ReadId vertexId=0; vertexId<vertexCount; vertexId++) {
        const ReadId vertexCluster = cluster[vertexId];
        const auto it = clusterMap.find(vertexCluster);
        if(it == clusterMap.end()) {
            clusterMap.insert(make_pair(vertexCluster, 1));
        } else {
            ++(it->second);
        }
    }

    {
        ofstream csv("Clusters.csv");
        csv << "Cluster,Size\n";
        for(const auto& p: clusterMap) {
            csv << p.first << ",";
            csv << p.second << "\n";
        }
    }



    // Summarize edges by cluster.
    // Key: pair(cluster0, cluster1) with cluster0<=cluster1;
    // Value: number of edges.
    std::map< pair<ReadId, ReadId>, ReadId> edgeTable;
    for(const ReadGraphEdge& edge: edges) {
        const ReadId vertexId0 = edge.orientedReadIds[0].getValue();
        const ReadId vertexId1 = edge.orientedReadIds[1].getValue();
        ReadId cluster0 = cluster[vertexId0];
        ReadId cluster1 = cluster[vertexId1];
        if(cluster1 < cluster0) {
            swap(cluster0, cluster1);
        }
        const auto it = edgeTable.find(make_pair(cluster0, cluster1));
        if(it == edgeTable.end()) {
            edgeTable.insert(make_pair(make_pair(cluster0, cluster1), 1));
        } else {
            ++(it->second);
        }
    }

    {
        ofstream csv("EdgeByClusters.csv");
        csv << "Cluster0,Cluster1,Size0,Size1,Edges\n";
        for(const auto& p: edgeTable) {
            const ReadId cluster0 = p.first.first;
            const ReadId cluster1 = p.first.second;
            csv << cluster0 << ",";
            csv << cluster1 << ",";
            csv << clusterMap[cluster0] << ",";
            csv << clusterMap[cluster1] << ",";
            csv << p.second << "\n";
        }
    }


    ofstream graphOut("ReadGraph.dot");
    graphOut << "graph ReadGraph {\n"
        "tooltip=\" \"";
    for(ReadId vertexId=0; vertexId<vertexCount; vertexId++) {
        const OrientedReadId orientedReadId = OrientedReadId(vertexId);
        const ReadId vertexLabel = cluster[vertexId];
        graphOut << "\"" << orientedReadId << "\"[" <<
            " tooltip=\"" << orientedReadId << " " << vertexLabel << "\""
            " color = \"" << 0.05 * (vertexLabel%20) << " 1 1\"];\n";
    }
    for(const ReadGraphEdge& edge: edges) {
        graphOut << "\"" << edge.orientedReadIds[0] << "\"--\"" <<
            edge.orientedReadIds[1] << "\"";
        /*
        const ReadId vertexId0 = edge.orientedReadIds[0].getValue();
        const ReadId vertexId1 = edge.orientedReadIds[1].getValue();
        const ReadId cluster0 = cluster[vertexId0];
        const ReadId cluster1 = cluster[vertexId1];
        if(cluster0 != cluster1) {
            graphOut << " [color=red]";
        }
        */
        graphOut << ";\n";
    }
    graphOut << "}\n";
}
