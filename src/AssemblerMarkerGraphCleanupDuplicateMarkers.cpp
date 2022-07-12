// Shasta.
#include "Assembler.hpp"
#include "MarkerConnectivityGraph.hpp"
#include "timestamp.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/connected_components.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "fstream.hpp"




/*******************************************************************************

During creation of marker graph vertices via the disjoint sets procedure,
it occasionally happens that a vertex has two or more markers from the same
oriented read. This is caused by small errors in marker alignments,
especially in the presence of multiple nearby copies of a marker.

These vertices, sometimes referred to in Shasta code and documentation as
"bad vertices", cause artifact in the marker graph which eventually lead
to assembly errors. As a result, "bad vertices" are forbidden by default -
they are simply not generated. Command line option --MarkerGraph.allowDuplicateMarkers
can be used to allow generation of these vertices.

However those missing vertices can cause other problems, especially
in problematic regions such as centromeres. It may be preferable to
"clean up" those bad vertices, rather than not generating them.
If the following command line options are used, "bad vertices" are allowed
to be generated during Assembler::createMarkerGraphVertices,
but then the code in this file (Assembler::cleanupDuplicateMarkers) is called to
"clean them up" in a couple of different ways:

--MarkerGraph.allowDuplicateMarkers --MarkerGraph.cleanupDuplicateMarkers

The clean up logic uses two common patterns that occur in "bad vertices".
The two patterns are defined in terms of the number of "duplicate markers" -
that is, markers for which another marker on the same oriented read
also exists in the same vertex. Note that the term "duplicate markers"
is misleading, as these markers are not really duplicate -
it's only the OrientedReadId that appears more than once in the same vertex.

- Pattern 1: If the number of duplicate markers is small (ratio of the
number of duplicate markers over total markers in the vertex is less
than pattern1Threshold), we simply remove the duplicate markers from the vertex.
If pattern1CreateNewVertices is true, each of the removed duplicate markers
is used to create a new vertex containing only that marker.
Otherwise, no vertex remains assigned to the duplicate markers.

- Pattern 2: If the vertex has too many duplicate markers for pattern 1,
we make an attempt to process it via pattern 2. We compute
connected components of the marker connectivity graph for the vertex,
but considering only duplicate markers.
It often happens that all of the large connected components
obtained in this way have no duplicate markers within each component.
If this happens, each connected component is used to generate a new vertex.
The non-duplicate markers each generate a one-marker vertex
if pattern2CreateNewVertices is true. Otherwise, they are
assigned no vertex.

In all cases, for both patterns, newly generated vertices must
have a number of markers at least equal to minCoverage,
and a number of markers on each strand at least equal to minCoveragePerStrand.
If this is not the case, no vertex is assigned to those markers.

*******************************************************************************/

void Assembler::cleanupDuplicateMarkers(
    uint64_t threadCount,
    uint64_t minCoverage,
    uint64_t minCoveragePerStrand,
    double pattern1Threshold,
    bool pattern1CreateNewVertices,
    bool pattern2CreateNewVertices)
{
    const bool debug = false;

    // Check that we have what we need.
    SHASTA_ASSERT(markers.isOpen());
    using CompressedVertexId = MarkerGraph::CompressedVertexId;
    MemoryMapped::Vector<CompressedVertexId>& vertexTable = markerGraph.vertexTable;
    SHASTA_ASSERT(vertexTable.isOpenWithWriteAccess);
    const MemoryMapped::VectorOfVectors<MarkerId, CompressedVertexId>& vertices = markerGraph.vertices();
    SHASTA_ASSERT(vertices.isOpen());
    const uint64_t vertexCount = vertices.size();
    SHASTA_ASSERT(markerGraph.reverseComplementVertex.isOpen);
    SHASTA_ASSERT(markerGraph.reverseComplementVertex.size() == vertexCount);

    cout << timestamp << "Cleaning up duplicate markers for " << vertexCount << " marker graph vertices." << endl;

    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store information that needs to be visible to the threads.
    cleanupDuplicateMarkersData.minCoverage = minCoverage;
    cleanupDuplicateMarkersData.minCoveragePerStrand = minCoveragePerStrand;
    cleanupDuplicateMarkersData.pattern1Threshold = pattern1Threshold;
    cleanupDuplicateMarkersData.pattern1CreateNewVertices = pattern1CreateNewVertices;
    cleanupDuplicateMarkersData.pattern2CreateNewVertices = pattern2CreateNewVertices;
    cleanupDuplicateMarkersData.badVertexCount = 0;
    cleanupDuplicateMarkersData.pattern1Count = 0;
    cleanupDuplicateMarkersData.pattern2Count = 0;
    cleanupDuplicateMarkersData.nextVertexId = vertexCount;

    // Process each vertex in multithreaded code.
    // For each vertex, we possibly change the vertexTable for the markers in that vertex (only).
    // Some vertices can disappear completely.
    const uint64_t batchSize = 100;
    setupLoadBalancing(vertexCount, batchSize);
    runThreads(&Assembler::cleanupDuplicateMarkersThreadFunction, threadCount);

    cout << "Found " << cleanupDuplicateMarkersData.badVertexCount <<
        " vertices with duplicate markers." << endl;
    cout << "Pattern 1 vertex count: " << cleanupDuplicateMarkersData.pattern1Count << endl;
    cout << "Pattern 2 vertex count: " << cleanupDuplicateMarkersData.pattern2Count << endl;

    // Renumber the vertex table to make sure vertices are numbered contiguously starting at 0.
    if(debug) {
        cout << "Maximum vertex id before renumbering of the vertex table " << cleanupDuplicateMarkersData.nextVertexId - 1 << endl;
    }
    const MarkerGraph::VertexId maxVertexId =
        markerGraph.renumberVertexTable(threadCount, cleanupDuplicateMarkersData.nextVertexId - 1);
    if(debug) {
        cout << "Maximum vertex id after renumbering of the vertex table " << maxVertexId << endl;
    }

    // Now we can recreate the vertices in the marker graph.
    markerGraph.createVerticesFromVertexTable(
        threadCount, maxVertexId);
    if(debug) {
        cout << "New number of vertices is " << markerGraph.vertices().size() << endl;
    }



    // Sanity check.
    if(debug) {
        for(MarkerGraph::VertexId vertexId=0; vertexId<markerGraph.vertices().size(); vertexId++) {
            if(markerGraph.vertices().size(vertexId) == 0) {
                cout << "Failing vertex id " << vertexId << endl;
            }
            SHASTA_ASSERT(markerGraph.vertices().size(vertexId) > 0);
        }
    }



    // Finally, recreate the reverse complement vertices.
    findMarkerGraphReverseComplementVertices(threadCount);


    cout << timestamp << "Cleaning up duplicate markers completed." << endl;
    cout << "Number of marker graph vertices is now " << markerGraph.vertices().size() << endl;
}



void Assembler::cleanupDuplicateMarkersThreadFunction(size_t threadId)
{
    const bool debug = true;
    ofstream out;
    if(debug) {
        out.open("cleanupDuplicateMarkers-" + to_string(threadId) + ".threadLog");
    }

    const uint64_t minCoverage = cleanupDuplicateMarkersData.minCoverage;
    const uint64_t minCoveragePerStrand = cleanupDuplicateMarkersData.minCoveragePerStrand;
    const double pattern1Threshold = cleanupDuplicateMarkersData.pattern1Threshold;
    const bool pattern1CreateNewVertices = cleanupDuplicateMarkersData.pattern1CreateNewVertices;
    const bool pattern2CreateNewVertices = cleanupDuplicateMarkersData.pattern2CreateNewVertices;

    uint64_t badVertexCount = 0;
    uint64_t pattern1Count = 0;
    uint64_t pattern2Count = 0;

    // The pairs (orientedReadId, marker ordinal) for the current vertex.
    vector<MarkerDescriptor> markerDescriptors;

    // Vector of flags that tells us which MarkerDescriptor are duplicate (duplicate in OrientedReadId only).
    vector<bool> isDuplicateOrientedReadId;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over all vertices in this batch.
        for(MarkerGraph::VertexId vertexId=begin; vertexId!=end; ++vertexId) {

            // Process one pair of reverse complemented vertices at a time.
            const MarkerGraph::VertexId vertexIdRc = markerGraph.reverseComplementVertex[vertexId];
            if(vertexIdRc < vertexId) {
                continue;
            }

            // If this vertex does not have duplicate markers, skip it.
            if(not isBadMarkerGraphVertex(vertexId)) {
                continue;
            }

            // This vertex has duplicate markers (more than one marker on the
            // same oriented read).
            if(vertexId == vertexIdRc) {
                ++badVertexCount;   // Unusual/exceptional case.
            } else {
                badVertexCount += 2;
            }
            if(debug) {
                if(vertexId == vertexIdRc) {
                    out << "Working on self-complementary vertex " <<
                        vertexId << endl;
                } else {
                    out << "Working on vertex " <<
                        vertexId << " and its reverse complement " << vertexIdRc << endl;
                }
            }

            // Get the pairs (orientedReadId, marker ordinal) for this vertex.
            const span<MarkerId> markerIds = markerGraph.getVertexMarkerIds(vertexId);
            const uint64_t markerCount = markerIds.size();
            SHASTA_ASSERT(markerCount > 1);
            markerDescriptors.clear();
            for(const MarkerId markerId: markerIds) {
                markerDescriptors.push_back(findMarkerId(markerId));
            }

            // Find the ones that are duplicate.
            // We take advantage of the fact that the markerDescriptors are sorted by OrientedReadId.
            isDuplicateOrientedReadId.resize(markerDescriptors.size());
            fill(isDuplicateOrientedReadId.begin(), isDuplicateOrientedReadId.end(), false);
            for(uint64_t i=1; i<markerCount; i++) {
                if(markerDescriptors[i-1].first == markerDescriptors[i].first) {
                    isDuplicateOrientedReadId[i-1] = true;
                    isDuplicateOrientedReadId[i] = true;
                }
            }
            const uint64_t duplicateCount =
                std::count(isDuplicateOrientedReadId.begin(), isDuplicateOrientedReadId.end(), true);
            if(debug) {
                out << duplicateCount << " duplicate markers out of " << markerCount << endl;
                for(uint64_t i=0; i<markerCount; i++) {
                    const auto& p = markerDescriptors[i];
                    out << p.first << " " << p.second;
                    if(isDuplicateOrientedReadId[i]) {
                        out << " duplicate";
                    }
                    out << endl;
                }
            }
            SHASTA_ASSERT(duplicateCount > 0);



            // Pattern 1: the number of duplicate markers is small.
            const double duplicateRatio = double(duplicateCount) / double(markerCount);
            if(duplicateRatio < pattern1Threshold) {
                if(debug) {
                    out << "Vertex " << vertexId << " processed as pattern 1 vertex." << endl;
                }
                SHASTA_ASSERT(duplicateCount < markerCount);
                if(vertexId == vertexIdRc) {
                    ++pattern1Count;   // Unusual/exceptional case.
                } else {
                    pattern1Count += 2;
                }
                cleanupDuplicateMarkersPattern1(vertexId,
                    minCoverage, minCoveragePerStrand,
                    pattern1CreateNewVertices, markerDescriptors, isDuplicateOrientedReadId,
                    debug, out);
                continue;
            }



            // Pattern 2: the duplicate vertices are in connected components,
            // and there are no duplications within each connected component.
            cleanupDuplicateMarkersPattern2(vertexId,
                minCoverage, minCoveragePerStrand,
                pattern2CreateNewVertices, markerDescriptors, isDuplicateOrientedReadId,
                debug, out);
            if(vertexId == vertexIdRc) {
                ++pattern2Count;   // Unusual/exceptional case.
            } else {
                pattern2Count += 2;
            }
        }
    }

    // Increment global counts.
    __sync_fetch_and_add(&cleanupDuplicateMarkersData.badVertexCount, badVertexCount);
    __sync_fetch_and_add(&cleanupDuplicateMarkersData.pattern1Count, pattern1Count);
    __sync_fetch_and_add(&cleanupDuplicateMarkersData.pattern2Count, pattern2Count);
}



void Assembler::cleanupDuplicateMarkersPattern1(
    MarkerGraph::VertexId vertexId,
    uint64_t minCoverage,
    uint64_t minCoveragePerStrand,
    bool createNewVertices,
    vector<MarkerDescriptor>& markerDescriptors,
    vector<bool>& isDuplicateOrientedReadId,
    bool debug,
    ostream& out)
{
    if(debug) {
        out << "Processing pattern 1 vertex " << vertexId << endl;
    }
    const uint64_t markerCount = markerDescriptors.size();
    SHASTA_ASSERT(isDuplicateOrientedReadId.size() == markerCount);
    array<uint64_t, 2> strandCoverage = {0, 0};

    // Loop over markers on this vertex.
    for(uint64_t i=0; i<markerCount; i++) {
        const pair<OrientedReadId, uint32_t>& p = markerDescriptors[i];

        // If not duplicate, just count coverage.
        if(not isDuplicateOrientedReadId[i]) {
            strandCoverage[p.first.getStrand()]++;
            continue;
        }

        // This is a duplicate marker.
        const MarkerId markerId = getMarkerId(p.first, p.second);
        const MarkerId markerIdRc = getReverseComplementMarkerId(p.first, p.second);

        if(createNewVertices and minCoverage<=1 and minCoveragePerStrand==0) {

            // Assign it to a new vertex.
            markerGraph.vertexTable[markerId] =
                cleanupDuplicateMarkersData.getAndIncrementNextVertexId();
            if(markerIdRc != markerId) {
                markerGraph.vertexTable[markerIdRc] =
                    cleanupDuplicateMarkersData.getAndIncrementNextVertexId();
            }
        } else {

            // Take this marker out of the current vertex, without
            // assigning it to a new vertex.
            markerGraph.vertexTable[markerId] = MarkerGraph::invalidCompressedVertexId;
            markerGraph.vertexTable[markerIdRc] = MarkerGraph::invalidCompressedVertexId;
        }
    }


    // Check if the remaining, non-duplicate vertices
    // satisfy our coverage criteria.
    if(
        strandCoverage[0]>=minCoveragePerStrand and
        strandCoverage[1]>=minCoveragePerStrand and
        (strandCoverage[0] + strandCoverage[1]) >= minCoverage
        ) {
        // We have enough coverage. We are done.
        return;
    }

    // If getting here, the vertex with the remaining non-duplicate markers
    // does not have enough coverage.
    // Assign all those markers to no vertex.
    for(uint64_t i=0; i<markerCount; i++) {
        if(isDuplicateOrientedReadId[i]) {
            continue;
        }

        const pair<OrientedReadId, uint32_t>& p = markerDescriptors[i];
        const MarkerId markerId = getMarkerId(p.first, p.second);
        const MarkerId markerIdRc = getReverseComplementMarkerId(p.first, p.second);
        markerGraph.vertexTable[markerId] = MarkerGraph::invalidCompressedVertexId;
        markerGraph.vertexTable[markerIdRc] = MarkerGraph::invalidCompressedVertexId;
    }

}



// Pattern 2: the duplicate vertices are in connected components,
// and there are no duplications within each connected component.
void Assembler::cleanupDuplicateMarkersPattern2(
    MarkerGraph::VertexId vertexId,
    uint64_t minCoverage,
    uint64_t minCoveragePerStrand,
    bool createNewVertices,
    vector<MarkerDescriptor>& markerDescriptors,
    vector<bool>& isDuplicateOrientedReadId,
    bool debug,
    ostream& out)
{
    using vertex_descriptor = MarkerConnectivityGraph::vertex_descriptor;
    using edge_descriptor = MarkerConnectivityGraph::edge_descriptor;

    if(debug) {
        out << "Processing pattern 2 vertex " << vertexId << endl;
    }
    const uint64_t markerCount = markerDescriptors.size();
    SHASTA_ASSERT(markerCount > 0);
    SHASTA_ASSERT(isDuplicateOrientedReadId.size() == markerCount);

    // Create the marker connectivity graph for this vertex.
    MarkerConnectivityGraph graph;
    MarkerConnectivityGraphVertexMap vertexMap;
    createMarkerConnectivityGraph(
        markerDescriptors.front().first, markerDescriptors.front().second, true, graph, vertexMap);
    SHASTA_ASSERT(num_vertices(graph) == markerCount);

    if(debug) {
        out << "The initial marker connectivity graph has " <<
            num_vertices(graph) << " vertices and " <<
            num_edges(graph) << " edges." << endl;
    }

    // Find the vertices that correspond to duplicate markers.
    std::set<vertex_descriptor> duplicateMarkerVertices;
    for(uint64_t i=0; i<markerCount; i++) {
        if(isDuplicateOrientedReadId[i]) {
            const MarkerDescriptor markerDescriptor = markerDescriptors[i];
            const auto it = vertexMap.find(markerDescriptor);
            SHASTA_ASSERT(it != vertexMap.end());
            const vertex_descriptor v = it->second;
            duplicateMarkerVertices.insert(v);
        }
    }

    // Only keep edges between duplicate markers.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, MarkerConnectivityGraph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        if(
            (duplicateMarkerVertices.find(v0) == duplicateMarkerVertices.end()) or
            (duplicateMarkerVertices.find(v1) == duplicateMarkerVertices.end())) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        remove_edge(e, graph);
    }

    if(debug) {
        out << "After edges involving non-duplicate vertices, the  marker connectivity graph has " <<
            num_vertices(graph) << " vertices and " <<
            num_edges(graph) << " edges." << endl;
    }

    // Compute connected components of the graph with the remaining edges.
    vector<uint64_t> component(markerCount);
    boost::connected_components(graph, &component[0]);
    std::map<uint64_t, vector<vertex_descriptor> > components;
    BGL_FORALL_VERTICES(v, graph, MarkerConnectivityGraph) {
        components[component[v]].push_back(v);
    }

    // Process the connected components one at a time.
    for(auto& p: components) {
        vector<vertex_descriptor>& component = p.second;

        // Get marker descriptors for this connected component and sort them.
        vector<MarkerDescriptor> componentDescriptors;
        for(const vertex_descriptor v: component) {
            componentDescriptors.push_back(graph[v]);
        }
        sort(componentDescriptors.begin(), componentDescriptors.end());

        if(debug) {
            out << "Found a connected component with " << componentDescriptors.size() << " markers:" << endl;
            for(const MarkerDescriptor& markerDescriptor: componentDescriptors) {
                out << markerDescriptor.first << " " << markerDescriptor.second << endl;
            }
        }

        // See if there are any duplicate markers (duplicate OrientedReadId's)
        // within this component.
        bool duplicatesFoundInThisComponent = false;
        for(uint64_t i=1; i<componentDescriptors.size(); i++) {
            const OrientedReadId thisOrientedReadId = componentDescriptors[i].first;
            const OrientedReadId previousOrientedReadId = componentDescriptors[i-1].first;
            if(thisOrientedReadId == previousOrientedReadId) {
                duplicatesFoundInThisComponent = true;
                break;
            }
        }
        if(debug) {
            out << "Duplicates found in this component: " << int(duplicatesFoundInThisComponent) << endl;
        }

        // Compute coverage for each strand.
        array<uint64_t, 2> strandCoverage = {0, 0};
        for(const MarkerDescriptor& markerDescriptor: componentDescriptors) {
            const OrientedReadId orientedReadId = markerDescriptor.first;
            strandCoverage[orientedReadId.getStrand()]++;
        }
        if(debug) {
            out << "Strand coverage for this component:" <<
                strandCoverage[0] << " " << strandCoverage[1] << endl;
        }



        // If there are no duplicates and coverage is sufficient,
        // create a new vertex with this component.
        if((not duplicatesFoundInThisComponent) and
            strandCoverage[0] >= minCoveragePerStrand and
            strandCoverage[1] >= minCoveragePerStrand and
            (strandCoverage[0] + strandCoverage[1]) >= minCoverage) {

            // Create a new vertex for this component, plus a second one for the reverse complement.
            const MarkerGraph::VertexId vertexId = cleanupDuplicateMarkersData.getAndIncrementNextVertexId();
            const MarkerGraph::VertexId vertexIdRc = cleanupDuplicateMarkersData.getAndIncrementNextVertexId();

            for(const MarkerDescriptor& markerDescriptor: componentDescriptors) {
                const MarkerId markerId = getMarkerId(markerDescriptor);
                const MarkerId markerIdRc = getReverseComplementMarkerId(markerDescriptor);
                markerGraph.vertexTable[markerId] = vertexId;
                if(markerIdRc != markerId) {
                    markerGraph.vertexTable[markerIdRc] = vertexIdRc;
                }
            }
            continue;
        }


        // We could not make a new vertex out of this component.
        // Each of this markers either becomes a new vertex on its own,
        // or gets assigned to no vertex.
        if(createNewVertices and minCoverage<=1 and minCoveragePerStrand==0) {
            for(const MarkerDescriptor& markerDescriptor: componentDescriptors) {
                const MarkerId markerId = getMarkerId(markerDescriptor);
                const MarkerId markerIdRc = getReverseComplementMarkerId(markerDescriptor);
                markerGraph.vertexTable[markerId] =
                    cleanupDuplicateMarkersData.getAndIncrementNextVertexId();
                if(markerIdRc != markerId) {
                    markerGraph.vertexTable[markerIdRc] =
                        cleanupDuplicateMarkersData.getAndIncrementNextVertexId();
                }
            }
        } else {
            for(const MarkerDescriptor& markerDescriptor: componentDescriptors) {
                const MarkerId markerId = getMarkerId(markerDescriptor);
                const MarkerId markerIdRc = getReverseComplementMarkerId(markerDescriptor);
                markerGraph.vertexTable[markerId] = MarkerGraph::invalidCompressedVertexId;
                markerGraph.vertexTable[markerIdRc] = MarkerGraph::invalidCompressedVertexId;
            }
        }
    }
}


