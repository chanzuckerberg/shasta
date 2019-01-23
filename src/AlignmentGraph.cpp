// shasta.
#include "AlignmentGraph.hpp"
#include "Alignment.hpp"
using namespace ChanZuckerberg;
using namespace shasta;

// Boost libraries.
// The boost gil library includes png.h,
// then uses int_p_NULL which is not defined in
// all versions of png (see Boost bug 3908,
// flaged as fixed but it is not obvious that that
// is the case). To deal with this, we defensively
// include pngh., then define int_p_NULL if necessary.
#include <png.h>
#ifndef int_p_NULL
#define int_p_NULL (int *)NULL
#endif
#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png_dynamic_io.hpp>

// Standard library.
#include "algorithm.hpp"
#include "fstream.hpp"


// Ccompute an alignment of the markers of two oriented reads.
void ChanZuckerberg::shasta::align(

    // Markers of the two oriented reads to be aligned, sorted by KmerId.
    const array<vector<MarkerWithOrdinal>, 2>& markers,

    // The maximum ordinal skip to be tolerated between successive markers
    // in the alignment.
    size_t maxSkip,

    // Marker frequency threshold.
    // When computing an alignment between two oriented reads,
    // marker kmers that appear more than this number of times
    // in either of the two oriented reads are discarded
    // (in both oriented reads).
    // Change to size_t when conversion completed.
    uint32_t maxMarkerFrequency,

    // Flag to control various types of debug output.
    bool debug,

    // The AlignmentGraph can be reused.
    // For performance, it should be reused when doing many alignments.
    AlignmentGraph& graph,

    // The computed alignment.
    // This should also be reused when performance is important.
    Alignment& alignment,

    // Also create alignment summary information.
    AlignmentInfo& alignmentInfo
    )
{
    graph.create(markers, maxMarkerFrequency, maxSkip, debug,
        alignment, alignmentInfo);
}




void AlignmentGraph::create(
    const array<vector<MarkerWithOrdinal>, 2>& markers,
    uint32_t maxMarkerFrequency,
    size_t maxSkip,
    bool debug,
    Alignment& alignment,
    AlignmentInfo& alignmentInfo)
{

    // Start with an empty graph.
    clear();

    // Write out the markers.
    if(debug) {
        writeMarkers(markers[0], "Markers-ByKmerId-0.csv");
        writeMarkers(markers[1], "Markers-ByKmerId-1.csv");
    }

    // Create the vertices - one for each pair of common markers.
    createVertices(markers, maxMarkerFrequency);
    sortVertices();

    // Add the start and finish vertices.
    // This must be done after sorting the remaining vertices,
    // because sorting alters vertex descriptors.
    vStart = addVertex();
    vFinish = addVertex();
    doneAddingVertices();
    if(debug) {
        writeVertices("AlignmentGraphVertices.csv");
    }

    // Create the edges.
    createEdges(uint32_t(markers[0].size()), uint32_t(markers[1].size()), maxSkip);
    if(debug) {
        writeEdges("AlignmentGraphEdges.csv");
    }
    doneAddingEdges();

    if(debug) {
        cout << "The alignment graph has " << vertexCount()-2;
        cout << " marker vertices and " << edgeCount() << " edges." << endl;
        writeGraphviz("AlignmentGraph.dot");
    }

    // Look for the shortest path between vStart and vFinish.
    findShortestPath(*this, vStart, vFinish, shortestPath, queue);
    if(shortestPath.empty()) {
        alignment.ordinals.clear();
        if(debug) {
            cout << "The shortest path is empty." << endl;
        }
        return;
    }
    if(debug) {
        cout << "The shortest path has " << shortestPath.size()-2;
        cout << " k-mer vertices." << endl;
        writeShortestPath("ShortestPath.csv");
    }

    // Store the alignment.
    alignment.ordinals.clear();
    for(const vertex_descriptor v: shortestPath) {
        if(v==vStart || v==vFinish) {
            continue;
        }
        const auto& vertex = (*this)[v];
        alignment.ordinals.push_back(
            array<uint32_t, 2>({uint32_t(vertex.ordinals[0]), uint32_t(vertex.ordinals[1])}));
    }

    // Store the alignment info.
    alignmentInfo.create(alignment, uint32_t(markers[0].size()), uint32_t(markers[1].size()));

    if(debug) {
        writeImage(markers[0], markers[1], alignment, "Alignment.png");
    }
}


void AlignmentGraph::writeMarkers(
    const vector<MarkerWithOrdinal>& markers,
    const string& fileName
    )
{
    ofstream csv(fileName);
    csv << "Index,KmerId,Ordinal,Position\n";

    for(size_t i=0; i<markers.size(); i++) {
        const MarkerWithOrdinal& marker = markers[i];
        csv << i << "," << marker.kmerId << "," << marker.ordinal << "," << marker.position << "\n";
    }

}



void AlignmentGraph::createVertices(
    const array<vector<MarkerWithOrdinal>, 2>& markers,
    uint32_t maxMarkerFrequency)
{
    // Some shorthands for readability.
    const vector<MarkerWithOrdinal>& markers0 = markers[0];
    const vector<MarkerWithOrdinal>& markers1 = markers[1];

    // Some iterators we will need.
    using MarkerIterator = vector<MarkerWithOrdinal>::const_iterator;
    const MarkerIterator begin0 = markers0.begin();
    const MarkerIterator end0   = markers0.end();
    const MarkerIterator begin1 = markers1.begin();
    const MarkerIterator end1   = markers1.end();

    // Initialize isLowFrequencyMarker flags to all true.
    // We will set to false the ones that need it,
    // when we encounter long streaks of the same marker.
    for(size_t i=0; i<2; i++) {
        isLowFrequencyMarker[i].clear();
        isLowFrequencyMarker[i].resize(markers[i].size(), true);
    }

    // Joint loop over the markers, looking for common k-mer ids.
    auto it0 = begin0;
    auto it1 = begin1;
    while(it0!=end0 && it1!=end1) {
        if(it0->kmerId < it1->kmerId) {
            ++it0;
        } else if(it1->kmerId < it0->kmerId) {
            ++it1;
        } else {

            // We found a common k-mer id.
            const KmerId kmerId = it0->kmerId;


            // This k-mer could appear more than once in each of the oriented reads,
            // so we need to find the streak of this k-mer in kmers0 and kmers1.
            MarkerIterator it0Begin = it0;
            MarkerIterator it1Begin = it1;
            MarkerIterator it0End = it0Begin;
            MarkerIterator it1End = it1Begin;
            while(it0End!=end0 && it0End->kmerId==kmerId) {
                ++it0End;
            }
            while(it1End!=end1 && it1End->kmerId==kmerId) {
                ++it1End;
            }
            const size_t streakLength0 = it0End - it0Begin;
            const size_t streakLength1 = it1End - it1Begin;


            if(streakLength0>maxMarkerFrequency || streakLength1>maxMarkerFrequency) {

                // At least one of these streaks is too long.
                // Flag these markers as high frequency markers.
                for(MarkerIterator jt0=it0Begin; jt0!=it0End; ++jt0) {
                    isLowFrequencyMarker[0][jt0->ordinal]= false;
                }
                for(MarkerIterator jt1=it1Begin; jt1!=it1End; ++jt1) {
                    isLowFrequencyMarker[1][jt1->ordinal]= false;
                }

            } else {

                // Both streaks are short enough.
                // Generate vertices in the alignment graph.

                // Loop over pairs in the streaks.
                for(MarkerIterator jt0=it0Begin; jt0!=it0End; ++jt0) {
                    for(MarkerIterator jt1=it1Begin; jt1!=it1End; ++jt1) {

                        // Generate a vertex corresponding to this pair
                        // of occurrences of this common k-mer.
                        AlignmentGraphVertex vertex;
                        vertex.kmerId = kmerId;
                        vertex.indexes[0] = jt0 - begin0;
                        vertex.indexes[1] = jt1 - begin1;
                        vertex.positions[0] = jt0->position;
                        vertex.positions[1] = jt1->position;
                        vertex.ordinals[0] = jt0->ordinal;
                        vertex.ordinals[1] = jt1->ordinal;
                        addVertex(vertex);
                    }
                }

            }

            // Continue joint loop over k-mers.
            it0 = it0End;
            it1 = it1End;
        }
    }


    // Compute correctedOrdinals, the ordinals keeping into account
    // only low frequency markers.
    for(size_t i=0; i<2; i++) {
        correctedOrdinals[i].resize(markers[i].size());
        uint32_t correctedOrdinal = 0;
        for(size_t j=0; j<markers[i].size(); j++) {
            if(isLowFrequencyMarker[i][j]) {
                correctedOrdinals[i][j] = correctedOrdinal++;
            } else {
                correctedOrdinals[i][j] =  std::numeric_limits<uint32_t>::max();
            }
        }
    }
}



void AlignmentGraph::writeVertices(const string& fileName) const
{
    ofstream csv(fileName);
    csv << "Vertex,KmerId,Index0,Index1,Ordinal0,Ordinal1,Position0,Position1\n";

    BGL_FORALL_VERTICES(v, *this, AlignmentGraph) {
        if(v==vStart || v==vFinish) {
            continue;
        }
        const AlignmentGraphVertex& vertex = (*this)[v];
        csv << v.v << ",";
        csv << vertex.kmerId << ",";
        csv << vertex.indexes[0] << ",";
        csv << vertex.indexes[1] << ",";
        csv << vertex.ordinals[0] << ",";
        csv << vertex.ordinals[1] << ",";
        csv << vertex.positions[0] << ",";
        csv << vertex.positions[1] << "\n";

    }

}



void AlignmentGraph::createEdges(
    uint32_t markerCount0,
    uint32_t markerCount1,
    size_t maxSkip)
{
    AlignmentGraph& graph = *this;

    const vertex_iterator itBegin = verticesBegin();
    const vertex_iterator itEnd = verticesEnd();

    // Loop over pairs of vertices with a gap no more
    // than maxSkip in sequence 0.
    // The vertices are sorted by position in sequence 0.
    for(vertex_iterator itA=itBegin; itA!=itEnd; ++itA) {
        const vertex_descriptor vA = *itA;
        if(vA==vStart || vA==vFinish) {
            continue;
        }
        const auto& vertexA = graph[vA];
        const int ordinalA0 = int(vertexA.ordinals[0]);
        const int ordinalA1 = int(vertexA.ordinals[1]);
        const int correctedOrdinalA0 = int(correctedOrdinals[0][ordinalA0]);
        const int correctedOrdinalA1 = int(correctedOrdinals[1][ordinalA1]);
        CZI_ASSERT(correctedOrdinalA0 < int(markerCount0));
        CZI_ASSERT(correctedOrdinalA1 < int(markerCount1));

        vertex_iterator itB = itA;
        ++itB;
        for(; itB!=itEnd; ++itB) {
            const vertex_descriptor vB = *itB;
            if(vB==vStart || vB==vFinish) {
                continue;
            }
            const auto& vertexB = graph[vB];
            const int ordinalB0 = int(vertexB.ordinals[0]);
            CZI_ASSERT(ordinalB0 >= ordinalA0);
            const int correctedOrdinalB0 = int(correctedOrdinals[0][ordinalB0]);
            CZI_ASSERT(correctedOrdinalB0 < int(markerCount0));

            // If we got too far, we can end the inner loop,
            // because vertices are sorted by position in sequence 0.
            if(correctedOrdinalB0 > correctedOrdinalA0 + int(maxSkip)) {
                break;
            }
            const int ordinalB1 = int(vertexB.ordinals[1]);
            const int correctedOrdinalB1 = int(correctedOrdinals[1][ordinalB1]);
            CZI_ASSERT(correctedOrdinalB1 < int(markerCount1));

            // Check that the skip in the 1 direction is less than maxSkip.
            if(abs(correctedOrdinalB1 - correctedOrdinalA1) > maxSkip) {
                continue;
            }

            // If getting here, we will add an edge.
            // We need to compute the weight.
            const int delta0 = correctedOrdinalB0 - correctedOrdinalA0;
            const int delta1 = correctedOrdinalB1 - correctedOrdinalA1;
            // const size_t weight = delta0*delta0 + delta1*delta1;
            const size_t weight = abs(delta0-1) + abs(delta1-1);

            // Add the edge.
            addEdge(vA, vB, AlignmentGraphEdge(weight));
        }
    }


    // Create edges from vStart and vFinish to all other vertices.
    for(vertex_iterator it=itBegin; it!=itEnd; ++it) {
        const vertex_descriptor v = *it;
        if(v==vStart || v==vFinish) {
            continue;
        }
        const auto& vertex = graph[v];
        const int ordinal0 = int(vertex.ordinals[0]);
        const int ordinal1 = int(vertex.ordinals[1]);
        const int correctedOrdinal0 = int(correctedOrdinals[0][ordinal0]);
        const int correctedOrdinal1 = int(correctedOrdinals[1][ordinal1]);
        const int deltaFinish0 = int(markerCount0) - correctedOrdinal0;
        const int deltaFinish1 = int(markerCount1) - correctedOrdinal1;

        /*
        addEdge(v, vStart,  AlignmentGraphEdge(
            ordinal0 * ordinal0 +
            ordinal1 * ordinal1));
        addEdge(v, vFinish, AlignmentGraphEdge(
            deltaFinish0 * deltaFinish0 +
            deltaFinish1 * deltaFinish1));
        */
        addEdge(v, vStart,  AlignmentGraphEdge(
            abs(correctedOrdinal0) +
            abs(correctedOrdinal1)));
        addEdge(v, vFinish, AlignmentGraphEdge(
            abs(deltaFinish0) +
            abs(deltaFinish1)));
    }


}



void AlignmentGraph::writeEdges(const string& fileName) const
{
    ofstream csv(fileName);
    csv << "V0,V1,Weight\n";

    BGL_FORALL_EDGES(e, *this, AlignmentGraph) {
        const vertex_descriptor v0 = source(e);
        const vertex_descriptor v1 = target(e);
        const AlignmentGraphEdge& edge = (*this)[e];
        csv << v0.v << ",";
        csv << v1.v << ",";
        csv << edge.weight << "\n";
    }

}



// Write in graphviz format, without the start and finish vertices.
void AlignmentGraph::writeGraphviz(const string& fileName) const
{
    // Open the output file and write the beginning of the graph.
    ofstream file(fileName);
    file << "graph G {\n";

    // Write the vertices.
    BGL_FORALL_VERTICES(v, *this, AlignmentGraph) {
        if(v==vStart || v==vFinish) {
            continue;
        }
        const AlignmentGraphVertex& vertex = (*this)[v];
        file << v.v;

        file << " [label=\"";
        file << "Vertex " << v.v << "\\n";
        file << "Kmer id " << vertex.kmerId << "\\n";
        file << "Ordinals " << vertex.ordinals[0] << " " << vertex.ordinals[1] << "\\n";
        file << "\"]";

        file << ";\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, *this, AlignmentGraph) {
        const vertex_descriptor v0 = source(e);
        if(v0==vStart || v0==vFinish) {
            continue;
        }
        const vertex_descriptor v1 = target(e);
        if(v1==vStart || v1==vFinish) {
            continue;
        }
        const AlignmentGraphEdge& edge = (*this)[e];
        file << v0.v << "--";
        file << v1.v;

        file << " [label=\"";
        file << edge.weight;
        file << "\"]";

        file << ";\n";
    }

    // Write the end of the graph to the output file.
    file << "}\n";
}


void AlignmentGraph::writeShortestPath(const string& fileName) const
{
    ofstream csv(fileName);
    csv << "Vertex,KmerId,Index0,Index1,Ordinal0,Ordinal1,Position0,Position1\n";

    for(const vertex_descriptor v: shortestPath) {
        if(v==vStart || v==vFinish) {
            continue;
        }

        const AlignmentGraphVertex& vertex = (*this)[v];
        csv << v.v << ",";
        csv << vertex.kmerId << ",";
        csv << vertex.indexes[0] << ",";
        csv << vertex.indexes[1] << ",";
        csv << vertex.ordinals[0] << ",";
        csv << vertex.ordinals[1] << ",";
        csv << vertex.positions[0] << ",";
        csv << vertex.positions[1] << "\n";
    }

}



// Write an image representing the markers and the computed alignment
// in 2-D ordinal space.
void AlignmentGraph::writeImage(
    const vector<MarkerWithOrdinal>& markers0,
    const vector<MarkerWithOrdinal>& markers1,
    const Alignment& alignment,
    const string& fileName) const
{
    using namespace boost::gil;

    // Create the image and the view.
    const size_t n0 = markers0.size();
    const size_t n1 = markers1.size();
    rgb8_image_t image(n0, n1);
    rgb8_image_t::view_t imageView = view(image);

    // Initialize it to black.
    const rgb8_pixel_t black(0, 0, 0);
    for(size_t i0=0; i0<n0; i0++) {
        for(size_t i1=0; i1<n1; i1++) {
            imageView(i0, i1) = black;
        }
    }


    // Write a grid.
    const size_t smallGridSpacing = 10;
    const rgb8_pixel_t lightGrey(12, 12, 12);
    for(size_t i0=0; i0<n0; i0+=smallGridSpacing) {
        for(size_t i1=0; i1<n1; i1++) {
            imageView(i0, i1) = lightGrey;
        }
    }
    for(size_t i1=0; i1<n1; i1+=smallGridSpacing) {
        for(size_t i0=0; i0<n0; i0++) {
            imageView(i0, i1) = lightGrey;
        }
    }
    const size_t largeGridSpacing = 50;
    const rgb8_pixel_t grey(32, 32, 32);
    for(size_t i0=0; i0<n0; i0+=largeGridSpacing) {
        for(size_t i1=0; i1<n1; i1++) {
            imageView(i0, i1) = grey;
        }
    }
    for(size_t i1=0; i1<n1; i1+=largeGridSpacing) {
        for(size_t i0=0; i0<n0; i0++) {
            imageView(i0, i1) = grey;
        }
    }

    // Write the markers.
    const rgb8_pixel_t red(255, 0, 0);
    for(size_t i0=0; i0<n0; i0++) {
        const MarkerWithOrdinal& marker0 = markers0[i0];
        for(size_t i1=0; i1<n1; i1++) {
            const MarkerWithOrdinal& marker1 = markers1[i1];
            if(marker0.kmerId == marker1.kmerId) {
                imageView(marker0.ordinal, marker1.ordinal) = red;
            }
        }
    }

    // Write the alignment.
    const rgb8_pixel_t green(0, 255, 0);
    for(const auto& p: alignment.ordinals) {
        imageView(p[0], p[1]) = green;
    }

    // Write it out.
    png_write_view(fileName, imageView);
}
