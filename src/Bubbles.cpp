#include "Bubbles.hpp"
#include "Assembler.hpp"
#include "computeLayout.hpp"
#include "orderPairs.hpp"
#include "prefixLength.hpp"
#include "writeGraph.hpp"
using namespace shasta;

#include <boost/graph/connected_components.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <queue>


Bubbles::Bubbles(
    const Assembler& assembler) :
    assembler(assembler)
{
    // Bubbles with discordant ratio greater than this value
    // are considered bad and will be removed.
    const double discordantRatioThreshold = 0.1;

    // Number of phasing iterations.
    const uint64_t iterationCount = 1;

    // Seed for random number generator.
    uint64_t randomSeed = 231;

    // Only used for graphics.
    const double minRelativePhase = 0.;

    findBubbles();
    cout << "Found " << bubbles.size() << " bubbles." << endl;
    cout << "Found "<< 2*assembler.getReads().readCount() << " oriented reads." << endl;

    fillOrientedReadsTable();
    writeOrientedReadsTable();

    createBubbleGraph();
    flagBadBubbles();
    removeBadBubbles(discordantRatioThreshold);
    writeBubbles();
    writeBubbleGraphGraphviz();
    cout << "The bubble graph has " << num_vertices(bubbleGraph) <<
        " vertices and " << num_edges(bubbleGraph) << " edges." << endl;

    phase(iterationCount, randomSeed, minRelativePhase);

}





// For now we only consider diploid bubbles, defined using the
// following strict criteria:
// - Source vertex v0 has out-degree 2.
// - Target vertex v1 has in-degree 2.
// - There are two parallel edges eA and eB, both v0->v1.

void Bubbles::findBubbles()
{
    const bool debug = false;

    // Loop over possible choices for v0 (source vertex of the bubble).
    const AssemblyGraph& assemblyGraph = *assembler.assemblyGraphPointer;
    const AssemblyGraph::VertexId vertexCount = assemblyGraph.vertices.size();
    uint64_t suppressedDueToRepeats = 0;
    for(AssemblyGraph::VertexId v0=0; v0<vertexCount; v0++) {

        // Check the out-degree.
        const auto outEdges0 = assemblyGraph.edgesBySource[v0];
        if(outEdges0.size() != 2) {
            continue;
        }

        const AssemblyGraph::EdgeId eA = outEdges0[0];
        const AssemblyGraph::EdgeId eB = outEdges0[1];
        const AssemblyGraph::Edge& edgeA = assemblyGraph.edges[eA];
        const AssemblyGraph::Edge& edgeB = assemblyGraph.edges[eB];

        // Check the target vertex.
        const AssemblyGraph::VertexId v1 = edgeA.target;
        if(edgeB.target != v1) {
            continue;
        }
        const auto inEdges1 = assemblyGraph.edgesByTarget[v1];
        if(inEdges1.size() != 2) {
            continue;
        }

        // Assemble the RLE sequence for the two edges.
        vector<Base> sequenceA;
        vector<Base> sequenceB;
        assembler.assembleAssemblyGraphEdgeRleStrict(eA, sequenceA);
        assembler.assembleAssemblyGraphEdgeRleStrict(eB, sequenceB);


        // If the RLE sequences are identical, skip it.
        // This can happen due to missing alignments.
        if(sequenceA == sequenceB) {
            continue;
        }

        if(debug) {
            const MarkerGraph::VertexId u0 = assemblyGraph.vertices[v0];
            const MarkerGraph::VertexId u1 = assemblyGraph.vertices[v1];
            cout << "Bubble at marker graph " << u0 << "->" << u1 << "\n";
            copy(sequenceA.begin(), sequenceA.end(), ostream_iterator<Base>(cout));
            cout << "\n";
            copy(sequenceB.begin(), sequenceB.end(), ostream_iterator<Base>(cout));
            cout << "\n";
        }

        // If the two sequences differ by a 2- of 3-base repeat,
        // this is probably an erroneous bubble.
        if(isShortRepeatCopyNumberDifference(sequenceA, sequenceB)) {
            ++suppressedDueToRepeats;
            if(debug) {
                cout << "Is short repeat difference - skipping.\n";
            }
            continue;
        }

        // We have a diploid bubble.
        bubbles.push_back(Bubble(v0, v1, eA, eB, assembler.markerGraph, assemblyGraph));

    }

    cout << suppressedDueToRepeats << " bubbles were suppressed due to 2- or 3-base repeats." << endl;
}



void Bubbles::writeBubbles()
{
    ofstream csv("Bubbles.csv");
    csv << "BubbleId,VertexId0,VertexId1,CoverageA,CoverageB,ConcordantSum,DiscordantSum,DiscordantRatio\n";

    for(uint64_t bubbleId=0; bubbleId<bubbles.size(); bubbleId++) {
        const Bubble& bubble = bubbles[bubbleId];
        csv << bubbleId << ",";
        csv << bubble.mv0 << ",";
        csv << bubble.mv1 << ",";
        csv << bubble.orientedReadIds[0].size() << ",";
        csv << bubble.orientedReadIds[1].size() << ",";
        csv << bubble.concordantSum << ",";
        csv << bubble.discordantSum << ",";
        csv << bubble.discordantRatio() << "\n";
    }
}


// Figure out if two sequences differ only by copy numbers in
// a 2- or 3-base repeat.
bool Bubbles::isShortRepeatCopyNumberDifference(
    const vector<Base>& x,
    const vector<Base>& y)
{
    const bool debug = false;

    // Get the lengths and their differences.
    const uint64_t nx = x.size();
    const uint64_t ny = y.size();

    // When they have the same length, we can certainly return false.
    if(nx == ny) {
        return false;
    }

    // Recursive call so x is shorter than y.
    if(ny < nx) {
        return isShortRepeatCopyNumberDifference(y, x);
    }
    SHASTA_ASSERT(nx < ny);

    // If the length difference is not a multiple of 2 or 3,
    // we can certainly return false.
    const uint64_t dn = ny - nx;
    if( ((dn%2) != 0) and
        ((dn%3) != 0) ) {
        if(debug) {
            cout << "Length difference is not a multiple of 2 or 3\n";
            cout << "nx " << nx << ", ny " << ny << ", dn " << dn << "\n";
        }
        return false;
    }

    const uint64_t commonPrefixLength = shasta::commonPrefixLength(x, y);
    const uint64_t commonSuffixLength = shasta::commonSuffixLength(x, y);

    // Find the portion of y that is not in x.
    uint64_t ix = commonPrefixLength;
    uint64_t iy = commonPrefixLength;
    uint64_t jx = nx - commonSuffixLength;
    uint64_t jy = ny - commonSuffixLength;
    while((jx<ix) or (jy<iy)) {
        ++jx;
        ++jy;
    }

    if(debug) {
        cout << "commonPrefixLength " << commonPrefixLength << "\n";
        cout << "commonSuffixLength " << commonSuffixLength << "\n";
        cout << "nx " << nx << "\n";
        cout << "ny " << ny << "\n";
        cout << "dn " << dn << "\n";
        cout << "ix " << ix << "\n";
        cout << "jx " << jx << "\n";
        cout << "iy " << iy << "\n";
        cout << "jy " << jy << "\n";
    }


    if(ix != jx) {
        // There is more than just an insertion.
        return false;
    }


    // x and y differ by an insertion in iy of range [iy, jy).
    SHASTA_ASSERT(ix == jx);
    SHASTA_ASSERT(jy - iy == dn);



    // Check for k base repeat.
    // We kept the entire common prefix, so we can check just to the left of the insertion.
    for(uint64_t k=2; k<=3; k++) {
        if((dn % k) != 0) {
            continue;
        }
        if(debug) {
            cout << "Trying k = " << k << "\n";
        }

        // Check that the inserted bases are a repeat with period k.
        const uint64_t m = dn / k;
        bool repeatViolationFound = false;
        for(uint64_t i=0; i<m; i++) {
            for(uint64_t j=0; j<k; j++) {
                if(y[iy + i*k + j] != y[iy + j]) {
                    repeatViolationFound = true;
                    break;
                }
            }
        }
        if(repeatViolationFound) {
            // The inserted portion is not a repeat of period k.
            continue;
        }

        // Check the previous k bases in both x and y.
        if(ix < k) {
            continue;
        }
        if(iy < k) {
            continue;
        }
        for(uint64_t j=0; j<k; j++) {
            if(y[iy - k + j] != y[ix + j]) {
                repeatViolationFound = true;
                break;
            }
            if(x[ix - k + j] != y[ix + j]) {
                repeatViolationFound = true;
                break;
            }
        }
        if(repeatViolationFound) {
            continue;
        }

        // It is an insertion in y of a repeat of length k.
        return true;
    }



    // None of the k values we tried worked.
    return false;
}



Bubbles::Bubble::Bubble(
    AssemblyGraph::VertexId av0,
    AssemblyGraph::VertexId av1,
    AssemblyGraph::EdgeId eA,
    AssemblyGraph::EdgeId eB,
    const MarkerGraph& markerGraph,
    const AssemblyGraph& assemblyGraph) :
    av0(av0),
    av1(av1),
    aEdgeIds({eA, eB}),
    mv0(assemblyGraph.vertices[av0]),
    mv1(assemblyGraph.vertices[av1])
{
    fillInOrientedReadIds(markerGraph, assemblyGraph);
}



void Bubbles::Bubble::fillInOrientedReadIds(
    const MarkerGraph& markerGraph,
    const AssemblyGraph& assemblyGraph)
{

    // Loop over both sides of the bubble.
    array<std::set<OrientedReadId>, 2> orientedReadIdSets;
    for(uint64_t i=0; i<2; i++) {
        const AssemblyGraph::EdgeId aEdgeId = aEdgeIds[i];

        // Loop over marker graph edges of this side.
        const span<const MarkerGraph::EdgeId> mEdgeIds = assemblyGraph.edgeLists[aEdgeId];
        for(const MarkerGraph::EdgeId mEdgeId : mEdgeIds) {

            // Loop over MarkerIntervals of this marker graph edge.
            const span<const MarkerInterval> markerIntervals =
                markerGraph.edgeMarkerIntervals[mEdgeId];
            for(const MarkerInterval& markerInterval: markerIntervals) {
                orientedReadIdSets[i].insert(markerInterval.orientedReadId);
            }
        }
    }

    // Now store them, discarding the ones that appear on both sides.
    for(uint64_t i=0; i<2; i++) {
        const std::set<OrientedReadId>& x = orientedReadIdSets[i];
        const std::set<OrientedReadId>& y = orientedReadIdSets[1-i];
        for(const OrientedReadId orientedReadId: x) {
            if(y.find(orientedReadId) == y.end()) {
                orientedReadIds[i].push_back(orientedReadId);
            }
        }
    }
}



void Bubbles::fillOrientedReadsTable()
{
    orientedReadsTable.resize(2* assembler.getReads().readCount());
    for(uint64_t bubbleId=0; bubbleId<bubbles.size(); bubbleId++) {
        const Bubble& bubble = bubbles[bubbleId];
        for(uint64_t side=0; side<2; side++) {
            for(const OrientedReadId orientedReadId: bubble.orientedReadIds[side]) {
                orientedReadsTable[orientedReadId.getValue()].push_back(make_pair(bubbleId, side));
            }
        }
    }
}



void Bubbles::writeOrientedReadsTable()
{
    ofstream csv("OrientedReadsTable.csv");
    csv << "OrientedReadId,Bubble,Side\n";
    for(ReadId readId=0; readId<assembler.getReads().readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            for(const auto& v: orientedReadsTable[orientedReadId.getValue()]) {
                csv << orientedReadId << ",";
                csv << v.first << ",";
                csv << v.second << "\n";
            }
        }
    }
}



void Bubbles::createBubbleGraph()
{
    // Create the vertices.
    for(uint64_t bubbleId=0; bubbleId<bubbles.size(); bubbleId++) {
        bubbleGraph.vertexTable.push_back(add_vertex(BubbleGraphVertex(bubbleId), bubbleGraph));
    }


    // Create the edges.
    for(uint64_t bubbleIdA=0; bubbleIdA<bubbles.size(); bubbleIdA++) {
        const Bubble& bubbleA = bubbles[bubbleIdA];
        const BubbleGraph::vertex_descriptor vA = bubbleGraph.vertexTable[bubbleIdA];

        for(uint64_t sideA=0; sideA<2; sideA++) {
            const vector<OrientedReadId>& orientedReadIds = bubbleA.orientedReadIds[sideA];

            for(const OrientedReadId orientedReadId: orientedReadIds) {
                const auto& v = orientedReadsTable[orientedReadId.getValue()];

                for(const auto& p: v) {
                    const uint64_t bubbleIdB = p.first;
                    if(bubbleIdB <= bubbleIdA) {
                        continue;
                    }
                    const uint64_t sideB = p.second;
                    const BubbleGraph::vertex_descriptor vB = bubbleGraph.vertexTable[bubbleIdB];

                    // Locate the edge between these two bubbles,
                    // creating it if necessary.
                    BubbleGraph::edge_descriptor e;
                    bool edgeWasFound = false;
                    tie(e, edgeWasFound) = edge(vA, vB, bubbleGraph);
                    if(not edgeWasFound) {
                        tie(e, edgeWasFound) = add_edge(vA, vB, bubbleGraph);
                    }
                    SHASTA_ASSERT(edgeWasFound);

                    // Update the matrix.
                    BubbleGraphEdge& edge = bubbleGraph[e];
                    ++edge.matrix[sideA][sideB];
                }
            }
        }
    }



    ofstream csv("BubbleGraphEdges.csv");
    csv << "BubbleIdA,BubbleIdB,m00,m11,m01,m10\n";
    BGL_FORALL_EDGES(e, bubbleGraph, BubbleGraph) {
        const BubbleGraphEdge& edge = bubbleGraph[e];
        csv << bubbleGraph[source(e, bubbleGraph)].bubbleId << ",";
        csv << bubbleGraph[target(e, bubbleGraph)].bubbleId << ",";
        csv << edge.matrix[0][0] << ",";
        csv << edge.matrix[1][1] << ",";
        csv << edge.matrix[0][1] << ",";
        csv << edge.matrix[1][0] << "\n";
    }
}



void Bubbles::writeBubbleGraphGraphviz() const
{

    ofstream out("BubbleGraph.dot");
    out << "graph G{\n"
        "node [shape=point];\n";

    BGL_FORALL_VERTICES(v, bubbleGraph, BubbleGraph) {
        const BubbleGraphVertex& vertex = bubbleGraph[v];
        const uint64_t bubbleId = vertex.bubbleId;
        const Bubble& bubble = bubbles[bubbleId];
        const double discordantRatio = bubble.discordantRatio();
        SHASTA_ASSERT(discordantRatio <= 0.5);
        const double hue = (1. - 2. * discordantRatio) * 0.333333;  //  Good = green, bad = red, via yellow
        out << bubbleId <<
            " [color=\"" << hue << ",1.,0.8\""
            " tooltip=\""  << bubbleId << " " << discordantRatio << "\""
            "];\n";
    }


    BGL_FORALL_EDGES(e, bubbleGraph, BubbleGraph) {
        const BubbleGraphEdge& edge = bubbleGraph[e];
        const uint64_t bubbleIdA = bubbleGraph[source(e, bubbleGraph)].bubbleId;
        const uint64_t bubbleIdB = bubbleGraph[target(e, bubbleGraph)].bubbleId;

        const double hue = (1. - edge.ambiguity()) * 0.333333; //  Good = green, bad = red, via yellow

        out << bubbleIdA << "--" << bubbleIdB <<
            " ["
            // "penwidth=" << 0.1*double(total) <<
            " color=\"" << hue << ",1.,0.8\""
            "];\n";
    }

    out << "}\n";
}



// Write a single component of the BubbleGraph in svg format.
// To compute sfdp layout, only consider edges
// for which relativePhase() >= minRelativePhase.
void Bubbles::writeBubbleGraphComponentSvg(
    uint64_t componentId,
    const vector<BubbleGraph::vertex_descriptor>& component,
    double minRelativePhase) const
{
    using vertex_descriptor = BubbleGraph::vertex_descriptor;
    using edge_descriptor = BubbleGraph::edge_descriptor;

    // Map the vertices of this component to integers.
    std::map<vertex_descriptor, uint64_t> componentVertexMap;
    for(uint64_t i=0; i<component.size(); i++) {
        componentVertexMap.insert(make_pair(component[i], i));
    }

    // Create a graph with only this connected component and only
    // edges with relativePhase >= minRelativePhase.
    // This will be used to compute the sfdp layout.
    ComponentGraph componentGraph(component.size());
    for(uint64_t i0=0; i0<component.size(); i0++) {
        const vertex_descriptor v0 = component[i0];
        BGL_FORALL_OUTEDGES(v0, e, bubbleGraph, BubbleGraph) {
            if(bubbleGraph[e].relativePhase() < minRelativePhase) {
                continue;
            }
            const vertex_descriptor v1 = target(e, bubbleGraph);
            auto it = componentVertexMap.find(v1);
            SHASTA_ASSERT(it != componentVertexMap.end());
            const uint64_t i1 = it->second;
            if(i0 < i1) {
                add_edge(i0, i1, componentGraph);
            }
        }
    }

    // Compute the layout of this graph.
    std::map<ComponentGraph::vertex_descriptor, array<double, 2> > positionMap;
    SHASTA_ASSERT(computeLayout(componentGraph, "sfdp", 600., positionMap) == ComputeLayoutReturnCode::Success);
    for(uint64_t i=0; i<component.size(); i++) {
        componentGraph[i].position = positionMap[i];
        // cout << componentGraph[i].position[0] << " " << componentGraph[i].position[1] << "\n";
    }

    // Add all the remaining edges.
    for(uint64_t i0=0; i0<component.size(); i0++) {
        const vertex_descriptor v0 = component[i0];
        BGL_FORALL_OUTEDGES(v0, e, bubbleGraph, BubbleGraph) {
            if(bubbleGraph[e].relativePhase() >= minRelativePhase) {
                // We already added this edge.
                continue;
            }
            const vertex_descriptor v1 = target(e, bubbleGraph);
            auto it = componentVertexMap.find(v1);
            SHASTA_ASSERT(it != componentVertexMap.end());
            const uint64_t i1 = it->second;
            if(i0 < i1) {
                add_edge(i0, i1, componentGraph);
            }
        }
    }



    // Graphics scaling.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = std::numeric_limits<double>::max();
    double yMax = std::numeric_limits<double>::min();
    BGL_FORALL_VERTICES_T(v, componentGraph, ComponentGraph) {
        const auto& position = componentGraph[v].position;
        xMin = min(xMin, position[0]);
        xMax = max(xMax, position[0]);
        yMin = min(yMin, position[1]);
        yMax = max(yMax, position[1]);
    }
    const double xyRange = max(xMax-xMin, yMax-yMin);
    const int svgSize = 1024;
    const double vertexRadiusPixels = 3.;
    const double vertexRadius = vertexRadiusPixels * xyRange / double(svgSize);
    const double edgeThicknessPixels = 0.4;
    const double edgeThickness = edgeThicknessPixels * xyRange / double(svgSize);



    // Write the componentGraph in svg format.

    // Vertex attributes. Color by phase.
    std::map<ComponentGraph::vertex_descriptor, WriteGraph::VertexAttributes> vertexAttributes;
    BGL_FORALL_VERTICES(v, componentGraph, ComponentGraph) {
        const int64_t phase = bubbleGraph[component[v]].phase;
        auto& attributes = vertexAttributes[v];
        attributes.radius = vertexRadius;
        if(phase == +1) {
            attributes.color = "hsl(240,50%,50%)";
        }
        if(phase == -1) {
            attributes.color = "hsl(300,50%,50%)";
        }
    }

    // Edge attributes. Color by relative phase.
    std::map<ComponentGraph::edge_descriptor, WriteGraph::EdgeAttributes> edgeAttributes;
    BGL_FORALL_EDGES(e, componentGraph, ComponentGraph) {
        const vertex_descriptor v0 = component[source(e, componentGraph)];
        const vertex_descriptor v1 = component[target(e, componentGraph)];
        edge_descriptor eb;
        bool edgeWasFound = false;
        tie(eb, edgeWasFound) = edge(v0, v1, bubbleGraph);
        const BubbleGraphEdge& edge = bubbleGraph[eb];
        const double relativePhase = edge.relativePhase();
        const double hue = (1. + relativePhase) * 60.; /// Goes from 0 (red) to 120 (green).
        auto& attributes = edgeAttributes[e];
        attributes.thickness = edgeThickness;
        attributes.color = "hsla(" + to_string(int(hue)) + ",50%,50%,20%)";
    }

    // Draw the svg.
    ofstream out("BubbleGraph-Component-" + to_string(componentId) + ".html");
    out << "<html><body>";
    WriteGraph::writeSvg(
        componentGraph,
        "component" + to_string(componentId),
        svgSize, svgSize,
        vertexAttributes,
        edgeAttributes,
        out);
    out << "</body></html>";
}



// Use the BubbleGraph to flag bad bubbles.
void Bubbles::flagBadBubbles()
{
    for(uint64_t bubbleId=0; bubbleId<bubbles.size(); bubbleId++) {
        Bubble& bubble = bubbles[bubbleId];
        bubble.concordantSum = 0;
        bubble.discordantSum = 0;
        BGL_FORALL_OUTEDGES(bubbleGraph.vertexDescriptor(bubbleId), e, bubbleGraph, BubbleGraph) {
            const BubbleGraphEdge& edge = bubbleGraph[e];
            const uint64_t diagonal = edge.matrix[0][0] + edge.matrix[1][1];
            const uint64_t offDiagonal = edge.matrix[0][1] + edge.matrix[1][0];
            const uint64_t concordant = max(diagonal, offDiagonal);
            const uint64_t discordant = min(diagonal, offDiagonal);
            bubble.concordantSum += concordant;
            bubble.discordantSum += discordant;
        }
    }
}



void Bubbles::removeBadBubbles(double discordantRatioThreshold)
{
    uint64_t removedCount = 0;
    for(uint64_t bubbleId=0; bubbleId<bubbles.size(); bubbleId++) {
        if(bubbles[bubbleId].discordantRatio() > discordantRatioThreshold) {
            const BubbleGraph::vertex_descriptor v = bubbleGraph.vertexDescriptor(bubbleId);
            clear_vertex(v, bubbleGraph);
            remove_vertex(v, bubbleGraph);
            bubbleGraph.vertexTable[bubbleId] = BubbleGraph::null_vertex();
            ++removedCount;
        }
    }

    cout << "Removed " << removedCount << " bad bubbles out of " <<
        bubbles.size() << " total." << endl;



}



// Phase bubbles and reads.
// This should be called after bad bubbles were already removed
// from the bubble graph.
void Bubbles::phase(
    uint64_t iterationCount,
    uint64_t randomSeed,
    double minRelativePhase)
{
    bubbleGraph.computeConnectedComponents();
    cout << "Found " << bubbleGraph.connectedComponents.size() <<
        " connected components of the bubble graph." << endl;

    std::mt19937 randomSource(randomSeed);
    for(uint64_t componentId=0; componentId<bubbleGraph.connectedComponents.size(); componentId++) {
        phaseComponent(
            componentId,
            bubbleGraph.connectedComponents[componentId],
            iterationCount,
            randomSource,
            minRelativePhase);
    }
}



void Bubbles::phaseComponent(
    uint64_t componentId,
    const vector<BubbleGraph::vertex_descriptor>& component,
    uint64_t iterationCount,
    std::mt19937& randomSource,
    double minRelativePhase)
{
    cout << "Phasing a connected component with " << component.size() <<
        " bubbles." << endl;

    // Iterate.
    for(uint64_t iteration=0; iteration<iterationCount; iteration++) {
        cout << "Iteration " << iteration << endl;
        phaseComponentIteration(component, randomSource);

    }

    // Write out this component in svg format.
    writeBubbleGraphComponentSvg(componentId, component, minRelativePhase);

}



void Bubbles::phaseComponentIteration(
    const vector<BubbleGraph::vertex_descriptor>& component,
    std::mt19937& randomSource)
{
    using G = BubbleGraph;
    G& g = bubbleGraph;
    using vertex_descriptor = G::vertex_descriptor;

    const bool debug = false;

    // Set the phase to zero for all vertices in this component.
    for(const vertex_descriptor v: component) {
        g[v].phase = 0;
    }

    // A priority queue ordered by the number of concordant reads
    // that lead to each vertex.
    using QueuedPair = pair<uint64_t, vertex_descriptor>;
    std::priority_queue<
        QueuedPair,
        vector<QueuedPair>,
        OrderPairsByFirstOnly<uint64_t, vertex_descriptor> > q;

    // Pick a random vertex to start and queue it.
    std::uniform_int_distribution<uint64_t>
        uniformDistribution(0, component.size() - 1);
    const vertex_descriptor v = component[uniformDistribution(randomSource)];
    g[v].phase = 1;
    q.push(make_pair(0, v));


    // Main search loop.
    while(not q.empty()) {

        // Dequeue the vertex with the top priority.
        const QueuedPair p = q.top();
        const vertex_descriptor v0 = p.second;
        q.pop();
        BubbleGraphVertex& vertex0 = g[v0];

        // Check all edges incident on this vertex.
        BGL_FORALL_OUTEDGES(v0, e, g, G) {
            SHASTA_ASSERT(source(e, g) == v0);
            const BubbleGraphEdge& edge = g[e];
            const vertex_descriptor v1 = target(e, g);
            BubbleGraphVertex& vertex1 = g[v1];

            // Figure out the optimal phase of v1 based on the phase of v0 and
            // this edge.
            const uint64_t diagonalCount= edge.diagonalCount();
            const uint64_t offDiagonalCount= edge.offDiagonalCount();
            int64_t phase1 =  vertex0.phase;
            if(diagonalCount < offDiagonalCount) {
                phase1 = -phase1;
            }


            if(vertex1.phase == 0) {
                // We have not been here before.
                vertex1.phase = phase1;
                q.push(make_pair(edge.concordantCount(), v1));
            } else {
                // We have been here before. If disagreement, keep track of it.
                if(vertex1.phase != phase1) {
                    if(debug) {
                        cout << "Disagreement " << vertex1.phase << " " << phase1 << " " <<
                            diagonalCount << " " << offDiagonalCount << endl;
                    }
                }
            }
        }
    }


    // When getting here, every bubble in this component was phased.
    for(const vertex_descriptor v: component) {
        SHASTA_ASSERT(bubbleGraph[v].phase != 0);
    }



    // Use the phasing of the bubbles to phase the reads.
    // First, gather all the reads in this component.
    std::set<OrientedReadId> orientedReadIds;
    for(const vertex_descriptor v: component) {
        const uint64_t bubbleId = bubbleGraph[v].bubbleId;
        const Bubble& bubble = bubbles[bubbleId];
        for(uint64_t side=0; side<2; side++) {
            for(const OrientedReadId orientedReadId: bubble.orientedReadIds[side]) {
                orientedReadIds.insert(orientedReadId);
            }
        }
    }
    cout << "Found " << orientedReadIds.size() << " oriented reads in this component." << endl;



    // For each oriented read, count the bubbles on each orientation.
    for(const OrientedReadId orientedReadId: orientedReadIds) {

        // Loop over pairs(bubbleId, side) where this orientedRead appears.
        const vector <pair <uint64_t, uint64_t> >& v = orientedReadsTable[orientedReadId.getValue()];
        uint64_t phasePlusCount = 0;
        uint64_t phaseMinusCount = 0;
        for(const auto& p: v) {
            const uint64_t bubbleId = p.first;
            const uint64_t side = p.second;
            const vertex_descriptor v = bubbleGraph.vertexTable[bubbleId];
            if(v == BubbleGraph::null_vertex()) {
                continue;
            }
            const int64_t bubblePhase = bubbleGraph[v].phase;
            int64_t orientedReadPhase = bubblePhase;
            if(side == 1) {
                orientedReadPhase = -orientedReadPhase;
            }

            if(orientedReadPhase == +1) {
                ++phasePlusCount;
            } else {
                ++phaseMinusCount;
            }
        }
        if(debug) {
            cout << orientedReadId << " " << phasePlusCount << " " << phaseMinusCount << "\n";
        }
    }
}



void Bubbles::BubbleGraph::computeConnectedComponents()
{
    using boost::connected_components;
    using boost::get;
    using boost::color_map;

    using G = BubbleGraph;
    G& g = *this;

    connected_components(g,
        get(&BubbleGraphVertex::component, g),
        color_map(get(&BubbleGraphVertex::color, g)));

    // Gather the vertices in each connected component.
    std::map<uint64_t, vector<vertex_descriptor> > componentMap;
    BGL_FORALL_VERTICES(v, g, G) {
        componentMap[g[v].component].push_back(v);
    }
    connectedComponents.clear();
    for(const auto& p: componentMap) {
        connectedComponents.push_back(p.second);
    }
}



