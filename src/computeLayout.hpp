#ifndef SHASTA_COMPUTE_SFDP_LAYOUT
#define SHASTA_COMPUTE_SFDP_LAYOUT

// Shasta.
#include "filesystem.hpp"
#include "platformDependent.hpp"
#include "runCommandWithTimeout.hpp"
#include "SHASTA_ASSERT.hpp"

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "array.hpp"
#include "fstream.hpp"
#include "iostream.hpp"
#include <map>
#include "utility.hpp"


namespace shasta {

    enum class ComputeLayoutReturnCode {
        Success,
        Error,
        Timeout,
        Signal
    };

    // Use Graphviz to compute the sfdp layout of a Boost graph.
    template<class Graph> ComputeLayoutReturnCode computeLayout(
        const Graph&,
        const string& layoutMethod,
        double timeout,
        std::map<typename Graph::vertex_descriptor, array<double, 2> >& positionMap,
        const string& additionalOptions="",
        const std::map<typename Graph::edge_descriptor, double>* edgeLengthMap = 0);

}


// The edge length map is only effective with neato and sfdp layouts.
template<class Graph> shasta::ComputeLayoutReturnCode shasta::computeLayout(
    const Graph& graph,
    const string& layoutMethod,
    double timeout,
    std::map<typename Graph::vertex_descriptor, array<double, 2> >& positionMap,
    const string& additionalOptions,
    const std::map<typename Graph::edge_descriptor, double>* edgeLengthMap)
{
    using vertex_descriptor = typename Graph::vertex_descriptor;

    // Create a vector of vertex descriptors and
    // a map from vertex descriptors to vertex indices.
    uint64_t i = 0;
    vector<vertex_descriptor> vertexVector;
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    BGL_FORALL_VERTICES_T(v, graph, Graph) {
        vertexVector.push_back(v);
        vertexIndexMap.insert(make_pair(v, i++));
    }
    const uint64_t vertexCount = i;

    // Create a dot file to contain the graph.
    const string uuid = to_string(boost::uuids::random_generator()());
    const string dotFileName = tmpDirectory() + uuid + ".dot";
    ofstream dotFile(dotFileName);
    dotFile <<
        "graph G {\n"
        "layout=" << layoutMethod << ";\n"
        "smoothing=triangle;\n"
        "node [shape=point];\n";

    // Write the vertices.
    for(uint64_t i=0; i<vertexCount; i++) {
        dotFile << i << ";\n";
    }

    // Write the edges.
    BGL_FORALL_EDGES_T(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        dotFile << vertexIndexMap[v0] << "--" << vertexIndexMap[v1];
        if(edgeLengthMap) {
            const auto it = edgeLengthMap->find(e);
            if(it != edgeLengthMap->end()) {
                dotFile << " [len=" << it->second << "]";
            }
        }
        dotFile<< ";\n";
    }

    dotFile << "}\n";
    dotFile.close();



    // Now use graphviz to compute the layout.
    // Use plain format output described here
    // https://www.graphviz.org/doc/info/output.html#d:plain
    const string plainFileName = dotFileName + ".txt";
    const string command = layoutMethod + " -T plain " + dotFileName + " -o " + plainFileName + " " + additionalOptions;
    bool timeoutTriggered = false;
    bool signalOccurred = false;
    int returnCode = 0;
    runCommandWithTimeout(command, timeout, timeoutTriggered, signalOccurred, returnCode);
    cout << "Dot file " << dotFileName << endl;
    // filesystem::remove(dotFileName);
    if(signalOccurred) {
        return ComputeLayoutReturnCode::Signal;
    }
    if(timeoutTriggered) {
        return ComputeLayoutReturnCode::Timeout;
    }
    if(returnCode!=0 ) {
        return ComputeLayoutReturnCode::Error;
    }



    // Extract vertex coordinates and store them in positionMap.
    positionMap.clear();
    ifstream plainFile(plainFileName);
    string line;
    vector<string> tokens;
    while(true) {

        // Read the next line.
        std::getline(plainFile, line);
        if( not plainFile) {
            break;
        }

        // Parse it.
        boost::algorithm::split(tokens, line, boost::algorithm::is_any_of(" "));

        // If not a line describing a vertex, skip it.
        if(tokens.front() != "node") {
            continue;
        }
        SHASTA_ASSERT(tokens.size() >= 4);

        // Get the vertex id.
        const string& vertexName = tokens[1];
        SHASTA_ASSERT(not vertexName.empty());
        const uint64_t vertexId = boost::lexical_cast<uint64_t>(vertexName);

        // Get the corresponding vertex descriptor.
        SHASTA_ASSERT(vertexId < vertexVector.size());
        const vertex_descriptor v = vertexVector[vertexId];

        // Store it in the layout.
        array<double, 2> x;
        x[0] = boost::lexical_cast<double>(tokens[2]);
        x[1] = boost::lexical_cast<double>(tokens[3]);
        positionMap.insert(make_pair(v, x));
    }
    plainFile.close();
    cout << plainFileName << endl;
    // filesystem::remove(plainFileName);

    return ComputeLayoutReturnCode::Success;

}

#endif

