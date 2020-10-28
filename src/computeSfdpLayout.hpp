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
#include <map>
#include "utility.hpp"


namespace shasta {

    enum class ComputeSfdpLayoutReturnCode {
        Success,
        SfdpError,
        Timeout,
        Signal
    };

    // Use Graphviz to compute the sfdp layout of a Boost graph.
    template<class Graph> ComputeSfdpLayoutReturnCode computeSfdpLayout(
        const Graph&,
        double timeout,
        std::map<typename Graph::vertex_descriptor, array<double, 2> >& positionMap);

}



template<class Graph> shasta::ComputeSfdpLayoutReturnCode shasta::computeSfdpLayout(
    const Graph& graph,
    double timeout,
    std::map<typename Graph::vertex_descriptor, array<double, 2> >& positionMap)
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
        "layout=sfdp;\n"
        "smoothing=triangle;\n"
        "node [shape=point];\n";

    // Write the vertices.
    for(uint64_t i=0; i<vertexCount; i++) {
        dotFile << i << ";\n";
    }

    // Write the edges.
    BGL_FORALL_EDGES_T(v, graph, Graph) {
        const vertex_descriptor v0 = source(v, graph);
        const vertex_descriptor v1 = target(v, graph);
        dotFile << vertexIndexMap[v0] << "--" << vertexIndexMap[v1] << ";\n";
    }

    dotFile << "}\n";
    dotFile.close();



    // Now use sfdp to compute the layout.
    // Use plain format output described here
    // https://www.graphviz.org/doc/info/output.html#d:plain
    const string plainFileName = dotFileName + ".txt";
    const string command = "sfdp -T plain " + dotFileName + " -o " + plainFileName;
    bool timeoutTriggered = false;
    bool signalOccurred = false;
    int returnCode = 0;
    runCommandWithTimeout(command, timeout, timeoutTriggered, signalOccurred, returnCode);
    filesystem::remove(dotFileName);
    if(signalOccurred) {
        return ComputeSfdpLayoutReturnCode::Signal;
    }
    if(timeoutTriggered) {
        return ComputeSfdpLayoutReturnCode::Timeout;
    }
    if(returnCode!=0 ) {
        return ComputeSfdpLayoutReturnCode::SfdpError;
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
    filesystem::remove(plainFileName);

    return ComputeSfdpLayoutReturnCode::Success;

}

#endif

