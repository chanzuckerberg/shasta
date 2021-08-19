#ifndef SHASTA_GFA_ASSEMBLY_GRAPH_HPP
#define SHASTA_GFA_ASSEMBLY_GRAPH_HPP

// A class used to facilitate output of an assembly graph in GFA 1 format.
// https://gfa-spec.github.io/GFA-spec/GFA1.html

// It is implemented as a Boost directed graph in which
// each edge corresponds to a GFA segment and vertices are used to
// generate GFA links.

// The Vertex can be any type that can be used as a key in an std map.

// Shasta.
#include "Base.hpp"
#include "SHASTA_ASSERT.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

// Standard library.
#include "cstdint.hpp"
#include "fstream.hpp"
#include <map>
#include "string.hpp"
#include "vector.hpp"



namespace shasta {

    class GfaAssemblyGraphEdge;

    template<class Vertex> using GfaAssemblyGraphBaseClass =
        boost::adjacency_list<
        boost::listS, boost::listS, boost::bidirectionalS,
        Vertex, GfaAssemblyGraphEdge>;

    template <class Vertex> class GfaAssemblyGraph;

}


// Each edge corresponds to a GFA segment.
class shasta::GfaAssemblyGraphEdge {
public:
    string name;

    bool sequenceIsAvailable = false;;
    vector<Base> sequence;          // Only valid if the above is true.

    bool sequenceLengthIsAvailable = false;;
    uint64_t sequenceLength = 0;    // Only valid if the above is true.

    GfaAssemblyGraphEdge(
        const string& name
        ) :
        name(name)
    {}

    GfaAssemblyGraphEdge(
        const string& name,
        uint64_t sequenceLength
        ) :
        name(name),
        sequenceLengthIsAvailable(true),
        sequenceLength(sequenceLength)
    {}

    GfaAssemblyGraphEdge(
        const string& name,
        const vector<Base>& sequence
        ) :
        name(name),
        sequenceIsAvailable(true),
        sequence(sequence),
        sequenceLengthIsAvailable(true),
        sequenceLength(uint64_t(sequence.size()))
    {}
};



template <class Vertex> class shasta::GfaAssemblyGraph :
    public GfaAssemblyGraphBaseClass<Vertex> {
public:

    // Add a segment with known sequence.
    void addSegment(
        const string& name,
        Vertex vertex0,
        Vertex vertex1,
        const vector<Base>& sequence
    )
    {
        const vertex_descriptor v0 = getVertex(vertex0);
        const vertex_descriptor v1 = getVertex(vertex1);
        boost::add_edge(v0, v1, GfaAssemblyGraphEdge(name, sequence), *this);
    }

    // Add a segment with unknown sequence,
    // but with known sequence length.
    void addSegment(
        const string& name,
        Vertex vertex0,
        Vertex vertex1,
        uint64_t sequenceLength
    )
    {
        const vertex_descriptor v0 = getVertex(vertex0);
        const vertex_descriptor v1 = getVertex(vertex1);
        boost::add_edge(v0, v1, GfaAssemblyGraphEdge(name, sequenceLength), *this);
    }

    // Add a segment with unknown sequence and sequence length.
    void addSegment(
        const string& name,
        Vertex vertex0,
        Vertex vertex1
    )
    {
        const vertex_descriptor v0 = getVertex(vertex0);
        const vertex_descriptor v1 = getVertex(vertex1);
        boost::add_edge(v0, v1, GfaAssemblyGraphEdge(name), *this);
    }

    // Add a path, specified as a sequence of segment ids
    // (without the "+" sign which is added automatically).
    void addPath(const string& name, const vector<string>& path)
    {
        paths.push_back(Path{name, path});
    }

    // Write out in GFA format.
    void write(const string& fileName) const
    {
        ofstream gfa(fileName);
        write(gfa);
    }
    void write(ostream& gfa) const
    {
        writeHeader(gfa);
        writeSegments(gfa);
        writeLinks(gfa);
        writePaths(gfa);
    }

    using BaseClass = GfaAssemblyGraphBaseClass<Vertex>;
    using vertex_descriptor = typename BaseClass::vertex_descriptor;
    using G = GfaAssemblyGraph<Vertex>;

private:

    // The paths to be written in the gfa file.
    // Each path consists of a sequence of segment ids
    // (without "+" signs, which are added automatically).
    class Path {
    public:
        string name;
        vector<string> segmentNames;
    };
    vector<Path> paths;

    static void writeHeader(ostream& gfa)
    {
        gfa << "H\tVN:Z:1.0\n";
    }

    // Write a segment for each edge.
    void writeSegments(ostream& gfa) const
    {
        const G& g = *this;

        // Write one segment for each edge.
        BGL_FORALL_EDGES_T(e, g, G) {
            const GfaAssemblyGraphEdge& edge = g[e];

            gfa << "S\t" << edge.name << "\t";

            if(edge.sequenceIsAvailable) {
                copy(edge.sequence.begin(), edge.sequence.end(),
                    ostream_iterator<Base>(gfa));
                gfa << "\tLN:i:" << edge.sequenceLength << "\n";
            } else if (edge.sequenceLengthIsAvailable) {
                gfa << "*\tLN:i:" << edge.sequenceLength << "\n";
            } else {
                gfa << "*\n";
            }
        }
    }


    // Write the GFA links.
    // For each vertex, we write a link for each pair of
    // incoming/outgoing edges.
    // This uses "0M" as the Cigar string.
    void writeLinks(ostream& gfa) const
    {
        const G& g = *this;

        BGL_FORALL_VERTICES_T(v, g, G) {
            BGL_FORALL_INEDGES_T(v, e0, g, G) {
                const string& name0 = g[e0].name;
                BGL_FORALL_OUTEDGES_T(v, e1, g, G) {
                    const string& name1 = g[e1].name;

                    gfa << "L\t" <<
                        name0 << "\t+\t" <<
                        name1 << "\t+\t0M\n";
                }
            }
        }
    }

    // Write the paths.
    void writePaths(ostream& gfa) const
    {
        for(const Path& path: paths) {
            gfa << "P\t" << path.name << "\t";
            for(uint64_t i=0; i<uint64_t(path.segmentNames.size()); i++) {
                if(i != 0) {
                    gfa << ",";
                }
                gfa << path.segmentNames[i] << "+";
            }
            gfa << "\t";
            for(uint64_t i=0; i<uint64_t(path.segmentNames.size()-1); i++) {
                if(i != 0) {
                    gfa << ",";
                }
                gfa << "0M";
            }
            gfa << "\n";

        }
    }

    // Map a user-provided Vertex to the corresponding vertex_descriptor.
    std::map<Vertex, vertex_descriptor> vertexMap;

    // Return the vertex_descriptor corresponding to a Vertex,
    // adding it if necessary.
    vertex_descriptor getVertex(const Vertex& vertex)
    {
        auto it = vertexMap.find(vertex);
        if(it == vertexMap.end()) {
            const vertex_descriptor v = boost::add_vertex(*this);
            bool wasInserted = false;
            tie(it, wasInserted) = vertexMap.insert(make_pair(vertex, v));
            SHASTA_ASSERT(wasInserted);
        }
        return it->second;
    }


};





#endif

