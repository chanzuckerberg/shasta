#include "Assembler.hpp"
#include "deduplicate.hpp"
using namespace shasta;



// Analyze oriented read paths in the marker graph and in the assembly graph.
// In the code and comments below, "segment" is synonym for
// "assembly graph edge".
void Assembler::analyzeOrientedReadPaths() const
{
    const AssemblyGraph& assemblyGraph = *assemblyGraphPointer;
    const uint64_t segmentCount = assemblyGraph.edges.size();
    cout << "The assembly graph has " << assemblyGraph.vertices.size() <<
        " vertices and " << segmentCount << " edges (segments)." << endl;



    // Each oriented read is guaranteed to correspond to a path in the marker graph.
    // But because of transitive reduction, that path does not necessarily correspond
    // to a path in the assembly graph.
    // Here we find, for each oriented read, the sequence of assembly graph edges
    // encountered along its marker graph path.
    // This vector is indexed by OrientedReadId::getValue().
    vector< vector<AssemblyGraph::EdgeId> > orientedReadSegments(2*readCount());
    vector<MarkerGraph::EdgeId> markerGraphPath;
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            vector<AssemblyGraph::EdgeId>& segments =
                orientedReadSegments[orientedReadId.getValue()];

            // Find the marker graph path of this oriented read.
            computeOrientedReadMarkerGraphPath(
                orientedReadId,
                0, uint32_t(markers.size(orientedReadId.getValue())-1),
                markerGraphPath);

            // Loop over the path.
            AssemblyGraph::EdgeId previousSegmentId =
                std::numeric_limits<AssemblyGraph::EdgeId>::max();
            for(const MarkerGraph::EdgeId markerGraphEdgeId: markerGraphPath) {

                // Get the corresponding segments.
                const span<const pair<AssemblyGraph::EdgeId, uint32_t> > v =
                    assemblyGraph.markerToAssemblyTable[markerGraphEdgeId];

                // If no segments, skip.
                if(v.size() == 0) {
                    continue;
                }

                // If detangling was used, there can be more than one,
                // and we don't want this here.
                SHASTA_ASSERT(v.size() == 1);

                // There is only one segment.
                const AssemblyGraph::EdgeId segmentId = v.front().first;

                // If same as the previous, slip.
                if(segmentId == previousSegmentId) {
                    continue;
                }

                // This is the next key assembly graph edge encountered
                // by this oriented read along its marker graph path. Store it.
                segments.push_back(segmentId);
                previousSegmentId = segmentId;
            }
        }
    }



    // Write a csv file with the sequence of assembly graph edges (segments)
    // for each oriented read.
    {
        ofstream csv("OrientedReadSegments.csv");
        for(ReadId readId=0; readId<readCount(); readId++) {
            for(Strand strand=0; strand<2; strand++) {
                const OrientedReadId orientedReadId(readId, strand);
                csv << orientedReadId << ",";

                vector<AssemblyGraph::EdgeId>& segments =
                    orientedReadSegments[orientedReadId.getValue()];
                for(const AssemblyGraph::EdgeId segmentId: segments) {
                    csv << segmentId << ",";
                }
                csv << "\n";
            }
        }
    }



    // Find segments that are encountered more than once in the sequence of any
    // oriented reads.
    vector<bool> isDuplicate(assemblyGraph.edges.size(), false);
    vector<AssemblyGraph::EdgeId> deduplicatedSegments;
    vector<uint64_t> frequency;
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const vector<AssemblyGraph::EdgeId>& segments =
                orientedReadSegments[orientedReadId.getValue()];

            // See if there are any duplicates.
            deduplicatedSegments = segments;
            deduplicateAndCount(deduplicatedSegments, frequency);
            for(uint64_t i=0; i<deduplicatedSegments.size(); i++) {
                if(frequency[i] > 1) {
                    isDuplicate[deduplicatedSegments[i]] = true;
                    /*
                    cout << "Segment " << deduplicatedSegments[i] << " encountered " <<
                        frequency[i] << " times by oriented read " <<
                        orientedReadId << endl;
                    */
                }
            }
        }
    }
    cout << "The following assembly graph edges (segments) are encountered "
        "more than once by one or more oriented reads:" << endl;
    uint64_t duplicateCount = 0;
    for(AssemblyGraph::EdgeId segmentId=0; segmentId<segmentCount; segmentId++) {
        if(isDuplicate[segmentId]) {
            cout << segmentId << " ";
            ++duplicateCount;
        }
    }
    cout << endl;
    cout << duplicateCount << " such segments found." << endl;



    // Table to contain, for each segment, its occurrences in
    // orientedReadSegments.
    // For each segmentId, we store a vector of pairs (orientedReadId, index)
    // such that
    // orientedReadSegments[orientedReadId.getValue()][index] == segmentId
    vector< vector< pair<OrientedReadId, uint64_t> > >
        segmentTable(segmentCount);
    for(ReadId readId=0; readId<readCount(); readId++) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const vector<AssemblyGraph::EdgeId>& segments =
                orientedReadSegments[orientedReadId.getValue()];

            for(uint64_t index=0; index<segments.size(); index++) {
                const AssemblyGraph::EdgeId segmentId = segments[index];
                segmentTable[segmentId].push_back(make_pair(orientedReadId, index));
            }
        }
    }



    // Write out the segment table.
    {
        ofstream csv("SegmentTable.csv");
        csv << "Segment,OrientedRead,Position\n";
        for(AssemblyGraph::EdgeId segmentId=0; segmentId<segmentCount; segmentId++) {
            const vector< pair<OrientedReadId, uint64_t> >& v = segmentTable[segmentId];
            for(const auto& p: v) {
                csv << segmentId << ",";
                csv << p.first << ",";
                csv << p.second << "\n";
            }

        }
    }
}
