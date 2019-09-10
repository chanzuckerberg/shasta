#include "Assembler.hpp"
#include "deduplicate.hpp"
using namespace shasta;


void Assembler::createPhasingData(
    size_t threadCount,
    double phasingSimilarityThreshold,
    int maxNeighborCount)
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    // Store information used by the phasing graph to create
    // binary data.
    phasingData.dataFileNamePrefix = largeDataFileNamePrefix;
    phasingData.dataPageSize = largeDataPageSize;

    // Find the oriented reads internal to each assembly graph edge.
    phasingGatherOrientedReads(threadCount);

    // Find forks in the assembly graph.
    assemblyGraph.forks.createNew(
        largeDataName("PhasingForks"), largeDataPageSize);
    assemblyGraph.createForks();



    // Count the number of times each pair of oriented
    // reads appears on the same branch or different branches of a fork.
    // The bool is true if the two reads are on the same branch.
    vector< pair< array<OrientedReadId, 2>, bool> > readPairs;
    for(uint64_t forkId=0; forkId<assemblyGraph.forks.size(); forkId++) {

        // Get the branches (edges) of this fork.
        MemoryAsContainer<AssemblyGraph::EdgeId> edgeIds = assemblyGraph.getForkEdges(forkId);
        const uint64_t branchCount = edgeIds.size();
        SHASTA_ASSERT(branchCount);
        // cout << "Working on a fork with " << branchCount << " branches." << endl;

        // Find oriented reads that are internal to exactly one branch.
        // For each oriented read we also store the branch index.
        vector< pair<OrientedReadId, uint64_t> > orientedReadTable;
        for(uint64_t branchId=0; branchId<branchCount; branchId++) {
            const AssemblyGraph::EdgeId edgeId = edgeIds[branchId];

            // Access the oriented reads internal to this branch.
            const MemoryAsContainer<OrientedReadId> orientedReadIds =
                phasingData.orientedReads[edgeId];
            /*
            cout << "Branch " << branchId << " (edge " << edgeId << ")"
                << " has " << orientedReadIds.size() <<
                " oriented reads." << endl;
            */

            // Store them in the oriented read table.
            for(const OrientedReadId orientedReadId: orientedReadIds) {
                orientedReadTable.push_back(make_pair(orientedReadId, branchId));
            }
        }


        // Sort, then remove the ones that appear in more than one branch.
        sort(orientedReadTable.begin(), orientedReadTable.end());
        /*
        cout << "Initial size of oriented read table " << orientedReadTable.size() << endl;
        for(const auto& p: orientedReadTable) {
            cout << p.first << " " << p.second << endl;
        }
        */

        const auto begin = orientedReadTable.begin();
        const auto end = orientedReadTable.end();
        auto input = begin;
        auto output = begin;
        while(input != end) {

            if((input+1==end) || ((input+1)->first != input->first)) {
                // cout << "Storing " << input->first << " " << input->second << endl;
                *output++ = *input++;
            } else {
                auto it = input;
                while(it!=end && it->first==input->first) {
                    ++it;
                }
                input = it;
                // cout << "Skipped to " << input->first << " " << input->second << endl;
            }
        }
        orientedReadTable.resize(output - begin);
        /*
        cout << "Final size of oriented read table " << orientedReadTable.size() << endl;
        for(const auto& p: orientedReadTable) {
            cout << p.first << " " << p.second << endl;
        }
        */


        // Sanity check.
        for(uint64_t i=1; i<orientedReadTable.size(); i++) {
            SHASTA_ASSERT(orientedReadTable[i-1].first != orientedReadTable[i].first);
        }


        // Now we can loop over pairs of reads.
        array<OrientedReadId, 2> orientedReadIdPair;
        for(uint64_t i0=0; i0<orientedReadTable.size()-1; i0++) {
            const auto& p0 = orientedReadTable[i0];
            orientedReadIdPair[0] = p0.first;
            const uint64_t branchId0 = p0.second;
            for(uint64_t i1=i0+1; i1<orientedReadTable.size(); i1++) {
                const auto& p1 = orientedReadTable[i1];
                orientedReadIdPair[1] = p1.first;
                const uint64_t branchId1 = p1.second;
                SHASTA_ASSERT(orientedReadIdPair[0] != orientedReadIdPair[1]);

                readPairs.push_back(make_pair(orientedReadIdPair, branchId0==branchId1));
            }
        }
    }


    // Deduplicate and count the oriented read id pairs.
    vector<uint64_t> frequency;
    deduplicateAndCount(readPairs, frequency);
    /*
    {
        ofstream csv("ReadPairs.csv");
        for(size_t i=0; i<readPairs.size(); i++) {
            const auto& p = readPairs[i];
            csv << p.first[0] << ",";
            csv << p.first[1] << ",";
            csv << int(p.second) << ",";
            csv << frequency[i] << "\n";
        }
    }
    */



    // Gather same branch/different branch frequencies for each pair.
    phasingData.forkStatistics.createNew(
        largeDataName("PhasingForkStatistics"), largeDataPageSize);
    for(size_t i=0; i<readPairs.size(); i++) {
        const auto& p = readPairs[i];
        const uint32_t f = uint32_t(frequency[i]);

        // See if we can add it to the last one.
        if(phasingData.forkStatistics.size()>0) {
            auto& last = phasingData.forkStatistics[phasingData.forkStatistics.size()-1];
            if(p.first == last.orientedReadIds) {
                // cout << "Added to " << last.orientedReadIds[0] << " " << last.orientedReadIds[1] << endl;
                if(p.second) {
                    last.sameBranchFrequency += f;
                } else {
                    last.differentBranchFrequency += f;
                }
                continue;
            }
        }

        // Create a new one.
        PhasingData::ForkStatistics statistics;
        statistics.orientedReadIds = p.first;
        if(p.second) {
            statistics.sameBranchFrequency += f;
        } else {
            statistics.differentBranchFrequency += f;
        }
        phasingData.forkStatistics.push_back(statistics);
        // cout << "Stored " << statistics.orientedReadIds[0] << " " << statistics.orientedReadIds[1] << endl;
    }
    {
        ofstream csv("ForkStatistics.csv");
        for(const PhasingData::ForkStatistics& s: phasingData.forkStatistics) {
            csv << s.orientedReadIds[0] << ",";
            csv << s.orientedReadIds[1] << ",";
            csv << s.sameBranchFrequency << ",";
            csv << s.differentBranchFrequency << "\n";
        }
    }



    // Write a graph.
    {
        ofstream dot("Graph.dot");
        dot << "graph G {\n";
        for(const PhasingData::ForkStatistics& s: phasingData.forkStatistics) {
            if(s.sameBranchFrequency>=4 && s.differentBranchFrequency<3) {
                dot << s.orientedReadIds[0].getValue() << "--";
                dot << s.orientedReadIds[1].getValue() << ";\n";
            }
        }
        dot << "}\n";
    }



#if 0
    // Find the assembly graph edges that each oriented read is internal to.
    phasingGatherAssemblyGraphEdges(threadCount);
    phasingSortAssemblyGraphEdges(threadCount);

    // Find the assembly graph edges related to each assembly graph edge.
    // Two assembly graph edges are related if they share at
    // least one internal oriented read.
    phasingData.gatherRelatedAssemblyGraphEdges();

    // Find oriented read pairs with phasing similarity greater than the threshold.
    phasingData.findSimilarPairs(threadCount, phasingSimilarityThreshold);

    // Only keep up to maxNeighborCount neighbors.
    phasingData.keepBestSimilarPairs(maxNeighborCount);

    // Write out the global phasing graph in graphviz format.
    phasingData.writeGraphviz();
#endif
}



void Assembler::accessPhasingData()
{
    phasingData.orientedReads.accessExistingReadOnly(
        largeDataName("PhasingGraphOrientedReads"));
    phasingData.assemblyGraphEdges.accessExistingReadOnly(
        largeDataName("PhasingGraphAssemblyGraphEdges"));
    phasingData.relatedAssemblyGraphEdges.accessExistingReadOnly(
        largeDataName("PhasingRelatedAssemblyGraphEdges"));
}



// Find the oriented reads internal to each assembly graph edge.
void Assembler::phasingGatherOrientedReads(size_t threadCount)
{
    phasingData.orientedReads.createNew(
        largeDataName("PhasingGraphOrientedReads"), largeDataPageSize);
    phasingData.orientedReads.beginPass1(assemblyGraph.edges.size());
    setupLoadBalancing(assemblyGraph.edges.size(), 1000);
    runThreads(&Assembler::phasingGatherOrientedReadsPass1, threadCount);
    phasingData.orientedReads.beginPass2();
    setupLoadBalancing(assemblyGraph.edges.size(), 1000);
    runThreads(&Assembler::phasingGatherOrientedReadsPass2, threadCount);
    phasingData.orientedReads.endPass2();
}



void Assembler::phasingGatherOrientedReadsPass1(size_t threadId)
{
    phasingGatherOrientedReadsPass(1);
}
void Assembler::phasingGatherOrientedReadsPass2(size_t threadId)
{
    phasingGatherOrientedReadsPass(2);
}
void Assembler::phasingGatherOrientedReadsPass(int pass)
{

    // Define this here to reduce memory allocation activity.
    vector<OrientedReadId> orientedReadIds;

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while (getNextBatch(begin, end)) {

        // Loop over assembly graph edges assigned to this batch.
        for (AssemblyGraph::EdgeId assemblyGraphEdgeId = begin;
            assemblyGraphEdgeId != end; ++assemblyGraphEdgeId) {

            // Access the marker graph edges corresponding
            // to this assembly graph edge.
            const MemoryAsContainer<MarkerGraph::EdgeId> markerGraphEdges =
                assemblyGraph.edgeLists[assemblyGraphEdgeId];
            const uint64_t n = markerGraphEdges.size();
            SHASTA_ASSERT(n > 0);

            // Gather the oriented read ids.
            orientedReadIds.clear();
            if (n == 1) {

                // There is only one marker graph edge.
                // The OrientedReadId's are the ones in that one edge.
                const MarkerGraph::EdgeId markerGraphEdgeId =
                    markerGraphEdges[0];
                const MemoryAsContainer<MarkerInterval> markerIntervals =
                    markerGraph.edgeMarkerIntervals[markerGraphEdgeId];
                for (const MarkerInterval markerInterval : markerIntervals) {
                    orientedReadIds.push_back(markerInterval.orientedReadId);
                }

            } else {

                // The OrientedReadId's are the ones from all of the
                // marker graph vertices internal to this chain.
                // These are the target vertices of every edge except the last in the chain.
                for (uint64_t i = 0; i < n - 1; i++) {
                    const MarkerGraph::EdgeId markerGraphEdgeId =
                        markerGraphEdges[i];
                    const MarkerGraph::Edge &markedGraphEdge =
                        markerGraph.edges[markerGraphEdgeId];
                    const MarkerGraph::VertexId markerGraphVertexId =
                        markedGraphEdge.target;

                    // Loop over the markers in this marker graph vertex.
                    const MemoryAsContainer<MarkerId> markerIds =
                        markerGraph.vertices[markerGraphVertexId];
                    for (const MarkerId markerId : markerIds) {
                        OrientedReadId orientedReadId;
                        tie(orientedReadId, ignore) = findMarkerId(markerId);
                        orientedReadIds.push_back(orientedReadId);
                    }

                }
            }

            // Deduplicate.
            sort(orientedReadIds.begin(), orientedReadIds.end(),
                std::greater<OrientedReadId>());
            orientedReadIds.resize(
                unique(orientedReadIds.begin(), orientedReadIds.end())
                    - orientedReadIds.begin());

            // Store.
            if (pass == 1) {
                phasingData.orientedReads.incrementCount(assemblyGraphEdgeId,
                    orientedReadIds.size());
            } else {
                for (const OrientedReadId orientedReadId : orientedReadIds) {
                    phasingData.orientedReads.store(assemblyGraphEdgeId,
                        orientedReadId);
                }
            }

        }

    }

}



// Find the assembly graph edges that each oriented read is internal to..
void Assembler::phasingGatherAssemblyGraphEdges(size_t threadCount)
{
    const uint64_t orientedReadCount = 2 * reads.size();

    phasingData.assemblyGraphEdges.createNew(
        largeDataName("PhasingGraphAssemblyGraphEdges"), largeDataPageSize);
    phasingData.assemblyGraphEdges.beginPass1(orientedReadCount);
    setupLoadBalancing(assemblyGraph.edges.size(), 1000);
    runThreads(&Assembler::phasingGatherAssemblyGraphEdgesPass1, threadCount);
    phasingData.assemblyGraphEdges.beginPass2();
    setupLoadBalancing(assemblyGraph.edges.size(), 1000);
    runThreads(&Assembler::phasingGatherAssemblyGraphEdgesPass2, threadCount);
    phasingData.assemblyGraphEdges.endPass2();
}

void Assembler::phasingGatherAssemblyGraphEdgesPass1(size_t threadId)
{
    phasingGatherAssemblyGraphEdgesPass(1);
}
void Assembler::phasingGatherAssemblyGraphEdgesPass2(size_t threadId)
{
    phasingGatherAssemblyGraphEdgesPass(2);
}



void Assembler::phasingGatherAssemblyGraphEdgesPass(int pass)
{

    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while (getNextBatch(begin, end)) {

        // Loop over assembly graph edges assigned to this batch.
        for (AssemblyGraph::EdgeId assemblyGraphEdgeId = begin;
            assemblyGraphEdgeId != end; ++assemblyGraphEdgeId) {

            // Access the oriented reads internal to this assembly graph edge.
            const MemoryAsContainer<OrientedReadId> orientedReadIds =
                phasingData.orientedReads[assemblyGraphEdgeId];

            // Loop over these oriented reads.
            for (const OrientedReadId orientedReadId : orientedReadIds) {
                if (pass == 1) {
                    phasingData.assemblyGraphEdges.incrementCountMultithreaded(
                        orientedReadId.getValue());
                } else {
                    phasingData.assemblyGraphEdges.storeMultithreaded(
                        orientedReadId.getValue(), assemblyGraphEdgeId);
                }
            }

        }

    }

}


void Assembler::phasingSortAssemblyGraphEdges(size_t threadCount)
{
    const uint64_t orientedReadCount = 2 * reads.size();
    setupLoadBalancing(orientedReadCount, 1000);
    runThreads(&Assembler::phasingSortAssemblyGraphEdgesThreadFunction,
        threadCount);

}
void Assembler::phasingSortAssemblyGraphEdgesThreadFunction(size_t threadId)
{
    // Loop over all batches assigned to this thread.
    uint64_t begin, end;
    while (getNextBatch(begin, end)) {

        // Loop over oriented reads assigned to this batch.
        for (uint64_t orientedReadId = begin; orientedReadId != end;
            ++orientedReadId) {

            //  Sort the assembly graph edges that this oriented read is internal to.
            MemoryAsContainer<AssemblyGraph::EdgeId> edges =
                phasingData.assemblyGraphEdges[orientedReadId];
            sort(edges.begin(), edges.end());
        }
    }

}



#if 0
double Assembler::computePhasingSimilarity(
    ReadId readId0, Strand strand0,
    ReadId readId1, Strand strand1)
{
    return computePhasingSimilarity(
        OrientedReadId(readId0, strand0),
        OrientedReadId(readId1, strand1)
        );
}
double Assembler::computePhasingSimilarity(
    OrientedReadId orientedReadId0,
    OrientedReadId orientedReadId1)
{
    return phasingData.computePhasingSimilarity(orientedReadId0, orientedReadId1);
}
#endif



double Assembler::computePhasingSimilarity(
    AssemblyGraph::EdgeId edgeId0,
    AssemblyGraph::EdgeId edgeId1)
{
    return phasingData.computePhasingSimilarity(edgeId0, edgeId1);
}



uint64_t Assembler::countCommonInternalOrientedReads(
    AssemblyGraph::EdgeId edgeId0,
    AssemblyGraph::EdgeId edgeId1)
{
    return phasingData.countCommonInternalOrientedReads(edgeId0, edgeId1);
}

