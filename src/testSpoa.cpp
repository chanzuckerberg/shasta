#include "testSpoa.hpp"
using namespace shasta;

#include "spoa/spoa.hpp"

void shasta::testSpoa()
{


    std::vector<std::string> sequences = {
        "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
        "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
        "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
        "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
        "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
        "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
    };

    const spoa::AlignmentType alignmentType = spoa::AlignmentType::kNW;
    const int8_t match = 5;
    const int8_t mismatch = -4;
    const int8_t gap = -8;
    auto alignment_engine = spoa::AlignmentEngine::Create(alignmentType, match, mismatch, gap);

    spoa::Graph graph;

    for (const auto& it: sequences) {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }

    /*
    std::string consensus = graph->generate_consensus();
    fprintf(stderr, "Consensus (%zu)\n", consensus.size());
    fprintf(stderr, "%s\n", consensus.c_str());
    */

    const std::vector<std::string> msa = graph.GenerateMultipleSequenceAlignment();

    fprintf(stderr, "Multiple sequence alignment\n");
    for (const auto& it: msa) {
        fprintf(stderr, "%s\n", it.c_str());
    }

}
