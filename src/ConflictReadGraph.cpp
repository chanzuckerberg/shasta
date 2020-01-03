#include "ConflictReadGraph.hpp"
using namespace shasta;



// See this Wikipedia article for greedy coloring in general:
// https://en.wikipedia.org/wiki/Greedy_coloring

// This code uses the Brelaz (1979) DSatur
// algorithm described there in this section
// https://en.wikipedia.org/wiki/Greedy_coloring#Adaptive

// The main Wikipedia article for this algorithm is here:
// https://en.wikipedia.org/wiki/DSatur

// Brelaz, Daniel (April 1979),
// "New methods to color the vertices of a graph",
// Communications of the ACM, 22 (4): 251–256, doi:10.1145/359094.359101

// Many alternative approximate coloring methods are available, for
// example
// https://www.gerad.ca/~alainh/RLFPaper.pdf

void ConflictReadGraph::colorConnectedComponent(const vector<VertexId>& component)
{

}
