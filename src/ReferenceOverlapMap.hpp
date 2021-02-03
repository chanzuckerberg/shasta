#ifndef SHASTA_REFERENCEOVERLAPMAP_HPP
#define SHASTA_REFERENCEOVERLAPMAP_HPP

/// The overlap map is based on a boost interval map. The interval map performs an "aggregation" operation whenever
/// multiple intervals share space on the number line, combining values for those key:value pairs, and splitting the
/// intervals at all boundaries. The overlap map uses this data structure for each chromosome/contig in the reference
/// alignment to infer overlap between reads


#include "ReadId.hpp"

#include <boost/icl/split_interval_map.hpp>
#include <boost/icl/interval_map.hpp>
#include <boost/icl/interval.hpp>

using boost::icl::total_enricher;
using boost::icl::interval_map;
using boost::icl::interval;

#include <unordered_map>
#include <string>
#include <set>

using std::unordered_map;
using std::string;
using std::set;


namespace shasta{
    class ReferenceOverlapMap;
}

class shasta::ReferenceOverlapMap {
public:
    unordered_map <string, interval_map <uint32_t, set<OrientedReadId>, total_enricher> > intervals;
    size_t size;

    void insert(string& region_name, uint32_t start, uint32_t stop, OrientedReadId id);

    ReferenceOverlapMap();

    void print(ostream& out);
};



#endif //SHASTA_REFERENCEOVERLAPMAP_HPP
