#include "ReferenceOverlapMap.hpp"


void shasta::ReferenceOverlapMap::insert(string& region_name, uint32_t start, uint32_t stop, OrientedReadId id) {
    // Check if this reference contig has been encountered before, and initialize it if needed
    if (intervals.count(region_name) == 0){
        interval_map <uint32_t, set<OrientedReadId>, total_enricher> empty_interval_map;
        intervals.emplace(region_name, empty_interval_map);
    }

    // Construct a Boost interval for this numeric range/coord
    auto a = interval<uint32_t>::right_open(start, stop);
    set<OrientedReadId> b = {id};

    // Within this reference contig, add the numeric interval
    intervals.at(region_name).add({a, b});
    size++;
}


void shasta::ReferenceOverlapMap::print(ostream& out){
    for (auto& item: intervals){
        out << item.first << '\n';
        for (auto& item1: item.second){
            out << item1.first << " -> ";
            for (auto& item2: item1.second){
                out << item2 << " " ;
            }
            out << '\n';
        }
    }
}


shasta::ReferenceOverlapMap::ReferenceOverlapMap():
        intervals(),
        size(0)
{}
