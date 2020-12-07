#ifndef SHASTA_MAP_HPP
#define SHASTA_MAP_HPP

// A quick and dirty hash map with open addressing and linear probing.
// See https://en.wikipedia.org/wiki/Open_addressing
// It has several limitations and deviations from
// the behavior of standard container, including:
// - Keys and values in unoccupied buckets are default constructed.
// - We hash the bits of the key.
// - Deletion is not supported.
// - 32-bit hashes are used.

#include "chrono.hpp"
#include "vector.hpp"

namespace shasta {
    template<class K, class V> class Map;
    void testMap();
}



template<class K, class V> class shasta::Map {
public:

    Map(uint32_t log2BucketCount = 8) :
        log2BucketCount(log2BucketCount)
    {
        SHASTA_ASSERT(log2BucketCount < 32);
        const uint32_t bucketCount = (1 << log2BucketCount);
        buckets.resize(bucketCount);
        mask = bucketCount - 1;
    }

    bool insert(const K& key, const V& value)
    {
        if(size > buckets.size()/2) {
            rehash(log2BucketCount+2);
        }
        uint32_t h = hash(key);
        while(true) {
            const uint32_t bucketIndex = (h++ & mask);
            Bucket& b = buckets[bucketIndex];
            if(not b.isOccupied) {
                b.isOccupied = true;
                b.p.first = key;
                b.p.second = value;
                ++size;
                return true;
            } else {
                if(b.p.first == key) {
                    // We already have this key.
                    return false;
                }
            }
        }
    }

    pair<K, V>* find(const K& key)
    {
        uint32_t h = hash(key);
        while(true) {
            const uint32_t bucketIndex = (h++ & mask);
            Bucket& b = buckets[bucketIndex];
            if(not b.isOccupied) {
                return 0;
            } else {
                if(b.p.first == key) {
                    return &b.p;
                }
            }
        }
    }

    void rehash(uint32_t log2BucketCount)
    {
        SHASTA_ASSERT(log2BucketCount < 32);
        Map<K, V> newMap(log2BucketCount);
        for(const Bucket& bucket: buckets) {
            if(bucket.isOccupied) {
                newMap.insert(bucket.p.first, bucket.p.second);
            }
        }
        swap(newMap);
    }

    void swap(Map<K, V>& that)
    {
        buckets.swap(that.buckets);
        std::swap(size, that.size);
        std::swap(mask, that.mask);
    }



    // Note that this default constructs objects
    // in unoccupied buckets.
    class Bucket {
    public:
        pair<K, V>p;
        bool isOccupied = false;
    };
    vector<Bucket> buckets;
    uint32_t log2BucketCount;
    uint32_t mask;
    uint32_t size = 0;

    uint32_t hash(const K& k) const
    {
        return MurmurHash2(&k, sizeof(K), 45341);
    }
};



inline void shasta::testMap()
{
    int n;
    std::cin >> n;

    for(int i=0; i<3; i++) {
        Map<int, int> m(20);
        const auto t0 = steady_clock::now();
        for(int i=0; i<n; i++) {
            m.insert(i, i);
        }
        const auto t1 = steady_clock::now();
        const auto dt = seconds(t1-t0);

        cout << dt << " " << m.size << " " << dt/m.size << " " << m.buckets.size() << endl;
    }
    /*
    for(uint64_t i=0; i<m.buckets.size(); i++) {
        if(m.buckets[i].isOccupied) {
            cout << m.buckets[i].p.first <<" " << m.buckets[i].p.second << endl;
        }
    }
    */


    cout << "All good." << endl;
}

#endif
