
#ifndef RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEREADER_HPP
#define RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEREADER_HPP

#include "BinaryIO.hpp"
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <unordered_map>

using std::pair;
using std::string;
using std::to_string;
using std::cout;
using std::ofstream;
using std::vector;
using std::runtime_error;
using std::unordered_map;


class CompressedRunnieIndex {
public:
    /// Attributes ///
    string name;
    uint64_t name_length;
    uint64_t sequence_byte_index;
    uint64_t sequence_length;

    /// Methods ///
};

ostream& operator<<(ostream& s, CompressedRunnieIndex& index);


class CompressedRunnieSequence {
public:
    /// Attributes ///
    string name;
    string sequence;
    vector <uint8_t> encoding;

    /// Methods ///
    void print_encoding();
};


class CompressedRunnieReader{
public:
    /// Attributes ///
    string sequence_file_path;
    int sequence_file_descriptor;

    uint64_t indexes_start_position;
    uint64_t channel_metadata_start_position;
    off_t file_length;

    // How many accessory channels will be paired 1:1 with each nucleotide sequence
    uint64_t n_channels;

    // What is the unit size of each channel
    vector<uint64_t> channel_sizes;

    vector<CompressedRunnieIndex> indexes;
    unordered_map<string,size_t> index_map;

    /// Methods ///
    CompressedRunnieReader(string file_path);
    void read_footer();
    void read_channel_metadata();
    void read_indexes();
    void read_index_entry(CompressedRunnieIndex& index_element, off_t& byte_index);
    void read_sequence(CompressedRunnieSequence& sequence, uint64_t read_index);
};

#endif //RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEREADER_HPP
