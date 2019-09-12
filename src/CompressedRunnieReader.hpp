
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
    uint64_t nameLength;
    uint64_t sequenceByteIndex;
    uint64_t sequenceLength;

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
    void printEncoding();
};


class CompressedRunnieReader{
public:
    /// Attributes ///
    string sequenceFilePath;
    int sequenceFileDescriptor;

    uint64_t indexesStartPosition;
    uint64_t channelMetadataStartPosition;
    off_t fileLength;

    // How many accessory channels will be paired 1:1 with each nucleotide sequence
    uint64_t nChannels;

    // What is the unit size of each channel
    vector<uint64_t> channelSizes;

    vector<CompressedRunnieIndex> indexes;
    unordered_map<string,size_t> indexMap;

    /// Methods ///
    CompressedRunnieReader(string filePath);
    void readFooter();
    void readChannelMetadata();
    void readIndexes();
    void readIndexEntry(CompressedRunnieIndex& indexElement, off_t& byteIndex);
    void readSequence(CompressedRunnieSequence& sequence, uint64_t readIndex);
    size_t countReads();
};

#endif //RUNLENGTH_ANALYSIS_COMPRESSEDRUNNIEREADER_HPP
