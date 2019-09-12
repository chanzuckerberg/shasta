
#include "CompressedRunnieReader.hpp"


ostream& operator<<(ostream& s, CompressedRunnieIndex& index) {
    s << "sequenceLength:\t" << index.sequenceLength << '\n';
    s << "sequenceByteIndex:\t" << index.sequenceByteIndex << '\n';
    s << "nameLength:\t" << index.nameLength << '\n';
    s << "name:\t" << index.name;

    return s;
}


void CompressedRunnieSequence::printEncoding(){
    for (auto& element: this->encoding) {
        cout << int(element) << ',';
    }
    cout << '\n';
}


CompressedRunnieReader::CompressedRunnieReader(string filePath) {
    this->sequenceFilePath = filePath;

    // Open the input file.
    this->sequenceFileDescriptor = ::open(filePath.c_str(), O_RDONLY);

    // Verify it is working
    if(this->sequenceFileDescriptor == -1) {
        throw runtime_error("ERROR: could not read " + filePath);
    }

    // Find file size in bytes
    this->fileLength = lseek(this->sequenceFileDescriptor, 0, SEEK_END);

    // Initialize remaining parameters using the file footer data
    this->readFooter();

    // Read table of contents, needed for indexed reading
    this->readIndexes();
}


size_t CompressedRunnieReader::countReads(){
    return this->indexes.size();
}


void CompressedRunnieReader::readSequence(CompressedRunnieSequence& sequence, uint64_t readIndex){
    off_t byteIndex = off_t(this->indexes[readIndex].sequenceByteIndex);
    preadStringFromBinary(this->sequenceFileDescriptor, sequence.sequence, this->indexes[readIndex].sequenceLength, byteIndex);
    preadVectorFromBinary(this->sequenceFileDescriptor, sequence.encoding, this->indexes[readIndex].sequenceLength, byteIndex);
}


void CompressedRunnieReader::readIndexEntry(CompressedRunnieIndex& indexElement, off_t& byteIndex){
    preadValueFromBinary(this->sequenceFileDescriptor, indexElement.sequenceByteIndex, byteIndex);
    preadValueFromBinary(this->sequenceFileDescriptor, indexElement.sequenceLength, byteIndex);
    preadValueFromBinary(this->sequenceFileDescriptor, indexElement.nameLength, byteIndex);
    preadStringFromBinary(this->sequenceFileDescriptor, indexElement.name, indexElement.nameLength, byteIndex);
}


void CompressedRunnieReader::readIndexes(){
    off_t byteIndex = off_t(this->indexesStartPosition);

    while (byteIndex > 0 and uint64_t(byteIndex) < this->channelMetadataStartPosition){
        CompressedRunnieIndex indexElement;
        this->readIndexEntry(indexElement, byteIndex);
        this->indexes.emplace_back(indexElement);

        // Update the mapping of read names to their places in the vector of indexes
        auto element = make_pair(indexElement.name, this->indexes.size() - 1);
        auto success = this->indexMap.insert(move(element)).second;
        if (not success){
            throw runtime_error("ERROR: possible duplicate read name (" + indexElement.name + ") found in runnie file: " + this->sequenceFilePath);
        }
    }
}


void CompressedRunnieReader::readChannelMetadata(){
    ///
    /// Read the description of the channels that accompany the nucleotide sequence. Not currently used for initializing
    /// data vectors, since CompressedRunnieSequence.encodings specifies the data type.
    ///
    off_t byteIndex = off_t(this->channelMetadataStartPosition);
    preadValueFromBinary(this->sequenceFileDescriptor, this->nChannels, byteIndex);
    this->channelSizes.resize(nChannels);

    for (uint64_t i=0; i<this->nChannels; i++) {
        preadValueFromBinary(this->sequenceFileDescriptor, this->channelSizes[i], byteIndex);
    }
}


void CompressedRunnieReader::readFooter(){
    off_t byteIndex = off_t(this->fileLength - 2*sizeof(uint64_t));
    preadValueFromBinary(this->sequenceFileDescriptor, this->indexesStartPosition, byteIndex);
    preadValueFromBinary(this->sequenceFileDescriptor, this->channelMetadataStartPosition, byteIndex);

    this-> readChannelMetadata();
}
