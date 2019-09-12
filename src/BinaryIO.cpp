
#ifndef RUNLENGTH_ANALYSIS_BINARYIO_CPP_H
#define RUNLENGTH_ANALYSIS_BINARYIO_CPP_H

#include "BinaryIO.hpp"


void writeStringToBinary(ostream& file, string& s){
    ///
    /// Without worrying about size conversions, write any string to a file using ostream.write
    ///

    file.write(reinterpret_cast<const char*>(s.data()), s.size());
}


void readStringFromBinary(istream& file, string& s, uint64_t length){
    ///
    /// Without worrying about size conversions, read any value to a file using ostream.write
    ///

    s.resize(length);
    file.read(reinterpret_cast<char*>(&s[0]), length);
}


void preadBytes(int fileDescriptor, char* bufferPointer, size_t bytesToRead, off_t& byteIndex){
    ///
    /// Reimplementation of binary readBytes(), but with Linux pread, which is threadsafe
    ///

    while (bytesToRead) {
        const ssize_t byteCount = ::pread(fileDescriptor, bufferPointer, bytesToRead, byteIndex);
        if (byteCount <= 0) {
            throw runtime_error("ERROR " + std::to_string(errno) + " while reading: " + string(::strerror(errno)));
        }
        bytesToRead -= byteCount;
        bufferPointer += byteCount;
        byteIndex += byteCount;
    }
}


void preadStringFromBinary(int fileDescriptor, string& s, uint64_t length, off_t& byteIndex){
    ///
    /// Reimplementation of binary readStringFromBinary(), but with Linux pread, which is threadsafe
    ///

    s.resize(length);

    size_t bytesToRead = length;
    char* bufferPointer = reinterpret_cast<char*>(&s[0]);

    preadBytes(fileDescriptor, bufferPointer, bytesToRead, byteIndex);
}

#endif //RUNLENGTH_ANALYSIS_BINARYIO_CPP_H
