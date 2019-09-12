
#ifndef RUNLENGTH_ANALYSIS_BINARYIO_HPP
#define RUNLENGTH_ANALYSIS_BINARYIO_HPP

#include <iostream>
#include <ostream>
#include <istream>
#include <string>
#include <vector>
#include <stdexcept>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>

using std::ostream;
using std::istream;
using std::string;
using std::cout;
using std::vector;
using std::runtime_error;


void writeStringToBinary(ostream& file, string& s);

void readStringFromBinary(istream& file, string& v, uint64_t length);


template<class T> void writeVectorToBinary(ostream& s, const vector<T>& v){
    ///
    /// Without worrying about size conversions, write any vector to a file using ostream.write
    ///

    s.write(reinterpret_cast<const char*>(v.data()), v.size()*sizeof(T));
}


template<class T> void writeValueToBinary(ostream& s, T v){
    ///
    /// Without worrying about size conversions, write any value to a file using ostream.write
    ///

    auto vTemp = v;
    s.write(reinterpret_cast<const char*>(&vTemp), sizeof(T));
}


template<class T> void readValueFromBinary(istream& s, T& v){
    ///
    /// Without worrying about size conversions, read any value from a file using istream.read
    ///
    cout << "Reading value size of: " << sizeof(T) << " at position: " << s.tellg() << '\n';
    s.read(reinterpret_cast<char*>(&v), sizeof(T));
}


template<class T> void readVectorFromBinary(istream& s, vector<T>& v, uint64_t length){
    ///
    /// Without worrying about size conversions, read any vector from a file using istream.read
    ///
    cout << "Reading vector of size: " << sizeof(T)*length << " at position: " << s.tellg() << '\n';

    v.resize(length);
    s.read(reinterpret_cast<char*>(v.data()), sizeof(T)*length);
}


void preadBytes(int fileDescriptor, char* bufferPointer, size_t bytesToRead, off_t& byteIndex);


void preadStringFromBinary(int fileDescriptor, string& s, uint64_t length, off_t& byteIndex);


template<class T> void preadValueFromBinary(int fileDescriptor,  T& v, off_t& byteIndex){
    ///
    /// Reimplementation of binary readValueFromBinary(), but with Linux pread, which is threadsafe
    ///

    size_t bytesToRead = sizeof(T);
    char* bufferPointer = reinterpret_cast<char*>(&v);

    preadBytes(fileDescriptor, bufferPointer, bytesToRead, byteIndex);
}


template<class T> void preadVectorFromBinary(int fileDescriptor, vector<T>& v, uint64_t length, off_t& byteIndex){
    ///
    /// Reimplementation of binary readVectorFromBinary(), but with Linux pread, which is threadsafe
    ///

    v.resize(length);

    size_t bytesToRead = sizeof(T)*length;
    char* bufferPointer = reinterpret_cast<char*>(v.data());

    preadBytes(fileDescriptor, bufferPointer, bytesToRead, byteIndex);
}


#endif //RUNLENGTH_ANALYSIS_BINARYIO_HPP
