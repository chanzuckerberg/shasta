
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


void write_string_to_binary(ostream& file, string& s);

void read_string_from_binary(istream& file, string& v, uint64_t length);


template<class T> void write_vector_to_binary(ostream& s, const vector<T>& v){
    ///
    /// Without worrying about size conversions, write any vector to a file using ostream.write
    ///

    s.write(reinterpret_cast<const char*>(v.data()), v.size()*sizeof(T));
}


template<class T> void write_value_to_binary(ostream& s, T v){
    ///
    /// Without worrying about size conversions, write any value to a file using ostream.write
    ///

    auto v_temp = v;
    s.write(reinterpret_cast<const char*>(&v_temp), sizeof(T));
}


template<class T> void read_value_from_binary(istream& s, T& v){
    ///
    /// Without worrying about size conversions, read any value from a file using istream.read
    ///
    cout << "Reading value size of: " << sizeof(T) << " at position: " << s.tellg() << '\n';
    s.read(reinterpret_cast<char*>(&v), sizeof(T));
}


template<class T> void read_vector_from_binary(istream& s, vector<T>& v, uint64_t length){
    ///
    /// Without worrying about size conversions, read any vector from a file using istream.read
    ///
    cout << "Reading vector of size: " << sizeof(T)*length << " at position: " << s.tellg() << '\n';

    v.resize(length);
    s.read(reinterpret_cast<char*>(v.data()), sizeof(T)*length);
}


void pread_bytes(int file_descriptor, char* buffer_pointer, size_t bytes_to_read, off_t& byte_index);


void pread_string_from_binary(int file_descriptor, string& s, uint64_t length, off_t& byte_index);


template<class T> void pread_value_from_binary(int file_descriptor,  T& v, off_t& byte_index){
    ///
    /// Reimplementation of binary read_value_from_binary(), but with Linux pread, which is threadsafe
    ///

    size_t bytes_to_read = sizeof(T);
    char* buffer_pointer = reinterpret_cast<char*>(&v);

    pread_bytes(file_descriptor, buffer_pointer, bytes_to_read, byte_index);
}


template<class T> void pread_vector_from_binary(int file_descriptor, vector<T>& v, uint64_t length, off_t& byte_index){
    ///
    /// Reimplementation of binary read_vector_from_binary(), but with Linux pread, which is threadsafe
    ///

    v.resize(length);

    size_t bytes_to_read = sizeof(T)*length;
    char* buffer_pointer = reinterpret_cast<char*>(v.data());

    pread_bytes(file_descriptor, buffer_pointer, bytes_to_read, byte_index);
}


#endif //RUNLENGTH_ANALYSIS_BINARYIO_HPP
