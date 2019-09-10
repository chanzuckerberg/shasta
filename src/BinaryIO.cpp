
#ifndef RUNLENGTH_ANALYSIS_BINARYIO_CPP_H
#define RUNLENGTH_ANALYSIS_BINARYIO_CPP_H

#include "BinaryIO.hpp"


void write_string_to_binary(ostream& file, string& s){
    ///
    /// Without worrying about size conversions, write any string to a file using ostream.write
    ///

    file.write(reinterpret_cast<const char*>(s.data()), s.size());
}


void read_string_from_binary(istream& file, string& s, uint64_t length){
    ///
    /// Without worrying about size conversions, read any value to a file using ostream.write
    ///

    s.resize(length);
    file.read(reinterpret_cast<char*>(&s[0]), length);
}


void pread_bytes(int file_descriptor, char* buffer_pointer, size_t bytes_to_read, off_t& byte_index){
    ///
    /// Reimplementation of binary read_bytes(), but with Linux pread, which is threadsafe
    ///

    while (bytes_to_read) {
        const ssize_t byte_count = ::pread(file_descriptor, buffer_pointer, bytes_to_read, byte_index);
        if (byte_count <= 0) {
            throw runtime_error("ERROR " + std::to_string(errno) + " while reading: " + string(::strerror(errno)));
        }
        bytes_to_read -= byte_count;
        buffer_pointer += byte_count;
        byte_index += byte_count;
    }
}


void pread_string_from_binary(int file_descriptor, string& s, uint64_t length, off_t& byte_index){
    ///
    /// Reimplementation of binary read_string_from_binary(), but with Linux pread, which is threadsafe
    ///

    s.resize(length);

    size_t bytes_to_read = length;
    char* buffer_pointer = reinterpret_cast<char*>(&s[0]);

    pread_bytes(file_descriptor, buffer_pointer, bytes_to_read, byte_index);
}

#endif //RUNLENGTH_ANALYSIS_BINARYIO_CPP_H
