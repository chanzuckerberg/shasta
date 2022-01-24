#ifndef SHASTA_TEE_HPP
#define SHASTA_TEE_HPP

// Class to "tee" an ostream to another ostream.
// That is, all output to the first ostream is duplicated
// to go to both ostreams.


#include <ios>


namespace shasta {
    class Tee;
}



class shasta::Tee : public std::streambuf {
public:
    std::streambuf* original = 0;
    std::streambuf* copy = 0;

    void duplicate(std::ostream& originalStream, std::ostream& copyStream)
    {
        original = originalStream.rdbuf(this);
        copy = copyStream.rdbuf();
    }

    int overflow(int c = EOF)
    {
        copy->sputc(char(c));
        return original->sputc(char(c));
    }

    int sync()
    {
        copy->pubsync();
        return original->pubsync();
    }
};


#endif

