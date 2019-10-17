#if 1
//#ifdef __linux__

#ifndef SHASTA_PNG_IMAGE_HPP
#define SHASTA_PNG_IMAGE_HPP

#include <png.h>
#include "string.hpp"
#include "vector.hpp"

namespace shasta {
    class PngImage;
}

// A simple class to write a png image.
class shasta::PngImage {
public:
    PngImage(int width, int height);
    void setPixel(int x, int y, int r, int g, int b);
    void write(const string& fileName) const;
private:
    int width;
    int height;
    vector<::png_byte> data;
};

#endif
#endif

