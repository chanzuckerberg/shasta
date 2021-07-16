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

    // Construct a magnified version of another PngImage.
    PngImage(const PngImage&, int magnifyFactor);

    // Write a square grid.
    void writeGrid(int spacing, int red, int green, int blue);

    // Magnify this image.
    void magnify(int magnifyFactor);

    // Swap the content of this image with the content of another.
    void swap(PngImage&);
private:
    int width;
    int height;
    vector<::png_byte> data;
};

#endif

