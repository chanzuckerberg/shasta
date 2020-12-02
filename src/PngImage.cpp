#include "PngImage.hpp"
using namespace shasta;

#include "stdexcept.hpp"

PngImage::PngImage(int width, int height) :
    width(width),
    height(height),
    data(3ULL*uint64_t(width)*uint64_t(height), 0)
{
}



void PngImage::setPixel(int x, int y, int r, int g, int b)
{
    ::png_byte* pointer = data.data() + 3* (width*uint64_t(y) + x);
    *pointer++ = png_byte(r);
    *pointer++ = png_byte(g);
    *pointer++ = png_byte(b);
}



void PngImage::write(const std::string& fileName) const
{
    // Open the png file.
    ::FILE* fp = std::fopen (fileName.c_str(), "w");
    if(!fp) {
        throw runtime_error("Error opening " + fileName);
    }


    
    // Initialize the PNG data structures.
    ::png_structp pngPointer = ::png_create_write_struct (PNG_LIBPNG_VER_STRING, 0, 0, 0);
    if(pngPointer == 0) {
        throw runtime_error("Error writing " + fileName);
    }
    
    ::png_infop infoPointer = ::png_create_info_struct (pngPointer);
    if(infoPointer == 0) {
        throw runtime_error("Error writing " + fileName);
    }

   ::png_set_IHDR (
        pngPointer,
        infoPointer,
        width,
        height,
        8,
        PNG_COLOR_TYPE_RGB,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT);
        

        
    // Create pointers to rows
    vector<png_byte*> rowPointers(height);
    for (int y=0; y<height; y++) {
        rowPointers[y] = const_cast<png_byte*>(&data[0] + 3*width*uint64_t(y));
    }
    
    // Write.
    ::png_init_io (pngPointer, fp);
    ::png_set_rows (pngPointer, infoPointer, rowPointers.data());
    ::png_write_png (pngPointer, infoPointer, PNG_TRANSFORM_IDENTITY, NULL);

    // Clean up.
    ::png_destroy_write_struct(&pngPointer, &infoPointer);
    std::fclose (fp);
}




// Construct a magnified version of another PngImage.
PngImage::PngImage(const PngImage& that, int magnifyFactor) :
    width(magnifyFactor * that.width),
    height(magnifyFactor * that.height),
    data(3 * width * height)
{
    for(int y=0; y<that.height; y++) {
        for(int dy=0; dy<magnifyFactor; dy++) {
            const int Y = magnifyFactor * y + dy;
            for(int x=0; x<that.width; x++) {
                const auto pixel = that.data.begin() + 3* (that.width * y + x);
                for(int dx=0; dx<magnifyFactor; dx++) {
                    const int X = magnifyFactor * x + dx;
                    setPixel(X, Y, pixel[0], pixel[1], pixel[2]);
                }
            }
        }
    }

}



// Magnify this image.
void PngImage::magnify(int magnifyFactor)
{
    PngImage magnified(*this, magnifyFactor);
    swap(magnified);
}



// Swap the content of this image with the content of another.
void PngImage::swap(PngImage& that)
{
    std::swap(width, that.width);
    std::swap(height, that.height);
    data.swap(that.data);
}


// Write a square grid.
void PngImage::writeGrid(int spacing, int red, int green, int blue)
{
    for(int x=0; x<width; x+=spacing) {
        for(int y=0; y<height; y++) {
            setPixel(x, y, red, green, blue);
        }
    }
    for(int x=0; x<width; x++) {
        for(int y=0; y<height; y+=spacing) {
            setPixel(x, y, red, green, blue);
        }
    }
}
