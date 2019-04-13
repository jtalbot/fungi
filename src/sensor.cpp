
#include "sensor.h"

/*#include <OpenEXR/ImfRgbaFile.h>
#include <OpenEXR/ImfChromaticities.h>
#include <OpenEXR/ImfStandardAttributes.h>

using namespace Imf;
using namespace Imath;

void Sensor::outputEXR(std::string const& filename) {
    std::vector<Rgba> exr_pixels(xres*yres);
    float s = ((float)xres*yres)/samples;
    for(int i = 0; i < xres*yres; i++) {
        exr_pixels[i].r = array[i].sum.r * s;
        exr_pixels[i].g = array[i].sum.g * s;
        exr_pixels[i].b = array[i].sum.b * s;
        exr_pixels[i].a = 1;
    }
    Header header(xres, yres);
    Chromaticities chroma(V2f(1,0), V2f(0,1), V2f(0,0), V2f(1./3, 1./3));
    addChromaticities(header, chroma);
    RgbaOutputFile file(file_name.c_str(), header, WRITE_RGBA);
    file.setFrameBuffer(&exr_pixels[0], 1, xres);
    file.writePixels(yres);
}*/

#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"

// See `examples/rgbe2exr/` for more details.
void Sensor::outputEXR(std::string const& filename)
{
    EXRHeader header;
    InitEXRHeader(&header);

    EXRImage image;
    InitEXRImage(&image);

    image.num_channels = 3;

    int64_t pixels = xres*yres;

    std::vector<float> images[3];
    images[0].resize(pixels);
    images[1].resize(pixels);
    images[2].resize(pixels);

    // Split RGBRGBRGB... into R, G and B layer
    float s = ((float)xres*yres)/samples;
    for (int i = 0; i < pixels; i++) {
      images[0][i] = array[i].sum.r * s;
      images[1][i] = array[i].sum.g * s;
      images[2][i] = array[i].sum.b * s;
    }

    float* image_ptr[3];
    image_ptr[0] = &(images[2].at(0)); // B
    image_ptr[1] = &(images[1].at(0)); // G
    image_ptr[2] = &(images[0].at(0)); // R

    image.images = (unsigned char**)image_ptr;
    image.width = xres;
    image.height = yres;

    header.num_channels = 3;
    header.channels = (EXRChannelInfo *)malloc(sizeof(EXRChannelInfo) * header.num_channels); 
    // Must be (A)BGR order, since most of EXR viewers expect this channel order.
    strncpy(header.channels[0].name, "B", 255); header.channels[0].name[strlen("B")] = '\0';
    strncpy(header.channels[1].name, "G", 255); header.channels[1].name[strlen("G")] = '\0';
    strncpy(header.channels[2].name, "R", 255); header.channels[2].name[strlen("R")] = '\0';

    header.pixel_types = (int *)malloc(sizeof(int) * header.num_channels); 
    header.requested_pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
    for (int i = 0; i < header.num_channels; i++) {
      header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
      header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF; // pixel type of output image to be stored in .EXR
    }

    const char* err = NULL; // or nullptr in C++11 or later.
    int ret = SaveEXRImageToFile(&image, &header, filename.c_str(), &err);
    if (ret != TINYEXR_SUCCESS) {
      fprintf(stderr, "Save EXR err: %s\n", err);
      FreeEXRErrorMessage(err); // free's buffer for an error message 
    }
    //printf("Saved exr file. [ %s ] \n", filename.c_str());

    free(header.channels);
    free(header.pixel_types);
    free(header.requested_pixel_types);
}

