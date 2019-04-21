
#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "color.h"
#include "la.h"
#include "sample.h"

inline int min(int a, int b) { return a < b ? a : b; }

inline int max(int a, int b) { return a < b ? b : a; }

struct Sensor {
    struct Pixel {
        rgba sum;
        float weight;
        unsigned int samples;

        Pixel() : sum(0, 0, 0, 1), weight(0), samples(0) {}
    };

    int xres, yres;
    float half;
    float inv_size;
    std::vector<Pixel> array;
    int samples;
    float scale;

    Sensor(int xres, int yres, float filtersize)
        : xres(xres),
          yres(yres),
          half(filtersize * 0.5f),
          inv_size(1.f / filtersize),
          array(xres * yres, Pixel()),
          samples(0),
          scale(minf(xres, yres)) {}

    void detect(P2 const& p, float wavelength, float radiance) {
        samples++;
        rgba c = Spectrum::torgb(wavelength, radiance);

        auto x = p.u;
        auto y = p.v;

        int left = max(0, ceil(x - half)),
            right = min(xres - 1, floor(x + half));
        int bottom = max(0, ceil(y - half)),
            top = min(yres - 1, floor(y + half));

        for (int iy = bottom; iy <= top; iy++) {
            for (int ix = left; ix <= right; ix++) {
                float weight = mitchell((x - ix) * inv_size) *
                               mitchell((y - iy) * inv_size) * inv_size *
                               inv_size;
                pixel(ix, iy).sum += c * weight;
                pixel(ix, iy).weight += weight;
                pixel(ix, iy).samples++;
            }
        }
    }

    void detect(int x, int y, rgba const& c) {
        samples++;

        float weight = 1;
        pixel(x, y).sum += c * weight;
        pixel(x, y).weight += weight;
        pixel(x, y).samples++;

        /*int left = max(0,ceil(x-half)), right = min(xres-1,floor(x+half));
        int bottom = max(0,ceil(y-half)), top = min(yres-1,floor(y+half));

        for(int iy = bottom; iy <= top; iy++) {
                for(int ix = left; ix <= right; ix++) {
                        float weight =
        mitchell((x-ix)*inv_size)*mitchell((y-iy)*inv_size)*inv_size*inv_size;
                        pixel(ix, iy).sum += c*weight;
                        pixel(ix, iy).weight += weight;
                        pixel(ix, iy).samples++;
                }
        }*/
    }

    float mitchell(float d) {
        static const float B = 1.f / 3;
        static const float C = 1.f / 3;
        float x = 4 * abs(d);
        if (x >= 2)
            return 0;
        else if (x >= 1)
            return ((-B - 6 * C) * x * x * x + (6 * B + 30 * C) * x * x +
                    (-12 * B - 48 * C) * x + (8 * B + 24 * C)) *
                   (1.f / 6);
        else
            return ((12 - 9 * B - 6 * C) * x * x * x +
                    (-18 + 12 * B + 6 * C) * x * x + (6 - 2 * B)) *
                   (1.f / 6);
    }

    void outputEXR(std::string const& filename);

    Pixel& pixel(int x, int y) {
        return array[y * xres + x];
    }

    P2 sample(int64_t i, int64_t samples, int& x, int& y) const {
        auto samplesPerPixel = std::div(samples, (int64_t)(xres * yres));
        if (i < samples - samplesPerPixel.rem) {
            i = i / samplesPerPixel.quot;
        } else {
            auto blocksize = (float)(xres * yres) / samplesPerPixel.rem;
            i = i - (samples - samplesPerPixel.rem);
            i = (int64_t)((i + gi_random()) * blocksize);
        }

        auto s = std::div(i, (int64_t)xres);
        x = s.rem;
        y = s.quot;
        auto p = sampleUnitSquare();
        auto uv = P2(s.rem + p.u, s.quot + p.v);

        auto u = (2 * uv.u - xres) / scale;
        auto v = (2 * uv.v - yres) / scale;

        return P2(u, -v);
    }

    /*
            bool outputPFM(std::string file_name) {
                    rgba c;

                    std::ofstream out(file_name.c_str(), std::ios::out |
       std::ios::binary);

                    if( !out.is_open() ) return false;

                    out << "PF\n";

                    out << xres << " " << yres << "\n";
                    out << -1.0 << "\n";

                    for ( int i = 0; i < yres; i++ ) {

                            for ( int j = 0; j < xres; j++ ) {

                                    c = pixel(j, i).sum;

                                    out.write( (char*)&c.r, sizeof( float ) );
                                    out.write( (char*)&c.g, sizeof( float ) );
                                    out.write( (char*)&c.b, sizeof( float ) );
                            }
                    }

                    out.close();

                    return true;
            }
    */
};
