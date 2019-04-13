
#include "material.h"

#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

Transform makeLocal(Dip const& d) {
    auto z = normalize(d.n);

    V3 x = normalize(d.didu);
    V3 y = cross(z, x);
    x = cross(y, z);

    Transform local(y, x, z, Point(0, 0, 0, 1));

    return local;
}

class ImageTexture::Impl {
   public:
    Impl(std::string const& filename);
    rgba eval(P2 const& uv) const;

    ~Impl() {
        for (int y = 0; y < height; y++) {
            free(row_pointers[y]);
        }
        free(row_pointers);
    }

   private:
    void readPNG(std::string const& filename);
    rgba value(int x, int y) const;

    int width, height;
    png_byte color_type;
    png_byte bit_depth;
    png_byte** row_pointers;
};

ImageTexture::Impl::Impl(std::string const& filename) { readPNG(filename); }

rgba ImageTexture::Impl::eval(P2 const& uv) const {
    float u = width * uv.u;
    float v = height * (1.f - uv.v);

    int x = floor(u - 0.5);
    int y = floor(v - 0.5);

    // (0.5, 0.5) should get 100% from (0,0)

    u = u - x - 0.5;
    v = v - y - 0.5;

    rgba a = value(x, y);
    rgba b = value(x + 1, y);
    rgba c = value(x, y + 1);
    rgba d = value(x + 1, y + 1);

    rgba result = a * ((1 - u) * (1 - v)) + b * (u * (1 - v)) +
                  c * ((1 - u) * v) + d * (u * v);

    return result;
}

float correct(float v) {
    if (v <= 0.04045f)
        return v * 1.f / 12.92f;
    else
        return std::pow((v + 0.055f) * 1.f / 1.055f, 2.4f);
}

rgba ImageTexture::Impl::value(int x, int y) const {
    while (x < 0) x += width;
    while (x >= width) x -= width;
    while (y < 0) y += height;
    while (y >= height) y -= height;

    auto row = row_pointers[y];
    png_bytep p = &(row[x * 4]);
    return rgba(correct(p[0] / 255.0f), correct(p[1] / 255.0f),
                correct(p[2] / 255.0f), correct(p[3] / 255.0f));
}

void ImageTexture::Impl::readPNG(std::string const& filename) {
    /*
     * A simple libpng example program
     * http://zarb.org/~gc/html/libpng.html
     *
     * Modified by Yoshimasa Niwa to make it much simpler
     * and support all defined color_type.
     *
     * To build, use the next instruction on OS X.
     * $ brew install libpng
     * $ clang -lz -lpng15 libpng_test.c
     *
     * Copyright 2002-2010 Guillaume Cottenceau.
     *
     * This software may be freely redistributed under the terms
     * of the X11 license.
     */

    // std::cout << "Loading png: " << filename << std::endl;

    FILE* fp = fopen(filename.c_str(), "rb");

    png_structp png =
        png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) abort();

    png_infop info = png_create_info_struct(png);
    if (!info) abort();

    if (setjmp(png_jmpbuf(png))) abort();

    png_init_io(png, fp);

    png_read_info(png, info);

    width = png_get_image_width(png, info);
    height = png_get_image_height(png, info);
    color_type = png_get_color_type(png, info);
    bit_depth = png_get_bit_depth(png, info);

    // Read any color_type into 8bit depth, RGBA format.
    // See http://www.libpng.org/pub/png/libpng-manual.txt

    if (bit_depth == 16) png_set_strip_16(png);

    if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);

    // PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth.
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
        png_set_expand_gray_1_2_4_to_8(png);

    if (png_get_valid(png, info, PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png);

    // These color_type don't have an alpha channel then fill it with 0xff.
    if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_GRAY ||
        color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_filler(png, 0xFF, PNG_FILLER_AFTER);

    if (color_type == PNG_COLOR_TYPE_GRAY ||
        color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_gray_to_rgb(png);

    png_read_update_info(png, info);

    row_pointers = (png_byte**)malloc(sizeof(png_byte*) * height);
    for (int y = 0; y < height; y++) {
        row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png, info));
    }

    png_read_image(png, row_pointers);

    fclose(fp);
}

ImageTexture::ImageTexture(std::string const& filename)
    : a(std::make_shared<ImageTexture::Impl const>(filename)) {}

rgba ImageTexture::eval(P2 const& uv) const { return a->eval(uv); }

rgba BumpMaterial::eval(P2 const& uv) const {
    // return bumpMap->eval(uv);
    return m->eval(uv);
}

std::pair<rgba, V3> BumpMaterial::sample(Dip const& dip, V3 const& i) const {
    auto d0 = bumpMap->eval(dip.uv).g;
    auto d1 = bumpMap->eval(dip.uv + dip.dtdu * 0.001).g;
    auto d2 = bumpMap->eval(dip.uv + dip.dtdv * 0.001).g;

    V3 du = dip.didu + dip.n * (d1 - d0) * dip.dtdu.length() / 0.01;
    V3 dv = dip.didv + dip.n * (d2 - d0) * dip.dtdv.length() / 0.01;

    auto n = cross(du, dv);

    auto newDip = Dip{dip.i, dip.p, n, dip.uv, du, dv, dip.dtdu, dip.dtdv};

    return m->sample(newDip, i);
}

std::pair<rgba, V3> DiffuseMaterial::sample(Dip const& dip,
                                            V3 const& in) const {
    auto out = V3(~makeLocal(dip) * sampleHemisphere(1));
    return std::make_pair(eval(dip.uv), out);
}

std::pair<rgba, V3> GlassMaterial::sample(Dip const& dip, V3 const& in) const {
    auto local = makeLocal(dip);

    auto i = local * in;
    // If the input vector length is proportional to the
    // IOR of the material through which it is traveling,
    // this code will maintain that through reflection (easy)
    // and transmission (a little harder). This invariant
    // simplifies the Fresnel and Snell computations quite
    // a bit.

    // Reflectance vector
    // Staying in same medium, so vector length doesn't change.
    auto r = V3(i.x, i.y, -i.z);

    // Ratio of IORs depends on whether ray is entering or leaving object
    auto ratio = i.z < 0 ? eta / 1.f : 1.f / eta;

    // Squared length of the incoming vector parallel to the surface
    auto p = (i.x * i.x + i.y * i.y);

    // Squared length of the output vector scaled by IOR ratio
    auto l = (p + i.z * i.z) * ratio * ratio;

    // If the output length is greater than the parallel part,
    // transmission can occur, otherwise total internal reflection occurs.
    if (l > p) {
        // Compute the z component for the transmitted vector
        // with the IOR scaled output length.
        auto t = V3(i.x, i.y, std::copysign(std::sqrt(l - p), i.z));

        // Compute Fresnel reflectance, then sample accordingly
        auto r1 = (i.z * ratio - t.z / ratio) / (i.z * ratio + t.z / ratio);
        auto r2 = (i.z - t.z) / (i.z + t.z);
        auto fr = 0.5f * (r1 * r1 + r2 * r2);

        if (gi_random() < fr)
            return std::make_pair(Kr->eval(dip.uv), V3(~local * r));
        else
            return std::make_pair(Kt->eval(dip.uv), V3(~local * t));
    } else {
        // Total internal reflection
        return std::make_pair(Kr->eval(dip.uv), V3(~local * r));
    }
}

/*std::pair<rgba, V3> GlassMaterial::sample(Dip const& dip, V3 const& in) const
{
    // If the input vector length is proportional to the
    // IOR of the material through which it is traveling,
    // this code will maintain that through reflection (easy)
    // and transmission (a little harder). This invariant
    // simplifies the Fresnel and Snell computations quite
    // a bit.

    auto i = in*dip.n;          // cos(in)
    auto p = Trans(dip.n*i);
    auto q = in-p;

    // Ratio of IORs depends on whether ray is entering or leaving object
    auto r = i > 0 ? eta/1.f : i < 0 ? 1.f/eta : 1.f;

    // Squared length of the incoming vector parallel to the surface
    auto q2 = q/q;

    // Squared length of the output vector scaled by IOR ratio
    auto l2 = (in/in)*r*r;

    // If the output length is greater than the parallel part,
    // transmission can occur, otherwise total internal reflection occurs.
    if(l2 > q2)
    {
        // Compute the z component for the transmitted vector
        // with the IOR scaled output length.
        //auto t = V3(i.x, i.y, std::copysign(std::sqrt(l2-q2),i.z));
        auto o = std::sqrt(l2-q2);       // cos(out)
        auto t = q + o*p;

        // Compute Fresnel reflectance, then sample accordingly
        auto r1 = (i*r-o/r) / (i*r+o/r);
        auto r2 = (i-o) / (i+o);
        auto fr = 0.5f*(r1*r1+r2*r2);

        if(gi_random() >= fr)
            return std::make_pair(Kt->eval(dip.uv), t);
    }

    // Reflectance vector
    // Staying in same medium, so vector length doesn't change.
    auto r = q + -1*p;

    // Total internal reflection or Fresnel reflectance
    return std::make_pair(Kr->eval(dip.uv), r);
}*/
