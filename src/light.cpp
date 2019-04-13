
#include "light.h"

#include <cmath>
#include <iostream>

#include "la.h"
#include "sample.h"
#include "tinyexr.h"

/*void test(InfiniteAreaLight const* l, int width, int height)
{
    // reproject spherical coordinates to concentric mapping
    rgba* out = new rgba[width*height];

    for(int y = 0; y < height; ++y)
    {
        for(int x = 0; x < width; ++x)
        {
            P2 p =
            V3 v = sampleHemisphere(P2(), 0);
            out[y*width+x] = l->eval(v);
        }
    }
}*/

InfiniteAreaLight::InfiniteAreaLight(rgba l, std::string const& filename)
    : l(std::move(l)), data(nullptr), width(0), height(0) {
    if (filename != "") {
        readExr(filename);
        buildDistribution();
    }
}

InfiniteAreaLight::~InfiniteAreaLight() { free(data); }

void InfiniteAreaLight::readExr(std::string const& filename) {
    // std::cout << "Reading EXR: " << filename << std::endl;
    const char* err = nullptr;

    int ret = LoadEXR(&data, &width, &height, filename.c_str(), &err);

    if (ret != TINYEXR_SUCCESS) {
        if (err) {
            fprintf(stderr, "ERR : %s\n", err);
            FreeEXRErrorMessage(err);  // release memory of error message.
        }
    }
}

void InfiniteAreaLight::buildDistribution() {
    float ccdf = 0;
    for (int y = 0; y < height; ++y) {
        std::vector<float> row;
        row.reserve(width);

        float cdf = 0;
        for (int x = 0; x < width; ++x) {
            int index = (y * width + x) * 4;
            auto color =
                rgba(data[index + 0], data[index + 1], data[index + 2], 1);
            float w = (color.r + color.g + color.b) / 3.0 + 0.00001;
            cdf += w;
            row.push_back(cdf);
        }

        cdf_columns.push_back(std::move(row));
        ccdf += cdf;
        cdf_rows.push_back(ccdf);
    }
}

std::pair<rgba, Point> InfiniteAreaLight::sample(Dip const& dip) const {
    // sample row
    auto v = std::lower_bound(cdf_rows.begin(), cdf_rows.end(),
                              cdf_rows.back() * gi_random()) -
             cdf_rows.begin();
    auto const& row = cdf_columns.at(v);

    // sample column
    auto u =
        std::lower_bound(row.begin(), row.end(), row.back() * gi_random()) -
        row.begin();

    int index = (v * width + u) * 4;

    auto l = rgba(data[index + 0], data[index + 1], data[index + 2], 1);

    auto scale = (l.r + l.g + l.b) / 3.0;

    l *= rgba(1 / scale, 1 / scale, 1 / scale, 1 / scale);

    auto phi = ((float)u / width) * 2 * M_PI;
    auto theta = ((float)v / height) * M_PI;

    auto p = Point(std::sin(theta) * std::cos(phi),
                   std::sin(theta) * std::sin(phi), std::cos(theta), 0);

    return std::make_pair(l, p);
}

rgba InfiniteAreaLight::eval(Point const& p) const {
    auto q = normalize(V3(p));

    auto phi = std::atan2(q.x, q.y);
    if (phi < 0) phi = phi + 2 * M_PI;
    auto theta = std::acos(std::max(std::min(q.z, 1.0f), -1.0f));

    int u = (phi / (2 * M_PI)) * width;
    int v = (theta / M_PI) * height;

    int index = (v * width + u) * 4;

    return rgba(data[index + 0], data[index + 1], data[index + 2], 1);
}
