
#pragma once

#include <cmath>

#include "la.h"

inline float gi_random() { return ((float)rand()) / ((float)RAND_MAX + 1); }

inline P2 sampleUnitSquare() { return P2(gi_random(), gi_random()); }

inline P2 sampleTwoSquare() {
    return P2(2 * gi_random() - 1, 2 * gi_random() - 1);
}

inline P2 sampleDisk(P2 const& r) {
    const auto u = r.u * 2 - 1;
    const auto v = r.v * 2 - 1;

    if (std::abs(u) > std::abs(v)) {
        const float t = M_PI / 4 * (v / u);
        return P2(u * std::cos(t), u * std::sin(t));
    } else if (std::abs(u) < std::abs(v)) {
        const float t = M_PI / 4 * (u / v);
        return P2(v * std::sin(t), v * std::cos(t));
    } else {
        const float t = M_PI / 4;
        return P2(u * std::cos(t), u * std::cos(t));
    }
}

inline V3 sampleHemisphere(P2 const& r, float cosPower) {
    auto d = sampleDisk(r);

    float r2 = d.u * d.u + d.v * d.v;
    float z = std::pow(1 - r2, 1.f / (cosPower + 1));
    float w = r2 > 0 ? sqrt((1 - z * z) / r2) : 0;

    return V3(d.u * w, d.v * w, z);
}

inline V3 sampleHemisphere(float cosPower) {
    return sampleHemisphere(sampleUnitSquare(), cosPower);
}
