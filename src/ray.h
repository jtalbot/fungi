
#pragma once

#include <algorithm>
#include <vector>

#include "la.h"

struct Ray {
    P3 o;
    V3 d;
};

inline Ray operator*(Transform const& t, Ray const& r) {
    return {P3(t * r.o), V3(t * r.d)};
}

struct Dip {
    P3 i;
    V3 didu, didv;
    P2 uv;
    P2 dtdu, dtdv;
};

inline Dip operator*(Transform const& t, Dip const& d) {
    return {P3(t * d.i), V3(t * d.didu), V3(t * d.didv),
                d.uv, d.dtdu, d.dtdv};
}
