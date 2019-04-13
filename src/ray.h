
#pragma once

#include <vector>
#include <algorithm>

#include "la.h"

struct Ray {
	Point o, d;
} __attribute__ ((aligned));

inline Ray operator*(Transform const& t, Ray const& r) {
    return (Ray){t*r.o, t*r.d};
}

struct Dip {
	Point i;
	Plane p;
    V3 n;
    P2 uv;
    V3 didu, didv;
    P2 dtdu, dtdv;
};

inline Dip operator*(Transform const& t, Dip const& d) {
    return (Dip){t*d.i, d.p*t, d.n, d.uv, V3(t*d.didu), V3(t*d.didv), d.dtdu, d.dtdv};
}

