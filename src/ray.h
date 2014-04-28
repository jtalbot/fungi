
#ifndef _RAY_H_
#define _RAY_H_

#include <vector>
#include <algorithm>

#include "la.h"

struct Ray {
	Point o, d;
} __attribute__ ((aligned));

Ray operator*(Transform const& t, Ray const& r) {
    return (Ray){t*r.o, t*r.d};
}

struct Dip {
	Point i;
	Plane p;
};

Dip operator*(Transform const& t, Dip const& d) {
    return (Dip){(~t)*d.i, d.p*(~t)};
}

#endif

