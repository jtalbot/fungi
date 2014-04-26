
#ifndef _RAY_H_
#define _RAY_H_

#include <vector>
#include <algorithm>

#include "la.h"

struct Ray {
	Point o, d;
} __attribute__ ((aligned));

struct Dip {
	Point i;
	Plane p;
};

#endif

