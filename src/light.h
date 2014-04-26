
#ifndef LIGHT_H
#define LIGHT_H

#include "ray.h"

struct PointLight {

	Point p;

	PointLight(Point p) : p(p) {}

	Point sample() const { return p; /*float u = gi_random(); float v = gi_random(); return p + u*Point(50,0,0,0) + v*Point(0,0,50,0);*/ }
};

#endif
