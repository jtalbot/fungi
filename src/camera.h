
#pragma once

#include "la.h"
#include "ray.h"
#include "sample.h"

struct Pinhole {
    Point l;
    Plane ip;
    float f;

    Pinhole(Point l, float f) : l(l), f(f) { ip = Plane(0, 0, 1, l.z); }

    Dip sample(P2 const& sensor) const {
        Dip dip = {P3(l), V3(1,0,0), V3(0,1,0),
                   sensor, P2(1,0), P2(0,1)};
        return dip;
    };

    void sample(P2 const& sensor, Dip& dip, V3& out) const {
        dip = sample(sensor);
        out = V3(sensor.u, sensor.v, f);
    }

    void sampleDir(Dip& dip, V3& out) const {
        float u = gi_random() * 2 - 1;
        float v = gi_random() * 2 - 1;

        sample(P2(u,v), dip, out);
    };

};

