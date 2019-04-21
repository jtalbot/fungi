
#pragma once

#include <limits>
#include "la.h"
#include "ray.h"

// Ray with inverse direction for faster intersection tests with boxes
struct FastRay : public Ray {
    V3 i;

    FastRay(Ray const& r) : Ray(r), i(1.f / r.d) {}
};

struct Box {
    P3 m, n;

    Box()
        : m(Infinity, Infinity, Infinity), n(-Infinity, -Infinity, -Infinity) {}

    bool intersect(FastRay const& r, const float t) const {
        const auto a = (m - r.o) * r.i;
        const auto b = (n - r.o) * r.i;

        float lmin = 0, lmax = t;
        lmin = maxf(lmin, minf(a.x, b.x));
        lmax = minf(lmax, maxf(a.x, b.x));

        lmin = maxf(lmin, minf(a.y, b.y));
        lmax = minf(lmax, maxf(a.y, b.y));

        lmin = maxf(lmin, minf(a.z, b.z));
        lmax = minf(lmax, maxf(a.z, b.z));

        return lmin <= lmax;
    }

    Box& insert(P3 const& p) {
        m = min(m, p);
        n = max(n, p);
        return *this;
    }

    Box& insert(Box const& b) {
        m = min(m, b.m);
        n = max(n, b.n);
        return *this;
    }
} __attribute__((aligned));

inline Box operator*(Transform const& t, Box const& b) {
    return Box()
        .insert(P3(t * b.m))
        .insert(P3(t * P3(b.n.x, b.m.y, b.m.z)))
        .insert(P3(t * P3(b.m.x, b.n.y, b.m.z)))
        .insert(P3(t * P3(b.m.x, b.m.y, b.n.z)))
        .insert(P3(t * P3(b.n.x, b.n.y, b.m.z)))
        .insert(P3(t * P3(b.m.x, b.n.y, b.n.z)))
        .insert(P3(t * P3(b.n.x, b.m.y, b.n.z)))
        .insert(P3(t * b.n));
}

class BVH {
public:
    struct Node {
        Box box;
        int begin, end;
    } __attribute__((aligned));

    const std::vector<int> o;
    const std::vector<Node> n;

    Node const& root() const { return n[0]; }

    static BVH construct(std::vector<Box> const& b);
};
