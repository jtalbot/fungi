
#include "shape.h"

#include <algorithm>

#include "sample.h"

BVH constructBVH(std::vector<P3> const& v,
                 std::vector<Mesh::Triangle> const& m) {
    std::vector<Box> b;
    b.reserve(m.size());

    for (auto const& t : m) {
        b.push_back(Box().insert(v[t.a]).insert(v[t.b]).insert(v[t.c]));
    }

    return BVH::construct(b);
}

std::vector<float> constructECDF(std::vector<P3> const& v,
                                 std::vector<Mesh::Triangle> const& m) {
    std::vector<float> ecdf;
    ecdf.reserve(m.size());

    float sum = 0;
    for (auto const& t : m) {
        sum += (v[t.a] * v[t.b] * v[t.c]).q();
        ecdf.push_back(sum);
    }

    return ecdf;
}

struct Mesh::Impl {
    const std::vector<P3> v;        // vertices
    const std::vector<V3> n;        // might be empty
    const std::vector<P2> uv;       // might be empty
    const std::vector<Triangle> m;  // triangles

    const BVH bvh;
    const std::vector<float> ecdf;

    Impl(std::vector<P3> vertexes, std::vector<V3> normals, std::vector<P2> uvs,
         std::vector<Triangle> triangles)
        : v(std::move(vertexes)),
          n(std::move(normals)),
          uv(std::move(uvs)),
          m(std::move(triangles)),
          bvh(constructBVH(this->v, this->m)),
          ecdf(constructECDF(this->v, this->m)) {}

    // Mesh functions
    bool min(BVH::Node const& node, FastRay const& r, Dip& dip,
             float& t) const {
        // metrics.bbtests += 1;
        if (!node.box.intersect(r, t)) return false;

        if (node.end > 0) {
            bool hit = false;
            for (int i = node.begin; i < node.end; i++) {
                auto intersection = intersect(m[bvh.o[i]], r);
                if (isInside(intersection, t)) {
                    t = intersection.z;
                    dip = sample(m[bvh.o[i]], intersection);
                    hit = true;
                }
            }
            // metrics.tritests += node.end - node.begin;
            return hit;
        } else {
            auto left = min(bvh.n[node.begin], r, dip, t);
            auto right = min(bvh.n[node.begin + 1], r, dip, t);
            return left || right;
        }
    }

    bool any(BVH::Node const& node, FastRay const& r, float t) const {
        // metrics.bbtests += 1;
        if (!node.box.intersect(r, t)) return false;

        if (node.end > 0) {
            for (int i = node.begin; i < node.end; i++) {
                auto intersection = intersect(m[bvh.o[i]], r);
                if (isInside(intersection, t)) {
                    // metrics.tritests += (i + 1) - node.begin;
                    return true;
                }
            }
            // metrics.tritests += node.end - node.begin;
            return false;
        } else {
            return any(bvh.n[node.begin], r, t) ||
                   any(bvh.n[node.begin + 1], r, t);
        }
    }

    // Triangle functions
    V3 intersect(Mesh::Triangle const& t, Ray const& r) const {
        auto x = Transform{v[t.b] - v[t.a], v[t.c] - v[t.a], -r.d,
                           Point(0, 0, 0, 1)};

        return V3(x * (r.o - v[t.a]));
    }

    static bool isInside(V3 const& i, float const t) {
        return i.x >= 0 && i.y >= 0 && i.x + i.y < 1 && i.z > 0.00001f &&
               i.z < t;
    }

    Dip sample(Triangle const& t, Point const& p) const {
        auto i = interpolate(v[t.a], v[t.b], v[t.c], p.x, p.y);

        auto didu = normalize(v[t.b] - v[t.a]);
        auto didv = normalize(v[t.c] - v[t.a]);

        auto tuv = uv.size() > 0
                       ? interpolate(uv[t.a], uv[t.b], uv[t.c], p.x, p.y)
                       : P2(p.x, p.y);

        auto dtdu = uv.size() > 0
                        ? P2(uv[t.b].u - uv[t.a].u, uv[t.b].v - uv[t.a].v)
                        : P2(1, 0);

        auto dtdv = uv.size() > 0
                        ? P2(uv[t.c].u - uv[t.a].u, uv[t.c].v - uv[t.a].v)
                        : P2(0, 1);

        return Dip{i, didu, didv, tuv, dtdu, dtdv};
    }

    Dip r(Triangle const& t) const {
        // Use transform that maintains stratification
        float u = gi_random(), v = gi_random();
        Point i((v > u) ? 0.5 * u : u - 0.5 * v,
                (v > u) ? v - 0.5 * u : 0.5 * v, 0, 1);

        return sample(t, i);
    }
};

Mesh::Mesh(std::vector<P3> vertexes, std::vector<V3> normals,
           std::vector<P2> uvs, std::vector<Triangle> triangles)
    : impl(std::make_unique<Impl>(vertexes, normals, uvs, triangles)) {}

Mesh::~Mesh() {}

Dip Mesh::r() const {
    float s = gi_random() * impl->ecdf.back();
    size_t i = std::upper_bound(impl->ecdf.begin(), impl->ecdf.end(), s) -
               impl->ecdf.begin();
    Dip dip = impl->r(impl->m[i]);
    // dip.p = dip.p *
    //        (ecdf.back() / ((i == 0) ? ecdf[0] : (ecdf[i] - ecdf[i - 1])));
    return dip;
}

float Mesh::w() const { return impl->ecdf.back(); }

Box Mesh::bb() const { return impl->bvh.root().box; }

bool Mesh::min(Ray const& r, Dip& dip, float& t) const {
    // metrics.rays++;
    return impl->min(impl->bvh.root(), r, dip, t);
}

bool Mesh::any(Ray const& r, float t) const {
    // metrics.rays++;
    return impl->any(impl->bvh.root(), r, t);
}
