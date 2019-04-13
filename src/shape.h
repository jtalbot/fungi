
#pragma once

#include <algorithm>
#include <vector>

#include "box.h"
#include "la.h"
#include "ray.h"
#include "sample.h"

class Shape {
   public:
    Shape() {}

    virtual ~Shape() {}

    virtual bool min(Ray const& r, Dip& dip, float& t) const = 0;
    virtual bool any(Ray const& r, float t) const = 0;

    virtual Dip r() const = 0;
    virtual float w() const = 0;

    virtual Box bb() const = 0;
};

class Mesh : public Shape {
   public:
    struct Triangle {
        int a, b, c;

        Triangle() {}

        Triangle(int a, int b, int c) : a(a), b(b), c(c) {}

        Point intersect(Mesh const& m, Ray const& r) const {
            auto x = Transform{Point(m.v[b]) - m.v[a], Point(m.v[c]) - m.v[a],
                               -r.d, Point(0, 0, 0, 1)};

            return x * (r.o - m.v[a]);
        }

        static bool isInside(Point const& i, float const t) {
            return i.x >= 0 && i.y >= 0 && i.x + i.y < 1 && i.z > 0.00001f &&
                   i.z < t;
        }

        Dip sample(Mesh const& m, Point const& i) const {
            Plane p = m.v[a] * m.v[b] * m.v[c];

            V3 n = m.n.size() > 0
                       ? interpolate(m.n[a], m.n[b], m.n[c], i.x, i.y)
                       : V3(p.a, p.b, p.c);
            P2 uv = m.uv.size() > 0
                        ? interpolate(m.uv[a], m.uv[b], m.uv[c], i.x, i.y)
                        : P2(i.x, i.y);

            auto dtdu = m.uv.size() > 0
                            ? P2(m.uv[b].u - m.uv[a].u, m.uv[b].v - m.uv[a].v)
                            : P2(1, 0);

            auto dtdv = m.uv.size() > 0
                            ? P2(m.uv[c].u - m.uv[a].u, m.uv[c].v - m.uv[a].v)
                            : P2(0, 1);

            auto didu = V3(m.v[b].x - m.v[a].x, m.v[b].y - m.v[a].y,
                           m.v[b].z - m.v[a].z);
            auto didv = V3(m.v[c].x - m.v[a].x, m.v[c].y - m.v[a].y,
                           m.v[c].z - m.v[a].z);

            return Dip{interpolate(m.v[a], m.v[b], m.v[c], i.x, i.y),
                       p,
                       n,
                       uv,
                       didu,
                       didv,
                       dtdu,
                       dtdv};
        }

        Dip r(Mesh const& m) const {
            // Use transform that maintains stratification
            float u = gi_random(), v = gi_random();
            Point i((v > u) ? 0.5 * u : u - 0.5 * v,
                    (v > u) ? v - 0.5 * u : 0.5 * v, 0, 1);

            return sample(m, i);
        }
    };

    Mesh(std::vector<P3> vertexes, std::vector<V3> normals, std::vector<P2> uvs,
         std::vector<Triangle> triangles)
        : Shape(),
          v(std::move(vertexes)),
          n(std::move(normals)),
          uv(std::move(uvs)),
          m(std::move(triangles)),
          rays(0),
          tritests(0),
          bbtests(0) {
        std::vector<Box> b;
        b.reserve(m.size());
        ecdf.reserve(m.size());

        float sum = 0;
        for (auto const& t : m) {
            sum += (v[t.a] * v[t.b] * v[t.c]).q();
            ecdf.push_back(sum);

            b.push_back(Box().insert(v[t.a]).insert(v[t.b]).insert(v[t.c]));
        }

        // printf("Starting BVH\n");
        bvh.construct(b);
        // printf("Ending BVH\n");
    }

    //~Mesh() {
    //	printf("tri/ray: %f\n", (float)tritests/rays);
    //	printf("bb/ray:  %f\n", (float)bbtests/rays);
    //}

    bool min(Ray const& r, Dip& dip, float& t) const {
        rays++;
        Point id(1.f / r.d.x, 1.f / r.d.y, 1.f / r.d.z, 0);
        if (bvh.root().intersect(r.o, id, t)) {
            return minBVH(bvh.n[0], r, id, dip, t);
        }
        return false;
    }

    bool any(Ray const& r, float t) const {
        rays++;
        Point id(1.f / r.d.x, 1.f / r.d.y, 1.f / r.d.z, 0);
        if (bvh.root().intersect(r.o, id, t)) {
            return anyBVH(bvh.n[0], r, id, t);
        }
        return false;
    }

    Dip r() const {
        float s = gi_random() * ecdf.back();
        size_t i = std::upper_bound(ecdf.begin(), ecdf.end(), s) - ecdf.begin();
        Dip dip = m[i].r(*this);
        dip.p = dip.p *
                (ecdf.back() / ((i == 0) ? ecdf[0] : (ecdf[i] - ecdf[i - 1])));
        return dip;
    }

    float w() const { return ecdf.back(); }

    Box bb() const { return bvh.root(); }

   private:
    bool minBVH(BVH::Node const& node, Ray const& r, Point const& id, Dip& dip,
                float& t) const {
        if (node.end > 0) {
            bool any = false;
            for (int i = node.begin; i < node.end; i++) {
                tritests++;
                auto intersection = m[bvh.o[i]].intersect(*this, r);
                if (Triangle::isInside(intersection, t)) {
                    t = intersection.z;
                    dip = m[bvh.o[i]].sample(*this, intersection);
                    any = true;
                }
            }
            return any;
        } else {
            bbtests += 2;
            bool left = (bvh.n[node.begin].box.intersect(r.o, id, t) &&
                         minBVH(bvh.n[node.begin], r, id, dip, t));
            bool right = (bvh.n[node.begin + 1].box.intersect(r.o, id, t) &&
                          minBVH(bvh.n[node.begin + 1], r, id, dip, t));
            return left || right;
        }
    }

    bool anyBVH(BVH::Node const& node, Ray const& r, Point const& id,
                float t) const {
        if (node.end > 0) {
            for (int i = node.begin; i < node.end; i++) {
                tritests++;
                if (Triangle::isInside(m[bvh.o[i]].intersect(*this, r), t))
                    return true;
            }
            return false;
        } else {
            bbtests += 2;
            return (bvh.n[node.begin].box.intersect(r.o, id, t) &&
                    anyBVH(bvh.n[node.begin], r, id, t)) ||
                   (bvh.n[node.begin + 1].box.intersect(r.o, id, t) &&
                    anyBVH(bvh.n[node.begin + 1], r, id, t));
        }
    }

    const std::vector<P3> v;        // vertices
    const std::vector<V3> n;        // might be empty
    const std::vector<P2> uv;       // might be empty
    const std::vector<Triangle> m;  // triangles
    BVH bvh;

    std::vector<float> ecdf;

    mutable int64_t rays;
    mutable int64_t tritests;
    mutable int64_t bbtests;
};

/*class Transformation : public Shape {
public:
    Transformation(Transform transform, Shape const* shape)
        : transform(transform)
        , shape(shape) {
        light = shape->light;
    }

    ~Transformation() {
        delete shape;
    }

    bool min(Ray const& r, Dip& dip, float& t) const {
        bool out = shape->min(transform*r, dip, t);
        if(out) {
            dip = (~transform)*dip;
        }
        return out;
    }

    bool any(Ray const& r, float t) const {
        return shape->any(transform*r, t);
    }

    Dip r() const {
        return (~transform) * shape->r();
    }

    float w() const {
        // this is only true if the transform only uniformly scales...
        // otherwise, it's only an approximation, which is fine,
        // since it's only guiding our importance sampling.
        // If the det is 0, then we need to up the weight to avoid
        // biasing the results...
        return maxf((~transform).det, 0.00001) * shape->w();
    }

    Box bb() const {
        return (~transform) * shape->bb();
    }

private:
    Transform transform;
    Shape const* shape;
};

class Group : public Shape {
public:
        Group(std::vector<Shape const*> shapes)
        : shapes(shapes) {

        double t = 0;
                std::vector<Box> b;
                b.reserve(shapes.size());

        for(auto i = shapes.begin(); i != shapes.end(); ++i) {
            // area really should be scaled by transform...
                    t += (*i)->light*(*i)->w();
                    ecdf.push_back(t);

            light |= (*i)->light;

            Box const& box = (*i)->bb();
                        b.push_back(box);
        }

                bvh.construct(b);
    }

        ~Group() {
        for(auto i = shapes.begin(); i != shapes.end(); ++i)
            delete *i;
        }

        bool min(Ray const& r, Dip& dip, float& t) const {
                Point id(1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z, 0);
                if(bvh.root().intersect(r.o, id, t)) {
                        return minBVH(bvh.n[0], r, id, dip, t);
                }
                return false;
        }

        bool any(Ray const& r, float t) const {
                Point id(1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z, 0);
                if(bvh.root().intersect(r.o, id, t)) {
                        return anyBVH(bvh.n[0], r, id, t);
                }
                return false;
        }

        Dip r() const {
                float s = gi_random() * ecdf.back();
                size_t i = std::upper_bound(ecdf.begin(), ecdf.end(), s) -
ecdf.begin(); Dip dip = shapes[i]->r(); dip.p = dip.p * (ecdf.back() / (i == 0 ?
ecdf[0] : (ecdf[i] - ecdf[i-1]))); return dip;
        }

    float w() const {
        return ecdf.back();
    }

    Box bb() const {
        return bvh.root();
    }

        bool minBVH(BVH::Node const& node, Ray const& r, Point const& id, Dip&
dip, float& t) const { if(node.end > 0) { bool any = false; for(int i =
node.begin; i < node.end; i++) { if(shapes[bvh.o[i]]->min(r, dip, t)) any =
true;
                        }
                        return any;
                }
                else {
                        bool left = (bvh.n[node.begin].box.intersect(r.o, id, t)
&& minBVH(bvh.n[node.begin], r, id, dip, t)); bool right =
(bvh.n[node.begin+1].box.intersect(r.o, id, t) && minBVH(bvh.n[node.begin+1], r,
id, dip, t)); return left || right;
                }
        }

        bool anyBVH(BVH::Node const& node, Ray const& r, Point const& id, float
t) const { if(node.end > 0) { for(int i = node.begin; i < node.end; i++) {
                                if(shapes[bvh.o[i]]->any(r, t)) return true;
                        }
                        return false;
                }
                else {
                        return 	(bvh.n[node.begin].box.intersect(r.o, id, t) &&
                                        anyBVH(bvh.n[node.begin], r, id, t)) ||
                                (bvh.n[node.begin+1].box.intersect(r.o, id, t)
&& anyBVH(bvh.n[node.begin+1], r, id, t));
                }
        }

private:
        std::vector<Shape const*> shapes;
        BVH bvh;

        std::vector<float> ecdf;
};*/
