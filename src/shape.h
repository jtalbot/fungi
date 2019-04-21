
#pragma once

#include <algorithm>
#include <vector>

#include "box.h"
#include "la.h"
#include "ray.h"

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
        uint32_t a, b, c;
    };

    struct Metrics {
        int64_t rays;
        int64_t tritests;
        int64_t bbtests;
    };

    Mesh(std::vector<P3> vertexes,
         std::vector<V3> normals,
         std::vector<P2> uvs,
         std::vector<Triangle> triangles);
    ~Mesh();

    bool min(Ray const& r, Dip& dip, float& t) const override final;
    bool any(Ray const& r, float t) const override final;

    Dip r() const override final;
    float w() const override final;
    Box bb() const override final;

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

/*class Sphere : public Shape {
public:
    Sphere(float radius)
        : radius(radius) {}

    bool min(Ray const& r, Dip& dip, float& t) const override final {
        auto v = P3(0,0,0)-P3(r.o);
        // project
        auto s = dot(v,V3(r.d));
        // d is the vector from the origin to the nearest point on the ray
        auto d = v-V3(r.d)*s;
        // the distance from the origin to the ray
        auto d2 = dot(d,d);

        // is the distance from the origin to the ray less than the radius?
        if(d2 > radius*radius) {
            return false;
        }

        // chord is the distance from the
        auto chord = sqrt(radius*radius - d2);
        // I think add/sub the chord and s is assuming a normalized ray.
        auto t1 = s-chord;
        auto t2 = s+chord;

        if(t1 > 0 && t1 < t) {
            auto hit = r.o+r.d*t1;
            auto normal = ;
            frame = ;
            t = t1;
            return true;
        }
        else if(t2 > 0 && t2 < t) {
            auto hit = r.o+r.d*t2;
            auto normal = ;
            frame = ;
            t = t2;
            return true;
        }

        return false;
    }

    bool any(Ray const& r, float t) const override final {
        auto v = P3(0,0,0)-P3(r.o);
        auto s = v * V3(r.d);
        auto d = v-V3(r.d)*s;
        auto d2 = d*d;

        if(d2 > radius*radius) {
            return false;
        }

        auto chord = sqrt(radius*radius - d2);
        auto t1 = s-chord;
        auto t2 = s+chord;

        if(t1 > 0 && t1 < t) {
            return true;
        }
        else if(t2 > 0 && t2 < t) {
            return true;
        }

        return false;
    }

    Dip r() const override final {
    }

    float w() const override final {
    }

    Box bb() const override final {

    }

private:
    float radius;
};*/

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
