
#pragma once

#include <algorithm>
#include <vector>

#include "box.h"
#include "la.h"
#include "light.h"
#include "material.h"
#include "ray.h"
#include "shape.h"

struct Instance {
    std::shared_ptr<const Shape> shape;
    std::shared_ptr<const Material> material;
    Transform transform;

    Instance const* min(Ray const& r, Dip& dip, float& t) const {
        if (shape->min(transform * r, dip, t)) {
            dip = (~transform) * dip;
            return this;
        }
        return nullptr;
    }

    bool any(Ray const& r, float t) const {
        return shape->any(transform * r, t);
    }

    Box bb() const { return (~transform) * shape->bb(); }
};

struct LightInstance {
    std::shared_ptr<const Light> light;
    Transform transform;

    rgba eval(Point const& p) const { return light->eval(transform * p); }

    std::pair<rgba, Point> sample(Dip const& d) const {
        auto r = light->sample(d);
        return std::make_pair(r.first, (~transform) * r.second);
    }
};

class Scene {
public:
    Scene(std::vector<Instance> instances, std::vector<LightInstance> lights)
        : instances(std::move(instances)), lights(std::move(lights)) {
        // double t = 0;
        std::vector<Box> b;
        b.reserve(this->instances.size());

        for (auto const& i : this->instances) {
            // area really should be scaled by transform...
            // t += (*i)->light*(*i)->w();
            // ecdf.push_back(t);

            // light |= (*i)->light;

            Box const& box = i.bb();
            b.push_back(box);
        }

        bvh.construct(b);
    }

    Instance const* min(Ray const& r, Dip& dip, float& t) const {
        return min(bvh.n[0], r, 1.f/r.d, dip, t);
    }

    bool any(Ray const& r, float t) const {
        return any(bvh.n[0], r, 1.f/r.d, t);
    }

    /*Dip r() const
{
            float s = gi_random() * ecdf.back();
            size_t i = std::upper_bound(ecdf.begin(), ecdf.end(), s) -
ecdf.begin(); Dip dip = shapes[i]->r(); dip.p = dip.p * (ecdf.back() / (i == 0 ?
ecdf[0] : (ecdf[i] - ecdf[i-1]))); return dip;
    }

float w() const {
    return ecdf.back();
}*/

    Box bb() const { return bvh.root(); }

    std::pair<rgba, Point> sampleL(Dip const& d) const {
        return lights.at(0).sample(d);
    }

    rgba L(Point const& p) const {
        rgba r(0, 0, 0, 1);
        for (auto const& l : lights) {
            r += l.eval(p);
        }
        return r;
    }

private:
    Instance const* min(BVH::Node const& node, Ray const& r, V3 const& id,
                           Dip& dip, float& t) const {
        if(!node.box.intersect(r.o, id, t))
            return nullptr;

        if (node.end > 0) {
            Instance const* out = nullptr;
            for (int i = node.begin; i < node.end; i++) {
                if (auto result = instances[bvh.o[i]].min(r, dip, t))
                    out = result;
            }
            return out;
        } else {
            auto left = min(bvh.n[node.begin], r, id, dip, t);
            auto right = min(bvh.n[node.begin + 1], r, id, dip, t);
            return right ? right : left;
        }
    }

    bool any(BVH::Node const& node, Ray const& r, V3 const& id, float t) const {
        if(!node.box.intersect(r.o, id, t))
            return false;

        if (node.end > 0) {
            for (int i = node.begin; i < node.end; i++) {
                if (instances[bvh.o[i]].any(r, t)) return true;
            }
            return false;
        } else {
            return any(bvh.n[node.begin], r, id, t) ||
                   any(bvh.n[node.begin + 1], r, id, t);
        }
    }

    const std::vector<Instance> instances;
    const std::vector<LightInstance> lights;
    BVH bvh;

    // std::vector<float> ecdf;
};

