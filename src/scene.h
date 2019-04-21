
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

inline BVH constructSceneBVH(std::vector<Instance> const& instances) {
    std::vector<Box> b;
    b.reserve(instances.size());

    for (auto const& i : instances) {
        b.push_back(i.bb());
    }

    return BVH::construct(b);
}

class Scene {
   public:
    Scene(std::vector<Instance> instances, std::vector<LightInstance> lights)
        : instances(std::move(instances)),
          lights(std::move(lights)),
          bvh(constructSceneBVH(this->instances)) {
        // double t = 0;
        // for (auto const& i : this->instances) {
        // area really should be scaled by transform...
        // t += (*i)->light*(*i)->w();
        // ecdf.push_back(t);

        // light |= (*i)->light;
        //}
    }

    Instance const* min(Ray const& r, Dip& dip, float& t) const {
        return min(bvh.root(), r, dip, t);
    }

    bool any(Ray const& r, float t) const { return any(bvh.root(), r, t); }

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

    Box bb() const { return bvh.root().box; }

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
    Instance const* min(BVH::Node const& node, FastRay const& r, Dip& dip,
                        float& t) const {
        if (!node.box.intersect(r, t)) return nullptr;

        if (node.end > 0) {
            Instance const* out = nullptr;
            for (int i = node.begin; i < node.end; i++) {
                if (auto result = instances[bvh.o[i]].min(r, dip, t))
                    out = result;
            }
            return out;
        } else {
            auto left = min(bvh.n[node.begin], r, dip, t);
            auto right = min(bvh.n[node.begin + 1], r, dip, t);
            return right ? right : left;
        }
    }

    bool any(BVH::Node const& node, FastRay const& r, float t) const {
        if (!node.box.intersect(r, t)) return false;

        if (node.end > 0) {
            for (int i = node.begin; i < node.end; i++) {
                if (instances[bvh.o[i]].any(r, t)) return true;
            }
            return false;
        } else {
            return any(bvh.n[node.begin], r, t) ||
                   any(bvh.n[node.begin + 1], r, t);
        }
    }

    const std::vector<Instance> instances;
    const std::vector<LightInstance> lights;
    const BVH bvh;

    // std::vector<float> ecdf;
};
