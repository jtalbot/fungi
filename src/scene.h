
#pragma once

#include <vector>
#include <algorithm>

#include "la.h"
#include "ray.h"
#include "box.h"
#include "shape.h"
#include "material.h"
#include "light.h"

struct Instance {
    std::shared_ptr<const Shape> shape;
    std::shared_ptr<const Material> material;
    Transform transform;
    
    Instance const* min(Ray const& r, Dip& dip, float& t) const {
        if(shape->min(transform*r, dip, t))
        {
            dip = (~transform)*dip;
            return this;
        }
        return nullptr;
    }

    bool any(Ray const& r, float t) const {
        return shape->any(transform*r, t);
    }

    Box bb() const {
        return (~transform) * shape->bb();
    }
};

struct LightInstance {
    std::shared_ptr<const Light> light;
    Transform transform;

    rgba eval(Point const& p) const {
        return light->eval(transform*p);
    }

    std::pair<rgba, Point> sample(Dip const& d) const {
        auto r = light->sample(d);
        return std::make_pair(r.first, (~transform) * r.second);
    }
};

class Scene {
public:
    Scene(std::vector<Instance> instances, std::vector<LightInstance> lights)
        : instances(std::move(instances))
        , lights(std::move(lights))
    {
        //double t = 0;
		std::vector<Box> b;
		b.reserve(this->instances.size());
        
        for(auto const& i : this->instances)
        {
            // area really should be scaled by transform...
		    //t += (*i)->light*(*i)->w();
		    //ecdf.push_back(t);

            //light |= (*i)->light;
			
            Box const& box = i.bb();
			b.push_back(box);
        }
		
		bvh.construct(b);
    }

	Instance const* min(Ray const& r, Dip& dip, float& t) const
    {
		Point id(1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z, 0);
		if(bvh.root().intersect(r.o, id, t)) {
			return minBVH(bvh.n[0], r, id, dip, t);
		}
		return nullptr;
	}

	bool any(Ray const& r, float t) const
    {
		Point id(1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z, 0);
		if(bvh.root().intersect(r.o, id, t)) {
			return anyBVH(bvh.n[0], r, id, t);
		}
		return false;
	}

	/*Dip r() const
    {
		float s = gi_random() * ecdf.back();
		size_t i = std::upper_bound(ecdf.begin(), ecdf.end(), s) - ecdf.begin();
		Dip dip = shapes[i]->r();
		dip.p = dip.p * (ecdf.back() / (i == 0 ? ecdf[0] : (ecdf[i] - ecdf[i-1])));
		return dip;
	}

    float w() const {
        return ecdf.back();
    }*/

    Box bb() const
    {
        return bvh.root();
    }

	Instance const* minBVH(BVH::Node const& node, Ray const& r, Point const& id, Dip& dip, float& t) const
    {
		if(node.end > 0) {
			Instance const* out = nullptr;
			for(int i = node.begin; i < node.end; i++) {
				if(auto result = instances[bvh.o[i]].min(r, dip, t)) out = result;
			}
			return out;
		}
		else {
			auto left = (bvh.n[node.begin].box.intersect(r.o, id, t) ? 
					minBVH(bvh.n[node.begin], r, id, dip, t) : nullptr);
			auto right = (bvh.n[node.begin+1].box.intersect(r.o, id, t) ? 
					minBVH(bvh.n[node.begin+1], r, id, dip, t) : nullptr);
			return right ? right : left;
		}
	}

	bool anyBVH(BVH::Node const& node, Ray const& r, Point const& id, float t) const
    {
		if(node.end > 0) {
			for(int i = node.begin; i < node.end; i++) {
				if(instances[bvh.o[i]].any(r, t)) return true;
			}
            return false;
		}
		else {
            return  (bvh.n[node.begin].box.intersect(r.o, id, t) && 
					anyBVH(bvh.n[node.begin], r, id, t)) || 
                    (bvh.n[node.begin+1].box.intersect(r.o, id, t) && 
                    anyBVH(bvh.n[node.begin+1], r, id, t));
		}
	}

    std::pair<rgba, Point> sampleL(Dip const& d) const
    {
        return lights.at(0).sample(d);
    }

    rgba L(Point const& p) const {
        rgba r(0,0,0,1);
        for(auto const& l : lights)
        {
            r += l.eval(p);
        }
        return r;
    }

private:
	std::vector<Instance> instances;
    std::vector<LightInstance> lights;
	BVH bvh;
	
	//std::vector<float> ecdf;
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
		size_t i = std::upper_bound(ecdf.begin(), ecdf.end(), s) - ecdf.begin();
		Dip dip = shapes[i]->r();
		dip.p = dip.p * (ecdf.back() / (i == 0 ? ecdf[0] : (ecdf[i] - ecdf[i-1])));
		return dip;
	}

    float w() const {
        return ecdf.back();
    }

    Box bb() const {
        return bvh.root();
    }

	bool minBVH(BVH::Node const& node, Ray const& r, Point const& id, Dip& dip, float& t) const {
		if(node.end > 0) {
			bool any = false;
			for(int i = node.begin; i < node.end; i++) {
				if(shapes[bvh.o[i]]->min(r, dip, t)) any = true;
			}
			return any;
		}
		else {
			bool left = (bvh.n[node.begin].box.intersect(r.o, id, t) && 
					minBVH(bvh.n[node.begin], r, id, dip, t));
			bool right = (bvh.n[node.begin+1].box.intersect(r.o, id, t) && 
					minBVH(bvh.n[node.begin+1], r, id, dip, t));
			return left || right;
		}
	}

	bool anyBVH(BVH::Node const& node, Ray const& r, Point const& id, float t) const {
		if(node.end > 0) {
			for(int i = node.begin; i < node.end; i++) {
				if(shapes[bvh.o[i]]->any(r, t)) return true;
			}
			return false;
		}
		else {
			return 	(bvh.n[node.begin].box.intersect(r.o, id, t) && 
					anyBVH(bvh.n[node.begin], r, id, t)) || 
				(bvh.n[node.begin+1].box.intersect(r.o, id, t) && 
					anyBVH(bvh.n[node.begin+1], r, id, t));
		}
	}

private:
	std::vector<Shape const*> shapes;
	BVH bvh;
	
	std::vector<float> ecdf;
};
*/

