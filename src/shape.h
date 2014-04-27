
#ifndef _SHAPE_H_
#define _SHAPE_H_

#include <vector>
#include <algorithm>

#include "la.h"
#include "ray.h"
#include "box.h"

inline float gi_random() {
	return ((float) rand()) / ((float)RAND_MAX+1);
}

class Shape {
public:
    virtual ~Shape() {}

    virtual bool min(Ray const& r, Dip& dip, float& t) const = 0;
    virtual bool any(Ray const& r, float t) const = 0;
    
    virtual Dip r() const = 0;
    virtual float w() const = 0;

    virtual Box bb() const = 0;
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
		    t += (*i)->w();
		    ecdf.push_back(t);
			
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


class Mesh : public Shape {
public:
	struct Triangle {
		int a, b, c;

        Triangle(int a, int b, int c)
            : a(a), b(b), c(c) {}

		bool any(Mesh const& m, Ray const& r, float const t) const {
			Point i = Transform(Point(m.v[b])-m.v[a], Point(m.v[c])-m.v[a], -r.d, Point(0,0,0,1)) * (r.o-m.v[a]);
			return i.x >= 0 && i.y >= 0 && i.x+i.y < 1 &&
                   i.z > 0.00001f && i.z < t;
		}

		bool min(Mesh const& m, Ray const& r, Dip& dip, float& t) const {
			Transform x(Point(m.v[b])-m.v[a], Point(m.v[c])-m.v[a], -r.d, Point(0,0,0,1));
			Point i = x * (r.o-m.v[a]);
			if(i.x >= 0 && i.y >= 0 && i.x+i.y < 1 &&
               i.z > 0.00001f && i.z < t) {
                t = i.z;
                dip = (Dip){
                    m.v[a]*(1-i.x-i.y) + m.v[b]*i.x + m.v[c]*i.y,
				    m.v[a]*m.v[b]*m.v[c] };
				return true;
			}
			return false;
		}

		Dip r(Mesh const& m) const {
			// Use transform that maintains stratification
            float u = gi_random(), v = gi_random();
            Point i( (v>u) ? 0.5*u : u-0.5*v,
                     (v>u) ? v-0.5*u : 0.5*v,
                     0, 1 );

			Dip dip = {
                m.v[a]*(1-i.x-i.y) + m.v[b]*i.x + m.v[c]*i.y,
                m.v[a]*m.v[b]*m.v[c] };
			return dip;
		} 
	};


	Mesh(std::vector<P3> vertexes, std::vector<Triangle> triangles)
        : v(vertexes)
        , m(triangles)
        , rays(0)
        , tritests(0)
        , bbtests(0) {
    
        float sum = 0;
		std::vector<Box> b;
		b.reserve(m.size());

		for(auto i = m.begin(); i != m.end(); ++i) {
		
		    sum += (v[i->a]*v[i->b]*v[i->c]).q();
		    ecdf.push_back(sum);

            b.push_back(Box()
                .insert(v[i->a])
                .insert(v[i->b])
                .insert(v[i->c]));
		}

        printf("Starting BVH\n");
		bvh.construct(b);
		printf("Ending BVH\n");
	}

	~Mesh() {
		printf("tri/ray: %f\n", (float)tritests/rays);
		printf("bb/ray:  %f\n", (float)bbtests/rays);
	}
	
	bool min(Ray const& r, Dip& dip, float& t) const {
		rays++;
		Point id(1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z, 0);
		if(bvh.root().intersect(r.o, id, t)) {
			return minBVH(bvh.n[0], r, id, dip, t);
		}
		return false;
	}

	bool any(Ray const& r, float t) const {
		rays++;
		Point id(1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z, 0);
		if(bvh.root().intersect(r.o, id, t)) {
			return anyBVH(bvh.n[0], r, id, t);
		}
		return false;
	}

	Dip r() const {
		float s = gi_random() * ecdf.back();
        size_t i = std::upper_bound(ecdf.begin(), ecdf.end(), s) - ecdf.begin();
		Dip dip = m[i].r(*this);
		dip.p = dip.p * (ecdf.back() / ((i == 0) ? ecdf[0] : (ecdf[i]-ecdf[i-1])));
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
				tritests++;
				if(m[bvh.o[i]].min(*this, r, dip, t)) any = true;
			}
			return any;
		}
		else {
			bbtests+=2;
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
				tritests++;
				if(m[bvh.o[i]].any(*this, r, t)) return true;
			}
			return false;
		}
		else {
			bbtests+=2;
			return 	(bvh.n[node.begin].box.intersect(r.o, id, t) &&
					anyBVH(bvh.n[node.begin], r, id, t)) || 
				(bvh.n[node.begin+1].box.intersect(r.o, id, t) &&
					anyBVH(bvh.n[node.begin+1], r, id, t));
		}
	}

private:

	std::vector<P3> v;        // vertices
	std::vector<Triangle> m;  // triangles
	BVH bvh;
	
	std::vector<float> ecdf;

	mutable int rays;
	mutable int tritests;
	mutable int bbtests;
};

#endif

