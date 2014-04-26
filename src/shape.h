
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
		std::vector<Point> c;
		b.reserve(shapes.size());
		c.reserve(shapes.size());
        
        for(auto i = shapes.begin(); i != shapes.end(); ++i) {
            // area really should be scaled by transform...
		    t += (*i)->w();
		    ecdf.push_back(t);
			
            Box const& box = (*i)->bb();
			b.push_back(box);
			c.push_back((box.m + box.n)*0.5);
        }
		
		bvh.construct(b,c);
    }

	~Group() {
        for(auto i = shapes.begin(); i != shapes.end(); ++i)
            delete *i;
	}
	
	bool min(Ray const& r, Dip& dip, float& t) const {
		Point id(1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z, 0);
		if(bvh.n[0].box().intersect(r.o, id, t)) {
			return minBVH(bvh.n[0], r, id, dip, t);
		}
		return false;
	}

	bool any(Ray const& r, float t) const {
		Point id(1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z, 0);
		if(bvh.n[0].box().intersect(r.o, id, t)) {
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
        return bvh.n[0].box();
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
			bool left = (bvh.n[node.begin].box().intersect(r.o, id, t) && 
					minBVH(bvh.n[node.begin], r, id, dip, t));
			bool right = (bvh.n[node.begin+1].box().intersect(r.o, id, t) && 
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
			return 	(bvh.n[node.begin].box().intersect(r.o, id, t) && 
					anyBVH(bvh.n[node.begin], r, id, t)) || 
				(bvh.n[node.begin+1].box().intersect(r.o, id, t) && 
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
		unsigned int i0, i1, i2;

        Triangle(unsigned int i0, unsigned int i1, unsigned int i2)
            : i0(i0), i1(i1), i2(i2) {}

		bool any(Mesh const* mesh, Ray const& r, float& t) const {
			Point const& p0 = mesh->vertexes[i0];
			Point const& p1 = mesh->vertexes[i1];
			Point const& p2 = mesh->vertexes[i2];
			// subtracting p0 here is the same as asserting that m is affine!
			// perhaps I should have an Affine Transform class?
			Transform m(p1-p0, p2-p0, -r.d, Point(0,0,0,1));
			Point i = m * (r.o-p0);
			if(i.x >= 0 && i.y >= 0 && i.x+i.y < 1 && i.z > 0.00001f && i.z < t) {
				t = i.z;
				return true;
			}
			return false;
		}

		bool min(Mesh const* mesh, Ray const& r, Dip& dip, float& t) const {
			Point const& p0 = mesh->vertexes[i0];
			Point const& p1 = mesh->vertexes[i1];
			Point const& p2 = mesh->vertexes[i2];
			// subtracting p0 here is the same as asserting that m is affine!
			// perhaps I should have an Affine Transform class?
			Transform m(p1-p0, p2-p0, -r.d, Point(0,0,0,1));
			Point i = m * (r.o-p0);
			if(i.x >= 0 && i.y >= 0 && i.x+i.y < 1 && i.z > 0.00001f && i.z < t) {
				t = i.z;
				float x = i.x, y = i.y;
				Plane p = mesh->vertexes[i0]*mesh->vertexes[i1]*mesh->vertexes[i2];
				Point i = mesh->vertexes[i0]*(1-x-y) + mesh->vertexes[i1]*x + mesh->vertexes[i2]*y;
				dip = (Dip){i, p};
				return true;
			}
			return false;
		}

		Dip r(Mesh const* mesh) const {
			float u = gi_random(), v = gi_random();
			Plane p = mesh->vertexes[i0]*mesh->vertexes[i1]*mesh->vertexes[i2];

			// Transform that maintains stratification
			float x = (v>u) ? 0.5*u : u-0.5*v;
			float y = (v>u) ? v-0.5*u : 0.5*v;

			//float x = (u+v > 1) ? 1-u : u;
			//float y = (u+v > 1) ? 1-v : v;

			Point i = mesh->vertexes[i0]*(1-x-y) + mesh->vertexes[i1]*x + mesh->vertexes[i2]*y;
			Dip dip = {i, p};
			return dip;
		} 
	};


	Mesh(std::vector<Point> points, std::vector<Triangle> triangles)
        : vertexes(points)
        , triangles(triangles)
        , rays(0)
        , tritests(0)
        , bbtests(0) {
    
        float t = 0;
		std::vector<Box> b;
		std::vector<Point> c;
		b.reserve(triangles.size());
		c.reserve(triangles.size());
		for(auto i = triangles.begin(); i != triangles.end(); ++i) {
		
            Point const& v0 = vertexes[i->i0];
		    Point const& v1 = vertexes[i->i1];
		    Point const& v2 = vertexes[i->i2];

		    Plane p = v0*v1*v2;		// plane proportional to area.
		    float a = sqrt(p.a*p.a + p.b*p.b + p.c*p.c);
		
		    t += a;
		    ecdf.push_back(t);
		
			Box box;
			box.insert(v0);
			box.insert(v1);
			box.insert(v2);
			b.push_back(box);
			c.push_back((box.m + box.n)*0.5);
		}

        printf("Starting BVH\n");
		bvh.construct(b,c);
		printf("Ending BVH\n");
	}

	~Mesh() {
		printf("tri/ray: %f\n", (float)tritests/rays);
		printf("bb/ray:  %f\n", (float)bbtests/rays);
	}
	
	bool min(Ray const& r, Dip& dip, float& t) const {
		rays++;
		Point id(1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z, 0);
		if(bvh.n[0].box().intersect(r.o, id, t)) {
			return minBVH(bvh.n[0], r, id, dip, t);
		}
		return false;
	}

	bool any(Ray const& r, float t) const {
		rays++;
		Point id(1.f/r.d.x, 1.f/r.d.y, 1.f/r.d.z, 0);
		if(bvh.n[0].box().intersect(r.o, id, t)) {
			return anyBVH(bvh.n[0], r, id, t);
		}
		return false;
	}

	Dip r() const {
		float s = gi_random() * ecdf.back();
        size_t i = std::upper_bound(ecdf.begin(), ecdf.end(), s) - ecdf.begin();
		Dip dip = triangles[i].r(this);
		dip.p = dip.p * (ecdf.back() / ((i == 0) ? ecdf[0] : (ecdf[i]-ecdf[i-1])));
		return dip;
	}

    float w() const {
        return ecdf.back();
    }

    Box bb() const {
        return bvh.n[0].box();
    }

	bool minBVH(BVH::Node const& node, Ray const& r, Point const& id, Dip& dip, float& t) const {
		if(node.end > 0) {
			bool any = false;
			for(int i = node.begin; i < node.end; i++) {
				tritests++;
				if(triangles[bvh.o[i]].min(this, r, dip, t)) any = true;
			}
			return any;
		}
		else {
			bbtests+=2;
			bool left = (bvh.n[node.begin].box().intersect(r.o, id, t) && 
					minBVH(bvh.n[node.begin], r, id, dip, t));
			bool right = (bvh.n[node.begin+1].box().intersect(r.o, id, t) && 
					minBVH(bvh.n[node.begin+1], r, id, dip, t));
			return left || right;
		}
	}

	bool anyBVH(BVH::Node const& node, Ray const& r, Point const& id, float t) const {
		if(node.end > 0) {
			for(int i = node.begin; i < node.end; i++) {
				tritests++;
				if(triangles[bvh.o[i]].any(this, r, t)) return true;
			}
			return false;
		}
		else {
			bbtests+=2;
			return 	(bvh.n[node.begin].box().intersect(r.o, id, t) && 
					anyBVH(bvh.n[node.begin], r, id, t)) || 
				(bvh.n[node.begin+1].box().intersect(r.o, id, t) && 
					anyBVH(bvh.n[node.begin+1], r, id, t));
		}
	}

private:

	std::vector<Point> vertexes;
	std::vector<Triangle> triangles;
	BVH bvh;
	
	std::vector<float> ecdf;

	mutable int rays;
	mutable int tritests;
	mutable int bbtests;
};

#endif

