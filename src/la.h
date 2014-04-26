
#ifndef _LA_H_
#define _LA_H_

#include <math.h>		// for sin, cos, INFINITY
#include <stdlib.h>

struct Point {
	float x, y, z, w;

	Point() {}
	Point(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}

	Point operator+(Point const& o) const {
		return Point(x+o.x, y+o.y, z+o.z, w+o.w);
	}

	// gcc can't turn a+-1*b into a sub even with -ffast-math.
	Point operator-(Point const& o) const {
		return Point(x-o.x, y-o.y, z-o.z, w-o.w);
	}

	Point operator-() const {
		return Point(-x, -y, -z, -w);
	}

	Point operator~() const {
		float s = 1./(x*x+y*y+z*z);
		return Point(x*s, y*s, z*s, w*s);
	}
	
	void print() const { printf("c(%f, %f, %f, %f)", x, y, z, w); }
} __attribute__ ((aligned));

struct Line {
	float i, j, k, l, m, n;
	Line(float i, float j, float k, float l, float m, float n) :
		i(i), j(j), k(k), l(l), m(m), n(n) {}
};


struct Plane {
	float a, b, c, d;

	Plane() {}
	Plane(float a, float b, float c, float d) : a(a), b(b), c(c), d(d) {}

	Plane operator+(Plane const& o) const {
		return Plane(a+o.a, b+o.b, c+o.c, d+o.d);
	}

	Plane operator~() const {
		float s = a*a+b*b+c*c;
		return Plane(a/s, b/s, c/s, d/s);
	}

	void print() const { printf("c(%f, %f, %f, %f)", a, b, c, d); }
} __attribute__ ((aligned));


inline Point operator*(Point const& p, float s) { return Point(p.x*s, p.y*s, p.z*s, p.w*s); }
inline Point operator*(float s, Point const& p) { return p*s; }

inline Plane operator*(Plane const& p, float s) { return Plane(p.a*s, p.b*s, p.c*s, p.d*s); }
inline Plane operator*(float s, Plane const& p) { return p*s; }

inline float operator*(Plane const& p, Point const& v) { return p.a*v.x + p.b*v.y + p.c*v.z + p.d*v.w; }

inline Line operator*(Point const& j, Point const& k) {
	return Line(	
		j.z*k.w-j.w*k.z, 
		j.w*k.y-j.y*k.w,
		j.w*k.x-j.x*k.w,
		j.x*k.z-j.z*k.x,
		j.y*k.z-j.z*k.y,
		j.x*k.y-j.y*k.x);
}

inline Plane operator*(Line const& b, Point const& a) {
	return Plane( 	 
		 a.y*b.i + a.z*b.j + a.w*b.m,
		-a.z*b.k - a.w*b.l - a.x*b.i,
		 a.w*b.n - a.x*b.j + a.y*b.k,
		-a.x*b.m + a.y*b.l - a.z*b.n );
}
inline Plane operator*(Point const& a, Line const& b) { return b*a; }


// Transform from global space to local space 
struct Transform {
	float det;			// the determinant
	Plane m0, m1, m2, m3;		// the transformation matrix (stored by rows)
	Point i0, i1, i2, i3;		// its inverse (stored by columns)

	Transform() {}

	Transform(Point const& i0, Point const& i1, Point const& i2, Point const& i3) :
		det(1.0f / (i1*i2*i3*i0)),
		// this is the matrix inverse...faster to negate the determinant 
		// than to reorder the points (compiler can find more common subexpressions)
		m0(i1*i2*i3*det), 
		m1(i0*i2*i3*-det), 
		m2(i0*i1*i3*det), 
		m3(i0*i1*i2*-det), 
		i0(i0), i1(i1), i2(i2), i3(i3) {}

	// some standard transforms (from global to local, the inverse of what you're used to!!) 
	// Reasoning: eliminates a transform when storing.
	static Transform Identity() {
		return Transform(Point(1,0,0,0), Point(0,1,0,0), Point(0,0,1,0), Point(0,0,0,1));
	}

	static Transform Scale(float sx, float sy, float sz, float sw=1) {
		return Transform(Point(sx,0,0,0), Point(0,sy,0,0), Point(0,0,sz,0), Point(0,0,0,sw));
	}

	static Transform Translate(float tx, float ty, float tz) {
		return Transform(Point(1,0,0,0), Point(0,1,0,0), Point(0,0,1,0), Point(tx,ty,tz,1));
	}

	static Transform  RotateX(float t) {
		return Transform(Point(1,0,0,0), Point(0,cos(t),sin(t),0), Point(0,-sin(t),cos(t),0), Point(0,0,0,1));
	}
	
	static Transform  RotateY(float t) {
		return Transform(Point(cos(t),0,-sin(t),0), Point(0,1,0,0), Point(sin(t),0,cos(t),0), Point(0,0,0,1));
	}
	
	static Transform  RotateZ(float t) {
		return Transform(Point(cos(t),sin(t),0,0), Point(-sin(t),cos(t),0,0), Point(0,0,1,0), Point(0,0,0,1));
	}

	// matrix inverse, swap m and i, reordering entries 
	Transform operator~() const {
		Transform t;
		t.det = 1.f/det;
		t.m0 = Plane(i0.x, i1.x, i2.x, i3.x);
		t.m1 = Plane(t.i0.y, t.i1.y, t.i2.y, t.i3.y);
		t.m2 = Plane(t.i0.z, t.i1.z, t.i2.z, t.i3.z);
		t.m3 = Plane(t.i0.w, t.i1.w, t.i2.w, t.i3.w);
		t.i0 = Point(t.m0.a, t.m1.a, t.m2.a, t.m3.a);
		t.i1 = Point(t.m0.b, t.m1.b, t.m2.b, t.m3.b);
		t.i2 = Point(t.m0.c, t.m1.c, t.m2.c, t.m3.c);
		t.i3 = Point(t.m0.d, t.m1.d, t.m2.d, t.m3.d);
		return t;
	}

	void print() const {
		printf("m <- rbind(\n");
		m0.print(); printf(",\n");
		m1.print(); printf(",\n");
		m2.print(); printf(",\n");
		m3.print(); printf(")\n");
		printf("i <- cbind(\n");
		i0.print(); printf(",\n");
		i1.print(); printf(",\n");
		i2.print(); printf(",\n");
		i3.print(); printf(")\n");
	}
} __attribute__ ((aligned));

inline Point operator*(Transform const& t, Point const& o) {		// matrix * point
	return Point(t.m0*o,t.m1*o,t.m2*o,t.m3*o);
}
	
inline Plane operator*(Plane const& p, Transform const& t) { 		// plane * matrix
	return Plane(p*t.i0,p*t.i1,p*t.i2,p*t.i3);
}

inline Transform operator*(Transform const& m, Transform const& n) {
	Transform p = ~n;
	Transform t;
	t.m0 = m.m0*p; t.m1 = m.m1*p; t.m2 = m.m2*p; t.m3 = m.m3*p;
	t.i0 = p*m.i0; t.i1 = p*m.i1; t.i2 = p*m.i2; t.i3 = p*m.i3;
	t.det = m.det*n.det;
	return t;
}

#endif

