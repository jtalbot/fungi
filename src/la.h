
#pragma once

#include <math.h>  // for sin, cos, INFINITY
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#define Infinity (std::numeric_limits<float>::infinity())

constexpr float minf(const float a, const float b) { return a < b ? a : b; }

constexpr float maxf(const float a, const float b) { return a > b ? a : b; }

struct Scalar;
struct Point;
struct Line;
struct Plane;
struct Pseudo;

struct P2;
struct P3;
struct V3;

struct Scalar {
    float s;

    Scalar() {}
    constexpr Scalar(float s) : s(s) {}
    explicit operator float() const { return s; }

    Scalar operator+(Scalar o) const { return s + o.s; }

    Scalar operator-(Scalar o) const { return s - o.s; }

    Scalar operator*(Scalar o) const { return s * o.s; }

    Scalar operator/(Scalar o) const { return s / o.s; }

    Scalar operator-() const { return -s; }
};

constexpr Scalar operator*(float f, Scalar s) { return f * s.s; }
constexpr Scalar operator*(Scalar s, float f) { return s.s * f; }
constexpr Scalar operator/(float f, Scalar s) { return f / s.s; }
constexpr Scalar operator/(Scalar s, float f) { return s.s / f; }
constexpr Scalar minf(Scalar s, Scalar t) { return minf(s.s, t.s); }
constexpr Scalar maxf(Scalar s, Scalar t) { return maxf(s.s, t.s); }

struct Pseudo {
    float s;

    Pseudo() {}
    constexpr explicit Pseudo(float s) : s(s) {}
    explicit operator float() const { return s; }

    Pseudo operator+(Pseudo o) const { return Pseudo(s + o.s); }

    Pseudo operator-(Pseudo o) const { return Pseudo(s - o.s); }

    Scalar operator*(Pseudo o) const { return Scalar(s * o.s); }

    Scalar operator/(Pseudo o) const { return Scalar(s / o.s); }

    Pseudo operator-() const { return Pseudo(-s); }
};

constexpr Pseudo operator*(float f, Pseudo p) { return Pseudo(f * p.s); }
constexpr Pseudo operator*(Pseudo p, float f) { return Pseudo(p.s * f); }
constexpr Pseudo operator/(float f, Pseudo p) { return Pseudo(f / p.s); }
constexpr Pseudo operator/(Pseudo p, float f) { return Pseudo(p.s / f); }
constexpr Pseudo minf(Pseudo s, Pseudo t) { return Pseudo(minf(s.s, t.s)); }
constexpr Pseudo maxf(Pseudo s, Pseudo t) { return Pseudo(maxf(s.s, t.s)); }

constexpr Pseudo operator*(Pseudo p, Scalar s) { return Pseudo(p.s * s.s); }
constexpr Pseudo operator*(Scalar s, Pseudo p) { return Pseudo(s.s * p.s); }

constexpr Pseudo operator/(Pseudo p, Scalar s) { return Pseudo(p.s / s.s); }
constexpr Pseudo operator/(Scalar s, Pseudo p) { return Pseudo(s.s / p.s); }

struct Point {
    // e0, e1, e2, e3
    float x, y, z, w;

    Point() {}
    constexpr Point(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}

    Point operator+(Point const& o) const {
        return Point(x + o.x, y + o.y, z + o.z, w + o.w);
    }

    Point operator-(Point const& o) const {
        return Point(x - o.x, y - o.y, z - o.z, w - o.w);
    }

    Point operator-() const { return Point(-x, -y, -z, -w); }

    // Plane operator~() const;

    void print() const { printf("c(%f, %f, %f, %f)", x, y, z, w); }
} __attribute__((aligned));

struct Plane {
    // e123, e023, e013, e012
    float a, b, c, d;

    Plane() {}
    constexpr Plane(float a, float b, float c, float d) : a(a), b(b), c(c), d(d) {}

    Plane operator+(Plane const& o) const {
        return Plane(a + o.a, b + o.b, c + o.c, d + o.d);
    }

    Plane operator-(Plane const& o) const {
        return Plane(a - o.a, b - o.b, c - o.c, d - o.d);
    }

    Plane operator-() const { return Plane(-a, -b, -c, -d); }

    // Point operator~() const;

    float q() const { return std::sqrt(a * a + b * b + c * c); }

    void print() const { printf("c(%f, %f, %f, %f)", a, b, c, d); }
} __attribute__((aligned));

/*constexpr Plane Point::operator~() const {
    return Plane(1.f/x, 1.f/y, 1.f/z, 1.f/w);
}

constexpr Point Plane::operator~() const {
    return Point(1.f/a, 1.f/b, 1.f/c, 1.f/d);
}*/

constexpr Point operator*(Point const& p, Scalar s) {
    return Point(p.x * s.s, p.y * s.s, p.z * s.s, p.w * s.s);
}

constexpr Point operator*(Scalar s, Point const& p) { return p * s; }

constexpr Plane operator*(Plane const& p, Scalar s) {
    return Plane(p.a * s.s, p.b * s.s, p.c * s.s, p.d * s.s);
}

constexpr Plane operator*(Scalar s, Plane const& p) { return p * s; }

constexpr Point operator*(Plane const& p, Pseudo q) {
    return Point(p.a * q.s, p.b * q.s, p.c * q.s, p.d * q.s);
}

constexpr Point operator*(Pseudo q, Plane const& p) {
    return Point(q.s * p.a, q.s * p.b, q.s * p.c, q.s * p.d);
}

constexpr Plane operator*(Point const& p, Pseudo q) {
    return Plane(p.x * q.s, p.y * q.s, p.z * q.s, p.w * q.s);
}

constexpr Plane operator*(Pseudo q, Point const& p) {
    return Plane(q.s * p.x, q.s * p.y, q.s * p.z, q.s * p.w);
}

constexpr Pseudo operator*(Plane const& p, Point const& v) {
    return Pseudo(p.a * v.x + p.b * v.y + p.c * v.z + p.d * v.w);
}

constexpr Pseudo operator*(Point const& v, Plane const& p) {
    return Pseudo(p.a * v.x + p.b * v.y + p.c * v.z + p.d * v.w);
}

// Line*Pseudo = Line

// Point*Point = Scalar + Line + 0
// Plane*Plane = Scalar + Line + 0
// Point*Plane =      0 + Line + Pseudo
// Plane*Point =      0 + Line + Pseudo

// Line*Line   = Scalar + Line + Pseudo

// Line*Point  =     Point + Plane
// Line*Plane  =     Point + Plane

struct Line {
    // e23, e31, e30
    // e02, e12, e01
    float i, j, k;  // The vector part of plucker coordinates
    float l, m, n;  // The plane part of plucker coordinates
    constexpr Line(float i, float j, float k, float l, float m, float n)
        : i(i), j(j), k(k), l(l), m(m), n(n) {}
};

constexpr Line operator*(Line const& l, Scalar s) {
    return Line(l.i * s.s, l.j * s.s, l.k * s.s, l.l * s.s, l.m * s.s,
                l.n * s.s);
}

constexpr Line operator*(Scalar s, Line const& l) { return l * s; }

constexpr Line operator*(Point const& j, Point const& k) {
    return Line(j.z * k.w - j.w * k.z, j.w * k.y - j.y * k.w,
                j.w * k.x - j.x * k.w, j.x * k.z - j.z * k.x,
                j.y * k.z - j.z * k.y, j.x * k.y - j.y * k.x);
}

constexpr Plane operator*(Line const& b, Point const& a) {
    return Plane(
        a.y * b.i + a.z * b.j + a.w * b.m, -a.z * b.k - a.w * b.l - a.x * b.i,
        a.w * b.n - a.x * b.j + a.y * b.k, -a.x * b.m + a.y * b.l - a.z * b.n);
}

constexpr Plane operator*(Point const& a, Line const& b) { return b * a; }

constexpr Point operator*(Line const& b, Plane const& a) {
    return Point(
        a.b * b.i + a.c * b.j + a.d * b.m, -a.c * b.k - a.d * b.l - a.a * b.i,
        a.d * b.n - a.a * b.j + a.b * b.k, -a.a * b.m + a.b * b.l - a.c * b.n);
}

constexpr Point operator*(Plane const& a, Line const& b) { return b * a; }

// Transforms
struct Transform {
    Pseudo det;            // the determinant
    Scalar idet;           // the determinant of the inverse
    Plane m0, m1, m2, m3;  // the matrix (stored by rows)
    Point i0, i1, i2, i3;  // the inverse (stored by columns)

    Transform() {}

    constexpr Transform(Point const& i0, Point const& i1, Point const& i2,
              Point const& i3)
        : det(i1 * i2 * i3 * i0),
          idet(1.f / det.s)
          // this is the matrix inverse, faster to negate the determinant than
          // to reorder the points (compiler can find more common
          // subexpressions)
          ,
          m0(i1 * i2 * i3 * idet),
          m1(i0 * i2 * i3 * -idet),
          m2(i0 * i1 * i3 * idet),
          m3(i0 * i1 * i2 * -idet),
          i0(i0),
          i1(i1),
          i2(i2),
          i3(i3) {}

    // some standard transforms
    static Transform Identity() {
        return {Point(1, 0, 0, 0), Point(0, 1, 0, 0), Point(0, 0, 1, 0),
                Point(0, 0, 0, 1)};
    }

    static Transform Scale(float sx, float sy, float sz, float sw = 1) {
        return {Point(sx, 0, 0, 0), Point(0, sy, 0, 0), Point(0, 0, sz, 0),
                Point(0, 0, 0, sw)};
    }

    static Transform Translate(float tx, float ty, float tz) {
        return {Point(1, 0, 0, 0), Point(0, 1, 0, 0), Point(0, 0, 1, 0),
                Point(tx, ty, tz, 1)};
    }

    static Transform RotateX(float t) {
        return {Point(1, 0, 0, 0), Point(0, cos(t), sin(t), 0),
                Point(0, -sin(t), cos(t), 0), Point(0, 0, 0, 1)};
    }

    static Transform RotateY(float t) {
        return {Point(cos(t), 0, -sin(t), 0), Point(0, 1, 0, 0),
                Point(sin(t), 0, cos(t), 0), Point(0, 0, 0, 1)};
    }

    static Transform RotateZ(float t) {
        return {Point(cos(t), sin(t), 0, 0), Point(-sin(t), cos(t), 0, 0),
                Point(0, 0, 1, 0), Point(0, 0, 0, 1)};
    }

    static Transform Rotate(V3 const& axis, float t);
    static Transform LookAt(P3 const& eye, P3 const& at, V3 const& up);

    // matrix inverse, swap m and i, reordering entries
    Transform operator~() const {
        Transform t;

        t.idet = det.s;
        t.det = Pseudo(idet.s);

        t.i0 = Point(m0.a, m1.a, m2.a, m3.a);
        t.i1 = Point(m0.b, m1.b, m2.b, m3.b);
        t.i2 = Point(m0.c, m1.c, m2.c, m3.c);
        t.i3 = Point(m0.d, m1.d, m2.d, m3.d);

        t.m0 = Plane(i0.x, i1.x, i2.x, i3.x);
        t.m1 = Plane(i0.y, i1.y, i2.y, i3.y);
        t.m2 = Plane(i0.z, i1.z, i2.z, i3.z);
        t.m3 = Plane(i0.w, i1.w, i2.w, i3.w);

        return t;
    }

    void print() const {
        printf("m <- rbind(\n");
        m0.print();
        printf(",\n");
        m1.print();
        printf(",\n");
        m2.print();
        printf(",\n");
        m3.print();
        printf(")\n");
        printf("i <- cbind(\n");
        i0.print();
        printf(",\n");
        i1.print();
        printf(",\n");
        i2.print();
        printf(",\n");
        i3.print();
        printf(")\n");
    }

} __attribute__((aligned));

// matrix * point
constexpr Point operator*(Transform const& t, Point const& o) {
    return Point(float(t.m0 * o), float(t.m1 * o), float(t.m2 * o),
                 float(t.m3 * o));
}

// plane * matrix
constexpr Plane operator*(Plane const& p, Transform const& t) {
    return Plane(float(p * t.i0), float(p * t.i1), float(p * t.i2),
                 float(p * t.i3));
}

inline Transform operator*(Transform const& m, Transform const& n) {
    Transform p = ~n;
    Transform t;
    t.m0 = m.m0 * p;
    t.m1 = m.m1 * p;
    t.m2 = m.m2 * p;
    t.m3 = m.m3 * p;
    t.i0 = p * m.i0;
    t.i1 = p * m.i1;
    t.i2 = p * m.i2;
    t.i3 = p * m.i3;
    t.det = Pseudo(m.det.s * n.det.s);
    t.idet = m.idet * n.idet;
    return t;
}

// Subspaces for performance and storage and ease of use.

struct P2 {
    float u, v;
    P2() {}
    constexpr P2(float u, float v) : u(u), v(v) {}

    P2 operator-() const { return P2(-u, -v); }

    P2 operator+(P2 const& o) const { return P2(u + o.u, v + o.v); }

    P2 operator-(P2 const& o) const { return P2(u - o.u, v - o.v); }

    P2 operator*(Scalar f) const { return P2(u * f.s, v * f.s); }

    P2 operator/(Scalar f) const { return P2(u / f.s, v / f.s); }

    Scalar length() const { return std::sqrt(float(u * u + v * v)); }
};

constexpr P2 interpolate(P2 const& a, P2 const& b, P2 const& c, float u, float v) {
    return P2(a.u * (1 - u - v) + b.u * u + c.u * v,
              a.v * (1 - u - v) + b.v * u + c.v * v);
}

struct P3 {
    float x, y, z;

    P3() {}
    constexpr P3(float x, float y, float z) : x(x), y(y), z(z) {}
    explicit P3(Point const& p);
    operator Point() const;
};

constexpr P3 interpolate(P3 const& a, P3 const& b, P3 const& c, float u, float v) {
    return P3(a.x * (1 - u - v) + b.x * u + c.x * v,
              a.y * (1 - u - v) + b.y * u + c.y * v,
              a.z * (1 - u - v) + b.z * u + c.z * v);
}

constexpr P3 min(P3 const& p0, P3 const& p1) {
    return P3(minf(p0.x, p1.x), minf(p0.y, p1.y), minf(p0.z, p1.z));
}

constexpr P3 max(P3 const& p0, P3 const& p1) {
    return P3(maxf(p0.x, p1.x), maxf(p0.y, p1.y), maxf(p0.z, p1.z));
}

struct V3 {
    float x, y, z;

    V3() {}
    constexpr V3(float x, float y, float z) : x(x), y(y), z(z) {}
    explicit V3(Point const& p);
    operator Point() const;

    V3 operator-() const { return V3(-x, -y, -z); }

    V3 operator+(V3 const& o) const { return V3(x + o.x, y + o.y, z + o.z); }

    V3 operator-(V3 const& o) const { return V3(x - o.x, y - o.y, z - o.z); }

    V3 operator*(Scalar f) const { return V3(x * f.s, y * f.s, z * f.s); }

    V3 operator/(Scalar f) const { return V3(x / f.s, y / f.s, z / f.s); }
};

constexpr V3 operator/(float f, V3 const& v) { return V3(f/v.x, f/v.y, f/v.z); }

constexpr V3 interpolate(V3 const& a, V3 const& b, V3 const& c, float u, float v) {
    return V3(a.x * (1 - u - v) + b.x * u + c.x * v,
              a.y * (1 - u - v) + b.y * u + c.y * v,
              a.z * (1 - u - v) + b.z * u + c.z * v);
}

inline V3 normalize(V3 const& a) {
    auto w = std::sqrt(float(a.x * a.x + a.y * a.y + a.z * a.z));
    return V3(a.x / w, a.y / w, a.z / w);
}

constexpr V3 cross(V3 const& a, V3 const& b) {
    return V3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
              a.x * b.y - a.y * b.x);
}

constexpr float dot(V3 const& a, V3 const& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

constexpr V3 operator-(P3 const& a, P3 const& b) {
    return V3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline P3::P3(Point const& p) : x(p.x / p.w), y(p.y / p.w), z(p.z / p.w) {
    //if (p.w == 0) throw "Can't convert Point with w = 0 to a P3";
}

inline P3::operator Point() const { return Point(x, y, z, 1.f); }

inline V3::V3(Point const& p) : x(p.x), y(p.y), z(p.z) {
    //if (p.w != 0) throw "Can't convert Point with w != 0 to a V3";
}

inline V3::operator Point() const { return Point(x, y, z, 0.f); }

inline Transform Transform::Rotate(V3 const& axis, float t) {
    V3 a = normalize(axis);
    float c = cos(t);
    float s = sin(t);
    float x = float(a.x), y = float(a.y), z = float(a.z);

    return Transform(Point(x * x + (1 - x * x) * c, y * x * (1 - c) + z * s,
                           z * x * (1 - c) - y * s, 0),
                     Point(x * y * (1 - c) - z * s, y * y + (1 - y * y) * c,
                           z * y * (1 - c) + x * s, 0),
                     Point(x * z * (1 - c) + y * s, y * z * (1 - c) - x * s,
                           z * z + (1 - z * z) * c, 0),
                     Point(0, 0, 0, 1));
}

inline Transform Transform::LookAt(P3 const& eye, P3 const& at, V3 const& up) {
    V3 dir = normalize(at - eye);
    V3 left = normalize(cross(up, dir));
    V3 newUp = cross(dir, left);
    return ~Transform(left, newUp, dir, eye);
}
