
#include <cfloat>

#include "catch.hpp"

#include "../src/la.h"

bool AlmostEqualRelative(float A, float B,
                         float maxRelDiff = FLT_EPSILON)
{
    // Calculate the difference.
    float diff = fabs(A - B);
    A = fabs(A);
    B = fabs(B);
    // Find the largest
    float largest = (B > A) ? B : A;
 
    if (diff <= largest * maxRelDiff)
        return true;
    return false;
}

TEST_CASE( "Transform Construction", "[Transform]" ) {

    auto t = Transform(
        Point(1,2,3,4),
        Point(4,3,2,1),
        Point(2,4,1,3),
        Point(3,2,4,1)
        );

    REQUIRE(float(t.det) == -50.f);
    REQUIRE(AlmostEqualRelative(float(t.idet), -1.f/50.f));

    REQUIRE(AlmostEqualRelative(float(t.m0.a), 3.f/10.f));
    REQUIRE(AlmostEqualRelative(float(t.m0.d), 5.f/10.f));
    REQUIRE(AlmostEqualRelative(float(t.m3.a), -5.f/10.f));
    REQUIRE(AlmostEqualRelative(float(t.m3.d), -5.f/10.f));
    
    REQUIRE(float(t.i0.x) == 1.f);
    REQUIRE(float(t.i1.y) == 3.f);
    REQUIRE(float(t.i2.z) == 1.f);
    REQUIRE(float(t.i3.w) == 1.f);
    REQUIRE(float(t.i0.w) == 4.f);
    REQUIRE(float(t.i3.x) == 3.f);
}


TEST_CASE( "Transform Inverse", "[Transform]" ) {

    auto t = Transform(
        Point(1,2,3,4),
        Point(4,3,2,1),
        Point(2,4,1,3),
        Point(3,2,4,1)
        );

    t = ~t;

    REQUIRE(AlmostEqualRelative(float(t.det), -1.f/50.f));
    REQUIRE(AlmostEqualRelative(float(t.idet), -50));
 
    REQUIRE(AlmostEqualRelative(float(t.m0.a), 1.f));
    REQUIRE(AlmostEqualRelative(float(t.m0.d), 3.f));
    REQUIRE(AlmostEqualRelative(float(t.m3.a), 4.f));
    REQUIRE(AlmostEqualRelative(float(t.m3.d), 1.f));
    
    REQUIRE(AlmostEqualRelative(float(t.i0.x), 3.f/10.f));
    REQUIRE(AlmostEqualRelative(float(t.i3.x), 5.f/10.f));
    REQUIRE(AlmostEqualRelative(float(t.i0.w), -5.f/10.f));
    REQUIRE(AlmostEqualRelative(float(t.i3.w), -5.f/10.f));
}


TEST_CASE( "Transform Multiplication", "[Transform]" ) {

    auto s = Transform::Translate(3,4,5);
    auto t = Transform::Scale(2,3,4,5);

    t = t*s;

    REQUIRE(AlmostEqualRelative(float(t.det), 120.f));
    REQUIRE(AlmostEqualRelative(float(t.idet), 1.f/120.f));
    
    REQUIRE(AlmostEqualRelative(float(t.m0.a), 1.f/2.f));
    REQUIRE(AlmostEqualRelative(float(t.m1.b), 1.f/3.f));
    REQUIRE(AlmostEqualRelative(float(t.m2.c), 1.f/4.f));
    REQUIRE(AlmostEqualRelative(float(t.m0.d), -3.f/2.f));
    REQUIRE(AlmostEqualRelative(float(t.m1.d), -4.f/3.f));
    REQUIRE(AlmostEqualRelative(float(t.m2.d), -5.f/4.f));
    REQUIRE(AlmostEqualRelative(float(t.m3.d), 1.f/5.f));
    
    REQUIRE(AlmostEqualRelative(float(t.i0.x), 2.f));
    REQUIRE(AlmostEqualRelative(float(t.i1.y), 3.f));
    REQUIRE(AlmostEqualRelative(float(t.i2.z), 4.f));
    REQUIRE(AlmostEqualRelative(float(t.i3.x), 15.f));
    REQUIRE(AlmostEqualRelative(float(t.i3.y), 20.f));
    REQUIRE(AlmostEqualRelative(float(t.i3.z), 25.f));
    REQUIRE(AlmostEqualRelative(float(t.i3.w), 5.f));
}

