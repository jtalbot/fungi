
#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include <stack>
#include <string>

#include <lua.hpp>

#include "la.h"

int fungiIdentity(lua_State* L) {
    int argc = lua_gettop(L);
    if (argc != 0) std::cerr << "fungiIdentity takes 0 arguments" << std::endl;

    auto t = new Transform(Point(1, 0, 0, 0), Point(0, 1, 0, 0),
                           Point(0, 0, 1, 0), Point(0, 0, 0, 1));
    lua_pushlightuserdata(L, t);

    return 1;
}

int fungiTranslate(lua_State* L) {
    int argc = lua_gettop(L);
    if (argc != 3) std::cerr << "fungiTranslate takes 3 arguments" << std::endl;

    double tx = lua_tonumber(L, 1);
    double ty = lua_tonumber(L, 2);
    double tz = lua_tonumber(L, 3);

    auto t = new Transform(Point(1, 0, 0, 0), Point(0, 1, 0, 0),
                           Point(0, 0, 1, 0), Point(tx, ty, tz, 1));
    lua_pushlightuserdata(L, t);

    return 1;
}

int fungiScale(lua_State* L) {
    int argc = lua_gettop(L);
    if (argc != 3) std::cerr << "fungiScale takes 3 arguments" << std::endl;

    double dx = lua_tonumber(L, 1);
    double dy = lua_tonumber(L, 2);
    double dz = lua_tonumber(L, 3);

    auto t = new Transform(Point(dx, 0, 0, 0), Point(0, dy, 0, 0),
                           Point(0, 0, dz, 0), Point(0, 0, 0, 1));
    lua_pushlightuserdata(L, t);

    return 1;
}

int fungiRotate(lua_State* L) {
    int argc = lua_gettop(L);
    if (argc != 4) std::cerr << "fungiRotate takes 4 arguments" << std::endl;

    double a = lua_tonumber(L, 1);
    V3 axis = normalize(V3((float)lua_tonumber(L, 2), (float)lua_tonumber(L, 3),
                           (float)lua_tonumber(L, 4)));

    float x = (float)axis.x, y = (float)axis.y, z = (float)axis.z;
    float c = cos(a);
    float s = sin(a);

    auto t =
        new Transform(Point(x * x + (1 - x * x) * c, y * x * (1 - c) + z * s,
                            z * x * (1 - c) - y * s, 0),
                      Point(x * y + (1 - c) - z * s, y * y * (1 - y * y) * c,
                            z * y * (1 - c) + x * s, 0),
                      Point(x * z + (1 - c) + y * s, y * z * (1 - c) - x * s,
                            z * z * (1 - z * z) * c, 0),
                      Point(0, 0, 0, 1));
    lua_pushlightuserdata(L, t);

    return 1;
}

int fungiLookAt(lua_State* L) {
    int argc = lua_gettop(L);
    if (argc != 9) std::cerr << "fungiLookAt takes 9 arguments" << std::endl;

    P3 eye = P3((float)lua_tonumber(L, 1), (float)lua_tonumber(L, 2),
                (float)lua_tonumber(L, 3));

    P3 at = P3((float)lua_tonumber(L, 4), (float)lua_tonumber(L, 5),
               (float)lua_tonumber(L, 6));

    V3 up = V3((float)lua_tonumber(L, 7), (float)lua_tonumber(L, 8),
               (float)lua_tonumber(L, 9));

    V3 dir = normalize(at - eye);
    V3 left = normalize(cross(up, dir));
    V3 newUp = cross(dir, left);

    auto t = new Transform(left, newUp, dir, eye);

    *t = ~*t;
    lua_pushlightuserdata(L, t);

    return 1;
}

int fungiTransform(lua_State* L) {
    int argc = lua_gettop(L);
    if (argc != 16)
        std::cerr << "fungiTransform takes 16 arguments" << std::endl;

    double a = lua_tonumber(L, 1);
    double b = lua_tonumber(L, 2);
    double c = lua_tonumber(L, 3);
    double d = lua_tonumber(L, 4);

    double e = lua_tonumber(L, 5);
    double f = lua_tonumber(L, 6);
    double g = lua_tonumber(L, 7);
    double h = lua_tonumber(L, 8);

    double i = lua_tonumber(L, 9);
    double j = lua_tonumber(L, 10);
    double k = lua_tonumber(L, 11);
    double l = lua_tonumber(L, 12);

    double m = lua_tonumber(L, 13);
    double n = lua_tonumber(L, 14);
    double o = lua_tonumber(L, 15);
    double p = lua_tonumber(L, 16);

    auto t = new Transform(Point(a, e, i, m), Point(b, f, j, n),
                           Point(c, g, k, o), Point(d, h, l, p));
    lua_pushlightuserdata(L, t);

    return 1;
}

int fungiMultiply(lua_State* L) {
    int argc = lua_gettop(L);
    if (argc != 2) std::cerr << "fungiMultiply takes 2 arguments" << std::endl;

    auto a = (Transform*)lua_touserdata(L, 1);
    auto b = (Transform*)lua_touserdata(L, 2);

    Transform* t = new Transform();
    *t = (*a) * (*b);
    lua_pushlightuserdata(L, t);

    return 1;
}

void registerAPI(lua_State* L) {
    lua_register(L, "fungiIdentity", fungiIdentity);
    lua_register(L, "fungiTranslate", fungiTranslate);
    lua_register(L, "fungiScale", fungiScale);
    lua_register(L, "fungiRotate", fungiRotate);
    lua_register(L, "fungiLookAt", fungiLookAt);
    lua_register(L, "fungiTransform", fungiTransform);
    lua_register(L, "fungiMultiply", fungiMultiply);
}
