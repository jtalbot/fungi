
#include <stdio.h>
#include <string>
#include <iostream>
#include <sys/time.h>
#include <stack>

#include <lua.hpp>

#include "la.h"
#include "ray.h"
#include "shape.h"
#include "color.h"
#include "sensor.h"
#include "camera.h"
#include "light.h"

#include "tiny_obj_loader.cc"

int time(lua_State* L)
{
  timeval time_tt;
  gettimeofday(&time_tt, NULL);
  double r = (double)time_tt.tv_sec * 1000 * 1000 + (double)time_tt.tv_usec;
  lua_pushnumber(L, r);
  return 1;
}


std::vector<Shape const*> all;
std::stack<Transform> transforms;
bool islight = false;


int push(lua_State* L) {
    transforms.push(transforms.top());
    return 0;
}

int pop(lua_State* L) {
    transforms.pop();
    return 0;
}

int scale(lua_State* L) {
    int argc = lua_gettop(L);
    if(argc != 3)
        std::cerr << "scale takes 3 arguments" << std::endl;

    double x = lua_tonumber(L, 1);
    double y = lua_tonumber(L, 2);
    double z = lua_tonumber(L, 3);

    transforms.top() = Transform::Scale(x,y,z) * transforms.top();

    return 0;
}

int translate(lua_State* L) {
    int argc = lua_gettop(L);
    if(argc != 3)
        std::cerr << "translate takes 3 arguments" << std::endl;

    double x = lua_tonumber(L, 1);
    double y = lua_tonumber(L, 2);
    double z = lua_tonumber(L, 3);

    transforms.top() = Transform::Translate(x,y,z) * transforms.top();

    return 0;
}

int rotatex(lua_State* L) {
    int argc = lua_gettop(L);
    if(argc != 1)
        std::cerr << "rotatex takes 1 argument" << std::endl;

    double r = lua_tonumber(L, 1);
    transforms.top() = Transform::RotateX(r) * transforms.top();
    return 0;
}

int rotatey(lua_State* L) {
    int argc = lua_gettop(L);
    if(argc != 1)
        std::cerr << "rotatey takes 1 argument" << std::endl;

    double r = lua_tonumber(L, 1);
    transforms.top() = Transform::RotateY(r) * transforms.top();
    return 0;
}

int rotatez(lua_State* L) {
    int argc = lua_gettop(L);
    if(argc != 1)
        std::cerr << "rotatez takes 1 argument" << std::endl;

    double r = lua_tonumber(L, 1);
    transforms.top() = Transform::RotateZ(r) * transforms.top();
    return 0;
}

std::vector<P3> makeVertexes(tinyobj::shape_t const& shape) {
    std::vector<P3> vertexes;
    for(size_t f = 0; f < shape.mesh.positions.size(); f += 3) {
        vertexes.push_back(P3(
            shape.mesh.positions[f+0],
            shape.mesh.positions[f+1],
            shape.mesh.positions[f+2]));
    }
    return vertexes;
}

std::vector<Mesh::Triangle> makeTriangles(tinyobj::shape_t const& shape) {
    std::vector<Mesh::Triangle> triangles;
    for(size_t f = 0; f < shape.mesh.indices.size(); f += 3) {
        triangles.push_back(Mesh::Triangle(
            shape.mesh.indices[f+0],
            shape.mesh.indices[f+1],
            shape.mesh.indices[f+2]));
    }
    return triangles;
}

std::vector<Shape const*> makeShapes(std::vector<tinyobj::shape_t> const& shapes) {
    std::cout << "# of shapes: " << shapes.size() << std::endl;

    std::vector<Shape const*> out;

    for(size_t i = 0; i < shapes.size(); ++i) {
        Mesh* m = new Mesh(makeVertexes(shapes[i]), makeTriangles(shapes[i]));
        m->light = islight;
        out.push_back(m);
    }
	
    return out;
}

int object(lua_State* L) {
    int argc = lua_gettop(L);
    if(argc != 1)
        std::cerr << "obj(filename) -- takes a single argument" << std::endl;

    std::string filename = lua_tostring(L, 1);

    std::vector<tinyobj::shape_t> shapes;
    std::string err = tinyobj::LoadObj(shapes, filename.c_str());

    if(!err.empty()) {
        std::cerr << err << std::endl;
        exit(1);
    }

    all.push_back(new Transformation( 
        transforms.top(), new Group( makeShapes(shapes) ) ));
    
    return 0;
}

int light(lua_State* L) {
    islight = true;
    return 0;
}

int render(lua_State* L) {
    int argc = lua_gettop(L);
    if(argc != 3)
        std::cerr << "render(filename) -- takes three arguments" << std::endl;

    const int64_t samples = (int64_t)lua_tonumber(L, 1);
    std::string filename = lua_tostring(L, 2);
   
    Shape const* shape = new Group(all);
    Pinhole c(Point(0,-30,-200,1), 1);
	const int width = 750, height = 500;
	Sensor s(width,height,2);
	PointLight light(Point(150,150,-150,1));

	for(int64_t i = 0; i < samples; i++) {
		//Dip p0 = c.sample();
		//Dip p1 = mesh.sample();
		//Point p2 = light.sample();

		//Ray ray = {p0.i, (p1.i+-1*p0.i)};
		//Ray ray2 = {p1.i, (p2+-1*p1.i)};	
		//float x,y;
		//c.project(ray, x, y);

		//if(!mesh.hit(ray, 0.999999) && !mesh.hit(ray2, 0.999999)) {
		//	float g = maxf(p0.p*(~ray.d),0)*maxf((p1.p*(-~ray.d)),0)*1/PI*maxf((p1.p*~ray2.d),0);
		//	s.detect(x, 1-y, 500, g);
		//}
		
		Dip p0, p1, p2;
		Point d0;
		c.sampleDir(p0, d0);
		Ray ray = {p0.i, d0};
        float t = Infinity;
		bool hit = shape->min(ray, p1, t);
        
        if(hit) {
            float x,y;
		    c.project(ray, x, y);
        	
		    //Point p2 = light.sample();
            p2 = shape->r();
            Ray ray2 = {p1.i, p2.i-p1.i};	
			// sample circle
			//float u = gi_random()*2*PI; 
			//float v = gi_random();
			//float x = cos(u)*v;
			//float y = sin(y)*v;
			//float z = sqrt(1-v);
			//out = f*Point(0,-1,0,0) + u*Point(1,0,0,0) + v*Point(0,0,1,0);
            if(/*p2.light &&*/ !shape->any(ray2, 0.9999)) {
                float g = maxf(((~p1.p)*(~ray2.d)),0)
                        * maxf(((~p2.p)*-(~ray2.d)),0);
			    s.detect(x,y, 625, 10000*g);
            }
		}
        else {
            float x,y;
		    c.project(ray, x, y);
            s.detect(x,y,500,300);
        }
	
        if(i%100000==0) {
            lua_pushvalue(L, 3);
            lua_pushnumber(L, (lua_Number)i);
            if( lua_pcall(L, 1, 1, 0) != 0 )
                std::cerr << "Error calling status function: " << lua_tostring(L, -1) << std::endl;
            bool cont = lua_toboolean(L, -1);
            lua_pop(L, 1);
            if(!cont)
                break;
        }
		//if(i%1000000==0) printf("%d\n",i);
	}

    {
        lua_pushvalue(L, 3);
        lua_pushnumber(L, (lua_Number)samples);
        if( lua_pcall(L, 1, 1, 0) != 0 )
            std::cerr << "Error calling status function: " << lua_tostring(L, -1) << std::endl;
        lua_pop(L, 1);
    }
	
	s.outputEXR(filename);

    printf("\n");
    delete shape;

    return 0;
}

void report_errors(lua_State* L, int status)
{
    if ( status!=0 ) {
        std::cerr << "-- " << lua_tostring(L, -1) << std::endl;
        lua_pop(L, 1); // remove error message
    }
}

int main(int argc, char** argv) {

    srand(1);

    transforms.push(Transform::Identity());

    for ( int n=1; n<argc; ++n ) {
        const char* file = argv[n];

        lua_State *L = lua_open();

        luaL_openlibs(L);

        lua_register(L, "time", time);
        lua_register(L, "scale", scale);
        lua_register(L, "push", push);
        lua_register(L, "pop", pop);
        lua_register(L, "translate", translate);
        lua_register(L, "rotatex", rotatex);
        lua_register(L, "rotatey", rotatey);
        lua_register(L, "rotatez", rotatez);
        lua_register(L, "obj", object);
        lua_register(L, "light", light);
        lua_register(L, "render", render);

        std::cerr << "Running: " << file << std::endl;

        int s = luaL_loadfile(L, file);

        if ( s==0 ) {
            // execute Lua program
            s = lua_pcall(L, 0, LUA_MULTRET, 0);
        }

        report_errors(L, s);
        lua_close(L);
    }

	return 0;
}
