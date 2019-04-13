
#pragma once

struct Pinhole {
	Point l;
	Plane ip;
	float f;

	Pinhole(Point l, float f) : l(l), f(f) {
		ip = Plane(0,0,1,l.z);
	}
	
	Dip sample() const {
		Dip dip = {l, ip};
		return dip;
	};

	void sampleDir(Dip& dip, Point& out) const {
		dip = (Dip){l, ip};
		
		float u = gi_random()*2-1; 
		float v = gi_random()*2-1; 
		out = f*Point(0,0,1,0) + u*Point(1,0,0,0) + v*Point(0,1,0,0);
	};

    void sample(P2 const& sensor, Dip& dip, Point& out) const {
        dip = (Dip){l, ip};
		out = f*Point(0,0,1,0) + sensor.u*Point(1,0,0,0) + sensor.v*Point(0,1,0,0);
    }

	/*void project(Ray const& r, float& x, float&y) const {
		x = r.d.x/r.d.z * f;
		y = r.d.y/r.d.z * f;
	}*/	
};

