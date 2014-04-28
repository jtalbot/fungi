
#ifndef _BOX_H
#define _BOX_H

#include <limits>
#include "la.h"

#define Infinity (std::numeric_limits<float>::infinity())

inline float minf(const float a, const float b) {
	return a < b ? a : b;
}

inline float maxf(const float a, const float b) {
	return a > b ? a : b;
}

inline P3 min(P3 const& p0, P3 const& p1) {
	return P3(minf(p0.x, p1.x), minf(p0.y, p1.y), minf(p0.z, p1.z));
}

inline P3 max(P3 const& p0, P3 const& p1) {
	return P3(maxf(p0.x, p1.x), maxf(p0.y, p1.y), maxf(p0.z, p1.z));
}

struct Box {

	P3 m,n;

	Box()
        : m( Infinity,  Infinity,  Infinity)
        , n(-Infinity, -Infinity, -Infinity)
    {}

	bool intersect(Point const& o, Point const& id, const float t) const {
		float lmin=0, lmax=t;
	
		const float xa = (m.x-o.x)*id.x;
		const float xb = (n.x-o.x)*id.x;
		lmin = maxf(lmin, minf(xa,xb));
		lmax = minf(lmax, maxf(xa,xb));

		const float ya = (m.y-o.y)*id.y;
		const float yb = (n.y-o.y)*id.y;
		lmin = maxf(lmin, minf(ya,yb));
		lmax = minf(lmax, maxf(ya,yb));

		const float za = (m.z-o.z)*id.z;
		const float zb = (n.z-o.z)*id.z;
		lmin = maxf(lmin, minf(za,zb));
		lmax = minf(lmax, maxf(za,zb));

		return lmin <= lmax;
	}

	Box& insert(P3 const& p) {
		m = min(m, p);
		n = max(n, p);
        return *this;
	}

	Box& insert(Box const& b) {
		m = min(m, b.m);
		n = max(n, b.n);
        return *this;
	}
} __attribute__ ((aligned));

inline Box operator*(Transform const& t, Box const& b) {
    return Box()
        .insert(P3(t * b.m))
        .insert(P3(t * P3(b.n.x,b.m.y,b.m.z)))
        .insert(P3(t * P3(b.m.x,b.n.y,b.m.z)))
        .insert(P3(t * P3(b.m.x,b.m.y,b.n.z)))
        .insert(P3(t * P3(b.n.x,b.n.y,b.m.z)))
        .insert(P3(t * P3(b.m.x,b.n.y,b.n.z)))
        .insert(P3(t * P3(b.n.x,b.m.y,b.n.z)))
        .insert(P3(t * b.n));
}


struct BVH {
	
	static const int K = 8;

	struct Node {
        Box box;
		int begin, end;
		Node(Box const& box, int begin, int end)
            : box(box), begin(begin), end(end) {}
	} __attribute__ ((aligned));

	std::vector<int> o;
	std::vector<Node> n;

	struct State {
		std::vector<Box> const& b;
		std::vector<P3> c;
		std::vector<int> temp;
		State(std::vector<Box> const& b, std::vector<P3> const& c) : b(b), c(c), temp(b.size()) {}
	};

    Box const& root() const {
        return n[0].box;
    }

	void construct(std::vector<Box> const& b) {
		o.reserve(b.size());
		n.reserve(2*b.size()-1);

        std::vector<P3> c;
		Box vb, cb;
		for(size_t i = 0; i < b.size(); i++) {
            P3 s( (b[i].m.x+b[i].n.x)*0.5,
                  (b[i].m.y+b[i].n.y)*0.5,
                  (b[i].m.z+b[i].n.z)*0.5 );

            o.push_back(i);
            c.push_back(s);

            vb.insert(b[i]);
            cb.insert(s);
        }

		State state(b, c);
		
		n.push_back(Node(vb, 0, o.size()));
		construct(state, 0, vb, cb, 0, o.size());
	}

	void construct(State& state, int nindex, Box const& vb, Box const& cb, int begin, int end) {
		int d;	
		// check termination conditions
		if(end-begin <= 2 || (d = dim(vb, cb)) < 0) {
			return;
		}

		int ni[K] = {0};
		Box vbi[K];
		Box cbi[K];

		float k1 = (K*0.9999) /
                (d == 0 ? cb.n.x-cb.m.x :
				(d == 1 ? cb.n.y-cb.m.y :
					      cb.n.z-cb.m.z));
			
		for(int i = begin; i < end; ++i) {
			int bin = k1 *
                  (d == 0 ? state.c[o[i]].x-cb.m.x :
				  (d == 1 ? state.c[o[i]].y-cb.m.y :
					        state.c[o[i]].z-cb.m.z));
			ni[bin]++;
			vbi[bin].insert(state.b[o[i]]);
			cbi[bin].insert(state.c[o[i]]);
		}

		float bestscore = Infinity; int best = 0;
		Box lvb, rvb, lcb, rcb;

		// L-to-R scan
		Box lvi[K];
		Box lci[K];
		int lni[K] = {0};
		for(int i = 1; i < K; ++i) {
			lvi[i] = lvi[i-1]; lvi[i].insert(vbi[i-1]);
			lci[i] = lci[i-1]; lci[i].insert(cbi[i-1]);
			lni[i] = lni[i-1] + ni[i-1];
		}

		// R-to-L scan
		Box rvi, rci;
		int rni = 0;
		for(int i = K-1; i > 0; --i) {
			rvi.insert(vbi[i]);
			rci.insert(cbi[i]);
			rni += ni[i];
			
			float score = lni[i] * area(lvi[i]) + rni * area(rvi);
			if(lni[i] > 0 && rni > 0 && score < bestscore) { 
				bestscore = score; 
				best = i;
				lvb = lvi[i];
				lcb = lci[i];
				rvb = rvi;
				rcb = rci;
			}
		}

		// partition left and right
		int bindex = begin;
		int eindex = end-1;
		for(int i = begin; i < end; ++i) {
			int bin = d == 0 ? (state.c[o[i]].x-cb.m.x)*k1 :
				  (d == 1 ? (state.c[o[i]].y-cb.m.y)*k1 :
					   (state.c[o[i]].z-cb.m.z)*k1) ;
			if(bin < best) 	{ state.temp[bindex++] = o[i]; }
			else 		{ state.temp[eindex--] = o[i]; }
		}
		for(int i = begin; i < end; ++i) {
			o[i] = state.temp[i];
		}

		// push back child pair...area seems to be a good ordering (better than ray direction?)
		if(area(lvb) < area(rvb)) {
			n.push_back(Node(lvb, begin, bindex));
			n.push_back(Node(rvb, eindex+1, end));
			n[nindex].begin = n.size()-2;
			n[nindex].end = -d;

			construct(state, n[nindex].begin, lvb, lcb, begin, bindex);
			construct(state, n[nindex].begin+1, rvb, rcb, eindex+1, end);
		} else {
			n.push_back(Node(rvb, eindex+1, end));
			n.push_back(Node(lvb, begin, bindex));
			n[nindex].begin = n.size()-2;
			n[nindex].end = -d;

			construct(state, n[nindex].begin, rvb, rcb, eindex+1, end);
			construct(state, n[nindex].begin+1, lvb, lcb, begin, bindex);
		}
	}

	int dim(Box const& vb, Box const& cb) {
		float x = cb.n.x-cb.m.x;
		float y = cb.n.y-cb.m.y;
		float z = cb.n.z-cb.m.z;
		if(x > y && x > z && x > 0) return 0;
		else if(y > z && y > 0) return 1;
		else if(z > 0) return 2;
		else return -1;
	}

	float area(Box const& b) {
		return 	(b.n.x-b.m.x)*(b.n.y-b.m.y) +
			(b.n.x-b.m.x)*(b.n.z-b.m.z) +
			(b.n.z-b.m.z)*(b.n.y-b.m.y);
	}
};

#endif
