
#include "box.h"

#include <limits>
#include "la.h"
#include "ray.h"

int dim(Box const& vb, Box const& cb) {
    auto t = cb.n-cb.m;
    if (t.x > t.y && t.x > t.z && t.x > 0)
        return 0;
    else if (t.y > t.z && t.y > 0)
        return 1;
    else if (t.z > 0)
        return 2;
    else
        return -1;
}

float area(Box const& b) {
    auto t = b.n-b.m;
    return t.x*t.y + t.x*t.z + t.y*t.z;
}

struct State {
    std::vector<Box> const& b;  // Bounding boxes
    std::vector<P3> c;          // Centers of bounding boxes

    std::vector<BVH::Node> n;  // BVH nodes
    std::vector<int> o;        // Indexes into bounding box list

    State(std::vector<Box> const& b) : b(b) {
        c.reserve(b.size());
        o.reserve(b.size());
        n.reserve(2 * b.size() - 1);
    }
};

void constructImpl(State& state, int nindex, Box const& cb) {
    auto& node = state.n[nindex];

    // check termination conditions

    // too few children to split
    if (node.end - node.begin <= 2) return;

    // all children coincide, so splitting won't help
    int d = dim(node.box, cb);
    if (d < 0) return;


    // bucket children into K equally spaced bins by their bounding box center
    static const int K = 8;

    int ni[K] = {0};
    Box vbi[K];
    Box cbi[K];

    auto cdelta = (K * 0.9999) / (cb.n-cb.m);
    float k1 = (d == 0 ? cdelta.x : (d == 1 ? cdelta.y : cdelta.z));

    for (int i = node.begin; i < node.end; ++i) {
        auto cdelta = k1 * (state.c[state.o[i]] - cb.m);
        int bin = (d == 0 ? cdelta.x : (d == 1 ? cdelta.y : cdelta.z));
        ni[bin]++;
        vbi[bin].insert(state.b[state.o[i]]);
        cbi[bin].insert(state.c[state.o[i]]);
    }

    // pick a split between the K bins using SAH
    float bestscore = Infinity;
    int best = 0;
    Box lvb, rvb, lcb, rcb;
    float lscore, rscore;

    // L-to-R scan
    // compute bounds for everything to the left of each potential
    // split location.
    Box lvi[K];
    Box lci[K];
    int lni[K] = {0};
    for (int i = 1; i < K; ++i) {
        lvi[i] = lvi[i - 1];
        lvi[i].insert(vbi[i - 1]);
        lci[i] = lci[i - 1];
        lci[i].insert(cbi[i - 1]);
        lni[i] = lni[i - 1] + ni[i - 1];
    }

    // R-to-L scan
    // compute bounds for everything to the right of each potential 
    // split location. And find optimal split based off SAH.
    Box rvi, rci;
    int rni = 0;
    for (int i = K - 1; i > 0; --i) {
        rvi.insert(vbi[i]);
        rci.insert(cbi[i]);
        rni += ni[i];

        float score = lni[i] * area(lvi[i]) + rni * area(rvi);
        if (lni[i] > 0 && rni > 0 && score < bestscore) {
            bestscore = score;
            best = i;
            lvb = lvi[i];
            lcb = lci[i];
            rvb = rvi;
            rcb = rci;
            lscore = lni[i] * area(lvi[i]);
            rscore = rni * area(rvi);
        }
    }

    // partition input boxes to left and right children (in place)
    int bindex = node.begin;
    int eindex = node.end;
    while(bindex < eindex) {
        auto cdelta = k1 * (state.c[state.o[bindex]] - cb.m);
        int bin = (d == 0 ? cdelta.x : (d == 1 ? cdelta.y : cdelta.z));
        if (bin < best) {
            bindex++;
        } else {
            std::swap(state.o[bindex], state.o[eindex-1]);
            eindex--;
        }
    }

    auto lchild = BVH::Node{lvb, node.begin, bindex};
    auto rchild = BVH::Node{rvb, eindex, node.end};

    // push back child pair...use SAH cost to order our traversal.
    // perhaps counterintuitively, we put the higher cost child first.
    // I think this is because we're much more likely to get a hit there,
    // which may end our ray traversal early.
    // at one point I did some testing and found this was better than
    // traversing based on ray direction, but I wouldn't trust those results.
    // regardless, it does have the advantage of having the same traversal
    // order for all rays which makes things simpler.
    if (lscore < rscore) {
        std::swap(lchild, rchild);
        std::swap(lcb, rcb);
    }

    auto index = state.n.size();    
    node.begin = index;  // index of first child node
    node.end = -d;       // split dimension (negative = inner node)

    state.n.push_back(lchild);
    state.n.push_back(rchild);
    constructImpl(state, index, lcb);
    constructImpl(state, index + 1, rcb);
}

BVH BVH::construct(std::vector<Box> const& b) {
    auto state = State(b);

    Box vb, cb;
    for (size_t i = 0; i < b.size(); i++) {
        P3 s((b[i].m.x + b[i].n.x) * 0.5, (b[i].m.y + b[i].n.y) * 0.5,
             (b[i].m.z + b[i].n.z) * 0.5);

        state.o.push_back(i);
        state.c.push_back(s);

        vb.insert(b[i]);
        cb.insert(s);
    }

    state.n.push_back({vb, 0, int(b.size())});
    constructImpl(state, 0, cb);

    return BVH{std::move(state.o), std::move(state.n)};
}
