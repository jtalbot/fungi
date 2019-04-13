
#pragma once

#include "color.h"
#include "la.h"
#include "ray.h"

class Light {
   public:
    virtual ~Light() {}

    // Point sample() const;
    virtual rgba eval(Point const&) const = 0;

    virtual std::pair<rgba, Point> sample(Dip const&) const = 0;
};

class PointLight : public Light {
   public:
    PointLight(Point p) : p(p) {}

    rgba eval(Point const&) const override final { return rgba(1, 1, 1, 1); }

    std::pair<rgba, Point> sample(Dip const&) const override final {
        return std::make_pair(rgba(1, 1, 1, 1), p);
    }

   private:
    Point p;
};

class InfiniteAreaLight : public Light {
   public:
    InfiniteAreaLight(rgba l, std::string const& filename);

    ~InfiniteAreaLight() override final;

    rgba eval(Point const& p) const override final;

    std::pair<rgba, Point> sample(Dip const&) const override final;

   private:
    const rgba l;
    float* data;
    int width;
    int height;

    std::vector<std::vector<float>> cdf_columns;
    std::vector<float> cdf_rows;

    void readExr(std::string const& filename);
    void buildDistribution();
};
