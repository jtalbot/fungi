
#pragma once

#include <memory>

#include "color.h"
#include "la.h"
#include "ray.h"
#include "sample.h"

class Texture {
   public:
    virtual ~Texture() {}
    virtual rgba eval(P2 const& uv) const = 0;
};

class ConstantTexture : public Texture {
   public:
    ConstantTexture(rgba v) : v(v) {}

    rgba eval(P2 const& uv) const override final { return v; }

   private:
    const rgba v;
};

class ScaleTexture : public Texture {
   public:
    ScaleTexture(std::shared_ptr<Texture const> a,
                 std::shared_ptr<Texture const> b)
        : a(std::move(a)), b(std::move(b)) {}

    rgba eval(P2 const& uv) const override final {
        return a->eval(uv) * b->eval(uv);
    }

   private:
    const std::shared_ptr<Texture const> a, b;
};

class MixTexture : public Texture {
   public:
    MixTexture(std::shared_ptr<Texture const> a,
               std::shared_ptr<Texture const> b,
               std::shared_ptr<Texture const> m)
        : a(std::move(a)), b(std::move(b)), m(std::move(m)) {}

    rgba eval(P2 const& uv) const override final {
        auto f = m->eval(uv);
        return (rgba(1, 1, 1, 1) + f * -1) * a->eval(uv) + f * b->eval(uv);
    }

   private:
    const std::shared_ptr<Texture const> a, b, m;
};

class ImageTexture : public Texture {
   public:
    ImageTexture(std::string const& filename);
    rgba eval(P2 const& uv) const override final;

   private:
    class Impl;
    const std::shared_ptr<Impl const> a;
};

class Material {
   public:
    virtual ~Material() {}
    virtual rgba eval(P2 const& uv) const = 0;
    virtual std::pair<rgba, V3> sample(Dip const& dip, V3 const& in) const = 0;
};

class MixMaterial : public Material {
   public:
    MixMaterial(std::shared_ptr<Texture const> amount,
                std::shared_ptr<Material const> a,
                std::shared_ptr<Material const> b)
        : amount(std::move(amount)), a(std::move(a)), b(std::move(b)) {}

    rgba eval(P2 const& uv) const override final {
        rgba f = amount->eval(uv);
        return (rgba(1, 1, 1, 1) + f * -1) * a->eval(uv) + f * b->eval(uv);
    }

    std::pair<rgba, V3> sample(Dip const& dip,
                               V3 const& in) const override final {
        rgba f = amount->eval(dip.uv);
        return (f.r < gi_random()) ? a->sample(dip, in) : b->sample(dip, in);
    }

   private:
    const std::shared_ptr<Texture const> amount;
    const std::shared_ptr<Material const> a, b;
};

class BumpMaterial : public Material {
   public:
    BumpMaterial(std::shared_ptr<Texture const> bumpMap,
                 std::shared_ptr<Material const> m)
        : bumpMap(std::move(bumpMap)), m(std::move(m)) {}

    rgba eval(P2 const& uv) const override final;

    std::pair<rgba, V3> sample(Dip const& dip,
                               V3 const& in) const override final;

   private:
    const std::shared_ptr<Texture const> bumpMap;
    const std::shared_ptr<Material const> m;
};

class DiffuseMaterial : public Material {
   public:
    DiffuseMaterial(std::shared_ptr<Texture const> kd) : kd(std::move(kd)) {}

    rgba eval(P2 const& uv) const override final { return kd->eval(uv); }

    std::pair<rgba, V3> sample(Dip const& dip,
                               V3 const& in) const override final;

   private:
    const std::shared_ptr<Texture const> kd;
};

class GlassMaterial : public Material {
   public:
    GlassMaterial(std::shared_ptr<Texture const> Kr,
                  std::shared_ptr<Texture const> Kt)
        : Kr(std::move(Kr)), Kt(std::move(Kt)), eta(1.5) {}

    rgba eval(P2 const& uv) const override final { return rgba(0, 0, 0, 0); }

    std::pair<rgba, V3> sample(Dip const& dip,
                               V3 const& in) const override final;

   private:
    const std::shared_ptr<Texture const> Kr, Kt;
    const float eta;
};

class MetalMaterial : public Material {
   public:
    MetalMaterial(std::shared_ptr<Texture const> eta,
                  std::shared_ptr<Texture const> k)
        : eta(std::move(eta)), k(std::move(k)) {}

    rgba eval(P2 const& uv) const override final { return rgba(0, 0, 0, 0); }

    std::pair<rgba, V3> sample(Dip const& dip,
                               V3 const& in) const override final;

   private:
    const std::shared_ptr<Texture const> eta, k;
};
