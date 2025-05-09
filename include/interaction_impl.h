#pragma once

#include <particles.h>

namespace cfe
{

template <class Tfp, class Tep, class... Args>
class basic_interaction
{
   public:
    using FP = Tfp;
    using EP = Tep;

    virtual double sound_speed(const FP& p, Args... args) const = 0;
    virtual double time_step(const FP& p, Args... args) const = 0;
    virtual vector force(const Tep& pi, const Tep& pj, Args... args) const = 0;
};

template <class... Args>
using spring_interaction =
    basic_interaction<particles::full_ptcl, particles::essential::spring, Args...>;

class spring : public spring_interaction<>
{
   public:
    double
    sound_speed(const FP& p) const override
    {
        return std::sqrt(p.spring_const[1] / p.mass) * p.smth;
    }

    double
    time_step(const FP& p) const override
    {
        return p.smth / sound_speed(p);
    }

    vector
    force(const EP& pi, const EP& pj) const override
    {
        auto diff = pi.pos - pj.pos;
        auto dist = diff.length();
        if(dist == 0) return {};
        auto diff_init = pi.init_pos - pj.init_pos;
        auto rest_len = diff_init.length();
        if(pi.r_search < rest_len) return {};
        auto stretch = dist - rest_len;
        auto type = pi.spring_const[0];
        vector force;
        if(type == 0)
        {
            auto factor = pi.spring_const[1] * stretch / dist;
            force = factor * diff;
        }
        else if(type == 1)
        {
            auto factor = pi.spring_const[1] * stretch / dist + pi.spring_const[2] * stretch;
            force = factor * diff;
        }
        else if(type == 2)
        {
            auto factor =
                pi.spring_const[1] * stretch / dist + pi.spring_const[2] * stretch * stretch;
            force = factor * diff;
        }
        else
        {
            throw std::runtime_error("Unknown spring type: " + std::to_string(type));
        }
        return force;
    }
};

class density_dependent_spring_with_init : public spring_interaction<double>
{
   public:
    double
    sound_speed(const FP& p, double dens_diff_std) const override
    {
        auto dens_diff = p.dens - p.std_dens;
        auto k = p.spring_const[1] + p.spring_const[2] * std::pow(dens_diff / dens_diff_std, 2);
        return std::sqrt(k / p.mass) * p.smth;
    }

    double
    time_step(const FP& p, double dens_diff_std) const override
    {
        return p.smth / sound_speed(p, dens_diff_std);
    }
};

}  // namespace cfe