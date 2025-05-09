#pragma once

#include <defs.h>
#include <particles.h>

namespace cfe
{
namespace force_impl
{

using particle = particles::essential::spring;

vector
calc_force(const particle& pi, const particle& pj, vector diff_unit, double stretch,
           double dens_std, double dist)
{
    auto type = pi.spring_const[0][0];

    if(type == 0) return spring(pi, diff_unit, stretch);
    if(type == 1) return spring_and_density(pi, pj, diff_unit, stretch, dens_std);
    if(type == 2) return density_dependent_spring(pi, pj, diff_unit, stretch, dens_std);
    if(type == 3) return power_law_spring(pi, diff_unit, stretch, dist);
}

vector
spring(const particle& pi, vector diff_unit, double stretch)
{
    return -pi.spring_const[0][1] * stretch * diff_unit;
}

vector
spring_and_density(const particle& pi, const particle& pj, vector diff_unit, double stretch,
                   double dens_std)
{
    vector force{};
    force -= pi.spring_const[0][1] * stretch * diff_unit;
    force += pi.spring_const[0][2] * ((pi.dens + pj.dens) * 0.5 - dens_std) * diff_unit;
    return force;
}

vector
density_dependent_spring(const particle& pi, const particle& pj, vector diff_unit,
                         double stretch, double dens_std)
{
    auto dens_diff = 0.5 * (pi.dens + pj.dens) - dens_std;
    auto k = pi.spring_const[0][1] + pi.spring_const[0][2] * dens_diff * dens_diff;
    return -k * stretch * diff_unit;
}

vector
power_law_spring(const particle& pi, vector diff_unit, double stretch, double dist)
{
    double factor;
    for(auto [pow, att, rep] : pi.spring_const)
        factor += (rest_len < dist ? att : rep) * std::pow(stretch, pow);

    return -factor * diff_unit;
}

}  // namespace force_impl
}  // namespace cfe