#pragma once

namespace cfe
{
namespace EoS
{

class ideal_fluid_fn
{
   public:
    __attribute__((always_inline)) inline double
    operator()(const double gamma, const double dens, const double eng) const
    {
        return (gamma - 1) * dens * eng;
    }
};
constexpr inline ideal_fluid_fn ideal_fluid;

class cold_term_fn
{
   public:
    __attribute__((always_inline)) inline double
    operator()(double bulk_sound_sq, double dens, double init_dens) const
    {
        return bulk_sound_sq * (dens - init_dens);
    }
};
constexpr inline cold_term_fn cold_term;

}  // namespace EoS
}  // namespace cfe