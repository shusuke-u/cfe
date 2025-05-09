#pragma once

#include <cmath>
#include <numbers>

#include <defs.h>

namespace cfe
{
namespace kernel
{

class gauss_fn
{
   private:
    static constexpr auto sqrt_pi = 1.772453850905516;

   public:
    __attribute__((always_inline)) inline constexpr double
    operator()(const vector diff, const double smth, const std::size_t dim) const
    {
        // Negative operator on `unsigned int` is dangerous because it assigns a negative number
        // to the return value of `unsigned int` type.
        return std::pow(smth * sqrt_pi, -(double)dim) *
               std::exp(-diff.square() / (smth * smth));
    }
};
constexpr inline gauss_fn gauss;

class grad_gauss_fn
{
   public:
    __attribute__((always_inline)) inline constexpr vector
    operator()(const vector diff, const double smth, const std::size_t dim) const
    {
        return -2 * diff * gauss(diff, smth, dim) / (smth * smth);
    }
};
constexpr inline grad_gauss_fn grad_gauss;

}  // namespace kernel
}  // namespace cfe