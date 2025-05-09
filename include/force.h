#pragma once

#include <sh/secured.h>
#include <sh/stack_map.h>

#include <defs.h>

namespace cfe
{

using spring_constant = std::array<double, 3>;

struct spring_info
{
    spring_constant constant;
    double length;
};

namespace interaction
{
namespace result
{

// class make_spring
// {
//    public:
//     sh::stack_map<std::size_t, spring_info, n_max_neighbors> springs;

//     void
//     clear()
//     {
//         springs.clear();
//     }
// };

class spring
{
   public:
    vector force;
    double eng_dot;

    std::size_t n_ngb;

    void
    clear()
    {
        force.clear();
        eng_dot = 0;
        n_ngb = 0;
    }
};

class repulsion
{
   public:
    vector force;

    std::size_t n_ngb;

    void
    clear()
    {
        force.clear();
        n_ngb = 0;
    }
};

class density
{
   public:
    double dens;

    // The number of neighbors interacted.
    std::size_t n_ngb;

    void
    clear()
    {
        dens = 0;
        n_ngb = 0;
    }
};

}  // namespace result
}  // namespace interaction

}  // namespace cfe