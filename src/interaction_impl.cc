#include <interaction_impl.h>

namespace cfe
{

double
spring::sound_speed(const FP& p) const
{
    return std::sqrt(p.spring_const[1] / p.mass) * p.smth;
}
double
spring::time_step(const FP& p) const
{
    return p.smth / sound_speed(p);
}
vector
spring::force(const EP& pi, const EP& pj) const
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
        throw sh::runtime_error{} << "Unknown type: " << type;
    }
    return force;
}

}  // namespace cfe