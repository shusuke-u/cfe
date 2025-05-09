#pragma once

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <sh/stack_map.h>

#include <force.h>

namespace cfe
{

namespace particles
{

template <class Tp>
concept particle =
    std::is_copy_assignable_v<Tp> and std::is_copy_constructible_v<Tp> and requires(Tp p) {
        { p.mass } -> std::same_as<double>;

        { p.pos } -> std::same_as<vector>;
        { p.vel } -> std::same_as<vector>;
        { p.mom } -> std::same_as<vector>;

        { p.belongs_to } -> std::same_as<std::array<bool, n_max_views>>;

        { p.make_belong_to(std::declval<std::size_t>()) };
        { p.make_leave_from(std::declval<std::size_t>()) };
    };

class full_ptcl
{
   public:
    // bool is_ghost;

    vector_t<bool> is_fixed;

    std::array<bool, n_max_views> belongs_to;

    std::size_t id;
    std::size_t owner_id;
    std::optional<std::size_t> unit_id;

    vector pos;
    vector vel;
    vector force;

    vector init_pos;
    vector init_vel;

    double std_dens;

    double mass;
    double pres;
    double dens;
    double eng;
    double eng_dot;
    double r_search_for_spring;
    double r_search_for_dens;
    double smth;

    double dissipate;

    double gamma = 5.0 / 3;

    // sh::stack_map<std::size_t, spring_info, n_max_neighbors> springs;

    // spring_constant spring_const;
    double spring_const;
    double bulk_sound_sq;

    std::size_t n_ngb_spr;
    std::size_t n_ngb_dens;
    std::size_t n_ngb_rep;

    full_ptcl();

    __attribute__((always_inline)) inline vector acc() const;
    __attribute__((always_inline)) inline vector mom() const;
    __attribute__((always_inline)) inline vector angular_vel(vector center = vector(0)) const;
    __attribute__((always_inline)) inline vector angular_mom(vector center = vector(0)) const;
    __attribute__((always_inline)) inline vector pos_diff() const;
    __attribute__((always_inline)) inline double kinetic_eng() const;
    __attribute__((always_inline)) inline double total_eng() const;

    double sound_speed() const;
    double time_step() const;

    struct header_t
    {
        friend std::ostream& operator<<(std::ostream& os, header_t);
    };
    static constexpr header_t header{};
    static const std::vector<std::string> header_v;

    friend std::ostream& operator<<(std::ostream& os, const full_ptcl& src);
    friend std::istream& operator>>(std::istream& is, full_ptcl& tar);

    // void copyFromForce(const interaction::result::make_spring& spr);
    void copyFromForce(const interaction::result::spring& spr);
    void copyFromForce(const interaction::result::repulsion& rep);
    void copyFromForce(const interaction::result::density& dens);

    PS::F64vec getPos() const;
    void setPos(const PS::F64vec& src);
};

inline vector
full_ptcl::acc() const
{
    return force / mass;
}
inline vector
full_ptcl::mom() const
{
    return mass * vel;
}
inline vector
full_ptcl::pos_diff() const
{
    return pos - init_pos;
}
inline vector
full_ptcl::angular_vel(vector center) const
{
    return (pos - center) ^ vel;
}
inline vector
full_ptcl::angular_mom(vector center) const
{
    return (pos - center) ^ mom();
}
inline double
full_ptcl::kinetic_eng() const
{
    return 0.5 * mass * vel.square();
}
inline double
full_ptcl::total_eng() const
{
    return eng + kinetic_eng();
}

namespace essential
{
// class make_spring
// {
//    public:
//     sh::stack_map<std::size_t, double, n_max_neighbors> springs;

//     double r_search;
//     double smth;

//     std::size_t id;
//     std::size_t owner_id;

//     vector pos;

//     spring_constant constant;

//     void copyFromFP(const full_ptcl& src);

//     void setPos(const PS::F64vec& src);

//     PS::F64vec getPos() const;

//     PS::F64 getRSearch() const;
// };
class spring
{
   public:
    // sh::stack_map<std::size_t, spring_info, n_max_neighbors> springs;

    std::size_t id;
    std::optional<std::size_t> unit_id;

    double mass;
    double pres;
    double std_dens;
    double dens;
    double smth;
    double r_search;

    // spring_constant spring_const;
    double spring_const;

    vector vel;
    vector pos;
    vector init_pos;

    void copyFromFP(const full_ptcl& src);

    void setPos(const PS::F64vec& src);

    PS::F64vec getPos() const;

    PS::F64 getRSearch() const;
};
class repulsion
{
   public:
    std::size_t id;
    std::optional<std::size_t> unit_id;
    double r_search;
    double dens;
    vector pos;
    void copyFromFP(const full_ptcl& src);
    void setPos(const PS::F64vec& src);
    PS::F64vec getPos() const;
    PS::F64 getRSearch() const;
};
class density
{
   public:
    std::size_t id;
    double mass;
    double dens;
    double smth;
    double r_search;
    vector pos;
    void copyFromFP(const full_ptcl& src);
    void setPos(const PS::F64vec& src);
    PS::F64vec getPos() const;
    PS::F64 getRSearch() const;
};
}  // namespace essential

}  // namespace particles
}  // namespace cfe