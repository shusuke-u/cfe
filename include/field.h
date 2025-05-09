#pragma once

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <sh/syntax.h>
#include <sh/min_max.h>

#include <flags.h>
#include <force.h>
#include <material.h>
#include <interaction.h>

namespace cfe
{

// output data type
struct _key_distribution
{
};
inline constexpr _key_distribution distr{};
struct param
{
};
struct bonds
{
};

// interaction type
struct repulsion_flag
{
};
constexpr inline repulsion_flag repulsion{};

class field : public particle_view<particle_range<particles::full_ptcl>>,
              public std::enable_shared_from_this<field>
{
   private:
    void enable_lazy_evaluation() = delete;
    void disable_lazy_evaluation() = delete;
    particle_view& recompute() = delete;
    filter_type filter() const = delete;
    void remove() = delete;

    using ptcl_t = particles::full_ptcl;

    friend class material;

    double time_;
    double time_step_;
    double CFL_factor_;
    double kernel_support_ratio_;
    // double spring_length_ratio_;
    // double spring_const_rest_length_ratio_;
    // double mean_ptcl_sep_;
    // double init_dens_;
    // double dens_diff_std_;

    /// @brief Box for periodic boundary conditions. The component in the direction that is not
    /// a periodic boundary is `std::nullopt`.
    vector_t<std::optional<sh::min_max<double>>> periodic_domain_;

#ifdef USE_ORG_PERIODIC_BOUNDARY
    PS::TreeForForceShort<interaction::result::spring, particles::essential::spring,
                          particles::essential::spring>::Gather spring_tree_;

    PS::TreeForForceShort<interaction::result::repulsion, particles::essential::repulsion,
                          particles::essential::repulsion>::Gather repulsion_tree_;

    PS::TreeForForceShort<interaction::result::density, particles::essential::density,
                          particles::essential::density>::Gather dens_tree_;
#else
    PS::TreeForForce<PS::SEARCH_MODE_GATHER, interaction::result::spring,
                     particles::essential::spring, particles::essential::spring,
                     PS::MomentShort, PS::MomentShort, PS::SuperParticleBase,
                     PS::CALC_DISTANCE_TYPE_NEAREST_XYZ>
        spring_tree_;

    PS::TreeForForce<PS::SEARCH_MODE_GATHER, interaction::result::repulsion,
                     particles::essential::repulsion, particles::essential::repulsion,
                     PS::MomentShort, PS::MomentShort, PS::SuperParticleBase,
                     PS::CALC_DISTANCE_TYPE_NEAREST_XYZ>
        repulsion_tree_;

    PS::TreeForForce<PS::SEARCH_MODE_GATHER, interaction::result::density,
                     particles::essential::density, particles::essential::density,
                     PS::MomentShort, PS::MomentShort, PS::SuperParticleBase,
                     PS::CALC_DISTANCE_TYPE_NEAREST_XYZ>
        dens_tree_;
#endif

    PS::DomainInfo dinfo_;

    std::vector<std::shared_ptr<material>> materials_;

    template <class Pred>
    void for_each(Pred pred);
    template <class Pred>
    void for_each(Pred pred) const;

    template <typename Pred>
    void for_each_material(Pred pred) const;

    void kick(const double time_step);
    void drift(const double time_step);
    void update_eng(const double time_step);
    void scatter();
    void clear_force_for_fixed_ptcl();
    void update_time_step();
    void equalize_density() const;
    void project_to_actual_dim();
    void apply_periodic_boundary();
    void set_periodic_domain();
    template <class Proj>
    void check_large_value(Proj&& proj, double bound) const;

    bool first_call_of_compute_force_ = true;
    bool first_call_of_compute_force_repulsion_ = true;
    bool first_call_of_compute_density_ = true;

   public:
    field() = delete;
    field(sh::syntax& conf);

    // ~Field();

    field(const field&) = delete;
    field(field&&) = delete;
    field& operator=(const field&) = delete;
    field& operator=(field&&) = delete;

    void initialize_particle_info();
    void finish_setup();
    void next_step();
    void compute_force();
    void compute_force(repulsion_flag, double factor);
    void compute_pressure();
    void compute_density();
    void next_step(repulsion_flag, double factor);
    void relax(sh::syntax& conf);
    void report_n_ptcl() const;
    void write(std::ostream& os) const;
    void write(const std::string& path) const;
    void write(std::ostream&& os, bonds);
    void write(std::ostream& os, param);
    void write(sh::syntax& conf, param);
    void check_NaN_and_inf() const;
    void write_springs(std::ostream&& os);

    double time() const;
    double time_step() const;
    double mean_square_acc() const;
    double total_eng() const;
    double kinetic_eng() const;
    double sound_speed() const;

    orthotope periodic_domain() const;

    /// @return the average of the distance to nearest neighbors
    double nearest_spacing() const;

    /// @return the mean particle spacing
    double mean_spacing() const;

    std::size_t n_ptcl() const;
    std::size_t n_ptcl_loc() const;

    std::shared_ptr<material> make_material(sh::syntax& conf);
};

}  // namespace cfe