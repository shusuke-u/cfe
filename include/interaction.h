#pragma once

#include <mutex>
#include <shared_mutex>

#include <defs.h>
#include <kernel.h>
#include <particles.h>

namespace cfe
{
namespace interaction
{
namespace functor
{

// class make_spring
// {
//    private:
//     const std::size_t M_n_dim;
//     const std::optional<double> M_mean_separation;
//     const std::optional<double> M_spr_len_rat;
//     const std::optional<double> M_const_rest_len_ratio;

//    public:
//     /// @brief Instance constructed with this never does nothing in `operator()()`.
//     make_spring();

//     make_spring(std::size_t n_dim, double mean_spacing, double spring_length_ratio, double
//     spring_const_rest_length_ratio);

//     void operator()(const particles::essential::make_spring* const list_ptcl_i, const PS::S32
//     n_ptcl_i,
//                     const particles::essential::make_spring* const list_ptcl_j, const PS::S32
//                     n_ptcl_j, interaction::result::make_spring* results) const;
// };

class spring_temp
{
    const std::size_t n_dim_;
    const vector_t<std::optional<sh::min_max<double>>> periodic_domain_;

   public:
    spring_temp(std::size_t n_dim, vector_t<std::optional<sh::min_max<double>>> domain)
        : n_dim_{ n_dim }, periodic_domain_{ domain }
    {
    }
    void operator()(const particles::essential::spring* const list_ptcl_i,
                    const PS::S32 n_ptcl_i,
                    const particles::essential::spring* const list_ptcl_j,
                    const PS::S32 n_ptcl_j, interaction::result::spring* results) const;
    void calc_force(const particles::essential::spring& pi,
                    const particles::essential::spring& pj,
                    interaction::result::spring& result) const;
};

class repulsion
{
    vector_t<std::optional<sh::min_max<double>>> periodic_domain_;
    // double mean_separation_;
    double repulsion_factor_;
    // double dens_std_;

   public:
    repulsion(vector_t<std::optional<sh::min_max<double>>> domain, double repulsion_constant)
        : periodic_domain_{ domain }, repulsion_factor_{ repulsion_constant }
    {
    }
    void operator()(const particles::essential::repulsion* const list_ptcl_i,
                    const PS::S32 n_ptcl_i,
                    const particles::essential::repulsion* const list_ptcl_j,
                    const PS::S32 n_ptcl_j, interaction::result::repulsion* results) const;
    void calc_force(const particles::essential::repulsion& pi,
                    const particles::essential::repulsion& pj,
                    interaction::result::repulsion& result) const;
};

class get_springs
{
   public:
    struct spring
    {
        std::pair<std::size_t, std::size_t> ptcl_id;
        double rest_length;
        double stretch;
        double constant;
        friend std::ostream&
        operator<<(std::ostream& os, const spring& s)
        {
            os << s.ptcl_id.first << '\t';
            os << s.ptcl_id.second << '\t';
            os << s.rest_length << '\t';
            os << s.constant << '\t';
            os << s.stretch << '\t';
            return os;
        }
        struct header_type
        {
        };
        static constexpr header_type header{};
        friend std::ostream&
        operator<<(std::ostream& os, header_type)
        {
            os << "ptcl.id[0]" << '\t';
            os << "ptcl.id[1]" << '\t';
            os << "rest_length" << '\t';
            os << "constant" << '\t';
            os << "stretch" << '\t';
            return os;
        }
    };
    using springs = std::vector<spring>;

   private:
    std::shared_ptr<springs> springs_;

    // Cannot have `std::mutex` object as is. Since `std::mutex` cannot be copied but
    // `TreeForForce::calcForceAll()` does copy this functor.
    mutable std::shared_mutex springs_mutex_;

    vector_t<std::optional<sh::min_max<double>>> periodic_domain_;

   public:
    get_springs() = delete;
    get_springs(vector_t<std::optional<sh::min_max<double>>> domain,
                std::shared_ptr<springs> springs)
        : periodic_domain_{ domain }, springs_{ std::move(springs) }
    {
    }

    get_springs(const get_springs& src)
    {
        std::shared_lock<std::shared_mutex> lock(src.springs_mutex_);
        springs_ = src.springs_;
        periodic_domain_ = src.periodic_domain_;
    }

    get_springs&
    operator=(const get_springs& src)
    {
        if(this != &src)
        {
            std::unique_lock<std::shared_mutex> lock1(springs_mutex_, std::defer_lock);
            std::shared_lock<std::shared_mutex> lock2(src.springs_mutex_, std::defer_lock);
            std::lock(lock1, lock2);
            springs_ = src.springs_;
            periodic_domain_ = src.periodic_domain_;
        }
        return *this;
    }

    void operator()(const particles::essential::spring* const list_ptcl_i,
                    const PS::S32 n_ptcl_i,
                    const particles::essential::spring* const list_ptcl_j,
                    const PS::S32 n_ptcl_j, interaction::result::spring* results);
};

class spring
{
   public:
    void operator()(const particles::essential::spring* const list_ptcl_i,
                    const PS::S32 n_ptcl_i,
                    const particles::essential::spring* const list_ptcl_j,
                    const PS::S32 n_ptcl_j, interaction::result::spring* results) const;
};

class SPH_like
{
   public:
    void operator()(const particles::essential::spring* const list_ptcl_i,
                    const PS::S32 n_ptcl_i,
                    const particles::essential::spring* const list_ptcl_j,
                    const PS::S32 n_ptcl_j, interaction::result::spring* results) const;
};

class Lanvegin
{
   private:
    const double _coef, _r_int, _r_int_sq;

   public:
    Lanvegin(double coef, double interaction_radius);
    void operator()(const particles::full_ptcl* list_ptcl_i, const PS::S32 n_ptcl_i,
                    const particles::full_ptcl* list_ptcl_j, const PS::S32 n_ptcl_j,
                    interaction::result::spring* const forces) const;
};

class visualize
{
   public:
    using list_type =
        std::vector<std::pair<particles::essential::spring, particles::essential::spring>>;

   private:
    const std::shared_ptr<list_type> _bond_list;

   public:
    visualize() = delete;
    visualize(std::shared_ptr<list_type> bond_list);
    void operator()(const particles::essential::spring* list_ptcl_i, const PS::S32 n_ptcl_i,
                    const particles::essential::spring* list_ptcl_j, const PS::S32 n_ptcl_j,
                    interaction::result::spring* const forces);
};

class density
{
   private:
    vector_t<std::optional<sh::min_max<double>>> periodic_domain_;
    std::size_t n_dim_;

   public:
    density(std::size_t n_dim, vector_t<std::optional<sh::min_max<double>>> domain);

    void operator()(const particles::essential::density* const list_ptcl_i,
                    const PS::S32 n_ptcl_i,
                    const particles::essential::density* const list_ptcl_j,
                    const PS::S32 n_ptcl_j, interaction::result::density* const results) const;
};

class clear
{
   public:
    void operator()(const particles::essential::spring* list_ptcl_i, const PS::S32 n_ptcl_i,
                    const particles::essential::spring* list_ptcl_j, const PS::S32 n_ptcl_j,
                    interaction::result::spring* const forces) const;
};

class do_nothing
{
   public:
    void operator()(const particles::full_ptcl* list_ptcl_i, const PS::S32 n_ptcl_i,
                    const particles::full_ptcl* list_ptcl_j, const PS::S32 n_ptcl_j,
                    interaction::result::spring* const forces) const;
};

}  // namespace functor
}  // namespace interaction
}  // namespace cfe
