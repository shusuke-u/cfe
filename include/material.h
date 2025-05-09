#pragma once

#include <set>
#include <vector>
#include <functional>
#include <unordered_map>
#include <optional>

#include <sh/syntax.h>

#include <particle_range.h>
#include <particle_view.h>

namespace cfe
{

class field;

class material : public recursive_particle_view<particle_range<particles::full_ptcl>>
{
   private:
    // using base = recursive_particle_view<particle_range<particles::full_ptcl>>;

    std::weak_ptr<field> field_;
    std::size_t id_;
    orthotope init_domain_;
    double init_mean_spacing_;
    double init_dens_;

   private:
    util::arrangement read_file(sh::syntax& conf, std::vector<particles::full_ptcl>& ptcl_buf);
    util::arrangement make_crystal_distribution(sh::syntax& conf,
                                                std::vector<particles::full_ptcl>& ptcl_buf);
    void initialize(sh::syntax& conf);
    void make_basic_subsets(sh::syntax& conf);

   public:
    /// @brief Prohibit copying.
    material(const material&) = delete;
    material& operator=(const material&) = delete;

    /// @brief Allow moving.
    material(material&&) = default;
    material& operator=(material&& src) = default;

    material(sh::syntax& conf, std::shared_ptr<particle_range<particles::full_ptcl>> ptcl_range,
             std::shared_ptr<geometry> geo, std::weak_ptr<field> field);

    std::size_t index() const;
    double init_mean_spacing() const;
    double equalize_density() const;
    orthotope init_domain() const;
    orthotope domain_by_edge_ptcls() const;
};

}  // namespace cfe