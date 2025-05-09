#pragma once

#include <vector>
#include <execution>
#include <functional>
#include <unordered_map>

#include <sh/finite.h>

#include <particles.h>
#include <util.h>

namespace cfe
{

/// @brief Contains the geometrical infomation.
struct geometry
{
    std::size_t n_dim;
};

/**
 * @fn
 * Base class for particle aliases.
 * It holds views to particles belonging to itself and provides methods for computing physical
 * values.
 *  */
template <std::ranges::range Range>
class particle_view
{
   public:
    using range_type = Range;
    using value_type = std::ranges::range_value_t<Range>;

    // class filter_type : public std::function<bool(const value_type&)>
    // {
    //    public:
    //     using std::function<bool(const value_type&)>::function;

    //     bool
    //     operator()(const value_type& ptcl) const
    //     {
    //         if(!(*this)) return true;  // avoid throwing std::bad_function_call

    //         // without calling the base function explicitly, infinite recursion occurs
    //         return std::function<bool(const value_type&)>::operator()(ptcl);
    //     }
    // };
    using filter_func_type = std::function<bool(const value_type&)>;

    class filter_type : public filter_func_type
    {
        using base = filter_func_type;

       public:
        using base::base;

        bool
        operator()(const value_type& ptcl) const
        {
            if(!base::operator bool()) return true;  // avoid throwing std::bad_function_call
            return base::operator()(ptcl);
        }
    };

   protected:
    /// @brief An index, which particles with the same index are to belong to `*this`.
    /// If this has not been allocated, lazy evaluation is disabled.
    sh::finite<std::size_t> id_;

    /// @brief An actual filter that defines particles belong to `*this`.
    /// Do not modify this.
    filter_type filter_;

    /// @brief A shared pointer to whole particle range.
    std::shared_ptr<range_type> ptcl_range_;

    std::shared_ptr<geometry> geom_;

   private:
    particle_view() = delete;

    particle_view(const particle_view&) = delete;
    particle_view& operator=(const particle_view&) = delete;

    /// @brief Make a `filter_view` of the particles that had been marked as belonging to
    /// `*this`. `filter_` does not considered now.
    /// @return A `filter_view` made.
    auto
    make_view() const
    {
        if(id_)
        {
            filter_type f = [id = *id_](const value_type& p) { return p.belongs_to[id]; };

            return *ptcl_range_ | std::views::filter(f);
        }
        else
        {
            return *ptcl_range_ | std::views::filter(filter_);
        }
    }

   public:
    /// @brief Construct without filtering and lazy evaluation.
    /// @param ptcl_range `shared_ptr` to the whole particle system.
    explicit particle_view(std::shared_ptr<range_type>&& ptcl_range,
                           std::shared_ptr<geometry>&& geom)
        : ptcl_range_{ std::move(ptcl_range) },
          geom_{ std::move(geom) },
          filter_{ [](const auto&) { return true; } },
          id_{}
    {
    }

    /// @brief Construct with a filter and lazy evaluation, then compute the particles belong to
    /// `*this` as an initialization.
    explicit particle_view(std::shared_ptr<range_type>&& ptcl_range,
                           std::shared_ptr<geometry>&& geom, filter_type filter,
                           sh::finite<std::size_t>&& id_for_lazy_eval)
        : ptcl_range_{ std::move(ptcl_range) },
          geom_{ std::move(geom) },
          filter_{ std::move(filter) },
          id_{ std::move(id_for_lazy_eval) }
    {
        enable_lazy_evaluation();  // this should be called after `id_` is initialized
        recompute();
    }

    particle_view(particle_view&& rhs) = default;
    particle_view& operator=(particle_view&&) = default;

    ~particle_view()
    {
        id_.release();
    }

    // const ptcl_alias& subset( filter_type filter ) const { return ptcl_alias{
    // _ptr_to_ptcl_sys, filter }; }

    void
    enable_lazy_evaluation()
    {
        id_.secure();
    }

    void
    disable_lazy_evaluation()
    {
        id_.release();
    }

    /// @brief Reapply the filter to the whole particle and mark the particles belonging to
    /// `*this`.
    void
    recompute()
    {
        if(id_)
        {
            auto set_belonging = [f = filter_, id = *id_](auto& p) { p.belongs_to[id] = f(p); };
            util::for_each(util::parallel, *ptcl_range_, set_belonging);
        }
        else
        {
            cerr << "Warning: Particle filter is recomputed while lazy evaluation is not "
                    "enabled, which dose not make sense. "
                    "Is this an expected call?\n";
        }
    }

    filter_type
    filter() const
    {
        return filter_;
    }

    void
    remove()
    {
        assert(ptcl_range_->size() < std::numeric_limits<int>::max());

        std::vector<int> indices;

        for(const auto& [i, ptcl] : *ptcl_range_ | sh::views::enumerate)
        {
            if(filter_(ptcl)) indices.push_back(static_cast<int>(i));
        }

        ptcl_range_->removeParticle(indices.data(), indices.size());
    }

    /// @brief Perform a read-only operation on all particles belonging to `*this`.
    /// @param pred Functor with `operator(const value_type&)`.
    template <class Pred, class... Args>
    Pred
    for_each(Pred pred, Args&&... args) const
    {
        util::for_each(make_view(), pred, std::forward<Args>(args)...);
        return pred;
    }

    /// @brief Perform a modifying operation on all particles belonging to `*this`.
    /// @param pred Functor with `operator(value_type&)`.
    template <class Pred, class... Args>
    Pred
    for_each(Pred pred, Args&&... args)
    {
        util::for_each(make_view(), pred, std::forward<Args>(args)...);
        return pred;
    }

    template <class Proj, class Value>
        requires std::assignable_from<std::invoke_result_t<Proj, value_type&>, Value>
    void
    fill(Proj proj, const Value& value)
    {
        util::fill(make_view(), value, proj);
    }

    template <class Proj>
    void
    clear(Proj proj)
    {
        util::clear(make_view(), proj);
    }

    template <class Proj = std::identity, class Weight = double>
    decltype(auto)
    sum(Proj proj = {}, Weight weight = 1.0) const
    {
        return util::sum(make_view(), proj, weight);
    }

    template <class Proj, class Weight = double>
    decltype(auto)
    mean(Proj proj, Weight weight = 1.0) const
    {
        return util::mean(make_view(), proj, weight);
    }

    template <class Proj>
    decltype(auto)
    min(Proj proj) const
    {
        return util::min(make_view(), {}, proj);
    }

    template <class Proj>
    decltype(auto)
    max(Proj proj) const
    {
        return util::max(make_view(), {}, proj);
    }

    template <class Proj>
    decltype(auto)
    min_element(Proj proj) const
    {
        return util::min_element(make_view(), {}, proj);
    }

    template <class Proj>
    decltype(auto)
    max_element(Proj proj) const
    {
        return util::max_element(make_view(), {}, proj);
    }

    util::arrangement
    box_info() const
    {
        return util::box_info(make_view(), geom_->n_dim);
    }

    /// @return the mean pressure of `*this`
    double
    pres() const
    {
        return mean(&particles::full_ptcl::pres);
    }

    /// @return the mean density of `*this`
    double
    dens() const
    {
        return mean(&particles::full_ptcl::dens);
    }

    /// @return the total mass of `*this`
    double
    mass() const
    {
        return sum(&particles::full_ptcl::mass);
    }

    double
    ptcl_mass() const
    {
        return mean(&particles::full_ptcl::mass);
    }

    /// @return the center of mass of `*this`
    vector
    pos() const
    {
        return mean(&particles::full_ptcl::pos, &particles::full_ptcl::mass);
    }

    /// @return the total momentum over total mass of `*this`
    vector
    vel() const
    {
        return mean(&particles::full_ptcl::vel, &particles::full_ptcl::mass);
    }

    /// @return the total momentum of `*this`
    vector
    mom() const
    {
        return sum(&particles::full_ptcl::mom);
    }

    vector
    angular_vel() const
    {
        return mean(&particles::full_ptcl::angular_vel, &particles::full_ptcl::mass);
    }

    vector
    angular_mom(vector center = vector::zero()) const
    {
        return sum([center](const auto& p) { return p.angular_mom(center); });
    }

    /// @return the number of particles belonging to `*this`
    std::size_t
    size() const
    {
        return util::size(make_view());
    }

    /// @brief Shift positions of all particles belonging to `*this`.
    /// @param shift the distance to be traveled.
    /// @return A reference to `*this`.
    void
    drift(vector shift) const
    {
        for_each([shift](auto& ptcl) { ptcl.pos += shift; });
    }

    /// @brief Add uniform velocity to those of all particles belonging to `*this`.
    /// @param vel the velocity to be added.
    void
    kick(vector vel) const
    {
        for_each([vel](auto& ptcl) { ptcl.vel += vel; });
    }

    /// @brief Set the center of mass of `*this`.
    /// @param src the position to be set.
    void
    set_pos(vector src) const
    {
        drift(src - pos());
    }

    /// @brief Set uniform velocity to those of all particles belonging to `*this`.
    /// @param vel the velocity to be set.
    void
    set_vel(vector vel) const
    {
        fill(&value_type::vel, vel);
    }

    /// @brief Clear velocity (set zero) to those of all particles belonging to `*this`.
    void
    clear_vel()
    {
        clear(&value_type::vel);
    }

    /// @brief Fix all particles belonging to `*this`.
    template <class... Args>
    void
    fix(Args... args)
    {
        for_each([=](auto& ptcl) { (ptcl.is_fixed[args] = ... = true); });
    }

    /// @brief Unfix all particles belonging to `*this`.
    template <class... Args>
    void
    unfix(Args... args)
    {
        for_each([=](auto& ptcl) { (ptcl.is_fixed[args] = ... = false); });
    }

    /// @return A boolean value for whether `*this` contains at least one particle or not.
    bool
    empty() const
    {
        return size() == 0;
    }

    void
    write(const std::string& path, double time) const
    {
        std::ofstream file(path);
        if(!file) throw sh::runtime_error{} << "Failed to open file: " << path;
        write(file, time);
    }

    void
    write(std::ostream& os, double time) const
    {
        util::write(os, make_view(), { size(), time });
    }

    double
    spacing() const
    {
        return util::box_info(make_view(), geom_->n_dim).spacing;
    }

    auto
    ptcl(std::size_t id) const
    {
        auto ptcl = util::find(make_view(), id, [](const auto& p) { return p.id; });
        auto found = ptcl.has_value();
        auto n_found = PS::Comm::getSum((int)found);
        if(0 == n_found) throw sh::INVALID_VALUE(id) << "particle not found";
        if(2 <= n_found) throw sh::exception{} << "duplicated id";
        auto rank_found = PS::Comm::getSum(found ? rank : 0UL);
        value_type ptcl_buf;
        if(found) ptcl_buf = ptcl.value();
        PS::Comm::broadcast(&ptcl_buf, 1, rank_found);
        return ptcl_buf;
    }

    template <class Pred, class Proj = std::identity>
    std::vector<value_type>
    find(Pred pred, Proj proj = {}) const
    {
        return util::find_if_all(make_view(), pred, proj);
    }

    void
    scale(double factor, std::size_t dim)
    {
        auto center = pos();
        auto shift = [&](auto& p) { p.pos[dim] += (p.pos[dim] - center[dim]) * (factor - 1); };
        for_each(shift);
    }

    template <class Proj = vector particles::full_ptcl::*>
    void
    make_random_shift(double magnitude, Proj proj = &particles::full_ptcl::pos) const
    {
        util::make_random_shift(make_view(), proj, magnitude);
    }
};

template <std::ranges::range Range>
class recursive_particle_view : public particle_view<Range>
{
   public:
    using base = particle_view<Range>;
    using typename base::filter_type;
    using typename base::range_type;
    using typename base::value_type;
    using subset_type = recursive_particle_view;

    /// @brief A map of stored sub aliases belong to `*this`.
    std::unordered_map<std::string, std::unique_ptr<subset_type>> subset_list_;

   public:
    using base::base;

    /// @brief Create a sub alias of a subset of `*this` which satisfies the given conditions.
    /// @param filter Functor with `bool operator()(const value_type&)`
    /// @return A reference to the part.
    subset_type&
    make_subset(const std::string& key, filter_type filter)
    {
        auto copied_range = base::ptcl_range_;
        auto copied_geom = base::geom_;

        filter_type new_filter = [this_filter = base::filter_, filter](const value_type& ptcl)
        { return this_filter(ptcl) && filter(ptcl); };

        auto new_subset = std::make_unique<subset_type>(
            std::move(copied_range), std::move(copied_geom), new_filter, base::id_.take());

        auto [it, success] = subset_list_.emplace(key, std::move(new_subset));

        if(!success)
            throw sh::runtime_error{}
                << "Failed to emplace a new object of subset with the key: " << key;

        return *it->second;
    }

    subset_type&
    make_subset(filter_type filter)
    {
        std::string key;
        do
        {
            static std::random_device rd;
            static std::mt19937 gen(rd());
            static std::uniform_int_distribution<std::size_t> dis(0, 1000000);
            key = "__auto" + std::to_string(dis(gen));
        } while(subset_list_.find(key) != subset_list_.end());

        return make_subset(key, filter);
    }

    subset_type&
    subset(const std::string& key)
    {
        auto it = subset_list_.find(key);

        if(it == subset_list_.end())
            throw sh::runtime_error{} << "No subset found with the key: " << key;

        return *it->second;
    }

    const subset_type&
    subset(const std::string& key) const
    {
        auto it = subset_list_.find(key);

        if(it == subset_list_.end())
            throw sh::runtime_error{} << "No subset found with the key: " << key;

        return *it->second;
    }

    const subset_type&
    operator[](const std::string& key) const
    {
        return subset(key);
    };

    subset_type&
    operator[](const std::string& key)
    {
        return subset(key);
    };

    subset_type
    operator&&(const subset_type& rhs) const
    {
        auto copied_range = base::ptcl_range_;
        auto copied_geom = base::geom_;

        auto new_filter = [this_filter = base::filter_, rhs_filter = rhs.base::filter_](
                              const value_type& p) { return this_filter(p) && rhs_filter(p); };

        return subset_type{ std::move(copied_range), std::move(copied_geom), new_filter,
                            base::id_.take() };
    }

    subset_type
    operator||(const subset_type& rhs) const
    {
        auto copied_range = base::ptcl_range_;
        auto copied_geom = base::geom_;

        auto new_filter = [this_filter = base::filter_, rhs_filter = rhs.base::filter_](
                              const value_type& p) { return this_filter(p) || rhs_filter(p); };

        return subset_type{ std::move(copied_range), std::move(copied_geom), new_filter,
                            base::id_.take() };
    }

    subset_type
    operator!() const
    {
        auto copied_range = base::ptcl_range_;
        auto copied_geom = base::geom_;

        // The `filter_` object at the time this function is called is used.
        // That is, any future changes to `this->filter_` will not be
        // reflected in the return object of this function. But it is ok because
        // `this->filter_` is constant.
        auto new_filter = [filter = base::filter_](const value_type& p) { return !filter(p); };

        return subset_type{ std::move(copied_range), std::move(copied_geom), new_filter,
                            base::id_.take() };
    }

    /// @brief Release the indexes reserved to mark particles belonging to `this*` and subsets
    /// of `*this`.
    void
    disable_lazy_evaluation_all()
    {
        for(const auto& [key, sub] : subset_list_) sub->disable_lazy_evaluation_all();
        base::disable_lazy_evaluation();
    }

    void
    enable_lazy_evaluation_all()
    {
        for(const auto& [key, sub] : subset_list_) sub->enable_lazy_evaluation_all();
        base::enable_lazy_evaluation_all();
    }
};

}  // namespace cfe