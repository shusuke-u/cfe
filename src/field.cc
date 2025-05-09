#include <ranges>
#include <random>
#include <numeric>
#include <algorithm>
#include <filesystem>

#include <sh/syntax.h>
#include <sh/ranges.h>

#include <eos.h>
#include <field.h>
#include <parallel.h>

namespace cfe
{

field::field(sh::syntax& conf)
    : particle_view{ std::make_shared<particle_range<particles::full_ptcl>>(),
                     std::make_shared<geometry>(
                         (std::size_t)conf["n_dim"].value_or(vector::size())) },
      time_{ conf["time_begin"].value_or(0) },
      CFL_factor_{ conf["CFL_factor"].value_or(0.1) },
      kernel_support_ratio_{ conf["kernel_support_ratio"].value_or(3) }
{
    if(geom_->n_dim < 1 or vector::size() < geom_->n_dim)
    {
        throw sh::INVALID_VALUE(geom_->n_dim)
            << "Out of the supported range [1, " << vector::size() << "].";
    }
    if(vector::size() < geom_->n_dim)
    {
        throw sh::INVALID_VALUE(geom_->n_dim)
            << "Exceeded the capacity of the vector class: " << vector::size() << ".";
    }

    // Initialize `PS::DomainInfo`.
    dinfo_.initialize();
    vector_t<bool> periodic = conf["periodic"].value_or(vector_t<bool>(false));
    conf["periodic"] = periodic.truncate(geom_->n_dim, false);
#ifdef USE_ORG_PERIODIC_BOUNDARY
    dinfo_.setBoundaryCondition(util::convert(periodic));
#else
    dinfo_.setBoundaryCondition(PS::BOUNDARY_CONDITION_OPEN);
#endif

    // Set the periodic domain manually.
    const auto given_domain = conf.optional<orthotope>("periodic_domain");
    for(auto a : std::views::iota(0UL, geom_->n_dim))
    {
        if(!periodic[a]) continue;

        if(given_domain)
        {
            periodic_domain_[a].emplace(given_domain->lower()[a], given_domain->upper()[a]);
        }
        else
        {
            periodic_domain_[a].emplace(std::numeric_limits<double>::max(),
                                        std::numeric_limits<double>::lowest());
        }

        dinfo_.setPosRootDomain1D(given_domain->lower()[a], given_domain->upper()[a], a);
    }

    // Initialize `PS::ParticleSystem`.
    ptcl_range_->initialize();
}

std::shared_ptr<material>
field::make_material(sh::syntax& conf)
{
    auto mat = std::make_shared<material>(conf, ptcl_range_, geom_, shared_from_this());
    auto domain = mat->init_domain();
    // Extend the periodic domain to include the new region
    for(auto [i, lw, up] : sh::views::zip(domain.lower(), domain.upper()) |
                               std::views::take(geom_->n_dim) | sh::views::enumerate)
        cout << std::string(2, '\t') << sh::dir(i) << " : " << lw << ' ' << up << '\n';

    return mat;
}

orthotope
field::periodic_domain() const
{
    orthotope domain;
    for(auto [src, lw, up] : sh::views::zip(periodic_domain_, domain.lower(), domain.upper()) |
                                 std::views::take(geom_->n_dim))
    {
        lw = src ? src->min : std::numeric_limits<double>::lowest();
        up = src ? src->max : std::numeric_limits<double>::max();
    }
    return domain;
}

// double
// field::nearest_spacing() const
// {
//     auto get_len = [](const sh::stack_map<std::size_t, spring_info,
//     n_max_neighbors>::iterator::value_type& p) { return p.second.length; };

//     auto min_spacing = [&](const particles::full_ptcl& ptcl)
//     {
//         if(ptcl.springs.size() == 0) throw sh::EXCEPTION << "Neighbor not found.";

//         auto it = std::ranges::min_element(ptcl.springs, std::less{}, get_len);

//         return get_len(*it);
//     };

//     auto nearest_spacing = [&](const particles::full_ptcl& ptcl)
//     {
//         auto limit = min_spacing(ptcl) * 1.1;

//         auto nearests =
//             ptcl.springs | std::views::transform(get_len) | std::views::filter([&](double
//             len) { return len < limit; });

//         return std::reduce(std::execution::par_unseq, nearests.begin(), nearests.end()) /
//         std::ranges::distance(nearests);
//     };

//     return mean(nearest_spacing);
// }

double
field::mean_spacing() const
{
    return std::pow(util::volume(box_info().domain, geom_->n_dim) / n_ptcl(),
                    1.0 / geom_->n_dim);
}

void
field::project_to_actual_dim()
{
    auto trunc = [this](particles::full_ptcl& p)
    {
        p.pos.truncate(geom_->n_dim);
        p.vel.truncate(geom_->n_dim);
    };
    for_each(trunc);
}

template <class Pred>
void
field::for_each(Pred pred)
{
    util::for_each(*ptcl_range_, pred);
}

template <class Pred>
void
field::for_each(Pred pred) const
{
    util::for_each(*ptcl_range_, pred);
}

template <typename Pred>
void
field::for_each_material(Pred pred) const
{
    constexpr auto deref = [](const auto& ptr) -> material& { return *ptr; };
    util::for_each(materials_, pred, deref);
}

void
field::kick(const double time_step)
{
    for_each([time_step](particles::full_ptcl& ptcl) { ptcl.vel += time_step * ptcl.acc(); });
}

void
field::drift(const double time_step)
{
    for_each([time_step](particles::full_ptcl& ptcl) { ptcl.pos += time_step * ptcl.vel; });
}

void
field::update_eng(const double time_step)
{
    for_each([time_step](particles::full_ptcl& ptcl) { ptcl.eng += time_step * ptcl.eng_dot; });
}

void
field::scatter()
{
    auto pred = [](particles::full_ptcl& p) { p.force -= p.dissipate * p.vel; };
    for_each(pred);
}

void
field::clear_force_for_fixed_ptcl()
{
    for(auto dim : vector::indices())
    {
        auto pred = [&](auto& p)
        {
            if(p.is_fixed[dim]) p.force[dim] = 0;
        };

        for_each(pred);
    }
}

void
field::finish_setup()
{
    project_to_actual_dim();
    set_periodic_domain();

    ptcl_range_->adjustPositionIntoRootDomain(dinfo_);
    dinfo_.decomposeDomainAll(*ptcl_range_);
    ptcl_range_->exchangeParticle(dinfo_);

    // call this before `initialize_particle_info()` and `flatten_density()`
    compute_density();

#ifdef USE_VOLUME_BASED_SPH
    // call this before `initialize_particle_info()`
    equalize_density();
#endif

    // call this before `compute_force()`
    initialize_particle_info();

    // call this before `compute_force()`
    compute_pressure();

    // call this after `initialize_particle_info()`
    clear(&particles::full_ptcl::force);
    compute_force();  // for making tree

    update_time_step();
}

void
field::compute_pressure()
{
    constexpr auto pred = [](particles::full_ptcl& p)
    { p.pres = EoS::cold_term(p.bulk_sound_sq, p.dens, p.std_dens); };
    for_each(pred);
}

void
field::set_periodic_domain()
{
    // auto periodic = util::convert(dinfo_.getBoundaryCondition());
    if(!std::ranges::empty(periodic_domain_ | sh::views::keep_true)) cout << "periodic axes:\n";

    for(auto [dim, axis] : periodic_domain_ | sh::views::enumerate)
    {
        if(!axis) continue;
        // assert(periodic[dim] == (bool)M_periodic_domain[dim]);
        dinfo_.setPosRootDomain1D(axis->min, axis->max, dim);
        cout << std::string(1, '\t') << sh::dir(dim) << " : " << *axis << '\n';
    }
}

template <class Proj>
void
field::check_large_value(Proj&& proj, double bound) const
{
    using value_type = std::decay_t<std::invoke_result_t<Proj, particles::full_ptcl>>;

    auto out_of_bound = [proj, bound](const particles::full_ptcl& ptcl)
    {
        auto var = std::invoke(proj, ptcl);

        if constexpr(std::ranges::range<value_type>)
        {
            auto pred = [bound](auto&& v) { return bound < std::abs(v); };
            return std::ranges::any_of(var, pred);
        }
        else
        {
            return bound < std::abs(var);
        }

        return false;
    };
    if(const auto& ptcl = util::find_if(*ptcl_range_, out_of_bound))
    {
        sh::runtime_error err;
        err << "Particle with too large value found:\n";
        err << particles::full_ptcl::header << '\n';
        err << *ptcl << '\n';
        throw err;
    }
}

void
field::equalize_density() const
{
    for_each_material(&material::equalize_density);
}

void
field::initialize_particle_info()
{
    auto init = [](particles::full_ptcl& p)
    {
        p.init_pos = p.pos;
        p.init_vel = p.vel;
        p.std_dens = p.dens;
    };
    for_each(init);
}

void
field::update_time_step()
{
    auto time_step = [&](const auto& p) { return p.time_step(); };
    time_step_ = min(time_step) * CFL_factor_;
}

void
field::next_step()
{
    update_time_step();

    clear(&particles::full_ptcl::force);
    // calculate spring force
    compute_force();
    clear_force_for_fixed_ptcl();

    // 1st order symplectic scheme
    scatter();
    // set_vel();
    kick(time_step_);
    drift(time_step_);
    update_eng(time_step_);
    time_ += time_step_;

    check_large_value(&particles::full_ptcl::pos, 1.5);

    ptcl_range_->adjustPositionIntoRootDomain(dinfo_);

    // Note: Do not call following two methods in order to reuse interaction list
#ifdef USE_ORG_PERIODIC_BOUNDARY
    dinfo_.decomposeDomainAll(*ptcl_range_);
    ptcl_range_->exchangeParticle(dinfo_);
#endif

    // compute density and update smoothing length
    compute_density();
    compute_pressure();
}

void
field::compute_force()
{
    // Make spring tree.
    if(first_call_of_compute_force_)
    {
        // // Here `r_search` must be same as the max spring length.
        // PS::TreeForForceShort<interaction::result::make_spring,
        // particles::essential::make_spring,
        //                       particles::essential::make_spring>::Gather make_spr;
        // make_spr.initialize(n_ptcl());
        // // This call does not set spring infomation not at all, since a value of the mean
        // particle spacing is required
        // // which is computed using the spring list.
        // make_spr.calcForceAllAndWriteBack(interaction::functor::make_spring{}, *ptcl_range_,
        // dinfo_, true,
        //                                   PS::MAKE_LIST_FOR_REUSE);
        // // Now one can compute the spacing.
        // const auto spacing = nearest_spacing();
        // // Set the accurate spring coefficients.
        // make_spr.calcForceAllAndWriteBack(
        //     interaction::functor::make_spring{ M_geom->n_dim, spacing, M_spring_length_ratio,
        //     spring_const_rest_length_ratio_
        //     }, *ptcl_range_, dinfo_, true, PS::REUSE_LIST);
        // // Enlarge `r_search` by `margin` so that all neighbors are certainly included when
        // the tree is updated in each step. auto margin = spacing * 2; auto sides =
        // M_periodic_domain | sh::views::keep_true | std::views::transform([](const auto& e) {
        // return e->bound();
        // }); std::optional<double> bound; if(not sides.empty()) bound = 0.5 *
        // std::ranges::min(sides); for_each(
        //     [&](particles::full_ptcl& p)
        //     {
        //         p.r_search_for_spring += margin;
        //         if(bound and *bound < p.r_search_for_spring) p.r_search_for_spring = *bound;
        //     });

        spring_tree_.initialize(n_ptcl());
        spring_tree_.calcForceAllAndWriteBack(
            interaction::functor::spring_temp{ geom_->n_dim, periodic_domain_ }, *ptcl_range_,
            dinfo_, true, PS::MAKE_LIST_FOR_REUSE);

        first_call_of_compute_force_ = false;
    }
    else
    {
#ifdef USE_ORG_PERIODIC_BOUNDARY
        spring_tree_.calcForceAllAndWriteBack(
            interaction::functor::spring_temp{ geom_->n_dim, periodic_domain_ }, *ptcl_range_,
            dinfo_, true, PS::MAKE_LIST);
#else
        spring_tree_.calcForceAllAndWriteBack(
            interaction::functor::spring_temp{ geom_->n_dim, periodic_domain_,
                                               spring_const_rest_length_ratio_, mean_ptcl_sep_,
                                               init_dens_, dens_diff_std_ },
            *ptcl_range_, dinfo_, true, PS::REUSE_LIST);
#endif
    }
}

void
field::next_step(repulsion_flag, double factor)
{
    // update_time_step();
    // compute_force();
    clear(&particles::full_ptcl::force);
    compute_force(repulsion, factor);

    clear_force_for_fixed_ptcl();
    scatter();
    // set_vel();
    kick(time_step_);
    drift(time_step_);
    update_eng(time_step_);
    time_ += time_step_;

    check_large_value(&particles::full_ptcl::pos, 1e2);

    ptcl_range_->adjustPositionIntoRootDomain(dinfo_);
}

void
field::compute_force(repulsion_flag, double factor)
{
    if(first_call_of_compute_force_repulsion_)
    {
        repulsion_tree_.initialize(n_ptcl());
        first_call_of_compute_force_repulsion_ = false;
    }

    repulsion_tree_.calcForceAllAndWriteBack(
        interaction::functor::repulsion{ periodic_domain_, factor }, *ptcl_range_, dinfo_, true,
        PS::MAKE_LIST);
}

void
field::compute_density()
{
    if(first_call_of_compute_density_)
    {
        dens_tree_.initialize(n_ptcl());
        first_call_of_compute_density_ = false;
    }

    auto update_smth = [this](particles::full_ptcl& p)
    {
        p.smth = std::pow(p.mass / p.dens, 1.0 / geom_->n_dim);
        p.r_search_for_dens = p.smth * kernel_support_ratio_;
    };

    dens_tree_.calcForceAllAndWriteBack(
        interaction::functor::density{ geom_->n_dim, periodic_domain_ }, *ptcl_range_, dinfo_,
        true, PS::MAKE_LIST);

    // This must be called after the above calculation, since, at the first call of this method,
    // the densities are not set. The smoothing lengths are set in the particles'
    // initialization.
    // for_each(update_smth);
}

void
field::apply_periodic_boundary()
{
    auto apply_pb = [this](particles::full_ptcl& p)
    {
        for(auto [ax, x] :
            sh::views::zip(periodic_domain_, p.pos) | std::views::take(geom_->n_dim))
        {
            while(x < ax->min) x += ax->min;
            while(ax->max < x) x -= ax->max;
        }
    };
    for_each(apply_pb);

    // for(auto a : std::views::iota(0ul, geom_->n_dim))
    // {
    //     const auto axis = periodic_domain_[a];

    //     for_each(
    //         [&](particles::full_ptcl& ptcl)
    //         {
    //             auto pos = ptcl.pos;
    //             while(pos[a] < axis->min) pos[a] += axis->min;
    //             while(axis->max < pos[a]) pos[a] -= axis->max;
    //             ptcl.pos = pos;
    //         });
    // }
}

void
field::check_NaN_and_inf() const
{
    auto nan_or_inf = [](vector v) { return v.is_nan() or v.is_inf(); };
    auto ptcls = find(nan_or_inf, &particles::full_ptcl::acc);

    if(not ptcls.empty())
    {
        cerr << "Particles with NaN or Inf accelerations were found.\n";
        cerr << particles::full_ptcl::header << '\n';
        for(const auto& p : ptcls) std::cerr << p << '\n';
    }
}

void
field::write_springs(std::ostream&& os)
{
    using interaction::functor::get_springs;

    auto springs = std::make_shared<get_springs::springs>();

    spring_tree_.calcForceAllAndWriteBack(get_springs{ periodic_domain_, springs },
                                          *ptcl_range_, dinfo_, true, PS::REUSE_LIST);

    assert(springs);
    auto springs_all = par::all_gather(*springs);
    if(in_root_proc)
    {
        os << get_springs::spring::header << '\n';
        for(const auto& s : springs_all) os << s << '\n';
    }
}

void
field::write(const std::string& path) const
{
    particle_view::write(path, time_);
}

void
field::write(std::ostream& os) const
{
    particle_view::write(os, time_);
}

#if 0
void field::dump_to( std::ostream& os, param ) {
    if( os.bad() ) PARTICLE_SIMULATOR_PRINT_ERROR( "Given `ostream` is bad." )
    if( _init_side == 0 ) PARTICLE_SIMULATOR_PRINT_ERROR( "Initial side of material has not been set." )
    const auto curr_side = side_length();
    const auto strain = ApplyEach( []( auto a, auto b ) { return a / b - 1; }, curr_side, _init_side );
    const auto poisson = -strain / strain[0];  // x-component has no meaning
    const auto stress = side_stress();
    const auto young = stress[0] / strain[0];
    if( root_proc ) os << strain << "   " << stress << "   " << poisson << "   " << young << "   " << '\n';
}

void field::dump_to( sh::syntax& conf, param ) {
    if( _init_side == 0 ) PARTICLE_SIMULATOR_PRINT_ERROR( "Initial side of material has not been set." )
    const auto curr_side = side_length();
    const auto strain = ApplyEach( []( auto a, auto b ) { return a / b - 1; }, curr_side, _init_side );
    const auto poisson = -strain / strain[0];  // x-component has no meaning
    const auto stress = side_stress();
    const auto young = stress[0] / strain[0];
    const auto msa = mean_square_acc();
    conf.set( "poisson_y", poisson[1] ).set( "poisson_z", poisson[2] ).set( "young", young ).set( "mean_sq_acc", msa );
}
#endif

void
field::write(std::ostream&& os, bonds)
{
    const auto pair_list = std::make_shared<interaction::functor::visualize::list_type>();
    spring_tree_.calcForceAllAndWriteBack(interaction::functor::visualize{ pair_list },
                                          *ptcl_range_, dinfo_, true, PS::MAKE_LIST);
    // dump
    if(in_root_proc) os << file_header{ n_ptcl(), time_ };
    std::stringstream ss;
    for(const auto& bond : *pair_list)
        ss << rank << sh::del << bond.first.pos << bond.second.pos << '\n';
    util::write(os, ss);
}

void
field::report_n_ptcl() const
{
    PS::Comm::barrier();
    std::cout << "(" << rank << "/" << n_procs << ") n_ptcl: " << n_ptcl_loc() << '\n';
    std::cout << std::endl;
}

double
field::time() const
{
    return time_;
}

double
field::kinetic_eng() const
{
    return sum(&particles::full_ptcl::kinetic_eng);
}

double
field::sound_speed() const
{
    auto proj = [&](const auto& p) { return p.sound_speed(); };
    return mean(proj);
}

double
field::total_eng() const
{
    return sum(&particles::full_ptcl::total_eng);
}

double
field::mean_square_acc() const
{
    auto square_acc = [](const particles::full_ptcl& p) { return p.acc().square(); };
    return mean(square_acc);
}

double
field::time_step() const
{
    return time_step_;
}

std::size_t
field::n_ptcl() const
{
    return ptcl_range_->getNumberOfParticleGlobal();
}

std::size_t
field::n_ptcl_loc() const
{
    return ptcl_range_->getNumberOfParticleLocal();
}

}  // namespace cfe