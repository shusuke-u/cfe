#include <field.h>
#include <material.h>
#include <crystals.h>

namespace cfe
{

util::arrangement
material::read_file(sh::syntax& conf, std::vector<particles::full_ptcl>& ptcl_buf)
{
    file_header header;
    std::string path = conf["path"];
    util::read(path, ptcl_buf, header);

    // Check the number of the data.
    auto n_ptcl_read_glob = util::size(ptcl_buf);
    if(in_root_proc && n_ptcl_read_glob != header.n_ptcl)
    {
        throw sh::runtime_error{}
            << "The number of the read particles is inconsistent with the file header." << '\n'
            << "n_ptcl_read_glob = " << n_ptcl_read_glob << '\n'
            << "header.n_ptcl    = " << header.n_ptcl << '\n';
    }

    // update n_ptcl
    conf["n_ptcl"] = n_ptcl_read_glob;

    // reuse read time if needed
    if(conf["reuse_time"].value_or(false)) field_.lock()->time_ = header.time;

    if(std::optional<double> eng = conf["eng"])
    {
        cout << std::string(2, '\t') << "eng     : " << *eng << '\n';
        for(auto& p : ptcl_buf) p.eng = *eng;
    }
    if(std::optional<double> dens = conf["dens"])
    {
        cout << std::string(2, '\t') << "dens    : " << *dens << '\n';
        for(auto& p : ptcl_buf) p.dens = *dens;
    }
    auto arr = util::box_info(ptcl_buf, geom_->n_dim);

    cout << std::string(2, '\t') << "path    : " << path << '\n';
    cout << std::string(2, '\t') << "n_ptcl  : " << arr.n_ptcl << '\n';
    cout << std::string(2, '\t') << "domain  : " << arr.domain << '\n';
    cout << std::string(2, '\t') << "spacing : " << arr.spacing << '\n';

    return arr;
}

util::arrangement
material::make_crystal_distribution(sh::syntax& conf,
                                    std::vector<particles::full_ptcl>& ptcl_buf)
{
    std::size_t n_ptcl_new = conf["n_ptcl"];
    orthotope domain = conf["domain"];

    crystals::arrange_fn::result_type arrange;
    const char* structure = conf["structure"];
    if(std::set{ "Cartesian" }.contains(structure))
    {
        arrange = crystals::arrange(crystals::Cartesian(geom_->n_dim), domain, n_ptcl_new,
                                    geom_->n_dim);
    }
    else if(std::set{ "tetrahedron", "densest" }.contains(structure))
    {
        arrange = crystals::arrange(crystals::densest(geom_->n_dim), domain, n_ptcl_new,
                                    geom_->n_dim);
    }
    else if(std::set{ "struct", "stratified" }.contains(structure))
    {
        arrange = crystals::arrange(crystals::tetrahedron(geom_->n_dim), domain, n_ptcl_new,
                                    geom_->n_dim);
    }
    else
    {
        throw sh::INVALID_VALUE(structure);
    }

    auto [distribution, nearest_spacing] = arrange;

    conf["domain"] = domain;
    conf["n_ptcl"] = n_ptcl_new;

    if(in_root_proc)
    {
        double ptcl_mass = conf["dens"] * util::volume(domain, geom_->n_dim) / n_ptcl_new;
        double eng = conf["eng"];
        vector vel = conf["vel"].value_or(vector::zero());
        std::size_t owner_id = conf["index"];

        ptcl_buf.reserve(n_ptcl_new);

        for(auto [id, dstr] : distribution | sh::views::enumerate)
        {
            auto& p = ptcl_buf.emplace_back();
            p.id = id;
            p.unit_id = dstr.unit_id;
            p.pos = domain.lower() + dstr.pos;
            p.vel = vel;
            p.eng = eng;
            p.mass = ptcl_mass;
            p.owner_id = owner_id;
        }
    }

    auto mean_spacing =
        std::pow(util::volume(domain, geom_->n_dim) / n_ptcl_new, 1.0 / geom_->n_dim);

    cout << std::string(2, '\t') << "unit structure : " << structure << '\n';
    cout << std::string(2, '\t') << "# of particles : " << n_ptcl_new << '\n';
    cout << std::string(2, '\t') << "min spacing    : " << nearest_spacing << '\n';
    cout << std::string(2, '\t') << "mean spacing   : " << mean_spacing << '\n';

    return { n_ptcl_new, domain, mean_spacing };
}

void
material::initialize(sh::syntax& conf)
{
    std::size_t owner_id = conf["index"];
    std::string mode = conf["mode"];

    cout << "Constructing material:" << '\n';
    cout << std::string(1, '\t') << "id   : " << owner_id << '\n';
    cout << std::string(1, '\t') << "mode : " << mode << '\n';

    util::arrangement arr;
    std::vector<particles::full_ptcl> ptcl_buf;

    if(mode == "read")
    {
        arr = read_file(conf, ptcl_buf);
    }
    else if(mode == "make")
    {
        arr = make_crystal_distribution(conf, ptcl_buf);
    }
    else
    {
        throw sh::INVALID_VALUE(mode);
    }

    // if(auto read_domain = conf.optional<orthotope>("domain"))
    // {
    //     if(not read_domain->contains(arr.domain))
    //         throw sh::invalid_value{}
    //             << "Manually set domain does not contain the whole of actual distribution.";
    //     arr.domain = *read_domain;
    // }
    cout << std::string(1, '\t') << "domain:\n";
    for(auto [dim, lw, up] : sh::views::zip(arr.domain.lower(), arr.domain.upper()) |
                                 std::views::take(geom_->n_dim) | sh::views::enumerate)
    {
        cout << std::string(2, '\t') << sh::dir(dim) << " : " << lw << ' ' << up << '\n';
    }

    double dissipate = conf["dissipate"].value_or(0);

    double spring_const, bulk_sound_sq;
    if(conf["set_both"].value_or(false))
    {
        bulk_sound_sq = conf["bulk_sound_sq"];
        spring_const = conf["spring_const"];
    }
    else
    {
        spring_const = std::pow(ptcl_buf.size() * 1e-4, 2.0 / geom_->n_dim - 2);
        conf["spring_const"] = spring_const;
        double spr_sph_ratio;
        if(conf["set_alpha"].value_or(false))
        {
            spr_sph_ratio = conf["spring_to_sph_ratio"];
        }
        else
        {
            auto calc_spr_sph_ratio = [](double poisson)
            {
                auto [a, b, c] = std::make_tuple(-0.11658747, 1.8306559, 0.37682884);
                return -1 / b * (poisson - c - a) / (poisson - c + a);
            };
            spr_sph_ratio = conf["spring_to_sph_ratio"] = calc_spr_sph_ratio(conf["poisson"]);
        }
        auto ptcl_mass = util::mean(ptcl_buf, &particles::full_ptcl::mass);
        bulk_sound_sq = 8e2 / spr_sph_ratio;
        bulk_sound_sq *= std::pow(arr.n_ptcl * 1e-4, 1);
        bulk_sound_sq *= spring_const * util::sq(arr.spacing) / ptcl_mass;
        conf["bulk_sound_sq"] = bulk_sound_sq;
    }

#if 0
    // adjust the magnitude of spring constant
    double stretch = arr.spacing * conf["spring_const_stretch_factor"].value_or(1);
    for(auto& [pow, att, rep] : spring_const)
    {
        att *= std::pow(stretch, -pow + 1);
        rep *= std::pow(stretch, -pow + 1);
    }
    conf["spring_const"] = spring_const;
#endif

#if 0
    // CFL time step criterion
    auto ptcl_mass = util::mean(ptcl_buf, &particles::full_ptcl::mass);
    spring_const_rest_length_ratio_ = conf["spring_const_rest_length_ratio"].value_or(1);
    time_step_ = std::numeric_limits<double>::max();
    auto snd_spd_max = std::numeric_limits<double>::lowest();
    for(auto dist : sh::views::linspace(0.0, arr.spacing * spring_length_ratio_, 100))
    {
        auto factor = (spring_const_rest_length_ratio_ - 1) / (spring_length_ratio_ - 1) *
                            (dist / arr.spacing - 1) +
                        1;
        if(factor < 0) factor = 0;
        double snd_spd_tmp{};
        for(auto [pow, att, rep] : spring_const)
        {
            auto spr_const = std::max(att, rep) * factor;
            snd_spd_tmp += std::sqrt(spr_const / ptcl_mass) * std::pow(dist, pow);
        }
        snd_spd_max = std::max(snd_spd_max, snd_spd_tmp);
        auto time_step_tmp = dist / snd_spd_tmp * CFL_factor_;
        time_step_ = std::min(time_step_tmp, time_step_);
    }
    conf["sound_speed_estimated"] = snd_spd_max;
    conf["time_step"] = time_step_;
#endif

    // Prevents some particles entering the interaction and others not, due to rounding errors,
    // when particles are evenly distributed.
    auto factor = 1.001;

    auto init = [&](particles::full_ptcl& p)
    {
        p.owner_id = owner_id;
        p.dissipate = dissipate;
        p.spring_const = spring_const;
        p.bulk_sound_sq = bulk_sound_sq;
        p.smth = factor * arr.spacing;
        p.r_search_for_spring = p.smth * field_.lock()->kernel_support_ratio_;
        p.r_search_for_dens = p.smth * field_.lock()->kernel_support_ratio_;
    };
    util::for_each(ptcl_buf, init);

    if(auto shift = conf.optional<double>("random_shift"))
        util::make_random_shift(ptcl_buf, &particles::full_ptcl::pos, *shift * arr.spacing);

    auto resister = [this](const particles::full_ptcl& p) { ptcl_range_->addOneParticle(p); };
    util::for_each(ptcl_buf, resister);

    init_dens_ = conf["dens"].value_or(util::mean(ptcl_buf, &particles::full_ptcl::dens));
    init_domain_ = arr.domain;
    init_mean_spacing_ = arr.spacing;
    id_ = owner_id;
}

void
material::make_basic_subsets(sh::syntax& conf)
{
    double side_width_ratio = conf["side_width_ratio"].value_or(1.5);

    auto side_width = vector(1) * init_mean_spacing_ * side_width_ratio;
    auto shrunk = init_domain_.make_extended(-side_width);

    make_subset("envelop", [&](const auto& p) { return !shrunk.contains(p.pos); });
    make_subset("lower.x", [&](const auto& p) { return p.pos[0] < shrunk.lower()[0]; });
    make_subset("upper.x", [&](const auto& p) { return shrunk.upper()[0] < p.pos[0]; });
    if(1 < geom_->n_dim)
    {
        make_subset("lower.y", [&](const auto& p) { return p.pos[1] < shrunk.lower()[1]; });
        make_subset("upper.y", [&](const auto& p) { return shrunk.upper()[1] < p.pos[1]; });
    }
    if(2 < geom_->n_dim)
    {
        make_subset("lower.z", [&](const auto& p) { return p.pos[2] < shrunk.lower()[2]; });
        make_subset("upper.z", [&](const auto& p) { return shrunk.upper()[2] < p.pos[2]; });
    }

    for(const auto& subset : subset_list_)
        if(subset.second->empty())
            throw sh::runtime_error{} << "The subset '" << subset.first << "' has no particle.";
}

material::material(sh::syntax& conf,
                   std::shared_ptr<particle_range<particles::full_ptcl>> ptcl_range,
                   std::shared_ptr<geometry> geo, std::weak_ptr<field> field)
    : recursive_particle_view{ std::move(ptcl_range), std::move(geo),
                               static_cast<filter_type>(
                                   [id = (size_t)conf["index"]](const auto& p)
                                   { return p.owner_id == id; }),
                               sh::finite<std::size_t>{ std::views::iota(0UL, n_max_views) } },
      field_(std::move(field))
{
    disable_lazy_evaluation();
    initialize(conf);
    make_basic_subsets(conf);
    // cerr << "constructed material: " << owner_id << '\n' << "n_ptcl = " <<
    // _ptr_to_ptcl_sys.lock()->getNumberOfParticleGlobal() << '\n' << "box = " <<
    // domain_init_ << '\n';
}

std::size_t
material::index() const
{
    return id_;
}

orthotope
material::init_domain() const
{
    return init_domain_;
}

double
material::init_mean_spacing() const
{
    return init_mean_spacing_;
}

double
material::equalize_density() const
{
    auto tune_mass = [this](particles::full_ptcl& p)
    {
        p.mass *= init_dens_ / p.dens;
        p.dens = init_dens_;  // Eliminates the need to recalculate densities
    };
    for_each(tune_mass);
    return init_dens_;
}

orthotope
material::domain_by_edge_ptcls() const
{
    orthotope domain;
    for(auto i : std::views::iota(0UL, geom_->n_dim))
    {
        std::string ax{ sh::dir(i) };
        domain.lower()[i] = (*this)["lower." + ax].pos()[i];
        domain.upper()[i] = (*this)["upper." + ax].pos()[i];
    }
    return domain;
}

}  // namespace cfe