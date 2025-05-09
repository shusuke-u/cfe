#include <chrono>
#include <filesystem>

#include <cfe.h>

constexpr inline auto pi = M_PI;

// Make a linear perturbation of one wave i.e. k = 2 \pi / L
std::tuple<double, cfe::vector, cfe::vector, double>
make_perturbation(sh::syntax& conf, cfe::field& field, cfe::material& plane)
{
    cfe::vector wave_num = conf["wave_num"];

    auto sound_speed = field.sound_speed();
    conf["estimated_sound_speed"] = sound_speed;

    auto angl_vel = wave_num.length() * sound_speed;
    auto period = 2 * pi / angl_vel;
    auto ampl = (cfe::vector)conf["amplitude"] * plane.init_mean_spacing();

    cfe::orthotope domain = conf["domain"];
    auto side = domain.diag() * wave_num.make_normalized();

    bool perturb_pos = conf["perturb_pos"];
    bool perturb_vel = conf["perturb_vel"];

    auto perturb = [&](cfe::particles::full_ptcl& p)
    {
        if(perturb_pos) p.pos += std::sin(wave_num * p.pos / side) * ampl;
        if(perturb_vel) p.vel += std::cos(wave_num * p.pos / side) * ampl * -angl_vel;
    };
    plane.for_each(perturb);

    auto ampl_norm = ampl.make_normalized();
    auto get_diff = [ampl_norm](const auto& p) { return p.pos_diff() * ampl_norm; };
    auto pos_max_id = plane.max_element(get_diff).id;

    field.compute_density();

    conf["estimated_period"] = period;
    conf["estimated_angular_vel"] = angl_vel;

    return { period, wave_num, ampl.make_normalized(), pos_max_id };
}

int
main(int argc, char* argv[])
{
    cfe::initialize(argc, argv);
    sh::syntax conf(2 <= argc ? argv[1] : "linear.py");

    cfe::field field{ conf };
    auto& plane = field.make_material(conf);
    field.finish_setup();
    field.write("debug/0.tsv");

    auto [period, wave_num, ampl_norm, sample_id] = make_perturbation(conf, field, plane);

    size_t step_out = conf["output_interval"];
    auto n_period = conf.optional<size_t>("n_period");
    auto step_end = conf.optional<size_t>("step_end");

    // sh::min_max<double> v_max{ +1e5, -1e5 };
    double dx_sam_prv{}, dx_sam_prv_prv{};
    double time_prv{}, time_std_prv{};
    std::vector<std::pair<double, double>> periods;

    auto cwd = std::filesystem::current_path();

    for(auto step : std::views::iota(0UL))
    {
        auto dx_sam = plane.ptcl(sample_id).pos_diff() * ampl_norm;
        if(dx_sam_prv_prv < dx_sam_prv && dx_sam < dx_sam_prv)
        {
            if(time_std_prv != 0) periods.emplace_back(time_prv, time_prv - time_std_prv);
            if(n_period && periods.size() == *n_period) break;
            time_std_prv = time_prv;
        }
        dx_sam_prv_prv = dx_sam_prv;
        dx_sam_prv = dx_sam;
        time_prv = field.time();

        if(step % step_out == 0)
        {
            if(!std::filesystem::exists(cwd)) PS::Abort();

            cfe::cout << "step:\t" << step << '\n';
            cfe::cout << "time:\t" << field.time() << '\n';
            cfe::cout << "total_eng:\t" << field.total_eng() << '\n';
            cfe::cout << "kinetic_eng:\t" << field.kinetic_eng() << '\n';
            cfe::cout << "total_mom:\t" << field.mom() << '\n';
            cfe::cout << "time_step:\t" << field.time_step() << '\n';
            cfe::cout << std::endl;

            plane.write("dump/" + std::to_string(step) + ".tsv", field.time());
        }
        if(step_end && *step_end < step) break;

        field.next_step();
    }

    if(cfe::in_root_proc)
    {
        std::cerr << "\ntime\tperiod\n";
        for(auto [t, p] : periods) std::cerr << t << '\t' << p << '\n';
        double period{};
        for(auto [t, p] : periods) period += p;
        period /= periods.size();
        auto angl_vel = 2 * pi / period;
        auto sound_speed = angl_vel / wave_num.length();

        std::cerr << '\n';
        std::cerr << "period\t" << period << '\n';
        std::cerr << "wave_num\t" << wave_num << '\n';
        std::cerr << "angular vel\t" << angl_vel << '\n';
        std::cerr << "sound speed\t" << sound_speed << '\n';

        conf["period"] = period;
        conf["angular_vel"] = angl_vel;
        conf["sound_speed"] = sound_speed;
        conf.write("result.yml", sh::syntax::yml);
        conf.write("result.tsv", sh::syntax::tsv, true);
    }

    cfe::finalize();
}