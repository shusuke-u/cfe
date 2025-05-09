#include <chrono>
#include <filesystem>

#include <cfe.h>

constexpr inline auto pi = M_PI;

// Make a linear perturbation of one wave i.e. k = 2 \pi / L
auto
make_perturbation(sh::syntax& conf, cfe::field& field, cfe::material& plane)
{
    double wave_num = conf["wave_num"];
    double ampl = conf["ampl"] * plane.init_mean_spacing();
    cfe::orthotope domain = conf["domain"];
    auto side = domain.diag().x();

    auto perturb = [&](cfe::particles::full_ptcl& p)
    {
        p.pos.y() += ampl * std::sin(wave_num * p.pos.x() / side);
        p.pos.z() += ampl * std::sin(wave_num * p.pos.x() / side + 0.5 * pi);
    };
    plane.for_each(perturb);

    cfe::vector_t<size_t> pos_max_id;
    for(auto [i, id] : pos_max_id | sh::views::enumerate)
    {
        auto diff = [i](const auto& p) { return p.pos_diff()[i]; };
        id = plane.max_element(diff).id;
    }

    field.compute_density();

    return std::tuple{ wave_num, ampl, pos_max_id };
}

struct period_measure
{
    double prev_val;
    double prev_prev_val;
    size_t count;
    double prev_time;
    struct period
    {
        double start_time, period;
    };
    std::vector<period> periods;

    auto
    measure(double curr_val, double curr_time)
    {
        if(count++ < 2) return false;

        bool ret = false;
        if(prev_prev_val < prev_val && prev_val > curr_val)
        {
            ret = true;
            if(periods.empty())
            {
                periods.emplace_back(curr_time, 0);
            }
            else
            {
                periods.back().period = curr_time - periods.back().start_time;
                periods.emplace_back(prev_time, curr_time);
            }
        }
        prev_prev_val = prev_val;
        prev_val = curr_val;
        prev_time = curr_time;
        return ret;
    }

    auto
    write(std::ostream& os, double wave_num)
    {
        os << "start_time\tperiod\n";
        using namespace std::views;
        for(auto [time, period] : periods | reverse | drop(1) | reverse)
            os << time << '\t' << period << '\n';

        auto period_view = periods | reverse | drop(1) | transform(&period::period);
        auto mean_period =
            std::reduce(period_view.begin(), period_view.end(), 0.0) / period_view.size();
        auto angl_vel = 2 * pi / mean_period;
        auto sound_speed = angl_vel / wave_num;

        os << "n_periods:\t" << period_view.size() << '\n';
        os << "average_period:\t" << mean_period << '\n';
        os << "angular_vel:\t" << angl_vel << '\n';
        os << "wave_num:\t" << wave_num << '\n';
        os << "sound_speed:\t" << sound_speed << '\n';
        return std::tuple(mean_period, angl_vel, sound_speed);
    }
};

int
main(int argc, char* argv[])
{
    cfe::initialize(argc, argv);
    sh::syntax conf(2 <= argc ? argv[1] : "circ_linear.py");

    cfe::field field{ conf };
    auto& plane = field.make_material(conf);

    field.finish_setup();

    auto [wave_num, ampl, sample_id] = make_perturbation(conf, field, plane);

    // field.write(std::ofstream{ "debug/dump/1.csv" }, cfe::distr);

    size_t step_out = conf["output_interval"];
    auto n_period = conf["n_period"];
    auto step_end = conf["step_end"];
    cfe::vector_t<period_measure> maximum_measures;

    auto cwd = std::filesystem::current_path();

    for(auto step : std::views::iota(0UL))
    {
        // auto v_max_now = field.max( [ampl_norm]( const auto& p ) { return p.vel * ampl_norm;
        // } ); v_max.max = std::max( v_max.max, v_max_now ); v_max.min = std::min( v_max.min,
        // v_max_now );

        cfe::vector_t<double> sample_disp;
        for(auto [measure, id, i] :
            sh::views::zip(maximum_measures, sample_id, std::views::iota(0)))
        {
            auto disp = plane.ptcl(id).pos_diff()[i];
            measure.measure(disp, field.time());
        }

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
        if(step_end.has_value() && (size_t)step_end < step) break;
        if(std::ranges::all_of(maximum_measures,
                               [](const auto& m) { return 10 <= m.periods.size(); }))
            break;

        field.next_step();
    }

    if(cfe::in_root_proc)
    {
        for(auto [i, measure] : maximum_measures | sh::views::enumerate)
        {
            if(i == 0) continue;
            std::string dir{ sh::dir(i) };
            std::cout << "wave in  " << dir << '\n';
            auto [period, angl_vel, sound] = measure.write(std::cout, wave_num);
            std::cout << '\n';
            conf["wave." + dir + ".period"] = period;
            conf["wave." + dir + ".angular_vel"] = angl_vel;
            conf["wave." + dir + ".sound_speed"] = angl_vel / wave_num;
        }
        conf.write("result.yml", sh::syntax::yml);
        conf.write("result.tsv", sh::syntax::tsv, true);
    }

    cfe::finalize();
}