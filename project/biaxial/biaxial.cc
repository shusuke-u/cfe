#include <chrono>
#include <filesystem>

#include <cfe.h>

int
main(int argc, char* argv[])
{
    cfe::initialize(argc, argv);
    sh::syntax conf(2 <= argc ? argv[1] : "tensile.py");

    cfe::field field(conf);
    auto material = field.make_material(conf);
    field.clear_vel();  // set `vel_init` zero
    field.finish_setup();

    auto& lower_side = material["lower.z"];
    auto& upper_side = material["upper.z"];
    if(lower_side.empty()) throw sh::runtime_error{} << "The 'lower_side' has no particle.";
    if(upper_side.empty()) throw sh::runtime_error{} << "The 'upper_side' has no particle.";

    lower_side.fix(0, 1, 2);
    upper_side.fix(0, 1, 2);

    double converge_bound = conf["poisson_converge_bound"];
    // double tensile_begin = conf["tensile_begin"];
    // double tensile_end = conf["tensile_end"];
    // cfe::vector tensile_shift;
    // tensile_shift[tensile_axis] = conf["tensile_length_ratio"] *
    //                               material.domain_by_edge_ptcls().diag()[tensile_axis] /
    //                               (tensile_end - tensile_begin);
    auto msa_converge_bound = conf.optional<double>("msa_converge_bound");
    auto clear_vel_int = conf.optional<size_t>("clear_vel_interval");
    size_t clear_vel_low = conf["clear_vel_lower"].value_or(0);

    auto init_side = material.domain_by_edge_ptcls().diag();

    for(auto i : cfe::vector::indices())
        field.scale(i == tensile_axis ? 1 + conf["tensile_length_ratio"]
                                      : 1 - conf["tensile_length_ratio"] *
                                                conf["poisson_estimated"].value_or(0.3),
                    i);

    size_t step_out = conf["step_out"];

    auto cwd = std::filesystem::current_path();

    for(auto step : std::views::iota(0UL))
    {
        if(step % step_out == 0)
        {
            if(not std::filesystem::exists(cwd)) cfe::abort();

            // static cfe::vector init_side;
            // static auto side_uninit = true;
            // if(side_uninit and tensile_begin < step)
            // {
            //     init_side = material.domain_by_edge_ptcls().diag();
            //     side_uninit = false;
            // }

            // while the forces of the side particles are cleared
            double msa{};
            static double msa_prev{};
            msa_prev = msa;
            conf["msa"] = msa = field.mean_square_acc();

            static cfe::vector poisson_prev;
            auto domain = material.domain_by_edge_ptcls();
            auto side = domain.diag();
            auto strain = (cfe::vector)sh::call([](auto a, auto b) { return a / b - 1; }, side,
                                                init_side);
            auto poisson = -strain / strain[tensile_axis];

            field.compute_force();
            auto stress = lower_side.mean(&cfe::particles::full_ptcl::force);
            auto young = stress[tensile_axis] / strain[tensile_axis];

            if(cfe::in_root_proc)
            {
                conf["side_init"] = init_side;
                conf["side"] = side;
                conf["domain"] = domain;
                conf["strain"] = strain;
                conf["stress"] = stress;
                conf["poisson"] = poisson;
                conf["young"] = young;
                conf.write(std::ofstream("result.yml"), sh::syntax::yml);
                conf.write(std::ofstream("result.tsv"), sh::syntax::tsv, true);
            }
            cfe::cout << "step:\t" << step << '\n';
            cfe::cout << "time:\t" << field.time() << '\n';
            cfe::cout << "msa:\t" << msa << '\n';
            cfe::cout << "E_tot:\t" << field.total_eng() << '\n';
            cfe::cout << "E_kin:\t" << field.kinetic_eng() << '\n';
            cfe::cout << "ang_mom:\t" << field.angular_mom() << '\n';
            cfe::cout << "Poisson:\t" << poisson << '\n';
            cfe::cout << "Young:\t" << young << '\n';
            cfe::cout << '\n';

            material.write(std::ofstream{ "dump/" + std::to_string(step) + ".tsv" },
                           field.time());

            auto poisson_converged = true;
            for(auto i : cfe::vector::indices())
                if(i != tensile_axis and
                   converge_bound < std::abs((poisson[i] - poisson_prev[i]) / poisson_prev[i]))
                    poisson_converged = false;
            poisson_prev = poisson;
            // auto stretch_ended = tensile_end < step;
            auto msa_converged = not msa_converge_bound or
                                 std::abs((msa - msa_prev) / msa_prev) < *msa_converge_bound;
            // if(not side_uninit and poisson_converged and stretch_ended and msa_converged)
            // break;
            if(poisson_converged and msa_converged) break;
        }

        // if(tensile_begin < step and step < tensile_end) upper_side.drift(tensile_shift);

        if(clear_vel_int and step % *clear_vel_int == 0 and clear_vel_low < step)
            material.clear_vel();

        field.next_step();
    }

    cfe::finalize();
}