#include <chrono>
#include <filesystem>

#include <cfe.h>

int
main(int argc, char* argv[])
{
    cfe::initialize(argc, argv);
    sh::syntax conf(2 <= argc ? argv[1] : "tensile.py");

    auto field = std::make_shared<cfe::field>(conf);
    auto material = field->make_material(conf);
    field->clear_vel();  // set `vel_init` zero
    field->finish_setup();

    std::string tensile_axis_str = conf["direction"];
    auto tensile_axis = sh::dir(tensile_axis_str);
    auto& lower_side = material->subset("lower." + tensile_axis_str);
    auto& upper_side = material->subset("upper." + tensile_axis_str);

    lower_side.fix(tensile_axis);
    upper_side.fix(tensile_axis);

    // double tensile_begin = conf["tensile_begin"];
    // double tensile_end = conf["tensile_end"];
    // cfe::vector tensile_shift;
    // tensile_shift[tensile_axis] = conf["tensile_length_ratio"] *
    //                               material->domain_by_edge_ptcls().diag()[tensile_axis] /
    //                               (tensile_end - tensile_begin);
    auto msa_converge_bound = conf["msa_converge_bound"];
    auto clear_vel_int = conf["clear_vel_interval"];
    auto min_step_bound = conf["min_step_bound"];
    size_t clear_vel_low = conf["clear_vel_lower"].value_or(0);

    auto init_side = material->domain_by_edge_ptcls().diag();

    auto [a, b, c] = std::tuple(-0.11658747, 1.8306559, 0.37682884);
    auto alpha = conf["spring_to_sph_ratio"];
    auto poisson_estimated = a * (b * alpha - 1) / (b * alpha + 1) + c;
    conf["poisson_estimated"] = poisson_estimated;
    for(auto i : cfe::vector::indices())
        field->scale(i == tensile_axis ? 1 + conf["tensile_length_ratio"]
                                       : 1 - conf["tensile_length_ratio"] * poisson_estimated,
                     i);

    size_t step_out = conf["step_out"];

    auto cwd = std::filesystem::current_path();

    for(auto step : std::views::iota(0UL))
    {
        if(step % step_out == 0)
        {
            if(!std::filesystem::exists(cwd)) cfe::abort();

            // static cfe::vector init_side;
            // static auto side_uninit = true;
            // if(side_uninit && tensile_begin < step)
            // {
            //     init_side = material->domain_by_edge_ptcls().diag();
            //     side_uninit = false;
            // }

            // while the forces of the side particles are cleared
            double msa{};
            static double msa_prev{};
            msa_prev = msa;
            conf["msa"] = msa = field->mean_square_acc();

            static cfe::vector poisson_prev;
            auto domain = material->domain_by_edge_ptcls();
            auto side = domain.diag();
            auto strain = (cfe::vector)sh::call([](auto a, auto b) { return a / b - 1; }, side,
                                                init_side);
            auto poisson = -strain / strain[tensile_axis];

            field->compute_force();
            auto stress = -upper_side.sum(&cfe::particles::full_ptcl::force,
                                          [](const auto& p) { return p.mass / p.dens; });
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
                conf.write("result.yml", sh::syntax::yml);
                conf.write("result.tsv", sh::syntax::tsv, true);
            }
            cfe::cout << "step:\t" << step << '\n';
            cfe::cout << "time:\t" << field->time() << '\n';
            cfe::cout << "msa:\t" << msa << '\n';
            cfe::cout << "E_tot:\t" << field->total_eng() << '\n';
            cfe::cout << "E_kin:\t" << field->kinetic_eng() << '\n';
            cfe::cout << "ang_mom:\t" << field->angular_mom() << '\n';
            cfe::cout << "Poisson:\t" << poisson << '\n';
            cfe::cout << "Young:\t" << young << '\n';
            cfe::cout << '\n';

            material->write("dump/" + std::to_string(step) + ".tsv", field->time());

            auto poisson_converged = true;
            for(auto [p, pp] : sh::views::zip(poisson, poisson_prev))
                if(p != pp && (double)conf["poisson_converge_bound"] < std::abs((p - pp) / pp))
                    poisson_converged = false;
            poisson_prev = poisson;

            auto msa_converged =
                msa_converge_bound.empty() ||
                std::abs((msa - msa_prev) / msa_prev) < (double)msa_converge_bound;
            auto step_passed = min_step_bound.empty() || (double)min_step_bound < step;
            if(poisson_converged && msa_converged && step_passed) break;
        }

        if(clear_vel_int.has_value() && step % (size_t)clear_vel_int == 0 &&
           clear_vel_low < step)
            material->clear_vel();

        field->next_step();
    }

    cfe::finalize();
}