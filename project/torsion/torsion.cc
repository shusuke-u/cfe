#include <chrono>
#include <filesystem>

#include <cfe.h>

int
main(int argc, char* argv[])
{
    cfe::initialize(argc, argv);
    sh::syntax conf(2 <= argc ? argv[1] : "torsion.py");

    conf["n_ptcl"] *= 4 / M_PI;
    double angle = conf["angle"];
    double radius = conf["radius"];
    double height = conf["height"];
    conf["domain"] = cfe::orthotope{ { -radius * 1.1, -radius * 1.1, 0 },
                                     { radius * 1.1, radius * 1.1, height } };
    cfe::vector center = { 0, 0, 0 };

    cfe::field field(conf);
    auto& lod = field.make_material(conf);
    field.clear_vel();  // set `vel_init` zero

    auto out_of_circle = [&](const cfe::particles::full_ptcl& p)
    { return radius < (p.pos - center).set(2, 0).length(); };
    lod.make_subset(out_of_circle).remove();

    field.finish_setup();

    auto& lower_side = lod["lower.z"];
    auto& upper_side = lod["upper.z"];
    lower_side.fix(0, 1, 2);
    upper_side.fix(0, 1, 2);

    auto torsion = [&](cfe::particles::full_ptcl& p)
    { p.pos = (p.pos - center).rotate(angle * p.pos.z() / height) + center; };
    lod.for_each(torsion);

    auto msa_converge_bound = conf.optional<double>("msa_converge_bound");
    auto torque_converge_bound = conf.optional<double>("torque_converge_bound");
    auto clear_vel_int = conf.optional<size_t>("clear_vel_interval");
    size_t clear_vel_low = conf["clear_vel_lower"].value_or(0);

    size_t step_out = conf["step_out"];

    auto cwd = std::filesystem::current_path();

    for(auto step : std::views::iota(0UL))
    {
        if(step % step_out == 0)
        {
            if(!std::filesystem::exists(cwd)) cfe::abort();

            // while the forces of the side particles are cleared
            double msa{};
            static double msa_prev;
            msa_prev = msa;
            conf["msa"] = msa = field.mean_square_acc();

            auto get_torque = [center](const cfe::particles::full_ptcl& p)
            { return p.mass / p.dens * (p.pos - center).set_z(0) ^ p.force; };
            double torque{};
            static double torque_prev;
            torque_prev = torque;
            torque = lod.mean(get_torque).z();
            conf["torque"] = torque;

            if(cfe::in_root_proc)
            {
                conf.write(std::ofstream("result.yml"), sh::syntax::yml);
                conf.write(std::ofstream("result.tsv"), sh::syntax::tsv, true);
            }
            cfe::cout << "step:\t" << step << '\n';
            cfe::cout << "time:\t" << field.time() << '\n';
            cfe::cout << "msa:\t" << msa << '\n';
            cfe::cout << "E_tot:\t" << field.total_eng() << '\n';
            cfe::cout << "E_kin:\t" << field.kinetic_eng() << '\n';
            cfe::cout << "ang_mom:\t" << field.angular_mom() << '\n';
            cfe::cout << "torque:\t" << torque << '\n';
            cfe::cout << '\n';

            lod.write("dump/" + std::to_string(step) + ".tsv", field.time());

            auto msa_converged = !msa_converge_bound ||
                                 std::abs((msa - msa_prev) / msa_prev) < *msa_converge_bound;
            auto torque_converged =
                !torque_converge_bound ||
                std::abs((torque - torque_prev) / torque_prev) < *torque_converge_bound;
            if(step != 0 && msa_converged && torque_converged) break;
        }

        if(clear_vel_int && step % *clear_vel_int == 0 && clear_vel_low < step) lod.clear_vel();

        field.next_step();
    }

    cfe::finalize();
}