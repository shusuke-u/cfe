#include <cmath>
#include <chrono>
#include <filesystem>

#include <cfe.h>

int
main(int argc, char* argv[])
{
    cfe::initialize(argc, argv);

    sh::syntax conf(2 <= argc ? argv[1] : "plate.py");

    cfe::field field(conf);
    auto plate = field.make_material(conf);
    field.finish_setup();

    field.write_springs(std::ofstream("debug/springs.tsv"));

    double init_vel = conf["init_vel"];
    std::size_t step_out = conf["step_out"];
    std::size_t step_end = conf["step_end"];

    auto k = 1.875;
    auto L = 1;
    auto M = std::sin(k * L) + std::sinh(k * L);
    auto N = std::cos(k * L) + std::cosh(k * L);
    auto Q = 2 * (std::cos(k * L) * std::sinh(k * L) - std::sin(k * L) * std::cosh(k * L));
    auto init_domain = plate.domain_by_edge_ptcls();
    auto set_vel = [&](auto& p)
    {
        auto x = p.pos[0] - init_domain.lower()[0];
        p.vel[1] = init_vel / Q *
                   (M * (std::cos(k * x) - std::cosh(k * x)) -
                    N * (std::sin(k * x) - std::sinh(k * x)));
    };
    plate.for_each(set_vel);
    plate["lower.x"].clear_vel();
    plate["lower.x"].fix(0, 1, 2);

    for(auto step : std::views::iota(0UL))
    {
        if(step % step_out == 0)
        {
            cfe::cout << "step:\t" << step << '\n';
            cfe::cout << "time:\t" << field.time() << '\n';
            cfe::cout << "msa:\t" << field.mean_square_acc() << '\n';
            cfe::cout << "E_tot:\t" << field.total_eng() << '\n';
            cfe::cout << "E_kin:\t" << field.kinetic_eng() << '\n';
            cfe::cout << "ang_mom:\t" << field.angular_mom() << '\n';
            cfe::cout << '\n';

            plate.write(std::ofstream("dump/" + std::to_string(step) + ".tsv"), field.time());
        }

        if(step == step_end) break;

        field.next_step();
    }

    conf.write(std::ofstream("result.yml"), sh::syntax::yml);
    conf.write(std::ofstream("result.tsv"), sh::syntax::tsv, true);

    cfe::finalize();
    return 0;
}