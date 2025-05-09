#include <chrono>
#include <filesystem>

#include <elasticity_simulator.h>

int
main(int argc, char* argv[])
{
    cfe::initialize(argc, argv);

    sh::syntax conf{ 2 <= argc ? argv[1] : "speed.py" };

    cfe::field field{ conf };
    auto whole = field.make_material(conf);

    cfe::clear_dir("./dump");
    cfe::clear_dir("./debug/distr");

    field.finish_setup();

    auto box = whole.box_info().make_scaled(1. / 3);
    auto& plane = whole.make_subset("plane", [box](const auto& p) { return box.contains(p.pos); });

    box = plane.box_info().make_scaled(1 - 0.05);
    plane.make_subset("low.x", [box](const auto& p) { return p.pos()[0] < box.min()[0]; });
    plane.make_subset("high.x", [box](const auto& p) { return box.top()[0] < p.pos()[0]; });
    plane["low.x"].drift(whole.spacing() * conf["shift"]);

    field.compute_density();
    field.write(std::ofstream{ "debug/distr/0.csv" }, cfe::distr);

    std::ofstream reach{ "reach.csv" };
    reach << std::scientific << "time,vel.x,\n";

    auto step = conf.get("step", 0UL);
    auto time_out = 0.;
    const auto time_end = conf.get("time_end", -1.);
    const auto time_dif = time_end / conf.get("n_out", 0UL);

    for(; field.time() <= time_end; ++step)
    {
        reach << field.time() << '\t' << plane["high.x"].vel()[0] << '\n';

        if(1e-5 < plane["high.x"].vel().length())
        {
            const auto speed = (plane["high.x"].pos() - plane["low.x"].pos())[0] / field.time();
            cfe::cout << "speed = " << speed << '\n';
            conf.set("speed", speed);
            break;
        }

        if(time_out <= field.time())
        {
            time_out += time_dif;

            cfe::cout << "step = " << step << '\n'
                      << "time = " << field.time() << '\n'
                      << "msa  = " << field.mean_square_acc(true) << '\n'
                      << "kin  = " << field.kinetic_energy() << '\n'
                      << '-' << '\n'
                      << std::flush;

            whole.write(std::ofstream{ "./dump/" + std::to_string(step) + ".csv" }, field.time());
        }

        field.next_step();
    }

    conf.write(std::ofstream{ "./conf_res.py" });

    cfe::finalize();
}
