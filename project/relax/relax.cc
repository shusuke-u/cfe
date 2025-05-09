#include <filesystem>

#include <cfe.h>

int
main(int argc, char* argv[])
{
    cfe::initialize(argc, argv);
    sh::syntax conf(2 <= argc ? argv[1] : "relax.py");

    cfe::field field(conf);
    auto& material = field.make_material(conf);
    field.finish_setup();

    // // Thin out particles
    // auto random_sample = [](const cfe::particles::full_ptcl& p)
    // {
    //     static std::mt19937 gen(std::random_device{}());
    //     static std::uniform_int_distribution<> dis(0, 9);
    //     return dis(gen) != 0;
    // };
    // material.make_subset(random_sample).remove();
    // material.write("final.tsv", 0);
    // PS::Comm::barrier();
    // cfe::abort();

    double msa_max = conf["max_msa_bound"];
    double rep_factor = conf["repulsion_factor"];

    std::optional<std::size_t> clear_vel_int = conf["clear_vel_interval"];
    size_t clear_vel_low = conf["clear_vel_lower"];

    size_t step_out = conf["step_out"];
    auto step_end = conf.optional<size_t>("step_end");

    auto cwd = std::filesystem::current_path();

    for(auto step : std::views::iota(0UL))
    {
        if(step % step_out == 0)
        {
            if(!std::filesystem::exists(cwd)) PS::Abort();

            auto msa = field.mean_square_acc();
            conf["msa"] = msa;

            cfe::cout << "step  : " << step << '\n';
            cfe::cout << "time  : " << field.time() << '\n';
            cfe::cout << "msa   : " << msa << '\n';
            cfe::cout << "E_tot : " << field.total_eng() << '\n';
            cfe::cout << "E_kin : " << field.kinetic_eng() << '\n';
            cfe::cout << "dt    : " << field.time_step() << '\n';
            cfe::cout << '\n';

            static size_t i{};
            field.write("dump/" + std::to_string(i++) + ".tsv");

            if(0 < step && msa < msa_max) break;
        }

        if(step_end && *step_end <= step) break;

        if(clear_vel_int && step % *clear_vel_int == 0 && clear_vel_low < step)
            material.clear_vel();

        field.next_step(cfe::repulsion, rep_factor);
    }

    conf.write("result.yml", sh::syntax::yml);
    conf.write("result.tsv", sh::syntax::tsv, true);

    field.clear_vel();
    field.write("final.tsv");

    cfe::finalize();
}