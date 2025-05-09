#include <eigen-3.4.0/Eigen/Eigen>

#include <elasticity_simulator.h>

int
main(int argc, char* argv[])
{
    cfe::initialize(argc, argv);

    sh::syntax conf{ 2 <= argc ? argv[1] : "godunov.py" };

    cfe::field field{ conf };
    auto plane = field.make_material(conf);

    field.finish_setup();

    field.write(std::ofstream{ "./dump/0.csv" }, cfe::distr);

    const auto size = plane.n_ptcl();

    Eigen::MatrixXd m(size, size);

    cfe::cerr << "size = " << size << '\n';

    std::size_t n_dim = conf["n_dim"];

    auto i = 0UL;

    plane.for_each(
        [&](const auto& p_i)
        {
            auto j = 0UL;

            plane.for_each(
                [&](const auto& p_j)
                {
                    m(i, j) = cfe::kernel::gaussian(p_i.pos - p_j.pos, p_i.smth, n_dim);

                    ++j;
                });

            ++i;
        });

    cfe::cout << "det( M ) = " << m.determinant() << '\n';

    // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver( m );

    // if( solver.info() != Eigen::Success ) PS::Abort();

    // std::cout << "The eigenvalues of A are:\n" << solver.eigenvalues() << std::endl;
    // std::cout << "Here's a matrix whose columns are eigenvectors of A \n"
    //           << "corresponding to these eigenvalues:\n"
    //           << solver.eigenvectors() << std::endl;

    cfe::finalize();
}
