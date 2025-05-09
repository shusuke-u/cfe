#include <cmath>
#include <array>
#include <chrono>
#include <numbers>
#include <fstream>
#include <sstream>
#include <iostream>
#include <optional>
#include <iomanip>

#include <sh/syntax.h>
#include <sh/ranges.h>
#include <sh/min_max.h>

struct spring
{
    double k_rep, k_att;

    spring( const sh::syntax& init ) : k_rep{ init["rep"] }, k_att{ init["att"] }
    {
    }

    double
    operator()( double grad_f ) const
    {
        return grad_f < 0 ? k_rep : k_att;
    }
};

struct parameter
{
    parameter( const sh::syntax& init )
        : omega{ 0 },
          scaling{ init["scaling"] },
          eigen_begin{ init["eigen_begin"] },
          range{ (double)init["bound_lower"], (double)init["bound_upper"] },
          n_div_range{ init["bound_n_div"] },
          spr{ init }
    {
    }

    double omega;
    double scaling;
    double eigen_begin;
    sh::min_max<double> range;
    size_t n_div_range;
    spring spr;
};

struct variable
{
    double f, grad_f;

    variable() = default;

    variable( const sh::syntax& init ) : f{ init["f_init"] }, grad_f{ init["grad_f_init"] }
    {
    }
};

struct derivative
{
    double grad_f, grad_grad_f;

    void
    calc( const parameter& param, variable var )
    {
        this->grad_grad_f = -param.omega * param.omega / param.spr( var.grad_f ) * var.f;
        this->grad_f = var.grad_f;
    }
};

// Solving ODEs
std::stringstream
solve_ODEs( const parameter& param, variable& var, bool output = false )
{
    std::array<derivative, 4UL> coef;
    variable temp;

    const auto delta_x = param.range.bound() / param.n_div_range;

    std::stringstream buf;
    if( output ) buf << "x\tf\tdfdx\n";

    for( auto x : sh::views::linspace( param.range.min, param.n_div_range, delta_x ) )
    {
        if( output ) buf << x << '\t' << var.f << '\t' << var.grad_f << '\n';

        // 1st
        coef[0].calc( param, var );

        // 2nd
        temp.f = var.f + 0.5 * delta_x * coef[0].grad_f;
        temp.grad_f = var.grad_f + 0.5 * delta_x * coef[0].grad_grad_f;
        coef[1].calc( param, temp );

        // 3rd
        temp.f = var.f + 0.5 * delta_x * coef[1].grad_f;
        temp.grad_f = var.grad_f + 0.5 * delta_x * coef[1].grad_grad_f;
        coef[2].calc( param, temp );

        // 4th
        temp.f = var.f + delta_x * coef[2].grad_f;
        temp.grad_f = var.grad_f + delta_x * coef[2].grad_grad_f;
        coef[3].calc( param, temp );

        var.f += ( coef[0].grad_f + 2 * coef[1].grad_f + 2 * coef[2].grad_f + coef[3].grad_f ) * delta_x / 6;
        var.grad_f +=
            ( coef[0].grad_grad_f + 2 * coef[1].grad_grad_f + 2 * coef[2].grad_grad_f + coef[3].grad_grad_f ) * delta_x / 6;
    }

    var.f *= param.scaling;

    std::cout << "f: " << var.f << '\t' << "grad.f: " << var.grad_f << '\n';

    return buf;
}

// Shooting method
double
find_omega( parameter& param, variable& var )
{

    // where the sign of determinant change?
    auto eigen_value = param.eigen_begin;
    auto eigen_prev = eigen_value;
    const auto eigen_diff = param.eigen_begin * 1e-3;
    eigen_value -= eigen_diff;
    solve_ODEs( param, var );
    auto det_prev = var.f;
    auto det = var.f;

    for( auto count : std::views::iota( 0UL ) )
    {
        if( std::fabs( eigen_value / eigen_diff - 1.0 ) < 1e0 ) eigen_value -= 1.5 * eigen_diff;

        det_prev = det;
        param.omega = eigen_value;
        solve_ODEs( param, var );
        det = var.f;
        // printf("ω/(C_ensL) -> %e\t det -> %e\n", eigen_value, det);

        if( det_prev * det < 0 ) break;

        if( count == 1e3 ) throw std::runtime_error( "No eigen value was found." );

        eigen_prev = eigen_value;
        eigen_value -= eigen_diff;
    }

    // Solve eigen value by dichotomy
    double eigen_mid;
    auto eigen_low = eigen_prev;
    auto eigen_high = eigen_value;
    auto det_low = det_prev;
    auto det_mid = det_prev;

    for( auto count : std::views::iota( 0UL ) )
    {
        eigen_mid = 0.5 * ( eigen_low + eigen_high );
        param.omega = eigen_value;
        solve_ODEs( param, var );
        det_mid = var.f;

        if( det_low * det_mid > 0 )
        {
            eigen_low = eigen_mid;
        }
        else
        {
            eigen_high = eigen_mid;
        }

        if( std::fabs( eigen_mid ) < 1e-10 )
        {
            eigen_mid = 0;
            break;
        }

        if( std::fabs( det_mid ) < 1e-5 ) break;
        // printf("kH -> %e\t σ√(H/g) -> %e\t det -> %e\n",k,eigen_mid, det_mid);

        if( count == 1e4 ) throw std::runtime_error( "Cannot converged by dichotomy." );
    }
    // Terminal check
    std::cout << "omega/C_ens: " << std::sqrt( eigen_mid ) << '\n';

    return eigen_mid;
}

void
find_eigen_funcs( const sh::syntax& init )
{
    parameter param( init );
    variable var( init );

    const double omega = find_omega( param, var );

    // preparation for output
    std::stringstream path;
    const auto now = std::chrono::time_point_cast<std::chrono::seconds>( std::chrono::system_clock::now() );
    const auto now_c = std::chrono::system_clock::to_time_t( now );
    path << "./result/date:" << std::put_time( std::localtime( &now_c ), "%y%m%d%H%M%S" ) << "-k_rep:" << param.spr.k_rep
         << "-k_att:" << param.spr.k_att << ".tsv";
    std::ofstream file( path.str(), std::ios::out );
    if( not file.is_open() ) throw std::runtime_error( "Cannot open the file with the provided path: " + path.str() );

    const auto buf = solve_ODEs( param, var, true );
    file << buf.str();
}

int
main()
{
    try
    {
        sh::syntax init( "eigen_fn.py" );

        find_eigen_funcs( init );
    }
    catch( const std::exception& e )
    {
        std::cerr << "Error: " << e.what() << '\n';
    }
}