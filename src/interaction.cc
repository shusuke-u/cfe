#include <set>
#include <span>

#include <sh/ranges/concat_view.h>

#include <interaction.h>
#include <flags.h>

namespace cfe
{
namespace interaction
{
namespace functor
{

#if 0
make_spring::make_spring() : M_n_dim{}, M_min_spacing{}, M_spr_len_rat{},
M_const_rest_len_ratio{}
{
}

make_spring::make_spring(std::size_t n_dim, double min_spacing, double spring_length_ratio,
                         double spring_const_rest_length_ratio)
    : M_n_dim{ n_dim },
      M_min_spacing{ min_spacing },
      M_spr_len_rat{ spring_length_ratio },
      M_const_rest_len_ratio{ spring_const_rest_length_ratio }
{
}

void
make_spring::operator()(const particles::essential::make_spring* const arr_pi, const PS::S32
n_pi,
                        const particles::essential::make_spring* const arr_pj, const PS::S32
                        n_pj, interaction::result::make_spring* results) const
{
    for(auto& pi : std::span{ arr_pi, (std::size_t)n_pi })
    {
        auto& result = *results++;

        for(const auto& pj : std::span{ arr_pj, (std::size_t)n_pj })
        {
            if(pi.owner_id != pj.owner_id or pi.id == pj.id) continue;

            const auto rel_pos = pi.pos - pj.pos;
            const auto dist = rel_pos.length();

            if(M_min_spacing) std::cout << "(" << pi.id << "," << pj.id << ") r_search: " << pi.r_search << " dist: " << dist << std::endl;

            if(pi.r_search < dist)
            {
                // std::cout << "Out of bound.\tr_search: " << pi.r_search << "\tdist: " <<
                dist << std::endl; continue;
            }

            auto constant = pi.constant;

            // The particle spacing is not known at the first call, which the spacing is
            derived using
            // the spring list. The spring constants must be updated (multiplied) once which
            is in the second call. if(M_min_spacing)
            {
                // gaussian
                // coef.att.multiply( kernel::gaussian( rel_pos, *M_min_spacing *
                *M_spr_len_rat, M_n_dim ) );
                // coef.rep.multiply( kernel::gaussian( rel_pos, *M_min_spacing *
                *M_spr_len_rat, M_n_dim ) );

                // linear
                auto factor = (*M_const_rest_len_ratio - 1) / (pi.r_search - *M_min_spacing)
                * (dist - *M_min_spacing) + 1; if(factor < 0) factor = 0; for(auto& [pow,
                att, rep] : constant)
                {
                    att *= factor;
                    rep *= factor;
                }
            }

            auto succeed = result.springs.insert(pj.id, { .constant = constant, .length =
            dist });

            if(not succeed) throw exception{} << "Neighbor list is full.";

            // std::cerr << "id  :  " << pi.id << '\t' << pj.id << '\n'     //
            //           << "pos : " << pi.pos << pj.pos << '\n'            //
            //           << "coef: " << coef.rep << '\t' << coef.att << '\n'  //
            //           << "-" << '\n';

            if(M_min_spacing and true)
            {
                static std::ofstream f{ "./debug/spring.csv" };
                // auto& f = std::cout;
                f << std::setprecision(16) << std::scientific << std::showpos;
                // f << pi.id << '\t' << pj.id << '\t' << dist << '\t' << coef.att << '\t' <<
                coef.rep << '\n'; f << pi.id << '\t' << pj.id << '\t' << dist << '\n';
            }
        }
    }
}
#endif

void
spring_temp::operator()(const particles::essential::spring* const arr_pi, const PS::S32 n_pi,
                        const particles::essential::spring* const arr_pj, const PS::S32 n_pj,
                        interaction::result::spring* results) const
{
    auto view_pi = std::span{ arr_pi, (std::size_t)n_pi };
    auto view_pj = std::span{ arr_pj, (std::size_t)n_pj };

#ifdef SELF_IMPL_PERIODIC_BOUNDARY
    // for periodic boundary, the following is needed
    // and use view_pj_all

    // make unique view of pj
    std::vector view_pj_unique(view_pj.begin(), view_pj.end());
    std::ranges::sort(view_pj_unique, {}, &particles::essential::spring::id);
    auto last =
        std::ranges::unique(view_pj_unique, {}, &particles::essential::spring::id).begin();
    view_pj_unique.erase(last, view_pj_unique.end());

    // make image pjs
    std::vector<particles::essential::spring> view_pj_all;
    for(const auto& pj : view_pj_unique)
    {
        enum shift
        {
            minus,
            none,
            plus,
        };
        vector_t<std::vector<shift>> image_shifts;
        for(auto d : vector::indices())
        {
            // A real particle-j always exists
            image_shifts[d].push_back(none);

            if(not periodic_domain_[d]) continue;
            auto side = periodic_domain_[d].value();
            if(side.max < pj.pos[d] + pj.r_search) image_shifts[d].push_back(minus);
            if(pj.pos[d] - pj.r_search < side.min) image_shifts[d].push_back(plus);
        }
        for(const auto& shift : sh::views::direct_product(image_shifts))
        {
            assert(shift.size() == vector::size());
            auto& p = view_pj_all.emplace_back(pj);
            for(auto d : vector::indices())
            {
                if(not periodic_domain_[d]) continue;
                auto bound = periodic_domain_[d]->bound();
                if(shift[d] == minus) p.pos[d] -= bound;
                if(shift[d] == plus) p.pos[d] += bound;
            }
        }
    }
#else
    // assert(!periodic_domain_[0] && !periodic_domain_[1] && !periodic_domain_[2]);
#endif

    for(const auto& pi : view_pi)
    {
        auto& result = *results++;
        result.clear();

#ifdef SELF_IMPL_PERIODIC_BOUNDARY
        for(const auto& pj : view_pj_all) calc_force(pi, pj, result);
#else
        for(const auto& pj : view_pj) calc_force(pi, pj, result);
#endif
        // assert(result.n_ngb == 122);
    }
}

void
spring_temp::calc_force(const particles::essential::spring& pi,
                        const particles::essential::spring& pj,
                        interaction::result::spring& result) const
{
    auto diff = pi.pos - pj.pos;
    auto dist = diff.length();
    if(dist == 0) return;
    // if(pi.r_search < dist) return;

    auto diff_init = pi.init_pos - pj.init_pos;

    for(auto [domain, diff, diff_init] : sh::views::zip(periodic_domain_, diff, diff_init))
    {
        if(!domain) continue;
        auto bound = domain->bound();

        // Ensure that pj does not overlap between real and image
        if(bound - pi.r_search < std::abs(diff)) return;

        if(diff_init - pi.r_search < -bound)
        {
            diff_init += bound;
        }
        else if(bound < diff_init + pi.r_search)
        {
            diff_init -= bound;
        }
    }

    auto rest_len = diff_init.length();
    if(pi.r_search < rest_len) return;
    auto stretch = dist - rest_len;

#if 0
    auto type = pi.spring_const[0];
    vector force;
    //--------------------------------------------------------------------------------
    // 0. Normal
    //--------------------------------------------------------------------------------
    if(type == 0.)
    {
        force = -pi.spring_const[1] * stretch * diff / dist;
    }
    //--------------------------------------------------------------------------------
    // 1. Asymmetric, rest length dependent, and power-law spring force
    //--------------------------------------------------------------------------------
    else if(type == 1.)
    {
        double factor{};
        for(auto i = 1; i < pi.spring_const.size() - 2; i += 3)
        {
            auto pow = pi.spring_const[i];
            auto att = pi.spring_const[i + 1];
            auto rep = pi.spring_const[i + 2];
            auto f = (0 < stretch ? att : rep) * std::pow(stretch, pow);
            factor += f;
            if(static size_t j{}; j++ % 10000000 == 0) cerr << f << '\n';
        }
        auto rest_len_factor = (spr_const_rest_len_ratio_ - 1) /
                                   (pi.r_search - mean_ptcl_sep_) *
                                   (rest_len - mean_ptcl_sep_) +
                               1;
        // auto rest_len_factor = rest_len < spr_const_rest_len_ratio_ * mean_ptcl_sep_ ? 0 : 1;
        if(rest_len_factor < 0) rest_len_factor = 0;
        force = -factor * rest_len_factor * diff / dist;
    }
    //--------------------------------------------------------------------------------
    // 2. Intra-unit and other
    //--------------------------------------------------------------------------------
    else if(type == 2.)
    {
        auto factor = pi.unit_id and pj.unit_id and pi.unit_id.value() == pj.unit_id.value()
                          ? pi.spring_const[1]
                          : pi.spring_const[2];

        force = -factor * stretch * diff / dist;
    }
    //--------------------------------------------------------------------------------
    // 3.1 Spring and density term
    //--------------------------------------------------------------------------------
    else if(type == 3.1)
    {
        auto f0 = -pi.spring_const[1] * stretch;
        auto f1 = +pi.spring_const[2] * 0.5 *
                  (pi.dens + pj.dens - pi.init_dens - pj.init_dens) / dens_diff_std_ *
                  mean_ptcl_sep_;
        // auto f1 = +pi.spring_const[2] * ((pi.dens + pj.dens) * 0.5 - dens_std_);
        force = (f0 + f1) * diff / dist;
        if(static size_t i{}; omp_get_thread_num() == 0 and i++ % 1000000 == 0)
            cerr << f0 << '\t' << f1 << '\t' << (f0 + f1) << '\n';
    }
    //--------------------------------------------------------------------------------
    // 3.2 Spring and density term without init_dens
    //--------------------------------------------------------------------------------
    else if(type == 3.2)
    {
        auto f0 = -pi.spring_const[1] * stretch;
        auto f1 = +pi.spring_const[2] * 0.5 * (pi.dens + pj.dens) / dens_std_ * mean_ptcl_sep_;
        force = (f0 + f1) * diff / dist;
        if(static size_t i{}; omp_get_thread_num() == 0 and i++ % 1000000 == 0)
            cerr << f0 << '\t' << f1 << '\t' << (f0 + f1) << '\n';
    }
    //--------------------------------------------------------------------------------
    // 4.1 Density-dependent spring constant
    //--------------------------------------------------------------------------------
    else if(type == 4.1)
    {
        // auto dens_diff = 0.5 * (pi.dens + pj.dens) - dens_std_;
        auto dens_diff = 0.5 * (pi.dens + pj.dens - pi.init_dens - pj.init_dens);
        auto k0 = pi.spring_const[1];
        auto k1 = pi.spring_const[2] * std::pow(dens_diff / dens_diff_std_, 2);
        force -= (k0 + k1) * stretch * diff / dist;
        if(static size_t i{}; omp_get_thread_num() == 0 and i++ % 1000000 == 0)
            cerr << -k0 * stretch << '\t' << -k1 * stretch << '\n';
    }
    //--------------------------------------------------------------------------------
    // 4.2 Density-dependent spring constant
    //--------------------------------------------------------------------------------
    else if(type == 4.2)
    {
        auto dens = 0.5 * (pi.dens + pj.dens);
        auto k0 = pi.spring_const[1];
        auto k1 = pi.spring_const[2] * std::pow(dens / dens_std_, 2);
        force -= (k0 + k1) * stretch * diff / dist;
        if(static size_t i{}; omp_get_thread_num() == 0 and i++ % 1000000 == 0)
            cerr << -k0 * stretch << '\t' << -k1 * stretch << '\n';
    }
    //--------------------------------------------------------------------------------
    // 5. Spring + SPH force
    //--------------------------------------------------------------------------------
    else if(type == 5.)
    {
        // spring
        auto k = pi.spring_const[1] * kernel::gauss(diff_init, pi.smth, n_dim_);
        auto f0 = -k * stretch * diff / dist;
        // SPH
        auto f1 = -pi.mass * pj.mass *
                  (pi.pres / (pi.dens * pi.dens) + pj.pres / (pj.dens * pj.dens)) *
                  kernel::grad_gauss(diff, pi.smth, n_dim_);
        force = f0 + f1;
        // if(static size_t i; omp_get_thread_num() == 0 and i++ % 1000000 == 0)
        //     cerr << f0.length() << '\t' << f1.length() << '\t' << force.length() << '\n';
    }
    //--------------------------------------------------------------------------------
    else
    {
        throw sh::runtime_error{} << "Unknown type: " << type;
    }
    //--------------------------------------------------------------------------------
#endif

    // spring
    auto k = pi.spring_const * kernel::gauss(diff_init, pi.smth, n_dim_);
    auto f0 = -k * stretch * diff / dist;
    // SPH
#ifdef USE_VOLUME_BASED_SPH
    auto f1 = -pi.pres * (pi.mass / pi.dens) * (pi.mass / pi.dens) *
                  kernel::grad_gauss(diff, pi.smth, n_dim_) -
              pj.pres * (pj.mass / pj.dens) * (pj.mass / pj.dens) *
                  kernel::grad_gauss(diff, pj.smth, n_dim_);
#else
    auto f1 = -pi.mass * pj.mass *
              (pi.pres / (pi.dens * pi.dens) + pj.pres / (pj.dens * pj.dens)) *
              kernel::grad_gauss(diff, pi.smth, n_dim_);
#endif
    auto force = f0 + f1;
    if(static size_t i; omp_get_thread_num() == 0 && i++ % 1000000 == 0)
        cerr << f0.length() << '\t' << f1.length() << '\t' << force.length() << '\n';

    result.force += force;
    result.eng_dot -= force * pi.vel;

    if(force.is_nan() || force.is_inf())
    {
        sh::runtime_error err;
        err << std::scientific << std::showpos;
        err << "Detected NaN or Inf:\n";
        err << "id:\t" << pi.id << '\t' << pj.id << '\n';
        err << "dens:\t" << pi.dens << '\t' << pj.dens << '\n';
        err << "init_dens:\t" << pi.dens << '\t' << pj.dens << '\n';
        err << "pos:\t" << pi.pos << pj.pos << '\n';
        err << "init_pos:\t" << pi.init_pos << pj.init_pos << '\n';
        err << "vel:\t" << pi.vel << pj.vel << '\n';
        err << "force:\t" << result.force << '\n';
        err << "eng_dot:\t" << result.eng_dot << '\n';
        err << "rest_len:\t" << rest_len << '\n';
        throw err;
    }

    ++result.n_ngb;
}

void
repulsion::operator()(const particles::essential::repulsion* const arr_pi, const PS::S32 n_pi,
                      const particles::essential::repulsion* const arr_pj, const PS::S32 n_pj,
                      interaction::result::repulsion* results) const
{
    auto view_pi = std::span{ arr_pi, (std::size_t)n_pi };
    auto view_pj = std::span{ arr_pj, (std::size_t)n_pj };

#ifdef SELF_IMPL_PERIODIC_BOUNDARY
    std::vector view_pj_unique(view_pj.begin(), view_pj.end());
    std::ranges::sort(view_pj_unique, {}, &particles::essential::repulsion::id);
    auto last =
        std::ranges::unique(view_pj_unique, {}, &particles::essential::repulsion::id).begin();
    view_pj_unique.erase(last, view_pj_unique.end());

    std::vector<particles::essential::repulsion> view_pj_all;
    for(const auto& pj : view_pj_unique)
    {
        enum shift
        {
            minus,
            none,
            plus,
        };
        vector_t<std::vector<shift>> image_shifts;
        for(auto d : vector::indices())
        {
            // A real particle-j always exists
            image_shifts[d].push_back(none);

            if(not periodic_domain_[d]) continue;
            auto side = periodic_domain_[d].value();
            if(side.max < pj.pos[d] + pj.r_search) image_shifts[d].push_back(minus);
            if(pj.pos[d] - pj.r_search < side.min) image_shifts[d].push_back(plus);
        }
        for(const auto& shift : sh::views::direct_product(image_shifts))
        {
            assert(shift.size() == vector::size());
            auto& p = view_pj_all.emplace_back(pj);
            for(auto d : vector::indices())
            {
                if(not periodic_domain_[d]) continue;
                auto bound = periodic_domain_[d]->bound();
                if(shift[d] == minus) p.pos[d] -= bound;
                if(shift[d] == plus) p.pos[d] += bound;
            }
        }
    }
#endif

    for(const auto& pi : view_pi)
    {
        auto& result = *results++;
        result.clear();

#ifdef SELF_IMPL_PERIODIC_BOUNDARY
        for(const auto& pj : view_pj_all) calc_force(pi, pj, result);
#else
        for(const auto& pj : view_pj) calc_force(pi, pj, result);
#endif
    }
}

void
repulsion::calc_force(const particles::essential::repulsion& pi,
                      const particles::essential::repulsion& pj,
                      interaction::result::repulsion& result) const
{
    if(pi.unit_id and pj.unit_id and pi.unit_id.value() == pj.unit_id.value()) return;

    auto diff = pi.pos - pj.pos;
    auto dist = diff.length();
    if(dist == 0 or pi.r_search < dist) return;

    vector force;
    //---------------------------------------------------------------------
    // Linear
    //---------------------------------------------------------------------
    force = -repulsion_factor_ * (dist - pi.r_search) * diff / dist;
    //---------------------------------------------------------------------
    // Density-dependent
    //---------------------------------------------------------------------
    // auto f = repulsion_factor_ * ((pi.dens + pj.dens) * 0.5 - dens_std_);
    // force = f * diff / dist;
    // if(static size_t i{}; i++ % 10000000 == 0) cerr << f << '\n';
    //---------------------------------------------------------------------
    result.force += force;

    if(result.force.is_nan() || result.force.is_inf())
    {
        sh::exception err;
        err << std::scientific << std::showpos;
        err << "Detected NaN or Inf:\n";
        err << "id:    " << pi.id << '\t' << pj.id << '\n';
        err << "pos:   " << pi.pos << pj.pos << '\n';
        err << "force: " << result.force << '\n';
        throw err;
    }

    ++result.n_ngb;
}

void
get_springs::operator()(const particles::essential::spring* const arr_pi, const PS::S32 n_pi,
                        const particles::essential::spring* const arr_pj, const PS::S32 n_pj,
                        interaction::result::spring* results)
{
    auto view_pi = std::span{ arr_pi, (std::size_t)n_pi };
    auto view_pj = std::span{ arr_pj, (std::size_t)n_pj };

    std::vector<particles::essential::spring> image_ptcls;
    for(auto i : vector::indices())
    {
        if(not periodic_domain_[i]) continue;
        auto side = periodic_domain_[i].value();
        auto bound = side.bound();

        for(const auto& pj : view_pj)
        {
            if(side.max <= pj.pos[i] + pj.r_search)
            {
                auto& p = image_ptcls.emplace_back(pj);
                p.pos[i] -= bound;
            }
            if(pj.pos[i] - pj.r_search <= side.min)
            {
                auto& p = image_ptcls.emplace_back(pj);
                p.pos[i] += bound;
            }
        }
    }

    auto view_pj_all = sh::views::concat(view_pj, image_ptcls);

    for(const auto& pi : view_pi)
    {
        auto& result = *results++;
        result.clear();

        for(const auto& pj : view_pj_all)
        {
            auto diff = pi.pos - pj.pos;
            auto dist = diff.length();
            if(dist == 0) return;

            auto diff_init = pi.init_pos - pj.init_pos;
            for(auto i : vector::indices())
            {
                if(not periodic_domain_[i]) continue;

                auto bound = periodic_domain_[i]->bound();

                // Ensure that pj does not overlap between real and image
                if(bound - pi.r_search < std::fabs(diff[i])) continue;

                if(diff_init[i] - pi.r_search < -bound)
                {
                    diff_init[i] += bound;
                }
                else if(bound < diff_init[i] + pi.r_search)
                {
                    diff_init[i] -= bound;
                }
            }
            auto rest_len = diff_init.length();
            if(pi.r_search < rest_len) continue;

            auto stretch = dist - rest_len;
            // auto constant =
            //     pi.unit_id and pj.unit_id and pi.unit_id.value() == pj.unit_id.value()
            //         ? pi.spring_const[1]
            //         : pi.spring_const[2];
            auto constant = pi.spring_const;

            assert(springs_);
            std::unique_lock<std::shared_mutex> lg(springs_mutex_);
            springs_->emplace_back(std::make_pair(pi.id, pj.id), rest_len, stretch, constant);

            ++result.n_ngb;
        }
    }
}

// void
// spring::operator()(const particles::essential::spring* const arr_pi, const PS::S32 n_pi,
//                    const particles::essential::spring* const arr_pj, const PS::S32 n_pj,
//                    interaction::result::spring* results) const
// {
//     for(const auto& pi : std::span{ arr_pi, (std::size_t)n_pi })
//     {
//         auto& result = *results++;
//         result.clear();

//         // int n_ngb{};

//         const auto r_search = pi.r_search;

//         // std::set<std::size_t> interacted;

//         // std::cerr << "n_pj = " << n_pj << '\n';

//         for(const auto& pj : std::span{ arr_pj, (std::size_t)n_pj })
//         {
//             const auto [found, spring] = pi.springs.find(pj.id);

//             if(not found) continue;

//             const auto rel_pos = pi.pos - pj.pos;
//             const auto dist = rel_pos.length();

//             if(dist == 0 or r_search < dist)
//             {
//                 std::cout << "skipped (" << pi.id << "," << pj.id << ")\n";
//                 std::cout << "pi.pos  : " << pi.pos << '\n';
//                 std::cout << "pj.pos  : " << pj.pos << '\n';
//                 std::cout << "distance : " << dist << '\n';
//                 std::cout << std::endl;
//                 continue;
//             }

//             // if( interacted.find( pj.id ) != interacted.end() ) continue;
//             // interacted.insert( pj.id );

//             // if( r_search < dist and false )
//             // {
//             //     std::cerr << "Cutoff" << '\n'
//             //               << "id:  " << pi.id << '\t' << pj.id << '\n'
//             //               << "pos: " << pi.pos << pj.pos << '\n'
//             //               << "-" << '\n';
//             //     continue;
//             // }
//             // std::cout << pi.id << ' ' << pj.id << ' ' << r_search << ' ' << dist <<
//             std::endl;

//             //----------------------------------------------------
//             auto stretch = dist - spring.length;
//             auto factor = 0.0;
//             for(auto [pow, att, rep] : spring.constant) factor += (spring.length < dist ? att
//             : rep) * std::pow(stretch, pow);
//             //----------------------------------------------------
//             static auto count{ 0UL };
//             if(count++ % 10000000UL == 0 and in_root_proc)
//             {
//                 for(auto [pow, att, rep] : spring.constant)
//                     cerr << (spring.length < dist ? att : rep) * std::pow(stretch, pow) <<
//                     '\t';
//                 cerr << std::endl;
//             }
//             //----------------------------------------------------

//             auto force = -factor * rel_pos / dist;

//             result.force += force;

//             result.eng_dot -= force * pi.vel;

//             // std::cout << stretch << '\t' << factor << '\t' << result.force << std::endl;

//             if(result.force.is_nan() || result.force.is_inf())
//             {
//                 std::cout << "-\nError in interaction functor: Detected NaN or Inf:\n-\n";
//                 std::cout << "id:      " << pi.id << '\t' << pj.id << '\n';
//                 std::cout << "pos:     " << pi.pos << pj.pos << '\n';
//                 std::cout << "vel:     " << pi.vel << pj.vel << '\n';
//                 std::cout << "mass:    " << pi.mass << '\t' << pj.mass << '\n';
//                 std::cout << "dens:    " << pi.dens << '\t' << pj.dens << '\n';
//                 std::cout << "pres:    " << pi.pres << '\t' << pj.pres << '\n';
//                 std::cout << "force:   " << result.force << '\n';
//                 std::cout << "eng_dot: " << result.eng_dot << '\n';
//                 std::cout << "spr_len: " << spring.length << '\n';
//                 std::cout << "-" << std::endl;
//                 PS::Abort();
//             }

//             // ++n_ngb;
//             ++result.n_ngb;

//             // if( pi.id == 0 ) std::cout << pj.id << '\t' << dist << '\t' << spr_frc <<
//             '\n';
//         }

//         if(result.n_ngb != pi.springs.size())
//         {
//             std::cout << "r_search: " << pi.r_search << '\n';
//             std::cout << "pos: " << pi.pos << '\n';
//             std::cout << "id: " << pi.id << '\n';
//             std::cout << std::flush;
//             throw exception{} << "# of particles interacted is inconsistent with the size
//             of spring info map. size: "
//                                 << pi.springs.size() << ", interacted: " << result.n_ngb;
//         }
//     }
// }

// void
// SPH_like::operator()(const particles::essential::spring* const arr_pi, const PS::S32 n_pi,
//                      const particles::essential::spring* const arr_pj, const PS::S32 n_pj,
//                      interaction::result::spring* results) const
// {
//     for(const auto& pi : std::span{ arr_pi, (std::size_t)n_pi })
//     {
//         auto& result = *results++;
//         result.clear();

//         const auto r_search = pi.getRSearch();

//         for(const auto& pj : std::span{ arr_pj, (std::size_t)n_pj })
//         {
//             const auto [found, spr] = pi.springs.find(pj.id);

//             if(not found)
//             {
//                 std::cerr << "neighbor info not found\n";
//                 continue;
//             }

//             const auto rel_pos = pi.pos - pj.pos;
//             const auto dist = rel_pos.length();

//             if(dist == 0) continue;

//             auto stretch = dist - spr.length;
//             auto factor_spr = 0.0;
//             for(auto [pow, att, rep] : spr.constant) factor_spr += (spr.length < dist ? att :
//             rep) * std::pow(stretch, pow); auto spring_force = factor_spr * rel_pos / dist;

//             double coef = 0;

//             coef += pi.pres / (pi.dens * pi.dens);
//             coef += pj.pres / (pj.dens * pj.dens);

//             result.force -= pi.mass * pj.mass * coef * spring_force;

//             result.eng_dot -= 0.5 * pj.mass * (pj.vel - pi.vel) * spring_force;

//             if(result.force.is_nan() || result.force.is_inf())
//             {
//                 std::cerr << '\n'
//                           << "-" << '\n'
//                           << "Error in interaction functor: Detected NaN: " << '\n'
//                           << "id:      " << pi.id << '\t' << pj.id << '\n'
//                           << "pos:     " << pi.pos << pj.pos << '\n'
//                           << "vel:     " << pi.vel << pj.vel << '\n'
//                           << "mass:    " << pi.mass << '\t' << pj.mass << '\n'
//                           << "dens:    " << pi.dens << '\t' << pj.dens << '\n'
//                           << "pres:    " << pi.pres << '\t' << pj.pres << '\n'
//                           << "force:   " << result.force << '\n'
//                           << "spr_len: " << spr.length << '\n'
//                           << "eng_dot: " << result.eng_dot << '\n'
//                           << "-" << '\n';
//                 PS::Abort();
//             }
//         }
//     }
// }

Lanvegin::Lanvegin(double coef, double interaction_radius)
    : _coef(coef),
      _r_int(interaction_radius),
      _r_int_sq(interaction_radius * interaction_radius)
{
}

void
Lanvegin::operator()(const particles::full_ptcl* arr_pi, const PS::S32 n_pi,
                     const particles::full_ptcl* arr_pj, const PS::S32 n_pj,
                     interaction::result::spring* const forces) const
{
    for(int i = 0; i < n_pi; ++i)
    {
        auto& force = forces[i];
        force.clear();
        const auto& ptcl = arr_pi[i];
        const double inv_mass = 1 / ptcl.mass;

        for(int j = 0; j < n_pj; ++j)
        {
            const auto& cont = arr_pj[j];

            if(ptcl.id == cont.id) continue;

            const auto rel_pos = ptcl.pos - cont.pos;
            const auto dist_sq = rel_pos.square();
            const auto dist = std::sqrt(dist_sq);
            const auto ratio = dist / _r_int;

            const auto acc_coef = (std::pow(ratio, -13) - std::pow(ratio, -7)) * _coef;

            force.force += acc_coef * rel_pos / dist;
        }
    }
}

visualize::visualize(std::shared_ptr<list_type> bond_list) : _bond_list(bond_list)
{
}

void
visualize::operator()(const particles::essential::spring* arr_pi, const PS::S32 n_pi,
                      const particles::essential::spring* arr_pj, const PS::S32 n_pj,
                      interaction::result::spring* const forces)
{
    for(int i = 0; i < n_pi; ++i)
    {
        forces[i].clear();
        const auto& ptcl = arr_pi[i];
        const auto r_search_sq = std::pow(ptcl.getRSearch(), 2);
        for(int j = 0; j < n_pj; ++j)
        {
            const auto& cont = arr_pj[j];
            if(&ptcl == &cont) continue;

            _bond_list->emplace_back(ptcl, cont);
        }
    }
}

density::density(std::size_t n_dim, vector_t<std::optional<sh::min_max<double>>> domain)
    : n_dim_{ n_dim }, periodic_domain_{ domain }
{
}

void
density::operator()(const particles::essential::density* const arr_pi, const PS::S32 n_pi,
                    const particles::essential::density* const arr_pj, const PS::S32 n_pj,
                    interaction::result::density* results) const
{
    auto view_pi = std::span{ arr_pi, (std::size_t)n_pi };
    auto view_pj = std::span{ arr_pj, (std::size_t)n_pj };

#ifdef SELF_IMPL_PERIODIC_BOUNDARY
    // for periodic boundary, the following is needed
    // make unique view of pj
    std::vector view_pj_unique(view_pj.begin(), view_pj.end());
    std::ranges::sort(view_pj_unique, {}, &particles::essential::density::id);
    auto last =
        std::ranges::unique(view_pj_unique, {}, &particles::essential::density::id).begin();
    view_pj_unique.erase(last, view_pj_unique.end());
    // image particles
    std::vector<particles::essential::density> view_pj_all;
    for(const auto& pj : view_pj_unique)
    {
        enum shift
        {
            minus,
            none,
            plus,
        };
        vector_t<std::vector<shift>> image_shifts;
        for(auto d : vector::indices())
        {
            // A real particle-j always exists
            image_shifts[d].push_back(none);

            if(not periodic_domain_[d]) continue;
            auto side = periodic_domain_[d].value();
            if(side.max < pj.pos[d] + pj.r_search) image_shifts[d].push_back(minus);
            if(pj.pos[d] - pj.r_search < side.min) image_shifts[d].push_back(plus);
        }
        for(const auto& shift : sh::views::direct_product(image_shifts))
        {
            assert(shift.size() == vector::size());
            auto& p = view_pj_all.emplace_back(pj);
            for(auto d : vector::indices())
            {
                if(not periodic_domain_[d]) continue;
                auto bound = periodic_domain_[d]->bound();
                if(shift[d] == minus) p.pos[d] -= bound;
                if(shift[d] == plus) p.pos[d] += bound;
            }
        }
    }
#else
    // assert(!periodic_domain_[0] && !periodic_domain_[1] && !periodic_domain_[2]);
#endif

    for(const auto& pi : view_pi)
    {
        auto& result = *results++;
        result.clear();

#ifdef SELF_IMPL_PERIODIC_BOUNDARY
        for(const auto& pj : view_pj_all)
#else
        for(const auto& pj : view_pj)
#endif
        {
            auto diff = pi.pos - pj.pos;

            if(pi.r_search * pi.r_search < diff.square()) continue;

#ifdef USE_VOLUME_BASED_SPH
            result.dens += kernel::gauss(diff, pi.smth, n_dim_);
#else
            result.dens += pj.mass * kernel::gauss(diff, pi.smth, n_dim_);
#endif

            ++result.n_ngb;
        }

#ifdef USE_VOLUME_BASED_SPH
        result.dens *= pi.mass;
#endif

#if 0
        // debug
        if(static auto r_search = pi.r_search; r_search != pi.r_search) throw sh::exception{};
        if(result.n_ngb != 5)
        {
            sh::exception err;
            err << std::scientific << std::showpos;
            err << '\n';
            err << "result.n_ngb:\t" << result.n_ngb << "\n";
            err << "r_search:\t" << pi.r_search << '\n';
            err << "pi:\n";
            err << pi.id << '\t' << pi.pos << '\n';
            err << "view_pj_unique:\n";
            for(const auto& pj : view_pj_unique)
            {
                // if(pi.r_search < (pi.pos - pj.pos).length()) continue;
                err << pj.id << '\t';
                err << pj.pos;
                err << (pi.r_search < (pi.pos - pj.pos).length() ? "" : "ok");
                err << '\n';
            }
            err << "view_pj_all:\n";
            for(const auto& pj : view_pj_all)
            {
                // if(pi.r_search < (pi.pos - pj.pos).length()) continue;
                err << pj.id << '\t';
                err << pj.pos;
                err << (pi.r_search < (pi.pos - pj.pos).length() ? "" : "ok");
                err << '\n';
            }
            throw err;
        }
#endif
    }
}

void
clear::operator()(const particles::essential::spring* arr_pi, const PS::S32 n_pi,
                  const particles::essential::spring* arr_pj, const PS::S32 n_pj,
                  interaction::result::spring* const forces) const
{
    for(int i = 0; i < n_pi; ++i) forces[i].clear();
}

void
do_nothing::operator()(const particles::full_ptcl* arr_pi, const PS::S32 n_pi,
                       const particles::full_ptcl* arr_pj, const PS::S32 n_pj,
                       interaction::result::spring* const forces) const
{
    // nothing is here
}

}  // namespace functor
}  // namespace interaction
}  // namespace cfe