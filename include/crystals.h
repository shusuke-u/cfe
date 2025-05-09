#include <random>

#include <sh/ranges.h>

#include <util.h>

namespace cfe
{
namespace crystals
{

class unit_base
{
   protected:
    const std::size_t n_dim_;

   public:
    vector side;
    double min_sep;
    bool use_id;

   public:
    unit_base(std::size_t n_dim) : n_dim_{ n_dim }
    {
    }
    unit_base(std::size_t n_dim, vector side, double min_sep, bool use_id)
        : n_dim_{ n_dim }, side{ side }, min_sep{ min_sep }, use_id{ use_id }
    {
    }

    void
    rescale(double factor)
    {
        side *= factor;
        min_sep *= factor;
    }
    virtual std::size_t n_positions() const = 0;
    virtual std::vector<vector> positions(double factor) const = 0;
};
class Cartesian : public unit_base
{
   public:
    Cartesian(std::size_t n_dim) : unit_base{ n_dim, vector(1), 1, false }
    {
    }
    std::vector<vector>
    positions(double) const override
    {
        return { vector(0) };
    }
    std::size_t
    n_positions() const override
    {
        return 1;
    }
};
class densest : public unit_base
{
   private:
    std::vector<vector> pos_buf;

   public:
    densest(std::size_t n_dim) : unit_base{ n_dim }
    {
        use_id = false;

        vector n1{ 0.5, 0.8660254038, 0.0 };
        vector n2{ 0.0, 0.5773502692, 0.8164965809 };
        n1.truncate(n_dim);
        n2.truncate(n_dim);

        if(1UL == n_dim)
        {
            pos_buf.emplace_back();
            side = 2 * n1;
        }
        if(2UL == n_dim)
        {
            pos_buf.push_back(+0.5 * n1);
            pos_buf.push_back(-0.5 * n1);
            side = 2 * n1;
        }
        if(3UL == n_dim)
        {
            pos_buf.push_back(+0.5 * n1 + 0.5 * n2);
            pos_buf.push_back(-0.5 * n1 + 0.5 * n2);
            pos_buf.push_back(+0.5 * n1 - 0.5 * n2);
            pos_buf.push_back(-0.5 * n1 - 0.5 * n2);
            side[0] = n1[0] + n2[0];
            side[1] = n1[1];
            side[2] = n1[2] + n2[2];
            side *= 2;
        }
        min_sep = 1;
    }
    std::vector<vector>
    positions(double factor) const override
    {
        std::vector<vector> ret;
        for(const auto& pos : pos_buf) ret.push_back(pos * factor);
        return ret;
    }
    std::size_t
    n_positions() const override
    {
        switch(n_dim_)
        {
            case 1: return 1;
            case 2: return 2;
            case 3: return 4;
            default: return 0;
        };
    }
};
class tetrahedron : public unit_base
{
   public:
    tetrahedron(std::size_t n_dim) : unit_base{ n_dim }
    {
        use_id = true;
        side = vector(1, sqrt(3) / 2, sqrt(6) / 3);
        min_sep = 1;
    }
    std::vector<vector>
    positions(double factor) const override
    {
        using std::sqrt;

        std::vector<vector> pos_buf;
        pos_buf.emplace_back();
        if(1 <= n_dim_) pos_buf.emplace_back(1.0, 0.0, 0.0);
        if(2 <= n_dim_) pos_buf.emplace_back(0.5, sqrt(3.0) / 2.0, 0.0);
        if(3 <= n_dim_) pos_buf.emplace_back(0.5, sqrt(3.0) / 3.0, sqrt(6.0) / 3.0);
        vector center;
        for(const auto& pos : pos_buf) center += pos;
        center /= pos_buf.size();
        for(auto& pos : pos_buf) pos -= center;

        using std::numbers::pi;
        static std::random_device rd;
        static std::mt19937 gen{ rd() };
        static std::uniform_real_distribution<double> uni_dis{ 0, 2 * pi };

        double psi{}, theta{}, phi{};
        if(2 <= n_dim_) psi = uni_dis(gen);
        if(3 <= n_dim_) theta = uni_dis(gen), phi = uni_dis(gen);
        for(auto& pos : pos_buf) pos.rotate(psi, theta, phi);

        static std::normal_distribution<double> nom_dis{ 0, 0.5 * side.length() };
        auto shift = vector::unit();
        shift.normalize(nom_dis(gen));
        if(2 <= n_dim_) psi = uni_dis(gen);
        if(3 <= n_dim_) theta = uni_dis(gen), phi = uni_dis(gen);
        shift.rotate(psi, theta, phi);
        for(auto& pos : pos_buf) pos += shift;

        for(auto& pos : pos_buf) pos *= factor;

        return pos_buf;
    }
    std::size_t
    n_positions() const override
    {
        return 1 + n_dim_;
    }
};

class arrange_fn
{
   public:
    struct pos_and_unit_id
    {
        std::optional<std::size_t> unit_id;
        vector pos;
    };
    using result_type = std::tuple<std::vector<pos_and_unit_id>, double>;

    result_type
    operator()(std::derived_from<unit_base> auto unit, orthotope& domain, std::size_t& n_ptcl,
               std::size_t n_dim) const
    {
        assert(n_ptcl != 0);
        assert(1 <= n_dim and n_dim <= 3);

        const auto indices = std::views::iota(0UL, n_dim);

        /// \# of particles in an unit
        const auto n_in_unit = unit.n_positions();
        /// side length of the domain
        const auto side = domain.diag();

        auto n_unit = n_ptcl / n_in_unit;

        if(n_unit == 0)
            throw sh::INVALID_VALUE(n_ptcl)
                << "The number of particles is too small. n_in_unit: " << n_in_unit;

        const auto factor =
            std::pow(side.prod_n(n_dim) / unit.side.prod_n(n_dim) / n_unit, 1. / n_dim);
        unit.rescale(factor);

        using namespace std::views;

        // \# of particles along the axes
        vector_t<std::size_t> n_unit_along;
        for(auto i : indices) n_unit_along[i] = std::trunc(side[i] / unit.side[i]);
        // vector_t<std::size_t> n_unit_along{ sh::call( std::trunc, sh::call( std::divides<>{},
        // side, unit_side ) ) };
        n_unit_along.truncate(n_dim, 1UL);
        // vector_t<std::size_t> n_unit_along =
        //     sh::call( &std::trunc, sh::call( std::divides<>{}, side, unit_side )
        //     ).make_truncated();

        // Update # of total units and particles.
        n_unit = n_unit_along.prod_n(n_dim);
        n_ptcl = n_unit * n_in_unit;

        // Shrunk the domain so that the length of the domain side is to be an integer multiple
        // of spacing. Use a custom function since `std::fmod` has bad behaviour. e.g.
        // `std::fmod( 1, 0.1 )` returns `0.1` rather than `0`. auto mod = []( double x, double
        // y ) { return x - std::floor( x / y ) * y; }; domain.extend( -(vector)sh::call( mod,
        // side, unit_side ) );

        // domain.lower() -= 0.5 * unit_side.make_truncated( n_dim );
        for(auto i : indices)
            domain.upper()[i] = domain.lower()[i] + n_unit_along[i] * unit.side[i];

        // buffer for the positions of the particles
        std::vector<pos_and_unit_id> ptcl_buf;
        ptcl_buf.reserve(n_ptcl);

        auto pos_init = 0.5 * unit.side.make_truncated(n_dim);

        auto unit_pos_set = sh::views::direct_product(
            sh::call(sh::views::linspace, pos_init, n_unit_along, unit.side));

        for(auto [unit_id, unit_pos] : unit_pos_set | sh::views::enumerate)
        {
            for(auto rel_pos : unit.positions(factor))
            {
                if(unit.use_id)
                    ptcl_buf.emplace_back(unit_id, vector(unit_pos) + rel_pos);
                else
                    ptcl_buf.emplace_back(std::nullopt, vector(unit_pos) + rel_pos);
            }
        }
        // This is inconsistent with the actual value of the mean spacing, while the value of
        // this is acceptable. const auto spacing = std::pow( unit_side.prod_n( n_dim ) /
        // unit.pos_buf.size(), 1.0 / n_dim );

        return { ptcl_buf, unit.min_sep };
    }
};
constexpr inline arrange_fn arrange{};

}  // namespace crystals
}  // namespace cfe