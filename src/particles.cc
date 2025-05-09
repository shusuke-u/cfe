#include <util.h>
#include <particles.h>

namespace cfe
{

namespace particles
{

full_ptcl::full_ptcl()
{
    is_fixed.fill(false);
    belongs_to.fill(false);
}

double
full_ptcl::sound_speed() const
{
    // auto type = spring_const[0];
    // if(type == 0.)
    // {
    //     return std::sqrt(spring_const[1] / mass) * smth;
    // }
    // else if(type == 3.1)
    // {
    //     // auto dens_diff = dens - init_dens;
    //     auto k = spring_const[1] + spring_const[2];
    //     return std::sqrt(k / mass) * smth;
    // }
    // else if(type == 3.2)
    // {
    //     auto k = spring_const[1] + spring_const[2] * dens / dens_std;
    //     return std::sqrt(k / mass) * smth;
    // }
    // else if(type == 4.1)
    // {
    //     auto dens_diff = dens - init_dens;
    //     auto k = spring_const[1] + spring_const[2] * std::pow(dens_diff / dens_diff_std, 2);
    //     return std::sqrt(k / mass) * smth;
    // }
    // else if(type == 4.2)
    // {
    //     auto k0 = spring_const[1];
    //     auto k1 = spring_const[2] * dens / dens_std;
    //     return std::sqrt((k0 + k1) / mass) * smth;
    // }
    // else if(type==5.)
    // {
    //     auto Cs0 = std::sqrt(spring_const[1] / mass) * smth;
    //     auto Cs1 = spring_const[2];
    //     return std::max(Cs0, Cs1);
    // }
    // else
    // {
    //     throw sh::runtime_error{} << "Unknown type: " << type;
    // }

    auto cs = std::sqrt(spring_const / mass) * smth;
    cs = std::max(cs, std::sqrt(bulk_sound_sq));
    return cs;
}

double
full_ptcl::time_step() const
{
    return smth / sound_speed();
}

PS::F64vec
full_ptcl::getPos() const
{
    return util::convert(pos);
}

void
full_ptcl::setPos(const PS::F64vec& src)
{
    pos = util::convert(src);
}

std::ostream&
operator<<(std::ostream& os, full_ptcl::header_t)
{
    using sh::del;
    os << "id" << del;
    os << "owner_id" << del;
    os << "unit_id" << del;
    os << "dens" << del;
    // os << "pres" << del;
    os << "eng" << del;
    os << "eng_dot" << del;
    os << "mass" << del;
    os << "n_ngb_spr" << del;
    // os << "n_ngb_spr_init" << del;
    os << "n_ngb_dens" << del;
    // os << "n_ngb_rep" << del;
    os << vector::suffix("pos");
    os << vector::suffix("vel");
    os << vector::suffix("pos_init");
    // os << vector::suffix("vel_init");
    os << vector::suffix("force");
    return os;
}

const std::vector<std::string> full_ptcl::header_v = []
{
    std::stringstream buf;
    buf << header;
    std::vector<std::string> header;
    header.insert(header.begin(), std::istream_iterator<std::string>(buf),
                  std::istream_iterator<std::string>());
    return header;
}();

std::ostream&
operator<<(std::ostream& os, const full_ptcl& src)
{
    using sh::del;

    os << src.id << del;
    os << src.owner_id << del;

    os << (src.unit_id ? std::to_string(*src.unit_id) : "None") << del;

    os << src.dens << del;
    // os << src.pres << del;
    os << src.eng << del;
    os << src.eng_dot << del;
    os << src.mass << del;
    os << src.n_ngb_spr << del;
    // os << src.springs.size() << del;
    os << src.n_ngb_dens << del;
    // os << src.n_ngb_rep << del;
    os << src.pos;
    os << src.vel;
    os << src.init_pos;
    // os << src.vel_init;
    os << src.force;

    return os;
}

std::istream&
operator>>(std::istream& is, full_ptcl& tar)
{
    is >> tar.id;
    is >> tar.owner_id;

    is >> std::ws;
    if(is.peek() == 'N')
    {
        is.ignore(4);
        tar.unit_id.reset();
    }
    else
    {
        tar.unit_id.emplace();
        is >> tar.unit_id.value();
    }

    is >> tar.dens;
    // is >> tar.pres;
    is >> tar.eng;
    is >> tar.eng_dot;
    is >> tar.mass;
    is >> tar.n_ngb_spr;
    // is >> tar.springs;
    is >> tar.n_ngb_dens;
    // is >> tar.n_ngb_rep;
    is >> tar.pos;
    is >> tar.vel;
    is >> tar.init_pos;
    // is >> tar.vel_init;
    is >> tar.force;

    return is;
}

// void
// full_ptcl::copyFromForce(const interaction::result::make_spring& spr)
// {
//     springs = spr.springs;
// }

void
full_ptcl::copyFromForce(const interaction::result::spring& src)
{
    force += src.force;
    eng_dot = src.eng_dot;
    n_ngb_spr = src.n_ngb;
}

void
full_ptcl::copyFromForce(const interaction::result::repulsion& src)
{
    force += src.force;
    n_ngb_rep = src.n_ngb;
}

void
full_ptcl::copyFromForce(const interaction::result::density& src)
{
    dens = src.dens;
    n_ngb_dens = src.n_ngb;
}

// void
// essential::make_spring::copyFromFP(const full_ptcl& src)
// {
//     id = src.id;
//     owner_id = src.owner_id;

//     pos = src.pos;

//     r_search = src.r_search_for_spring;
//     smth = src.smth;

//     constant = src.spring_const;
// }

void
essential::spring::copyFromFP(const full_ptcl& src)
{
    // springs = src.springs;

    id = src.id;
    unit_id = src.unit_id;

    mass = src.mass;
    pres = src.pres;
    std_dens = src.std_dens;
    dens = src.dens;
    smth = src.smth;

    spring_const = src.spring_const;

    vel = src.vel;
    pos = src.pos;
    init_pos = src.init_pos;

    r_search = src.r_search_for_spring;
}

void
essential::repulsion::copyFromFP(const full_ptcl& src)
{
    id = src.id;
    unit_id = src.unit_id;
    r_search = src.r_search_for_spring;
    dens = src.dens;
    pos = src.pos;
}
void
essential::density::copyFromFP(const full_ptcl& src)
{
    id = src.id;
    pos = src.pos;
    dens = src.dens;
    mass = src.mass;
    smth = src.smth;
    r_search = src.r_search_for_dens;
}

// void
// essential::make_spring::setPos(const PS::F64vec& src)
// {
//     pos = util::convert(src);
// }

void
essential::spring::setPos(const PS::F64vec& src)
{
    pos = util::convert(src);
}

void
essential::repulsion::setPos(const PS::F64vec& src)
{
    pos = util::convert(src);
}

void
essential::density::setPos(const PS::F64vec& src)
{
    pos = util::convert(src);
}

// PS::F64vec
// essential::make_spring::getPos() const
// {
//     return util::convert(pos);
// }

PS::F64vec
essential::spring::getPos() const
{
    return util::convert(pos);
}

PS::F64vec
essential::repulsion::getPos() const
{
    return util::convert(pos);
}

PS::F64vec
essential::density::getPos() const
{
    return util::convert(pos);
}

// PS::F64
// essential::make_spring::getRSearch() const
// {
//     return r_search;
// }

PS::F64
essential::spring::getRSearch() const
{
    return r_search;
}

PS::F64
essential::repulsion::getRSearch() const
{
    return r_search;
}

PS::F64
essential::density::getRSearch() const
{
    return r_search;
}

}  // namespace particles
}  // namespace cfe