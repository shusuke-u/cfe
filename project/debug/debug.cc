#include <span>
#include <ranges>

#include <particle_simulator.hpp>

class Force
{
   public:
    void
    clear()
    {
    }
};

class FP
{
   public:
    PS::S32 id;
    PS::F64vec pos;
    PS::F64 r_search;

    auto
    getPos() const
    {
        return pos;
    }

    void
    setPos(PS::F64vec src)
    {
        pos = src;
    }

    void
    copyFromForce(const Force &)
    {
    }
};

class EP
{
   public:
    PS::S32 id;
    PS::F64vec pos;
    PS::F64 r_search;

    auto
    getPos() const
    {
        return pos;
    }
    void
    setPos(const PS::F64vec &src)
    {
        pos = src;
    }
    auto
    getRSearch() const
    {
        return r_search;
    }
    void
    copyFromFP(const FP &fp)
    {
        id = fp.id;
        pos = fp.pos;
        r_search = fp.r_search;
    }
};

class CalcForce
{
   public:
    void
    operator()(const EP *epis, const PS::S32 ni, const EP *epjs, const PS::S32 nj, Force *force)
    {
        std::vector<EP> images;
        for(const auto &epj : std::span{ epjs, (size_t)nj })
        {
            if(10 < epj.pos.x + epj.r_search)
            {
                auto pos = epj.pos;
                pos.x -= 10;
                images.push_back(epj);
                images.back().pos = pos;
            }
            if(epj.pos.x - epj.r_search < 0)
            {
                auto pos = epj.pos;
                pos.x += 10;
                images.push_back(epj);
                images.back().pos = pos;
            }
        }

        for(const auto &epi : std::span{ epis, (size_t)ni })
        {
            auto r_search_sq = epi.r_search * epi.r_search;

            for(const auto &epj : std::span{ epjs, (size_t)nj })
            {
                if(epi.id == epj.id) continue;
                auto disp = epj.pos - epi.pos;
                if(r_search_sq < disp * disp) continue;
                std::cerr << "id: " << epi.id << '\t' << epj.id << '\n';
                std::cerr << "x : " << epi.pos.x << '\t' << epj.pos.x << '\n';
            }
            for(const auto &epj : images)
            {
                if(epi.id == epj.id) continue;
                auto disp = epj.pos - epi.pos;
                if(r_search_sq < disp * disp) continue;
                std::cerr << "id: " << epi.id << '\t' << epj.id << '\n';
                std::cerr << "x : " << epi.pos.x << '\t' << epj.pos.x << '\n';
            }
        }
    }
};

int
main(int argc, char **argv)
{
    PS::Initialize(argc, argv);

    auto rank = PS::Comm::getRank();
    auto procs = PS::Comm::getNumberOfProc();

    PS::ParticleSystem<FP> sys;
    sys.initialize();

    PS::DomainInfo dinfo;
    dinfo.initialize();
    //    dinfo.setBoundaryCondition( PS::BOUNDARY_CONDITION_PERIODIC_XYZ );
    dinfo.setPosRootDomain({ 0 }, { 10 });

    auto n_ptcl{ 3 };
    sys.setNumberOfParticleLocal(n_ptcl);

    for(auto i : std::views::iota(0, n_ptcl))
    {
        sys[i].id = i;
        sys[i].r_search = 2;
    }

    //    PS::TreeForForceShort<Force, EP, EP>::Gather tree;
    PS::TreeForForce<PS::SEARCH_MODE_GATHER, Force, EP, EP, PS::MomentShort, PS::MomentShort, PS::SuperParticleBase,
                     PS::CALC_DISTANCE_TYPE_NEAREST_XYZ>
        tree;

    tree.initialize(n_ptcl);

    auto step{ 0 };

    if(rank == 0) std::cerr << step++ << '\n';
    sys[0].pos.x = 8.5;
    sys[1].pos.x = 9.5;
    sys[2].pos.x = 5.0;
    sys.adjustPositionIntoRootDomain(dinfo);
    dinfo.decomposeDomainAll(sys);
    sys.exchangeParticle(dinfo);
    tree.calcForceAllAndWriteBack(CalcForce{}, sys, dinfo, true, PS::MAKE_LIST_FOR_REUSE);
    // tree.calcForceAllAndWriteBack( CalcForce{}, sys, dinfo, true, PS::MAKE_LIST );

    if(rank == 0) std::cerr << step++ << '\n';
    sys[0].pos.x += 1;
    sys[1].pos.x += 1;
    sys.adjustPositionIntoRootDomain(dinfo);
    dinfo.decomposeDomainAll(sys);
    // sys.exchangeParticle( dinfo );
    tree.calcForceAllAndWriteBack(CalcForce{}, sys, dinfo, true, PS::REUSE_LIST);
    // tree.calcForceAllAndWriteBack( CalcForce{}, sys, dinfo, true, PS::MAKE_LIST );

    if(rank == 0) std::cerr << step++ << '\n';
    sys[0].pos.x -= 2;
    sys[1].pos.x -= 2;
    sys.adjustPositionIntoRootDomain(dinfo);
    dinfo.decomposeDomainAll(sys);
    // sys.exchangeParticle( dinfo );
    tree.calcForceAllAndWriteBack(CalcForce{}, sys, dinfo, true, PS::REUSE_LIST);
    // tree.calcForceAllAndWriteBack( CalcForce{}, sys, dinfo, true, PS::MAKE_LIST );

    PS::Finalize();
}