#pragma once

#include <vector>
#include <cassert>
#include <unordered_map>
#include <set>
#include <unordered_set>

#include <mpi.h>

#include <sh/ranges.h>

#include <particles.h>

namespace cfe
{

class spring_set
{
   private:
    static constexpr auto M_n_intxn{ 2ul };

    struct spring
    {
        std::vector<std::size_t> ptcl_ids;
        double constant;
        double rest_length;
    };

    std::vector<spring> M_spring_set;

    auto
    calculate_force(std::array<particles::full_ptcl, M_n_intxn> &ptcl_list, const spring &spring)
    {
        assert(M_n_intxn == 2ul);

        auto &[p1, p2] = ptcl_list;

        const auto disp = p2.pos - p1.pos;
        const auto dist = disp.length();
        const auto stretch = dist - spring.rest_length;

        return spring.constant * (stretch / dist) * disp;
    }

   public:
    void
    calculate_force_all(std::vector<particles::full_ptcl> &ptcls_loc)
    {
        // ID-to-local-particle map
        std::unordered_map<std::size_t, std::reference_wrapper<particles::full_ptcl>> id_to_ptcl_loc;
        for(auto &ptcl : ptcls_loc) id_to_ptcl_loc[ptcl.id] = std::ref(ptcl);

        // [proc][id]
        std::vector<std::vector<std::size_t>> needed_ids(n_procs);
        std::vector<std::reference_wrapper<spring>> uncalculated_spring;
        // Collect particle IDs needed by each process
        for(auto &spring : M_spring_set)
        {
            // check if the particles to interact are in the same process
            auto lacked = false;
            std::set<std::reference_wrapper<particles::full_ptcl>> spring_ptcls;
            for(auto id : spring.ptcl_ids)
            {
                auto it = id_to_ptcl_loc.find(id);

                // found in the same process
                if(it != id_to_ptcl_loc.end())
                {
                    if(not lacked) spring_ptcls.insert(it->second);
                }
                // not found
                else
                {
                    needed_ids[rank].push_back(it->first);
                    lacked = true;
                }
            }
            if(not lacked)
            {
                calculate_force(spring_ptcls);
            }
            else
            {
                uncalculated_spring.push_back(std::ref(spring));
            }
        }

        // broadcast needed IDs in each process
        for(auto [src_rank, src_ids] : needed_ids | sh::views::enumerate)
        {
            // broadcast size
            std::size_t src_ids_size;
            if(rank == src_rank) src_ids_size = src_ids.size();
            MPI_Bcast(&src_ids_size, 1, MPI_UNSIGNED_LONG, src_rank, MPI_COMM_WORLD);
            // broadcast data
            if(rank != src_rank) src_ids.resize(src_ids_size);
            MPI_Bcast(src_ids.data(), src_ids_size, MPI_UNSIGNED_LONG, src_rank, MPI_COMM_WORLD);
        }

        // check if this process has needed ids
        for(auto)
        {
        }

        std::unordered_map<std::size_t, particles::full_ptcl> ptcls_recv;

        // Exchange required particle data using MPI
        std::vector<std::size_t> needed_ids_vector(needed_ids_not_here.begin(), needed_ids_not_here.end());
        std::vector<std::size_t> send_counts(n_procs, 0), recv_counts(n_procs);
        std::vector<std::size_t> send_displs(n_procs, 0), recv_displs(n_procs, 0);

        // Simulate an exchange plan (details depend on domain decomposition)
        // You should define send_counts, recv_counts, send_displs, recv_displs properly.
        // For simplicity, assume send and receive all needed IDs for now.

        // send and receive particle data
        for(auto rank_tar : std::views::iota(0ul, n_procs))
        {
            if(rank_tar == rank) continue;

            // send needed IDs
        }

        // Compute forces for interactions
        for(const auto &interaction : M_spring_set)
        {
            const Particle &p1 = id_to_ptcl_loc[interaction.id1];
            const Particle &p2 = id_to_ptcl_loc[interaction.id2];

            const auto force = calculate_force(p1, p2, interaction);
            id_to_ptcl_loc[interaction.id1].force += force;
            id_to_ptcl_loc[interaction.id2].force -= force;
        }

        // Reduce forces to the owning processes
        // Use MPI_Reduce or custom logic to combine forces
    }
};
}  // namespace cfe