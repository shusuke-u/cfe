#pragma once

#include <vector>
#include <ranges>
#include <optional>

#include <mpi.h>

#include <particle_simulator.hpp>

namespace cfe
{
namespace par
{

template <class Tp>
std::optional<std::vector<Tp>>
gather(const Tp& send_buf, int dst)
{
    std::vector<Tp> recv_buf(n_procs);
    PS::Comm::gather(const_cast<Tp*>(&send_buf), 1, recv_buf.data(), dst);
    if(rank != dst) return {};
    return recv_buf;
}

template <class Tp>
std::optional<std::vector<Tp>>
gather(const std::optional<Tp>& send_buf, int dst)
{
    auto send = send_buf.has_value();
    auto send_v = all_gather(send);  // Collect data presence across processes

    std::vector<int> recv_counts(n_procs, 0);
    std::vector<int> disp(n_procs, 0);
    int total_recv = 0;

    for(int i = 0; i < n_procs; ++i)
    {
        recv_counts[i] = send_v[i] ? 1 : 0;
        total_recv += recv_counts[i];
    }

    for(int i = 1; i < n_procs; ++i) disp[i] = disp[i - 1] + recv_counts[i - 1];

    Tp* send_ptr = send ? &(*send_buf) : nullptr;

    std::vector<Tp> recv_buf;
    if(rank == dst) recv_buf.resize(total_recv);

    PS::Comm::gatherV(send_ptr, send ? 1 : 0, rank == dst ? recv_buf.data() : nullptr,
                      recv_counts.data(), disp.data(), dst);

    return (rank == dst) ? std::optional{ recv_buf } : std::nullopt;
}

template <std::ranges::sized_range R>
std::optional<std::vector<std::ranges::range_value_t<R>>>
gather(R&& send_buf, int dst)
{
    int send_count = std::ranges::size(send_buf);

    // Gather the number of elements each process sends
    auto recv_counts_opt = gather(send_count, dst);

    // All processes call `gatherV`, but only root (`dst`) stores the result.
    std::vector<int> recv_counts(n_procs, 0);
    std::vector<int> disp(n_procs, 0);
    int total_recv_size = 0;

    if(rank == dst)
    {
        recv_counts = *recv_counts_opt;

        // Compute displacements for variable-sized receive buffer
        for(int i = 1; i < n_procs; ++i) disp[i] = disp[i - 1] + recv_counts[i - 1];

        total_recv_size = disp[n_procs - 1] + recv_counts[n_procs - 1];
    }

    // Convert non-contiguous ranges to contiguous buffer
    std::vector<std::ranges::range_value_t<R>> send_copy(send_buf.begin(), send_buf.end());

    // Allocate receive buffer only for root
    std::vector<std::ranges::range_value_t<R>> recv_buf;
    if(rank == dst) recv_buf.resize(total_recv_size);

    // All processes must call gatherV
    PS::Comm::gatherV(send_copy.data(), send_count,
                      rank == dst ? recv_buf.data() : nullptr,  // Non-root sends nullptr
                      recv_counts.data(), disp.data(), dst);

    // Only root returns the gathered result
    if(rank == dst) return recv_buf;
    return std::nullopt;
}
template <class Tp>
std::vector<Tp>
all_gather(const Tp& send_buf)
{
    std::vector<Tp> recv_data(n_procs);
    PS::Comm::allGather(&send_buf, 1, recv_data.data());
    return recv_data;
}

template <class Tp>
std::vector<Tp>
all_gather(const std::optional<Tp>& send_buf)
{
    auto send = send_buf.has_value();
    auto send_v = all_gather(send);
    auto total_recv = std::ranges::count(send_v, true);

    std::vector<int> disp(n_procs, 0);
    int count = 0;
    for(int i = 1; i < n_procs; ++i) disp[i] = disp[i - 1] + (int)send_v[i - 1];

    Tp* send_ptr = nullptr;
    if(send) send_ptr = &(*send_buf);

    std::vector<Tp> recv_buf(total_recv);
    PS::Comm::allGatherV(send_ptr, (int)send, recv_buf.data(), &total_recv, disp.data());

    return recv_buf;
}

template <std::ranges::sized_range R>
// using R = std::vector<double>;
std::vector<std::ranges::range_value_t<R>>
all_gather(const R& send_buf)
{
    int send_count = std::ranges::size(send_buf);

    auto recv_counts = all_gather(send_count);

    // Compute displacements
    std::vector<int> disp(n_procs, 0);
    for(int i = 1; i < n_procs; ++i) disp[i] = disp[i - 1] + recv_counts[i - 1];

    // Allocate receive buffer
    int total_recv_size = disp[n_procs - 1] + recv_counts[n_procs - 1];
    std::vector<std::ranges::range_value_t<R>> recv_buf(total_recv_size);

    // Convert non-contiguous ranges to contiguous buffer
    std::vector<std::ranges::range_value_t<R>> send_copy(send_buf.begin(), send_buf.end());

    PS::Comm::allGatherV(const_cast<std::ranges::range_value_t<R>*>(send_copy.data()),
                         send_count, recv_buf.data(), recv_counts.data(), disp.data());

    return recv_buf;
}

template <class Tp>
void
broadcast(Tp& data, int root)
{
    PS::Comm::broadcast(&data, 1, root);
}

template <std::ranges::sized_range Range>
void
broadcast(Range& data, int root)
{
    auto data_size = std::ranges::size(data);

    broadcast(data_size, root);

    if(PS::Comm::getRank() != root)
    {
        data.resize(data_size);
    }

    PS::Comm::broadcast(data.data(), data_size, root);
}

}  // namespace par
}  // namespace cfe