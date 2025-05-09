#pragma once

#include <complex>
#include <execution>

// #include <eigen-3.4.0/unsupported/Eigen/Polynomials>

#include <sh/callable.h>
#include <sh/exceptions.h>

#include <defs.h>
#include <parallel.h>

namespace cfe
{
namespace util
{

struct square_fn
{
    template <class Tp>
    __attribute__((always_inline)) inline constexpr decltype(auto)
    operator()(Tp a) const noexcept
    {
        return a * a;
    }
};
constexpr inline square_fn sq{};
constexpr inline square_fn square{};

struct cube_fn
{
    template <class Tp>
    __attribute__((always_inline)) inline constexpr decltype(auto)
    operator()(Tp a) const noexcept
    {
        return a * a * a;
    }
};
constexpr inline cube_fn cb{};
constexpr inline cube_fn cube{};

struct convert_fn
{
    template <class Tp>
    constexpr PS::Vector3<Tp>
    operator()(vector_t<Tp>&& src) const
    {
        return { std::move(src[0]), std::move(src[1]), std::move(src[2]) };
    }

    template <class Tp>
    constexpr PS::Vector3<Tp>
    operator()(const vector_t<Tp>& src) const
    {
        return { src[0], src[1], src[2] };
    }

    template <class Tp>
    constexpr vector_t<Tp>
    operator()(PS::Vector3<Tp>&& src) const
    {
        return { std::move(src.x), std::move(src.y), std::move(src.z) };
    }

    template <class Tp>
    constexpr vector_t<Tp>
    operator()(const PS::Vector3<Tp>& src) const
    {
        return { src.x, src.y, src.z };
    }

    template <class Tp>
    PS::Orthotope3<Tp>
    operator()(const orthotope_t<Tp>& src) const
    {
        return { { src.lower()[0], src.lower()[1], src.lower()[2] },
                 { src.upper()[0], src.upper()[1], src.upper()[2] } };
    }

    template <class Tp>
    orthotope_t<Tp>
    operator()(const PS::Orthotope3<Tp>& src) const
    {
        return { vector{ src.low_[0], src.low_[1], src.low_[2] },
                 vector{ src.high_[0], src.high_[1], src.high_[2] } };
    }

    constexpr PS::BOUNDARY_CONDITION
    operator()(const vector_t<bool>& bc) const
    {
        if(bc[0] == true && bc[1] == true && bc[2] == true)
            return PS::BOUNDARY_CONDITION_PERIODIC_XYZ;
        if(bc[0] == true && bc[1] == true && bc[2] == false)
            return PS::BOUNDARY_CONDITION_PERIODIC_XY;
        if(bc[0] == false && bc[1] == true && bc[2] == true)
            return PS::BOUNDARY_CONDITION_PERIODIC_YZ;
        if(bc[0] == true && bc[1] == false && bc[2] == true)
            return PS::BOUNDARY_CONDITION_PERIODIC_XZ;
        if(bc[0] == true && bc[1] == false && bc[2] == false)
            return PS::BOUNDARY_CONDITION_PERIODIC_X;
        if(bc[0] == false && bc[1] == true && bc[2] == false)
            return PS::BOUNDARY_CONDITION_PERIODIC_Y;
        if(bc[0] == false && bc[1] == false && bc[2] == true)
            return PS::BOUNDARY_CONDITION_PERIODIC_Z;

        return PS::BOUNDARY_CONDITION_OPEN;  // Default case for open boundary
    }

    constexpr vector_t<bool>
    operator()(int bc) const
    {
        vector_t<bool> result;

        switch(bc)
        {
            case PS::BOUNDARY_CONDITION_PERIODIC_XYZ:
                result = vector_t<bool>(true, true, true);
                break;
            case PS::BOUNDARY_CONDITION_PERIODIC_XY:
                result = vector_t<bool>{ true, true, false };
                break;
            case PS::BOUNDARY_CONDITION_PERIODIC_XZ:
                result = vector_t<bool>{ true, false, true };
                break;
            case PS::BOUNDARY_CONDITION_PERIODIC_X:
                result = vector_t<bool>{ true, false, false };
                break;
            case PS::BOUNDARY_CONDITION_PERIODIC_YZ:
                result = vector_t<bool>{ false, true, true };
                break;
            case PS::BOUNDARY_CONDITION_PERIODIC_Y:
                result = vector_t<bool>{ false, true, false };
                break;
            case PS::BOUNDARY_CONDITION_PERIODIC_Z:
                result = vector_t<bool>{ false, false, true };
                break;
            case PS::BOUNDARY_CONDITION_OPEN:
                result = vector_t<bool>{ false, false, false };
                break;
            default: throw std::invalid_argument("Unknown boundary condition value.");
        }

        return result;
    }
};
constexpr inline convert_fn convert{};

struct parallel_tag
{
};
constexpr inline parallel_tag parallel{};

struct for_each_fn
{
    template <std::ranges::range Range, class Pred, class Proj = std::identity, class... Args>
    constexpr Pred
    operator()(Range&& range, Pred pred, Proj proj = {}, Args... args) const
    {
        auto first = std::ranges::begin(range);
        auto last = std::ranges::end(range);
        for(; first != last; ++first) std::invoke(pred, std::invoke(proj, *first), args...);
        return pred;
    }

    /// @brief  Parallel version of `util::for_each`.
    /// @tparam Pred Unary function object type. This must NOT have any side
    ///         effects. That is, no element of `range` may be modified other
    ///         than the element itself.
    /// @tparam Range Range type.
    /// @param range Range to be iterated.
    /// @param pred Unary function object to be called for each element.
    /// @return Unary function object.
    template <std::ranges::random_access_range Range, class Pred, class... Args>
    constexpr Pred
    operator()(parallel_tag, Range&& range, Pred pred, Args... args) const
    {
        auto size = std::ranges::size(range);

#pragma omp parallel for
        for(std::size_t i = 0; i < size; ++i) std::invoke(pred, range[i], args...);
        return pred;
    }
};
constexpr inline for_each_fn for_each{};

struct fill_fn
{
    template <std::ranges::range Range, class Value, class Proj = std::identity>
        requires std::indirectly_writable<std::projected<std::ranges::iterator_t<Range>, Proj>,
                                          const Value&>
    constexpr void
    operator()(Range&& range, const Value& value, Proj proj = {}) const
    {
        auto assign = [proj, &value](auto&& elem) { std::invoke(proj, elem) = value; };
        for_each(std::forward<Range>(range), assign);
    }
};
constexpr inline fill_fn fill{};

struct clear_fn
{
    template <std::ranges::range Range, class Proj = std::identity>
    constexpr void
    operator()(Range&& range, Proj proj = {}) const
    {
        using value_type = std::invoke_result_t<Proj, std::ranges::range_value_t<Range>>;
        fill(std::forward<Range>(range), value_type{}, proj);
    }
};
constexpr inline clear_fn clear{};

struct reduce_fn
{
    template <std::ranges::range Range, class Binary, class Proj, class Value>
    constexpr auto
    operator()(Range&& range, Binary&& binary, Proj&& proj, Value&& init) const
    {
        return
            // `binary` must not have any side effects. That is, no element of
            // `range` may be modified other than the element itself.
            std::transform_reduce(
                // std::execution::par_unseq, //
                std::ranges::begin(range),     //
                std::ranges::end(range),       //
                std::forward<Value>(init),     //
                std::forward<Binary>(binary),  //
                [proj](auto&& elem)
                { return std::invoke(proj, std::forward<decltype(elem)>(elem)); }  //
            );
    }
};
constexpr inline reduce_fn reduce{};

struct sum_fn
{
    template <class Tp>
    constexpr Tp
    operator()(const Tp& val) const
    {
        return PS::Comm::getSum(val);
    }
    template <class Tp>
    constexpr vector_t<Tp>
    operator()(const vector_t<Tp>& val) const
    {
        return convert(PS::Comm::getSum(convert(val)));
    }
    template <std::ranges::range Range, class Proj, class Weight = double>
    constexpr auto
    operator()(Range&& range, Proj proj, Weight weight = 1.0) const
    {
        // A member object of `sh::secured` must be casted to a primitive type
        // since it cannot be assigned in usual way. This class wraps a
        // projection functor (e.g. member object pointer, member function
        // pointer, and generic lambda) so that its `operator()` returns a
        // primitive type even if the class has `sh::secured`.
        sh::callable proj_fn{ proj };
        sh::callable weight_fn{ weight };

        using value_type = decltype(proj_fn(std::declval<std::ranges::range_value_t<Range>>()));

        auto weighting = [proj_fn, weight_fn](const auto& elem)
        { return proj_fn(elem) * weight_fn(elem); };

        return operator()(reduce(range, std::plus<>{}, weighting, value_type{}));
    }
};
constexpr inline sum_fn sum{};

struct size_fn
{
    template <std::ranges::range Range>
    constexpr std::size_t
    operator()(Range&& range) const
    {
        return sum(std::ranges::distance(std::forward<Range>(range)));
    }
};
constexpr size_fn size{};

struct empty_fn
{
    template <std::ranges::range Range>
    constexpr bool
    operator()(Range&& range) const
    {
        return size(std::forward<Range>(range)) == 0;
    }
};
constexpr empty_fn empty{};

struct mean_fn
{
    template <std::ranges::range Range, class Proj = std::identity, class Weight = double>
    constexpr auto
    operator()(Range&& range, Proj proj = {}, Weight weight = 1.0) const
    {
        return sum(range, proj, weight) / sum(range, weight);
    }
};
constexpr inline mean_fn mean{};

struct find_fn
{
    template <std::ranges::range Range, class Value, class Proj = std::identity>
    constexpr std::optional<std::ranges::range_value_t<Range>>
    operator()(Range&& range, const Value& value, Proj proj = {}) const
    {
        // Do not forward `range` to avoid dangling
        auto it = std::ranges::find(range, value, proj);
        auto end = std::ranges::end(range);
        if(it == end) return {};
        return *it;
    }
};
constexpr inline find_fn find{};

struct find_all_fn
{
    template <std::ranges::range Range, class Value, class Proj = std::identity>
    constexpr std::vector<std::ranges::range_value_t<Range>>
    operator()(Range&& range, const Value& value, Proj proj = {}) const
    {
        using Element = std::ranges::range_value_t<Range>;

        std::vector<Element> local_results;
        for(auto&& elem : range)
            if(value == std::invoke(proj, elem)) local_results.push_back(elem);

        return par::all_gather(local_results);
    }
};
constexpr inline find_all_fn find_all{};

struct find_if_fn
{
    template <std::ranges::range Range, class Pred, class Proj = std::identity>
    constexpr std::optional<std::ranges::range_value_t<Range>>
    operator()(Range&& range, Pred pred, Proj proj = {}) const
    {
        auto it = std::ranges::find_if(range, pred, proj);
        auto end = std::ranges::end(range);
        if(it == end) return {};
        return *it;
    }
};
constexpr inline find_if_fn find_if{};

struct find_if_all_fn
{
    template <std::ranges::range Range, class Pred, class Proj = std::identity>
    constexpr std::vector<std::ranges::range_value_t<Range>>
    operator()(Range&& range, Pred pred, Proj proj = {}) const
    {
        using Element = std::ranges::range_value_t<Range>;

        std::vector<Element> local_results;
        for(auto&& elem : range)
            if(std::invoke(pred, std::invoke(proj, elem))) local_results.push_back(elem);

        return par::all_gather(local_results);
    }
};
constexpr inline find_if_all_fn find_if_all{};

struct min_fn
{
    template <class Tp>
    constexpr Tp
    operator()(const Tp& val) const
    {
        return PS::Comm::getMinValue(val);
    }
    template <class Tp>
    constexpr vector_t<Tp>
    operator()(const vector_t<Tp>& val) const
    {
        return convert(PS::Comm::getMinValue(convert(val)));
    }
    template <std::ranges::range R, class Comp = std::ranges::less, class Proj = std::identity>
    constexpr decltype(auto)
    operator()(R&& range, Comp comp = {}, Proj proj = {}) const
    {
        auto first = std::ranges::begin(range);
        auto last = std::ranges::end(range);

        if(first == last) throw sh::invalid_value{} << "Cannot compute min on an empty range";

        auto it = std::ranges::min_element(first, last, comp, proj);

        assert(it != last);
        return (*this)(std::invoke(proj, *it));
    }
};
constexpr inline min_fn min{};

struct max_fn
{
    template <class Tp>
    constexpr Tp
    operator()(const Tp& val) const
    {
        return PS::Comm::getMaxValue(val);
    }
    template <class Tp>
    constexpr vector_t<Tp>
    operator()(const vector_t<Tp>& val) const
    {
        return convert(PS::Comm::getMaxValue(convert(val)));
    }
    template <std::ranges::range R, class Comp = std::ranges::less, class Proj = std::identity>
    constexpr decltype(auto)
    operator()(R&& range, Comp comp = {}, Proj proj = {}) const
    {
        auto first = std::ranges::begin(range);
        auto last = std::ranges::end(range);

        if(first == last)
            throw sh::invalid_argument{} << "Cannot compute max on an empty range";

        auto it = std::ranges::max_element(first, last, comp, proj);

        assert(it != last);
        return (*this)(std::invoke(proj, *it));
    }
};
constexpr inline max_fn max{};

struct min_element_fn
{
    template <std::ranges::forward_range Range, class Proj = std::identity,
              class Comp = std::ranges::less>
        requires std::indirect_strict_weak_order<
            Comp, std::projected<std::ranges::iterator_t<Range>, Proj>>
    [[nodiscard]]
    constexpr auto
    operator()(Range&& range, Comp comp = {}, Proj proj = {}) const
    {
        if(empty(range))
        {
            throw sh::invalid_argument{} << "Cannot compute min on an empty range";
        }

        auto first = std::ranges::begin(range);
        auto last = std::ranges::end(range);

        // find local min iterator
        auto it = std::ranges::min_element(first, last, comp, proj);
        assert(it != last);

        using elm_t = std::ranges::range_value_t<Range>;
        using val_t = std::invoke_result_t<Proj, elm_t>;

        // materialize local element + its projected value
        elm_t local_elm = *it;
        val_t local_val = std::invoke(proj, local_elm);

        // --- Gather projected values to root and pick the global min + its rank ---
        val_t global_val{};
        int winner_rank{};
        int tag = 10;

        if(rank == root_rank)
        {
            bool first_recv = true;
            for(int i = 0; i < n_procs; ++i)
            {
                std::optional<val_t> tmp;
                if(i == root_rank)
                {
                    tmp = local_val;
                }
                else
                {
                    PS::Comm::recv(&tmp, 1, i, tag);
                }
                // pick the smaller value
                if(tmp && (first_recv || std::invoke(comp, *tmp, global_val)))
                {
                    global_val = *tmp;
                    winner_rank = i;
                    first_recv = false;
                }
            }
            assert(!first_recv);  // we must have seen at least one value
        }
        else
        {
            PS::Comm::send(&local_val, 1, root_rank, tag);
        }

        // tell everyone who had the min
        PS::Comm::broadcast(&winner_rank, 1, root_rank);

        // --- Broadcast the element itself from the winner process ---
        elm_t global_elm{};
        if(rank == winner_rank)
        {
            global_elm = local_elm;
        }
        PS::Comm::broadcast(&global_elm, 1, winner_rank);

        return global_elm;
    }
};
inline constexpr min_element_fn min_element{};

struct max_element_fn
{
    template <std::ranges::forward_range Range, class Proj = std::identity,
              class Comp = std::ranges::less>
        requires std::indirect_strict_weak_order<
            Comp, std::projected<std::ranges::iterator_t<Range>, Proj>>
    [[nodiscard]]
    constexpr auto
    operator()(Range&& range, Comp comp = {}, Proj proj = {}) const
    {
        if(empty(range))
        {
            throw sh::invalid_argument{} << "Cannot compute max on an empty range";
        }

        auto first = std::ranges::begin(range);
        auto last = std::ranges::end(range);

        // find local max iterator
        auto it = std::ranges::max_element(first, last, comp, proj);
        assert(it != last);

        using elm_t = std::ranges::range_value_t<Range>;
        using val_t = std::invoke_result_t<Proj, elm_t>;

        // materialize local element + its projected value
        elm_t local_elm = *it;
        val_t local_val = std::invoke(proj, local_elm);

        // --- Gather projected values to root and pick the global max + its rank ---
        val_t global_val{};
        int winner_rank{};
        int tag = 20;

        if(rank == root_rank)
        {
            bool first_recv = true;
            for(int i = 0; i < n_procs; ++i)
            {
                std::optional<val_t> tmp;
                if(i == root_rank)
                {
                    tmp = local_val;
                }
                else
                {
                    PS::Comm::recv(&tmp, 1, i, tag);
                }
                if(tmp && (first_recv || std::invoke(comp, global_val, *tmp)))
                {
                    global_val = *tmp;
                    winner_rank = i;
                    first_recv = false;
                }
            }
            assert(!first_recv);  // we must have seen at least one value
        }
        else
        {
            PS::Comm::send(&local_val, 1, root_rank, tag);
        }

        // tell everyone who had the max
        PS::Comm::broadcast(&winner_rank, 1, root_rank);

        // --- Broadcast the element itself from the winner process ---
        elm_t global_elm{};
        if(rank == winner_rank)
        {
            global_elm = local_elm;
        }
        PS::Comm::broadcast(&global_elm, 1, winner_rank);

        return global_elm;
    }
};
inline constexpr max_element_fn max_element{};

struct minmax_fn
{
    template <class Tp>
    struct result_type
    {
        Tp min, max;
    };

    template <typename T, typename Comp = std::ranges::less>
    constexpr result_type<T>
    operator()(const T& local_value, Comp comp = {}) const
    {
        T min_val_loc = local_value;
        T max_val_loc = local_value;
        std::size_t min_val_rank = rank;
        std::size_t max_val_rank = rank;

        T min_val_glb = min_val_loc;
        T max_val_glb = max_val_loc;

        if(rank == root_rank)
        {
            for(int i = 0; i < n_procs; ++i)
            {
                if(i == root_rank) continue;

                T min_rcv, max_rcv;
                PS::Comm::recv(&min_rcv, 1, i, 30);
                PS::Comm::recv(&max_rcv, 1, i, 31);

                if(std::invoke(comp, min_rcv, min_val_glb))
                {
                    min_val_glb = min_rcv;
                    min_val_rank = i;
                }
                if(std::invoke(comp, max_val_glb, max_rcv))
                {
                    max_val_glb = max_rcv;
                    max_val_rank = i;
                }
            }
        }
        else
        {
            PS::Comm::send(&min_val_loc, 1, root_rank, 30);
            PS::Comm::send(&max_val_loc, 1, root_rank, 31);
        }

        par::broadcast(min_val_rank, root_rank);
        par::broadcast(max_val_rank, root_rank);

        // Ensure values are properly set before broadcasting
        if(rank == root_rank)
        {
            min_val_glb = min_val_loc;
            max_val_glb = max_val_loc;
        }

        par::broadcast(min_val_glb, min_val_rank);
        par::broadcast(max_val_glb, max_val_rank);

        return { min_val_glb, max_val_glb };
    }

    template <std::ranges::range R, class Comp = std::ranges::less, class Proj = std::identity>
    constexpr result_type<std::ranges::range_value_t<R>>
    operator()(R&& range, Comp comp = {}, Proj proj = {}) const
    {
        using T = std::ranges::range_value_t<R>;

        if(empty(range)) throw sh::invalid_value{} << "Cannot compute minmax on an empty range";

        auto [min_it, max_it] = std::ranges::minmax_element(range, comp, proj);

        return { min(std::invoke(proj, *min_it)), max(std::invoke(proj, *max_it)) };
    }
};

constexpr inline minmax_fn minmax{};
struct product_fn
{
    template <std::ranges::range Range>
    constexpr auto
    operator()(const Range& range, std::optional<std::size_t> size = {}) const
    {
        using T = std::ranges::range_value_t<Range>;

        if(not size) size = util::size(range);

        return std::accumulate(range.begin(), range.begin() + *size, T{ 1 },
                               std::multiplies<T>{});
    }
};
constexpr inline product_fn product{};

struct project_fn
{
    template <std::ranges::viewable_range Range>
    constexpr void
    operator()(Range& range, std::size_t n_dim,
               std::ranges::range_value_t<Range> init = {}) const
    {
        for(auto& e : range | std::ranges::views::drop(n_dim)) e = init;
    }
};

struct volume_fn
{
    template <class Tp>
    constexpr double
    operator()(const vector_t<Tp>& diag, std::size_t n_dim) const
    {
        return diag.prod_n(n_dim);
    }
    template <class Tp>
    constexpr double
    operator()(const orthotope_t<Tp>& box, std::size_t n_dim) const
    {
        return operator()(box.diag(), n_dim);
    }
};
constexpr inline volume_fn volume{};

struct write_fn
{
    void
    check(const std::ostream& os) const
    {
        if(!os) throw sh::runtime_error{} << "Given output stream is in bad condition.";
    }
    template <std::ranges::range Range>
    void
    operator()(std::string& path, Range&& range, const file_header& header) const
    {
        std::ofstream file;
        if(in_root_proc)
        {
            file.open(path);
            if(!file) throw sh::runtime_error{} << "Cannot open the file with path: " << path;
        }
        (*this)(file, std::forward<Range>(range), header);
    }
    template <std::ranges::range Range>
    void
    operator()(std::ostream& os, Range&& range, const file_header& header) const
    {
        check(os);

        if(in_root_proc)
        {
            os << header;
            os << std::ranges::range_value_t<Range>::header << '\n';
        }

        std::stringstream ss;
        ss << std::scientific << std::showpos << std::setprecision(12);

        auto write = [&ss](const auto& e) { ss << e << '\n'; };
        for_each(std::forward<Range>(range), write);

        (*this)(os, ss);
    }

    void
    operator()(std::ostream& os, const std::stringstream& ss) const
    {
        check(os);

        auto data_loc = ss.str();
        auto size_loc = static_cast<int>(data_loc.size());
        // Gather sizes from all processes
        std::vector<int> recv_counts(n_procs);
        PS::Comm::gather(&size_loc, 1, recv_counts.data());
        // Calculate displacements for each process
        std::vector<int> displs(n_procs, 0);
        if(in_root_proc)
        {
            int offset = 0;
            for(auto i : std::views::iota(0UL, n_procs))
            {
                displs[i] = offset;
                offset += recv_counts[i];
            }
        }
        // Gather all data at root process
        std::vector<char> data_glob;
        if(in_root_proc)
        {
            data_glob.resize(displs.back() + recv_counts.back());
        }
        PS::Comm::gatherV(data_loc.data(), size_loc, data_glob.data(), recv_counts.data(),
                          displs.data(), 0);
        // Root process writes to file
        if(in_root_proc)
        {
            os.write(data_glob.data(), data_glob.size());
        }
    }
};
constexpr inline write_fn write{};

struct read_fn
{
    void
    check(const std::istream& is) const
    {
        if(!is) throw sh::exception{} << "Given input stream is in bad condition.";
    }
    template <std::ranges::range Range>
    void
    operator()(std::string& path, Range& range, file_header& header) const
    {
        std::ifstream file;
        if(in_root_proc)
        {
            file.open(path);
            if(!file) throw sh::runtime_error{} << "Cannot open the file with path: " << path;
        }
        (*this)(file, range, header);
    }
    template <std::ranges::range Range>
    void
    operator()(std::istream& is, Range& range, file_header& header) const
    {
        if(in_root_proc) check(is);

        using ptcl_type = std::ranges::range_value_t<Range>;

        // Root process reads the header.
        if(in_root_proc)
        {
            is >> header >> std::ws;  // Read the file header.

            std::vector<std::string> read_ptcl_header;
            std::string line;
            std::getline(is, line);
            std::istringstream iss(line);
            std::string key;
            while(iss >> key) read_ptcl_header.push_back(key);

            if(read_ptcl_header != ptcl_type::header_v)
                throw sh::runtime_error{}
                    << "The particle data of the read file is inconsistent"
                       " with that of the container class.";
        }

        // Broadcast the header to all processes.
        PS::Comm::broadcast(&header, 1, root_rank);

        // Vector for storing data.
        std::vector<ptcl_type> send_buf;

        if(in_root_proc)
        {
            // Root process reads data.
            ptcl_type val;
            while(is >> val) send_buf.push_back(val);
        }

        // Calculate and broadcast the size of data each process will receive.
        std::size_t n_send = 0;
        if(in_root_proc) n_send = send_buf.size() / n_procs;
        PS::Comm::broadcast(&n_send, 1, root_rank);

        // Resize the receiving range.
        range.resize(n_send);

        // Scatter the data to all processes.
        PS::Comm::scatter(send_buf.data(), n_send, range.data(), root_rank);

        // Handle remainder data for root process.
        if(in_root_proc)
        {
            auto beg = send_buf.begin() + n_send * n_procs;
            range.insert(range.end(), beg, send_buf.end());
        }
    }
    void
    operator()(std::istream& is, std::stringstream& ss) const
    {
        if(in_root_proc) check(is);

        // Read the entire file content into a stringstream on rank 0
        std::string data_glob;
        if(in_root_proc)
        {
            is.seekg(0, std::ios::end);
            std::streamsize size = is.tellg();
            data_glob.resize(size);
            is.seekg(0, std::ios::beg);
            is.read(&data_glob[0], size);
        }

        // Broadcast the size of the file data to all processes
        std::vector<int> recv_counts(n_procs);
        int size_loc = 0;
        if(in_root_proc)
        {
            size_loc = static_cast<int>(data_glob.size());
        }
        PS::Comm::broadcast(&size_loc, 1, 0);

        // Split the file content into chunks for each process
        std::vector<int> displs(n_procs, 0);
        if(in_root_proc)
        {
            int offset = 0;
            for(int i = 0; i < n_procs; ++i)
            {
                recv_counts[i] = size_loc / n_procs;
                displs[i] = offset;
                offset += recv_counts[i];
            }
        }

        // Allocate buffer for receiving local data
        std::string data_loc(recv_counts[rank], '\0');

        // Use sendIrecvV to distribute data among processes
        PS::Comm::sendIrecvV(data_glob.data(), nullptr, recv_counts.data(), displs.data(),
                             n_procs, data_loc.data(), nullptr, &recv_counts[rank],
                             &displs[rank], n_procs);

        // Store the local data into the stringstream
        ss.str(data_loc);
    }
};
constexpr inline read_fn read{};

struct make_random_shift_fn
{
    template <std::ranges::range R, class Proj = std::identity>
    // requires(std::assignable_from<std::invoke_result_t<std::ranges::range_value_t<R&&>,
    // Proj>,
    //                               vector>)
    void
    operator()(R&& range, Proj proj, double magnitude) const
    {
        static std::random_device rd;
        static std::mt19937 gen{ rd() };
        static std::normal_distribution<double> dis_nor{ 0, 1 };
        static std::uniform_real_distribution<double> dis_uni{ 0, 1 };

        auto shift_fn = [&](auto& p)
        {
            auto shift_vec = vector::unit();
            shift_vec.normalize(magnitude * dis_nor(gen));
            shift_vec.rotate(M_PI * dis_uni(gen), M_PI * dis_uni(gen));
            std::invoke(proj, p) += shift_vec;
        };
        util::for_each(range, shift_fn);
    }
};
constexpr inline make_random_shift_fn make_random_shift{};

struct arrangement
{
    std::size_t n_ptcl;
    orthotope domain;
    double spacing;
};

struct box_info_fn
{
   private:
    template <std::ranges::range Range>
    orthotope
    calc_tight_box(Range&& range, std::size_t n_dim) const
    {
        orthotope box{ vector::constant(std::numeric_limits<double>::max()),
                       vector::constant(std::numeric_limits<double>::min()) };

        for(auto d = 0UL; d < n_dim; ++d)
        {
            if(std::ranges::begin(range) != std::ranges::end(range))
            {  // `std::ranges::empty` is deleted for `filter_view`
                const auto [min, max] = std::minmax_element(
                    std::ranges::begin(range), std::ranges::end(range),
                    [d](const auto& a, const auto& b) { return a.pos[d] < b.pos[d]; });
                box.lower()[d] = min->pos[d];
                box.upper()[d] = max->pos[d];
            }
            box.lower()[d] = min(box.lower()[d]);
            box.upper()[d] = max(box.upper()[d]);
        }

        return box.truncate(n_dim);
    }

    // Shrunk the `box` so that it can be filled exactly homogeneously with
    // elements aligned with `spacing` in the periodic boundary condition, i.e.
    // the length of the `box` edge is to be an integer multiple of `spacing`.
    void
    shrunk_box_for_fitting(orthotope& box, double spacing, std::size_t n_dim) const
    {
        const auto diag = box.diag();

        vector shift;

        // Use a custom function since `std::fmod` has bad behaviour. e.g.
        // `std::fmod( 1, 0.1 )` returns `0.1` instead of `0`.
        auto mod = [](double x, double y) { return x - std::floor(x / y) * y; };

        for(auto i : std::ranges::views::iota(0UL, n_dim)) shift[i] = -mod(diag[i], spacing);
        // cerr << std::setprecision( 16 );
        // cerr << "box     : " << box << '\n';
        // cerr << "diag    : " << diag << '\n';
        // cerr << "spacing : " << spacing << '\n';
        // cerr << "fmod    : " << -shift << '\n';

        box.extend(shift);
    }

    double
    calc_spacing_in_tight_box(const orthotope& tight_box, std::size_t size,
                              std::size_t n_dim) const
    {
        if(size < 2UL)
            throw std::invalid_argument{ "Provided size must be or greater than 2." };

        auto l = tight_box.diag().make_truncated(n_dim);

        if(n_dim == 1)
        {
            return l[0] / (size - 1);
        }
        else if(n_dim == 2)
        {
            return 0.5 / (size - 1) *
                   (l[0] + l[1] +
                    std::sqrt(square(l[0] + l[1]) + 4 * (size - 1) * l[0] * l[1]));
        }
        else if(n_dim == 3)
        {
            return [](double a, double b, double c, double d)
            {
                b /= a;
                c /= a;
                d /= a;
                auto p = c - (b * b) / 3.0;
                auto q = (2.0 * b * b * b) / 27.0 - (b * c) / 3.0 + d;
                auto shift = b / 3.0;
                auto discriminant = std::pow(q / 2, 2) + std::pow(p / 3, 3);
                if(discriminant <= 0.0)
                    throw sh::exception{} << "The equation does not have one real root and two "
                                             "complex conjugate roots.";
                auto sqrt_discriminant = std::sqrt(discriminant);
                auto u = std::cbrt(-q / 2.0 + sqrt_discriminant);
                auto v = std::cbrt(-q / 2.0 - sqrt_discriminant);
                auto t = u + v;
                return t - shift;
            }(size - 1, -l[0] - l[1] - l[2], -l[0] * l[1] - l[1] * l[2] - l[2] * l[0],
              -l[0] * l[1] * l[2]);
        }
        else
        {
            throw sh::invalid_argument{} << "Invalid dimension: " << n_dim;
        }
    }

    double
    calc_spacing_in_box_with_margin(const orthotope& box_with_margin, std::size_t size,
                                    std::size_t n_dim, vector unit_side_ratio) const
    {
        return std::pow(volume(box_with_margin, n_dim) / size, 1. / n_dim);
    }

    vector_t<std::size_t>
    n_elms_along(const orthotope& box, double spacing, std::size_t n_dim) const
    {
        vector_t<std::size_t> n_along;
        for(auto i : std::ranges::views::iota(0UL, n_dim))
            n_along[i] = std::trunc(box.diag()[i] / spacing);

        // cerr << "org  : " << box.diagonal()[ 0 ] / spacing << '\n';
        // cerr << "round: " << std::round( box.diagonal()[ 0 ] / spacing ) <<
        // '\n'; cerr << "trunc: " << std::trunc( box.diagonal()[ 0 ] / spacing )
        // << '\n';

        return n_along.truncate(n_dim, 1UL);
    }

   public:
    /// @brief Assume the `box` has margin. Given `box` and `size` might be
    /// modified.
    /// @return `first`: mean element spacing
    /// @return `second`: # of elements along axes
    std::pair<double, vector_t<std::size_t>>
    operator()(orthotope& box, std::size_t& size, std::size_t n_dim,
               vector unit_side_ratio) const
    {
        box.truncate(n_dim);

        auto spacing = calc_spacing_in_box_with_margin(box, size, n_dim, unit_side_ratio);
        auto n_along = n_elms_along(box, spacing, n_dim);

        size = product(n_along, n_dim);

        // cerr << "box: " << box << '\n';
        shrunk_box_for_fitting(box, spacing, n_dim);

        // cerr << "box.diagonal() / spacing = " << box.diagonal() / spacing <<
        // '\n';

        return { spacing, n_along };
    }

    /// @brief Calculate the box size and mean spacing from given `range`, in
    /// `n_dim`. The resultant box is not tight, rather includes a margin, and
    /// can be used as is for periodic boundary.
    /// @return `first`: box size
    /// @return `second`: mean spacing of the elements
    template <std::ranges::range Range>
    arrangement
    operator()(Range&& range, std::size_t n_dim) const
    {
        auto n = size(range);
        auto tight_box = calc_tight_box(range, n_dim);
        auto spacing = calc_spacing_in_tight_box(tight_box, n, n_dim);
        // Add margin for the case of periodic boundary.
        auto box = tight_box.make_extended(vector(spacing).truncate(n_dim));

        return { n, box, spacing };
    }
};
constexpr inline box_info_fn box_info{};

}  // namespace util
}  // namespace cfe