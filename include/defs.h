#pragma once

#include <limits>
#include <ranges>
#include <string>
#include <sstream>
#include <filesystem>

#include <omp.h>

#include <particle_simulator.hpp>

#include <sh/orthotope.h>
#include <sh/exceptions.h>

namespace cfe
{

constexpr inline auto null_index = std::numeric_limits<PS::S32>::max();

inline bool in_root_proc;
inline std::size_t rank, n_procs;
inline constexpr std::size_t root_rank{ 0 };

inline void
initialize(int argc, char* argv[])
{
    static bool called = false;

    PS::Initialize(argc, argv);
    in_root_proc = (PS::Comm::getRank() == 0);
    rank = PS::Comm::getRank();
    n_procs = PS::Comm::getNumberOfProc();

    // std::set_terminate(
    //     []()
    //     {
    //         std::cerr << "Terminated due to uncaught exception. Aborting all processes...\n";
    //         PS::Abort();
    //     } );

    std::cout << std::scientific << std::showpos;
    std::cerr << std::scientific << std::showpos;

    called = true;
}

inline void
finalize()
{
    PS::Finalize();
}

inline void
abort(int err = -1)
{
    PS::Abort(err);
}

struct file_header
{
   public:
    std::size_t n_ptcl;
    double time;

    friend std::ostream&
    operator<<(std::ostream& os, const file_header& src)
    {
        using sh::del;

        os << std::scientific << std::showpos;
        os << "n_ptcl" << del;
        os << "time" << del;
        os << '\n';
        os << src.n_ptcl << del;
        os << src.time << del;
        os << '\n';

        return os;
    }

    friend std::istream&
    operator>>(std::istream& is, file_header& tar)
    {
        is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // skip header
        return is >> tar.n_ptcl >> tar.time >> std::ws;
    }

    size_t
    readAscii(FILE* fp)
    {
        fscanf(fp, "%zd\n%lf\n", &n_ptcl, &time);  // Must be read in one function call.
        return n_ptcl;
    }

    void
    writeAscii(FILE* fp) const
    {
        fprintf(fp, "%zd\n", n_ptcl);
        fprintf(fp, "%lf\n", time);
    }

    size_t
    readBinary(FILE* fp)
    {
        fread(&n_ptcl, sizeof(decltype(n_ptcl)), 1, fp);
        fread(&time, sizeof(decltype(time)), 1, fp);
        return n_ptcl;
    }

    void
    writeBinary(FILE* fp) const
    {
        fwrite(&n_ptcl, sizeof(decltype(n_ptcl)), 1, fp);
        fwrite(&time, sizeof(decltype(time)), 1, fp);
    }
};

class ostream : public std::ostream
{
   public:
    ostream(std::ostream& stream) : std::ostream(&_buffer), _buffer(stream)
    {
        *this << std::scientific << std::showpos;
    }

   private:
    class streambuf : public std::streambuf
    {
       public:
        streambuf(std::ostream& stream) : _buf(stream.rdbuf())
        {
        }

       protected:
        virtual int
        overflow(int c) override
        {
            if(PS::Comm::getRank() == 0)
            {
                if(c != EOF)
                {
                    return _buf->sputc(c);
                }
            }
            return c;
        }

        virtual std::streamsize
        xsputn(const char* s, std::streamsize n) override
        {
            if(PS::Comm::getRank() == 0)
            {
                return _buf->sputn(s, n);
            }
            return n;
        }

       private:
        std::streambuf* _buf;
    };

    streambuf _buffer;
};

inline ostream cout{ std::cout };
inline ostream cerr{ std::cerr };

// inline void
// clear_dir( const std::string& path_to_dir )
// {
//     if( not in_root_proc ) return;
//     using namespace std::filesystem;
//     for( auto& entry :
//          directory_iterator( path_to_dir ) | std::views::filter( []( const auto& etr ) {
//          return is_regular_file( etr ); } ) )
//     {
//         if( not remove( entry.path() ) ) std::cerr << "Cannot delete file: " << entry.path()
//         << '\n';
//     }
// }

inline void
clear_dir(const std::string& path_to_dir)
{
    if(not in_root_proc) return;

    using namespace std::ranges;
    namespace fs = std::filesystem;

    // Wrap the directory_iterator in a subrange to make it a viewable range
    auto directory = subrange(fs::directory_iterator(path_to_dir), fs::directory_iterator());

    // Now you can use the ranges pipeline with operator|
    for(const auto& entry :
        directory | views::filter([](const auto& e) { return fs::is_regular_file(e); }))
    {
        if(not fs::remove(entry.path()))
            std::cerr << "Cannot delete file: " << entry.path() << '\n';
    }
}

inline constexpr auto n_max_views{ 16ul };
// inline constexpr auto n_max_neighbors{ 64ul };

#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
inline constexpr auto size_of_vec{ 2UL };
template <class T = double>
using vector_t = sh::vector<T, 2>;
template <class T = double>
using orthotope_t = sh::orthotope<T, 2>;
#else
inline constexpr auto vector_size{ 3ul };
template <class Tp>
using vector_t = sh::vector<Tp, 3ul>;
template <class Tp>
using orthotope_t = sh::orthotope<Tp, 3ul>;
#endif

using vector = vector_t<double>;
using orthotope = orthotope_t<double>;

}  // namespace cfe