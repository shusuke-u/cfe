#pragma once

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <vector>

#include <force.h>

namespace cfe
{

namespace particles
{

template <class T_f>
concept particle = std::is_copy_assignable_v<T_f> && std::is_copy_constructible_v<T_f> && requires(T_f p) {
    { p.mass() } -> std::same_as<PS::F64>;
    { p.pos() } -> std::same_as<vector>;
    { p.vel() } -> std::same_as<vector>;
    { p.mom() } -> std::same_as<vector>;

    { p.set_pos(std::declval<vector>()) };
    { p.set_vel(std::declval<vector>()) };
    { p.add_pos(std::declval<vector>()) };
    { p.add_vel(std::declval<vector>()) };

    { p.belongs_to(std::declval<PS::U32>()) } -> std::same_as<bool>;
    { p.make_belong_to(std::declval<PS::U32>()) };
    { p.make_leave_from(std::declval<PS::U32>()) };
};

class basic_index
{
   protected:
    PS::U32 _index{};

   public:
    void set_index(PS::S32 index);

    PS::S32 index() const;
};

class basic_pos
{
   protected:
    vector _pos{};

   public:
    void set_pos(const vector src);

    void add_pos(const vector src);

    /// @return The position of `*this`.
    vector pos() const;

    /// @brief APIs for FDPS.
    PS::F64vec getPos() const;
    void setPos(const vector src);
};

class _base_mass
{
   protected:
    PS::F64 _mass{};

   public:
    PS::F64 mass() const;

    void set_mass(PS::F64 src);
};

class _base_r_search
{
   protected:
    PS::F64 _r_search{};

   public:
    /// @brief APIs for FDPS
    PS::F64 getRSearch() const;
};

namespace spring
{
struct coef
{
    PS::F64 att, rep;
};
class J : public virtual basic_index, public virtual basic_pos, public _base_r_search
{
   private:
    vector _init_pos;

   public:
    vector init_pos() const;

    void set_init_pos();

    void copyFromFP(const J& fp);
};
class I : public J
{
   private:
    coef _coef;

   public:
    I() : J{}
    {
    }

    I&
    operator=(const I& src)
    {
        if(this == &src) return *this;
        if(&src == nullptr)
        {
            std::cerr << "nullptr at spring::I::operator=\n";
            PS::Abort();
        }
        J::operator=(src);
        _coef = src._coef;
        return *this;
    }

    coef spring_coef() const;

    void set_spring_coef(coef src);

    void set_spr_len_lim(PS::F64 src);

    /// @brief APIs for FDPS
    void copyFromFP(const I& fp);
};
}  // namespace spring

class density : public virtual _base_mass, public virtual basic_pos, public _base_r_search
{  // Density should be the same for i and j. `r_search` is not virtual.
   private:
    PS::F64 _smth;

   public:
    density&
    operator=(const density& src)
    {
        if(this == &src) return *this;
        _base_mass::operator=(src);
        basic_pos::operator=(src);
        _base_r_search::operator=(src);
        _smth = src._smth;
        return *this;
    }

    PS::F64 smth() const;

    void set_smth_for_dens(PS::F64 smoothing_length, PS::F64 kernel_support_ratio);

    void copyFromFP(const density& src);
};

class full_ptcl : public spring::I, public density
{
   private:
    bool _is_ghost{ false };

    basic_vector<bool> _is_fixed{ false };

    std::array<bool, n_max_views> _belongs_to;

    PS::S32 _material_index;

    vector _vel;

    vector _force;

    PS::F64 _conserved_energy, _pot;

    PS::F64 _dens;

   public:
    full_ptcl() : spring::I{}, density{}
    {
    }

    basic_vector<bool> is_fixed() const;

    bool belongs_to(PS::U32 index) const;

    vector vel() const;
    vector acc() const;
    vector mom() const;
    PS::S32 material_index() const;
    PS::F64 spring_pot() const;

    PS::F64 internal_eng() const;
    PS::F64 kinetic_eng() const;

    void make_belong_to(PS::U32 index);
    void make_leave_from(PS::U32 index);

    void fix(basic_vector<bool> direction = true);
    void unfix(basic_vector<bool> direction = true);

    void set_vel(vector vel);
    void set_acc(vector acc);
    void set_pot(PS::F64 pot);
    void set_internal_eng(PS::F64 internal_eng);
    void set_material_index(PS::S32 index);

    void add_vel(vector vel);
    void add_pot(PS::F64 pot);
    void add_force(vector src);

    void clear_dens();
    void clear_force();

    bool is_ghost() const;

    friend std::ostream& operator<<(std::ostream& os, const full_ptcl& rhs);

    // APIs for FDPS
    void writeAscii(FILE* fp) const;
    void readAscii(FILE* fp);

    void copyFromForce(const interaction::result::spring& spr);
    void copyFromForce(const interaction::result::density& dens);
};

}  // namespace particles
}  // namespace cfe
