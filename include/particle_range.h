#pragma once

#include <particle_simulator.hpp>

namespace cfe
{

template <class Ptcl>
class particle_range : public ParticleSimulator::ParticleSystem<Ptcl>
{
   public:
    using base = ParticleSimulator::ParticleSystem<Ptcl>;

    using value_type = Ptcl;
    using iterator = Ptcl*;
    using const_iterator = const Ptcl*;

    iterator
    begin()
    {
        return &(*this)[0];
    }
    iterator
    end()
    {
        return &(*this)[this->getNumberOfParticleLocal()];
    }

    const_iterator
    begin() const
    {
        return &(*this)[0];
    }
    const_iterator
    end() const
    {
        return &(*this)[this->getNumberOfParticleLocal()];
    }

    const_iterator
    cbegin() const
    {
        return &(*this)[0];
    }
    const_iterator
    cend() const
    {
        return &(*this)[this->getNumberOfParticleLocal()];
    }

    value_type&
    front()
    {
        return (*this)[0];
    }
    value_type&
    back()
    {
        return (*this)[this->getNumberOfParticleLocal() - 1];
    }

    /// @brief LOCAL size of `*this`.
    std::size_t
    size() const
    {
        return this->getNumberOfParticleLocal();
    }

    bool
    empty() const
    {
        return is_locally_empty();
    }

    bool
    is_locally_empty() const
    {
        return this->getNumberOfParticleLocal() == 0;
    }
    bool
    is_globally_empty() const
    {
        return this->getNumberOfParticleGlobal() == 0;
    }
};

}  // namespace cfe