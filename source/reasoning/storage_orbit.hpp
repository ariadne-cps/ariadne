// Grid-backed storage orbits belong to reasoning, not evolution.
#ifndef ARIADNE_REASONING_STORAGE_ORBIT_HPP
#define ARIADNE_REASONING_STORAGE_ORBIT_HPP

#include "evolution/orbit.hpp"
#include "reasoning/storage.hpp"

namespace Ariadne {

template<>
class Orbit<Storage>
{
    struct Data;
  public:
    typedef Storage EnclosureListType;

    Orbit(const Storage& initial);
    Orbit(const Storage& initial, const Storage& reach,
          const Storage& intermediate, const Storage& final);

    Storage const& initial() const;
    Storage const& reach() const;
    Storage const& intermediate() const;
    Storage const& final() const;
  private:
    std::shared_ptr<Data> _data;
};

struct Orbit<Storage>::Data {
    Data(const Grid& grid, const EffectiveVectorMultivariateFunction& auxiliary_mapping)
            : initial(grid,auxiliary_mapping), reach(initial), intermediate(initial), final(initial) { }
    Storage initial;
    Storage reach;
    Storage intermediate;
    Storage final;
};

} // namespace Ariadne

#endif // ARIADNE_REASONING_STORAGE_ORBIT_HPP
