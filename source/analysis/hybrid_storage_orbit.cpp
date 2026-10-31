#include "analysis/hybrid_storage_orbit.hpp"

namespace Ariadne {

struct Orbit<HybridStorage>::Data {
    Data(const HybridGrid& grid, HybridAutomatonInterface const& sys)
        : initial(grid,sys), reach(grid,sys), intermediate(grid,sys), final(grid,sys) { }
    HybridStorage initial;
    HybridStorage reach;
    HybridStorage intermediate;
    HybridStorage final;
};

Orbit<HybridStorage>::
Orbit(const HybridStorage& initial_set)
    : _data(new Data(initial_set.grid(),initial_set.system()))
{
    this->_data->initial=initial_set;
}

Orbit<HybridStorage>::
Orbit(const HybridStorage& initial_set,
      const HybridStorage& reach_set,
      const HybridStorage& intermediate_set,
      const HybridStorage& final_set)
    : _data(new Data(initial_set.grid(),initial_set.system()))
{
    this->_data->initial=initial_set;
    this->_data->reach=reach_set;
    this->_data->intermediate=intermediate_set;
    this->_data->final=final_set;
}

HybridStorage const&
Orbit<HybridStorage>::
initial() const
{
    return this->_data->initial;
}

HybridStorage const&
Orbit<HybridStorage>::
reach() const
{
    return this->_data->reach;
}

HybridStorage const&
Orbit<HybridStorage>::
intermediate() const
{
    return this->_data->intermediate;
}

HybridStorage const&
Orbit<HybridStorage>::
final() const
{
    return this->_data->final;
}

Void Orbit<HybridStorage>::draw(CanvasInterface& c, const Set<DiscreteLocation>& l, const Variables2d& v) const {
    this->reach().draw(c,l,v);
}

} // namespace Ariadne
