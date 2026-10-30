// Operations coupling hybrid set geometry to dynamic systems.
#include "function/functional.hpp"
#include "config.hpp"
#include "hybrid/hybrid_paving.hpp"
#include "hybrid/hybrid_storage.hpp"
#include "hybrid/hybrid_orbit.hpp"
#include "system/hybrid_automaton_interface.hpp"
#include "hybrid/hybrid_graphics.hpp"
#include "drawing/grid_paving_draw.hpp"
#include "labelling/figure.hpp"
#include "numeric/rounding.hpp"
#include "symbolic/assignment.hpp"

namespace Ariadne {

HybridGridTreePaving extend_auxiliary(const HybridGridTreePaving& hybrid_paving,
                                      const HybridAutomatonInterface& hybrid_automaton) {
    HybridGrid hybrid_state_grid = hybrid_paving.grid();
    HybridSpace hybrid_state_space = hybrid_state_grid.space();
    HybridScalings scalings = hybrid_state_grid.scalings();

    HybridSpace hybrid_state_auxiliary_space = hybrid_automaton.state_auxiliary_space();
    HybridGrid hybrid_state_auxiliary_grid(hybrid_state_auxiliary_space,scalings);

    HybridGridTreePaving result(hybrid_state_auxiliary_grid);
    for (auto iter = hybrid_paving.locations_begin(); iter!=hybrid_paving.locations_end(); ++iter) {
        DiscreteLocation const& location=iter->first;
        GridTreePaving const& state_paving=iter->second;
        RealSpace auxiliary_space=hybrid_automaton.continuous_auxiliary_space(location);
        EffectiveVectorMultivariateFunction auxiliary_function=hybrid_automaton.auxiliary_function(location);
        Grid auxiliary_grid=scalings.grid(location,auxiliary_space);
        GridTreePaving state_auxiliary_paving=outer_skew_product(state_paving,auxiliary_grid,auxiliary_function);
        result[location]=state_auxiliary_paving;
    }
    return result;
}



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


} // namespace Ariadne
