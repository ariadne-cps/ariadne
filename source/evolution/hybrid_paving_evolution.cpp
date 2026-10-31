// Operations coupling hybrid set geometry to dynamic systems.
#include "function/functional.hpp"
#include "config.hpp"
#include "hybrid/hybrid_paving.hpp"
#include "system/hybrid_automaton_interface.hpp"
#include "hybrid/hybrid_graphics.hpp"
#include "drawing/grid_paving_draw.hpp"
#include "labelled/figure.hpp"
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

} // namespace Ariadne
