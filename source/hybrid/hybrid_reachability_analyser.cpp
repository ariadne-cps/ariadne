/***************************************************************************
 *            hybrid_reachability_analyser.cpp
 *
 *  Copyright  2006-11  Alberto Casagrande, Pieter Collins, Davide Bresolin
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "function/functional.hpp"
#include "config.h"

#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>


#include "utility/exceptions.hpp"

#include "numeric/numeric.hpp"

#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"

#include "geometry/box.hpp"
#include "geometry/list_set.hpp"
#include "geometry/grid_set.hpp"

#include "solvers/integrator.hpp"
#include "solvers/solver.hpp"

#include "hybrid/hybrid_automaton_interface.hpp"

#include "hybrid/hybrid_time.hpp"
#include "hybrid/hybrid_space.hpp"
#include "hybrid/hybrid_orbit.hpp"
#include "hybrid/hybrid_set.hpp"
#include "hybrid/hybrid_evolver.hpp"
#include "hybrid/hybrid_reachability_analyser.hpp"

#include "utility/logging.hpp"
#include "output/graphics.hpp"
#include "solvers/linear_programming.hpp"

#include "hybrid/hybrid_graphics.hpp"

namespace Ariadne {

static const double DEFAULT_MAXIMUM_ENCLOSURE_RADIUS=0.25;
static const double DEFAULT_GRID_LENGTH=0.125;

template class ReachabilityAnalyser<HybridAutomatonInterface>;
template class ReachabilityAnalyserConfiguration<HybridAutomatonInterface>;

inline Real operator-(Real const& r1, Rational const& q2) {
    return r1-Real(q2);
}

inline Float64Bounds operator*(Integer n1,  const Float64Value& x2) {
    ARIADNE_ASSERT(n1==n1.get_si());
    Float64Value x1((Int)n1.get_si());
    return x1*x2;
}

template<> inline Float64Approximation numeric_cast<Float64Approximation>(Real const& x) {
    return Float64Approximation(x,Precision64());
}

template<> inline Float64Value numeric_cast<Float64Value>(Real const& x) {
    return Float64Value(numeric_cast<Float64Approximation>(x).raw());
}

inline DiscreteTimeType div_floor(Real const& t, ExactNumericType h) {
    return integer_cast<Nat>(numeric_cast<Float64Approximation>(t)/h);
}

template<> Nat integer_cast<Nat,Real>(Real const& r);

template<> inline Nat integer_cast<Nat>(Float64Bounds const& x) {
    return integer_cast<Nat>(Float64Approximation(x).raw());
}

inline Nat compute_time_steps(HybridTime time, HybridTime lock_to_grid_time) {
    return integer_cast<Nat>(time.continuous_time()/lock_to_grid_time.continuous_time());
}

inline String itoa(Nat m) {
    std::stringstream ss; if(m<10) { ss<<"0"; } ss<<m; return ss.str();
}

inline const HybridExactBoxes& cast_exact(const HybridExactBoxes& boxes) {
    return boxes;
}

} // namespace Ariadne

#include "dynamics/reachability_analyser.tpl.hpp"

namespace Ariadne {

HybridReachabilityAnalyser::
HybridReachabilityAnalyser(
        const SystemType& system,
        const HybridEvolverInterface& evolver)
    : ReachabilityAnalyser<HybridAutomatonInterface>(system,evolver)
{
}


HybridReachabilityAnalyser*
HybridReachabilityAnalyser::
clone() const
{
    return new HybridReachabilityAnalyser(*this);
}




Void
HybridReachabilityAnalyser::_adjoin_upper_reach_evolve(HybridGridTreeSet& reach_cells,
                                                       HybridGridTreeSet& evolve_cells,
                                                       const HybridGridTreeSet& set,
                                                       const HybridTerminationCriterion& termination,
                                                       const Int accuracy,
                                                       const HybridEvolverInterface& evolver) const
{
    ARIADNE_LOG(6,"HybridReachabilityAnalyser::_adjoin_upper_reach_evolve(...)\n");
    HybridGrid grid=set.grid();
    HybridGridTreeSet cells=set;
    cells.mince(accuracy);

    ARIADNE_LOG(6,"Evolving "<<cells.size()<<" cells\n");
    for(HybridGridCell const& cell : cells) {
        ARIADNE_LOG(7,"Evolving cell = "<<cell<<"\n");
        HybridEnclosure initial_enclosure = evolver.enclosure(cell.box());
        ListSet<HybridEnclosure> reach_enclosures;
        ListSet<HybridEnclosure> final_enclosures;
        make_lpair(reach_enclosures,final_enclosures) = evolver.reach_evolve(initial_enclosure,termination,UPPER_SEMANTICS);
        ARIADNE_LOG(7,"  computed "<<reach_enclosures.size()<<" reach enclosures and "<<final_enclosures.size()<<" final enclosures.\n");
        ARIADNE_LOG(7,"  adjoining reach enclosures to grid... ");
        for(HybridEnclosure const& enclosure : reach_enclosures) {
            enclosure.adjoin_outer_approximation_to(reach_cells,accuracy);
        }
        ARIADNE_LOG(7,"done\n");
        ARIADNE_LOG(7,"  adjoining final enclosures to grid... ");
        for(HybridEnclosure const& enclosure : final_enclosures) {
            enclosure.adjoin_outer_approximation_to(evolve_cells,accuracy);
        }
        ARIADNE_LOG(7,"done.\n");
    }
    ARIADNE_LOG(6,"  final reach size = "<<reach_cells.size()<<"\n");
    ARIADNE_LOG(6,"  final evolve size = "<<evolve_cells.size()<<"\n");
    ARIADNE_LOG(6,"Done.\n");
}




HybridReachabilityAnalyserConfiguration::ReachabilityAnalyserConfiguration(ReachabilityAnalyser<HybridAutomatonInterface>& analyser)
    : _analyser(analyser)
{
    set_transient_time(0.0);
    set_transient_steps(0);
    set_lock_to_grid_time(1.0);
    set_lock_to_grid_steps(1);
    set_maximum_grid_depth(3);
    set_maximum_grid_height(16);
    set_grid(std::shared_ptr<HybridGrid>(new HybridGrid(_analyser.system().state_space(),SimpleHybridScaling())));
    set_outer_overspill_policy(ChainOverspillPolicy::OVERSPILL_ERROR);
}


OutputStream&
HybridReachabilityAnalyserConfiguration::write(OutputStream& os) const
{
    os << "HybridReachabilityAnalyserSettings"
       << "(\n  transient_time=" << transient_time()
       << ",\n  transient_steps=" << transient_steps()
       << ",\n  lock_to_grid_steps=" << lock_to_grid_steps()
       << ",\n  lock_to_grid_time=" << lock_to_grid_time()
       << ",\n  maximum_grid_depth=" << maximum_grid_depth()
       << ",\n  maximum_grid_height=" << maximum_grid_height()
       << ",\n  bounding_domain=" << (bounding_domain_ptr() ? "available" : "none")
       << ",\n  grid=" << grid()
       << ",\n  outer_overspill_policy=" << outer_overspill_policy()
       << "\n)\n";
    return os;
}


Void
HybridReachabilityAnalyserConfiguration::set_bounding_domain_ptr(const std::shared_ptr<HybridExactBoxes> value)
{
 //   ARIADNE_ASSERT_MSG(possibly(value_ptr->space() == _analyser.system().state_space()),
 //           "The bounding domain to set has a different hybrid space than the system.");

    _bounding_domain_ptr = value;
}


Void
HybridReachabilityAnalyserConfiguration::set_grid(const std::shared_ptr<HybridGrid> value_ptr)
{
    ARIADNE_ASSERT_MSG(possibly(value_ptr->space() == _analyser.system().state_space()),
            "The grid to set has a different hybrid space than the system.");

    _grid_ptr = value_ptr;
}

OutputStream& operator<<(OutputStream& os, const ChainOverspillPolicy& policy)
{
    switch(policy) {
        case ChainOverspillPolicy::OVERSPILL_IGNORE: os<<"ignore"; break;
        case ChainOverspillPolicy::OVERSPILL_WARNING: os<<"warning"; break;
        case ChainOverspillPolicy::OVERSPILL_ERROR: os<<"error"; break;
        default: abort();
    } return os;
}

} // namespace Ariadne
