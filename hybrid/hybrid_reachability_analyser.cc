/***************************************************************************
 *            hybrid_reachability_analyser.cc
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

#include "function/functional.h"
#include "config.h"

#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>


#include "utility/exceptions.h"

#include "numeric/numeric.h"

#include "algebra/vector.h"
#include "algebra/matrix.h"

#include "geometry/box.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"

#include "solvers/integrator.h"
#include "solvers/solver.h"

#include "hybrid/hybrid_automaton_interface.h"

#include "hybrid/hybrid_time.h"
#include "hybrid/hybrid_space.h"
#include "hybrid/hybrid_orbit.h"
#include "hybrid/hybrid_set.h"
#include "hybrid/hybrid_evolver.h"
#include "hybrid/hybrid_reachability_analyser.h"

#include "utility/logging.h"
#include "output/graphics.h"
#include "solvers/linear_programming.h"


namespace Ariadne {

static const double DEFAULT_MAXIMUM_ENCLOSURE_RADIUS=0.25;
static const double DEFAULT_GRID_LENGTH=0.125;

HybridReachabilityAnalyser::
~HybridReachabilityAnalyser()
{
}

inline Real operator-(Real const& r1, Rational const& q2) {
    return r1-Real(q2);
}

inline ValidatedFloat64 operator*(Integer n1,  const ExactFloat64& x2) {
    ARIADNE_ASSERT(n1==n1.get_si());
    ExactFloat64 x1((Int)n1.get_si());
    return x1*x2;
}

template<> inline ApproximateFloat64 numeric_cast<ApproximateFloat64>(Real const& x) {
    return ApproximateFloat64(x);
}

template<> inline ExactFloat64 numeric_cast<ExactFloat64>(Real const& x) {
    return ExactFloat64(ApproximateFloat64(x).raw());
}

inline DiscreteTimeType div_floor(Real const& t, ExactNumber h) {
    return integer_cast<Nat>(numeric_cast<ApproximateFloat64>(t)/h);
}

template<> inline Nat integer_cast<Nat>(ValidatedFloat64 const& x) {
    return integer_cast<Nat>(ApproximateFloat64(x).raw());
}


HybridReachabilityAnalyser::
HybridReachabilityAnalyser(
        const SystemType& system,
        const HybridEvolverInterface& evolver)
    : _system(system.clone())
    , _evolver(evolver.clone())
    , _configuration(new ConfigurationType(*this))
{
}


HybridReachabilityAnalyser*
HybridReachabilityAnalyser::
clone() const
{
    return new HybridReachabilityAnalyser(*this);
}


Pair<HybridGridTreeSet,HybridGridTreeSet>
HybridReachabilityAnalyser::_reach_evolve_resume(const ListSet<HybridEnclosure>& initial_enclosures,
                                                const HybridTime& time,
                                                const Int accuracy,
                                                ListSet<HybridEnclosure>& evolve_enclosures,
                                                Semantics semantics,
                                                const HybridEvolverInterface& evolver) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::_reach_evolve_resume(...)\n");
    const HybridGrid& grid=this->_configuration->grid();
    Pair<HybridGridTreeSet,HybridGridTreeSet> result=make_pair(HybridGridTreeSet(grid),HybridGridTreeSet(grid));
    HybridGridTreeSet& reach_cells=result.first; HybridGridTreeSet& evolve_cells=result.second;

    for(ListSet<HybridEnclosure>::ConstIterator encl_iter=initial_enclosures.begin(); encl_iter!=initial_enclosures.end(); ++encl_iter) {
        ListSet<HybridEnclosure> current_reach_enclosures;
        ListSet<HybridEnclosure> current_evolve_enclosures;
        make_lpair(current_reach_enclosures,current_evolve_enclosures) = evolver.reach_evolve(*encl_iter,time,semantics);

        for(ListSet<HybridEnclosure>::ConstIterator enclosure_iter=current_reach_enclosures.begin(); enclosure_iter!=current_reach_enclosures.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(reach_cells,accuracy);
        }
        for(ListSet<HybridEnclosure>::ConstIterator enclosure_iter=current_evolve_enclosures.begin(); enclosure_iter!=current_evolve_enclosures.end(); ++enclosure_iter) {
            const Enclosure& orig_encl = enclosure_iter->continuous_set();
            evolve_enclosures.adjoin(HybridEnclosure(enclosure_iter->location(),enclosure_iter->space(),
                                                        Enclosure(orig_encl.domain(),orig_encl.space_function(),
                                                                  ValidatedScalarFunction::zero(orig_encl.time_function().argument_size()),
                                                                  orig_encl.constraints(),orig_encl.function_factory())));
            enclosure_iter->adjoin_outer_approximation_to(evolve_cells,accuracy);
        }
    }

    ARIADNE_LOG(3,"  final reach size = "<<reach_cells.size()<<"\n");
    ARIADNE_LOG(3,"  final evolve size = "<<evolve_cells.size()<<"\n");
    ARIADNE_LOG(3,"  final evolve enclosures size = "<<evolve_enclosures.size()<<"\n");
    ARIADNE_LOG(2,"Done.\n");
    return result;
}


// Helper functions for operators on lists of sets.
HybridGridTreeSet
HybridReachabilityAnalyser::_upper_reach(const HybridGridTreeSet& set,
                                         const HybridTime& time,
                                         const Int accuracy,
                                         const HybridEvolverInterface& evolver) const
{
    HybridGrid grid=set.grid();
    HybridGridTreeSet result(grid);
    HybridGridTreeSet cells=set;
    cells.mince(accuracy);
    for(HybridGridTreeSet::ConstIterator cell_iter=cells.begin(); cell_iter!=cells.end(); ++cell_iter) {
        HybridEnclosure initial_enclosure = evolver.enclosure(cell_iter->box());
        ListSet<HybridEnclosure> reach = evolver.reach(initial_enclosure,time,UPPER_SEMANTICS);
        for(ListSet<HybridEnclosure>::ConstIterator enclosure_iter=reach.begin(); enclosure_iter!=reach.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(result,accuracy);
        }
    }
    return result;
}


HybridGridTreeSet
HybridReachabilityAnalyser::_upper_evolve(const HybridGridTreeSet& set,
                                          const HybridTime& time,
                                          const Int accuracy,
                                          const HybridEvolverInterface& evolver) const
{
    HybridGrid grid=set.grid();
    HybridGridTreeSet result(grid); HybridGridTreeSet cells=set; cells.mince(accuracy);
    for(HybridGridTreeSet::ConstIterator cell_iter=cells.begin(); cell_iter!=cells.end(); ++cell_iter) {
        ARIADNE_LOG(5,"Evolving cell = "<<*cell_iter<<"\n");
        HybridEnclosure initial_enclosure = evolver.enclosure(cell_iter->box());
        ListSet<HybridEnclosure> final = evolver.evolve(initial_enclosure,time,UPPER_SEMANTICS);
        for(ListSet<HybridEnclosure>::ConstIterator enclosure_iter=final.begin(); enclosure_iter!=final.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(result,accuracy);
        }
    }
    ARIADNE_LOG(4,"_upper_evolve result size = "<<result.size()<<"\n");
    return result;
}


Void
HybridReachabilityAnalyser::_adjoin_upper_reach_evolve(HybridGridTreeSet& reach_cells,
                                                       HybridGridTreeSet& evolve_cells,
                                                       const HybridGridTreeSet& set,
                                                       const HybridTerminationCriterion& termination,
                                                       const Int accuracy,
                                                       const HybridEvolverInterface& evolver) const
{
    ARIADNE_LOG(4,"HybridReachabilityAnalyser::_adjoin_upper_reach_evolve(...)\n");
    HybridGrid grid=set.grid();
    HybridGridTreeSet cells=set;
    cells.mince(accuracy);

    ARIADNE_LOG(3,"Evolving "<<cells.size()<<" cells\n");
    for(HybridGridCell const& cell : cells) {
        ARIADNE_LOG(5,"Evolving cell = "<<cell<<"\n");
        HybridEnclosure initial_enclosure = evolver.enclosure(cell.box());
        ListSet<HybridEnclosure> reach_enclosures;
        ListSet<HybridEnclosure> final_enclosures;
        make_lpair(reach_enclosures,final_enclosures) = evolver.reach_evolve(initial_enclosure,termination,UPPER_SEMANTICS);
        ARIADNE_LOG(5,"  computed "<<reach_enclosures.size()<<" reach enclosures and "<<final_enclosures.size()<<" final enclosures.\n");
        ARIADNE_LOG(5,"  adjoining reach enclosures to grid... ");
        for(HybridEnclosure const& enclosure : reach_enclosures) {
            enclosure.adjoin_outer_approximation_to(reach_cells,accuracy);
        }
        ARIADNE_LOG(5,"done\n");
        ARIADNE_LOG(5,"  adjoining final enclosures to grid... ");
        for(HybridEnclosure const& enclosure : final_enclosures) {
            enclosure.adjoin_outer_approximation_to(evolve_cells,accuracy);
        }
        ARIADNE_LOG(5,"done.\n");
    }
    ARIADNE_LOG(3,"  final reach size = "<<reach_cells.size()<<"\n");
    ARIADNE_LOG(3,"  final evolve size = "<<evolve_cells.size()<<"\n");
    ARIADNE_LOG(2,"Done.\n");
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_evolve(const OvertSetInterfaceType& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_evolve(...)\n");
    Int grid_depth = this->_configuration->maximum_grid_depth();
    Int grid_height = this->_configuration->maximum_grid_height();
    HybridGrid grid=this->_configuration->grid();
    HybridGridTreeSet initial_cells(grid); HybridGridTreeSet final_cells(grid);

    // Improve accuracy of initial set for lower computations
    initial_cells.adjoin_lower_approximation(initial_set,grid_height,grid_depth+4);
    ARIADNE_LOG(3,"initial_cells.size()="<<initial_cells.size()<<"\n");
    ARIADNE_LOG(3,"computing lower evolution.");
    for(HybridGridTreeSet::ConstIterator cell_iter=initial_cells.begin(); cell_iter!=initial_cells.end(); ++cell_iter) {
        ARIADNE_LOG(3,".");
        HybridGridCell cell=*cell_iter;
        HybridEnclosure initial_enclosure=_evolver->enclosure(cell_iter->box());
        ListSet<HybridEnclosure> final_enclosures=_evolver->evolve(initial_enclosure,time,LOWER_SEMANTICS);
        for(ListSet<HybridEnclosure>::ConstIterator enclosure_iter=final_enclosures.begin(); enclosure_iter!=final_enclosures.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(final_cells,grid_depth);
        }
    }
    ARIADNE_LOG(3,"\n");
    return final_cells;
}



HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_reach(const OvertSetInterfaceType& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach(set,time)\n");
    Int grid_depth = this->_configuration->maximum_grid_depth();
    Int grid_height = this->_configuration->maximum_grid_height();
    const HybridGrid& grid=this->_configuration->grid();
    HybridGridTreeSet initial_cells(grid); HybridGridTreeSet reach_cells(grid);

    ARIADNE_LOG(3,"Adjoining initial set to the grid...\n");
    // Improve accuracy of initial set for lower computations
    initial_cells.adjoin_lower_approximation(initial_set,grid_height,grid_depth+4);
    ARIADNE_LOG(3,"initial_cells.size()="<<initial_cells.size()<<"\n");
    ARIADNE_LOG(3,"Computing lower reach set...");
    for(HybridGridTreeSet::ConstIterator cell_iter=initial_cells.begin(); cell_iter!=initial_cells.end(); ++cell_iter) {
        ARIADNE_LOG(3,".");
        HybridGridCell cell=*cell_iter;
        HybridEnclosure initial_enclosure=_evolver->enclosure(cell_iter->box());
        ListSet<HybridEnclosure> reach_enclosures=_evolver->reach(initial_enclosure,time,LOWER_SEMANTICS);
        for(ListSet<HybridEnclosure>::ConstIterator enclosure_iter=reach_enclosures.begin(); enclosure_iter!=reach_enclosures.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(reach_cells,grid_depth);
        }
    }
    ARIADNE_LOG(3,"\n");
    return reach_cells;
}


Pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
lower_reach_evolve(const OvertSetInterfaceType& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach_evolve(...)\n");
    Int grid_depth = this->_configuration->maximum_grid_depth();
    Int grid_height = this->_configuration->maximum_grid_height();

    const HybridGrid& grid=this->_configuration->grid();

    HybridGridTreeSet initial_cells(grid);

    HybridGridTreeSet reach=(grid); HybridGridTreeSet evolve_cells(grid);

    // Improve accuracy of initial set for lower computations
    initial_cells.adjoin_lower_approximation(initial_set,grid_height,grid_depth+4);
    ARIADNE_LOG(3,"initial_cells.size()="<<initial_cells.size()<<"\n");
    ARIADNE_LOG(3,"computing lower evolution.");
    for(HybridGridTreeSet::ConstIterator cell_iter=initial_cells.begin(); cell_iter!=initial_cells.end(); ++cell_iter) {
        ARIADNE_LOG(3,".");
        HybridEnclosure initial_enclosure=_evolver->enclosure(cell_iter->box());
        ListSet<HybridEnclosure> reach_enclosures;
        ListSet<HybridEnclosure> final_enclosures;
        make_lpair(reach_enclosures,final_enclosures) = _evolver->reach_evolve(initial_enclosure,time,LOWER_SEMANTICS);
        for(ListSet<HybridEnclosure>::ConstIterator enclosure_iter=reach_enclosures.begin(); enclosure_iter!=reach_enclosures.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(reach,grid_depth);
        }
        for(ListSet<HybridEnclosure>::ConstIterator enclosure_iter=final_enclosures.begin(); enclosure_iter!=final_enclosures.end(); ++enclosure_iter) {
            enclosure_iter->adjoin_outer_approximation_to(evolve_cells,grid_depth);
        }
    }
    ARIADNE_LOG(3,"\n");
    return make_pair(reach,evolve_cells);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
lower_reach(const OvertSetInterfaceType& initial_set) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::lower_reach(set)\n");
    ExactNumber transient_time = this->_configuration->transient_time();
    Int transient_steps = this->_configuration->transient_steps();
    HybridTime hybrid_transient_time(transient_time, transient_steps);
    ExactNumber lock_to_grid_time=this->_configuration->lock_to_grid_time();
    Int lock_to_grid_steps=this->_configuration->lock_to_grid_steps();
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,lock_to_grid_steps);
    Int maximum_grid_depth = this->_configuration->maximum_grid_depth();
    Int maximum_grid_height = this->_configuration->maximum_grid_height();
    ARIADNE_LOG(3,"transient_time=("<<transient_time<<","<<transient_steps<<")\n");
    ARIADNE_LOG(3,"lock_to_grid_time=("<<lock_to_grid_time<<","<<lock_to_grid_steps<<")\n");
    ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");

    const HybridGrid& grid=this->_configuration->grid();

    Bool has_bounding_domain = static_cast<Bool>(this->_configuration->bounding_domain_ptr());

    HybridGridTreeSet bounding(grid);
    if (has_bounding_domain)
        bounding.adjoin_outer_approximation(*this->_configuration->bounding_domain_ptr(),maximum_grid_depth);

    ARIADNE_LOG(5,"maximum_grid_height="<<maximum_grid_height<<"\n");
    ARIADNE_LOG(5,"bounding_size="<<bounding.size()<<"\n");

    HybridGridTreeSet initial_cells(grid), evolve_cells(grid);
    initial_cells.adjoin_lower_approximation(initial_set,maximum_grid_height,maximum_grid_depth+4);
    if (has_bounding_domain) initial_cells.restrict(bounding);
    ARIADNE_LOG(5,"initial_size="<<initial_cells.size()<<"\n");

    ListSet<HybridEnclosure> starting_enclosures;

    HybridGridTreeSet reach_cells(grid);

    if(transient_time > 0.0_exact || transient_steps > 0) {
        ARIADNE_LOG(3,"Computing transient evolution...\n");

        ListSet<HybridEnclosure> initial_enclosures;
        for(HybridGridTreeSet::ConstIterator cell_iter=initial_cells.begin(); cell_iter!=initial_cells.end(); ++cell_iter)
            initial_enclosures.adjoin(_evolver->enclosure(cell_iter->box()));

        make_lpair(reach_cells,evolve_cells) = _reach_evolve_resume(initial_enclosures,hybrid_transient_time,
                maximum_grid_depth,starting_enclosures,LOWER_SEMANTICS,*_evolver);

        evolve_cells.restrict_to_height(maximum_grid_height);
        ARIADNE_LOG(3,"Completed restriction to height.");
        if (has_bounding_domain) reach_cells.restrict(bounding);

        ARIADNE_LOG(5,"reach size="<<reach_cells.size()<<"\n");
        ARIADNE_LOG(5,"evolve size="<<evolve_cells.size()<<"\n");
        ARIADNE_LOG(3,"  found "<<reach_cells.size()<<" cells.\n");
    }
    ARIADNE_LOG(3,"Computing recurrent evolution...\n");

    while(!evolve_cells.is_empty()) {

        HybridGridTreeSet new_reach_cells(grid);
        ListSet<HybridEnclosure> new_evolve_enclosures;
        make_lpair(new_reach_cells,evolve_cells) = _reach_evolve_resume(starting_enclosures,hybrid_lock_to_grid_time,
                maximum_grid_depth,new_evolve_enclosures,LOWER_SEMANTICS,*_evolver);

        ARIADNE_LOG(3,"Removing...\n");
        evolve_cells.remove(reach_cells);
        ARIADNE_LOG(3,"Restricting to height...\n");
        evolve_cells.restrict_to_height(maximum_grid_height);
        ARIADNE_LOG(3,"Restricting to bounding...\n");
        if (has_bounding_domain) evolve_cells.restrict(bounding);
        ARIADNE_LOG(3,"Adjoining...\n");
        reach_cells.adjoin(new_reach_cells);
        ARIADNE_LOG(3,"  found "<<new_reach_cells.size()<<" cells, of which "<<evolve_cells.size()<<" are new.\n");

        starting_enclosures = new_evolve_enclosures;
    }
    reach_cells.recombine();
    reach_cells.restrict_to_height(maximum_grid_height);
    if (has_bounding_domain) reach_cells.restrict(bounding);

    return reach_cells;
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_evolve(const CompactSetInterfaceType& initial_set,
             const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_evolve(...)\n");
    const HybridGrid& grid=this->_configuration->grid();
    HybridGridTreeSet evolve_cells(grid);
    Int grid_depth = this->_configuration->maximum_grid_depth();
    evolve_cells.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial_evolve.size()="<<evolve_cells.size()<<"\n");
    Real real_time=time.continuous_time();
    DiscreteTimeType discrete_steps=time.discrete_time();
    ExactNumber lock_to_grid_time=this->_configuration->lock_to_grid_time();
    DiscreteTimeType time_steps=integer_cast<Nat>(numeric_cast<ExactFloat64>(real_time)/lock_to_grid_time);
    Real remainder_time=real_time-Real(time_steps*Rational(lock_to_grid_time));
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");

    for(Nat i=0; i!=time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        evolve_cells=this->_upper_evolve(evolve_cells,hybrid_lock_to_grid_time,grid_depth,*_evolver);
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve_cells.is_empty() && possibly(remainder_time > 0)) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        evolve_cells=this->_upper_evolve(evolve_cells,hybrid_remainder_time,grid_depth,*_evolver);
    }
    evolve_cells.recombine();
    ARIADNE_LOG(4,"final_evolve.size()="<<evolve_cells.size()<<"\n");
    return evolve_cells;
}



HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
upper_reach(const CompactSetInterfaceType& initial_set,
            const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach(set,time)\n");
    ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n");
    const HybridGrid& grid=this->_configuration->grid();
    HybridGridTreeSet evolve_cells(grid);
    Int grid_depth = this->_configuration->maximum_grid_depth();
    ARIADNE_LOG(4,"grid_depth="<<grid_depth<<"\n");
    evolve_cells.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial size = "<<evolve_cells.size()<<"\n");
    HybridGridTreeSet reach_cells(evolve_cells);
    ARIADNE_LOG(4,"reach size ="<<reach_cells.size()<<"\n");
    Real real_time=time.continuous_time();
    DiscreteTimeType discrete_steps=time.discrete_time();
    ExactNumber lock_to_grid_time = this->_configuration->lock_to_grid_time();
    DiscreteTimeType time_steps=div_floor(real_time,lock_to_grid_time);
    Real remainder_time=real_time-time_steps*Rational(lock_to_grid_time);
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);

    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");
    ARIADNE_LOG(3,"discrete_steps="<<discrete_steps<<"\n");
    HybridGridTreeSet found_cells(grid), accumulated_evolve_cells(grid);
    for(Nat i=0; i!=time_steps; ++i) {
        accumulated_evolve_cells.adjoin(found_cells);
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        this->_adjoin_upper_reach_evolve(reach_cells,evolve_cells,evolve_cells,hybrid_lock_to_grid_time,grid_depth,*_evolver);
        ARIADNE_LOG(5,"found.size()="<<found_cells.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve_cells.size()<<"\n");
        evolve_cells.remove(accumulated_evolve_cells);
        reach_cells.adjoin(found_cells);
        accumulated_evolve_cells.adjoin(evolve_cells);
        ARIADNE_LOG(3,"  found "<<found_cells.size()<<" cells, with "<<evolve_cells.size()<<" new intermediate.\n");
        if(evolve_cells.is_empty()) break;
        ARIADNE_LOG(6,"evolve_cells="<<evolve_cells<<"\n");
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve_cells.is_empty() && possibly(remainder_time > 0)) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        this->_adjoin_upper_reach_evolve(found_cells,accumulated_evolve_cells,evolve_cells,hybrid_remainder_time,grid_depth,*_evolver);
        reach_cells.adjoin(found_cells);
    }
    // This last step is necessary to add the final set to the result.
    reach_cells.adjoin(evolve_cells);
    reach_cells.recombine();
    ARIADNE_LOG(4,"final_reach size = "<<reach_cells.size()<<"\n");
    return reach_cells;
}



Pair<HybridReachabilityAnalyser::SetApproximationType,HybridReachabilityAnalyser::SetApproximationType>
HybridReachabilityAnalyser::
upper_reach_evolve(const CompactSetInterfaceType& initial_set,
                   const TimeType& time) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::upper_reach_evolve(...)\n");
    ARIADNE_LOG(4,"initial_set="<<initial_set<<"\n");
    const HybridGrid& grid=this->_configuration->grid();
    HybridGridTreeSet evolve_cells(grid);
    Int grid_depth = this->_configuration->maximum_grid_depth();
    ARIADNE_LOG(4,"grid_depth="<<grid_depth<<"\n");
    evolve_cells.adjoin_outer_approximation(initial_set,grid_depth);
    ARIADNE_LOG(4,"initial_evolve"<<evolve_cells<<"\n");
    HybridGridTreeSet reach_cells(evolve_cells);
    ARIADNE_LOG(4,"reach="<<reach_cells<<"\n");
    Real real_time=time.continuous_time();
    DiscreteTimeType discrete_steps=time.discrete_time();
    ExactNumber lock_to_grid_time = this->_configuration->lock_to_grid_time();
    DiscreteTimeType time_steps=div_floor(real_time,lock_to_grid_time);
    Real remainder_time=real_time-time_steps*Rational(lock_to_grid_time);
    HybridTime hybrid_lock_to_grid_time(lock_to_grid_time,discrete_steps);
    HybridTime hybrid_remainder_time(remainder_time,discrete_steps);
    ARIADNE_LOG(3,"real_time="<<real_time<<"\n");
    ARIADNE_LOG(3,"time_steps="<<time_steps<<"  lock_to_grid_time="<<lock_to_grid_time<<"\n");

    HybridGridTreeSet found_cells(grid);
    for(Nat i=0; i!=time_steps; ++i) {
        ARIADNE_LOG(3,"computing "<<i+1<<"-th reachability step...\n");
        this->_adjoin_upper_reach_evolve(found_cells,evolve_cells,evolve_cells,hybrid_lock_to_grid_time,grid_depth,*_evolver);
        ARIADNE_LOG(5,"found.size()="<<found_cells.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve_cells.size()<<"\n");
        reach_cells.adjoin(found_cells);
        ARIADNE_LOG(3,"  found "<<found_cells.size()<<" cells.\n");
    }
    ARIADNE_LOG(3,"remainder_time="<<remainder_time<<"\n");
    if(!evolve_cells.is_empty() && possibly(remainder_time > 0)) {
        ARIADNE_LOG(3,"computing evolution for remainder time...\n");
        this->_adjoin_upper_reach_evolve(found_cells,evolve_cells,evolve_cells,hybrid_remainder_time,grid_depth,*_evolver);
        reach_cells.adjoin(found_cells);
    }
    reach_cells.recombine();
    ARIADNE_LOG(4,"reach="<<reach_cells<<"\n");
    evolve_cells.recombine();
    ARIADNE_LOG(4,"evolve="<<evolve_cells<<"\n");
    return std::make_pair(reach_cells,evolve_cells);
}


HybridReachabilityAnalyser::SetApproximationType
HybridReachabilityAnalyser::
outer_chain_reach(
        const CompactSetInterfaceType& initial_set) const
{
    ARIADNE_LOG(2,"HybridReachabilityAnalyser::outer_chain_reach(...)\n");
    Set<DiscreteEvent> lock_to_grid_events=this->_configuration->lock_to_grid_events();
    ExactNumber transient_time = this->_configuration->transient_time();
    DiscreteTimeType transient_steps = this->_configuration->transient_steps();
    HybridTerminationCriterion transient_termination(EffectiveNumber(transient_time), transient_steps, lock_to_grid_events);
    ExactNumber lock_to_grid_time=this->_configuration->lock_to_grid_time();
    DiscreteTimeType lock_to_grid_steps=this->_configuration->lock_to_grid_steps();
    HybridTerminationCriterion recurrent_termination(EffectiveNumber(lock_to_grid_time),lock_to_grid_steps,lock_to_grid_events);
    Int maximum_grid_depth = this->_configuration->maximum_grid_depth();
    Int maximum_grid_height = this->_configuration->maximum_grid_height();
    ARIADNE_LOG(3,"transient_time=("<<transient_time<<","<<transient_steps<<")\n");
    ARIADNE_LOG(3,"lock_to_grid_time=("<<lock_to_grid_time<<","<<lock_to_grid_steps<<")\n");
    ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");

    const HybridGrid& grid=this->_configuration->grid();

    Bool has_bounding_domain = static_cast<Bool>(this->_configuration->bounding_domain_ptr());

    HybridGridTreeSet bounding(grid);
    if (has_bounding_domain)
        bounding.adjoin_outer_approximation(*this->_configuration->bounding_domain_ptr(),maximum_grid_depth);

    ARIADNE_LOG(5,"maximum_grid_height="<<maximum_grid_height<<"\n");
    ARIADNE_LOG(5,"bounding_size="<<bounding.size()<<"\n");

    HybridGridTreeSet initial_cells(grid);
    initial_cells.adjoin_outer_approximation(initial_set,maximum_grid_depth);
    _checked_restriction(initial_cells,bounding);
    ARIADNE_LOG(5,"initial_size="<<initial_cells.size()<<"\n");


    HybridGridTreeSet reach_cells(grid);
    HybridGridTreeSet evolve_cells(grid);
    if(transient_time > 0.0_exact || transient_steps > 0) {
        ARIADNE_LOG(3,"Computing transient evolution...\n");
        this->_adjoin_upper_reach_evolve(reach_cells,evolve_cells,initial_cells,transient_termination,maximum_grid_depth,*_evolver);
        _checked_restriction(evolve_cells,bounding);
        evolve_cells.recombine();
        evolve_cells.mince(maximum_grid_depth);
        ARIADNE_LOG(5,"transient_reach_size="<<reach_cells.size()<<"\n");
        ARIADNE_LOG(5,"transient evolve_size="<<evolve_cells.size()<<"\n");
        ARIADNE_LOG(3,"  found "<<reach_cells.size()<<" cells.\n");
    } else {
        evolve_cells=initial_cells;
    }

    ARIADNE_LOG(3,"Computing recurrent evolution...\n");
    HybridGridTreeSet starting_cells = evolve_cells;
    HybridGridTreeSet current_evolve_cells(grid);


    while(!starting_cells.is_empty()) {
        current_evolve_cells = evolve_cells;
        this->_adjoin_upper_reach_evolve(reach_cells,evolve_cells,starting_cells,
                                         recurrent_termination,maximum_grid_depth,*_evolver);
        ARIADNE_LOG(5,"reach.size()="<<reach_cells.size()<<"\n");
        ARIADNE_LOG(5,"evolve.size()="<<evolve_cells.size()<<"\n");
        _checked_restriction(evolve_cells,bounding);
        evolve_cells.recombine();
        evolve_cells.mince(maximum_grid_depth);
        starting_cells = evolve_cells;
        starting_cells.remove(current_evolve_cells);
        starting_cells.mince(maximum_grid_depth);
        ARIADNE_LOG(3,"  evolved to "<<evolve_cells.size()<<" cells, of which "<<starting_cells.size()<<" are new.\n");
        // ARIADNE_LOG(5,"bounded_new_size="<<found_cells.size()<<"\n");
    }
    _checked_restriction(reach_cells,bounding);
    return reach_cells;
}


Void HybridReachabilityAnalyser::_checked_restriction(HybridGridTreeSet& set, const HybridGridTreeSet& bounding) const
{
    ChainOverspillPolicy policy = this->_configuration->outer_overspill_policy();
    HybridGridTreeSet set_copy(set.grid());

    if (policy != OVERSPILL_IGNORE) {
        set_copy = set;
    }
    set.restrict_to_height(this->_configuration->maximum_grid_height());
    if (!bounding.is_empty()) set.restrict(bounding);

    if (policy != OVERSPILL_IGNORE) {
        set_copy.remove(set);
        if (!set_copy.is_empty()) {
            if (policy == OVERSPILL_WARNING) {
                ARIADNE_WARN("The computed chain reach has been restricted, an outer approximation is .");
            }
            if (policy == OVERSPILL_ERROR) {
                ARIADNE_THROW(OuterChainOverspill,"outer_chain_reach","The computed chain reach has been restricted.");
            }
        }
    }
}



HybridReachabilityAnalyserConfiguration::HybridReachabilityAnalyserConfiguration(HybridReachabilityAnalyser& analyser)
    : _analyser(analyser)
{
    set_transient_time(0.0);
    set_transient_steps(0);
    set_lock_to_grid_time(1.0);
    set_lock_to_grid_steps(1);
    set_maximum_grid_depth(3);
    set_maximum_grid_height(16);
    set_grid(std::shared_ptr<HybridGrid>(new HybridGrid(_analyser.system().state_space(),SimpleHybridScaling())));
    set_outer_overspill_policy(OVERSPILL_ERROR);
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
HybridReachabilityAnalyserConfiguration::set_bounding_domain_ptr(const std::shared_ptr<HybridBoxes> value)
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
        case OVERSPILL_IGNORE: os<<"ignore"; break;
        case OVERSPILL_WARNING: os<<"warning"; break;
        case OVERSPILL_ERROR: os<<"error"; break;
        default: os << "unknown"; break;
    } return os;
}

} // namespace Ariadne
