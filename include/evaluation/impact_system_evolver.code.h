/***************************************************************************
 *            impact_system_evolver.code.h
 *
 *  Copyright  2006-8  Alberto Casagrande, Pieter Collins
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
 
#include "impact_system_evolver.h"

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "combinatoric/lattice_set.h"

#include "geometry/box.h"
#include "geometry/box_list_set.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/grid_approximation.h"

#include "system/impact_system.h"

#include "evaluation/evolution_profiler.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/standard_integrator.h"
#include "evaluation/standard_approximator.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/standard_reducer.h"

#include "output/logging.h"

namespace Ariadne {
  
template<class ES>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
Evolver()
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class ES>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
Evolver(const EvolutionParameters<R>&)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class ES>
Evolver<ImpactSystem<typename ES::real_type>,ES>::
Evolver(const EvolutionParameters<R>& parameters,
        const ApplicatorInterface<ES>& applicator, 
        const IntegratorInterface<ES>& integrator, 
        const SatisfierInterface<ES>& satisfier, 
        const SubdividerInterface<ES>& subdivider, 
        const ReducerInterface<ES>& reducer)
  : _parameters(parameters.clone()),
    _integrator(integrator.clone()),
    _subdivider(subdivider.clone()),
    _reducer(reducer.clone()),
    _profiler(new EvolutionProfiler)
{ }

template<class ES>
void
Evolver<ImpactSystem<typename ES::real_type>,ES>::
_evolution(ESL& final,
           ESL& intermediate, 
           const Sys& sys,
           const ES& initial,
           const T& time,
           Semantics semantics,
           bool reach) const
{
  typedef Box<R> Bx;
  VectorField<R> vf(sys.vector_field());

  uint verbosity=this->verbosity();
  ARIADNE_LOG(7,"ImpactSystemEvolver::_step(...)\n");

  TESL working;
  working.adjoin(TES(T(0),initial)); 

  while(working.size()!=0) {
    TES ts=working.pop();
    const ES& ws=ts.set();
    ARIADNE_LOG(7,"  ts="<<ts<<"\n");
    ARIADNE_ASSERT(ts.time()<=time);
    if(ts.time()==time) {
      final.adjoin(ts.set());
    } else if(radius(ts.set()) > maximum_basic_set_radius()) {
      if(semantics==upper_semantics) {
        this->adjoin_subdivision(working,ts);
        ++this->_profiler->subdivisions;
      } 
    } else {
      ++this->_profiler->time_steps;
      Bx bb; T h; ES rs; 
      make_lpair(h,bb)=this->flow_bounds(vf,this->bounding_box(ts.set()));
      this->_profiler->total_stepping_time+=h;
      this->_profiler->minimum_time_step=std::min(h,this->_profiler->minimum_time_step);
      h=std::min(h,T(time-ts.time()));
      {
        rs=this->reachability_step(vf,ws,h,bb);
        intermediate.adjoin(rs);
      }
      // Need to do reachability step first to avoid clobbering bs reference
      ts=this->integration_step(vf,ts,h,bb);
      ts=this->reduce(ts);
      working.adjoin(ts);
    }
  }
}



}  // namespace Ariadne
