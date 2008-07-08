/***************************************************************************
 *            map_evolver.code.h
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
 
#include "map_evolver.h"

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

#include "evaluation/standard_applicator.h"
#include "evaluation/standard_subdivider.h"
#include "evaluation/cascade_reducer.h"

#include "system/map.h"


#include "output/logging.h"

namespace Ariadne {
  

template<class ES>
Evolver< Map<typename ES::real_type>, ES>::
Evolver(const EvolutionParameters<R>& parameters,
        const ApplicatorInterface<ES>& applicator, 
        const SubdividerInterface<ES>& subdivider, 
        const ReducerInterface<ES>& reducer)
  : _parameters(parameters.clone()),
    _applicator(applicator.clone()),
    _subdivider(subdivider.clone()),
    _reducer(reducer.clone())
{ }

template<class ES>
Evolver< Map<typename ES::real_type>, ES>::
Evolver(const EvolutionParameters<R>& parameters,
        const ApplicatorInterface<ES>& applicator) 
{ 
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  Evolver(parameters,applicator,subdivider,reducer);
}

template<class ES>
Evolver< Map<typename ES::real_type>, ES>::
Evolver(const EvolutionParameters<R>& parameters) { 
  StandardApplicator< ES > applicator;
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  Evolver(parameters,applicator,subdivider,reducer);
}


template<class ES>
Evolver< Map<typename ES::real_type>, ES>::
Evolver() { 
  EvolutionParameters<R> parameters;
  StandardApplicator< ES > applicator;
  StandardSubdivider< ES > subdivider;
  CascadeReducer< ES > reducer(3);
  Evolver(parameters,applicator,subdivider,reducer);
}





template<class ES>
void
Evolver< Map<typename ES::real_type>, ES>::
_evolution(ESL& final,
           ESL& intermediate, 
           const Sys& system,
           const ES& initial,
           const T& time,
           Semantics semantics,
           bool reach) const
{
  uint verbosity=this->verbosity();

  TESL working;
  working.adjoin(TES(T(0),initial)); 
  ARIADNE_LOG(5,"  working.size()="<<working.size()<<"\n");
  if(reach) {
    intermediate.adjoin(initial);
  }
  while(working.size()!=0) {
    TES ts=working.pop();
    ARIADNE_LOG(5,"  ts="<<ts<<", r="<<this->radius(ts)<<"\n");
    if(this->radius(ts) > this->maximum_enclosure_radius()) {
      if(semantics==upper_semantics) {
        ARIADNE_LOG(5,"    subdivide...\n");
        //ARIADNE_LOG(7,"      "<<this->subdivide(ts)<<"\n");
        this->adjoin_subdivision(working,ts);
      } else if(semantics==lower_semantics) {
        ARIADNE_LOG(5,"    blocking...\n");
      }
    } else if(ts.time()==time) {
      ARIADNE_LOG(5,"    end...\n");
      final.adjoin(ts.set());
    } else {
      ARIADNE_LOG(5,"  apply... ");
      ts=this->apply(system,ts);
      ARIADNE_LOG(7,ts);
      ts=this->reduce(ts);
      ARIADNE_LOG(7,ts);
      ARIADNE_LOG(5,"\n");
      if(reach) {
        intermediate.adjoin(ts.set());
      }
      working.adjoin(ts);
    }
  }
}



}  // namespace Ariadne

