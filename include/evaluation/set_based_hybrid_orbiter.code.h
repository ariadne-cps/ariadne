/***************************************************************************
 *            set_based_hybrid_orbiter.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include "set_based_hybrid_orbiter.h"

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "base/tuple.h"

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "combinatoric/lattice_set.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/hybrid_set.h"

#include "system/hybrid_automaton.h"

#include "evaluation/bounder_interface.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/approximator_interface.h"

#include "evaluation/lohner_integrator.h"
#include "evaluation/standard_bounder.h"

#include "output/logging.h"

namespace Ariadne {
  
namespace Evaluation { static int& verbosity = applicator_verbosity; }



template<class BS>
Evaluation::SetBasedHybridOrbiter<BS>::SetBasedHybridOrbiter(const EvolutionParameters<R>& parameters, 
                                                             const IntegratorInterface<BS>& integrator, 
                                                             const ApplicatorInterface<BS>& applicator, 
                                                             const ApproximatorInterface<BS>& approximator, 
                                                             const SatisfierInterface<BS>& satisfier)
  : _bounder(new StandardBounder<R>()),
    _integrator(integrator.clone()),
    _applicator(applicator.clone()),
    _approximator(approximator.clone()),
    _satisfier(satisfier.clone()),
    _parameters(parameters.clone())
{
}



template<class BS>
Evaluation::SetBasedHybridOrbiter<BS>*
Evaluation::SetBasedHybridOrbiter<BS>::clone() const 
{
  return new SetBasedHybridOrbiter<BS>(*this);
}


template<class BS>
Geometry::HybridGrid<typename BS::real_type>
Evaluation::SetBasedHybridOrbiter<BS>::grid(const Geometry::HybridSpace& loc) const
{
  return Geometry::HybridGrid<R>(loc,this->_parameters->grid_length());
}

template<class BS>
Geometry::HybridGridCellListSet<typename BS::real_type>
Evaluation::SetBasedHybridOrbiter<BS>::outer_approximation(const Geometry::HybridListSet<BS>& hls, const Geometry::HybridGrid<R>& hg) const
{
  return Geometry::outer_approximation(hls,hg);
}




template<class BS>
Geometry::HybridListSet<BS> 
Evaluation::SetBasedHybridOrbiter<BS>::evolution(Geometry::HybridListSet<BS>& result,
                                                 std::vector< tuple<Q,Z,DS,BS> >& working,
                                                 const System::HybridAutomaton<R>& automaton, 
                                                 const Geometry::HybridListSet<BS>& initial_set, 
                                                 const Numeric::Rational& time,
                                                 const Semantics semantics,
                                                 const EvolutionType evolution_type) const
{
  typedef Geometry::DiscreteState DS;
  std::vector< tuple<Q,Z,DS,BS> > working;
  Geometry::HybridListSet<BS> result;

  Q h=0; 
  Q t=0; 
  Z n=0;
  DS ds(0);
  BS bs;
  Bx bb;

  // Make a list of working sets
  for(typename Geometry::HybridListSet<BS>::const_iterator loc_iter=initial_set.begin();
      loc_iter!=initial_set.end(); ++loc_iter)
  {
    ds=loc_iter->discrete_state();
    const Geometry::ListSet<BS>& ls=loc_iter->continuous_state_set();
    for(uint i=0; i!=ls.size(); ++i) {
      working.push_back(make_tuple(t,n,ds,ls[i]));
    }
  }

  // Process working sets until empty
  while(!working.empty()) {
    make_ltuple(t,n,ds,bs)=working.back();
    working.pop_back();
    const DM& mode=automaton.mode(ds);
    reference_vector<const DT> transitions=automaton.transitions(ds);
    const VF& vf=mode.dynamic();
    make_lpair(h,bb)=this->flow_bounds(vf,this->bounding_box(bs));
    BS ebs=this->integration_step(vf,bs,h,bb);
    BS rbs=this->reachability_step(vf,bs,h,bb);
    
    typedef typename reference_vector<const DT>::const_iterator transitions_const_iterator;

    for(transitions_const_iterator transition_iter=transitions.begin();
        transition_iter!=transitions.end(); ++transition_iter)
    {
      Z imn=n+1;
      const DT& transition=*transition_iter;
      DS imds=transition.destination().discrete_state();
      if(semantics==upper_semantics) {
        if(possibly(this->intersects(rbs,transition.activation()))) {
          BS imbs=this->apply(transition.reset(),rbs);
          working.push_back(make_tuple(t,imn,imds,imbs));
        }
      } else if(semantics==lower_semantics) {
        if(definitely(this->subset(bs, transition.activation()))) {
          BS imbs=this->apply(transition.reset(),bs);
          working.push_back(make_tuple(t,imn,imds,imbs));
        }
      }
    }

    t+=h;
    if(semantics==upper_semantics) {
      if(possibly(this->intersects(ebs,mode.invariant()))) {
        working.push_back(make_tuple(t,n,ds,ebs));
      }
    } else if(semantics==lower_semantics) {
      if(definitely(this->subset(rbs,mode.invariant()))) {
        working.push_back(make_tuple(t,n,ds,ebs));
      }
    }

    if( semantics==upper_semantics || definitely(this->subset(rbs,mode.invariant())) ) {
      if( evolution_type==evolve ) {
        if(t==time) {
          result.adjoin(make_pair(ds,ebs));
        }
      } else if( evolution_type==reach ) {
        result.adjoin(make_pair(ds,rbs));
      }
    }
  }

  return result;

}


template<class BS>
Geometry::HybridGridCellListSet<typename BS::real_type> 
Evaluation::SetBasedHybridOrbiter<BS>::upper_evolve(const System::HybridAutomaton<R>& f, 
                                                    const Geometry::HybridGridCell<R>& hgc, 
                                                    const Numeric::Rational& t) const
{
  typedef Geometry::DiscreteState DS;
  DS ds=hgc.discrete_state();
  BS bs=this->over_approximation(hgc.continuous_state_set());
  Geometry::HybridListSet<BS> hls(f.locations());
  hls.adjoin(ds,bs);
  
  hls=this->evolution(f,hls,t,upper_semantics,evolve);

  return outer_approximation(hls,this->grid(f.locations()));
}

template<class BS>
Geometry::HybridGridCellListSet<typename BS::real_type> 
Evaluation::SetBasedHybridOrbiter<BS>::upper_reach(const System::HybridAutomaton<R>& f, 
                                                const Geometry::HybridGridCell<R>& gc, 
                                                const Numeric::Rational& t) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  /*
  const Geometry::HybridGrid<R>& grid=gc.grid();
  BS bs=this->over_approximation(Geometry::HybridBox<R>(gc));
  std::vector<BSL> orbit=this->_orbit(f,bs,n);
  GCLS result(grid);
  for(size_type i=0; i!=orbit.size(); ++i) {
    result.adjoin(this->outer_approximation(orbit[n],grid));
  }
  result.unique_sort();
  return result;
  */
}


template<class BS>
Geometry::HybridBox<typename BS::real_type>
Evaluation::SetBasedHybridOrbiter<BS>::lower_evolve(const System::HybridAutomaton<R>& f, 
                                         const Geometry::HybridBox<R>& bx, 
                                         const Numeric::Rational& t) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  /*
  BS bs=this->over_approximation(bx);
  std::vector<BSL> orbit=this->_orbit(f,bs,n);
  return this->bounding_box(orbit[n]);
  */
}


template<class BS>
Geometry::HybridBoxListSet<typename BS::real_type>
Evaluation::SetBasedHybridOrbiter<BS>::lower_reach(const System::HybridAutomaton<R>& f, 
                                        const Geometry::HybridBox<R>& bx, 
                                        const Numeric::Rational& t) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  /*
  BS bs=this->over_approximation(bx);
  std::vector<BSL> orbit=this->_orbit(f,bs,n);
  Geometry::HybridBoxListSet<R> result;
  for(size_type i=0; i!=orbit.size(); ++i) {
    result.adjoin(this->bounding_box(orbit[i]));
  }
  return result;
  */
} 




}





