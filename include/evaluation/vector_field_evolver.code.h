/***************************************************************************
 *            vector_field_evolver.code.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>
#include <cstring>
#include <cassert>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "base/array.h"
#include "base/tuple.h"

#include "numeric/rational.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/box.h"
#include "geometry/box_list_set.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"

#include "system/vector_field.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/approximator_interface.h"

#include "evaluation/standard_approximator.h"
#include "evaluation/standard_subdivider.h"

#include "output/logging.h"

#include "vector_field_evolver.h"

namespace {
    
using namespace Ariadne;


template<class BS, class BST>
struct pair_first_less {
  bool operator()(const std::pair<BS,BST>& ts1, const std::pair<BS,BST>& ts2) {
    return ts1.first < ts2.first; 
  }
};

/*!\ brief A class representing pre-computed bounds for an integration step. */
template<class BS>
class IntegrationStepBound {
  typedef typename BS::real_type R;
 public:
  /*!\ brief Constructor. */
  IntegrationStepBound(const Geometry::Box<R>& bound, const BS& integration_time) 
    : _bound(bound), _integration_time(integration_time) { }
  /*!\ brief Constructor. */
  IntegrationStepBound(const Geometry::Box<R>& bound, const Numeric::Interval<BS>& integration_time) 
    : _bound(bound), _integration_time(integration_time.upper()) { }
  /*! The spacial bound for the integrations step. */
  const Geometry::Box<R>& bound() const { return _bound; }
  /*! The step size in time of the integrations step. */
  const BS& integration_time() const { return _integration_time; }
 private:
  Geometry::Box<R> _bound;
  BS _integration_time;
};

}


namespace Ariadne {
    
namespace Evaluation { static int& verbosity = integrator_verbosity; }
    

using namespace Geometry;
using namespace System;

template<class BS>
Evaluation::VectorFieldEvolver<BS>::~VectorFieldEvolver()
{
}

template<class BS> 
Evaluation::VectorFieldEvolver<BS>::VectorFieldEvolver(const EvolutionParameters<R>& parameters, 
                                                      const IntegratorInterface<BS>& integrator)
  : _parameters(parameters.clone()),
    _integrator(integrator.clone()),
    _approximator(new StandardApproximator<BS>()),
    _subdivider(new StandardSubdivider<BS>)
{
}

template<class BS> 
Evaluation::VectorFieldEvolver<BS>::VectorFieldEvolver(const EvolutionParameters<R>& parameters, 
                                                      const IntegratorInterface<BS>& integrator,
                                                      const ApproximatorInterface<BS>& approximator)
  : _parameters(parameters.clone()),
    _integrator(integrator.clone()),
    _approximator(approximator.clone()),
    _subdivider(new StandardSubdivider<BS>)
{
}


template<class BS>
Evaluation::VectorFieldEvolver<BS>*
Evaluation::VectorFieldEvolver<BS>::clone() const
{
  return new VectorFieldEvolver<BS>(*this);
}




template<class BS>
Evaluation::EvolutionParameters<typename BS::real_type>&
Evaluation::VectorFieldEvolver<BS>::parameters() 
{
  return *this->_parameters;
}

template<class BS>
const Evaluation::EvolutionParameters<typename BS::real_type>&
Evaluation::VectorFieldEvolver<BS>::parameters() const
{
  return *this->_parameters;
}





template<class BS>
void
Evaluation::VectorFieldEvolver<BS>::_step(BSL& evolve,
                                          BSL& reach, 
                                          TBSL& working,
                                          const VF& vf,
                                          const T& time,
                                          Semantics semantics) const
{
  TBS tbs=working.pop();
  BS const& bs=tbs.set();
  if(tbs.time()==time) {
    evolve.adjoin(bs);
  } else if(radius(bs) > maximum_basic_set_radius()) {
    if(semantics==upper_semantics) {
      working.adjoin(this->subdivide(tbs));
    } 
  } else {
    Bx bb; T h; BS rbs; 
    make_lpair(h,bb)=this->flow_bounds(vf,this->bounding_box(bs));
    tbs=this->integration_step(vf,tbs,h,bb);
    {
      rbs=this->reachability_step(vf,bs,h,bb);
      reach.adjoin(rbs);
    }
    working.adjoin(tbs);
  }
}



template<class BS>
Geometry::GridCellListSet<typename BS::real_type>
Evaluation::VectorFieldEvolver<BS>::_upper_evolve(const System::VectorField<R>& vector_field, 
                                                  const Geometry::GridCellListSet<R>& initial_set,
                                                  const Numeric::Rational& time) const
{
  BSL reach, evolve;
  TBSL working=this->timed_basic_set_list(initial_set);
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,vector_field,time,upper_semantics);
  }
  return this->outer_approximation(evolve,initial_set.grid());
}


template<class BS>
Geometry::GridCellListSet<typename BS::real_type>
Evaluation::VectorFieldEvolver<BS>::_upper_reach(const System::VectorField<R>& vector_field, 
                                                 const Geometry::GridCellListSet<R>& initial_set,
                                                 const Numeric::Rational& time) const
{
  BSL reach, evolve;
  TBSL working=this->timed_basic_set_list(initial_set);
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,vector_field,time,upper_semantics);
  }
  return this->outer_approximation(reach,initial_set.grid());
}








template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::VectorFieldEvolver<BS>::upper_evolve(const System::VectorField<R>& vector_field, 
                                                 const Geometry::SetInterface<R>& initial_set, 
                                                 const Numeric::Rational& time) const
{
  BSL reach, evolve;
  TBSL working=this->timed_basic_set_list(this->outer_approximation(initial_set));
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,vector_field,time,upper_semantics);
  }
  return new GMS(this->outer_approximation(evolve,this->grid(vector_field.dimension())));
}

template<class BS>
Geometry::SetInterface<typename BS::real_type>* 
Evaluation::VectorFieldEvolver<BS>::upper_reach(const System::VectorField<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_set,
                                               const Numeric::Rational& time) const
{
  std::cerr<<__FUNCTION__<<std::endl;
  BSL reach, evolve;
  std::cerr<<"outer_approximation(initial_set)="<<this->outer_approximation(initial_set)<<std::endl;
  TBSL working=this->timed_basic_set_list(this->outer_approximation(initial_set));
  std::cerr<<"working="<<working<<std::endl;
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,vector_field,time,upper_semantics);
  }
  std::cerr << "reach="<<reach<<"\nevolve="<<evolve<<std::endl;
  return new GMS(this->outer_approximation(reach,this->grid(vector_field.dimension())));
}


template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::VectorFieldEvolver<BS>::lower_evolve(const System::VectorField<R>& vector_field,
                                                 const Geometry::SetInterface<R>& initial_set,
                                                 const Numeric::Rational& time) const
{
  BSL reach, evolve;
  TBSL working=this->timed_basic_set_list(this->lower_approximation(initial_set));
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,vector_field,time,lower_semantics);
  }
  return new GMS(this->lower_approximation(evolve,this->grid(vector_field.dimension())));
}

template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::VectorFieldEvolver<BS>::lower_reach(const System::VectorField<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_set,
                                               const Numeric::Rational& time) const
{
  BSL reach, evolve;
  TBSL working=this->timed_basic_set_list(this->lower_approximation(initial_set));
  while(working.size()!=0) { 
    this->_step(evolve,reach,working,vector_field,time,lower_semantics);
  }
  return new GMS(this->lower_approximation(reach,this->grid(vector_field.dimension())));
}





template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::VectorFieldEvolver<BS>::lower_reach(const System::VectorField<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_sete) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::VectorFieldEvolver<BS>::chainreach(const System::VectorField<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_set) const
{
  Bx bb=this->bounding_domain(vector_field);
  Gr grid=this->grid(vector_field.dimension());
  T time=this->lock_to_grid_time();
  GMS* result=new GMS(grid,bb);
  GCLS found=this->outer_approximation(initial_set);
  found=this->_upper_reach(vector_field,found,time);
  while(!found.empty()) {
    result->adjoin(found);
    found=this->_upper_evolve(vector_field,found,time);
    found.remove(*result);
  }
  return result;
}


template<class BS>
Geometry::SetInterface<typename BS::real_type>*
Evaluation::VectorFieldEvolver<BS>::viable(const System::VectorField<R>& vector_field,
                                          const Geometry::SetInterface<R>& bounding_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class BS>
tribool
Evaluation::VectorFieldEvolver<BS>::verify(const System::VectorField<R>& vector_field,
                                          const Geometry::SetInterface<R>& initial_set,
                                          const Geometry::SetInterface<R>& safe_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}










}

