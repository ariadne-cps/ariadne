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

#include "numeric/rational.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/rectangular_set.h"

#include "system/vector_field.h"
#include "system/affine_vector_field.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/bounder.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/lohner_integrator.h"

#include "output/logging.h"

#include "vector_field_orbiter_interface.h"
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

template<class R>
Evaluation::VectorFieldEvolver<R>::~VectorFieldEvolver()
{
}

template<class R>
Evaluation::VectorFieldEvolver<R>::VectorFieldEvolver(const VectorFieldEvolver<R>& i)
  : _parameters(i._parameters->clone()),
    _orbiter(i._orbiter->clone())
{
}


template<class R>
Evaluation::VectorFieldEvolver<R>*
Evaluation::VectorFieldEvolver<R>::clone() const
{
  return new VectorFieldEvolver<R>(*this);
}

template<class R>
Evaluation::VectorFieldEvolver<R>::VectorFieldEvolver(const EvolutionParameters<R>& parameters)
  : _parameters(new EvolutionParameters<R>(parameters)),
    _orbiter()
{
}




template<class R>
Evaluation::EvolutionParameters<R>&
Evaluation::VectorFieldEvolver<R>::parameters() 
{
  return *this->_parameters;
}

template<class R>
const Evaluation::EvolutionParameters<R>&
Evaluation::VectorFieldEvolver<R>::parameters() const
{
  return *this->_parameters;
}






template<class R> inline
Geometry::GridCellListSet<R>
Evaluation::VectorFieldEvolver<R>::evolve(const System::VectorFieldInterface<R>& vector_field, 
                                          const Geometry::GridCell<R>& initial_set, 
                                          const Numeric::Rational& time) const
{
  return this->_orbiter->evolve(vector_field,initial_set,time);
}

template<class R> inline
Geometry::GridCellListSet<R>
Evaluation::VectorFieldEvolver<R>::reach(const System::VectorFieldInterface<R>& vector_field, 
                                         const Geometry::GridCell<R>& initial_set, 
                                         const Numeric::Rational& time) const
{
  return this->_orbiter->reach(vector_field,initial_set,time);
}




template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::evolve(const System::VectorFieldInterface<R>& vector_field, 
                                          const Geometry::SetInterface<R>& initial_set, 
                                          const Numeric::Rational& time) const
{
  Geometry::Grid<R> grid=this->grid();
  Geometry::GridMaskSet<R> result(grid);
  Geometry::GridMaskSet<R> initial(grid);
  for(typename GridMaskSet<R>::const_iterator cell_iter=initial.begin(); cell_iter!=initial.end(); ++cell_iter) {
    result.adjoin(this->evolve(vector_field,*cell_iter,time));
  }
  return new GridMaskSet<R>(result);
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::bounded_evolve(const System::VectorFieldInterface<R>& vector_field,
                                                  const Geometry::SetInterface<R>& initial_set,
                                                  const Geometry::SetInterface<R>& bounding_set,
                                                  const time_type& time) const
{
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::bounded_reach(const System::VectorFieldInterface<R>& vector_field,
                                                 const Geometry::SetInterface<R>& initial_set,
                                                 const Geometry::SetInterface<R>& bounding_set,
                                                 const time_type& time) const
{
}




template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::lower_reach(const System::VectorFieldInterface<R>& vector_field,
                                               const Geometry::SetInterface<R>& initial_set,
                                               const Geometry::SetInterface<R>& bounding_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::chainreach(const System::VectorFieldInterface<R>& vector_field, 
                                              const Geometry::SetInterface<R>& initial_set) const
{
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::chainreach(const System::VectorFieldInterface<R>& vector_field,
                                      const Geometry::SetInterface<R>& initial_set,
                                      const Geometry::Box<R>& bounding_box) const
{
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::viable(const System::VectorFieldInterface<R>& vector_field,
                                          const Geometry::SetInterface<R>& bounding_set) const
{
}


template<class R>
tribool
Evaluation::VectorFieldEvolver<R>::verify(const System::VectorFieldInterface<R>& vector_field,
                                  const Geometry::SetInterface<R>& initial_set,
                                  const Geometry::SetInterface<R>& safe_set) const
{
}










}

