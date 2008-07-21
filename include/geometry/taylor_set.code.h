/***************************************************************************
 *            taylor_set.code.h
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
 
#include <iostream>
#include <vector>
#include <algorithm>

#include "taylor_set.h"

#include "exceptions.h"
#include "base/array.h"
#include "base/tuple.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"

#include "output/logging.h"

namespace Ariadne {
    
extern int verbosity; 
   
template<class X>
std::ostream& operator<<(std::ostream& os, const SparseDifferentialVector<X>& sdv);
    
template<class R>
TaylorSet<R>::TaylorSet() 
  : _model()
{
}
    
template<class R>
TaylorSet<R>::TaylorSet(dimension_type d) 
  : _model(ApproximateTaylorModel<R>::constant(Vector<I>(0u,I(-1,1)),Vector<R>(0u,0),Vector<A>(d,A(0)),1u,0u))
{
}
    
template<class R>
TaylorSet<R>::TaylorSet(const ApproximateTaylorModel<R>& model) 
  : _model(model)
{
  uint as=model.argument_size();
  ARIADNE_ASSERT(model.centre()==Vector<R>(as,R(0)));
  ARIADNE_ASSERT(model.domain()==Vector<I>(as,I(-1,1)));
}
    
template<class R>
TaylorSet<R>::TaylorSet(const Box<R>& box) 
  : _model(ApproximateTaylorModel<R>::scaling(Vector<I>(box.dimension(),I(-1,1)),
                                              Vector<R>(box.dimension(),R(0)),
                                              box.position_vectors(),
                                              _max_degree,0u))
{
}
    

template<class R>
TaylorSet<R>::TaylorSet(const Zonotope<R>& z) 
  : _model(ApproximateTaylorModel<R>::affine(Vector<I>(z.number_of_generators(),I(-1,1)),
                                             Vector<R>(z.number_of_generators(),R(0)),
                                             z.centre().position_vector(),
                                             z.generators(),
                                             _max_degree,0u))
{
}
    

template<class R>
dimension_type
TaylorSet<R>::dimension() const
{
  return this->_model.result_size();
}

template<class R>
tribool
TaylorSet<R>::empty() const
{
  return this->_model.argument_size()==0 ? indeterminate : false;
}

template<class R>
tribool
TaylorSet<R>::bounded() const
{
  return true;
}


template<class R>
R
TaylorSet<R>::radius() const
{
  return this->bounding_box().radius();
}


template<class R>
Point<R>
TaylorSet<R>::centre() const
{
  return Point<R>(this->_model.evaluate(Vector<A>(this->_model.centre())));
}

template<class R>
Box<R>
TaylorSet<R>::bounding_box() const
{
  return Box<R>(this->_model.evaluate(this->_model.domain()));
}

    
template<class R> 
tribool
TaylorSet<R>::superset(const Box<R>& bx) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

  
template<class R> 
tribool
TaylorSet<R>::intersects(const Box<R>& bx) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

  
template<class R> 
tribool
TaylorSet<R>::disjoint(const Box<R>& bx) const
{
  return Ariadne::disjoint(zonotope_over_approximation(*this),bx);
}

  
template<class R> 
tribool
TaylorSet<R>::subset(const Box<R>& bx) const
{
  return Ariadne::subset(this->bounding_box(),bx);
}

  
    
template<class R> 
std::ostream&
TaylorSet<R>::write(std::ostream& os) const
{
  return os << "TaylorSet( expansion=" << this->_model.expansion() << " )";
}


template<class R>
void
TaylorSet<R>::_instantiate() 
{
  Box<R>* bx = 0;
  TaylorSet<R>* ts = 0;
  zonotope_over_approximation(*ts);
}
    

template<class R> 
Zonotope<R> 
zonotope_over_approximation(const TaylorSet<R>& ts)  
{
  typedef typename traits<R>::interval_type I;

  const ApproximateTaylorModel<R>& model=ts.model();
  uint n=model.result_size();
  uint m=model.argument_size();
  Vector<I> zc(n);
  Matrix<R> zg(n,m);
  make_lpair(zc,zg) = affine_model(model);
  Zonotope<R> z(Point<I>(zc),zg);
  return z;
}


} // namespace Ariadne
