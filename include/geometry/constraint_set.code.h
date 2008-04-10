/***************************************************************************
 *            constraint_set.code.h
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
 
#include "constraint_set.h"
#include "linear_algebra/vector.h"
#include "geometry/euclidean_space.h"
#include "geometry/point.h"
#include "geometry/box.h"

#include "function/identity_function.h"

namespace {

using namespace Ariadne;

template<class X1, class X2> inline
tribool less(const LinearAlgebra::Vector<X1>& v1, const LinearAlgebra::Vector<X2>& v2)
{
  tribool result=true;
  for(size_type i=0; i!=v1.size(); ++i) {
    result=result && (v1[0]<v2[0]);
    if(!result) { return false; }
  }
  return result;
}
  
}


namespace Ariadne {
    
template<class R>
Geometry::ConstraintSet<R>::ConstraintSet(const Geometry::Box<R>& bx)
  : _function_ptr(new Function::IdentityFunction<R>(bx.dimension()))
  , _codomain(bx) 
{ 
}

template<class R>
Geometry::ConstraintSet<R>::ConstraintSet(const Function::FunctionInterface<R>& f, 
                                          const Geometry::Box<R>& bx)
  : _function_ptr(f.clone())
  , _codomain(bx) 
{ 
}



template<class R>
Geometry::ConstraintSet<R>::~ConstraintSet() 
{ 
}


template<class R>
Geometry::ConstraintSet<R>* 
Geometry::ConstraintSet<R>::clone() const 
{
  return new ConstraintSet<R>(*this->_function_ptr,this->_codomain); 
}


template<class R>
dimension_type 
Geometry::ConstraintSet<R>::dimension() const 
{
  return this->_function_ptr->argument_size(); 
}

template<class R>
Geometry::EuclideanSpace
Geometry::ConstraintSet<R>::space() const 
{
  return EuclideanSpace(this->dimension());
}




template<class R>
tribool 
Geometry::ConstraintSet<R>::contains(const Point<R>& pt) const 
{
  const Box<R>& dom=this->_codomain;
  const Function::FunctionInterface<R>& f=*this->_function_ptr;
  return dom.contains(Point<A>(f(pt.position_vector())));
}


template<class R>
tribool 
Geometry::ConstraintSet<R>::contains(const Point<A>& pt) const 
{
  const Box<R>& dom=this->_codomain;
  const Function::FunctionInterface<R>& f=*this->_function_ptr;
  return dom.contains(Point<A>(f(pt.position_vector())));
}


template<class R>
tribool 
Geometry::ConstraintSet<R>::superset(const Box<R>& bx) const 
{
  const Box<R>& dom=this->_codomain;
  const Function::FunctionInterface<R>& f=*this->_function_ptr;
  return Geometry::subset(Box<R>(f(bx.position_vectors())),dom);
}


template<class R>
tribool 
Geometry::ConstraintSet<R>::intersects(const Box<R>& bx) const 
{
  return !this->disjoint(bx);
}


template<class R>
tribool 
Geometry::ConstraintSet<R>::disjoint(const Box<R>& bx) const 
{
  const Box<R>& dom=this->_codomain;
  const Function::FunctionInterface<R>& f=*this->_function_ptr;
  return Geometry::disjoint(dom,Box<R>(f(bx.position_vectors())));
}


template<class R>
tribool 
Geometry::ConstraintSet<R>::subset(const Box<R>& bx) const 
{
  return indeterminate;
}


template<class R>
tribool 
Geometry::ConstraintSet<R>::bounded() const 
{
  return indeterminate;
}      


template<class R>
Geometry::Box<R> 
Geometry::ConstraintSet<R>::bounding_box() const 
{
  if(dynamic_cast<const Function::IdentityFunction<R>*>(&*this->_function_ptr)) {
    return this->_codomain;
  } else {
    throw UnboundedSet("ConstraintSet::bounding_box(): cannot be computed in general case");
  }
}      


template<class R>
std::ostream& 
Geometry::ConstraintSet<R>::write(std::ostream& os) const 
{
  return os << "ConstraintSet( function=" << *this->_function_ptr
            << ", codomain=" << this->_codomain << ")";
}



template<class R>
size_type
Geometry::ConstraintSet<R>::number_of_constraints() const
{
  return this->_function_ptr->result_size();
}


template<class R>
const Function::FunctionInterface<R>&
Geometry::ConstraintSet<R>::function() const 
{
  return *this->_function_ptr;
}

template<class R>
const Geometry::Box<R>&
Geometry::ConstraintSet<R>::codomain() const 
{
  return this->_codomain;
}



}
