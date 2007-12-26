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
#include "geometry/point.h"
#include "geometry/box.h"

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
Geometry::ConstraintSet<R>::ConstraintSet(const Function::FunctionInterface<R>& f)
  : SetInterface<R>(), _function_ptr(f.clone()) 
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
  return new ConstraintSet<R>(*this->_function_ptr); 
}


template<class R>
dimension_type 
Geometry::ConstraintSet<R>::dimension() const 
{
  return this->_function_ptr->argument_size(); 
}



template<class R>
tribool 
Geometry::ConstraintSet<R>::contains(const Point<R>& pt) const 
{
  const Function::FunctionInterface<R>& f=*this->_function_ptr;
  LinearAlgebra::Vector<A> v=pt.position_vector();
  LinearAlgebra::Vector<R> z=LinearAlgebra::Vector<R>::zero(this->number_of_constraints());
  LinearAlgebra::Vector<A> r=f(v);
  return ::less(z,r);
}


template<class R>
tribool 
Geometry::ConstraintSet<R>::contains(const Point<A>& pt) const 
{
  const Function::FunctionInterface<R>& f=*this->_function_ptr;
  LinearAlgebra::Vector<A> v=pt.position_vector();
  LinearAlgebra::Vector<R> z=LinearAlgebra::Vector<R>::zero(this->number_of_constraints());
  LinearAlgebra::Vector<A> r=f(v);
  return ::less(z,r);
}


template<class R>
tribool 
Geometry::ConstraintSet<R>::superset(const Box<R>& r) const 
{
  const Function::FunctionInterface<R>& f=*this->_function_ptr;
  LinearAlgebra::Vector<A> v=r.position_vectors();
  LinearAlgebra::Vector<R> z=LinearAlgebra::Vector<R>::zero(this->number_of_constraints());
  LinearAlgebra::Vector<A> fv=f(v);
  return ::less(z,fv);
}


template<class R>
tribool 
Geometry::ConstraintSet<R>::intersects(const Box<R>& r) const 
{
  return !this->disjoint(r);
}


template<class R>
tribool 
Geometry::ConstraintSet<R>::disjoint(const Box<R>& r) const 
{
  const Function::FunctionInterface<R>& f=*this->_function_ptr;
  LinearAlgebra::Vector<A> v=r.position_vectors();
  LinearAlgebra::Vector<R> z=LinearAlgebra::Vector<R>::zero(this->number_of_constraints());
  return ::less(f(v),z);
}


template<class R>
tribool 
Geometry::ConstraintSet<R>::subset(const Box<R>& r) const 
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
  throw UnboundedSet("ConstraintSet::bounding_box(): cannot be computed in general case");
}      


template<class R>
std::ostream& 
Geometry::ConstraintSet<R>::write(std::ostream& os) const 
{
  return os << "ConstraintSet( function=" << *this->_function_ptr << " )";
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
Geometry::Point<typename Geometry::ConstraintSet<R>::A> 
Geometry::ConstraintSet<R>::function(const Point<A>& pt) const 
{
  const Function::FunctionInterface<R>& f=*this->_function_ptr;
  LinearAlgebra::Vector<A> v=pt.position_vector();
  return Point<A>(f(v));
}

}
