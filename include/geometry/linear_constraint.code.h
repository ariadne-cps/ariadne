/***************************************************************************
 *            linear_constraint.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  pieter.collins@cwi.nl
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
 
#include <limits>

#include "constraint.h"
#include "numeric/interval.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/point.h"
#include "geometry/box.h"


namespace Ariadne {

    




template<class R>
LinearConstraint<R>::~LinearConstraint()
{
}


template<class R>
LinearConstraint<R>::LinearConstraint(const Vector<R> a, Comparison cmp, const R& b)
  : _a(a), _b(b), _c(cmp)
{
}

template<class R>
LinearConstraint<R>::LinearConstraint(const Vector<R> a, const R& b)
  : _a(a), _b(b), _c(less)
{
}


template<class R>
LinearConstraint<R>* 
LinearConstraint<R>::clone() const
{
  return new LinearConstraint<R>(this->_a, this->_c, this->_b);
}


template<class R>
dimension_type
LinearConstraint<R>::dimension() const 
{
  return this->_a.size();
}


template<class R>
smoothness_type 
LinearConstraint<R>::smoothness() const 
{
  return std::numeric_limits<smoothness_type>::max();
}


template<class R>
Comparison
LinearConstraint<R>::comparison() const 
{
  return this->_c;
}


template<class R>
Polyhedron<R>
LinearConstraint<R>::polyhedron() const 
{
  if(this->_c==less) {
    return Polyhedron<R>(Matrix<R>(1,this->dimension(),_a.begin()),
                         Vector<R>(1,&_b));
  } else {
    return Polyhedron<R>(-Matrix<R>(1,this->dimension(),_a.begin()),
                         -Vector<R>(1,&_b));
  }
}





template<class R>
std::ostream& 
LinearConstraint<R>::write(std::ostream& os) const
{
  //return os << "LinearConstraint( a=" << this->_a << ", b=" << this->_b << ", c='" << (this->_c==less ? "<" : ">") << "' )";
  return os << this->_a << ".x"<< (this->_c==less ? "<" : ">") << "=" << this->_b;
}

      
template<class R>
typename LinearConstraint<R>::A
LinearConstraint<R>::value(const Point<A>& pt) const
{
  return inner_product(pt.position_vector(),Vector<A>(this->_a)) - this->_b;
}


template<class R>
Vector<typename LinearConstraint<R>::A>
LinearConstraint<R>::gradient(const Point<A>& pt) const
{
  return Vector<A>(this->_a);
}



} // namespace Ariadne
