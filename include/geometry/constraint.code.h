/***************************************************************************
 *            constraint.code.h
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
 
#include "constraint.h"
#include "../numeric/interval.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"

namespace {
 
using namespace Ariadne;

template<class R> inline 
tribool compare_zero(const Numeric::Interval<R>& ivl, Geometry::Comparison cmp) {
  if(ivl.upper()<0) {
    return (cmp==Geometry::less);
  } else if(ivl.lower()>0) {
    return (cmp==Geometry::greater);
  } else {
    return indeterminate;
  }
}

}

namespace Ariadne {
    
template<class R>
Geometry::Constraint<R>::Constraint(const System::FunctionInterface<R>& f, Comparison cmp)
  : _function_ptr(f.clone()), _comparison(cmp)
{ 
}


template<class R>
Geometry::Constraint<R>::~Constraint() 
{ 
}


template<class R>
Geometry::Constraint<R>* 
Geometry::Constraint<R>::clone() const 
{
  return new Constraint<R>(*this->_function_ptr); 
}


template<class R>
dimension_type 
Geometry::Constraint<R>::dimension() const 
{
  return this->_function_ptr->argument_size(); 
}



template<class R>
tribool 
Geometry::satisfies(const Rectangle<R>& r, const Constraint<R>& c)
{
  typedef typename Numeric::traits<R>::interval_type I;
  Point<I> pt(r);
  I v=c.value(pt);
  return ::compare_zero(v,c.comparison());
}

template<class R>
tribool 
Geometry::satisfies(const Zonotope<R,R>& z, const Constraint<R>& c)
{
  typedef Numeric::Interval<R> I;
  LinearAlgebra::Vector<I> e(z.number_of_generators(),I(-1,1));
  const Point<R>& zc=z.centre();
  const LinearAlgebra::Matrix<R>& zG=z.generators();
  Rectangle<R> bb=z.bounding_box();
  Point<I> bpt(bb);
  Numeric::Interval<R> v=c.value(zc)+LinearAlgebra::inner_product(c.gradient(bb)*zG,e);
  return ::compare_zero(v,c.comparison());
}
      
template<class R>
tribool 
Geometry::satisfies(const Zonotope<Numeric::Interval<R>,R>& z, const Constraint<R>& c)
{
  typedef Numeric::Interval<R> I;
  LinearAlgebra::Vector<I> e(z.number_of_generators(),I(-1,1));
  const Point<I>& zc=z.centre();
  const LinearAlgebra::Matrix<R>& zG=z.generators();
  Rectangle<R> bb=z.bounding_box();
  Point<I> bpt(bb);
  Numeric::Interval<R> v=c.value(zc)+LinearAlgebra::inner_product(c.gradient(bb)*zG,e);
  return ::compare_zero(v,c.comparison());
}
      

template<class R>
std::ostream& 
Geometry::Constraint<R>::write(std::ostream& os) const 
{
  return os << "Constraint( function=" << *this->_function_ptr << ", comparison=" << (this->_comparison==less ? "<" : ">") << " )";
}



template<class R>
const Geometry::Comparison&
Geometry::Constraint<R>::comparison() const 
{
  return this->_comparison;
}


template<class R>
const System::FunctionInterface<R>&
Geometry::Constraint<R>::function() const 
{
  return *this->_function_ptr;
}


template<class R>
typename Geometry::Constraint<R>::A 
Geometry::Constraint<R>::value(const Point<A>& pt) const
{
  return this->_function_ptr->operator()(pt.position_vector())[0];
}


template<class R>
LinearAlgebra::Vector<typename Geometry::Constraint<R>::A>
Geometry::Constraint<R>::gradient(const Point<A>& pt) const
{
  return LinearAlgebra::Vector<A>(this->dimension(),this->_function_ptr->jacobian(pt.position_vector()).begin());
}




template<class R>
void
Geometry::Constraint<R>::instantiate() 
{
  Rectangle<R>* r=0;
  Zonotope<R,R>* z=0;
  Zonotope<Numeric::Interval<R>,R>* ez=0;
  Constraint<R>* c=0;
  
  satisfies(*r,*c);
  satisfies(*z,*c);
  satisfies(*ez,*c);
}


}
