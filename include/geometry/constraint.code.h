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
 
#include <limits>

#include "constraint.h"
#include "numeric/interval.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/polyhedron.h"
#include "geometry/zonotope.h"


namespace Ariadne {

    



template<class R>
Constraint<R>::Constraint(const FunctionInterface<R>& f, Comparison cmp)
  : _function_ptr(f.clone()), _comparison(cmp)
{ 
}


template<class R>
Constraint<R>::~Constraint() 
{ 
}


template<class R>
Constraint<R>* 
Constraint<R>::clone() const 
{
  return new Constraint<R>(*this->_function_ptr); 
}


template<class R>
dimension_type 
Constraint<R>::dimension() const 
{
  return this->_function_ptr->argument_size(); 
}


template<class R>
smoothness_type 
Constraint<R>::smoothness() const 
{
  return this->_function_ptr->smoothness(); 
}





template<class R>
std::ostream& 
Constraint<R>::write(std::ostream& os) const 
{
  return os << "Constraint( ... )";
  return os << "Constraint( function=" << *this->_function_ptr << ", comparison=" << (this->_comparison==less ? "<" : ">") << " )";
}



template<class R>
Comparison
Constraint<R>::comparison() const 
{
  return this->_comparison;
}


template<class R>
const FunctionInterface<R>&
Constraint<R>::function() const 
{
  return *this->_function_ptr;
}


template<class R>
typename Constraint<R>::A 
Constraint<R>::value(const Point<A>& pt) const
{
  return this->_function_ptr->operator()(pt.position_vector())[0];
}


template<class R>
Vector<typename Constraint<R>::A>
Constraint<R>::gradient(const Point<A>& pt) const
{
  return Vector<A>(this->dimension(),this->_function_ptr->jacobian(pt.position_vector()).begin());
}




template<class R>
void
Constraint<R>::instantiate() 
{
  /*
  typedef Interval<R> I;
  Box<R>* r=0;
  Zonotope<R,R>* z=0;
  Zonotope<I,R>* ez=0;
  Zonotope<I,I>* iz=0;
  ConstraintInterface<R>* ci=0;
  Constraint<R>* c=0;
  
  satisfies(*r,*c);
  satisfies(*z,*c);
  satisfies(*ez,*c);
  satisfies(*r,*ci);
  satisfies(*iz,*ci);
  */
}







/*

template<class R>
Interval<R> 
value(const ConstraintInterface<R>& c, const Box<R>& r)
{
  typedef typename traits<R>::interval_type I;
  return c.value(Point<I>(r));
}


template<class R>
Interval<R> 
value(const ConstraintInterface<R>& c, const Zonotope< Interval<R> >& z)
{
  typedef typename traits<R>::interval_type I;
  Vector<I> e(z.number_of_generators(),I(-1,1));
  Point<I> bpt(z.bounding_box());
  return c.value(z.centre())+inner_product(c.gradient(bpt),z.generators()*e);
}


template<class R>
tribool 
satisfies(const Box<R>& r, const ConstraintInterface<R>& c)
{
  typedef typename traits<R>::interval_type I;
  Point<I> pt(r);
  I v=c.value(pt);
  return ::compare_zero(v,c.comparison());
}

template<class R>
tribool 
satisfies(const Box<R>& r, const Constraint<R>& c)
{
  typedef typename traits<R>::interval_type I;
  Point<I> pt(r);
  I v=c.value(pt);
  return ::compare_zero(v,c.comparison());
}

template<class R>
tribool 
satisfies(const Zonotope<R,R>& z, const Constraint<R>& c)
{
  typedef Interval<R> I;
  Vector<I> e(z.number_of_generators(),I(-1,1));
  const Point<R>& zc=z.centre();
  const Matrix<R>& zG=z.generators();
  Box<R> bb=z.bounding_box();
  Point<I> bpt(bb);
  Interval<R> v=c.value(zc)+inner_product(c.gradient(bb)*zG,e);
  return ::compare_zero(v,c.comparison());
}
      
template<class R>
tribool 
satisfies(const Zonotope<Interval<R>,R>& z, const Constraint<R>& c)
{
  typedef Interval<R> I;
  Vector<I> e(z.number_of_generators(),I(-1,1));
  const Point<I>& zc=z.centre();
  const Matrix<R>& zG=z.generators();
  Box<R> bb=z.bounding_box();
  Point<I> bpt(bb);
  Interval<R> v=c.value(zc)+inner_product(c.gradient(bb)*zG,e);
  return ::compare_zero(v,c.comparison());
}
      

template<class R>
tribool 
satisfies(const Zonotope< Interval<R>, Interval<R> >& z, 
                    const ConstraintInterface<R>& c)
{
  typedef Interval<R> I;
  Vector<I> e(z.number_of_generators(),I(-1,1));
  const Point<I>& zc=z.centre();
  const Matrix<I>& zG=z.generators();
  Box<R> bb=z.bounding_box();
  Point<I> bpt(bb);
  Interval<R> v=c.value(zc)+inner_product(c.gradient(bb)*zG,e);
  return ::compare_zero(v,c.comparison());
}
      
*/

} // namespace Ariadne
