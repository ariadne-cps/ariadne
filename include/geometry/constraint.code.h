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
#include "geometry/rectangle.h"
#include "geometry/polyhedron.h"
#include "geometry/zonotope.h"


namespace Ariadne {

    



template<class R>
Geometry::DifferentiableConstraint<R>::DifferentiableConstraint(const Function::DifferentiableFunctionInterface<R>& f, Comparison cmp)
  : _function_ptr(f.clone()), _comparison(cmp)
{ 
}


template<class R>
Geometry::DifferentiableConstraint<R>::~DifferentiableConstraint() 
{ 
}


template<class R>
Geometry::DifferentiableConstraint<R>* 
Geometry::DifferentiableConstraint<R>::clone() const 
{
  return new DifferentiableConstraint<R>(*this->_function_ptr); 
}


template<class R>
dimension_type 
Geometry::DifferentiableConstraint<R>::dimension() const 
{
  return this->_function_ptr->argument_size(); 
}


template<class R>
smoothness_type 
Geometry::DifferentiableConstraint<R>::smoothness() const 
{
  return this->_function_ptr->smoothness(); 
}





template<class R>
std::ostream& 
Geometry::DifferentiableConstraint<R>::write(std::ostream& os) const 
{
  return os << "DifferentiableConstraint( ... )";
  return os << "DifferentiableConstraint( function=" << *this->_function_ptr << ", comparison=" << (this->_comparison==less ? "<" : ">") << " )";
}



template<class R>
Geometry::Comparison
Geometry::DifferentiableConstraint<R>::comparison() const 
{
  return this->_comparison;
}


template<class R>
const Function::DifferentiableFunctionInterface<R>&
Geometry::DifferentiableConstraint<R>::function() const 
{
  return *this->_function_ptr;
}


template<class R>
typename Geometry::DifferentiableConstraint<R>::A 
Geometry::DifferentiableConstraint<R>::value(const Point<A>& pt) const
{
  return this->_function_ptr->operator()(pt.position_vector())[0];
}


template<class R>
LinearAlgebra::Vector<typename Geometry::DifferentiableConstraint<R>::A>
Geometry::DifferentiableConstraint<R>::gradient(const Point<A>& pt) const
{
  return LinearAlgebra::Vector<A>(this->dimension(),this->_function_ptr->jacobian(pt.position_vector()).begin());
}




template<class R>
void
Geometry::DifferentiableConstraint<R>::instantiate() 
{
  /*
  typedef Numeric::Interval<R> I;
  Box<R>* r=0;
  Zonotope<R,R>* z=0;
  Zonotope<I,R>* ez=0;
  Zonotope<I,I>* iz=0;
  ConstraintInterface<R>* ci=0;
  DifferentiableConstraint<R>* c=0;
  
  satisfies(*r,*c);
  satisfies(*z,*c);
  satisfies(*ez,*c);
  satisfies(*r,*ci);
  satisfies(*iz,*ci);
  */
}







/*

template<class R>
Numeric::Interval<R> 
Geometry::value(const DifferentiableConstraintInterface<R>& c, const Box<R>& r)
{
  typedef typename Numeric::traits<R>::interval_type I;
  return c.value(Point<I>(r));
}


template<class R>
Numeric::Interval<R> 
Geometry::value(const DifferentiableConstraintInterface<R>& c, const Zonotope< Numeric::Interval<R> >& z)
{
  typedef typename Numeric::traits<R>::interval_type I;
  LinearAlgebra::Vector<I> e(z.number_of_generators(),I(-1,1));
  Geometry::Point<I> bpt(z.bounding_box());
  return c.value(z.centre())+LinearAlgebra::inner_product(c.gradient(bpt),z.generators()*e);
}


template<class R>
tribool 
Geometry::satisfies(const Box<R>& r, const DifferentiableConstraintInterface<R>& c)
{
  typedef typename Numeric::traits<R>::interval_type I;
  Point<I> pt(r);
  I v=c.value(pt);
  return ::compare_zero(v,c.comparison());
}

template<class R>
tribool 
Geometry::satisfies(const Box<R>& r, const DifferentiableConstraint<R>& c)
{
  typedef typename Numeric::traits<R>::interval_type I;
  Point<I> pt(r);
  I v=c.value(pt);
  return ::compare_zero(v,c.comparison());
}

template<class R>
tribool 
Geometry::satisfies(const Zonotope<R,R>& z, const DifferentiableConstraint<R>& c)
{
  typedef Numeric::Interval<R> I;
  LinearAlgebra::Vector<I> e(z.number_of_generators(),I(-1,1));
  const Point<R>& zc=z.centre();
  const LinearAlgebra::Matrix<R>& zG=z.generators();
  Box<R> bb=z.bounding_box();
  Point<I> bpt(bb);
  Numeric::Interval<R> v=c.value(zc)+LinearAlgebra::inner_product(c.gradient(bb)*zG,e);
  return ::compare_zero(v,c.comparison());
}
      
template<class R>
tribool 
Geometry::satisfies(const Zonotope<Numeric::Interval<R>,R>& z, const DifferentiableConstraint<R>& c)
{
  typedef Numeric::Interval<R> I;
  LinearAlgebra::Vector<I> e(z.number_of_generators(),I(-1,1));
  const Point<I>& zc=z.centre();
  const LinearAlgebra::Matrix<R>& zG=z.generators();
  Box<R> bb=z.bounding_box();
  Point<I> bpt(bb);
  Numeric::Interval<R> v=c.value(zc)+LinearAlgebra::inner_product(c.gradient(bb)*zG,e);
  return ::compare_zero(v,c.comparison());
}
      

template<class R>
tribool 
Geometry::satisfies(const Zonotope< Numeric::Interval<R>, Numeric::Interval<R> >& z, 
                    const DifferentiableConstraintInterface<R>& c)
{
  typedef Numeric::Interval<R> I;
  LinearAlgebra::Vector<I> e(z.number_of_generators(),I(-1,1));
  const Point<I>& zc=z.centre();
  const LinearAlgebra::Matrix<I>& zG=z.generators();
  Box<R> bb=z.bounding_box();
  Point<I> bpt(bb);
  Numeric::Interval<R> v=c.value(zc)+LinearAlgebra::inner_product(c.gradient(bb)*zG,e);
  return ::compare_zero(v,c.comparison());
}
      
*/

} // namespace Ariadne
