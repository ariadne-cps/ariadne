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
#include "../numeric/interval.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/polyhedron.h"
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
Geometry::Constraint<R>::Constraint(const Function::FunctionInterface<R>& f, Comparison cmp)
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
smoothness_type 
Geometry::Constraint<R>::smoothness() const 
{
  return this->_function_ptr->smoothness(); 
}



template<class R>
Numeric::Interval<R> 
Geometry::value(const ConstraintInterface<R>& c, const Rectangle<R>& r)
{
  typedef typename Numeric::traits<R>::interval_type I;
  return c.value(Point<I>(r));
}


template<class R>
Numeric::Interval<R> 
Geometry::value(const ConstraintInterface<R>& c, const Zonotope< Numeric::Interval<R> >& z)
{
  typedef typename Numeric::traits<R>::interval_type I;
  LinearAlgebra::Vector<I> e(z.number_of_generators(),I(-1,1));
  Geometry::Point<I> bpt(z.bounding_box());
  return c.value(z.centre())+LinearAlgebra::inner_product(c.gradient(bpt),z.generators()*e);
}


template<class R>
tribool 
Geometry::satisfies(const Rectangle<R>& r, const ConstraintInterface<R>& c)
{
  typedef typename Numeric::traits<R>::interval_type I;
  Point<I> pt(r);
  I v=c.value(pt);
  return ::compare_zero(v,c.comparison());
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
tribool 
Geometry::satisfies(const Zonotope< Numeric::Interval<R>, Numeric::Interval<R> >& z, 
                    const ConstraintInterface<R>& c)
{
  typedef Numeric::Interval<R> I;
  LinearAlgebra::Vector<I> e(z.number_of_generators(),I(-1,1));
  const Point<I>& zc=z.centre();
  const LinearAlgebra::Matrix<I>& zG=z.generators();
  Rectangle<R> bb=z.bounding_box();
  Point<I> bpt(bb);
  Numeric::Interval<R> v=c.value(zc)+LinearAlgebra::inner_product(c.gradient(bb)*zG,e);
  return ::compare_zero(v,c.comparison());
}
      

template<class R>
std::ostream& 
Geometry::Constraint<R>::write(std::ostream& os) const 
{
  return os << "Constraint( ... )";
  return os << "Constraint( function=" << *this->_function_ptr << ", comparison=" << (this->_comparison==less ? "<" : ">") << " )";
}



template<class R>
Geometry::Comparison
Geometry::Constraint<R>::comparison() const 
{
  return this->_comparison;
}


template<class R>
const Function::FunctionInterface<R>&
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
  typedef Numeric::Interval<R> I;
  Rectangle<R>* r=0;
  Zonotope<R,R>* z=0;
  Zonotope<I,R>* ez=0;
  Zonotope<I,I>* iz=0;
  ConstraintInterface<R>* ci=0;
  Constraint<R>* c=0;
  
  Geometry::value(*c,*r);
  Geometry::value(*c,*iz);

  satisfies(*r,*c);
  satisfies(*z,*c);
  satisfies(*ez,*c);
  satisfies(*r,*ci);
  satisfies(*iz,*ci);
}




template<class R>
Geometry::LinearConstraint<R>::~LinearConstraint()
{
}


template<class R>
Geometry::LinearConstraint<R>::LinearConstraint(const LinearAlgebra::Vector<R> a, Comparison cmp, const R& b)
  : _a(a), _b(b), _c(cmp)
{
}


template<class R>
Geometry::LinearConstraint<R>* 
Geometry::LinearConstraint<R>::clone() const
{
  return new LinearConstraint<R>(this->_a, this->_c, this->_b);
}


template<class R>
dimension_type
Geometry::LinearConstraint<R>::dimension() const 
{
  return this->_a.size();
}


template<class R>
smoothness_type 
Geometry::LinearConstraint<R>::smoothness() const 
{
  return std::numeric_limits<smoothness_type>::max();
}


template<class R>
Geometry::Comparison
Geometry::LinearConstraint<R>::comparison() const 
{
  return this->_c;
}


template<class R>
Geometry::Polyhedron<R>
Geometry::LinearConstraint<R>::polyhedron() const 
{
  if(this->_c==less) {
    return Polyhedron<R>(LinearAlgebra::Matrix<R>(1,this->dimension(),_a.begin()),
                         LinearAlgebra::Vector<R>(1,&_b));
  } else {
    return Polyhedron<R>(-LinearAlgebra::Matrix<R>(1,this->dimension(),_a.begin()),
                         -LinearAlgebra::Vector<R>(1,&_b));
  }
}





template<class R>
std::ostream& 
Geometry::LinearConstraint<R>::write(std::ostream& os) const
{
  //return os << "LinearConstraint( a=" << this->_a << ", b=" << this->_b << ", c='" << (this->_c==less ? "<" : ">") << "' )";
  return os << this->_a << ".x"<< (this->_c==less ? "<" : ">") << "=" << this->_b;
}

      
template<class R>
typename Geometry::LinearConstraint<R>::A
Geometry::LinearConstraint<R>::value(const Point<A>& pt) const
{
  return LinearAlgebra::inner_product(pt.position_vector(),LinearAlgebra::Vector<A>(this->_a)) - this->_b;
}


template<class R>
LinearAlgebra::Vector<typename Geometry::LinearConstraint<R>::A>
Geometry::LinearConstraint<R>::gradient(const Point<A>& pt) const
{
  return LinearAlgebra::Vector<A>(this->_a);
}



} // namespace Ariadne
