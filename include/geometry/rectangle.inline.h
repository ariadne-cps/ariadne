/***************************************************************************
 *            rectangle.inline.h
 *
 *  Copyright 2005-6  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 


#include "../throw.h"
#include "../base/array.h"
#include "../base/iterator.h"
#include "../base/tribool.h"

#include "../numeric/arithmetic.h"
#include "../numeric/function.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"

#include "../geometry/exceptions.h"
#include "../geometry/point.h"
#include "../geometry/rectangle_expression.h"

namespace Ariadne {
  namespace Geometry {
  
    template<class R> inline
    Rectangle<R>::Rectangle(size_type d)
      : _data(2*d)
    { 
      if(d!=0) { this->_data[0]=1; this->_data[1]=0; }
    }
    
    template<class R> template<class ForwardIterator> inline
    Rectangle<R>::Rectangle(ForwardIterator b, ForwardIterator e)
      : _data(2*std::distance(b,e))
    {
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        this->set_lower_bound(i,b->lower());
        this->set_upper_bound(i,b->upper());
        ++b;
      }
    }
    
    template<class R> inline
    Rectangle<R>::Rectangle(const array< Numeric::Interval<R> >& a)
      : _data(2*a.size())
    {
      for(dimension_type i=0; i!=a.size(); ++i) {
        this->set_lower_bound(i,a[i].lower());
        this->set_upper_bound(i,a[i].upper());
      }
    }
    
    template<class R> inline
    Rectangle<R>::Rectangle(const std::vector< Numeric::Interval<R> >& v)
      : _data(2*v.size())
    {
      for(dimension_type i=0; i!=v.size(); ++i) {
        this->set_lower_bound(i,v[i].lower());
        this->set_upper_bound(i,v[i].upper());
      }
    }

    template<class R> inline
    Rectangle<R>::Rectangle(const Point<R>& pt)
      : _data(2*pt.dimension())
    {
      for(dimension_type i=0; i!=pt.dimension(); ++i) {
        this->set_lower_bound(i,pt[i]);
        this->set_upper_bound(i,pt[i]);
      }
    }
    
    template<class R> inline
    Rectangle<R>::Rectangle(const Point< Numeric::Interval<R> >& pt)
      : _data(2*pt.dimension())
    {
      for(dimension_type i=0; i!=pt.dimension(); ++i) {
        this->set_lower_bound(i,pt[i].lower());
        this->set_upper_bound(i,pt[i].upper());
      }
    }
    
    template<class R> inline
    Rectangle<R>::Rectangle(const R& lower_bound, const R& upper_bound)
      : _data(2)
    {
      this->set_lower_bound(0,lower_bound);
      this->set_upper_bound(0,upper_bound);
    }

    template<class R> inline
    Rectangle<R>::Rectangle(const Point<R>& pt1, const Point<R>& pt2) 
      : _data(2*pt1.dimension())
    {
      if(pt1.dimension()!=pt2.dimension()) {
        ARIADNE_THROW(IncompatibleDimensions,"Rectangle(Point pt1, Point pt2)","pt1=" << pt1 << ", pt2=" << pt2);
      }
      for (size_type i=0; i!=this->dimension(); ++i) {
        this->set_lower_bound(i,Numeric::min_exact(pt1[i],pt2[i]));
        this->set_upper_bound(i,Numeric::max_exact(pt1[i],pt2[i]));
      }
    }
    
    
    template<class R> inline
    Rectangle<R>::Rectangle(const LinearAlgebra::Vector< Numeric::Interval<R> >& iv)
      : _data(2*iv.size())
    {
      for (size_type i=0; i!=this->dimension(); ++i) {
        this->set_lower_bound(i,iv(i).lower());
        this->set_upper_bound(i,iv(i).upper());
      }
    }

    template<class R> template<class E> inline
    Rectangle<R>::Rectangle(const RectangleExpression<E>& original)
      : _data(2*original().dimension())
    {         
      const E& expression=original();
      for (size_type i=0; i!=this->dimension(); ++i) {
        this->_data[2*i]=expression.lower_bound(i);
        this->_data[2*i+1]=expression.upper_bound(i);
      }
    }
  
    template<class R> inline
    Rectangle<R>::Rectangle(const Rectangle<R>& original)
      : _data(original._data)
    { }
  
    template<class R> inline
    Rectangle<R>& Rectangle<R>::operator=(const Rectangle<R>& A) {
      if(this != &A) {
        this->_data = A._data;
      }
      return *this;
    }

    template<class R> template<class E> inline
    Rectangle<R>& Rectangle<R>::operator=(const RectangleExpression<E>& original)
    {         
      const E& expression=original();
      this->_data.resize(2*expression.dimension());
      for (size_type i=0; i!=this->dimension(); ++i) {
        this->set_lower_bound(i,expression.lower_bound(i));
        this->set_upper_bound(i,expression.upper_bound(i));
      }
      return *this;
    }
    
    
    // Conversion operators
    template<class R> inline
    Rectangle<R>::operator Point< Numeric::Interval<R> >() const 
    {
      return Point< Numeric::Interval<R> >(this->dimension(),reinterpret_cast<const Numeric::Interval<R>*>(this->_data.begin()));
    }    
    
    
    // Comparison operators
    template<class R> inline
    bool Rectangle<R>::operator==(const Rectangle<R>& A) const
    {
      if (A.empty() && this->empty()) { return true; }
      if (A.empty() || this->empty()) { return false; }
      if(this->dimension()!=A.dimension()) { return false; }
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        if (this->lower_bound(i)!=A.lower_bound(i)) { return false; }
        if (this->upper_bound(i)!=A.upper_bound(i)) { return false; }
      }
      return true;
    }
    
    template<class R> inline
    bool Rectangle<R>::operator!=(const Rectangle<R>& A) const 
    {
      return !(*this == A);
    }


    // Data access
    template<class R> inline
    array<R>& Rectangle<R>::data()
    {
      return this->_data;
    }
    
    template<class R> inline
    const array<R>& Rectangle<R>::data() const
    {
      return this->_data;
    }
    
    
   
    template<class R> inline
    const R& Rectangle<R>::lower_bound(dimension_type i) const 
    {
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle::lower_bound(dimension_type i) const","*this=" << *this << ", i=" << i);
      }
      return this->_data[2*i];
    }
    
    template<class R> inline
    R& Rectangle<R>::lower_bound(dimension_type i) 
    {
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle::lower_bound(dimension_type i)","self=" << *this << ", i=" << i);
      }
      return this->_data[2*i];
    }
    
    template<class R> inline
    const R& Rectangle<R>::upper_bound(dimension_type i) const 
    {
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle::upper_bound(dimension_type i) const","self=" << *this << ", i=" << i);
      }
      return this->_data[2*i+1];
    }
    
    template<class R> inline
    R& Rectangle<R>::upper_bound(dimension_type i) 
    {
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle::upper_bound(dimension_type i)","self=" << *this << ", i=" << i);
      }
      return this->_data[2*i+1];
    }
    
    
    template<class R> inline
    Numeric::Interval<R>& Rectangle<R>::operator[] (dimension_type i) 
    {
      return reinterpret_cast<Numeric::Interval<R>&>(this->_data[2*i]);
    }
 
    template<class R> inline
    const Numeric::Interval<R>& Rectangle<R>::operator[] (dimension_type i) const 
    {
      return reinterpret_cast<const Numeric::Interval<R>&>(this->_data[2*i]);
    }
    
    
    template<class R> inline
    const Numeric::Interval<R>& Rectangle<R>::interval(dimension_type i) const 
    {
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle::interval(dimension_type i) const","self=" << *this << ", i=" << i);
      }
      return reinterpret_cast<const Numeric::Interval<R>&>(this->_data[2*i]);
    }
    
    
    template<class R> inline
    Point<R> Rectangle<R>::lower_corner() const 
    {
      Point<R> result(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result[i]=this->lower_bound(i);
      }
      return result;
    }
    
    template<class R> inline
    Point<R> Rectangle<R>::upper_corner() const 
    {
      Point<R> result(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result[i]=this->upper_bound(i);
      }
      return result;
    }
    
    
    template<class R> inline
    LinearAlgebra::Vector< Numeric::Interval<R> > Rectangle<R>::position_vectors() const 
    {
      LinearAlgebra::Vector< Numeric::Interval<R> > result(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result(i)=this->interval(i);
      }
      return result;
    }


    // Modifying operations
    template<class R> inline
    void Rectangle<R>::clear()
    {
      if(this->_data.size()!=0) {
        this->_data[0]=1;
        this->_data[1]=0;
      }
    }
    
    template<class R> inline
    void Rectangle<R>::set_interval(dimension_type i, Numeric::Interval<R> x)
    {
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle::set_interval(dimension_type i, Interval x) const","self=" << *this << ", i=" << i);
      }
      this->set_lower_bound(i,x.lower());
      this->set_upper_bound(i,x.upper());
    }
    
    template<class R> inline
    void Rectangle<R>::set_lower_bound(dimension_type i, const R& l) 
    {
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle::set_lower_bound(dimension_type i, Real l) const","self=" << *this << ", i=" << i);
      }
      this->_data[2*i]=l;
    }
    
    template<class R> inline
    void Rectangle<R>::set_upper_bound(dimension_type i, const R& u) 
    {
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle::set_upper_bound(dimension_type i, Real u) const","self=" << *this << ", i=" << i);
      }
      this->_data[2*i+1]=u;
    }

    template<class R> inline
    Rectangle<R>& Rectangle<R>::expand_by(const real_type& delta) 
    {
      Numeric::Interval<R> expand(-delta,delta);
      for (size_type j=0; j< this->dimension(); ++j) {
        (*this)[j]+=expand;
      }
        
      return *this;
    }
      
    template<class R> inline
    Rectangle<R> Rectangle<R>::expand(const real_type& delta) const
    {
      Rectangle<R> result(*this);
      result.expand_by(delta);
      return result;
    }
      
    template<class R> inline
    Rectangle<R> Rectangle<R>::neighbourhood(const real_type& delta) const
    {
      Rectangle<R> result(this->dimension());
      for (size_type j=0; j!=this->dimension(); ++j) {
        result.set_lower_bound(j,sub_approx(this->lower_bound(j),delta));
        result.set_upper_bound(j,add_approx(this->upper_bound(j),delta));
      }
      return result;
    }
      


    // Rectangle geometric operations
    template<class R> inline
    size_type Rectangle<R>::dimension() const 
    {
      return this->_data.size()/2;
    }
    
    template<class R> inline
    tribool Rectangle<R>::empty() const 
    {
      tribool result=false;
      if(this->dimension()==0) {
        return true;
      }
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        if(this->lower_bound(i) > this->upper_bound(i)) {
          return true;
        }
        if(this->lower_bound(i)== this->upper_bound(i)) {
          result=indeterminate;
        }
      }
      return result;
    }
    
    template<class R> inline
    Point<R> Rectangle<R>::centre() const
    {
      Point<R> result(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result[i]=Numeric::div_approx(Numeric::add_approx(this->lower_bound(i),this->upper_bound(i)),R(2));
      }
      return result;
    }
    
    template<class R> inline
    R Rectangle<R>::radius() const 
    {
      R diameter=0;
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        diameter=Numeric::max_up(diameter,Numeric::sub_up(this->upper_bound(i),this->lower_bound(i)));
      }
      return div_up(diameter,R(2));
    }
    
    template<class R> inline
    R Rectangle<R>::volume() const 
    {
      R result=1;
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result=mul_approx(result,sub_approx(this->upper_bound(i),this->lower_bound(i)));
      }
      return result;
    }
    

    template<class R>
    inline
    tribool 
    Rectangle<R>::contains(const Point<R>& pt) const 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,pt,"Rectangle::contains(Point pt)");
      tribool result=true;
      const Rectangle<R>& self=*this;
      for (size_type i=0; i!=self.dimension(); ++i) {
        if(self.lower_bound(i)>pt[i] || pt[i]>self.upper_bound(i)) {
          return false;
        }
        if(self.lower_bound(i)==pt[i] || pt[i]==self.upper_bound(i)) { 
          result=indeterminate;
        }
      }
      return result;
    }


    template<class R> inline
    tribool Rectangle<R>::bounded() const 
    { 
      return true; 
    }
    
    template<class R> inline
    Rectangle<R> Rectangle<R>::bounding_box() const 
    {
      return *this;
    }

    
    
    
    
    template<class R> inline
    Rectangle< Numeric::Interval<R> >::Rectangle(dimension_type d)
      : _data(2*d) { }
    
    template<class R> inline
    Rectangle< Numeric::Interval<R> >::Rectangle(const I& x1, const I& x2)
      : _data(2)
    {
      if(x1.lower()<=x2.lower() && x1.upper()<=x2.upper()) {
        this->set_lower_bound(0,x1);
        this->set_upper_bound(0,x2);
      } else if(x1.lower()>=x2.lower() && x1.upper()>=x2.upper()) {
        this->set_lower_bound(0,x2);
        this->set_upper_bound(0,x1);
      } else {
        ARIADNE_THROW(InvalidValue,"Rectangle<Interval>(Interval x1, Interval x2)","x1=" << x1 << ", x2=" << x2);
      }
    }


    template<class R> inline
    Rectangle< Numeric::Interval<R> >::Rectangle(const Point<I>& pt1, const Point<I>& pt2)
      : _data(2*pt1.dimension())
    {
      if(pt1.dimension()!=pt2.dimension()) {
        ARIADNE_THROW(IncompatibleDimensions,"Rectangle<Interval>(Point pt1, Point pt2)","pt1=" << pt1 << ", pt2=" << pt2);
      }
      for (size_type i=0; i!=this->dimension(); ++i) {
        if(pt1[i].lower()<=pt2[i].lower() && pt1[i].upper()<=pt2[i].upper()) {
          this->set_lower_bound(i,pt1[i]);
          this->set_upper_bound(i,pt2[i]);
        } else if(pt1[i].lower()>=pt2[i].lower() && pt1[i].upper()>=pt2[i].upper()) {
          this->set_lower_bound(0,pt2[i]);
          this->set_upper_bound(0,pt1[i]);
        } else {
          ARIADNE_THROW(InvalidValue,"Rectangle<Interval>(Point pt1, Point pt2)","pt1=" << pt1 << ", pt2=" << pt2);
        }
      }
    }

    template<class R> template<class E> inline
    Rectangle< Numeric::Interval<R> >::Rectangle(const RectangleExpression<E>& e)
      : _data(2*e().dimension()) 
    { 
      this->assign(e()); 
    }
    
    template<class R> template<class E> inline
    Rectangle< Numeric::Interval<R> >& 
    Rectangle< Numeric::Interval<R> >::operator=(const RectangleExpression<E>& e)
    {
      this->_data.resize(e().dimension()); 
      this->assign(e); 
    }
    
    template<class R> inline
    dimension_type Rectangle< Numeric::Interval<R> >::dimension() const { 
      return this->_data.size()/2; 
    }
    
    template<class R> inline
    const Numeric::Interval<R>& 
    Rectangle< Numeric::Interval<R> >::lower_bound(const dimension_type& i) const 
    { 
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle<Interval>::lower_bound(dimension_type i) const","self=" << *this << ", i=" << i);
      }
      return _data[2*i]; 
    }
    
    template<class R> inline
    const Numeric::Interval<R>& 
    Rectangle< Numeric::Interval<R> >::upper_bound(const dimension_type& i) const
    { 
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle<Interval>::upper_bound(dimension_type i) const","self=" << *this << ", i=" << i);
      }
      return _data[2*i+1];
    }
    
    template<class R> inline
    void 
    Rectangle< Numeric::Interval<R> >::set_lower_bound(const dimension_type& i, const Numeric::Interval<R>& x)
    { 
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle<Interval>::set_lower_bound(dimension_type i, Interval x) const",
                      "self=" << *this << ", i=" << i << ", x=" << x);
      }
      _data[2*i]=x; 
    }
    
    template<class R> inline
    void 
    Rectangle< Numeric::Interval<R> >::set_upper_bound(const dimension_type& i, const Numeric::Interval<R>& x)
    { 
      if(i>=this->dimension()) {
        ARIADNE_THROW(InvalidCoordinate,"Rectangle<Interval>::set_upper_bound(dimension_type i, Interval x) const",
                      "self=" << *this << ", i=" << i << ", x=" << x);
      }
      _data[2*i+1]=x; 
    }
    
    template<class R> inline
    Rectangle<R> 
    Rectangle< Numeric::Interval<R> >::bounding_box() const
    { 
      Rectangle<R> result(this->dimension());
      for(dimension_type i=0; i!=result.dimension(); ++i) {
        result.set_lower_bound(i,this->lower_bound(i).lower());
        result.set_upper_bound(i,this->upper_bound(i).upper());
      }
      return result;
    }
    
    template<class R> template<class RE> inline
    void 
    Rectangle< Numeric::Interval<R> >::assign(const RE& re) 
    { 
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        this->_data[2*i]=re.lower_bound(i); this->_data[2*i+1]=re.upper_bound(i);
      }
    }

  
    
  
    template<class R> inline
    Rectangle<R> over_approximation(const Rectangle<R>& r) 
    {
      return r;
    }


    template<class R> inline
    Rectangle<R> under_approximation(const Rectangle<R>& r) 
    {
      return r;
    }


    template<class R> inline
    Rectangle<R> over_approximation(const Rectangle< Numeric::Interval<R> >& ir) 
    {
      Rectangle<R> result(ir.dimension());
      for(dimension_type i=0; i!=result.dimension(); ++i) {
        result.set_lower_bound(i,ir.lower_bound(i).lower());
        result.set_upper_bound(i,ir.upper_bound(i).upper());
      }
      return result;
    }
    
    template<class R> inline
    Rectangle<R> under_approximation(const Rectangle< Numeric::Interval<R> >& ir) 
    {
      Rectangle<R> result(ir.dimension());
      for(dimension_type i=0; i!=result.dimension(); ++i) {
        result.set_lower_bound(i,ir.lower_bound(i).upper());
        result.set_upper_bound(i,ir.upper_bound(i).lower());
      }
      return result;
    }
    
      
    template<class R> inline
    tribool 
    equal(const Rectangle<R>& r1, const Rectangle<R>& r2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"equal(Rectangle r1, Rectangle r2)")
      for(size_type i=0; i!=r1.dimension(); ++i) {
        if(r1.lower_bound(i)!=r2.lower_bound(i) || r1.upper_bound(i)!=r2.upper_bound(i)) {
          return false;
        }
      }
      return indeterminate;
    }
      
    template<class R> inline
    tribool 
    disjoint(const Rectangle<R>& r1, const Rectangle<R>& r2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"disjoint(Rectangle r1, Rectangle r2)")
      tribool result=false;
      for(size_type i=0; i!=r1.dimension(); ++i) {
        if(r1.lower_bound(i)>r2.upper_bound(i) || r1.upper_bound(i)<r2.lower_bound(i)) {
          return true;
        }
        if(r1.lower_bound(i)==r2.upper_bound(i) || r1.upper_bound(i)==r2.lower_bound(i)) {
          result=indeterminate;
        }
      }
      return result;
    }
       
  
    template<class R> inline
    tribool 
    subset(const Rectangle<R>& r1, const Rectangle<R>& r2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"subset(Rectangle r1, Rectangle r2)")
      tribool result=true;
      for (size_type i=0; i!=r1.dimension(); ++i) {
        if(r1.lower_bound(i)<r2.lower_bound(i) || r1.upper_bound(i)>r2.upper_bound(i)) {
          return false;
        }
        if(r1.lower_bound(i)==r2.lower_bound(i) || r1.upper_bound(i)==r2.upper_bound(i)) {
          result=indeterminate;
        }
      }
      return result;
    }
    
    
    template<class R> inline
    Rectangle<R> 
    closed_intersection(const Rectangle<R>& r1, const Rectangle<R>& r2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"Rectangle closed_intersection(Rectangle r1, Rectangle r2)");
      Rectangle<R> r3(r1.dimension());
      for(size_type i=0; i != r3.dimension(); ++i) {
        r3[i]=Numeric::intersection(r1[i],r2[i]);
      }
      return r3;
    }
  
    template<class R> inline
    Rectangle<R> 
    open_intersection(const Rectangle<R>& r1, const Rectangle<R>& r2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"Rectangle closed_intersection(Rectangle r1, Rectangle r2)");
      Rectangle<R> r3(r1.dimension());
      for(size_type i=0; i != r3.dimension(); ++i) {
        r3[i]=Numeric::intersection(r1[i],r2[i]);
        if(r3[i].lower()>=r3[i].upper()) {
          r3[i]=Numeric::Interval<R>();
        }
      }
      return r3;
    }
  
    template<class R> inline
    Rectangle<R>
    rectangular_hull(const Rectangle<R>& r1, const Rectangle<R>& r2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"Rectangle rectangular_hull(Rectangle r1, Rectangle r2)");
      Rectangle<R> r3(r1.dimension());
      for(size_type i=0; i != r3.dimension(); ++i) {
        r3[i]=Numeric::hull(r1[i],r2[i]);
      }
      return r3;
    }
  
    template<class R> inline
    Rectangle<R>
    rectangular_hull(const Rectangle<R>& r, const Point<R>& pt)
    {
      return rectangular_hull(r,Rectangle<R>(pt));
    }
  
    template<class R> inline
    Rectangle<R>
    rectangular_hull(const Point<R>& pt, const Rectangle<R>& r)
    {
      return rectangular_hull(Rectangle<R>(pt),r);
    }

    template<class R> inline
    Rectangle<R>
    rectangular_hull(const Point<R>& pt1, const Point<R>& pt2)
    {
      return Rectangle<R>(pt1,pt2);
    }

 
    template<class R1, class R2> inline
    Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> 
    minkowski_sum(const Rectangle<R1>& r1, const Rectangle<R2>& r2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"Rectangle minkowski_sum(Rectangle r1, Rectangle r2)");
      Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> r3(r1.dimension());
      for(dimension_type i=0; i!=r3.dimension(); ++i) {
        r3.set_lower_bound(i,r1.lower_bound(i)+r2.lower_bound(i));
        r3.set_lower_bound(i,r1.upper_bound(i)+r2.upper_bound(i));
      }
      return r3;
    }
  
    template<class R1, class R2> inline
    Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> 
    minkowski_difference(const Rectangle<R1>& r1, const Rectangle<R2>& r2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"Rectangle minkowski_difference(Rectangle r1, Rectangle r2)");
      Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> r3(r1.dimension());
      for(dimension_type i=0; i!=r3.dimension(); ++i) {
        r3.set_lower_bound(i,r1.lower_bound(i)-r2.lower_bound(i));
        r3.set_upper_bound(i,r1.upper_bound(i)-r2.upper_bound(i));
      }
      return r3;
    }
    
    
    template<class R> inline
    tribool 
    subset(const Rectangle<R>& r, ListSet< Geometry::Rectangle<R> >& ls);
        
  
    
    
    
    template<class R> inline
    LinearAlgebra::Vector< Numeric::Interval<R> > 
    operator-(const Geometry::Rectangle<R>& r1,
              const Geometry::Rectangle<R>& r2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r1,r2,"Vector operator-(Rectangle r1, Rectangle r2)");
      return r1.position_vectors()-r2.position_vectors();
    }
  
    template<class R, class E> inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::VectorExpression<E>& v)
    {
      const E& ve=v();
      ARIADNE_CHECK_DIMENSION_EQUALS_SIZE(r,ve,"Vector operator-(Rectangle r, VectorExpression ve)");
      LinearAlgebra::Vector< Numeric::Interval<R> > iv=ve;
      return r+iv; 
    }
      
    template<class R>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      ARIADNE_CHECK_DIMENSION_EQUALS_SIZE(r,v,"Vector operator+(Rectangle r, Vector v)");
       
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+v(i));
      }
      return result;
    }
  
    template<class R>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector< Numeric::Interval<R> >& iv)
    {
      Geometry::Rectangle<R> result(r.dimension());
      ARIADNE_CHECK_DIMENSION_EQUALS_SIZE(r,iv,"Vector operator+(Rectangle r, Vector<Interval> iv)");
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+iv(i));
      }
      return result;
    }
  
    
    template<class R> inline
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      ARIADNE_CHECK_DIMENSION_EQUALS_SIZE(r,v,"Vector operator-(Rectangle r, Vector v)");
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]-v(i));
      }
      return result;
    }
  
    template<class R> inline
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector< Numeric::Interval<R> >& iv)
    {
      Geometry::Rectangle<R> result(r.dimension());
      ARIADNE_CHECK_DIMENSION_EQUALS_SIZE(r,iv,"Vector operator-(Rectangle r, Vector<Interval> v)");
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]-iv(i));
      }
      return result;
    }
  
    template<class R> inline
    Geometry::Rectangle<R> 
    scale(const Geometry::Rectangle<R>& r, const R& scale_factor) 
    {
      Geometry::Rectangle<R> result(r.dimension());
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,scale_factor*r[i]);
      }
      return result;
    }
    
    


    template<class R> inline
    RectangleVerticesIterator<R>::RectangleVerticesIterator(const Rectangle<R>& r, const bool end)
      : _r(&r), _i(end==true ? (1<<(r.dimension()-1))*3 : 0), _parity(0), _pt(r.lower_corner()) { }
      
    template<class R> inline
    bool RectangleVerticesIterator<R>::equal(const RectangleVerticesIterator<R>& other) const {
      return this->_i==other._i && this->_r==other._r; }
      
    template<class R> inline
    const Point<R>& RectangleVerticesIterator<R>::dereference() const { 
      return this->_pt; }
      
    template<class R> inline
    void RectangleVerticesIterator<R>::increment() { 
      uint j=0; uint m=1; if(this->_parity) { while(!(m&(this->_i))) { ++j; m*=2u; } ++j; m*=2u; }
      this->_parity=!this->_parity;
      if(j==this->_r->dimension()) { this->_i+=m; return; }
      if(m&(this->_i)) { this->_pt[j]=this->_r->lower_bound(j); this->_i-=m; }
      else { this->_pt[j]=this->_r->upper_bound(j); this->_i+=m; }
    }
    
    template<class R> inline 
    std::ostream& operator<<(std::ostream& os, const Rectangle<R>& r) {
      return r.write(os);
    }
    
    template<class R> inline
    std::istream& operator>>(std::istream& is, Rectangle<R>& r) {
      return r.read(is);
    }
  
  }
}
