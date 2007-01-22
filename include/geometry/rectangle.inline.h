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
 


#include "../declarations.h"

#include "../base/array.h"
#include "../base/iterator.h"
#include "../base/tribool.h"
#include "../exceptions.h"

#include "../numeric/arithmetic.h"
#include "../numeric/function.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"

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
    Rectangle<R>::Rectangle(const array< Interval<R> >& a)
      : _data(2*a.size())
    {
      for(dimension_type i=0; i!=a.size(); ++i) {
        this->set_lower_bound(i,a[i].lower());
        this->set_upper_bound(i,a[i].upper());
      }
    }
    
    template<class R> inline
    Rectangle<R>::Rectangle(const std::vector< Interval<R> >& v)
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
    Rectangle<R>::Rectangle(const Point< Interval<R> >& pt)
      : _data(2*pt.dimension())
    {
      for(dimension_type i=0; i!=pt.dimension(); ++i) {
        this->set_lower_bound(i,pt[i].lower());
        this->set_upper_bound(i,pt[i].upper());
      }
    }
    
    template<class R> inline
    Rectangle<R>::Rectangle(const Point<R>& pt1, const Point<R>& pt2) 
      : _data(2*pt1.dimension())
    {
      check_equal_dimensions(pt1,pt2,__PRETTY_FUNCTION__);
      for (size_type i=0; i!=this->dimension(); ++i) {
        this->set_lower_bound(i,Numeric::min_exact(pt1[i],pt2[i]));
        this->set_upper_bound(i,Numeric::max_exact(pt1[i],pt2[i]));
      }
    }
    
    
    template<class R> inline
    Rectangle<R>::Rectangle(const LinearAlgebra::Vector< Interval<R> >& iv)
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
    Rectangle<R>::operator Point< Interval<R> >() const 
    {
      return Point< Interval<R> >(this->dimension(),reinterpret_cast<const Interval<R>*>(this->_data.begin()));
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
      check_coordinate(*this,i,__PRETTY_FUNCTION__);
      return this->_data[2*i];
    }
    
    template<class R> inline
    R& Rectangle<R>::lower_bound(dimension_type i) 
    {
      check_coordinate(*this,i,__PRETTY_FUNCTION__);
      return this->_data[2*i];
    }
    
    template<class R> inline
    const R& Rectangle<R>::upper_bound(dimension_type i) const 
    {
      check_coordinate(*this,i,__PRETTY_FUNCTION__);
      return this->_data[2*i+1];
    }
    
    template<class R> inline
    R& Rectangle<R>::upper_bound(dimension_type i) 
    {
      check_coordinate(*this,i,__PRETTY_FUNCTION__);
      return this->_data[2*i+1];
    }
    
    
    template<class R> inline
    Interval<R>& Rectangle<R>::operator[] (dimension_type i) 
    {
      check_coordinate(*this,i,__PRETTY_FUNCTION__);
      return reinterpret_cast<Interval<R>&>(this->_data[2*i]);
    }
 
    template<class R> inline
    const Interval<R>& Rectangle<R>::operator[] (dimension_type i) const 
    {
      check_coordinate(*this,i,__PRETTY_FUNCTION__);
      return reinterpret_cast<const Interval<R>&>(this->_data[2*i]);
    }
    
    
    template<class R> inline
    const Interval<R>& Rectangle<R>::interval(dimension_type i) const 
    {
      check_coordinate(*this,i,__PRETTY_FUNCTION__);
      return reinterpret_cast<const Interval<R>&>(this->_data[2*i]);
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
    LinearAlgebra::Vector< Interval<R> > Rectangle<R>::position_vectors() const 
    {
      LinearAlgebra::Vector< Interval<R> > result(this->dimension());
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
    void Rectangle<R>::set_interval(dimension_type i, Interval<R> x)
    {
      check_coordinate(*this,i,__PRETTY_FUNCTION__);
      this->set_lower_bound(i,x.lower());
      this->set_upper_bound(i,x.upper());
    }
    
    template<class R> inline
    void Rectangle<R>::set_lower_bound(dimension_type i, const R& l) 
    {
      check_coordinate(*this,i,__PRETTY_FUNCTION__);
      this->_data[2*i]=l;
    }
    
    template<class R> inline
    void Rectangle<R>::set_upper_bound(dimension_type i, const R& u) 
    {
      check_coordinate(*this,i,__PRETTY_FUNCTION__);
      this->_data[2*i+1]=u;
    }

    template<class R> inline
    Rectangle<R>& Rectangle<R>::expand_by(const real_type& delta) 
    {
      Interval<R> expand(-delta,delta);
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
    Rectangle<R>::contains(const Point<R>& p) const 
    {
      tribool result=true;
      check_equal_dimensions(*this,p,__PRETTY_FUNCTION__);
      const Rectangle<R>& self=*this;
      for (size_type i=0; i!=self.dimension(); ++i) {
        if(self.lower_bound(i)>p[i] || p[i]>self.upper_bound(i)) {
          return false;
        }
        if(self.lower_bound(i)==p[i] || p[i]==self.upper_bound(i)) { 
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
    Rectangle< Interval<R> >::Rectangle(dimension_type d)
      : _data(2*d) { }
    
    template<class R> template<class E> inline
    Rectangle< Interval<R> >::Rectangle(const RectangleExpression<E>& e)
      : _data(2*e().dimension()) 
    { 
      this->assign(e()); 
    }
    
    template<class R> template<class E> inline
    Rectangle< Interval<R> >& 
    Rectangle< Interval<R> >::operator=(const RectangleExpression<E>& e)
    {
      this->_data.resize(e().dimension()); 
      this->assign(e); 
    }
    
    template<class R> inline
    dimension_type Rectangle< Interval<R> >::dimension() const { 
      return this->_data.size()/2; 
    }
    
    template<class R> inline
    const Interval<R>& 
    Rectangle< Interval<R> >::lower_bound(const dimension_type& i) const 
    { 
      return _data[2*i]; 
    }
    
    template<class R> inline
    const Interval<R>& 
    Rectangle< Interval<R> >::upper_bound(const dimension_type& i) const
    { 
      return _data[2*i+1];
    }
    
    template<class R> inline
    void 
    Rectangle< Interval<R> >::set_lower_bound(const dimension_type& i, const Interval<R>& x)
    { 
      _data[2*i]=x; 
    }
    
    template<class R> inline
    void 
    Rectangle< Interval<R> >::set_upper_bound(const dimension_type& i, const Interval<R>& x)
    { 
      _data[2*i+1]=x; 
    }
    
    template<class R> template<class RE> inline
    void 
    Rectangle< Interval<R> >::assign(const RE& re) 
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
    Rectangle<R> over_approximation(const Rectangle< Interval<R> >& ir) 
    {
      Rectangle<R> result(ir.dimension());
      for(dimension_type i=0; i!=result.dimension(); ++i) {
        result.set_lower_bound(i,ir.lower_bound(i).lower());
        result.set_upper_bound(i,ir.upper_bound(i).upper());
      }
      return result;
    }
    
    template<class R> inline
    Rectangle<R> under_approximation(const Rectangle< Interval<R> >& ir) 
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
    equal(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(size_type i=0; i!=A.dimension(); ++i) {
        if(A.lower_bound(i)!=B.lower_bound(i) || A.upper_bound(i)!=B.upper_bound(i)) {
          return false;
        }
      }
      return indeterminate;
    }
      
    template<class R> inline
    tribool 
    disjoint(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      tribool result=false;
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(size_type i=0; i!=A.dimension(); ++i) {
        if(A.lower_bound(i)>B.upper_bound(i) || A.upper_bound(i)<B.lower_bound(i)) {
          return true;
        }
        if(A.lower_bound(i)==B.upper_bound(i) || A.upper_bound(i)==B.lower_bound(i)) {
          result=indeterminate;
        }
      }
      return result;
    }
       
  
    template<class R> inline
    tribool 
    subset(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      tribool result=true;
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for (size_type i=0; i!=A.dimension(); ++i) {
        if(A.lower_bound(i)<B.lower_bound(i) || A.upper_bound(i)>B.upper_bound(i)) {
          return false;
        }
        if(A.lower_bound(i)==B.lower_bound(i) || A.upper_bound(i)==B.upper_bound(i)) {
          result=indeterminate;
        }
      }
      return result;
    }
    
    
    template<class R> inline
    Rectangle<R> 
    closed_intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Numeric::intersection(A[i],B[i]);
      }
      return C;
    }
  
    template<class R> inline
    Rectangle<R> 
    open_intersection(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Numeric::intersection(A[i],B[i]);
        if(C[i].lower()>=C[i].upper()) {
          C[i]=Interval<R>();
        }
      }
      return C;
    }
  
    template<class R> inline
    Rectangle<R>
    rectangular_hull(const Rectangle<R>& A, const Rectangle<R>& B)
    {
      Rectangle<R> C(A.dimension());
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(size_type i=0; i != C.dimension(); ++i) {
        C[i]=Numeric::hull(A[i],B[i]);
      }
      return C;
    }
  
    template<class R> inline
    Rectangle<R>
    rectangular_hull(const Rectangle<R>& A, const Point<R>& B)
    {
      return rectangular_hull(A,Rectangle<R>(B));
    }
  
    template<class R> inline
    Rectangle<R>
    rectangular_hull(const Point<R>& A, const Rectangle<R>& B)
    {
      return rectangular_hull(B,A);
    }

    template<class R> inline
    Rectangle<R>
    rectangular_hull(const Point<R>& A, const Point<R>& B)
    {
      return Rectangle<R>(A,B);
    }

 
    template<class R1, class R2> inline
    Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> 
    minkowski_sum(const Rectangle<R1>& A, const Rectangle<R2>& B)
    {
      Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> C(A.dimension());
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(dimension_type i=0; i!=C.dimension(); ++i) {
        C.set_lower_bound(i,A.lower_bound(i)+B.lower_bound(i));
        C.set_lower_bound(i,A.upper_bound(i)+B.upper_bound(i));
      }
      return C;
    }
  
    template<class R1, class R2> inline
    Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> 
    minkowski_difference(const Rectangle<R1>& A, const Rectangle<R2>& B)
    {
      Rectangle<typename Numeric::traits<R1,R2>::arithmetic_type> C(A.dimension());
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      for(dimension_type i=0; i!=C.dimension(); ++i) {
        C.set_lower_bound(i,A.lower_bound(i)-B.lower_bound(i));
        C.set_upper_bound(i,A.upper_bound(i)-B.upper_bound(i));
      }
      return C;
    }
    
    
    template<class R> inline
    tribool 
    subset(const Rectangle<R>& A, ListSet<R,Geometry::Rectangle>& B);
        
  
    
    
    
    template<class R> inline
    LinearAlgebra::Vector< Interval<R> > 
    operator-(const Geometry::Rectangle<R>& r1,
              const Geometry::Rectangle<R>& r2)
    {
      check_equal_dimensions(r1,r2,__PRETTY_FUNCTION__);
      return r1.position_vectors()-r2.position_vectors();
    }
  
    template<class R, class E> inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::VectorExpression<E>& v)
    {
      const E& ev=v();
      check_dimension(r,ev.size(),__PRETTY_FUNCTION__);
      LinearAlgebra::Vector< Interval<R> > iv=ev;
      return r+iv; 
    }
      
    template<class R>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      check_dimension(r,v.size(),__PRETTY_FUNCTION__);
       
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+v(i));
      }
      return result;
    }
  
    template<class R>
    inline
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector< Interval<R> >& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      check_dimension(r,v.size(),__PRETTY_FUNCTION__);
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+v(i));
      }
      return result;
    }
  
    
    template<class R> inline
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector<R>& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      check_dimension(r,v.size(),__PRETTY_FUNCTION__);
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]-v(i));
      }
      return result;
    }
  
    template<class R> inline
    Geometry::Rectangle<R> 
    operator-(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::Vector< Interval<R> >& v)
    {
      Geometry::Rectangle<R> result(r.dimension());
      check_dimension(r,v.size(),__PRETTY_FUNCTION__);
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]-v(i));
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
