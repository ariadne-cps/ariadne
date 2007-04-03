/***************************************************************************
 *            point.inline.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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

#include <cassert>

namespace Ariadne {
  namespace Geometry {
   
    template<class R> inline
    Point<R>::Point()
      : _vector(0) 
    { 
    }

    template<class R> inline
    Point<R>::Point(dimension_type d)
      : _vector(d) 
    {
      for(size_type i=0; i!=dimension(); ++i) {
        _vector(i)=real_type(0);
        
      }
    }

    template<class R> template<class Rl> inline
    Point<R>::Point(dimension_type d, const Rl* data, size_type inc)
      : _vector(d,data,inc) 
    { 
    }
      
    template<class R> template<class ForwardIterator> inline
    Point<R>::Point(ForwardIterator b, ForwardIterator e)
      : _vector(std::distance(b,e))
    {
      for(size_type i=0; i!=dimension(); ++i) {
        _vector[i]=*b;
        ++b;
      }
    }

    template<class R> inline
    Point<R>::Point(const LinearAlgebra::Vector<R>& position)
      : _vector(position) 
    { 
    }

    template<class R> template<class R2> inline
    Point<R>::Point(const Point<R2>& original)
      : _vector(original.position_vector())
    {
    }
      
    template<class R> inline
    Point<R>::Point(const Point<R>& original)
      : _vector(original._vector)
    { 
    }

    template<class R> inline
    Point<R>& 
    Point<R>::operator=(const Point<R>& original) {
      if(this!=&original) { 
        this->_vector=original._vector; 
      }
      return *this; 
    }
      
 
    template<class R> inline
    bool 
    Point<R>::operator==(const Point<R>& A) const 
    {
      return this->_vector==A._vector;
    }
      
    template<class R> inline
    bool 
    Point<R>::operator!=(const Point<real_type>& A) const 
    {
      return !( *this == A );
    }

    template<class R> inline
    const array<R>&
    Point<R>::data() const 
    {
      return this->_vector.data();
    }

    template<class R> inline
    dimension_type 
    Point<R>::dimension() const 
    {
      return this->_vector.size();
    }

    template<class R> inline
    void
    Point<R>::resize(dimension_type d) 
    {
      this->_vector.resize(d);
    }

    template<class R> inline
    R& 
    Point<R>::operator[](dimension_type index) 
    {
      return  (this->_vector(index));
    }

    template<class R> inline
    const R& 
    Point<R>::operator[](dimension_type index) const 
    {
      return  (this->_vector(index));
    }

    template<class R> inline
    R& 
    Point<R>::at(dimension_type index) 
    {
      ARIADNE_CHECK_COORDINATE(*this,index,__PRETTY_FUNCTION__);
      return  (this->_vector(index));
    }

    template<class R> inline
    const R& 
    Point<R>::at(dimension_type index) const 
    {
      ARIADNE_CHECK_COORDINATE(*this,index,__PRETTY_FUNCTION__);
      return  (this->_vector(index));
    }

    template<class R> inline
    const LinearAlgebra::Vector<R>& 
    Point<R>::position_vector() const 
    {
      return this->_vector; 
    }






    template<class R> inline
    Point<R> approximation(const Point<R>& pt) 
    {
      return pt;
    }
    

    template<class R> inline
    bool contains_value(const Point< Numeric::Interval<R> >& ipt, const Point<R>& pt) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(ipt,pt,__PRETTY_FUNCTION__);
      for(dimension_type i=0; i!=ipt.dimension(); ++i) {
        if(!contains_value(ipt[i],pt[i])) {
          return false;
        }
      }
      return true;
    }
    
    
    template<class R> inline
    Point<R> approximation(const Point< Numeric::Interval<R> >& ipt) 
    {
      Point<R> result(ipt.dimension());
      for(dimension_type i=0; i!=ipt.dimension(); ++i) {
        result[i]=ipt[i].centre();
      }
      return result;
    }
    
    template<class R> inline
    R error_bound(const Point< Numeric::Interval<R> >& ipt) 
    {
      R result(0);
      for(dimension_type i=0; i!=ipt.dimension(); ++i) {
        result=Numeric::max(result,error_bound(ipt[i]));
      }
      return result;
    }
    

    template<class R> inline
    Point<typename Numeric::traits<R>::arithmetic_type>
    minkowski_sum(const Point<R>& pt1, const Point<R>& pt2) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(pt1,pt2,__PRETTY_FUNCTION__);
      return Point<typename Numeric::traits<R>::arithmetic_type>(pt1.position_vector()+pt2.position_vector());
    }
    
    template<class R> inline
    Point<typename Numeric::traits<R>::arithmetic_type>
    minkowski_difference(const Point<R>& pt1, const Point<R>& pt2) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(pt1,pt2,__PRETTY_FUNCTION__);
      return Point<typename Numeric::traits<R>::arithmetic_type>(pt1.position_vector()-pt2.position_vector());
    }
    

    template<class R1,class R2> inline
    LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator-(const Point<R1> pt1, const Point<R2>& pt2) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(pt1,pt2,__PRETTY_FUNCTION__);
      return pt1.position_vector()-pt2.position_vector();
    }
    
    template<class R1,class R2> inline
    Point<typename Numeric::traits<R1,R2>::arithmetic_type> 
    operator+(const Point<R1>& pt, const LinearAlgebra::Vector<R2>& v)
    {
      ARIADNE_CHECK_DIMENSION(pt,v.size(),__PRETTY_FUNCTION__);
      return Point<typename Numeric::traits<R1,R2>::arithmetic_type>(pt.position_vector() + v);
    }


    template<class R1,class R2> inline
    Point<typename Numeric::traits<R1,R2>::arithmetic_type> 
    operator-(const Point<R1>& pt, const LinearAlgebra::Vector<R2>& v)
    {
      ARIADNE_CHECK_DIMENSION(pt,v.size(),__PRETTY_FUNCTION__);
      return Point<typename Numeric::traits<R1,R2>::arithmetic_type>(pt.position_vector() - v);
    }

    template<class R> inline
    Point<R> 
    add_approx(const Point<R>& pt, const LinearAlgebra::Vector<R>& v)
    {
      ARIADNE_CHECK_DIMENSION(pt,v.size(),__PRETTY_FUNCTION__);
      return Point<R>(add_approx(pt.position_vector(),v));
    }

    template<class R> inline
    Point<R> 
    sub_approx(const Point<R>& pt, const LinearAlgebra::Vector<R>& v)
    {
      ARIADNE_CHECK_DIMENSION(pt,v.size(),__PRETTY_FUNCTION__);
      return Point<R>(sub_approx(pt.position_vector(),v));
    }


    
    template<class R> inline
    Point<R>
    approximate_value(const Point< Numeric::Interval<R> >& pt) 
    {
      Point<R> result(pt.dimension());
      for(dimension_type i=0; i!=result.dimension(); ++i) {
        result[i]=approximate_value(pt[i]);
      }
      return result;
    }
    

    template<class R> inline 
    Point<R> 
    project_on_dimensions(const Point<R> &A, const Base::array<bool>& dims) 
    {
      if (A.dimension()!=dims.size()) {
        throw std::runtime_error("project_on_dimensions(const Point & ,...): the two parameters have different dimension");
      }
      size_type new_dim=0;

      for (size_type i=0; i< A.dimension(); i++) {
        if (dims[i]) {
          new_dim++;
        }
      }
      Point<R> new_point(new_dim);
      
      size_type new_i=0;
      for (size_t i=0; i<dims.size(); i++) {
        if (dims[i]) {
           new_point.set(new_i,A[i]);
           new_i++;
        } 
      }

      return new_point;
    }

    template<class R> inline 
    Point<R> 
    project_on_dimensions(const Point<R> &A, 
                          const size_type& x, const size_type& y, const size_type& z) 
    {
      if ((A.dimension()<=x)||(A.dimension()<=y)||(A.dimension()<=z)) {
        throw std::runtime_error("project_on_dimensions(const Point & ,...): the two parameters have different dimension");
      }
      
      Point<R> new_point(3);
      
      new_point.set(0,A[x]);
      new_point.set(1,A[y]);
      new_point.set(2,A[z]);

      return new_point;
    }

    template<class R> inline 
    Point<R> 
    project_on_dimensions(const Point<R> &A, 
                          const size_type &x, const size_type&y) 
    {
      if ((A.dimension()<=x)||(A.dimension()<=y)) {
         throw "project_on_dimensions(const Point& ,...): one of the projection dimensions is greater than the Point dimension";
      }
      Point<R> new_point(2);
      
      new_point.set(0,A[x]);
      new_point.set(1,A[y]);

      return new_point;
    }
    



    template<class R> inline 
    std::ostream& 
    operator<<(std::ostream& os, const Point<R>& pt) 
    {
      return pt.write(os);
    }
    
    template<class R> inline
    std::istream& 
    operator>>(std::istream& is, Point<R>& pt) 
    {
      return pt.read(is);
    }


  }
}
