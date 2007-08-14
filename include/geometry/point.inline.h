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
      ARIADNE_CHECK_COORDINATE(*this,index,"R& Point at(dimension_type index)");
      return  (this->_vector(index));
    }

    template<class R> inline
    const R& 
    Point<R>::at(dimension_type index) const 
    {
      ARIADNE_CHECK_COORDINATE(*this,index,"const R& Point at(dimension_type index) const");
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
    Point<R> approximation(const Point< Numeric::Interval<R> >& ipt) 
    {
      Point<R> result(ipt.dimension());
      for(dimension_type i=0; i!=ipt.dimension(); ++i) {
        result[i]=midpoint(ipt[i]);
      }
      return result;
    }
    

    template<class R> inline
    Point<R>
    midpoint(const Point< Numeric::Interval<R> >& pt) 
    {
      Point<R> result(pt.dimension());
      for(dimension_type i=0; i!=result.dimension(); ++i) {
        result[i]=midpoint(pt[i]);
      }
      return result;
    }
    

    template<class R> inline
    R radius(const Point< Numeric::Interval<R> >& ipt) 
    {
      R result(0);
      for(dimension_type i=0; i!=ipt.dimension(); ++i) {
        result=Numeric::max(result,radius(ipt[i]));
      }
      return result;
    }


    template<class R> inline
    bool encloses(const Point< Numeric::Interval<R> >& ipt, const Point<R>& pt) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(ipt,pt,"bool encloses(Point<Interval>,Point<Real>)");
      for(dimension_type i=0; i!=ipt.dimension(); ++i) {
        if(!encloses(ipt[i],pt[i])) {
          return false;
        }
      }
      return true;
    }
    

    template<class R> inline
    bool refines(const Point< Numeric::Interval<R> >& ipt1, const Point< Numeric::Interval<R> >& ipt2) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(ipt1,ipt2,"bool refines(Point<Interval>,Point<Interval>)");
      for(dimension_type i=0; i!=ipt1.dimension(); ++i) {
        if(!refines(ipt1[i],ipt2[i])) {
          return false;
        }
      }
      return true;
    }
    

    template<class R> inline
    Point<typename Numeric::traits<R>::arithmetic_type>
    minkowski_sum(const Point<R>& pt1, const Point<R>& pt2) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(pt1,pt2,"Point minkowski_sum(Point,Point)");
      return Point<typename Numeric::traits<R>::arithmetic_type>(pt1.position_vector()+pt2.position_vector());
    }
    
    template<class R> inline
    Point<typename Numeric::traits<R>::arithmetic_type>
    minkowski_difference(const Point<R>& pt1, const Point<R>& pt2) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(pt1,pt2,"Point minkowski_difference(Point,Point)");
      return Point<typename Numeric::traits<R>::arithmetic_type>(pt1.position_vector()-pt2.position_vector());
    }
    

    template<class R1,class R2> inline
    LinearAlgebra::Vector<typename Numeric::traits<R1,R2>::arithmetic_type>
    operator-(const Point<R1>& pt1, const Point<R2>& pt2) 
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(pt1,pt2,"Vector operator-(Point,Point)");
      return pt1.position_vector()-pt2.position_vector();
    }
    
    template<class R1,class R2> inline
    Point<typename Numeric::traits<R1,R2>::arithmetic_type> 
    operator+(const Point<R1>& pt, const LinearAlgebra::Vector<R2>& v)
    {
      ARIADNE_CHECK_DIMENSION(pt,v.size(),"Point operator+(Point,Vector)");
      return Point<typename Numeric::traits<R1,R2>::arithmetic_type>(pt.position_vector() + v);
    }


    template<class R1,class R2> inline
    Point<typename Numeric::traits<R1,R2>::arithmetic_type> 
    operator-(const Point<R1>& pt, const LinearAlgebra::Vector<R2>& v)
    {
      ARIADNE_CHECK_DIMENSION(pt,v.size(),"Point operator-(Point,Vector)");
      return Point<typename Numeric::traits<R1,R2>::arithmetic_type>(pt.position_vector() - v);
    }

    template<class R> inline
    Point<R> 
    add_approx(const Point<R>& pt, const LinearAlgebra::Vector<R>& v)
    {
      ARIADNE_CHECK_DIMENSION(pt,v.size(),"Point add_approx(Point,Vector)");
      return Point<R>(add_approx(pt.position_vector(),v));
    }

    template<class R> inline
    Point<R> 
    sub_approx(const Point<R>& pt, const LinearAlgebra::Vector<R>& v)
    {
      ARIADNE_CHECK_DIMENSION(pt,v.size(),"Point sub_approx(Point,Vector)");
      return Point<R>(sub_approx(pt.position_vector(),v));
    }


    

    template<class R> inline 
    Point<R> 
    project_on_dimensions(const Point<R> &pt, const Base::array<bool>& dims) 
    {
      if (pt.dimension()!=dims.size()) {
        ARIADNE_THROW(IncompatibleDimensions,"Point project_on_dimensions(Point pt, BooleanArray dims)","pt.dimension()="<<pt<<", dims.size()="<<dims.size());
      }
      size_type new_dim=0;

      for (size_type i=0; i< pt.dimension(); i++) {
        if (dims[i]) {
          new_dim++;
        }
      }
      Point<R> new_point(new_dim);
      
      size_type new_i=0;
      for (size_t i=0; i<dims.size(); i++) {
        if (dims[i]) {
           new_point.set(new_i,pt[i]);
           new_i++;
        } 
      }

      return new_point;
    }

    template<class R> inline 
    Point<R> 
    project_on_dimensions(const Point<R> &pt, 
                          const size_type& x, const size_type& y, const size_type& z) 
    {
      if ((pt.dimension()<=x)||(pt.dimension()<=y)||(pt.dimension()<=z)) {
        ARIADNE_THROW(InvalidIndex,"Point project_on_dimensions(Point pt, size_type x, size_type y, size_type z)",
                      "pt.dimension()="<<pt<<", x="<<x<<", y="<<y<<", z="<<z);
      }
      
      Point<R> new_point(3);
      
      new_point.set(0,pt[x]);
      new_point.set(1,pt[y]);
      new_point.set(2,pt[z]);

      return new_point;
    }

    template<class R> inline 
    Point<R> 
    project_on_dimensions(const Point<R> &pt, 
                          const size_type &x, const size_type&y) 
    {
      if ((pt.dimension()<=x)||(pt.dimension()<=y)) {
        ARIADNE_THROW(InvalidIndex,"Point project_on_dimensions(Point pt, size_type x, size_type y)",
                      "pt.dimension()="<<pt<<", x="<<x<<", y="<<y);
      }
      Point<R> new_point(2);
      
      new_point.set(0,pt[x]);
      new_point.set(1,pt[y]);

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
