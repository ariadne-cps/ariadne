/***************************************************************************
 *            polyhedron.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 *  along with this program; if not, write to bouthe Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "polyhedron.h"

#include "../base/stlio.h"

#include "../numeric/interval.h"
#include "../numeric/arithmetic.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/linear_program.h"

#include "../geometry/ddconv.h"
#include "../geometry/ddconv.code.h"
#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"
#include "../geometry/polytope.h"

namespace Ariadne {
  namespace Geometry {

    
    template<class R1, class R2, template<class> class BS>
    tribool 
    _subset(const BS<R1>& A, const Polyhedron<R2>& B)
    {
      tribool result=true;
      for(typename BS<R1>::vertices_const_iterator v=A.vertices_begin();
          v!=A.vertices_end(); ++v)
      {
        for(typename Polyhedron<R2>::constraints_const_iterator c=B.constraints_begin();
            c!=B.constraints_end(); ++c)
        {
          result=result && (c->satisfied_by(*v));
          if(!result) { return result; }
        }
      }
      return result;
    }
    
    
  
    template<class R>
    Polyhedron<R>::Polyhedron(dimension_type d)
      : _dimension(d), 
        _number_of_constraints(0), 
        _data()
    {
    }
    
    
    template<class R>
    Polyhedron<R>::Polyhedron(dimension_type d, size_type nc, const R* data)
      : _dimension(d), 
        _number_of_constraints(nc), 
        _data(data,data+(d+1)*nc)
    { 
    }
   
   
    template<class R>
    Polyhedron<R>::Polyhedron(const LinearAlgebra::Matrix<R>& A,
                              const LinearAlgebra::Vector<R>& b) 
      : _dimension(A.number_of_columns()), 
        _number_of_constraints(A.number_of_rows()), 
        _data((A.number_of_columns()+1)*A.number_of_rows())
    {
      LinearAlgebra::check_size(b,A.number_of_rows(),__PRETTY_FUNCTION__);
      dimension_type d=this->dimension();
      dimension_type nc=this->number_of_constraints();
      LinearAlgebra::MatrixSlice<R>(nc,d,this->begin(),d+1u,1u)=-A;
      LinearAlgebra::VectorSlice<R>(nc,this->begin()+d,d+1u)=b;
    }
  
    template<class R>
    Polyhedron<R>::Polyhedron(const PointList<R>& pts)
      : _dimension(pts.dimension()), _number_of_constraints(0), _data()
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
  
    template<>
    Polyhedron<Numeric::Rational>::Polyhedron(const PointList<Numeric::Rational>& pts)
    {
      (*this)=Polyhedron<Numeric::Rational>(Polytope<Numeric::Rational>(pts));
    }
  

    template<class R>
    LinearAlgebra::MatrixSlice<R>
    Polyhedron<R>::_constraints()
    {
      return LinearAlgebra::MatrixSlice<R>(this->_number_of_constraints,
                                           this->_dimension+1u,
                                           this->begin(),
                                           this->_dimension+1u,
                                           1u);
    }


    template<class R>
    Polyhedron<R>::Polyhedron(const Rectangle<R>& r)
      : _dimension(r.dimension()), 
        _number_of_constraints(r.dimension()*2u),
        _data((r.dimension()+1u)*r.dimension()*2u,static_cast<R>(0))
    {
      dimension_type d=r.dimension();
      LinearAlgebra::MatrixSlice<R> constraints=this->_constraints();
      for(size_type i=0; i!=d; ++i) {
        constraints(i,i)=static_cast<R>(1);
        constraints(i,d)=-r.lower_bound(i);
        constraints(i+d,i)=static_cast<R>(-1); 
        constraints(i+d,d)=r.upper_bound(i);
      }
    }
   

    template<class R>
    Polyhedron<R>&
    Polyhedron<R>::operator=(const Polyhedron<R>& p)
      
    {
      if(this!=&p) {
        this->_dimension=p.dimension();
        this->_number_of_constraints=p.number_of_constraints();
        this->_data=p.data();
      }
      return *this;
    }
    
    template<class R>
    typename Polyhedron<R>::constraints_const_iterator
    Polyhedron<R>::constraints_begin() const
    {
      return constraints_const_iterator(*this,0u);
    }
 
    template<class R>
    typename Polyhedron<R>::constraints_const_iterator
    Polyhedron<R>::constraints_end() const
    {
      return constraints_const_iterator(*this,this->number_of_constraints());
    }
 
    template<class R>
    dimension_type
    Polyhedron<R>::dimension() const
    {
      return this->_dimension;
    }



    template<class R> inline
    Rectangle<R> 
    bounding_box(const Polyhedron<R>& p)
    {
      using Numeric::Rational;
      Rectangle<Rational> qbb=Polytope<Rational>(Polyhedron<Rational>(p)).bounding_box();
      Rectangle<R> bb(p.dimension());
      for(dimension_type i=0; i!=bb.dimension(); ++i) {
        bb.set_lower_bound(i,Numeric::conv_down<R>(qbb.lower_bound(i)));
        bb.set_upper_bound(i,Numeric::conv_up<R>(qbb.upper_bound(i)));
      }
      return bb;
    }
 
    template<class R> inline
    Rectangle< Numeric::Interval<R> > 
    bounding_box(const Polyhedron< Numeric::Interval<R> >& p)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    Rectangle<R> 
    Polyhedron<R>::bounding_box() const
    {
      return Geometry::bounding_box(*this);
    }



    template<class R>
    tribool 
    Polyhedron<R>::empty() const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool 
    Polyhedron<R>::bounded() const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool 
    Polyhedron<R>::contains(const Point<R>& pt) const
    {
      check_equal_dimensions(*this,pt,__PRETTY_FUNCTION__);
      tribool result=true;
      for(constraints_const_iterator i=this->constraints_begin(); 
          i!=this->constraints_end(); ++i) {
        result=result && i->satisfied_by(pt); 
        if(!result) { return result; }
      }
      return result;
    }
    
    
    
    template<>
    tribool 
    disjoint(const Polyhedron<Numeric::Rational>& p, const Rectangle<Numeric::Rational>& r)
    {
      typedef Numeric::Rational R;
      typedef Numeric::Rational F;

      if(verbosity>7) { std::cerr << __PRETTY_FUNCTION__ << std::endl; }
      if(verbosity>8) { std::cerr << p << " " << r << std::endl; }
      
      check_equal_dimensions(p,r,__PRETTY_FUNCTION__);
      dimension_type d=p.dimension();
      size_type nc=p.number_of_constraints();
      size_type nnc=0; // number of negative constraints Ax<=b
      
      // Translate x'=x-l;  Ax<=b A(x'+l)<=b Ax'<= b-Al b'=b-Al 
      LinearAlgebra::Vector<F> u=r.upper_corner()-r.lower_corner();
      LinearAlgebra::Matrix<R> A=p.A();
      LinearAlgebra::Vector<F> b=p.b()-A*r.lower_corner().position_vector();
      
      if(verbosity>8) { std::cerr << "u'=" << u << " A'=" << A << " b'=" << b << std::endl; }

      // count number of negative constraints
      for(size_type i=0; i!=nc; ++i) {
        if(b(i)<0) {
          ++nnc;
        }
      }

      LinearAlgebra::Matrix<F> T(nc+d+1,d+nnc+1);
      // set up constraints Ax<=b
      size_type k=0; // negative constraint number
      for(size_type i=0; i!=nc; ++i) {
        if(b(i)<0) {
          for(size_type j=0; j!=d; ++j) {
            T(i,j)=-A(i,j);
          }
          T(i,d+k)=-1;
          T(i,d+nnc)=-b(i);
          ++k;
        } else {
          for(size_type j=0; j!=d; ++j) {
            T(i,j)=A(i,j);
          }
          T(i,d+nnc)=b(i);
        }          
      }
      //
      for(size_type i=0; i!=d; ++i) {
        T(nc+i,i)=1;
        T(nc+i,d+nnc)=u(i);
      }
      
      // SetInterface value function for feasibility problem
      k=0;
      for(size_type i=0; i!=nc; ++i) {
        if(b(i)<0) {
          for(size_type j=0; j!=d+nnc; ++j) {
            T(nc+d,j)-=T(i,j);
          }
          T(nc+d,d+nnc)-=T(i,d+nnc);
        }
      }

      if(verbosity>8) { std::cerr << T << std::endl; }

      LinearAlgebra::LinearProgram<F> lp(T);
      tribool result=!lp.is_feasible();
      return result;
    }
    
    
    template<class R>
    tribool 
    disjoint(const Polyhedron< Numeric::Interval<R> >& p, const Rectangle< Numeric::Interval<R> >& r)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool 
    disjoint(const Polyhedron<R>& p, const Rectangle<R>& r)
    {
      return disjoint(Polyhedron<Numeric::Rational>(p),Rectangle<Numeric::Rational>(r));
    }


    template<class R>
    tribool 
    disjoint(const Rectangle<R>& r, const Polyhedron<R>& p)
    {
      return disjoint(p,r);
    }
    
    
    /*!Set up linear programming problem to solve A1*s1 <= b1; A2*s2 <= b2; A1*s1 - A2*s2 ==0 */
    template<class R>
    tribool 
    disjoint(const Polyhedron<R>& plyhd1, const Polyhedron<R>& plyhd2)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
    
    template<class R>
    tribool 
    subset(const Polyhedron< Numeric::Interval<R> >& A, const Rectangle< Numeric::Interval<R> >& B)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    tribool 
    subset(const Polyhedron< Numeric::Interval<R> >& A, const Polyhedron< Numeric::Interval<R> >& B)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    tribool 
    subset(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return Geometry::subset(Polytope<Numeric::Rational>(Polyhedron<Numeric::Rational>(A)),Polyhedron<Numeric::Rational>(B));
    }
    
    template<class R>
    tribool 
    subset(const Polyhedron<R>& A, const Rectangle<R>& B)
    {
      return Geometry::subset(Polytope<Numeric::Rational>(A),Rectangle<Numeric::Rational>(B));
    }
    
    
    
    
    template<class R1,class R2>
    tribool 
    subset(const Rectangle<R1>& A, const Polyhedron<R2>& B)
    {
      return Geometry::_subset(A,B);
    }
    
    template<class R1,class R2>
    tribool 
    subset(const Zonotope<R1>& A, const Polyhedron<R2>& B)
    {
      return Geometry::_subset(A,B);
    }
    
    template<class R1,class R2>
    tribool 
    subset(const Polytope<R1>& A, const Polyhedron<R2>& B)
    {
      return Geometry::_subset(A,B);
    }
    
    
    
    template<class R>
    tribool 
    equal(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return Geometry::subset(A,B) && Geometry::subset(B,A); 
    }
    

    template<class R>
    Polyhedron<R> 
    open_intersection(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      check_equal_dimensions(A,B,__PRETTY_FUNCTION__);
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
    template<class R>
    Polyhedron<R> 
    closed_intersection(const Polyhedron<R>& p1, const Polyhedron<R>& p2)
    {
      check_equal_dimensions(p1,p2,__PRETTY_FUNCTION__);
      dimension_type d=p1.dimension();
      size_type nc1=p1.number_of_constraints();
      size_type nc2=p2.number_of_constraints();
      LinearAlgebra::Matrix<R> A(nc1+nc2,d);
      LinearAlgebra::MatrixSlice<R>(nc1,d,A.begin(),d,1)=p1.A();
      LinearAlgebra::MatrixSlice<R>(nc2,d,A.begin()+nc1*d,d,1)=p2.A();
      LinearAlgebra::Vector<R> b=direct_sum(p1.b(),p2.b());
      return Polyhedron<R>(A,b);
    }
          
    
    template<class R>
    Polyhedron<R> 
    polyhedron(const Rectangle<R>& r)
    {
      return Polyhedron<R>(r);
    }
    
    
    template<class R>
    Polyhedron<typename Numeric::traits<R>::arithmetic_type>
    polyhedron(const Polytope<R>& pltp)
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      
      dimension_type d=pltp.dimension();
      size_type nv=pltp.number_of_vertices();
      
      const LinearAlgebra::Matrix<R> G=pltp.generators();
      
      std::vector< LinearAlgebra::Vector<F> > result;
      std::vector< LinearAlgebra::Vector<F> > argument;
      
      LinearAlgebra::Vector<F> tmp(d+1);
      
      for(size_type j=0; j!=nv; ++j) {
        for(size_type i=0; i!=d+1u; ++i) {
          tmp(i)=G(i,j);
        }
        argument.push_back(tmp);
      }
      ddconv(result,argument);     
      
      size_type nc=result.size();
      array<F> data((d+1)*nc);
      for(size_type i=0; i!=nc; ++i) {
        for(size_type j=0; j!=d+1u; ++j) {
          data[i*(d+1u)+j]=result[i](j);
        }
      }
      
      return Polyhedron<F>(d,nc,data.begin());
    }
    

    
    template<class R>  
    std::ostream& 
    Constraint<R>::write(std::ostream& os) const
    {
      LinearAlgebra::Vector<R> a=-LinearAlgebra::Vector<R>(this->_d,this->_a);
      R b=this->_a[this->_d];
      os << "Constraint(a=" << a << ",b=" << b << ")";
      return os;
    }
    
    template<class R>  
    std::ostream& 
    Polyhedron<R>::write(std::ostream& os) const
    {
      //return os << "Polyhedron( A=" << this->A() << ", b=" << this->b() << " )";
      os << "Polyhedron( constraints=";
      dimension_type d=this->dimension();
      size_type nc=this->number_of_constraints();
      for(size_type i=0; i!=nc; ++i) {
        os << ( i==0 ? "[" : "," );
        for(size_type j=0; j!=d; ++j) {
          os << ( j==0 ? "(" : ",");
          os << this->_data[i*(d+1)+j]; 
        }
        os << ":" << this->_data[i*(d+1)+d] << ")";
      }
      os << "] )";
      return os;
    }
    
    template<class R>  
    std::istream& 
    Polyhedron<R>::read(std::istream& is) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }


    template<class R>
    void
    Polyhedron<R>::_instantiate_geometry_operators() 
    {
      Rectangle<R> r;
      Zonotope<R> z;
      Polytope<R> c;
      Polyhedron<R> p;
      
      equal(p,p);
      disjoint(r,p);
      disjoint(p,r);
      disjoint(p,p);
      subset(r,p);
      subset(z,p);
      subset(c,p);
      subset(p,p);
      subset(p,r);
  
      closed_intersection(p,p);
      open_intersection(p,p);
      
      polyhedron(r);
      polyhedron(c);
    }


  }
}
