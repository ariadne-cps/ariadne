/***************************************************************************
 *            polyhedron.tpl
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

#include "../geometry/ddconv.h"
#include "../geometry/ddconv.tpl"
#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"
#include "../geometry/polytope.h"

namespace Ariadne {
  namespace Geometry {

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
    }
    
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
    Polyhedron<R>::Polyhedron(dimension_type n) : _A(0,n), _b(0)
    {
    }
    
    template<class R>
    Polyhedron<R>::Polyhedron(const LinearAlgebra::Matrix<R>& A,
                              const LinearAlgebra::Vector<R>& b) 
      : _A(A), _b(b)
    {
      assert(A.number_of_rows()==b.size());
    }
    
    template<class R>
    Polyhedron<R>::Polyhedron(const PointList<R>& ptl)
      : _A(), _b()
    {
      (*this)=Polyhedron<R>(Polytope<R>(ptl));
    }
    
    template<class R>
    Polyhedron<R>::Polyhedron(const Rectangle<R>& r)
      : _A(2*r.dimension(),r.dimension()), _b(r.dimension())
    {
      dimension_type n=r.dimension();
      for(size_type i=0; i!=n; ++i) {
        this->_A(i,i)=R(1); 
        this->_A(i+n,i)=R(-1); 
        this->_b(i)=r.upper_bound(i);
        this->_b(i+n)=r.lower_bound(i);
      }
    }
   
    template<class R>
    Polyhedron<R>::Polyhedron(const Polytope<R>& pltp)
      : _A(), _b()
    {   
      throw std::runtime_error("Polyhedron<R>::Polyhedron(const Polytope<R>& pltp) not implemented");
    }
    
    template<>
    Polyhedron<Rational>::Polyhedron(const Polytope<Rational>& pltp)
      : _A(), _b()
    {
      typedef Rational R;
      typedef Rational F;
      
      dimension_type d=pltp.dimension();
      size_type nv=pltp.number_of_vertices();
      
      LinearAlgebra::Matrix<R>& A=this->_A;
      LinearAlgebra::Vector<R>& b=this->_b;
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
      A.resize(nc,d);
      b.resize(nc);
      for(size_type i=0; i!=nc; ++i) {
        for(size_type j=0; j!=d; ++j) {
          A(i,j)=argument[i](j);
        }
        b(i)=argument[i](d);
      }
    }
    
    template<class R>
    Polyhedron<R>::Polyhedron(const Polyhedron<R>& p)
      : _A(p._A), _b(p._b)
    {
    }
   
    template<class R>
    Polyhedron<R>&
    Polyhedron<R>::operator=(const Polyhedron<R>& p)
      
    {
      if(this!=&p) { 
        this->_A=p._A;
        this->_b=p._b;
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
      return this->_A.number_of_columns(); 
    }

    template<class R>
    Rectangle<R> 
    Polyhedron<R>::bounding_box() const 
    {
      throw std::runtime_error("Polyhedron<R>::bounding_box() const not implemented");
    }
      
    template<class R>
    PointList<Rational>
    Polyhedron<R>::vertices() const
    {
      return Polytope<Rational>(*this).vertices();
    }

    template<class R>
    tribool 
    Polyhedron<R>::empty() const
    {
      throw std::runtime_error("Polyhedron<R>::empty() const");
    }

    template<class R>
    tribool 
    Polyhedron<R>::contains(const Point<R>& pt) const
    {
      tribool result=true;
      for(constraints_const_iterator i=this->constraints_begin(); 
          i!=this->constraints_end(); ++i) {
        result=result && i->satisfied_by(pt); 
        if(!result) { return result; }
      }
      return result;
    }
    
    template<class R>
    tribool 
    equal(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return Geometry::subset(A,B) && Geometry::subset(B,A); 
    }
    
    
    
    template<class R>
    tribool 
    disjoint(const Polyhedron<R>& p, const Rectangle<R>& r)
    {
      throw("disjoint(const Polyhedron<R>&, const Rectangle<R>&) not implemented");
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
      throw("disjoint(const Polyhedron<R>&, const Polyhedron<R>&) not implemented");
      
      typedef Rational F;
      LinearAlgebra::Matrix<F> qA(plyhd1.number_of_constraints() + plyhd2.number_of_constraints(), plyhd1.dimension());
      
      LinearAlgebra::Matrix<F> qA1=plyhd1.A();
      LinearAlgebra::Vector<F> qb1=plyhd1.b();
      LinearAlgebra::Matrix<F> qA2=plyhd2.A();
      LinearAlgebra::Vector<F> qb2=plyhd2.b();
      
      
      //if(Geometry::disjoint(C_Polyhedron(A),C_Polyhedron(B))) { return true; }
      //if(Geometry::interiors_intersect(C_Polyhedron(A),C_Polyhedron(B))) { return false; }
      //return indeterminate;
    }
    
    
    
    template<class R>
    tribool 
    subset(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      return Geometry::subset(Polytope<Rational>(Polyhedron<Rational>(A)),Polyhedron<Rational>(B));
    }
    
    template<class R>
    tribool 
    subset(const Polyhedron<R>& A, const Rectangle<R>& B)
    {
      return Geometry::subset(Polytope<Rational>(A),Rectangle<Rational>(B));
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
    Polyhedron<R> 
    open_intersection(const Polyhedron<R>& A, const Polyhedron<R>& B)
    {
      if(A.dimension()!=B.dimension()) {
        throw std::runtime_error("Incompatible dimensions in open_intersection(const Polyhedron<R>&, const Polyhedron<R>&)");
      }
      throw std::runtime_error("open_intersection(const Polyhedron<R>&, const Polyhedron<R>&) not implemented");
    }
    
    
    template<class R>
    Polyhedron<R> 
    closed_intersection(const Polyhedron<R>& p1, const Polyhedron<R>& p2)
    {
      if(p1.dimension()!=p1.dimension()) {
        throw std::runtime_error("Incompatible dimensions in closed_intersection(const Polyhedron<R>&, const Polyhedron<R>&)");
      }
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
    std::ostream& 
    Constraint<R>::write(std::ostream& os) const
    {
      os << "Constraint(a=";
      Utility::write_sequence(os,this->_a,this->_a+this->_d);
      os << ",b=" << *this->_b << ")";
      return os;
    }
    
    template<class R>  
    std::ostream& 
    Polyhedron<R>::write(std::ostream& os) const
    {
      return os << "Polyhedron( A=" << this->A() << ", b=" << this->b() << " )";
    }
    
    template<class R>  
    std::istream& 
    Polyhedron<R>::read(std::istream& is) 
    {
      throw std::runtime_error("std::istream& operator>>(std::istream&, Polyhedron<R>&) not implemented");
    }


  }
}
