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

#include "../debug.h"

#include "../base/stlio.h"

#include "../numeric/interval.h"
#include "../numeric/arithmetic.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_programming/linear_program.h"

#include "../geometry/ddconv.h"
#include "../geometry/ddconv.code.h"
#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"
#include "../geometry/polytope.h"

#include "../output/logging.h"


namespace Ariadne {
  namespace Geometry {

    template<class R> 
    tribool 
    empty(const Polyhedron<R>& plhd) 
    {
      return empty(Polyhedron<Numeric::Rational>(plhd));
    }

  
    template<class R> 
    tribool 
    empty(const Polyhedron< Numeric::Interval<R> >& plhd) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
  

    template<> 
    tribool 
    empty(const Polyhedron<Numeric::Rational>& plhd) {
      try {
        Polytope<Numeric::Rational> pltp(plhd);
        ARIADNE_LOG(7,"empty(plhd): plhd="<<plhd<<" pltp="<<pltp);
        return pltp.empty();
      }
      catch(UnboundedSet& e) {
        return false;
      }
    }
  
      
    
    template<class R> 
    tribool 
    bounded(const Polyhedron<R>& plhd) 
    {
      return bounded(Polyhedron<Numeric::Rational>(plhd));
    }

  
    template<class R> 
    tribool 
    bounded(const Polyhedron< Numeric::Interval<R> >& plhd) 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
  

    template<> 
    tribool 
    bounded(const Polyhedron<Numeric::Rational>& plhd)
    {
      try {
        Polytope<Numeric::Rational> pltp(plhd);
        return pltp.bounded();
      }
      catch(UnboundedSet& e) {
        return false;
      }
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
      using namespace LinearAlgebra;
      ARIADNE_CHECK_SIZE(b,A.number_of_rows(),"Polyhedron::Polyhedron(Matrix A, Vector b)");
      dimension_type d=this->dimension();
      dimension_type nc=this->number_of_constraints();
      LinearAlgebra::MatrixSlice<R>(nc,d,this->begin(),d+1u,1u)=-A;
      LinearAlgebra::VectorSlice<R>(nc,this->begin()+d,d+1u)=b;
    }
  

    template<class R>
    Polyhedron<R>::Polyhedron(const LinearAlgebra::Matrix<R>& C) 
      : _dimension(C.number_of_columns()-1), 
        _number_of_constraints(C.number_of_rows()), 
        _data(C.data())
    {
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
    void
    Polyhedron<R>::new_constraint(const Constraint<R>& c)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,c,"void Polyhedron::new_constraint(Constraint& c)");
      dimension_type d=this->_dimension;
      size_type sz=this->_data.size();
      this->_data.resize(sz+d+1u);
      for(dimension_type i=0; i!=d; ++i) {
        this->_data[sz+i]=c._a[i];
        this->_data[sz+d]=c._a[d];
      }
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
    Polyhedron<R>::operator=(const Polyhedron<R>& plhd)
      
    {
      if(this!=&plhd) {
        this->_dimension=plhd.dimension();
        this->_number_of_constraints=plhd.number_of_constraints();
        this->_data=plhd.data();
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
    bounding_box(const Polyhedron<R>& plhd)
    {
      using Numeric::Rational;
      Rectangle<Rational> qbb=Polytope<Rational>(Polyhedron<Rational>(plhd)).bounding_box();
      Rectangle<R> bb(plhd.dimension());
      for(dimension_type i=0; i!=bb.dimension(); ++i) {
        bb.set_lower_bound(i,Numeric::conv_down<R>(qbb.lower_bound(i)));
        bb.set_upper_bound(i,Numeric::conv_up<R>(qbb.upper_bound(i)));
      }
      return bb;
    }
 
    template<class R> inline
    Rectangle< Numeric::Interval<R> > 
    bounding_box(const Polyhedron< Numeric::Interval<R> >& plhd)
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
      return Geometry::empty(*this);
    }

    template<class R>
    tribool 
    Polyhedron<R>::bounded() const
    {
      return Geometry::bounded(*this);
    }

    template<class R>
    tribool 
    Polyhedron<R>::contains(const Point<R>& pt) const
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(*this,pt,"tribool Polyhedron::contains(Point pt)");
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
    disjoint(const Polyhedron<Numeric::Rational>& plhd, const Rectangle<Numeric::Rational>& r)
    {
      typedef Numeric::Rational R;
      typedef Numeric::Rational F;

      if(verbosity>7) { std::clog << "tribool disjoint(Polyhedron<Rational>,Polyhedron<Rational>)" << std::endl; }
      if(verbosity>8) { std::clog << plhd << " " << r << std::endl; }
      
      ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd,r,"tribool disjoint(Polyhedron<Rational> plhd, Rectangle<Rational> r)");
      dimension_type d=plhd.dimension();
      size_type nc=plhd.number_of_constraints();
      size_type nnc=0; // number of negative constraints Ax<=b
      
      // Translate x'=x-l;  Ax<=b A(x'+l)<=b Ax'<= b-Al b'=b-Al 
      LinearAlgebra::Vector<F> u=r.upper_corner()-r.lower_corner();
      LinearAlgebra::Matrix<R> A=plhd.A();
      LinearAlgebra::Vector<F> b=plhd.b()-A*r.lower_corner().position_vector();
      
      if(verbosity>8) { std::clog << "u'=" << u << " A'=" << A << " b'=" << b << std::endl; }

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

      if(verbosity>8) { std::clog << T << std::endl; }

      LinearProgramming::LinearProgram<F> lp(T);
      tribool result=!lp.is_feasible();
      return result;
    }
    
    
    template<class R>
    tribool 
    disjoint(const Polyhedron< Numeric::Interval<R> >& plhd, const Rectangle< Numeric::Interval<R> >& r)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool 
    disjoint(const Polyhedron<R>& plhd, const Rectangle<R>& r)
    {
      return disjoint(Polyhedron<Numeric::Rational>(plhd),Rectangle<Numeric::Rational>(r));
    }


    template<class R>
    tribool 
    disjoint(const Rectangle<R>& r, const Polyhedron<R>& plhd)
    {
      return disjoint(plhd,r);
    }
    
    
    /*!Set up linear programming problem to solve A1*s1 <= b1; A2*s2 <= b2; A1*s1 - A2*s2 ==0 */
    template<class R>
    tribool 
    disjoint(const Polyhedron<R>& plhd1, const Polyhedron<R>& plhd2)
    {
      return closed_intersection(plhd1,plhd2).empty();
    }
    
    
    
    template<class R>
    tribool 
    subset(const Polyhedron< Numeric::Interval<R> >& plhd, const Rectangle< Numeric::Interval<R> >& r)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    tribool 
    subset(const Polyhedron< Numeric::Interval<R> >& plhd1, const Polyhedron< Numeric::Interval<R> >& plhd2)
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    tribool 
    subset(const Polyhedron<R>& plhd1, const Polyhedron<R>& plhd2)
    {
      return Geometry::subset(Polytope<Numeric::Rational>(Polyhedron<Numeric::Rational>(plhd1)),Polyhedron<Numeric::Rational>(plhd2));
    }
    
    template<class R>
    tribool 
    subset(const Polyhedron<R>& plhd, const Rectangle<R>& r)
    {
      ARIADNE_LOG(3,"Geometry::subset(plhd,r): ""plhd="<<plhd<<", r="<<r<<"\n");
      tribool result=true;
      dimension_type d=plhd.dimension();
      array<R> data(d+1u);
      for(size_type i=0; i!=d+1u; ++i) {
        data[i]=0;
      }
      Polyhedron<R> halfspace;
      ARIADNE_LOG(9,"data="<<data<<"\n");
     
      for(size_type j=0; j!=d; ++j) {
        // Check if disjoint from lower halfspace
        data[j]=-1;
        data[d]=r.lower_bound(j);
        ARIADNE_LOG(9,"data="<<data<<"\n");
        halfspace=Polyhedron<R>(d,1u,data.begin());
        ARIADNE_LOG(8,"halfspace="<<halfspace<<"\n");
        result=result && disjoint(plhd,halfspace);
        // Check if disjoint from upper halfspace
        data[j]=1;
        data[d]=-r.upper_bound(j);
        halfspace=Polyhedron<R>(d,1u,data.begin());
        ARIADNE_LOG(8,"halfspace="<<halfspace<<"\n");
        result=result && disjoint(plhd,halfspace);
        // Early return
        if(result==false) {
          return result;
        }
        data[j]=0;
      }
      return result;
    }

    
    
    
    
    template<class R1,class R2>
    tribool 
    subset(const Rectangle<R1>& r, const Polyhedron<R2>& plhd)
    {
      return Geometry::_subset(r,plhd);
    }
    
    template<class R1,class R2>
    tribool 
    subset(const Zonotope<R1>& z, const Polyhedron<R2>& plhd)
    {
      return Geometry::_subset(z,plhd);
    }
    
    template<class R1,class R2>
    tribool 
    subset(const Polytope<R1>& pltp, const Polyhedron<R2>& plhd)
    {
      return Geometry::_subset(pltp,plhd);
    }
    
    
    
    template<class R>
    tribool 
    equal(const Polyhedron<R>& plhd1, const Polyhedron<R>& plhd2)
    {
      return Geometry::subset(plhd1,plhd2) && Geometry::subset(plhd2,plhd1); 
    }
    

    template<class R>
    Polyhedron<R> 
    open_intersection(const Polyhedron<R>& plhd1, const Polyhedron<R>& plhd2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd1,plhd2,"Polyhedron open_intersection(Polyhedron plhd1, Polyhedron plhd2)");
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    
    template<class R>
    Polyhedron<R> 
    closed_intersection(const Polyhedron<R>& plhd1, const Polyhedron<R>& plhd2)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd1,plhd2,"Polyhedron closed_intersection(Polyhedron plhd1, Polyhedron plhd2)");
      dimension_type d=plhd1.dimension();
      size_type nc1=plhd1.number_of_constraints();
      size_type nc2=plhd2.number_of_constraints();
      LinearAlgebra::Matrix<R> A(nc1+nc2,d);
      LinearAlgebra::MatrixSlice<R>(nc1,d,A.begin(),d,1)=plhd1.A();
      LinearAlgebra::MatrixSlice<R>(nc2,d,A.begin()+nc1*d,d,1)=plhd2.A();
      LinearAlgebra::Vector<R> b=direct_sum(plhd1.b(),plhd2.b());
      return Polyhedron<R>(A,b);
    }
          
    
    template<class R>
    Polyhedron<R> 
    closed_intersection(const Rectangle<R>& r, const Polyhedron<R>& plhd)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(r,plhd,"Polyhedron closed_intersection(Rectangle r, Polyhedron plhd)");
      return closed_intersection(Polyhedron<R>(r),plhd);
    }
          
    
    template<class R>
    Polyhedron<R> 
    closed_intersection(const Polyhedron<R>& plhd, const Rectangle<R>& r)
    {
      ARIADNE_CHECK_EQUAL_DIMENSIONS(plhd,r,"Polyhedron closed_intersection(Polyhedron plhd, Rectangle r)");
      return closed_intersection(plhd,Polyhedron<R>(r));
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
      
      if(nv==0) {
        // empty polytope; return empty polyhedron
        array<F> data(d+1,F(0));
        data[d]=-1;
        return Polyhedron<F>(d,1,data.begin());
      }

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
    std::string
    Polyhedron<R>::name()
    {
      return std::string("Polyhedron")+"<"+Numeric::name<R>()+">";
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
      closed_intersection(r,p);
      closed_intersection(p,r);
      open_intersection(p,p);
      
      polyhedron(r);
      polyhedron(c);
    }


  }
}

