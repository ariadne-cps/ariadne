/***************************************************************************
 *            polyhedron.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include "numeric/numerical_types.h"
#include "numeric/interval.h"
#include "geometry/ppl_polyhedron.h"

#include <ppl.hh>


namespace Ariadne {
  namespace Geometry {

    Parma_Polyhedra_Library::Constraint
    ppl_constraint(const LinearAlgebra::Vector<Rational>&  a, 
                   const Rational& b) 
    {
      Parma_Polyhedra_Library::Linear_Expression e;
      Parma_Polyhedra_Library::Linear_Expression c;
      Integer multiplier = denominator(b);
      
      for (uint i = 0; i!=a.size(); ++i) {
        multiplier=Ariadne::lcm(multiplier,denominator(a[i]));
      }
      for (uint i = 0; i!=a.size(); ++i) {
        Rational coefficient = a[i];
        Integer scaled_coefficient = numerator(coefficient)*(multiplier/denominator(coefficient));
        e += Parma_Polyhedra_Library::Variable(i) * scaled_coefficient;
      }
      Integer scaled_coefficient = numerator(b)*(multiplier/denominator(b));
      c = Parma_Polyhedra_Library::Linear_Expression(scaled_coefficient);
      
      return e<=c;
    }
    
    Parma_Polyhedra_Library::Constraint
    ppl_open_constraint(const LinearAlgebra::Vector<Rational>&  a, 
                        const Rational& b) 
    {
      Parma_Polyhedra_Library::Linear_Expression e;
      Parma_Polyhedra_Library::Linear_Expression c;
      Integer multiplier = denominator(b);
      
      for (uint i = 0; i!=a.size(); ++i) {
        multiplier=Ariadne::lcm(multiplier,denominator(a[i]));
      }
      for (uint i = 0; i!=a.size(); ++i) {
        Rational coefficient = a[i];
        Integer scaled_coefficient = numerator(coefficient)*(multiplier/denominator(coefficient));
        e += Parma_Polyhedra_Library::Variable(i) * scaled_coefficient;
      }
      Integer scaled_coefficient = numerator(b)*(multiplier/denominator(b));
      c = Parma_Polyhedra_Library::Linear_Expression(scaled_coefficient);
      
      return e<c;
    }
    
    
    /* Convert the vector \a p into a Parma Polyhedra Library generator. */
    Parma_Polyhedra_Library::Generator 
    ppl_generator(const LinearAlgebra::Vector<Rational>& p) 
    {
      Integer den=1;
      for (size_type j=0; j!=p.size(); ++j) {
        den=lcm(den,denominator(p[j]));
      }
        
      Parma_Polyhedra_Library::Linear_Expression ppl_lin_expr;
      for (size_type j=0; j!=p.size(); ++j) {
        ppl_lin_expr+=( (numerator(p[j])*(den/denominator(p[j]))) * Parma_Polyhedra_Library::Variable(j) );
      }
      return Parma_Polyhedra_Library::Generator::point(ppl_lin_expr,den);
    }


    
    Parma_Polyhedra_Library::NNC_Polyhedron 
    ppl_open_polyhedron(const LinearAlgebra::Matrix<Rational>& A, 
                        const LinearAlgebra::Vector<Rational>& b)
    {
      LinearAlgebra::Vector<Rational> rw;
      Rational cnst;
      Parma_Polyhedra_Library::Constraint_System cs;
      
      for(size_type i=0; i!=A.size1(); ++i) {
        rw=row(A,i);
        cnst=b[i];
        cs.insert(ppl_open_constraint(rw,cnst));
      }
      return Parma_Polyhedra_Library::NNC_Polyhedron(cs);
    }
    
    Parma_Polyhedra_Library::NNC_Polyhedron 
    ppl_nnc_polyhedron(const LinearAlgebra::Matrix<Rational>& G)
    {
      Parma_Polyhedra_Library::Generator_System ppl_gs;
      Parma_Polyhedra_Library::Linear_Expression ppl_lin_expr;
      
      LinearAlgebra::Vector<Rational> s;
      for(size_type j=0; j!=G.size2(); ++j) {
        s=column(G,j);
        ppl_gs.insert(ppl_generator(s));
      }
      return Parma_Polyhedra_Library::NNC_Polyhedron(ppl_gs);
    }
    
    Parma_Polyhedra_Library::C_Polyhedron 
    ppl_polyhedron(const LinearAlgebra::Matrix<Rational>& A, 
                   const LinearAlgebra::Vector<Rational>& b)
    {
      LinearAlgebra::Vector<Rational> rw;
      Rational cnst;
      Parma_Polyhedra_Library::Constraint_System cs;
      
      for(size_type i=0; i!=A.size1(); ++i) {
        rw=row(A,i);
        cnst=b[i];
        cs.insert(ppl_constraint(rw,cnst));
      }
      return Parma_Polyhedra_Library::C_Polyhedron(cs);
    }
    
    // FIXME: This is inefficient
    Parma_Polyhedra_Library::C_Polyhedron 
    ppl_polyhedron(const LinearAlgebra::Vector<Rational>& c, 
                   const LinearAlgebra::Matrix<Rational>& G)
    {
      Parma_Polyhedra_Library::Generator_System ppl_gs;
      LinearAlgebra::Vector<Rational> e(G.size2());
      LinearAlgebra::Vector<Rational> v(c);
      size_type m=G.size2();
      for(size_type i=0; i!=(1u<<m); ++i) {
        for(size_type j=0; j!=m; ++j) {
          e[j] = (i&(1u<<j)) ? +1 : -1;
        }
        v=c+G*e;
        ppl_gs.insert(ppl_generator(v));
      }
      return Parma_Polyhedra_Library::C_Polyhedron(ppl_gs);
    }
    
    Parma_Polyhedra_Library::C_Polyhedron  
    ppl_polyhedron(const LinearAlgebra::Matrix<Rational>& G)
    {
      Parma_Polyhedra_Library::Generator_System ppl_gs;
      Parma_Polyhedra_Library::Linear_Expression ppl_lin_expr;
      
      LinearAlgebra::Vector<Rational> s;
      for(size_type j=0; j!=G.size2(); ++j) {
        s=column(G,j);
        ppl_gs.insert(ppl_generator(s));
      }
      return Parma_Polyhedra_Library::C_Polyhedron(ppl_gs);
    }
   
    
    Parma_Polyhedra_Library::C_Polyhedron  
    ppl_polyhedron(const LinearAlgebra::IntervalVector<Rational>& r)
    {
      Parma_Polyhedra_Library::Constraint_System ppl_cs;
      LinearAlgebra::Vector<Rational> e(r.size());
      Rational bd;
      for (size_type i=0; i<r.size(); ++i) {
        e[i]=1;
        bd=r[i].upper();
        ppl_cs.insert(ppl_constraint(e,bd));
        e[i]=-1;
        bd=r[i].lower();
        ppl_cs.insert(ppl_constraint(e,bd));
      }
      return Parma_Polyhedra_Library::C_Polyhedron(ppl_cs);
    }
    
    
    Parma_Polyhedra_Library::C_Polyhedron  
    ppl_polyhedron(const LinearAlgebra::Vector<Rational>& p)
    {      
      Parma_Polyhedra_Library::Generator_System ppl_gs;
      ppl_gs.insert(ppl_generator(p));
      return Parma_Polyhedra_Library::C_Polyhedron(ppl_gs);
    }
    
    
  
    LinearAlgebra::Matrix<Rational>
    generators(const Parma_Polyhedra_Library::C_Polyhedron& A)
    {
      LinearAlgebra::Matrix<Rational> result;
      const Parma_Polyhedra_Library::Generator_System& gs = A.minimized_generators();
      size_type d=gs.space_dimension();
      size_type n=std::distance(gs.begin(),gs.end());
      result.resize(d,n);
     
      size_type j=0;
      for(Parma_Polyhedra_Library::Generator_System::const_iterator iter=gs.begin(); iter!=gs.end(); ++iter) {
        Parma_Polyhedra_Library::Generator gen = *iter;
        Parma_Polyhedra_Library::Coefficient_traits::const_reference denom = gen.divisor();
        for(size_type i=0; i!=d; ++i) {
          Parma_Polyhedra_Library::Coefficient_traits::const_reference num = gen.coefficient(Parma_Polyhedra_Library::Variable(i));
          result(i,j)=Rational(num)/denom;
        } 
        ++j;
      }        
      return result;
    }
    
    
    LinearAlgebra::Matrix<Rational> 
    constraints_matrix(const Parma_Polyhedra_Library::C_Polyhedron& A)
    {
      LinearAlgebra::Matrix<Rational> result;
      const Parma_Polyhedra_Library::Constraint_System& cs = A.minimized_constraints();
      size_type d=cs.space_dimension();
      size_type n=std::distance(cs.begin(),cs.end());
      result.resize(n,d);
     
      size_type i=0;
      for(Parma_Polyhedra_Library::Constraint_System::const_iterator iter=cs.begin(); iter!=cs.end(); ++iter) {
        Parma_Polyhedra_Library::Constraint c = *iter;
         for(size_type j=0; j!=d; ++j) {
           result(i,j)=-c.coefficient(Parma_Polyhedra_Library::Variable(j));
        } 
        ++i;
      }        
      return result;
    }
    
    LinearAlgebra::Vector<Rational> 
    constraints_vector(const Parma_Polyhedra_Library::C_Polyhedron& A)
    {
      LinearAlgebra::Vector<Rational> result;
      const Parma_Polyhedra_Library::Constraint_System& cs = A.minimized_constraints();
      size_type n=std::distance(cs.begin(),cs.end());
      result.resize(n);
     
      size_type i=0;
      for(Parma_Polyhedra_Library::Constraint_System::const_iterator iter=cs.begin(); iter!=cs.end(); ++iter) {
        const Parma_Polyhedra_Library::Constraint& c = *iter;
        result(i)=c.inhomogeneous_term();
        ++i;
      }        
      return result;
    }
    
     
    
    dimension_type dimension(const Parma_Polyhedra_Library::C_Polyhedron& A)
    {
      return A.space_dimension();
    }
    
                  
    bool empty(const Parma_Polyhedra_Library::C_Polyhedron& A)
    {
      return A.is_empty();
    }
    
                  
    bool empty_interior(const Parma_Polyhedra_Library::C_Polyhedron& A)
    {
      return A.space_dimension()!=A.affine_dimension();
    }
                  
    
    
    bool equal(const Parma_Polyhedra_Library::C_Polyhedron& A,
               const Parma_Polyhedra_Library::C_Polyhedron& B)
    {
      return A==B;
    }
                  
    bool disjoint(const Parma_Polyhedra_Library::C_Polyhedron& A,
                  const Parma_Polyhedra_Library::C_Polyhedron& B)
    {
      return A.is_disjoint_from(B);
    }
    
    bool interiors_intersect(const Parma_Polyhedra_Library::C_Polyhedron& A,
                             const Parma_Polyhedra_Library::C_Polyhedron& B)
    {
      return !empty_interior(intersection(A,B));
    }
    
    bool inner_subset(const Parma_Polyhedra_Library::C_Polyhedron& A,
                      const Parma_Polyhedra_Library::C_Polyhedron& B)
    {
      return interior(B).contains(closure(A));
    }
    
    bool subset(const Parma_Polyhedra_Library::C_Polyhedron& A,
                const Parma_Polyhedra_Library::C_Polyhedron& B)
    {
      return B.contains(A);
    }
     
    

    Parma_Polyhedra_Library::NNC_Polyhedron 
    interior(const Parma_Polyhedra_Library::C_Polyhedron& A)
    {     
      return ppl_open_polyhedron(constraints_matrix(A),constraints_vector(A));
    }
     
    Parma_Polyhedra_Library::NNC_Polyhedron 
    closure(const Parma_Polyhedra_Library::C_Polyhedron& A)
    {     
      return Parma_Polyhedra_Library::NNC_Polyhedron(A.generators());
    }
     
   
    
    Parma_Polyhedra_Library::C_Polyhedron
    intersection(const Parma_Polyhedra_Library::C_Polyhedron& A,
                 const Parma_Polyhedra_Library::C_Polyhedron& B)
    {
      Parma_Polyhedra_Library::C_Polyhedron result(A);
      result.intersection_assign(B);
      return result;
    }
                  
    Parma_Polyhedra_Library::C_Polyhedron
    regular_intersection(const Parma_Polyhedra_Library::C_Polyhedron& A,
                         const Parma_Polyhedra_Library::C_Polyhedron& B)
    {
      Parma_Polyhedra_Library::C_Polyhedron result=intersection(A,B);
      if(empty_interior(result)) { result=Parma_Polyhedra_Library::C_Polyhedron(A.space_dimension()); }
      return result;
    }
                
    Parma_Polyhedra_Library::C_Polyhedron
    convex_hull(const Parma_Polyhedra_Library::C_Polyhedron& A,
                const Parma_Polyhedra_Library::C_Polyhedron& B)
    {
      Parma_Polyhedra_Library::C_Polyhedron result(A); 
      result.poly_hull_assign_and_minimize(B);
      return result; 
    }
                       
    /* FIXME: Reimplement */
    Parma_Polyhedra_Library::C_Polyhedron
    minkowski_sum(const Parma_Polyhedra_Library::C_Polyhedron& A, 
                  const Parma_Polyhedra_Library::C_Polyhedron& B)
    {
      assert(A.space_dimension()==B.space_dimension());
      LinearAlgebra::Matrix<Rational> Av=generators(A);
      LinearAlgebra::Matrix<Rational> Bv=generators(B);
      LinearAlgebra::Matrix<Rational> Rv(Av.size1(),Av.size2()*Bv.size2());
      for(size_type i=0; i!=Rv.size1(); ++i) {
        for(size_type j=0; j!=Av.size2(); ++j) {
          for(size_type k=0; k!=Bv.size2(); ++k) {
            Rv(i,j*Bv.size2()+k)=Av(i,j)+Bv(i,k);
          }
        }
      }
      return ppl_polyhedron(Rv);
    }
  
    Parma_Polyhedra_Library::C_Polyhedron
    minkowski_difference(const Parma_Polyhedra_Library::C_Polyhedron& A, 
                         const Parma_Polyhedra_Library::C_Polyhedron& B)
    {
      assert(A.space_dimension()==B.space_dimension());
      LinearAlgebra::Matrix<Rational> Av=generators(A);
      LinearAlgebra::Matrix<Rational> Bv=generators(B);
      LinearAlgebra::Matrix<Rational> Rv(Av.size1(),Av.size2()*Bv.size2());
      for(size_type i=0; i!=Rv.size1(); ++i) {
        for(size_type j=0; j!=Av.size2(); ++j) {
          for(size_type k=0; k!=Bv.size2(); ++k) {
            Rv(i,j*Bv.size2()+k)=Av(i,j)-Bv(i,k);
          }
        }
      }
      return ppl_polyhedron(Rv);
    }
                   
  }
}
