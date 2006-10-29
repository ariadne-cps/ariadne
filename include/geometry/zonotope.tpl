/***************************************************************************
 *            zonotope.tpl
 *
 *  6 February 2006
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
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#include <iostream>
#include <vector>
#include <algorithm>

#include <ppl.hh>

#include "zonotope.h"

#include "../base/array.h"
#include "../base/exceptions.h"
#include "../numeric/conversion.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/lu_matrix.h"
#include "../linear_algebra/linear_program.h"

#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/polytope.h"
#include "../geometry/polyhedron.h"
#include "../geometry/parallelotope.h"
#include "../geometry/list_set.h"

// include for column operations
//#include "../linear_algebra/matrix.tpl" 

namespace Ariadne {
  namespace Geometry {
    
    template<class R>
    void
    Zonotope<R>::_instantiate_geometry_operators() 
    {
      Rectangle<R> r;
      Zonotope<R> z;
      Geometry::disjoint(r,z);
      Geometry::disjoint(z,r);
      Geometry::disjoint(z,z);
      Geometry::subset(r,z);
      Geometry::subset(z,r);
      Geometry::subset(z,z);
      Geometry::minkowski_sum(z,z);
      Geometry::minkowski_difference(z,z);
    }
    
    template<class R>
    void
    Zonotope< Interval<R> >::_instantiate_geometry_operators() 
    {
      Zonotope< Interval<R> > iz;
      Geometry::over_approximation(iz);
    }
      
    template<class R>
    Zonotope<R>::Zonotope(const Rectangle<R>& r)
      : _centre(r.dimension()), _generators(r.dimension(),r.dimension())
    {
      for(size_type i=0; i!=dimension(); ++i) {
        this->_centre[i] = r[i].centre();
        this->_generators(i,i) = r[i].radius();
      }
    }

    template<class R>
    Zonotope< Interval<R> >::Zonotope(const Rectangle<R>& r)
      : _centre(r.dimension()), _generators(r.dimension(),r.dimension())
    {
      for(size_type i=0; i!=dimension(); ++i) {
        this->_centre[i] = r[i].centre();
        this->_generators(i,i) = r[i].radius();
      }
    }

    
    template<class R>
    Zonotope<R>::Zonotope(const std::string& s)
      : _centre(), _generators()
    {
      std::stringstream ss(s);
      ss >> *this;
    }


 
    template<class R>
    typename Zonotope<R>::vertices_const_iterator
    Zonotope<R>::vertices_begin() const 
    {
      return ZonotopeVerticesIterator<R>(*this,false);
    }
    
    
    template<class R>
    typename Zonotope<R>::vertices_const_iterator
    Zonotope<R>::vertices_end() const 
    {
      return ZonotopeVerticesIterator<R>(*this,true);
    }
    
    
    template<class R>
    PointList<Rational>
    Zonotope<R>::vertices() const
    {
      //std::cerr << "Zonotope<R>::vertices()" << std::endl;
      Zonotope<Rational> qz(*this);
      PointList<Rational> qv(qz.dimension());
      for(Zonotope<Rational>::vertices_const_iterator v=qz.vertices_begin();
          v!=qz.vertices_end(); ++v)
      {
        qv.push_back(*v);
      }
      return qv;
    }
    
  
    
    //template<class R>
    //Zonotope<R>::operator Polytope<typename Numeric::traits<R>::arithmetic_type> () const
    //{
    //  return Polytope<F>(this->vertices());
    //}
    
    //template<class R>
    //Zonotope<R>::operator Polyhedron<typename Numeric::traits<R>::arithmetic_type> () const
    //{
    //  return Polyhedron<F>(this->vertices());
    //}
    
    template<class R>
    Zonotope<R>::operator Polytope<Rational> () const
    {
      return Polytope<Rational>(Zonotope<Rational>(*this).vertices());
    }
    
    
    template<class R>
    Zonotope<R>::operator Polyhedron<Rational> () const
    {
      return Polyhedron<Rational>(Zonotope<Rational>(*this).vertices());
    }
    
  
    
    template<class R>
    tribool 
    Zonotope<R>::empty() const 
    {
      LinearAlgebra::LUMatrix<F> lu(this->generators());
      return indeterminate;
      //return lu.full_row_rank();
    }
      
    template<class R>
    Zonotope<typename Zonotope<R>::F> 
    Zonotope<R>::scale(const Zonotope<R>& z, const R& scale_factor) 
    {
      Geometry::Point<F> new_centre(z.centre().position_vector()*F(scale_factor));
      LinearAlgebra::Matrix<F> new_generators=z.generators()*F(scale_factor);
      return Geometry::Zonotope<F>(new_centre, new_generators);
    }

    template<class R>
    tribool 
    equal(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      return subset(A,B) && subset(B,A);
    }
    
    template<class R>
    Rectangle<R>
    Zonotope<R>::bounding_box() const
    {
      LinearAlgebra::Vector< Interval<R> > v(this->number_of_generators());
      for(size_type j=0; j!=this->number_of_generators(); ++j) {
        v[j]=Interval<R>(R(-1),1);
      }
      v=this->generators()*v;
      return Rectangle<R>(this->centre().position_vector()+v);
    }

  
    template<class R>
    ListSet<R,Zonotope>
    Zonotope<R>::subdivide() const 
    {
      //size_type n=this->dimension();
      ListSet<R,Geometry::Zonotope> result(this->dimension());
      
      R two=2;
      size_type d=(this->generators()).number_of_rows();
      size_type m=(this->generators()).number_of_columns();
      
      matrix_type new_generators(d,m);
      for(size_type i=0; i!=d; ++i) {
        for(size_type j=0; i!=m; ++j) {
          new_generators(i,j)=div_up(this->generators()(i,j),two);
        }
      }

      state_type first_centre=this->centre();
      for(size_type i=0; i<m; i++) {
        first_centre=sub_approx(first_centre,LinearAlgebra::Vector<R>(column(new_generators,i)));
      }
      
      for(unsigned long k=0; k!=1u<<m; ++k) {
        state_type new_centre=first_centre;
        for(size_type i=0; i<m; i++) {
          if(k & 1u<<i) {
            new_centre=add_approx(new_centre,this->generator(i));
          }
        }
        result.adjoin(Zonotope(new_centre,new_generators));
      }
      return result;
    }
    
    template<class R>
    ListSet<R,Zonotope>
    Zonotope<R>::divide() const 
    {
      size_type n=this->dimension();
      size_type m=(this->generators()).number_of_columns();
      R two=2;
      ListSet<R,Geometry::Zonotope> result(this->dimension());
      
      matrix_type new_generators=this->generators();
      
      R max_norm=0;
      size_type max_column=0;
      for(size_type j=0; j<m; j++) {
        R norm = LinearAlgebra::norm(LinearAlgebra::Vector<R>(column(new_generators,j)));
        if(norm>max_norm) {
          max_norm=norm;
          max_column=j;
        }
      }
      
      size_type j=max_column;
      for(size_type i=0; i!=n; ++i) {
        div_up(new_generators(i,j),two);
      }
      
      state_type new_centre=sub_approx(this->centre(),LinearAlgebra::Vector<R>(column(new_generators,j)));
      result.adjoin(Zonotope<R>(new_centre,new_generators));
      new_centre=sub_approx(new_centre,LinearAlgebra::Vector<R>(column(new_generators,j)));
      result.adjoin(Zonotope(new_centre,new_generators));

      return result;
    }
    

  

   
   
    
   
    
    template<class R>
    void 
    Zonotope<R>::minimize_generators(void) 
    {
      return;
    }
    

    template<class R>
    inline
    tribool 
    norm_grtr(const LinearAlgebra::Vector<R>& v1, const LinearAlgebra::Vector<R>& v2) 
    {
      return LinearAlgebra::norm(v1)>LinearAlgebra::norm(v2);
    }
    
    template<class R>
    void 
    Zonotope<R>::sort_generators(void)
    {
      std::vector< LinearAlgebra::Vector<R> > generator_vectors;
      for(size_type j=0; j!=this->number_of_generators(); ++j) {
        generator_vectors.push_back(column(this->_generators,j));
      }
      
      std::stable_sort(generator_vectors.begin(),generator_vectors.end(),&norm_grtr<R>);
      
      for(size_type j=0; j!=this->number_of_generators(); ++j) {
        for(size_type i=0; i!=this->dimension(); ++i) {
          this->_generators(i,j)=generator_vectors[j](i);
        }
      }
    }
  
    
    
/*
    template<class R>
    void 
    Zonotope<R>::compute_linear_inequalities(matrix_type& A, vector_type& b) const
    {
      throw std::runtime_error("Zonotope<R>::compute_linear_inequalities(matrix_type&, vector_type&) const not implemented");
      using namespace Ariadne::LinearAlgebra;
     
      const matrix_type &gen=this->_generators;
      matrix_type Space=trans(gen);
      R b2;
      size_type ngen=this->number_of_generators();
      
      //compute complanar hyperplanes
      array<size_type> col,row;
      // FIXME: Can't use lu_decompose for general R
      //A=lu_decompose(Space,col,row); 

      if (col.size()==row.size()) {
        A=matrix_type(2*ngen,this->dimension());
        b=vector_type(2*ngen);
        Space=matrix_type();
      } 
      else { 
        remove_null_columns(A,row,col);
        Space=compute_space(A,row,col);

        A=matrix_type(2*ngen+2*Space.number_of_rows(),
                        this->dimension());
        b=vector_type(2*ngen+2*Space.number_of_rows());
      }
      // TO IMPROVE: we new compute both (col_i)^T(cols_j) and 
      // (col_j)^T(cols_i)
      for (size_type i=0; i<ngen; i++) {
        b(2*i)=0.0;
        for (size_type j=0; j<ngen; j++) {
          //b(2*i)-=abs(this->principle_direction(i)*this->principle_direction(j));
          b2=0.0;
          for (size_type k=0; k< this->dimension(); k++) {
            b2+=gen(k,i)*gen(k,j);
          }
          b(2*i)-=abs(b2);
        }
        
        b2=0.0;
        for (size_type k=0; k< this->dimension(); k++) {
          b2+=gen(k,i)*this->centre()[k];
        }

        b(2*i+1)=b(2*i)+b2;
        b(2*i)=b(2*i)-b2;        
        
        for (size_type k=0; k< this->dimension(); k++) {
           A(2*i,k)=-gen(k,i);
           A(2*i+1,k)=gen(k,i);
        }
      }

      for (size_type i=0; i< Space.number_of_rows(); i++) {
        for (size_type j=0; j< this->dimension(); j++) {
          A(2*(i+ngen),j)=Space(i,j);
          A(2*(i+ngen)+1,j)=-Space(i,j);
        }
      }
    }
*/    

    /*!Set up LP problem to solve \f$c+Ge=p\f$; \f$-1<=e<=1\f$.
     * Change variables so that the problem becomes \f$Ge=p-c-G1;\ 0\leq e\leq2\f$.
     * Change sign of \f$ Ge=p-c-G1\f$ to make right-hand side positive.
     */
    template<class R>
    tribool 
    Zonotope<R>::contains(const Point<R>& pt) const 
    {
      //std::cerr << "Zonotope<R>::contains(const Point<R>& )" << std::endl;
      //typedef typename Numeric::traits<R,R>::arithmetic_type F;
      typedef Rational F;
      const Zonotope<R>& z=*this;
      check_dimension(z,pt,"Zonotope<R>::contains(Point<R>)");
      dimension_type d=z.dimension();
      dimension_type m=z.number_of_generators();
      
      F zero=0;
      F one=1;
      F two=2;
      
      const Geometry::Point<F> qc=z.centre();
      const Geometry::Point<F> qp=pt;
      const LinearAlgebra::Matrix<F> qG=z.generators();
      const LinearAlgebra::Vector<F> qo(m,one);
      const LinearAlgebra::Vector<F> zv(m,zero);
      const LinearAlgebra::Vector<F> tv(m,two);

      LinearAlgebra::Vector<F> qrhs=qp-qc+qG*qo;
      
      
      
      
      LinearAlgebra::Matrix<F> T(d+m+1,m+1);
      
      // Set up constraints e+se=2
      for(dimension_type j=0; j!=m; ++j) {
        T(j,j)=1;
        T(j,m)=2;
      }
      
      // Set up constraints Ge \pm ax = p-c+G1
      for(dimension_type i=0; i!=d; ++i) {
        if(qrhs(i)>=F(0)) {
          for(dimension_type j=0; j!=m; ++j) {
            T(m+i,j)=qG(i,j); 
          }
          T(m+i,m)=qrhs(i);
        } else {
          for(dimension_type j=0; j!=m; ++j) {
            T(m+i,j)=-qG(i,j); 
          }
          T(m+i,m)=-qrhs(i);
        }
      }
      
      // Set up cost function ax^T 1
      for(dimension_type i=0; i!=d; ++i) {
        for(dimension_type j=0; j!=m; ++j) {
          T(m+d,j)-=T(m+i,j);
        }
        T(m+d,m)-=T(m+i,m);
      }
      
      LinearAlgebra::LinearProgram<F> lp(T);
      //std::cerr << lp.tableau() << std::endl;
      tribool result=lp.is_feasible();
      //std::cerr << lp.tableau() << std::endl;
      return result;
    }

  
    
    /*!Set up linear program to solve 
     *   \f\[x=c+Ge;\ l\leq x\leq u;\ -1\leq e\leq1\f\].
     *
     * Change variables to normalize \f$x\f$ and \f$e\f$
     *   \f\[x'=x-l,\ e'=e+1;   x'-Ge' = c-G1-l;  0\leq x\leq u-l; \ 0\leq e\leq 2.\f\] 
     * 
     * Introduce slack variables sx and se, and artificial variables ax. Problem in standard form
     *   \f\[ \matrix{I&0\\0&I\\\pm I&\mp G}\matrix{x'\\e'} + \matrix{I&&\\&I&\\&&I}\matrix{sx\\se\\ax} = \matrix{u-l\\2\\\pm(c-G1-l)}
     * 
     */
    template<class R>
    tribool
    disjoint(const Zonotope<R>& z, const Rectangle<R>& r)
    {
      //std::cerr << "Zonotope<R>::disjoint(const Zonotope<R>&, const Rectangle<R>&)" << std::endl;
      //std::cerr << r << "\n" << z << std::endl;
      check_dimension(z,r,"disjoint(Zonotope<R>,Rectangle<R>)");
      dimension_type d=z.dimension();
      size_type m=z.number_of_generators();
    
      //This method has some problems. It fails when the parameters are
      //[1,17/16]x[19/16,5/4] and 
      //  Zonotope(
      //    centre=(1/2, 1/10)
      //    directions=[ 1,1/2; 1/2,3/5 ]
      //  )
      
      // Construct tableau for testing intersection of zonotope and rectangle
      // Rectangle  l<=x<=u
      // Zonotope  x==c+Ge,  -1<=e<=1
      // 
      // Translate x'=x-l,  e'=e+1
      //   x'+l==c+G(e'-1) ->  x'-Ge' == c-l-G1
      //   0<=x'<=u-l      ->  x' + sx' == u-l
      //   0<=e'<=2        ->  e' + se' == 2
      //  
      // Change sign of RHS of first equality if necessary
      // Introduce slack variables for last two inequalities
      typedef Rational F;
      LinearAlgebra::Matrix<F> T(2*d+m+1,d+m+1);

      const Geometry::Point<R>& l=r.lower_corner();
      const Geometry::Point<R>& u=r.upper_corner();
      const Geometry::Point<R>& c=z.centre();
      const LinearAlgebra::Matrix<R>& G=z.generators();

      const LinearAlgebra::Vector<F> qo(m,F(1));
      const LinearAlgebra::Vector<F> ql=l.position_vector();
      const LinearAlgebra::Vector<F> qu=u.position_vector();
      const LinearAlgebra::Vector<F> qc=c.position_vector();
      const LinearAlgebra::Matrix<F> qG=G;
      const LinearAlgebra::Vector<F> qrhs=qc-ql-qG*qo;
      
      // Set up constraints x+sx=u-l
      for(size_type i=0; i!=d; ++i) {
        T(i,i)=1;
        T(i,d+m)=qu(i)-ql(i);
      }
      
      // Set up constraints e+se=2
      for(size_type j=0; j!=m; ++j) {
        T(d+j,d+j)=1;
        T(d+j,d+m)=2;
      }
      
      // Set up constraints x-Ge \pm ax=c-l-G1 
      for(size_type i=0; i!=d; ++i) {
        if(qrhs(i)>=F(0)) {
          T(i+d+m,i)=1;
          for(size_type j=0; j!=m; ++j) {
            T(i+d+m,m+j)=-qG(i,j);
          }
          T(i+d+m,d+m)=qrhs(i);
        }
        else {
          T(i+d+m,i)=-1;
          for(size_type j=0; j!=m; ++j) {
            T(i+d+m,m+j)=qG(i,j);
          }
          T(i+d+m,d+m)=-qrhs(i);
        }
      } 
      
      // Set up cost function ax^T1
      for(size_type i=0; i!=d; ++i) {
        T(2*d+m,i) -= T(i+d+m,i);
        for(size_type j=0; j!=m; ++j) {
          T(2*d+m,d+j) -= T(i+d+m,d+j);
        }
        T(2*d+m,d+m) -= T(i+d+m,d+m);
      }
      
      //std::cerr << "T=" << T << std::endl;
      LinearAlgebra::LinearProgram<F> lp(T);
      tribool result=!lp.is_feasible();
      //std::cerr << "T=" << lp.tableau() << std::endl;
      //std::cerr << "disjoint(" << z << "," << r << ")=" << result << std::endl;
      return result;
    }

   
    template<class R>
    tribool
    disjoint(const Rectangle<R>& r, const Zonotope<R>& z)
    {
      //std::cerr << "disjoint(const Rectangle<R>&, const Zonotope<R>&)" << std::endl;
      return disjoint(z,r);
    }
    
  

    template<class R>
    tribool
    disjoint(const Zonotope<R>& z1, const Zonotope<R>& z2)
    {
      typedef Rational F;
      if(z1.dimension()!=z2.dimension()) {
        throw std::runtime_error("Invalid dimension in disjoint(const Zonotope<R>&, const Zonotope<R>&)");
      }
      
      dimension_type d=z1.dimension();
      size_type m1=z1.number_of_generators();
      size_type m2=z2.number_of_generators();
      
      LinearAlgebra::Matrix<F> T(m1+m2+d+1,m1+m2+1);

      const LinearAlgebra::Vector<F> qo1(m1,F(1));
      const LinearAlgebra::Vector<F> qo2(m2,F(1));

      Geometry::Point<F> qc1=z1.centre();
      LinearAlgebra::Matrix<F> qG1=z1.generators();
      Geometry::Point<F> qc2=z2.centre();
      LinearAlgebra::Matrix<F> qG2=z2.generators();
      LinearAlgebra::Vector<F> qrhs = qG1*qo1 - qG2*qo2 + (qc2 - qc1);
            
      // Set up constraints e1 + se1 = 2
      for(size_type j1=0; j1!=m1; ++j1) {
        T(j1,j1)=1;
        T(j1,m1+m2)=2;
      }
      
      // Set up constraints e2 + se2 = 2
      for(size_type j2=0; j2!=m2; ++j2) {
        T(m1+j2,m1+j2)=1;
        T(m1+j2,m1+m2)=2;
      }
      
      // Set up constraints G1*e1 - G2*e2 = (c2 - G2*1) - (c1 - G1*1)
      for(size_type i=0; i!=d; ++i) {
        if(qrhs(i)>=F(0)) {
          for(size_type j1=0; j1!=m1; ++j1) {
            T(m1+m2+i,j1)=qG1(i,j1);
          }
          for(size_type j2=0; j2!=m2; ++j2) {
            T(m1+m2+i,m1+j2)=qG2(i,j2);
          }
          T(m1+m2+i,m1+m2)=qrhs(i);
        }
        else {
          for(size_type j1=0; j1!=m1; ++j1) {
            T(m1+m2+i,j1)=-qG1(i,j1);
          }
          for(size_type j2=0; j2!=m2; ++j2) {
            T(m1+m2+i,m1+j2)=-qG2(i,j2);
          }
          T(m1+m2+i,m1+m2)=-qrhs(i);
        }
      } 
      
      // Set up cost function ax^T1
      for(size_type i=0; i!=d; ++i) {
        for(size_type j1=0; j1!=m1; ++j1) {
          T(m1+m2+d,j1) -= T(m1+m2+i,j1);
        }
        for(size_type j2=0; j2!=m2; ++j2) {
          T(m1+m2+d,m1+j2) -= T(m1+m2+i,m1+j2);
        }
        T(m1+m2+d,m1+m2) -= T(m1+m2+i,m1+m2);
      }
      
      LinearAlgebra::LinearProgram<F> lp(T);
      
      tribool result=!lp.is_feasible();
      
      //std::cerr << "disjoint(" << z1 << "," << z2 << ")=" << result << std::endl;
      return result;
}
    
    


    template<class R>
    tribool 
    subset(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      throw std::runtime_error("subset(const Rectangle<R>& A, const Zonotope<R>& B) not implemented");
      typedef Rational F;
      return Geometry::subset(Rectangle<F>(A),B.operator Polyhedron<F>());
    }

    template<class R>
    tribool 
    subset(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Geometry::subset(A.bounding_box(),B);
    }

    template<class R>
    tribool 
    subset(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      throw std::runtime_error("subset(const Zonotope<R>& A, const Zonotope<R>& B) not implemented");
      typedef Rational F;
      return Geometry::subset(A.operator Polyhedron<F>(),B.operator Polyhedron<F>());
    }

  
    

    template<class R> 
    Zonotope<typename Numeric::traits<R>::arithmetic_type>
    minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
      
      if (A.dimension()!=B.dimension()) {
          throw std::domain_error(
              "minkowski_sum: the two zonotopes have different dimension.");
      }
      
      Geometry::Point<F> new_centre=Geometry::minkowski_sum(A.centre(),B.centre());
      LinearAlgebra::Matrix<R> new_generators=LinearAlgebra::concatenate_columns(A.generators(),B.generators());
      return Zonotope<F>(new_centre,new_generators);
    }
  
 
    
    template<class R> 
    Zonotope<typename Numeric::traits<R>::arithmetic_type> 
    minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;
     
      if (A.dimension()!=B.dimension()) {
          throw std::domain_error(
              "minkowski_difference: the two zonotopes have different dimension.");
      }
      
      return Zonotope<F>(Geometry::minkowski_difference(A.centre(),B.centre()),
                         LinearAlgebra::concatenate_columns(A.generators(),B.generators()));
    }


    
    template<class R> 
    Rectangle<R> 
    Zonotope< Interval<R> >::bounding_box() const
    {
      LinearAlgebra::Vector< Interval<R> > e(this->number_of_generators(),Interval<R>(-1,1));
      
      Point< Interval<R> > b=this->centre()+this->generators()*e;
      Rectangle<R> result(this->dimension());
      for(size_type i=0; i!=this->dimension(); ++i) {
        result.set_lower_bound(i,b[i].lower());
        result.set_lower_bound(i,b[i].upper());
      }
      return result;
    }
    
    template<class R> 
    Zonotope<R> 
    over_approximation(const Zonotope< Interval<R> >& iz)
    {
      return orthogonal_over_approximation(iz);
    }
    
    template<class R> 
    Zonotope< Interval<R> > 
    box_over_approximation(const Zonotope< Interval<R> >& iz)
    {
      LinearAlgebra::Vector< Interval<R> > c=iz.centre();
      LinearAlgebra::Matrix< Interval<R> > G=iz.generators();
      
      R x,r;
      for(dimension_type i=0; i!=G.number_of_rows(); ++i) {
        for(size_type j=0; j!=G.number_of_columns(); ++j) {
          x=approximate_value(G(i,j));
          r=error_bound(G(i,j));
          c(i)=c(i)+Interval<R>(-r,r);
          G(i,j)=x;
        }
      }
      return Zonotope< Interval<R> >(c,G);
    }
    
  

/*    
    template<class R>
    Zonotope<R>
    Zonotope<R>::over_approximation(const Geometry::Point< Interval<R> > &c, const LinearAlgebra::Matrix< Interval<R> >& A)
    {
      return Zonotope<R>(c+LinearAlgebra::radius_row_sum(A),LinearAlgebra::approximate_value(A));
    }
  */

    template<class R>
    std::ostream&
    Zonotope<R>::write(std::ostream& os) const 
    {
      const Zonotope<R>& z=*this;
      if(z.dimension() > 0) {
        if (z.empty()) {
          os << "Zonotope( )" << std::endl;
        } else {
          os << "Zonotope( centre=" << z.centre();
          os << ", directions=" << z.generators();
          os << ") ";
        } 
      }
      return os;
    }
    
    template<class R>
    std::ostream&
    Zonotope< Interval<R> >::write(std::ostream& os) const 
    {
      const Zonotope< Interval<R> >& z=*this;
      if(z.dimension() > 0) {
        os << "Zonotope( centre=" << z.centre();
        os << " directions=" << z.generators();
        os << ") ";
      }
      return os;
    }
    
    template<class R>
    std::istream& 
    Zonotope<R>::read(std::istream& is)
    {
      throw std::domain_error("Not implemented");
    }
      
  

  }
}
