/***************************************************************************
 *            zonotope.h
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

#include "../linear_algebra/vector.h"
#include "../geometry/zonotope.h"


namespace Ariadne {
  namespace Geometry {

    
    
    template<typename R>
    bool 
    Zonotope<R>::contains(const State& point) const 
    {
      if (point.dimension()!=this->dimension()) {
        throw std::domain_error("This object and parameter have different space dimensions");
      }  
      
      if (this->empty()) { return false; }
    
      Matrix A(0,0);
      Vector b(0);
              
      this->compute_linear_inequalities(A,b);

      Vector result(A*point.position_vector()-b);

      for (size_type i=0; i<result.size(); i++) {
         if (result(i)<0) return false;
      }
      
      return true;
    }
      
    template<typename R>
    bool 
    Zonotope<R>::interior_contains(const State& point) const 
    {
       using namespace Ariadne::LinearAlgebra;
      
      if (point.dimension()!=this->dimension()) {
       throw std::domain_error("This object and parameter have different space dimensions");
      } 
    
      if (this->empty_interior()) { return false; }
  
      Matrix A(0,0);
      Vector b(0);
              
      this->compute_linear_inequalities(A,b);

      Vector result(A*point.position_vector()-b);
      
      for (size_type i=0; i<result.size(); i++) {
         if (result(i)<=0) return false;
      }
      
      return true;
    }
      
    template<typename R>
    bool 
    Zonotope<R>::operator==(const Zonotope<Real>& A) const
    {
      using namespace LinearAlgebra;
            
      const Matrix &this_gen=this->_generators;
      const Matrix &A_gen=A._generators;
      size_type directions=this->number_of_generators();

      if (!have_same_dimensions(this_gen,A_gen)) {
        return false;
      }
      array<bool> not_found(directions,true);
      bool searching_for_equiv;

      size_type j2;
      for (size_type j=0; j< directions; j++) {
        j2=j; 
        searching_for_equiv=true;
        while ((j2< directions)&&(searching_for_equiv)) {
          if (not_found[j2]) {
            if (equivalent_columns(this_gen, j, A_gen, j2)) { 
              searching_for_equiv=false;
              not_found[j2]=false;
            }
          }
          j2++;
        }

        if (searching_for_equiv) {
          return false;
        }
      }

      return true;
    }
    
    template<typename R>
    Rectangle<R> 
    Zonotope<R>::bounding_box() const 
    {
      using namespace Ariadne::LinearAlgebra;
      Vector offset(this->dimension());
      
      for(size_type i=0; i!=this->dimension(); ++i) {
        for(size_type j=0; j!=this->_generators.size2(); ++j) {
           offset(i) += abs(this->_generators(i,j));
        }
      }
      return Rectangle<R>(this->centre()-offset, this->centre()-offset);
    }

    template<typename R>
    void 
    Zonotope<R>::minimize_generators(void) 
    {
      using namespace LinearAlgebra;
           
      size_type i,j,i2,j2;
      Matrix &gen=this->_generators;
      gen=remove_null_columns_but_one(gen);
      size_type rows=gen.size1();
      
      // if the first row is null the zonotope is a point and it is already
      // minimized
      if (find_first_not_null_in_col(gen,0)==rows) {
        return;
      }
      
      size_type cols=this->number_of_generators();
      size_type min_cols=cols;
      R coef,coef2;
     
      array<size_type> dependences(cols,cols);
      array<bool> same_sign(cols,true);
      
      for (j=0; j<cols; j++) {
        i=find_first_not_null_in_col(gen,j);
             
        // if the first row is null the zonotope is a point and it is already
        // minimized
        assert(i!=rows);

        coef=gen(i,j);
        
        for (j2=j+1; j2<cols; j2++) {
          i2=find_first_not_null_in_col(gen,j2);

          if (i==i2) {
            coef2=gen(i2,j2);

            //check whever the j-th and the j2-th columns are linear depended
            //or not.
            while ((i2<rows)&&(gen(i2,j2)*coef==gen(i2,j)*coef2))
              i2++;

            // if they are
            if (i2==rows) {
              if (dependences[j]==cols) {
                dependences[j2]=j;
                if (coef*coef2<0) same_sign[j2]=false;
                min_cols--;
              } else {
                dependences[j2]=dependences[j];
                if ((coef*coef2>0)^same_sign[j]) same_sign[j2]=false;
              }
            }
          }
        }
      }
      
      if (min_cols!= cols) {
        Matrix new_gen(rows,min_cols);

        j2=0;
        for (j=0; j< cols; j++) {
          if (dependences[j]==cols) {
            dependences[j]=j2;
            
            for (i=0; i<rows; i++)
              new_gen(i,j2)=gen(i,j);
            
            j2++;
          } else {
             
             if (same_sign[j]) {
               for (i=0; i<rows; i++)
                 new_gen(i,dependences[j])+=gen(i,j);
             } else {
               for (i=0; i<rows; i++)
                 new_gen(i, dependences[j])-=gen(i,j);
             }
          }
        }

        gen=new_gen;
      }
      sort_generators();
    }
    
    
    template<typename R>
    inline
    bool 
    norm_grtr(const LinearAlgebra::vector<R>& v1, const LinearAlgebra::vector<R>& v2) 
    {
      return LinearAlgebra::norm(v1)>LinearAlgebra::norm(v2);
    }
    
    
    template<typename R>
    void 
    Zonotope<R>::sort_generators(void)
    {
      std::vector< LinearAlgebra::vector<R> > generator_vectors;
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
    
    template<typename R> 
    Zonotope<R> 
    minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      using namespace Ariadne::LinearAlgebra;
     
      if (A.dimension()!=B.dimension()) {
          throw std::domain_error(
              "minkowski_sum: the two zonotopes have different dimension.");
      }
      
      dimension_type col_A=A._generators.size2(); 
      dimension_type col_B=B._generators.size2(); 
      
      matrix<R> gen(A.dimension(), col_A+col_B);

      for (size_type i=0; i!=col_A; ++i) {
         for (size_type j=0; j!=A.dimension(); ++j) {
            gen(j,i)=A._generators(j,i);
         }
      }

      
      for (size_type i=0; i!=col_B; ++i) {
        for (size_type j=0; j!=A.dimension(); ++j) {
          gen(j,i+col_A)=B._generators(j,i);
        }
      }
      
      return Zonotope<R>(A._centre+(B._centre).position_vector(),gen);
    }
   
    template<typename R> 
    Zonotope<R> minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      using namespace Ariadne::LinearAlgebra;
     
      if (A.dimension()!=B.dimension()) {
          throw std::domain_error(
              "minkowski_difference: the two zonotopes have different dimension.");
      }
      
      dimension_type col_A=A._generators.size2();
      dimension_type col_B=B._generators.size2(); 
      
      matrix<R> gen(A.dimension(), col_A+col_B);

      for (size_type i=0; i!=col_A; ++i) {
         for (size_type j=0; j!=A.dimension(); ++j) {
          gen(j,i)=A._generators(j,i);
        }
      }

      
      for (size_type i=0; i!=col_B; ++i) {
        for (size_type j=0; j!=A.dimension(); ++j) {
          gen(j,i+col_A)=B._generators(j,i);
         }
      }
      
      return Zonotope<R>(A._centre-(B._centre).position_vector(), 
                         remove_null_columns_but_one(gen));
    }
    template<typename R>
    void 
    Zonotope<R>::compute_linear_inequalities(Matrix& A, Vector& b) const
    {
      using namespace Ariadne::LinearAlgebra;
     
      const Matrix &gen=this->_generators;
      Matrix Space=trans(gen);
      R b2;
      size_type ngen=this->number_of_generators();
      
      //compute complanar hyperplanes
      array<size_type> col,row;
      A=lu_decompose(Space,col,row); 

      if (col.size()==row.size()) {
        A=Matrix(2*ngen,this->dimension());
        b=Vector(2*ngen);
        Space=Matrix();
      } 
      else { 
        remove_null_columns(A,row,col);
        Space=compute_space(A,row,col);

        A=Matrix(2*ngen+2*Space.size1(),
                        this->dimension());
        b=Vector(2*ngen+2*Space.size1());
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
          b2+=gen(k,i)*this->_centre[k];
        }

        b(2*i+1)=b(2*i)+b2;
        b(2*i)=b(2*i)-b2;        
        
        for (size_type k=0; k< this->dimension(); k++) {
           A(2*i,k)=-gen(k,i);
           A(2*i+1,k)=gen(k,i);
        }
      }

      for (size_type i=0; i< Space.size1(); i++) {
        for (size_type j=0; j< this->dimension(); j++) {
          A(2*(i+ngen),j)=Space(i,j);
          A(2*(i+ngen)+1,j)=-Space(i,j);
        }
      }
    }
      
    
    template <typename R>
    Zonotope<R>::operator Polyhedron<R>() const 
    {
      using namespace Ariadne::LinearAlgebra;
      
      matrix<R> A(0,0);
      vector<R> b(0);
      
      this->compute_linear_inequalities(A,b);
      
      return Polyhedron<R>(A,b);
    }
    
    
    
    
    template<typename R>
    Geometry::Parallelotope<R>
    over_approximating_parallelotope(const Geometry::Zonotope<R>& z)
    {
      std::cerr << "over_approximating_parallelotope(const Geometry::Zonotope<R>& z)" << std::endl;
      dimension_type n=z.dimension();
      LinearAlgebra::matrix<R> A(n,n);
      for(dimension_type i=0; i!=n; ++i) {
        for(dimension_type j=0; j!=n; ++j) {
          A(i,j)=z.generators()(i,j);
        }
      }
      LinearAlgebra::matrix<R> B=LinearAlgebra::inverse(A);
      const Point<R>& c=z.centre();
      
      Parallelotope<R> p(c,A);
      while(!subset(z,p)) {
        std::cerr << "A=" << A << std::endl;
        A*=2;
        p=Parallelotope<R>(c,A);
      }
      std::cerr << "A=" << A << std::endl;
      return p;
    }        
    
    template<typename R>
    Zonotope<R>
    operator+(const Rectangle<R>& r, const LinearAlgebra::zonotopic_vector<R>& v)
    {
      std::cerr << "operator+(const Rectangle<R>& r, const zonotopic_vector<R>& v)" << std::endl;
      dimension_type n=r.dimension();
      LinearAlgebra::matrix<R> r_generators(n,n);
      for(size_type i=0; i!=n; ++i) {
        r_generators(i,i)=(r[i].upper()-r[i].lower())/2;
      }
      std::cerr << LinearAlgebra::concatenate_columns(v.generators(),r_generators) << std::endl;
      return Zonotope<R>(r.centre()+v.centre(),LinearAlgebra::concatenate_columns(v.generators(),r_generators));
    }
    
    template<typename R>
    Zonotope<R>
    operator+(const Zonotope<R>& z, const LinearAlgebra::zonotopic_vector<R>& v)
    {
      return Zonotope<R>(z.centre()+v.centre(),LinearAlgebra::concatenate_columns(z.generators(),v.generators()));
    }
    
    
    
    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Zonotope<R>& z) 
    {
      if(z.dimension() > 0) {
        os << "Zonotope(\n  centre=" << z.centre();
        os << "\n  directions=" << z.generators();
        os << "\n) ";
      }

      return os;
    }
    
    template <typename R>
    std::istream& 
    operator>>(std::istream& is, Zonotope<R>& z)
    {
      throw std::domain_error("Not implemented");
    }
      
    

  }
}
