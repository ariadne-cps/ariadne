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

#include "zonotope.h"

#include "../base/array.h"
#include "../linear_algebra/vector.h"

#include "../geometry/list_set.h"
#include "../geometry/lattice_set.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Geometry {
    
    template<typename R>
    Zonotope<R>::Zonotope(const Rectangle<R>& r)
      : _central_block(r.dimension()), _generators(r.dimension(),r.dimension())
    {
      for(size_type i=0; i!=dimension(); ++i) {
        this->_central_block[i] = r[i].centre();
        this->_generators(i,i) = r[i].radius();
      }

      this->_generators=remove_null_columns_but_one(this->_generators);
    }

    template<typename R>
    Zonotope<R>::Zonotope(const Parallelotope<R>& p)
      : _central_block(p.centre()), _generators(p.generators())
    {
      this->minimize_generators();
    }

    template<typename R>
    Zonotope<R>::Zonotope(const std::string& s)
      : _central_block(), _generators()
    {
      std::stringstream ss(s);
      ss >> *this;
    }
         
    
    template<typename R>
    bool 
    Zonotope<R>::operator==(const Zonotope<real_type>& A) const
    {
      using namespace LinearAlgebra;
            
      const matrix_type &this_gen=this->_generators;
      const matrix_type &A_gen=A._generators;
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
    
    template <typename R>
    Rectangle<R>
    Zonotope<R>::bounding_box() const
    {
      LinearAlgebra::IntervalVector<R> e(this->number_of_generators());
      for(size_type j=0; j!=this->number_of_generators(); ++j) {
        e[j]=Interval<R>(-1,1);
      }
      return this->central_block()+this->generators()*e;
    }


    template <typename R>
    ListSet<R,Zonotope>
    Zonotope<R>::subdivide() const 
    {
      // FIXME: What to do about central block?
      //size_type n=this->dimension();
      size_type m=(this->generators()).size2();
      ListSet<R,Geometry::Zonotope> result(this->dimension());
      matrix_type new_generators=this->generators()/2;
      
      state_type first_centre=this->centre();
      for(size_type i=0; i<m; i++) {
        first_centre=first_centre-LinearAlgebra::Vector<R>((this->generator(i))/2);
      }
      
      array<index_type> lower(m,0);
      array<index_type> upper(m,2);
      array<index_type> finish(m,0);
      finish[m-1]=2;
      lattice_iterator end(finish,lower,upper);

      for(lattice_iterator iter(lower,lower,upper); iter!=end; ++iter) {
        array<index_type> ary=*iter;
        state_type new_centre=first_centre;
        for(size_type i=0; i<m; i++) {
          if(ary[i]==1) {
            new_centre=new_centre+this->generator(i);
          }
        }
        result.adjoin(Zonotope(new_centre,new_generators));
      }
      return result;
    }
    
    template <typename R>
    ListSet<R,Zonotope>
    Zonotope<R>::divide() const 
    {
      size_type n=this->dimension();
      size_type m=(this->generators()).size2();
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
        new_generators(i,j)/=2;
      }
      
      state_type new_centre=this->centre()-LinearAlgebra::Vector<R>(column(new_generators,j)/2);
      result.adjoin(Zonotope<R>(new_centre,new_generators));
      new_centre=new_centre+LinearAlgebra::Vector<R>(column(new_generators,j));
      result.adjoin(Zonotope(new_centre,new_generators));

      return result;
    }
    


    template<typename R>
    bool 
    Zonotope<R>::contains(const state_type& point) const 
    {
      return Polyhedron<Rational>(*this).contains(Point<Rational>(point));
    }

    
    template<typename R>
    bool 
    Zonotope<R>::interior_contains(const state_type& point) const 
    {
      // FIXME: Own routine
      Point<Rational> rational_point=convert_to< Point<Rational> >(point);
      return Polyhedron<Rational>(*this).interior_contains(rational_point);
    }


    template<typename R>
    bool 
    Zonotope<R>::superset(const Rectangle<R>& rect) const 
    {
      if (this->empty()) {
        return false;
      }
      for (size_type i=0; i< (size_type)(1<<rect.dimension()); i++) {
        if (!this->contains(rect.vertex(i))) {
          return false;
        }
      }

      return true;
    }
     


    template <typename R>
    std::vector< Point<Rational> > 
    Zonotope<R>::vertices() const
    {
      return Polyhedron<Rational>(*this).vertices();
    }

    template <typename R>
    std::vector< Point<R> > 
    Zonotope<R>::approximate_vertices() const
    {
      return Polyhedron<R>(this->centre(),this->_extended_generators()).vertices();
    }

    template <typename R>
    LinearAlgebra::Matrix<R>  
    Zonotope<R>::_extended_generators() const
    {
      size_type number_of_extra_generators=0;
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        if(this->_central_block[i].lower()==this->_central_block[i].upper()) {
          ++number_of_extra_generators;
        }
      }
      if(number_of_extra_generators==0) {
        return this->_generators;
      }
      else {
        LinearAlgebra::Matrix<R> result(this->dimension(),
            this->number_of_generators()+number_of_extra_generators);
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          for(size_type j=0; j!=this->number_of_generators(); ++j) {
            result(i,j)=this->_generators(i,j);
          }
        }
        size_type j=this->number_of_generators();
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          R radius=this->_central_block[i].radius();
          if(radius!=0) {
            result(i,j)=radius;
            ++j;
          }
        }
        return result;
      }
    }
      
    template <typename R>
    Zonotope<R>::operator Polyhedron<Rational>() const 
    {
      Point<Rational> c(this->centre());
      LinearAlgebra::Matrix<Rational> g(this->_extended_generators());
      return Polyhedron<Rational>(c,g);
      //return Polyhedron<Rational>(this->_possible_vertices());
    }
    
    
    
    template<typename R>
    Parallelotope<R>
    Zonotope<R>::over_approximating_parallelotope() const
    {
      typedef typename numerical_traits<R>::field_extension_type F;
      const Zonotope<R>& z=*this;
      
      dimension_type n=z.dimension();
      LinearAlgebra::Matrix<R> A(n,n);
      for(dimension_type i=0; i!=n; ++i) {
        for(dimension_type j=0; j!=n; ++j) {
          A(i,j)=z.generators()(i,j);
        }
      }
      LinearAlgebra::Matrix<F> B=LinearAlgebra::inverse(A);
      const Point<R>& c=z.centre();
      
      Parallelotope<R> p(c,A);
      while(!subset(z,p)) {
        A*=2;
        p=Parallelotope<R>(c,A);
      }
      return p;
    }        
    
    
    
    template<typename R>
    void 
    Zonotope<R>::minimize_generators(void) 
    {
      using namespace LinearAlgebra;
           
      size_type i,j,i2,j2;
      matrix_type &gen=this->_generators;
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
        matrix_type new_gen(rows,min_cols);

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
    norm_grtr(const LinearAlgebra::Vector<R>& v1, const LinearAlgebra::Vector<R>& v2) 
    {
      return LinearAlgebra::norm(v1)>LinearAlgebra::norm(v2);
    }
    
    template<typename R>
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
    
    
    
    template<typename R>
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

        A=matrix_type(2*ngen+2*Space.size1(),
                        this->dimension());
        b=vector_type(2*ngen+2*Space.size1());
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

      for (size_type i=0; i< Space.size1(); i++) {
        for (size_type j=0; j< this->dimension(); j++) {
          A(2*(i+ngen),j)=Space(i,j);
          A(2*(i+ngen)+1,j)=-Space(i,j);
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
      
      const Matrix<R> &A_gen=A.generators();
      const Matrix<R> &B_gen=B.generators();
      
      dimension_type col_A=(A_gen).size2(); 
      dimension_type col_B=(B_gen).size2(); 
      
      Matrix<R> gen(A.dimension(), col_A+col_B);

      for (size_type i=0; i!=col_A; ++i) {
         for (size_type j=0; j!=A.dimension(); ++j) {
            gen(j,i)=A_gen(j,i);
         }
      }

      for (size_type i=0; i!=col_B; ++i) {
        for (size_type j=0; j!=A.dimension(); ++j) {
          gen(j,i+col_A)=B_gen(j,i);
        }
      }
      
      return Zonotope<R>(A.centre()+(B.centre()).position_vector(),gen);
    }
   
    template<typename R> 
    Zonotope<R> 
    minkowski_sum(const Rectangle<R>& A, const Zonotope<R>& B)
    {
      return minkowski_sum(Zonotope<R>(A),B);
    }
    
    template<typename R> 
    Zonotope<R> 
    minkowski_sum(const Zonotope<R>& A, const Rectangle<R>& B)
    {
      return minkowski_sum(A,Zonotope<R>(B));
    }
    
    
    template<typename R> 
    Zonotope<R> minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      using namespace Ariadne::LinearAlgebra;
     
      if (A.dimension()!=B.dimension()) {
          throw std::domain_error(
              "minkowski_difference: the two zonotopes have different dimension.");
      }
      
      const Matrix<R> &A_gen=A.generators();
      const Matrix<R> &B_gen=B.generators();
      
      dimension_type col_A=(A_gen).size2(); 
      dimension_type col_B=(B_gen).size2(); 
      
      Matrix<R> gen(A.dimension(), col_A+col_B);

      for (size_type i=0; i!=col_A; ++i) {
         for (size_type j=0; j!=A.dimension(); ++j) {
          gen(j,i)=A_gen(j,i);
        }
      }

      
      for (size_type i=0; i!=col_B; ++i) {
        for (size_type j=0; j!=A.dimension(); ++j) {
          gen(j,i+col_A)=B_gen(j,i);
         }
      }
      
      return Zonotope<R>(A.centre()-(B.centre()).position_vector(), 
                         remove_null_columns_but_one(gen));
    }

    template<typename R> 
    Zonotope<R> 
    minkowski_difference(const Rectangle<R>& A, const Zonotope<R>& B)
    {
      return minkowski_difference(Zonotope<R>(A),B);
    }
    
    template<typename R> 
    Zonotope<R> 
    minkowski_difference(const Zonotope<R>& A, const Rectangle<R>& B)
    {
      return minkowski_difference(A,Zonotope<R>(B));
    }
    
    
    
    template<typename R>
    Zonotope<R>
    operator+(const Rectangle<R>& r, const LinearAlgebra::TransformationSystem<R>& v)
    {
      //std::cerr << "operator+(const Rectangle<R>& r, const TransformationSystem<R>& v)" << std::endl;
      dimension_type n=r.dimension();
      LinearAlgebra::Matrix<R> r_generators(n,n);
      for(size_type i=0; i!=n; ++i) {
        r_generators(i,i)=(r[i].upper()-r[i].lower())/2;
      }
      //std::cerr << LinearAlgebra::concatenate_columns(v.generators(),r_generators) << std::endl;
      return Zonotope<R>(r.centre()+v.centre(),LinearAlgebra::concatenate_columns(v.generators(),r_generators));
    }
    
    template<typename R>
    Zonotope<R>
    operator+(const Zonotope<R>& z, const LinearAlgebra::TransformationSystem<R>& v)
    {
      return Zonotope<R>(z.centre()+v.centre(),LinearAlgebra::concatenate_columns(z.generators(),v.generators()));
    }
    
    
    
    template<typename R>
    Zonotope<R>
    Zonotope<R>::over_approximation(const Rectangle<R> &c, const LinearAlgebra::IntervalMatrix<R>& A)
    {
      LinearAlgebra::Matrix<R> Amid=A.centre();
      return Zonotope<R>(c+A.radius_row_sum(),A.centre());
    }
    


    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Zonotope<R>& z) 
    {
      if(z.dimension() > 0) {
        if (!(z.empty())) {
          os << "Zonotope(\n  centre=" << z.centre();
          os << "\n  directions=" << z.generators();
          os << "\n) ";
        } else {
          os << "Zonotope( Empty )" << std::endl;
        }
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
