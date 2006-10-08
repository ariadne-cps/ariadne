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
#include "../numeric/conversion.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../combinatoric/lattice_set.h"

#include "../geometry/ppl_polyhedron.h"
#include "../geometry/point.h"
#include "../geometry/point_list.h"
#include "../geometry/rectangle.h"
#include "../geometry/polyhedron.h"
#include "../geometry/parallelotope.h"
#include "../geometry/list_set.h"

// include for column operations
#include "../linear_algebra/matrix.tpl" 

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

      //this->_generators=remove_null_columns_but_one(this->_generators);
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
    Zonotope<R>::empty_interior() const 
    {
      return !LinearAlgebra::independent_rows(this->_generators);
    }
      
    template<typename R>
    Zonotope<R> 
    Zonotope<R>::scale(const Zonotope<R>& z, const R& scale_factor) 
    {
      Geometry::Rectangle<R> new_central_block=Geometry::scale(z.central_block(),scale_factor);
      LinearAlgebra::Matrix< Interval<R> > new_interval_generators=z.generators()*Interval<R>(scale_factor);
      new_central_block=new_central_block+LinearAlgebra::radius_row_sum(new_interval_generators);
      LinearAlgebra::Matrix<R> new_generators=LinearAlgebra::centre(new_interval_generators);
      return Geometry::Zonotope<R>(new_central_block, new_generators);
    }

    template<typename R>
    bool 
    Zonotope<R>::equal(const Zonotope<real_type>& A, const Zonotope<real_type>& B)
    {
      using namespace LinearAlgebra;
            
      const matrix_type &A_gen=A._generators;
      const matrix_type &B_gen=B._generators;
      size_type directions=A.number_of_generators();

      if (!have_same_dimensions(A_gen,B_gen)) {
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
            if (equivalent_columns(A_gen, j, B_gen, j2)) { 
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
      LinearAlgebra::Vector< Interval<R> > e(this->number_of_generators());
      for(size_type j=0; j!=this->number_of_generators(); ++j) {
        e[j]=Interval<R>(-1,1);
      }
      return this->central_block()+(this->generators()*e);
    }


    template <typename R>
    ListSet<R,Zonotope>
    Zonotope<R>::subdivide() const 
    {
      // FIXME: What to do about central block?
      //size_type n=this->dimension();
      size_type m=(this->generators()).number_of_columns();
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
      size_type m=(this->generators()).number_of_columns();
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
      LinearAlgebra::Vector<Rational> v=point.position_vector();
      return Geometry::subset(ppl_polyhedron(v),
                              Parma_Polyhedra_Library::C_Polyhedron(*this));
    }

    
    template<typename R>
    bool 
    Zonotope<R>::interior_contains(const state_type& point) const 
    {
      // FIXME: Own routine
      LinearAlgebra::Vector<Rational> v=point.position_vector();
      return Geometry::inner_subset(ppl_polyhedron(v),
                                    Parma_Polyhedra_Library::C_Polyhedron(*this));
    }




    template <typename R>
    PointList<Rational>
    Zonotope<R>::vertices() const
    {
      return PointList<Rational>(Geometry::generators(Parma_Polyhedra_Library::C_Polyhedron(*this)));
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
    Zonotope<R>::operator Parma_Polyhedra_Library::C_Polyhedron() const 
    {
      LinearAlgebra::Vector<Rational> c(this->centre().position_vector());
      LinearAlgebra::Matrix<Rational> g(this->_extended_generators());
      return ppl_polyhedron(c,g);
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
      while(!Geometry::subset(z,p)) {
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
      size_type rows=gen.number_of_rows();
      
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
    

      /*
    template<typename R>
    bool
    Zonotope<R>::disjoint(const Rectangle<R>& r) const
    {
      const Zonotope<R>& z=*this;
      assert(z.dimension()==r.dimension());
      dimension_type n=z.dimension();
      dimension_type m=z.number_of_generators();
    
      This method has some problems. It fails when the parameters are
      [1,17/16]x[19/16,5/4] and 
        Zonotope(
          centre=(1/2, 1/10)
          directions=[ 1,1/2; 1/2,3/5 ]
        )
      // Construct tableau for testing intersection of rectangle and point
      // Rectangle  a<=x<=b
      // Parallelotope  x==c+Ae,  -1<=e<=1
      // 
      // Translate x'=x-a,  e'=e+1
      //   0<=x'<=b-a
      //   0<=e'<=2
      //   x'+a==c+A(e'-1) ->  x'-Ae' == c-a-A1
      //  
      // Introduce slack variables for first two inequalities
      // Introduce auxiliary variables for last equality, changing sign of RHS if necessary
      // 
      // Need to minimise sum of auxiliary variables -> add sum of last rows 
      // to get value function.
      LinearAlgebra::Matrix<Rational> T(3*n+1,2*m+1);

      const Geometry::Point<R>& a=r.lower_corner();
      const Geometry::Point<R>& b=r.upper_corner();
      const Geometry::Point<R>& c=z.centre();
      const LinearAlgebra::Matrix<R>& A=z.generators();

      for(size_type i=0; i!=n; ++i) {
        T(i,i)=1;
        T(i,2*m)=Rational(b[i])-Rational(a[i]);
        T(n+i,m+i)=1;
        T(n+i,2*m)=2;

        // Compute rhs = c[i]-a[i]-(A*1)[i]
        Rational rhs=Rational(c[i]) - Rational(a[i]);
        for(size_type j=0; j!=m; ++j) {
          rhs-=Rational(A(i,j));
        }
        
        if(rhs>=0) {
          T(2*n+i,i)=1;
          for(size_type j=0; j!=n; ++j) {
            T(2*n+i,m+j)=-A(i,j);
          }
          T(2*n+i,2*m)=rhs;
        }
        else {
          T(2*n+i,i)=-1;
          for(size_type j=0; j!=m; ++j) {
            T(2*n+i,m+j)=A(i,j);
          }
          T(2*n+i,2*m)=-rhs;
        }
        for(size_type j=0; j!=2*m; ++j) {
          T(3*n,j)-=T(2*n+i,j);
        }
        T(3*n,2*m)-=T(2*n+i,2*m);
      }
      
      LinearAlgebra::LinearProgram<Rational> lp(T);
      
      bool result=(lp.optimal_value()!=0);
      
      if(result!=Geometry::disjoint(Polyhedron<R>(*this), Polyhedron<R>(r))) {
        std::cerr << "Incorrect result for \n  " << r << "\nand\n" << *this << "\n";
        std::cerr << T << "\n" << lp.tableau() << "\n";
        std::cerr << convert_to<double>(lp.tableau()(3*n,2*m)) << "\n";
        assert(false);
      }
      return Geometry::disjoint(Polyhedron<Rational>(*this), Polyhedron<Rational>(r));
    }
    */
   
    template <typename R>
    bool 
    Zonotope<R>::disjoint(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return !minkowski_difference(A,B).contains(Point<R>(A.dimension()));
    }
    
    template <typename R>
    bool 
    Zonotope<R>::disjoint(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return disjoint(A,Zonotope<R>(B));
    }
    
    
    
    template <typename R>
    bool 
    Zonotope<R>::interiors_intersect(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return !minkowski_difference(A,B).interior_contains(Point<R>(A.dimension()));
    }
   
    template <typename R>
    bool 
    Zonotope<R>::interiors_intersect(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return interiors_intersect(A,Zonotope<R>(B));
    }
    
    
    
    template <typename R>
    bool 
    Zonotope<R>::inner_subset(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return Geometry::inner_subset(Parma_Polyhedra_Library::C_Polyhedron(A),
                                    Parma_Polyhedra_Library::C_Polyhedron(B));
    }

    template <typename R>
    bool 
    Zonotope<R>::inner_subset(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Geometry::inner_subset(Parma_Polyhedra_Library::C_Polyhedron(A),
                                    Parma_Polyhedra_Library::C_Polyhedron(B));
    }

    template <typename R>
    bool 
    Zonotope<R>::inner_subset(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Geometry::inner_subset(Parma_Polyhedra_Library::C_Polyhedron(A),
                                    Parma_Polyhedra_Library::C_Polyhedron(B));
    }

    

    template <typename R>
    bool 
    Zonotope<R>::subset(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return Geometry::subset(Parma_Polyhedra_Library::C_Polyhedron(A),
                              Parma_Polyhedra_Library::C_Polyhedron(B));
    }

    template <typename R>
    bool 
    Zonotope<R>::subset(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      return Geometry::subset(Parma_Polyhedra_Library::C_Polyhedron(A),
                              Parma_Polyhedra_Library::C_Polyhedron(B));
    }

    template <typename R>
    bool 
    Zonotope<R>::subset(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return Geometry::subset(Parma_Polyhedra_Library::C_Polyhedron(A),
                              Parma_Polyhedra_Library::C_Polyhedron(B));
    }

    
    

    template<typename R> 
    Zonotope<R> 
    Zonotope<R>::minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      using namespace Ariadne::LinearAlgebra;
     
      if (A.dimension()!=B.dimension()) {
          throw std::domain_error(
              "minkowski_sum: the two zonotopes have different dimension.");
      }
      
      const Matrix<R> &A_gen=A.generators();
      const Matrix<R> &B_gen=B.generators();
      
      dimension_type col_A=(A_gen).number_of_columns(); 
      dimension_type col_B=(B_gen).number_of_columns(); 
      
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
    Zonotope<R>::minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      using namespace Ariadne::LinearAlgebra;
     
      if (A.dimension()!=B.dimension()) {
          throw std::domain_error(
              "minkowski_difference: the two zonotopes have different dimension.");
      }
      
      const Matrix<R> &A_gen=A.generators();
      const Matrix<R> &B_gen=B.generators();
      
      dimension_type col_A=(A_gen).number_of_columns(); 
      dimension_type col_B=(B_gen).number_of_columns(); 
      
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
    Zonotope<R>::over_approximation(const Rectangle<R> &c, const LinearAlgebra::Matrix< Interval<R> >& A)
    {
      LinearAlgebra::Matrix<R> Amid=LinearAlgebra::centre(A);
      return Zonotope<R>(c+LinearAlgebra::radius_row_sum(A),LinearAlgebra::centre(A));
    }
    


    template <typename R>
    std::ostream&
    Zonotope<R>::write(std::ostream& os) const 
    {
      const Zonotope<R>& z=*this;
      if(z.dimension() > 0) {
        if (!(z.empty())) {
          os << "Zonotope( centre=" << z.centre();
          os << " directions=" << z.generators();
          os << ") ";
        } else {
          os << "Zonotope( )" << std::endl;
        }
      }

      return os;
    }
    
    template <typename R>
    std::istream& 
    Zonotope<R>::read(std::istream& is)
    {
      throw std::domain_error("Not implemented");
    }
      
    

  }
}
