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
 
/*! \file zonotope.h
 *  \brief Zonotopes (affine images of cuboids).
 */

#ifndef _ARIADNE_ZONOTOPE_H
#define _ARIADNE_ZONOTOPE_H

#include <iosfwd>
#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../base/utility.h"
#include "../base/interval.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/parallelotope.h"
#include "../geometry/polyhedron.h"

namespace Ariadne {
  namespace Geometry {
    template < typename R > class Zonotope;

    template < typename R > class Rectangle;
    template < typename R > class Polyhedron;
    template < typename R, template <typename> class BS > class ListSet;

    template<typename R> Zonotope<R> minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B);
    template<typename R> Zonotope<R> minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B);

    template<typename R> bool interiors_intersect(const Zonotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool interiors_intersect(const Zonotope<R>& A, const Rectangle<R>& B);
    template<typename R> bool interiors_intersect(const Rectangle<R>& A, const Zonotope<R>& B);
    template<typename R> bool interiors_intersect(const Zonotope<R>& A, const Parallelotope<R>& B);
    template<typename R> bool interiors_intersect(const Parallelotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool disjoint(const Zonotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool disjoint(const Zonotope<R>& A, const Rectangle<R>& B);
    template<typename R> bool disjoint(const Rectangle<R>& A, const Zonotope<R>& B);
    template<typename R> bool disjoint(const Parallelotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool disjoint(const Zonotope<R>& A, const Parallelotope<R>& B);
    template<typename R> bool inner_subset(const Zonotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool inner_subset(const Rectangle<R>& A, const Zonotope<R>& B);
    template<typename R> bool inner_subset(const Zonotope<R>& A, const Rectangle<R>& B);
    template<typename R> bool inner_subset(const Parallelotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool inner_subset(const Zonotope<R>& A, const Parallelotope<R>& B);
    template<typename R> bool subset(const Zonotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool subset(const Rectangle<R>& A, const Zonotope<R>& B);
    template<typename R> bool subset(const Zonotope<R>& A, const Rectangle<R>& B);
    template<typename R> bool subset(const Parallelotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool subset(const Zonotope<R>& A, const Parallelotope<R>& B);

    template<typename R> bool subset_of_open_cover(const Zonotope<R>& A, const ListSet<R, Zonotope >& U);
    template<typename R> bool inner_subset(const Zonotope<R>& A, const ListSet<R,Zonotope>& U);
    template<typename R> bool subset(const Zonotope<R>& A, const ListSet<R,Zonotope>& B);
    
    template<typename R> std::ostream& operator<<(std::ostream&, const Zonotope<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Zonotope<R>&);

    /*! \brief A zonotope of arbitrary dimension.
     * 
     * A zonotope is a set of the form \f$c+Ae\f$, where \f$||e||_{infty}\leq1\f$.
     * The intersection and membership tests are performed using algorithms from: <br/>
     * Guibas, Leonidas J.; Nguyen, An; Zhang, Li, "Zonotopes as bounding volumes."  <it>Proceedings of the Fourteenth Annual ACM-SIAM Symposium on Discrete Algorithms</it> (Baltimore, MD, 2003),  803--812, ACM, New York, 2003.
     */
    template <typename R>
    class Zonotope {
       /*! \brief Performs the Minkoswi difference of two zonotopes */
      friend Zonotope<R> minkowski_difference <> (const Zonotope<R>& A,
                                          const Zonotope<R>& B);

       /*! \brief Performs the Minkoswi sum of two zonotopes */
      friend Zonotope<R> minkowski_sum <> (const Zonotope<R>& A,
                                          const Zonotope<R>& B);
       
      /*! \brief Tests intersection of interiors. */
      /*friend bool interiors_intersect <> (const Zonotope<R>& A,
                                          const Zonotope<R>& B);
      */
      
       /*! \brief Tests disjointness */
      /*friend bool disjoint <> (const Zonotope<R>& A,
                               const Zonotope<R>& B);
        */
      
      /*! \brief Tests if \a A is a subset of the interior of \a B. */
      friend bool inner_subset <> (const Zonotope<R>& A,
                                   const Zonotope<R>& B);

      /*! \brief Tests inclusion. */
      friend bool subset <> (const Zonotope<R>& A,
                             const Zonotope<R>& B);

      /*! \brief Tests if \a A is a subset of the interior of \a B. */
      friend bool inner_subset <> (const Zonotope<R>& A,
                                   const ListSet<R,::Ariadne::Geometry::Zonotope>& B);

      /*! \brief Tests if \a A is a subset of \a B. */
      friend bool subset <> (const Zonotope<R>& A,
                             const ListSet<R,::Ariadne::Geometry::Zonotope>& B);


      /*! \brief Tests inclusion in an open cover, represented as a ListSet.
       */
      friend bool subset_of_open_cover <> (const Zonotope<R>& A,
                                           const ListSet<R,::Ariadne::Geometry::Zonotope>& B);

     public:
      /*! \brief The type of denotable real number used for the corners. */
      typedef R Real;
      /*! \brief The type of denotable point contained by the rectangle. */
      typedef Point<R> State;
      /*! \brief The type of vectors. */
      typedef Ariadne::LinearAlgebra::vector<R> Vector;
      /*! \brief The type of matrix giving principal directions. */
      typedef Ariadne::LinearAlgebra::matrix<R> Matrix;
     private:
      /* Zonotope's centre. */
      State _centre;
      
      /* Zonotope's principal directions. */
      Matrix _generators;
    
     public:
      /*! \brief Default constructor constructs an empty zonotope of dimension \a n. */
      explicit Zonotope(size_t n = 0)
        : _centre(n),  _generators(n,0) { }
      
      /*! \brief Construct from centre and directions. */
      explicit Zonotope(const State& c, const Matrix& m)
        : _centre(c), _generators(m)
      {
        using namespace Ariadne::LinearAlgebra;
        
        if (c.dimension()!=number_of_rows(m)) {
          throw std::domain_error(
              "The the matrix of principal directions does not have the same number of rows as the point dimension.");
        }

        this->minimize_generators();
      }
       
      /*! \brief Construct from a rectangle. */
      explicit Zonotope(const Rectangle<Real>& r)
        : _centre(r.dimension())
      {
        using namespace LinearAlgebra;

        if(r.lower_bound(0) > r.upper_bound(0)) {
          this->_generators=Matrix(r.dimension(),0);
        }
              
        this->_generators=Matrix(r.dimension(),r.dimension());
        for(size_t i=0; i!=dimension(); ++i) {
          this->_centre[i] = (r.lower_bound(i)+r.upper_bound(i))/2;
          this->_generators(i,i) = (r.upper_bound(i)-r.lower_bound(i))/2;
        }

        this->_generators=remove_null_columns_but_one(this->_generators);
      }

      /*! \brief Construct from a Parallelotope. */
      explicit Zonotope(const Parallelotope<Real>& p)
        : _centre(p.centre()), _generators(p.generators())
      {
        using namespace Ariadne::LinearAlgebra;

        this->minimize_generators();
      }

      /*! \brief Construct from a string literal. */
      explicit Zonotope(const std::string& s)
        : _centre(), _generators()
      {
        std::stringstream ss(s);
        ss >> *this;
      }
      
      /*! \brief Copy constructor. */
      Zonotope(const Zonotope<R>& original)
        : _centre(original._centre),
          _generators(original._generators)
      { }
      
      /*! \brief Copy assignment operator. */
      Zonotope<R>& operator=(const Zonotope<R>& original) {
        if(this != &original) {
          this->_centre = original._centre;
          this->_generators = original._generators;
        }
        return *this;
      }
      
      /*! \brief A rectangle containing the given zonotope. */
      inline Rectangle<R> bounding_box() const {
         using namespace Ariadne::LinearAlgebra;
        Vector offset(this->dimension());
        
        for(size_t i=0; i!=this->dimension(); ++i) {
          for(size_t j=0; j!=number_of_columns(this->_generators); ++j) {
             offset(i) += abs(this->_generators(i,j));
          }
        }
        
        return Rectangle<R>(this->centre()-offset, this->centre()-offset);
      }
      
      /*! \brief Convert to a polyhedron. */
      inline operator Polyhedron<R> () const;
      
      /*! \brief The dimension of the Euclidean space the zonotope lies in. */
      inline size_t dimension() const {
        return (this->_centre).dimension();
      }
      
      /*! \brief True if the zonotope is empty. */
      inline bool empty() const {
         using namespace Ariadne::LinearAlgebra;
        
              if (number_of_columns(this->_generators)==0) return true;
        return false;
      }
      
      /*! \brief True if the zonotope has empty interior. */
      inline bool empty_interior() const {
         using namespace Ariadne::LinearAlgebra;

        return !independent_rows(this->_generators);
      }
      
      /*! \brief The centre of the zonotope. */
      inline State centre() const {
        return this->_centre;
      }
      
      /*! \brief The \a n th of principle direction. */
      inline Vector principle_direction(size_t n) const {
        Vector out(this->dimension());

        for (size_t i=0; i<this->dimension(); i++) {
           out(i)=this->_generators(i,n);
        }
        
        return out;
      }
      
      /*! \brief The matrix of principle directions. */
      inline Matrix principle_directions() const {
        return this->_generators;
      }
     
      /*! \brief Tests if the zonotope contains \a point. */
      inline bool contains(const State& point) const {
        if (point.dimension()!=this->dimension()) {
          throw std::domain_error("This object and parameter have different space dimensions");
        }  
        
        if (this->empty()) { return false; }
      
        Matrix A(0,0);
        Vector b(0);
                
        this->compute_linear_inequalities(A,b);

        Vector result(A*point.position_vector()-b);

        for (size_t i=0; i<result.size(); i++) {
           if (result(i)<0) return false;
        }
        
        return true;
      }
      
      /*! \brief Tests if the interior of the zonotope contains \a point. */
      inline bool interior_contains(const State& point) const {
         using namespace Ariadne::LinearAlgebra;
        
        if (point.dimension()!=this->dimension()) {
         throw std::domain_error("This object and parameter have different space dimensions");
        } 
      
        if (this->empty_interior()) { return false; }
    
        Matrix A(0,0);
        Vector b(0);
                
        this->compute_linear_inequalities(A,b);

        Vector result(A*point.position_vector()-b);
        
        for (size_t i=0; i<result.size(); i++) {
           if (result(i)<=0) return false;
        }
        
        return true;
      }
      
      /*! \brief The equality operator. */
      inline bool operator==(const Zonotope<Real>& A) const
      {
        using namespace LinearAlgebra;
              
        const Matrix &this_gen=this->_generators;
        const Matrix &A_gen=A._generators;
        size_t directions=number_of_columns(A_gen);

        if (!have_same_dimensions(this_gen,A_gen))
          return false;
        
              array<bool> not_found(directions,true);
        bool searching_for_equiv;

        size_t j2;
        for (size_t j=0; j< directions; j++) {
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

          if (searching_for_equiv) 
            return false;
        }

        return true;
      }
      
      /*! \brief The inequality operator */
      inline bool operator!=(const Zonotope<Real>& A) const {
        return !(*this == A);
      }

      friend std::ostream&
      operator<< <> (std::ostream& os, 
                     const Zonotope<R>& r);
      
      friend std::istream&
      operator>> <> (std::istream& is, 
                     Zonotope<R>& r);
      
     private:
      // Minimize the generator matrix
      inline void minimize_generators(void);
      
      // The linear inequalities defining the zonotope.
      inline void compute_linear_inequalities(Matrix&, Vector&) const;
    };
  
    template<typename R>
    inline 
    void Zonotope<R>::minimize_generators(void) {

      using namespace LinearAlgebra;
           
      size_t i,j,i2,j2;
      Matrix &gen=this->_generators;
      gen=remove_null_columns_but_one(gen);
      size_t rows=number_of_rows(gen);
      
      // if the first row is null the zonotope is a point and it is already
      // minimized
      if (find_first_not_null_in_col(gen,0)==rows) {
        return;
      }
      
      size_t cols=number_of_columns(gen);
      size_t min_cols=cols;
      R coef,coef2;
     
      array<size_t> dependences(cols,cols);
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
    }
    
    template<typename R>
    inline 
    void Zonotope<R>::compute_linear_inequalities(Matrix& A, Vector& b) const
    {
      using namespace Ariadne::LinearAlgebra;
     
      const Matrix &gen=this->_generators;
      Matrix Space=trans(gen);
      size_t number_of_generators=number_of_columns(gen);
      R b2;
     
      //compute complanar hyperplanes
      array<size_t> col,row;
      A=lu_decompose(Space,col,row); 

      if (col.size()==row.size()) {
        A=Matrix(2*number_of_generators,this->dimension());
        b=Vector(2*number_of_generators);
        Space=Matrix();
      } else { 
        remove_null_columns(A,row,col);
        Space=compute_space(A,row,col);

        A=Matrix(2*number_of_generators+2*number_of_rows(Space),
                        this->dimension());
        b=Vector(2*number_of_generators+2*number_of_rows(Space));
      }
      // TO IMPROVE: we new compute both (col_i)^T(cols_j) and 
      // (col_j)^T(cols_i)
      for (size_t i=0; i< number_of_generators; i++) {
        b(2*i)=0.0;
        for (size_t j=0; j< number_of_generators; j++) {
          //b(2*i)-=abs(this->principle_direction(i)*this->principle_direction(j));
          b2=0.0;
          for (size_t k=0; k< this->dimension(); k++) {
            b2+=gen(k,i)*gen(k,j);
          }
          b(2*i)-=abs(b2);
        }
        
        b2=0.0;
        for (size_t k=0; k< this->dimension(); k++) {
          b2+=gen(k,i)*this->_centre[k];
        }

        b(2*i+1)=b(2*i)+b2;
        b(2*i)=b(2*i)-b2;        
        
        for (size_t k=0; k< this->dimension(); k++) {
           A(2*i,k)=-gen(k,i);
           A(2*i+1,k)=gen(k,i);
        }
      }

      for (size_t i=0; i< number_of_rows(Space); i++) {
        for (size_t j=0; j< this->dimension(); j++) {
          A(2*(i+number_of_generators),j)=Space(i,j);
          A(2*(i+number_of_generators)+1,j)=-Space(i,j);
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
    
    /*! \brief Performs the Minkoswi sum of two zonotopes */
    template<typename R> 
    inline
    Zonotope<R> minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      using namespace Ariadne::LinearAlgebra;
     
      if (A.dimension()!=B.dimension()) {
          throw std::domain_error(
              "minkowski_sum: the two zonotopes have different dimension.");
      }
      
      dimension_type col_A=number_of_columns(A._generators), 
                     col_B=number_of_columns(B._generators); 
      
      matrix<R> gen(A.dimension(), col_A+col_B);

      for (size_t i=0; i!=col_A; ++i) {
         for (size_t j=0; j!=A.dimension(); ++j) {
            gen(j,i)=A._generators(j,i);
         }
      }

      
      for (size_t i=0; i!=col_B; ++i) {
        for (size_t j=0; j!=A.dimension(); ++j) {
          gen(j,i+col_A)=B._generators(j,i);
        }
      }
      
      return Zonotope<R>(A._centre+(B._centre).position_vector(), 
                         remove_null_columns_but_one(gen));
    }
   
    /*! \brief Performs the Minkoswi difference of two zonotopes */
    template<typename R> 
    inline
    Zonotope<R> minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B)
    {
      using namespace Ariadne::LinearAlgebra;
     
      if (A.dimension()!=B.dimension()) {
          throw std::domain_error(
              "minkowski_difference: the two zonotopes have different dimension.");
      }
      
      dimension_type col_A=number_of_columns(A._generators), 
                     col_B=number_of_columns(B._generators); 
      
      matrix<R> gen(A.dimension(), col_A+col_B);

      for (size_t i=0; i!=col_A; ++i) {
         for (size_t j=0; j!=A.dimension(); ++j) {
          gen(j,i)=A._generators(j,i);
        }
      }

      
      for (size_t i=0; i!=col_B; ++i) {
        for (size_t j=0; j!=A.dimension(); ++j) {
          gen(j,i+col_A)=B._generators(j,i);
         }
      }
      
      return Zonotope<R>(A._centre-(B._centre).position_vector(), 
                         remove_null_columns_but_one(gen));
    }

    
    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool disjoint(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return !minkowski_difference(A,B).contains(Point<R>(A.dimension(),0));
    }
   
    
    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool disjoint(const Rectangle<R>& A, const Zonotope<R>& B) 
    {
      Zonotope<R> z_A(A);

      return disjoint(z_A,B);
    }

    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool disjoint(const Zonotope<R>& A, const Rectangle<R>& B) 
    {
      return disjoint(B,A);
    }

    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool disjoint(const Parallelotope<R>& A, const Zonotope<R>& B) 
    {
      Zonotope<R> z_A(A);

      return disjoint(z_A,B);
    }

    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool disjoint(const Zonotope<R>& A, const Parallelotope<R>& B) 
    {
      return disjoint(B,A);
    }

    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool interiors_intersect(const Zonotope<R>& A,
                             const Zonotope<R>& B) 
    {
      return !minkowski_sum(A,B).interior_contains(Point<R>(A.dimension(),0));
    }
   
    template <typename R>
    inline
    bool interiors_intersect(const Zonotope<R>& A,
                             const Rectangle<R>& B) 
    {
      Zonotope<R> z_B(B);
      return interiors_intersect(A,z_B);
    }
    
    template <typename R>
    inline
    bool interiors_intersect(const Rectangle<R>& A,
                             const Zonotope<R>& B) 
    {
      Zonotope<R> z_A(A);
      return interiors_intersect(z_A,B);
    }
   
    template <typename R>
    inline
    bool interiors_intersect(const Zonotope<R>& A,
                             const Parallelotope<R>& B) 
    {
      Zonotope<R> z_B(B);
      return interiors_intersect(A,z_B);
    }
    
    template <typename R>
    inline
    bool interiors_intersect(const Parallelotope<R>& A,
                             const Zonotope<R>& B) 
    {
      Zonotope<R> z_A(A);
      return interiors_intersect(z_A,B);
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool inner_subset(const Zonotope<R>& A,
                      const Zonotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool inner_subset(const Rectangle<R>& A,
                      const Zonotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool inner_subset(const Zonotope<R>& A,
                      const Rectangle<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    template <typename R>
    inline
    bool inner_subset(const Parallelotope<R>& A,
                      const Zonotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool inner_subset(const Zonotope<R>& A,
                      const Parallelotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }


    /*! \brief Tests inclusion of \a A in \a B. */
    template <typename R>
    inline
    bool subset(const Zonotope<R>& A,
                      const Zonotope<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in \a B. */
    template <typename R>
    inline
    bool subset(const Rectangle<R>& A,
                      const Zonotope<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool subset(const Zonotope<R>& A,
                      const Rectangle<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    template <typename R>
    inline
    bool subset(const Parallelotope<R>& A,
                      const Zonotope<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool subset(const Zonotope<R>& A,
                      const Parallelotope<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    
    /*! \brief Tests inclusion in an open cover.  */
    template <typename R>
    inline
    bool subset_of_open_cover(const Zonotope<R>& A,
                              const ListSet<R, Zonotope >& cover) 
    {
      throw std::domain_error("subset_of_open_cover(Zonotope, std::vector<Zonotope>) not implemented");
    }

    
    /*! \brief Tests inclusion of \a A om the interior of \a B. */
    template <typename R>
    inline
    bool inner_subset(const Zonotope<R>& A,
                      const ListSet<R,Zonotope>& B) 
    {
      throw std::domain_error("subset_of_closed_cover(Zonotope, std::vector<Zonotope>) not implemented");
    }

      
    

  }
}

#endif /* _ARIADNE_ZONOTOPE_H */
