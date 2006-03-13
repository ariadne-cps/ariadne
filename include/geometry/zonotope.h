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
    template<typename R> bool disjoint(const Zonotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool inner_subset(const Zonotope<R>& A, const Zonotope<R>& B);
    template<typename R> bool subset(const Zonotope<R>& A, const Zonotope<R>& B);

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
      friend bool interiors_intersect <> (const Zonotope<R>& A,
                                          const Zonotope<R>& B);

       /*! \brief Tests disjointness */
      friend bool disjoint <> (const Zonotope<R>& A,
                               const Zonotope<R>& B);

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
      typedef Point<R> Point;
      /*! \brief The type of vectors. */
      typedef Ariadne::LinearAlgebra::vector<R> Vector;
      /*! \brief The type of matrix giving principal directions. */
      typedef Ariadne::LinearAlgebra::matrix<R> Matrix;
     private:
      /* Zonotope's centre. */
      Point _centre;
      
      /* Zonotope's principal directions. */
      Matrix _generators;
     
     public:
      /*! \brief Default constructor constructs an empty zonotope of dimension \a n. */
      explicit Zonotope(size_t n = 0)
        : _centre(n),  _generators(n,0) { }
      
      /*! \brief Construct from centre and directions. */
      explicit Zonotope(const Point& c, const Matrix& m)
        : _centre(c)
      {
	using namespace Ariadne::LinearAlgebra;
	
	if (c.dimension()!=number_of_rows(m)) {
          throw std::domain_error(
              "The the matrix of principal directions does not have the same number of rows as the point dimension.");
        }

        this->_generators= remove_null_columns_but_one(m);
      }
       
      /*! \brief Construct from a rectangle. */
      explicit Zonotope(const Rectangle<Real>& r)
        : _centre(r.dimension())
      {
        if(r.lower_bound(0) > r.upper_bound(0)) {
	  this->_generators=Matrix(r.dimension(),0);
	}
	      
	this->_generators=Matrix(r.dimension(),r.dimension());
        for(size_t i=0; i!=dimension(); ++i) {
          this->_centre[i] = (r.lower_bound(i)+r.upper_bound(i))/2;
          this->_generators(i,i) = (r.upper_bound(i)-r.lower_bound(i))/2;
        }
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
      	if (number_of_columns(this->_generators)==0) return true;
	return false;
      }
      
      /*! \brief True if the zonotope has empty interior. */
      inline bool empty_interior() const {
 	using namespace Ariadne::LinearAlgebra;

        return !independent_rows(this->_generators);
      }
      
      /*! \brief The centre of the zonotope. */
      inline Point centre() const {
        return this->_centre;
      }
      
      /*! \brief The \a n th of principle direction. */
      inline Vector principle_direction(size_t n) const {
	Vector out(this->dimension());

	for (int i=0; i!=this->dimension(); ++i) {
	   out(i)=this->_generators(i,n);
        }
	
        return out;
      }
      
      /*! \brief The matrix of principle directions. */
      inline Matrix principle_directions() const {
        return this->_generators;
      }
      
      /*! \brief Tests if the zonotope contains \a point. */
      inline bool contains(const Point& point) const {
        if (point.dimension()!=this->dimension()) {
          throw std::domain_error("This object and parameter have different space dimensions");
        }  
        
	const Matrix &gen=this->_generators;
	
        if (this->empty()) { return false; }
	      
        throw std::domain_error("Zonotope::contains(...)  not implemented");
	
        return true;
      }
      
      /*! \brief Tests if the interior of the zonotope contains \a point. */
      inline bool interior_contains(const Point& point) const {
	if (point.dimension()!=this->dimension()) {
         throw std::domain_error("This object and parameter have different space dimensions");
        } 
      
        if (this->empty_interior()) { return false; }
     	
	throw std::domain_error("Zonotope::interior_contains(...)  not implemented");

        return true;
      }
      
      /*! \brief The equality operator (not implemented).
       *
       * Not currently implemented, since it requires matching the columns of 
       * the matrix of principal directions. 
       */
      inline bool operator==(const Zonotope<Real>& A) const
      {
        throw std::domain_error("Zonotope::operator==(...)  not implemented");
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
      // The linear inequalities defining the zonotope.
      inline void compute_linear_inequalities(Matrix&, Vector&, Vector&) const;
    };
  
    template<typename R>
    inline 
    void Zonotope<R>::compute_linear_inequalities(Matrix& A, Vector& o, Vector& b) const
    {
      /*
      using namespace Ariadne::LinearAlgebra;
      size_t n=this->dimension();
      const Vector &c=(this->centre()).position_vector();

      //FIXME: can't use matrix inverse here
      A=inverse(this->principle_directions());
      o=A*c;

      b=vector<R>(n);
      for(size_t i=0; i!=n; ++i) {
        b(i)=1;
      }
      */
      throw std::domain_error("Zonotope::compute_linear_inequalities(...)  not implemented");	    
    }
      
    template <typename R>
    Zonotope<R>::operator Polyhedron<R>() const 
    {
      using namespace Ariadne::LinearAlgebra;
      
      typedef typename Zonotope<R>::Real Real;
      typedef typename Zonotope<R>::Point Point;
      
      size_t n = this->dimension();
      
      /* Express in form invs * x - offst in [-bnds,+bnds] */
      matrix<R> invs;
      vector<R> offst;
      vector<R> bnds;
      this->compute_linear_inequalities(invs,offst,bnds);
      
     
      matrix<R> A(2*n,n);
      vector<R> b(2*n);
      
      for(size_t i=0; i!=n; ++i) {
        for(size_t j=0; j!=n; ++j) {
          A(i,j) = -invs(i,j);
          A(i+n,j) = invs(i,j);
        }
        b(i) = bnds(i)-offst(i);
        b(i+n) = bnds(i)+offst(i);
      }
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
      
      return Zonotope<R>(A._centre+B._centre, remove_null_columns_but_one(gen));
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
      
      return Zonotope<R>(A._centre-B._centre, remove_null_columns_but_one(gen));
    }

    
    /*! \brief Tests disjointness */
    template <typename R>
    bool disjoint(const Zonotope<R>& A, const Zonotope<R>& B) 
    {
      return !minkowski_difference(A,B).contains(Point<R>(A.dimension(),0));
    }
       
    /*! \brief Tests intersection of interiors */
    template <typename R>
    bool interiors_intersect(const Zonotope<R>& A,
                             const Zonotope<R>& B) 
    {
      return !minkowski_sum(A,B).interior_contains(Point<R>(A.dimension(),0));
    }
    
    
    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    bool inner_subset(const Zonotope<R>& A,
                      const Zonotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }


    /*! \brief Tests inclusion in an open cover.  */
    template <typename R>
    bool subset_of_open_cover(const Zonotope<R>& A,
                              const ListSet<R, Zonotope >& cover) 
    {
      throw std::domain_error("subset_of_open_cover(Zonotope, std::vector<Zonotope>) not implemented");
    }

    
    /*! \brief Tests inclusion of \a A om the interior of \a B. */
    template <typename R>
    bool inner_subset(const Zonotope<R>& A,
                      const ListSet<R,Zonotope>& B) 
    {
      throw std::domain_error("subset_of_closed_cover(Zonotope, std::vector<Zonotope>) not implemented");
    }



    template <typename R>
    std::ostream&
    operator<<(std::ostream& os, const Zonotope<R>& z) 
    {
      if(z.dimension() > 0) {
        os << "Zonotope(\n  centre=" << z.centre();
        os << "\n  directions=" << z.principle_directions();
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

#endif /* _ARIADNE_ZONOTOPE_H */
