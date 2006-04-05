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

#include "../linear_algebra/zonotopic_vector.h"

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
     * The intersection and membership tests are performed using algorithms from: <br>
     * Guibas, Leonidas J.; Nguyen, An; Zhang, Li, "Zonotopes as bounding volumes."  <i>Proceedings of the Fourteenth Annual ACM-SIAM Symposium on Discrete Algorithms</i> (Baltimore, MD, 2003),  803--812, ACM, New York, 2003.
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
      explicit Zonotope(size_type n = 0)
        : _centre(n),  _generators(n,0) { }
      
      /*! \brief Construct from centre and directions. */
      explicit Zonotope(const State& c, const Matrix& m)
        : _centre(c), _generators(m)
      {
        using namespace Ariadne::LinearAlgebra;
        
        if (c.dimension()!=m.size1()) {
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
        for(size_type i=0; i!=dimension(); ++i) {
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
      Rectangle<R> bounding_box() const;
      
      /*! \brief Convert to a polyhedron. */
      operator Polyhedron<R> () const;
      
      /*! \brief The dimension of the Euclidean space the zonotope lies in. */
      dimension_type dimension() const {
        return this->_centre.dimension();
      }
      
      /*! \brief The number of generators of the zonotope. */
      size_type number_of_generators() const {
        return this->_generators.size2();
      }
      
      /*! \brief True if the zonotope is empty. */
      bool empty() const {
        return false;
      }
      
      /*! \brief True if the zonotope has empty interior. */
      bool empty_interior() const {
        using namespace Ariadne::LinearAlgebra;

        return !independent_rows(this->_generators);
      }
      
      /*! \brief The centre of the zonotope. */
      State centre() const {
        return this->_centre;
      }
      
      /*! \brief The matrix of principle directions. */
      Matrix generators() const {
        return this->_generators;
      }
     
      /*! \brief Tests if the zonotope contains \a point. */
      bool contains(const State& point) const;
      
      /*! \brief Tests if the interior of the zonotope contains \a point. */
      bool interior_contains(const State& point) const;
      
      /*! \brief The equality operator. */
      bool operator==(const Zonotope<Real>& A) const;
      
      /*! \brief The inequality operator */
      bool operator!=(const Zonotope<Real>& A) const {
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
      void minimize_generators(void);
      
      // Order the generator matrix by norm.
      void sort_generators(void);
      
      // The linear inequalities defining the zonotope.
      void compute_linear_inequalities(Matrix&, Vector&) const;
    };
  
    
    

    
    /*! \brief Performs the Minkoswi sum of two zonotopes */
    template<typename R> 
    Zonotope<R> 
    minkowski_sum(const Zonotope<R>& A, const Zonotope<R>& B);
   
    /*! \brief Performs the Minkoswi difference of two zonotopes */
    template<typename R> 
    Zonotope<R> 
    minkowski_difference(const Zonotope<R>& A, const Zonotope<R>& B);
    
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

    
    /*! Compute an over-approximatiing parallelotope to \a z. */
    template <typename R>
    Geometry::Parallelotope<R>
    over_approximating_parallelotope(const Geometry::Zonotope<R>& z);      
    
    template<typename R>
    Zonotope<R>
    operator+(const Rectangle<R>& r, const LinearAlgebra::zonotopic_vector<R>& v);
    
    template<typename R>
    Zonotope<R>
    operator+(const Zonotope<R>& z, const LinearAlgebra::zonotopic_vector<R>& v);
    

  }
}

#endif /* _ARIADNE_ZONOTOPE_H */
