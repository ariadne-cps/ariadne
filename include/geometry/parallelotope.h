/***************************************************************************
 *            parallelotope.h
 *
 *  6 January 2006
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
 
/*! \file parallelotope.h
 *  \brief Parallelotopes.
 */

#ifndef _ARIADNE_PARALLELOTOPE_H
#define _ARIADNE_PARALLELOTOPE_H

#include <iosfwd>

#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../base/utility.h"
#include "../base/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/polyhedron.h"
#include "../geometry/zonotope.h"
#include "../geometry/geometry_declarations.h"

namespace Ariadne {
  namespace Geometry {
    template < typename R > class Parallelotope;

    template < typename R > class Rectangle;
    template < typename R > class Polyhedron;
    template < typename R > class Zonotope;

    template < typename R, template <typename> class BS > class ListSet;

    template <typename R> Parallelotope<R> intersection(const Parallelotope<R>& A, const Parallelotope<R>& B);
    template <typename R> Parallelotope<R> regular_intersection(const Parallelotope<R>& A, const Parallelotope<R>& B);

    template<typename R> bool interiors_intersect(const Parallelotope<R>& A, const Parallelotope<R>& B);
    template<typename R> bool disjoint(const Parallelotope<R>& A, const Parallelotope<R>& B);
    template<typename R> bool inner_subset(const Parallelotope<R>& A, const Parallelotope<R>& B);
    template<typename R> bool subset(const Parallelotope<R>& A, const Parallelotope<R>& B);

    template<typename R> bool subset_of_open_cover(const Parallelotope<R>& A, const ListSet<R, Parallelotope >& list);
    template<typename R> bool inner_subset(const Parallelotope<R>& rect, const ListSet<R,Parallelotope>& A);
    template<typename R> bool subset(const Parallelotope<R>& rect, const ListSet<R,Parallelotope>& A);
    
    template<typename R> std::ostream& operator<<(std::ostream&, const Parallelotope<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Parallelotope<R>&);

    /*! \brief A parallelotope of arbitrary dimension.
     */
    template <typename R>
    class Parallelotope {
      /*! \brief Makes intersection */
      friend Parallelotope<R> intersection <> (const Parallelotope<R>& A,
                                           const Parallelotope<R>& B);

      /*! \brief Makes intersection of interiors */
      friend Parallelotope<R> regular_intersection <> (const Parallelotope<R>& A,
                                                   const Parallelotope<R>& B);

       /*! \brief Tests intersection of interiors. */
      friend bool interiors_intersect <> (const Parallelotope<R>& A,
                                          const Parallelotope<R>& B);

       /*! \brief Tests disjointness */
      friend bool disjoint <> (const Parallelotope<R>& A,
                               const Parallelotope<R>& B);

      /*! \brief Tests if \a A is a subset of the interior of \a B. */
      friend bool inner_subset <> (const Parallelotope<R>& A,
                                   const Parallelotope<R>& B);

      /*! \brief Tests inclusion. */
      friend bool subset <> (const Parallelotope<R>& A,
                             const Parallelotope<R>& B);

      /*! \brief Tests if \a A is a subset of the interior of \a B. */
      friend bool inner_subset <> (const Parallelotope<R>& A,
                                   const ListSet<R,Ariadne::Geometry::Parallelotope>& B);

      /*! \brief Tests if \a A is a subset of \a B. */
      friend bool subset <> (const Parallelotope<R>& A,
                             const ListSet<R,Ariadne::Geometry::Parallelotope>& B);


      /*! \brief Tests inclusion in an open cover, represented as a ListSet.
       */
      friend bool subset_of_open_cover <> (const Parallelotope<R>& A,
                                           const ListSet<R,Ariadne::Geometry::Parallelotope>& B);

     public:
      /*! \brief The unsigned integer type used to denote temphe array positions. */
      typedef size_t size_type;
      /*! \brief The type of denotable real number used for the corners. */
      typedef R Real;
      /*! \brief The type of denotable point contained by the rectangle. */
      typedef Point<R> Point;
      /*! \brief The type of matrix giving principal directions. */
      typedef Ariadne::LinearAlgebra::vector<R> Vector;
      /*! \brief The type of matrix giving principal directions. */
      typedef Ariadne::LinearAlgebra::matrix<R> Matrix;
     private:
      /* Parallelotope's centre. */
      Point _centre;
      
      /* Parallelotope's principal directions. */
      Matrix _generators;
      
     public:
      /*! \brief Default constructor constructs an empty parallelotope of dimension \a n. */
      inline explicit Parallelotope(size_type n = 0)
        : _centre(n),  _generators(n,n) { }
      
      /*! \brief Construct from centre and directions. */
      inline explicit Parallelotope(const Point& c, const Matrix& m)
        : _centre(c), _generators(m)
      {
        if (LinearAlgebra::number_of_rows(m)!=LinearAlgebra::number_of_columns(m)) {
          throw std::domain_error(
              "The the matrix of principal directions is not a square matrix");
        }
        
        if (c.dimension()!=LinearAlgebra::number_of_rows(m)) {
          throw std::domain_error("The centre and directions have different dimensions.");
        }
      }
      
      /*! \brief Construct from a Rectangle. */
      inline explicit Parallelotope(const Rectangle<Real>& r)
        : _centre(r.dimension()), _generators(r.dimension(),r.dimension())
      {
        for(size_type i=0; i!=dimension(); ++i) {
          _centre[i] = (r.lower_bound(i)+r.upper_bound(i))/2;
          _generators(i,i) = (r.upper_bound(i)-r.lower_bound(i))/2;
        }
      }
      
      /*! \brief Construct from a string literal. */
      inline explicit Parallelotope(const std::string& s)
        : _centre(), _generators()
      {
        std::stringstream ss(s);
        ss >> *this;
      }
      
      /*! \brief Copy constructor. */
      Parallelotope(const Parallelotope<R>& original)
        : _centre(original._centre),
          _generators(original._generators)
      { }
      
      /*! \brief Copy assignment operator. */
      Parallelotope<R>& operator=(const Parallelotope<R>& original) {
        if(this != &original) {
          this->_centre = original._centre;
          this->_generators = original._generators;
        }
        return *this;
      }
      
      /*! \brief The dimension of the Euclidean space the parallelotope lies in. */
      inline size_type dimension() const {
        return (this->_centre).dimension();
      }
      
      /*! \brief True if the parallelotope is empty. */
      inline bool empty() const {
        // FIXME: This is probably ok since we're not talking about interior.
        return false;
      }
      
      /*! \brief True if the paralleltope has empty interior. */
      inline bool empty_interior() const {
        using namespace Ariadne::LinearAlgebra;

        return !independent_rows(this->_generators);     
      }
      
      /*! \brief The centre of the parallelotope. */
      inline Point centre() const {
        return this->_centre;
      }
      
      /*! \brief The \a n th of principle direction. */
      inline Vector generator(size_type n) const {
        return column(this->_generators,n);
      }
      
      /*! \brief The matrix of principle directions. */
      inline Matrix generators() const {
        return this->_generators;
      }
     
      /*! \brief The equality operator */
      bool operator==(const Parallelotope<R>& A) const;
      
      /*! \brief The inequality operator */
      bool operator!=(const Parallelotope<R>& A) const; 

      /*! \brief A rectangle containing the given parallelotope. */
      Rectangle<R> bounding_box() const;
      
      /*! \brief Convert to a zonotope. */
      operator Zonotope<R> () const;
      
      /*! \brief Convert to a polyhedron. */
      operator Polyhedron<R> () const;
       
      /*! \brief Tests if the parallelotope contains \a point. */
      bool contains(const Point& point) const;
      
      /*! \brief Tests if the interior of the parallelotope contains \a point. */
      bool interior_contains(const Point& point) const;

      /*! \brief Subdivide into smaller pieces. */
      ListSet<R,Ariadne::Geometry::Parallelotope> subdivide() const;
      
      /* \brief Tests disjointness with a rectangle */
      bool disjoint(const Rectangle<R>& r) const;
      
      friend std::ostream&
      operator<< <> (std::ostream& os, 
                     const Parallelotope<R>& r);
      
      friend std::istream&
      operator>> <> (std::istream& is, 
                     Parallelotope<R>& r);
      
     private:
      void compute_linear_inequalities(Matrix&, Vector&, Vector&) const;
      Vector coordinates(const Point& s) const;
    };
  
    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool 
    disjoint(const Parallelotope<R>& A, const Parallelotope<R>& B) 
    {
      return disjoint(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    
    /*! \brief Tests disjointness */
    template <typename R>
    inline 
    bool 
    disjoint(const Parallelotope<R>& A, const Rectangle<R>& B) 
    {
      return A.disjoint(B);
    }
    
    /*! \brief Tests disjointness */
    template <typename R>
    inline
    bool 
    disjoint(const Rectangle<R>& A, const Parallelotope<R>& B) 
    {
      return B.disjoint(A);
    }
    
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool 
    interiors_intersect(const Parallelotope<R>& A,
                        const Parallelotope<R>& B) 
    {
      return interiors_intersect(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool 
    interiors_intersect(const Parallelotope<R>& A,
                        const Rectangle<R>& B) 
    {
      return interiors_intersect(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    
    /*! \brief Tests intersection of interiors */
    template <typename R>
    inline
    bool
    interiors_intersect(const Rectangle<R>& A,
                        const Parallelotope<R>& B) 
    {
      return interiors_intersect(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    
    
    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool
    inner_subset(const Parallelotope<R>& A,
                 const Parallelotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool
    inner_subset(const Parallelotope<R>& A,
                 const Rectangle<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion of \a A in the interior of \a B. */
    template <typename R>
    inline
    bool
    inner_subset(const Rectangle<R>& A,
                 const Parallelotope<R>& B) 
    {
      return inner_subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }

    /*! \brief Tests inclusion */
    template <typename R>
    inline
    bool 
    subset(const Parallelotope<R>& A, 
           const Parallelotope<R>& B) 
    {
      return subset(Polyhedron<R>(A), Polyhedron<R>(B));
    }
    

    /*! \brief Tests inclusion in an open cover.  */
    template <typename R>
    inline
    bool 
    subset_of_open_cover(const Parallelotope<R>& A,
                         const ListSet<R, Parallelotope >& cover) 
    {
      throw std::domain_error("subset_of_open_cover(Parallelotope, std::vector<Parallelotope>) not implemented");
    }

    
    /*! \brief Tests inclusion of \a A om the interior of \a B. */
    template <typename R>
    inline
    bool 
    inner_subset(const Parallelotope<R>& A,
                 const ListSet<R,Parallelotope>& B) 
    {
      throw std::domain_error("subset_of_closed_cover(Parallelotope, std::vector<Parallelotope>) not implemented");
    }

  }
}

#endif /* _ARIADNE_PARALLELOTOPE_H */
