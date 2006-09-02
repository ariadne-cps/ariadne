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

#include "../declarations.h"

#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/zonotope.h"

namespace Ariadne {
  namespace Geometry {
 
    template<> 
    inline bool is_a<Parallelotope,Parallelotope>() { return true; }
    template<> 
    inline bool is_a<Parallelotope,Zonotope>() { return true; }
    template<>
    inline bool is_a<Parallelotope,Polyhedron>() { return true; }

    /* Forward declaration of friends. */
    template<typename R> std::ostream& operator<<(std::ostream&, const Parallelotope<R>&);
    template<typename R> std::istream& operator>>(std::istream&, Parallelotope<R>&);

    /*! \ingroup BasicSet
     *  \brief A parallelotope of arbitrary dimension.
     */
    template <typename R>
    class Parallelotope : public Zonotope<R> {
     private:
      typedef typename numerical_traits<R>::field_extension_type F;
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the parallelotope. */
      typedef Point<R> state_type;
      /*! \brief The type of vector in the underlying space. */
      typedef LinearAlgebra::Vector<R> vector_type;
      /*! \brief The type of matrix giving principal directions. */
      typedef LinearAlgebra::Matrix<R> matrix_type;
     public:
      //@{ 
      //! \name Constructors
      /*! \brief Default constructor constructs an empty parallelotope of dimension \a n. */
      explicit Parallelotope(dimension_type n) : Zonotope<R>(n,n) { }
      
      /*! \brief Construct from centre and directions. */
      explicit Parallelotope(const vector_type& c, const matrix_type& m) 
        : Zonotope<R>(state_type(c),m)
      {
        if (m.size1()!=m.size2()) {
          throw std::domain_error(
              "The the Matrix of principal directions is not a square Matrix");
        }
      }
      
      /*! \brief Construct from centre and directions. */
      explicit Parallelotope(const state_type& c, const matrix_type& m)
        : Zonotope<R>(c,m)
      {
        if (m.size1()!=m.size2()) {
          throw std::domain_error(
              "The the Matrix of principal directions is not a square Matrix");
        }
      }
       
      /*! \brief Construct from a rectangle. */
      explicit Parallelotope(const Rectangle<real_type>& r)
        : Zonotope<R>(r) { }
      
      /*! \brief Construct from a string literal. */
      explicit Parallelotope(const std::string& s) : Zonotope<R>(s) { }
        
      
      /*! \brief Copy constructor. */
      Parallelotope(const Parallelotope<R>& original) : Zonotope<R>(original) { }
      //@}
      
      
      //@{
      //! \name Geometric operations
      /*! \brief Tests if the parallelotope contains \a point. */
      bool contains(const state_type& point) const;
      
      /*! \brief Tests if the interior of the parallelotope contains \a point. */
      bool interior_contains(const state_type& point) const;

      /*! \brief The vertices of the parallelotope. */
      PointList<Rational> vertices() const;
      
      /*! \brief Subdivide into two smaller pieces. */
      ListSet<R,Geometry::Parallelotope> divide() const;
      
      /*! \brief Subdivide into smaller pieces in each dimension. */
      ListSet<R,Geometry::Parallelotope> subdivide() const;      
      //@}
      
      /*! \brief Computes an over approximation from an "interval parallelotope". */
      static Parallelotope<R> over_approximation(const Rectangle<R>& c, const LinearAlgebra::IntervalMatrix<R>& A);
      /*! \brief Scale the parallelotope by at least \a sf. */
      static Parallelotope<R> scale(const Parallelotope<R>& p, const R& sf);

      friend std::istream& operator>> <> (std::istream& is, Parallelotope<R>& r);
     private:
      static LinearAlgebra::Matrix<R> compute_generators(const Rectangle<R>& r);
      static void compute_linear_inequalities(matrix_type&, vector_type&, vector_type&);
      LinearAlgebra::Vector<F> coordinates(const state_type& s) const;
    };
 
    template<typename R>
    inline
    Parallelotope<R> 
    scale(const Parallelotope<R>& p, const R& scale_factor) 
    {
      return Parallelotope<R>::scale(p,scale_factor);
    }
    


    template <typename R>
    inline
    bool 
    subset_of_open_cover(const Parallelotope<R>& A,
                         const ListSet<R, Parallelotope >& cover) 
    {
      throw std::domain_error("subset_of_open_cover(Parallelotope, ListSet<Parallelotope>) not implemented");
    }

    
    template <typename R>
    inline
    bool 
    inner_subset(const Parallelotope<R>& A,
                 const ListSet<R,Parallelotope>& B) 
    {
      throw std::domain_error("inner_subset(Parallelotope, ListSet<Parallelotope>) not implemented");
    }

  }
}

#endif /* _ARIADNE_PARALLELOTOPE_H */
