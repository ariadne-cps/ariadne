/***************************************************************************
 *            sphere.h
 *
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
 
/*! \file sphere.h
 *  \brief Solid spheres in Euclidean space.
 */

#ifndef ARIADNE_SPHERE_H
#define ARIADNE_SPHERE_H

#include <iosfwd>

#include "../base/array.h"
#include "../geometry/exceptions.h"
#include "../geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    class basic_set_tag;

    /*! \ingroup BasicSet
     *  \brief A ball \f$||x-c||\leq r\f$ of arbitrary dimension.
     */
    template<class R>
    class Sphere {
     public:
      /*! \brief A tag describing the type of set. */
      typedef basic_set_tag set_category;
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the simplex. */
      typedef Point<R> state_type;
     
     private:
      /* Simplex's vertices. */
      Point<R> _centre;
      R _radius;
     
     public:
      //@{ 
      //! \name Constructors
      /*! \brief Default constructor constructs standard simplex of dimension \a n. */
      Sphere(size_type n = 0);
    
      /*! \brief Construct from centre and radius. */
      explicit Sphere(const Point<R>& c, const R& r);
     
      /*! \brief Construct from centre and radius. */
      explicit Sphere(const Point<R>& c, R r);
     
      /*! \brief Construct from a string literal. */
      explicit Sphere(const std::string& s);
      
      /*! \brief Copy constructor. */
      Sphere(const Sphere<R>& original);
      
      /*! \brief Copy assignment operator. */
      Sphere<R>& operator=(const Sphere<R>& original);
      //@}
      
      
      //@{ 
      //! \name Comparison operators
      /*! \brief The equality operator.
       */
      bool operator==(const Sphere<R>& other) const;
      
      /*! \brief The inequality operator */
      bool operator!=(const Sphere<R>& other) const;
      //@}
      
      
      //@{ 
      //! \name Data elements
      /*! \brief The centre of the ball. */
      const Point<R>& centre() const;
      
      /*! \brief The radius of the ball. */
      const R& radius() const;
      //@}
      
      
      //@{ 
      //! \name Geometric operations
      /*! \brief The dimension of the Euclidean space the ball lies in. */
      size_type dimension() const;
      
      /*! \brief True if the ball is empty. */
      bool empty() const;
      
      /*! \brief True if the simplex has empty interior. */
      bool empty_interior() const;
      
      /*! \brief Tests if \a point is contained in the ball. */
      tribool contains(const Point<R>& point) const;
      //@}
      
#ifdef DOXYGEN
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Sphere<R>& A, const Sphere<R>& B);
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Rectangle<R>& A, const Sphere<R>& B);
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Sphere<R>& A, const Rectangle<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Sphere<R>& A, const Sphere<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Rectangle<R>& A, const Sphere<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Sphere<R>& A, const Rectangle<R>& B);
      //@}
      
      //@{
      //! \name Geometric binary operations
      /*! \brief The Minkoswi sum of two spheres */
      friend Sphere<R> minkowski_sum(const Sphere<R>& A, const Sphere<R>& B);
      /*! \brief The Minkoswi difference of two spheres */
      friend Sphere<R> minkowski_difference(const Sphere<R>& A, const Sphere<R>& B);
      //@}
#endif      
      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
      
    };
         
    template<class R> 
    R euclidean_distance_square_down(const Point<R>&,const Point<R>&);
    
    template<class R> 
    R euclidean_distance_square_up(const Point<R>&,const Point<R>&);
      
  }
}

#include "sphere.inline.h"

#endif /* ARIADNE_SPHERE_H */
