/***************************************************************************
 *            ellipsoid.h
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
 
/*! \file ellipsoid.h
 *  \brief Solid ellipsoids in Euclidean space.
 */

#ifndef ARIADNE_ELLIPSOID_H
#define ARIADNE_ELLIPSOID_H

#include <iosfwd>

#include "base/array.h"
#include "linear_algebra/matrix.h"
#include "geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    class basic_set_tag;

    /*! \ingroup BasicSet
     *  \brief An ellipsoid \f$(x-c)^T A (x-c)\leq 1\f$ of arbitrary dimension.
     */
    template<class R>
    class Ellipsoid {
     public:
      /*! \brief A tag describing the type of set. */
      typedef basic_set_tag set_category;
      /*! \brief The type of denotable real number used for the corners. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the simplex. */
      typedef Point<R> state_type;
     
     private:
      Point<R> _centre;
      LinearAlgebra::Matrix<R> _bilinear_form;
     
     public:
      //! \name Constructors 
      //@{ 
      /*! \brief Default constructor constructs unit sphere centred at the origin. \a n. */
      Ellipsoid(size_type n = 0);
    
      /*! \brief Construct from centre and bilinear form. */
      explicit Ellipsoid(const Point<R>& c, const LinearAlgebra::Matrix<R>& r);
     
      /*! \brief Construct from a string literal. */
      explicit Ellipsoid(const std::string& s);
      
      /*! \brief Conversion constructor from a sphere. */
      /* FIXME: Conversion only works if R is a field; otherwise division need not be exact */
      Ellipsoid(const Sphere<R>& s);
      
      /*! \brief Copy constructor. */
      Ellipsoid(const Ellipsoid<R>& original);
      
      /*! \brief Copy assignment operator. */
      Ellipsoid<R>& operator=(const Ellipsoid<R>& original);
      //@}
      
      //@{
      //! \name Comparison operators
      /*! \brief The equality operator.
       */
      bool operator==(const Ellipsoid<R>& other) const;
      
      /*! \brief The inequality operator. */
      bool operator!=(const Ellipsoid<R>& other) const;
      //@}
      
      //! \name Data elements
      //@{
      /*! \brief The centre of the ellipsoid. */
      const Point<R>& centre() const;
      
      /*! \brief The bilinear form describing axes of the ellipsoid. */
      const LinearAlgebra::Matrix<R>& bilinear_form() const;
      //@}
      
      //! \name Geometric Operators
      //@{
      /*! \brief The dimension of the Euclidean space the ellipsoid lies in. */
      size_type dimension() const;
      
      /*! \brief True if the ellipsoid is empty. */
      bool empty() const;

      /*! \brief True if the ellipsoid has empty interior. */
      bool empty_interior() const;

      /*! \brief A rectangle containing the ellipsoid. */
      Rectangle<R> bounding_box() const;
      
      /*! \brief Tests if \a point is contained in the ellipsoid. */
      tribool contains(const Point<R>& pt) const;
      //@}
      
#ifdef DOXYGEN
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Ellipsoid<R>& A, const Ellipsoid<R>& B);
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Rectangle<R>& A, const Ellipsoid<R>& B);
      /*! \brief Tests disjointness */
      friend tribool disjoint(const Ellipsoid<R>& A, const Rectangle<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Ellipsoid<R>& A, const Ellipsoid<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Rectangle<R>& A, const Ellipsoid<R>& B);
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Ellipsoid<R>& A, const Rectangle<R>& B);
      //@}
      
      //@{
      //! \name Geometric binary operations
      /*! \brief The Minkoswi sum of two ellipsoids */
      friend Ellipsoid<R> minkowski_sum(const Ellipsoid<R>& A, const Ellipsoid<R>& B);
      /*! \brief The Minkoswi difference of two ellipsoids */
      friend Ellipsoid<R> minkowski_difference(const Ellipsoid<R>& A, const Ellipsoid<R>& B);
      //@}
#endif
           
      //@{ 
      //! \name Input/output operations
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}
      
     private:
      static void _instantiate_geometry_operators();
    };
    
    template<class R>  
    tribool disjoint(const Ellipsoid<R>& A, const Ellipsoid<R>& B); 
    
    template<class R>  
    tribool disjoint(const Ellipsoid<R>& A, const Rectangle<R>& B); 
    
    template<class R>  
    tribool disjoint(const Rectangle<R>& A, const Ellipsoid<R>& B); 

   
    template<class R> 
    tribool subset(const Ellipsoid<R>& A, const Ellipsoid<R>& B); 

    template<class R> 
    tribool subset(const Ellipsoid<R>& A, const Rectangle<R>& B); 
    
    template<class R> 
    tribool subset(const Rectangle<R>& A, const Ellipsoid<R>& B); 
    

    template<class R>
    Ellipsoid<R> minkowski_sum(const Ellipsoid<R>& A, const Ellipsoid<R>& B);

    template<class R>
    Ellipsoid<R> minkowski_difference(const Ellipsoid<R>& A, const Ellipsoid<R>& B);


    template<class R>
    Ellipsoid<R> scale(const Ellipsoid<R>& e, const R& scale_factor);

    
    template<class R> 
    std::ostream& operator<<(std::ostream& os, const Ellipsoid<R>& e);
    
    template<class R> 
    std::istream& operator>>(std::istream& is, Ellipsoid<R>& e);


  }
}

#include "ellipsoid.inline.h"

#endif /* ARIADNE_ELLIPSOID_H */
