/***************************************************************************
 *            taylor_set.h
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
 
/*! \file taylor_set.h
 *  \brief Taylor sets (polynomial images of cuboids).
 */

#ifndef _ARIADNE_TAYLOR_SET_H
#define _ARIADNE_TAYLOR_SET_H

#include <iosfwd>

#include "../base/iterator.h"
#include "../base/tribool.h"

#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/tensor.h"

namespace Ariadne {
  namespace Geometry {

    template<class R> class Point;
    template<class R> class Rectangle;
    template<class R> class Zonotope;

    /*!\ingroup BasicSet
     * \brief The image of a cuboid under a polynomial.
     *
     * \b Storage: A %Zonotope in dimension d with n generators is described by
     * d(n+1) real numbers. The ith component of the centre is given by the
     * ith element, and the ith component of the kth generator is given by the
     * (kd+i)th element.
     */
    template<class R>
    class TaylorSet {
      typedef typename Numeric::traits<R>::arithmetic_type F; 
      typedef typename Numeric::traits<R>::interval_type I; 
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the Taylor set. */
      typedef Point<R> state_type;

     private:
      Point<R> _centre;
      LinearAlgebra::Matrix<R> _generators;
     public:
      //@{
      //! \name Constructors
      /*! \brief Default constructor constructs an empty set. */
      explicit TaylorSet();
     
      /*! \brief Construct a Taylor set of dimension \a d with centre at the origin. */
      explicit TaylorSet(dimension_type d);
     
      /*! \brief Construct from a string literal. */
      explicit TaylorSet(const std::string& s);
      
      /*! \brief Copy constructor. */
      template<class R1> 
      TaylorSet(const TaylorSet<R1>& ts);
      
      /*! \brief Assign from a Rectangle. */
      TaylorSet<R>& operator=(const Rectangle<R>& r);
      
      /*! \brief Assign from a Zonotope. */
      TaylorSet<R>& operator=(const Zonotope<R>& z);
      
      /*! \brief Copy assignment operator. */
      TaylorSet<R>& operator=(const TaylorSet<R>& ts);
      //@}
      
      
      
      //@{ 
      //! \name Data access
      /*! \brief The centre. */
      Point<R> centre() const;
      //@}
      
      

      //@{
      //! \name Geometric operations.
      /*! \brief The dimension of the Euclidean space the Taylor set lies in. */
      dimension_type dimension() const;
      
      /*! \brief True if the Taylor set is empty. */
      tribool empty() const;
      
      /*! \brief Checks for boundedness. */
      tribool bounded() const;
      
      /*! \brief The radius of the Taylor set. */
      R radius() const;
      
      /*! \brief Subdivide into two smaller pieces. */
      ListSet<R,Geometry::TaylorSet> divide() const;
      
      /*! \brief Subdivide into smaller pieces in each dimension. */
      ListSet<R,Geometry::TaylorSet> subdivide() const;
      
      /*! \brief A rectangle containing the given Taylor set. */
      Rectangle<R> bounding_box() const;
      
      //@}
      
#ifdef DOXYGEN
      //@{
      //! \name Geometric approximation operations
      /*! \brief The Minkoswi sum of two zonotopes */
      friend Zonotope<F> over_approximation(const TaylorSet<R>& A);
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
  
    template<class R> Zonotope<R> over_approximation(const TaylorSet<R>&);
 
    template<class R>
    std::ostream& operator<<(std::ostream& os, const Zonotope<R>& z);
    
    template<class R> 
    std::istream& operator>>(std::istream& is, Zonotope<R>& z);


  }
}

#endif /* _ARIADNE_TAYLOR_SET_H */
