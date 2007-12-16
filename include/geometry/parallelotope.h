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

#ifndef ARIADNE_PARALLELOTOPE_H
#define ARIADNE_PARALLELOTOPE_H

#include <iosfwd>

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "geometry/exceptions.h"
#include "geometry/point.h"
#include "geometry/zonotope.h"

namespace Ariadne {
  namespace Geometry {
 
    class basic_set_tag;

    /*!\ingroup BasicSet
     * \brief A parallelotope of arbitrary dimension.
     *
     * A %Parallelotope is a type of Zonotope, and has the same data storage 
     * structure, so is implemented as a subclass. However, it is easier to convert a parallelotope to a
     * Polyhedron defined by its inequalities, since this can be performed by matrix inversion.
     * Hence it is easier to check if a parallelotope contains a point.
     *
     * A parallelotope is always bounded and nonempty. A nondegenerate 
     * parallelotope is regular, meaning that it is the closure of its interior.
     *
     * \b storage Since a parallelotope of dimension \a d has exactly \a d generators,
     * the space required for a parallelotope depends only on the dimension,
     * and is given by \f$d(d+1)\f$.
     *
     */
    template<class XC, class XG=XC>
    class Parallelotope
      : public Zonotope<XC,XG> {
     private:
      typedef typename Numeric::traits<XC>::number_type R;
      typedef typename Numeric::traits<R>::arithmetic_type F;
      typedef typename Numeric::traits<R>::interval_type I;
     public:
      /*! \brief The real number type. */
      typedef R real_type;
      /*! \brief The type of denotable point contained by the parallelotope. */
      typedef Point<R> state_type;
     public:
      //@{ 
      //! \name Constructors
      /*! \brief Default constructor constructs an empty parallelotope of dimension \a n. */
      explicit Parallelotope(dimension_type n=0);
      
      /*! \brief Construct from centre and directions. */
      template<class RC,class RG> explicit Parallelotope(const Point<RC>& c, const LinearAlgebra::Matrix<RG>& m);
       
      /*! \brief Construct from a rectangle. */
      template<class RR> explicit Parallelotope(const Rectangle<RR>& r);
      
      /*! \brief Construct from a string literal. */
      explicit Parallelotope(const std::string& str);
        
      
      /*! \brief Construct from a zonotope. */
      template<class RC,class RG> Parallelotope(const Zonotope<RC,RG>& z);

      /*! \brief Copy constructor. */
      template<class RC,class RG> Parallelotope(const Parallelotope<RC,RG>& original);

      /*! \brief Assign from a rectangle. */
      template<class RR> Parallelotope<XC,XG>& operator=(const Rectangle<RR>& r);
      //@}
      
      
      //@{
      //! \name Geometric operations

      /*! \brief Resize to \a d dimensions. */
      void resize(dimension_type d);
      
      /*! \brief Subdivide into two smaller pieces. */
      ListSet< Parallelotope<XC,XG> > divide() const;
      
      /*! \brief Subdivide into smaller pieces in each dimension. */
      ListSet< Parallelotope<XC,XG> > subdivide() const;      

      /*! \brief An approximation to the volume. */
      F volume() const;
      //@}
      
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
     private:
      static tribool _instantiate();
     private:
      static LinearAlgebra::Matrix<R> compute_generators(const Rectangle<R>& r);
      static void compute_linear_inequalities(LinearAlgebra::Matrix<R>&, LinearAlgebra::Vector<R>&, LinearAlgebra::Vector<R>&);
      LinearAlgebra::Vector<F> coordinates(const Point<R>& s) const;
      void _compute_generators_inverse() const;
     private:
      mutable LinearAlgebra::Matrix<F> _generators_inverse;
    };
    

    /*! \brief Computes an over approximation from a parallelotope (returns the argument). */
    template<class R> 
    Parallelotope<R,R> over_approximation(const Parallelotope<R,R>& p);
    
    /*! \brief Computes an over approximation from an interval parallelotope. */
    template<class R> 
    Parallelotope<R,R> over_approximation(const Parallelotope< Numeric::Interval<R>, R >& p);
    
    /*! \brief Computes an over approximation from an interval parallelotope. */
    template<class R> 
    Parallelotope<Numeric::Interval<R>,R> over_approximation(const Parallelotope< Numeric::Interval<R>, Numeric::Interval<R> >& p);
    
    /*! \brief Computes an over approximation from a zonotope using a qr factorization. */
    template<class R> 
    Parallelotope<R> orthogonal_over_approximation(const Zonotope<R>& z);
    
    /*! \brief Computes an over approximation from an uniform error zonotope using a qr factorization. */
    template<class R> 
    Parallelotope<Numeric::Interval<R>,R> orthogonal_over_approximation(const Zonotope<Numeric::Interval<R>,R>& ez);
    
    /*! \brief Computes an over approximation from an interval zonotope using a qr factorization. */
    template<class R> 
    Parallelotope< Numeric::Interval<R> > orthogonal_over_approximation(const Zonotope< Numeric::Interval<R> >& iz);
    
    
    
    template<class XC, class XG>
    std::istream& operator>>(std::ostream& is, Parallelotope<XC,XG>& p);
    
    template<class XC, class XG>
    std::ostream& operator<<(std::ostream& os, const Parallelotope<XC,XG>& p);
    
  }
}

#include "parallelotope.inline.h"

#endif /* ARIADNE_PARALLELOTOPE_H */
