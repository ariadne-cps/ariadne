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
#include "../exceptions.h"

#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/point.h"
#include "../geometry/zonotope.h"

namespace Ariadne {
  namespace Geometry {
 
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
    template<class R>
    class Parallelotope : public Zonotope<R> {
     private:
      typedef typename traits<R>::arithmetic_type F;
      typedef typename traits<R>::interval_type I;
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
      explicit Parallelotope(const LinearAlgebra::Vector<R>& c, const LinearAlgebra::Matrix<R>& m);
      
      /*! \brief Construct from centre and directions. */
      explicit Parallelotope(const Point<R>& c, const LinearAlgebra::Matrix<R>& m);
       
      /*! \brief Construct from a rectangle. */
      template<class Rl> explicit Parallelotope(const Rectangle<Rl>& r);
      
      /*! \brief Construct from a string literal. */
      explicit Parallelotope(const std::string& s);
        
      
      /*! \brief Construct from a zonotope. */
      Parallelotope(const Zonotope<R>& z);

      /*! \brief Copy constructor. */
      Parallelotope(const Parallelotope<R>& original);

      /*! \brief Assign from a rectangle. */
      template<class Rl> Parallelotope<R>& operator=(const Rectangle<Rl>& r);
      //@}
      
      
      //@{
      //! \name Geometric operations
      /*! \brief Tests if the parallelotope contains \a point. */
      tribool contains(const Point<R>& point) const;

#ifdef DOXYGEN
      /*! \brief Tests inclusion of \a A in \a B. */
      friend tribool subset(const Rectangle<R>& A, const Parallelotope<R>& B);
#endif
      
      /*! \brief The vertices of the parallelotope. */
      PointList<F> vertices() const;
      
      /*! \brief Subdivide into two smaller pieces. */
      ListSet< Parallelotope<R> > divide() const;
      
      /*! \brief Subdivide into smaller pieces in each dimension. */
      ListSet< Parallelotope<R> > subdivide() const;      

      /*! \brief An approximation to the volume. */
      R volume() const;
      //@}
      
      /*! \brief Scale the parallelotope by at least \a sf. */
      static Parallelotope<R> scale(const Parallelotope<R>& p, const R& sf);

      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
     private:
      static tribool _instantiate_geometry_operators();
     private:
      static tribool subset(const Rectangle<R>& r, const Parallelotope<R>& p);
      static LinearAlgebra::Matrix<R> compute_generators(const Rectangle<R>& r);
      static void compute_linear_inequalities(LinearAlgebra::Matrix<R>&, LinearAlgebra::Vector<R>&, LinearAlgebra::Vector<R>&);
      LinearAlgebra::Vector<F> coordinates(const Point<R>& s) const;
      void _compute_generators_inverse() const;
     private:
      mutable LinearAlgebra::Matrix<F> _generators_inverse;
    };
    
    
    
     
    template<class R>
    class Parallelotope< Interval<R> > 
      : public Zonotope< Interval<R> >
    {
      typedef typename traits<R>::arithmetic_type F;
      typedef typename traits<R>::interval_type I;
     public:
      /* The real number type. */
      typedef I real_type;
      /* The type of denotable point contained by the parallelotope. */
      typedef Point<I> state_type;
     
      Parallelotope(dimension_type d=0);
      
      template<class Rl1, class Rl2> 
      Parallelotope(const Point<Rl1>& c, const LinearAlgebra::Matrix<Rl2>& g);
      
      Parallelotope(const Rectangle<R>& r);
      
      Parallelotope(const Zonotope<R>& z);

      Parallelotope(const Zonotope<I>& z);
      
      /* Tests if the parallelotope contains \a point. */
      tribool contains(const Point<I>& point) const;

      /* The vertices of the parallelotope. */
      PointList<F> vertices() const;
      
      /* Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      
     private:
      void _compute_generators_inverse() const;
      static tribool _instantiate_geometry_operators();
     private:
      mutable LinearAlgebra::Matrix<I> _generators_inverse;
    };
    
    
    template<class R>
    tribool
    subset(const Rectangle<R>& r, const Parallelotope<R>& p);
    

    /*! \brief Computes an over approximation from a parallelotope (returns the argument). */
    template<class R> 
    Parallelotope<R> over_approximation(const Parallelotope<R>& p);
    
    /*! \brief Computes an over approximation from an interval parallelotope. */
    template<class R> 
    Parallelotope<R> over_approximation(const Parallelotope< Interval<R> >& p);
    
    /*! \brief Computes an over approximation from a zonotope using a qr factorization. */
    template<class R> 
    Parallelotope<R> orthogonal_over_approximation(const Zonotope<R>& z);
    
    /*! \brief Computes an over approximation from an interval zonotope using a qr factorization. */
    template<class R> 
    Parallelotope<R> orthogonal_over_approximation(const Zonotope< Interval<R> >& z);
    
    
    
    template<class R>
    std::istream& operator>>(std::ostream& is, Parallelotope<R>& p);
    
    template<class R>
    std::ostream& operator<<(std::ostream& os, const Parallelotope<R>& p);
    
    template<class R>
    std::ostream& operator<<(std::ostream& os, const Parallelotope< Interval<R> >& p) ;
    
  }
}

#include "parallelotope.inline.h"

#endif /* _ARIADNE_PARALLELOTOPE_H */
