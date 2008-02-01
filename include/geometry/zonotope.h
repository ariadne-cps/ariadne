/***************************************************************************
 *            zonotope.h
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
 
/*! \file zonotope.h
 *  \brief Zonotopes (affine images of cuboids).
 */

#ifndef ARIADNE_ZONOTOPE_H
#define ARIADNE_ZONOTOPE_H

#include <iosfwd>

#include "base/iterator.h"
#include "base/tribool.h"

#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "function/declarations.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "geometry/point.h"


namespace Ariadne {
  namespace Geometry {
    
    class basic_set_tag;
    template<class X> class Point;
    template<class R> class Box;
    template<class X> class Polyhedron;
    template<class BS> class ListSet;
    template<class R> class ConstraintSet;

    template<class R> class Zonotope;
  


    template<class R> tribool empty(const Zonotope<R>& z);
    template<class R> tribool bounded(const Zonotope<R>& z);
    template<class R> R radius(const Zonotope<R>&);
    template<class R> Box<R> bounding_box(const Zonotope<R>& z);
    

    template<class R> tribool contains(const Zonotope<R>& z, const Point<R>& pt);
    template<class R> tribool disjoint(const Zonotope<R>& z, const Box<R>& r);
    template<class R> tribool intersects(const Zonotope<R>& z, const Box<R>& r);
    template<class R> tribool superset(const Zonotope<R>& r, const Box<R>& z);
    template<class R> tribool subset(const Zonotope<R>& z, const Box<R>& r);
    
    template<class R> tribool subset(const Box<R>& z, const Zonotope<R>& r);
    template<class R> tribool disjoint(const Box<R>& r, const Zonotope<R>& z);

    template<class R> tribool subset(const Zonotope<R>& z, const Polyhedron<R>& p);

    template<class R> ListSet< Zonotope<R> > split(const Zonotope<R>& z);

    template<class R> Zonotope<R> approximation(const Zonotope<R>& z);
    template<class R> Zonotope<R> over_approximation(const Zonotope<R>& z);
    template<class R> Zonotope<R> error_free_over_approximation(const Zonotope<R>&);
    template<class R> Zonotope<R> orthogonal_over_approximation(const Zonotope<R>&);
    template<class R> Zonotope<R> nonsingular_over_approximation(const Zonotope<R>&);
    template<class R>Zonotope<R> cascade_over_approximation(const Zonotope<R>& z, size_type maximum_number_of_blocks);
    
    template<class R> tribool subset(const Zonotope<R>& z, const ConstraintSet<R>& cs);
    template<class R> tribool disjoint(const Zonotope<R>& z, const ConstraintSet<R>& cs);
    template<class R> tribool intersects(const Zonotope<R>& z, const ConstraintSet<R>& cs);
    
    template<class R> Zonotope<R> apply(const Function::AffineModel<R>& am, const Zonotope<R>& z);
        
    template<class R> std::ostream& operator<<(std::ostream& os, const Zonotope<R>& z);
    template<class R> std::istream& operator>>(std::istream& is, Zonotope<R>& z);


    /*!\ingroup BasicSet
     * \brief A zonotope of arbitrary dimension.
     * 
     * A zonotope is a set of the form \f$c+Ge\f$, where \f$||e||_{\infty}\leq1\f$.
     * The columns of the matrix \f$G\f$ are the <em>generators</em> of the 
     * zonotope. 
     *
     * Zonotopes are always bounded.
     * A zonotope always contains its centre point, so can never be empty.
     * However, it may not be regular.
     *
     *
     * The intersection and membership tests may be performed using algorithms from: <br>
     * Guibas, Leonidas J.; Nguyen, An; Zhang, Li, "Zonotopes as bounding volumes."  
     * <i>Proceedings of the Fourteenth Annual ACM-SIAM Symposium on Discrete Algorithms</i> 
     * (Baltimore, MD, 2003),  803--812, ACM, New York, 2003.
     *
     * \b Storage: A %Zonotope in dimension d with n generators is described by
     * d(n+1) real numbers. The ith component of the centre is given by the
     * ith element, and the ith component of the kth generator is given by the
     * (kd+i)th element.
     */
    template<class R>
    class Zonotope {
      typedef Numeric::Interval<R> I;
     private:
      Geometry::Point<R> _centre;
      LinearAlgebra::Matrix<R> _generators;
      LinearAlgebra::Vector<R> _error;
     public:
      //@{
      //! \name Typedefs 
      /*! \brief The type of real number used to describe the zonotope. */
      typedef R real_type;
      //@}

      //@{
      //! \name Constructors
      /*! \brief Default constructor yields a zonotope with dimension zero and no generators. */
      explicit Zonotope();
      /*! \brief Construct a zonotope of dimension \a d with no generators. */
      explicit Zonotope(dimension_type d);
      /*! \brief Construct a zonotope of dimension \a n with centre at the origin and \a m generators. */
      explicit Zonotope(dimension_type d, size_type m);

      /*! \brief Construct from centre, generators, and a uniform error term. */
      explicit Zonotope(const Point<R>& c, const LinearAlgebra::Matrix<R>& G, const LinearAlgebra::Vector<R>& e);
      /*! \brief Construct from centre and generators. */
      explicit Zonotope(const Point<R>& c, const LinearAlgebra::Matrix<R>& G);
      /*! \brief Construct from interval centre and a generator matrix. */
      explicit Zonotope(const Point<I>& c, const LinearAlgebra::Matrix<R>& G);
      /*! \brief Construct from centre and an interval generator matrix. */
      explicit Zonotope(const Point<R>& c, const LinearAlgebra::Matrix<I>& G);
      /*! \brief Construct from an interval centre and an interval generator matrix. */
      explicit Zonotope(const Point<I>& c, const LinearAlgebra::Matrix<I>& G);


      /*! \brief Construct a zonotope of dimension \a d with centre at the origin and \a m generators from the data beginning at \a ptr. */
      template<class XX> explicit Zonotope(dimension_type d, size_type m, const XX* ptr);
      /*! \brief Type conversion constructor. */
      template<class XX> Zonotope(const Zonotope<XX>& z);

      /*! \brief Convert from a box. */
      Zonotope(const Box<R>& r);

      /*! \brief Copy constructor. */
      Zonotope(const Zonotope<R>& z);
      /*! \brief Copy assignment operator. */
      Zonotope<R>& operator=(const Zonotope<R>& z);

      //@}
      
      //@{ 
      //! \name Data access
      /*! \brief The dimension of the zonotope. */
      dimension_type dimension() const;

      /*! \brief The number of generators of the zonotope. */
      size_type number_of_generators() const;

      /*! \brief The domain. */
      Box<R> domain() const;

      /*! \brief The centre. */
      const Point<R>& centre() const;

      /*! \brief The matrix of principle directions. */
      const LinearAlgebra::Matrix<R>& generators() const;
     
      /*! \brief The uniform error bound. */
      const LinearAlgebra::Vector<R>& error() const;
     
      /*! \brief A bounding box for the set. */
      Box<R> bounding_box() const;
     
      /*! \brief The radius of the set in the supremum norm. */
      R radius() const;
     
      /*! \brief Test if the set contains a point. */
      tribool contains(const Point<R>& pt) const;
     
      //@}
      
      
      //@{
      //! \name Geometric binary predicates
      /*! \brief Tests disjointness of \a z and \a r. */
      friend tribool disjoint<>(const Zonotope<R>& z, const Box<R>& r);
      /*! \brief Tests if \a z and \a r intersect. */
      friend tribool intersects<>(const Zonotope<R>& z, const Box<R>& r);
      /*! \brief Tests inclusion of \a r in \a z. */
      friend tribool superset<>(const Zonotope<R>& z, const Box<R>& r);
      /*! \brief Tests inclusion of \a z in \a r. */
      friend tribool subset<>(const Zonotope<R>& z, const Box<R>& r);
      /*! \brief Tests disjointness of \a r and \a z. */
      friend tribool disjoint<>(const Box<R>& z, const Zonotope<R>& r);
      /*! \brief Tests inclusion of \a r in \a z. */
      friend tribool subset<>(const Box<R>& z, const Zonotope<R>& r);

      /*! \brief Tests inclusion of \a z in \a p. */
      friend tribool subset<>(const Zonotope<R>& z, const Polyhedron<R>& p);

      /*! \brief Tests disjointness of \a z and \a cs. */
      friend tribool disjoint<>(const Zonotope<R>& z, const ConstraintSet<R>& cs);
      /*! \brief Tests intersection of \a z and \a cs. */
      friend tribool intersects<>(const Zonotope<R>& z, const ConstraintSet<R>& cs);
      /*! \brief Tests inclusion of \a z and \a cs. */
      friend tribool subset<>(const Zonotope<R>& z, const ConstraintSet<R>& cs);
      //@}

      //@{
      //! \name Approximation operations.
      /*! \brief Compute an simplified approximation of the zonotope \a z. */
      friend Zonotope<R> approximation<>(const Zonotope<R>& z);
      /*! \brief Compute an over-approximation of the zonotope \a z. */
      friend Zonotope<R> over_approximation<>(const Zonotope<R>& z);
      /*! \brief Compute an over-approximation of the zonotope \a z without a uniform error term. */
      friend Zonotope<R> error_free_over_approximation<>(const Zonotope<R>&);
      /*! \brief Compute an over-approximation of the zonotope \a z by a non-coordinate aligned orthotope. */
      friend Zonotope<R> orthogonal_over_approximation<>(const Zonotope<R>&);
      /*! \brief Compute an over-approximation of a zonotope \a z with nonsingular generator matrix. */
      friend Zonotope<R> nonsingular_over_approximation<>(const Zonotope<R>&);
      /*! \brief Compute a cascade-over-approximation of the zonotope \a z with \a maximum_number_of_blocks blocks of \a d generators. */
      friend Zonotope<R> cascade_over_approximation<>(const Zonotope<R>& z, size_type maximum_number_of_blocks);
      //@}
   
      //@{
      //! \name Function operations.
      /*! \brief Compute the image of \a z under a function given by the concrete model \a am. */
      friend Zonotope<R> apply<>(const Function::AffineModel<R>& am, const Zonotope<R>& z);
      //@}
     private:
      static void _instantiate();
    };
  


 
  }
}

#include "zonotope.inline.h"
#include "zonotope.template.h"

#endif /* ARIADNE_ZONOTOPE_H */
