/***************************************************************************
 *            zonotope.h
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
 * 
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
 *  \brief Zonotopes in Euclidean space.
 */

#ifndef ARIADNE_ZONOTOPE_H
#define ARIADNE_ZONOTOPE_H

#include <iosfwd>

#include "tribool.h"

#include "vector.h"
#include "matrix.h"

namespace Ariadne {

template<class X> class Vector;
template<class X> class Matrix;
 
class Point;
class Box;
class Zonotope;
template<class BS> class ListSet;

class AffineModel;



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

class Zonotope {
 private:
  Vector<Float> _centre;
  Matrix<Float> _generators;
  Vector<Float> _error;
 public:
  //@{
  //! \name Constructors
  /*! \brief Default constructor yields a zonotope with dimension zero and no generators. */
  explicit Zonotope();
  /*! \brief Construct a zonotope of dimension \a d with no generators. */
  explicit Zonotope(uint d);
  /*! \brief Construct a zonotope of dimension \a n with centre at the origin and \a m generators. */
  explicit Zonotope(uint d, uint m);
  
  /*! \brief Construct from centre, generators, and a uniform error term. */
  explicit Zonotope(const Vector<Float>& c, const Matrix<Float>& G, const Vector<Float>& e);
  /*! \brief Construct from centre and generators. */
  explicit Zonotope(const Vector<Float>& c, const Matrix<Float>& G);
  /*! \brief Construct from interval centre and a generator matrix. */
  explicit Zonotope(const Vector<Interval>& c, const Matrix<Float>& G);
  /*! \brief Construct from centre and an interval generator matrix. */
  explicit Zonotope(const Vector<Float>& c, const Matrix<Interval>& G);
  /*! \brief Construct from an interval centre and an interval generator matrix. */
  explicit Zonotope(const Vector<Interval>& c, const Matrix<Interval>& G);
  
  
  /*! \brief Construct a zonotope of dimension \a d with centre at the origin and \a m generators from the data beginning at \a ptr. */
  template<class XX> explicit Zonotope(uint d, uint m, const XX* ptr);
  
  /*! \brief Convert from a box. */
  Zonotope(const Box& r);
  /*! \brief Copy constructor. */
  Zonotope(const Zonotope& z);
  /*! \brief Copy assignment operator. */
  Zonotope& operator=(const Zonotope& z);
  
  //@}
  
  //@{ 
  //! \name Data access
  /*! \brief The dimension of the zonotope. */
  uint dimension() const;
  
  /*! \brief The number of generators of the zonotope. */
  uint number_of_generators() const;
  
  /*! \brief The domain. */
  Vector<Interval> domain() const;
  
  /*! \brief The centre. */
  const Vector<Float>& centre() const;
  
  /*! \brief The matrix of principle directions. */
  const Matrix<Float>& generators() const;
  
  /*! \brief The uniform error bound. */
  const Vector<Float>& error() const;
  
  /*! \brief A bounding box for the set. */
  Vector<Interval> bounding_box() const;
  
  /*! \brief The radius of the set in the supremum norm. */
  Float radius() const;
  
  /*! \brief Test if the set contains a point. */
  tribool contains(const Point& pt) const;
  
  //@}
  
  
  //@{
  //! \name Geometric binary predicates
  /*! \brief Tests disjointness of \a z and \a r. */
  friend tribool disjoint(const Zonotope& z, const Box& r);
  /*! \brief Tests if \a z and \a r intersect. */
  friend tribool intersects(const Zonotope& z, const Box& r);
  /*! \brief Tests inclusion of \a z in \a r. */
  friend tribool subset(const Zonotope& z, const Box& r);
  /*! \brief Tests disjointness of \a r and \a z. */
  friend tribool disjoint(const Box& r, const Zonotope& z);
  //@}
  
  //@{
  //! \name Approximation operations.
  /*! \brief Compute an simplified approximation of the zonotope \a z. */
  friend Zonotope approximation(const Zonotope& z);
  /*! \brief Compute an over-approximation of the zonotope \a z. */
  friend Zonotope over_approximation(const Zonotope& z);
  /*! \brief Compute an over-approximation of the zonotope \a z without a uniform error term. */
  friend Zonotope error_free_over_approximation(const Zonotope&);
  /*! \brief Compute an over-approximation of the zonotope \a z by a non-coordinate aligned orthotope. */
  friend Zonotope orthogonal_over_approximation(const Zonotope&);
  /*! \brief Compute an over-approximation of a zonotope \a z with nonsingular generator matrix. */
  friend Zonotope nonsingular_over_approximation(const Zonotope&);
  /*! \brief Compute a cascade-over-approximation of the zonotope \a z with \a b blocks of \a d generators. */
  friend Zonotope cascade_over_approximation(const Zonotope& z, uint b);
  //@}
  
  //@{
  //! \name Function operations.
  /*! \brief Compute the image of \a z under a function given by the concrete model \a am. */
  friend Zonotope apply(const AffineModel& am, const Zonotope& z);
  //@}
};


tribool empty(const Zonotope& z);
tribool bounded(const Zonotope& z);
Float radius(const Zonotope& z);
Box bounding_box(const Zonotope& z);


tribool contains(const Zonotope& z, const Point& pt);
tribool disjoint(const Zonotope& z, const Box& r);
tribool intersects(const Zonotope& z, const Box& r);
tribool subset(const Zonotope& z, const Box& r);


ListSet<Zonotope> split(const Zonotope& z);

Zonotope approximation(const Zonotope& z);
Zonotope over_approximation(const Zonotope& z);
Zonotope error_free_over_approximation(const Zonotope&);
Zonotope orthogonal_over_approximation(const Zonotope&);
Zonotope nonsingular_over_approximation(const Zonotope&);
Zonotope cascade_over_approximation(const Zonotope& z, uint b);

Zonotope apply(const AffineModel& am, const Zonotope& z);

std::ostream& operator<<(std::ostream& os, const Zonotope& z);
std::istream& operator>>(std::istream& is, Zonotope& z);


template<class X> inline
Zonotope::Zonotope(uint d, uint m, const X* ptr)
  : _centre(d,ptr), _generators(d,m,ptr+d), _error(d)
{
}


} // namespace Ariadne

#endif // ARIADNE_ZONOTOPE_H
