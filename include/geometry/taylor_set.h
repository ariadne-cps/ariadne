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

#ifndef ARIADNE_TAYLOR_SET_H
#define ARIADNE_TAYLOR_SET_H

#include <iosfwd>

#include "base/iterator.h"
#include "base/tribool.h"

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "function/approximate_taylor_model.h"

namespace Ariadne {
  

template<class R> class ApproximateTaylorModel;

class EuclideanSpace;
template<class R> class Point;
template<class R> class Box;
template<class R> class Zonotope;
template<class BS> class ListSet;


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
  typedef typename traits<R>::approximate_arithmetic_type A; 
  typedef typename traits<R>::arithmetic_type F; 
  typedef typename traits<R>::interval_type I; 
 public:
  /*! \brief The real number type. */
  typedef R real_type;
  /*! \brief The type of denotable point contained by the Taylor set. */
  typedef Point<R> state_type;
  /*! \brief The type of space the set lies in. */
  typedef EuclideanSpace space_type;
 private:
  static const uint _max_degree=7u;
  ApproximateTaylorModel<R> _model;
 public:
  //@{
  //! \name Constructors
  /*! \brief Default constructor constructs an empty set. */
  explicit TaylorSet();
  
  /*! \brief Construct a Taylor set of dimension \a d with centre at the origin. */
  explicit TaylorSet(dimension_type d);
  
  /*! \brief Construct a set from a model. */
  explicit TaylorSet(const ApproximateTaylorModel<R>& model);
  
  /*! \brief Construct from a Box. */
  TaylorSet(const Box<R>& r);
  
  /*! \brief Construct from a Zonotope. */
  TaylorSet(const Zonotope<R>& z);
  
  //@}
  
  
  
  //@{ 
  //! \name Data access
  /*! \brief The function model describing the set. */
  const ApproximateTaylorModel<R>& model() const;
  //@}
  
  
  
  //@{
  //! \name Geometric operations.
  /*! \brief The Euclidean space the Taylor set lies in. */
  EuclideanSpace space() const;
  
  /*! \brief The dimension of the Euclidean space the Taylor set lies in. */
  dimension_type dimension() const;
  
  /*! \brief True if the Taylor set is empty. */
  tribool empty() const;
  
  /*! \brief Checks for boundedness. */
  tribool bounded() const;
  
  /*! \brief The centre of the set. */
  Point<R> centre() const;
  
  /*! \brief The radius of the Taylor set. */
  R radius() const;
  
  /*! \brief Subdivide into two smaller pieces. */
  ListSet< TaylorSet<R> > split() const;
  
  /*! \brief Subdivide into smaller pieces in each dimension. */
  ListSet< TaylorSet<R> > subdivide() const;
  
  /*! \brief A rectangle containing the given Taylor set. */
  Box<R> bounding_box() const;

  /*! \brief Test if the set is a superset of a box. (Not currently implemented.) */
  tribool superset(const Box<R>& bx) const;

  /*! \brief Test if the set intersects a box. (Not currently implemented.) */
  tribool intersects(const Box<R>& bx) const;
 
  /*! \brief Test if the set is disjoint from a box. */
  tribool disjoint(const Box<R>& bx) const;

  /*! \brief Test if the set is a subset of a box. */
  tribool subset(const Box<R>& bx) const;

  //@}
  
  
  //@{ 
  //! \name Input/output operations
  /*! \brief Write to an output stream. */
  std::ostream& write(std::ostream& os) const;
  //@}
 private:
  static void _instantiate();
};


template<class R> 
Zonotope<R> zonotope_over_approximation(const TaylorSet<R>& ts);


template<class R> inline
const ApproximateTaylorModel<R>& TaylorSet<R>::model() const {
  return this->_model; }


template<class R> inline
tribool subset(const Box<R>& bx, const TaylorSet<R>& ts) {
  return ts.superset(bx); }

template<class R> inline
tribool superset(const TaylorSet<R>& ts, const Box<R>& bx) {
  return ts.superset(bx); }

template<class R> inline
tribool intersect(const TaylorSet<R>& ts, const Box<R>& bx) {
  return ts.intersects(bx); }

template<class R> inline
tribool disjoint(const TaylorSet<R>& ts, const Box<R>& bx) {
  return ts.disjoint(bx); }

template<class R> inline
tribool subset(const TaylorSet<R>& ts, const Box<R>& bx) {
  return ts.subset(bx); }

template<class R> inline
Box<R> bounding_box(const TaylorSet<R>& ts) {
  return ts.bounding_box(); }

template<class R> inline
ListSet< TaylorSet<R> > split(const TaylorSet<R>& ts) {
  return ts.split(); }

template<class R> inline
std::ostream& operator<<(std::ostream& os, const TaylorSet<R>& ts) {
  return ts.write(os); }




} // namespace Ariadne

#endif /* ARIADNE_TAYLOR_SET_H */
