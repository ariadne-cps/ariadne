/***************************************************************************
 *            approximate_taylor_model.h
 *
 *  Copyright 2008  Pieter Collins
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
 
/*! \file approximate_taylor_model.h
 *  \brief Models of functions based on approximate Taylor series.
 */
#ifndef ARIADNE_APPROXIMATE_TAYLOR_MODEL_H
#define ARIADNE_APPROXIMATE_TAYLOR_MODEL_H

#include <iosfwd>
#include "numeric.h"
#include "vector.h"


namespace Ariadne {

template<class X> class array;
template<class X> class Vector;
template<class X> class Matrix;

class MultiIndex;
template<class X> class SparseDifferential;
template<class D> class DifferentialVector;

class ApproximateTaylorModel;
class FunctionInterface;



ApproximateTaylorModel operator+(const ApproximateTaylorModel&, const ApproximateTaylorModel&);
ApproximateTaylorModel operator-(const ApproximateTaylorModel&, const ApproximateTaylorModel&);
ApproximateTaylorModel operator*(const ApproximateTaylorModel&, const ApproximateTaylorModel&);
ApproximateTaylorModel operator+(const ApproximateTaylorModel&, const Vector<Float>&);


ApproximateTaylorModel project(const ApproximateTaylorModel&, const Slice& slc);
ApproximateTaylorModel project(const ApproximateTaylorModel&, const Range& rng);
ApproximateTaylorModel recentre(const ApproximateTaylorModel&, const Vector<Interval>& bx, const Vector<Float>& pt);
ApproximateTaylorModel restrict(const ApproximateTaylorModel&, const Vector< Interval >& bx);
ApproximateTaylorModel truncate(const ApproximateTaylorModel&, const Vector<Float>&, uint, uint);
std::pair<Vector<Interval>, Matrix<Float> > affine(const ApproximateTaylorModel&);

ApproximateTaylorModel join(const ApproximateTaylorModel& f, const ApproximateTaylorModel& g);

ApproximateTaylorModel embed(const ApproximateTaylorModel&, const Vector<Interval>&, const Vector<Float>&, uint);
ApproximateTaylorModel compose(const ApproximateTaylorModel&, const ApproximateTaylorModel&);
ApproximateTaylorModel inverse(const ApproximateTaylorModel&);
ApproximateTaylorModel implicit(const ApproximateTaylorModel&);
ApproximateTaylorModel derivative(const ApproximateTaylorModel&, uint);
ApproximateTaylorModel antiderivative(const ApproximateTaylorModel&, uint);
ApproximateTaylorModel flow(const ApproximateTaylorModel&);
ApproximateTaylorModel integrate(const ApproximateTaylorModel&, const Float& time);
ApproximateTaylorModel hitting(const ApproximateTaylorModel& vf, const ApproximateTaylorModel& g);
tribool disjoint(const ApproximateTaylorModel& vf, const Vector<Interval>& g);
Vector< Interval > solve(const ApproximateTaylorModel&, const Vector<Float>&);

std::ostream& operator<<(std::ostream&, const ApproximateTaylorModel&);

/* \brief A taylor_model with multivalued output, using a den.
 *  \ingroup FunctionModel
 */
class ApproximateTaylorModel {
  typedef Interval I;
  typedef Float R;
  typedef Float A;
 public:
  typedef Vector<Float> PointType;
  typedef Vector<Interval> BoxType;
  typedef DifferentialVector< SparseDifferential<R> > DifferentialVectorType;
  typedef DifferentialVector< SparseDifferential<I> > IntervalDifferentialVectorType;

  //! \brief Destructor.
  ~ApproximateTaylorModel();
  //! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. 
  ApproximateTaylorModel();
  //! \brief The zero Taylor model in \a as variables with size \a rs image, order \a o and smoothness \a s, defined on the whole space with centre at the origin. 
  ApproximateTaylorModel(uint rs, uint as, ushort o, ushort s);
  
  //! \brief Construct the identity model on a domain. 
  ApproximateTaylorModel(const BoxType& domain);
  
  //! \brief Construct from a domain, centre, and a derivative expansion. 
  ApproximateTaylorModel(const BoxType& domain, const PointType& centre, 
                         const DifferentialVectorType& expansion);
  
  //! \brief Construct from a domain, centre, and two derivative expansions, one for the centre and one over the entire domain. Included for compatibility with TaylorModel class.
  ApproximateTaylorModel(const BoxType& domain, const PointType& centre, 
                         const IntervalDifferentialVectorType& centre_expansion,
                         const IntervalDifferentialVectorType& domain_expansion);
  
  //! \brief Construct an approximate model of \a function over \a domain, using a Taylor expansion of the given \a order. The \a smoothness parameter is not used.
  ApproximateTaylorModel(const BoxType& domain, const FunctionInterface& function,
                         ushort order, ushort smoothness);
  
  //! \brief Construct from a domain, centre, a function, a maximum order of the polynomial expansion and a dummy \a smoothness parameter. 
  ApproximateTaylorModel(const BoxType& domain, const PointType& centre,
                            const FunctionInterface& function,
                            ushort order, ushort smoothness);
  
 
  //! \brief Copy constructor.
  ApproximateTaylorModel(const ApproximateTaylorModel& atm);
  
  //! \brief Copy assignment operator.
  ApproximateTaylorModel& operator=(const ApproximateTaylorModel& atm);
  
  
  // Data access
  //! \brief The data used to define the centre of the Taylor model. 
  const DifferentialVectorType& expansion() const;
  
  // FIXME: A model should not have a dimension!! Needed because we're using 
  // a model directly as a set.
  uint dimension() const { return this->result_size(); }

  // Data access
  //! \brief The order of the Taylor model. 
  ushort order() const;
  //! \brief The smoothness of the function. 
  ushort smoothness() const;
  //! \brief The size of the argument. 
  uint argument_size() const;
  //! \brief The size of the result. 
  uint result_size() const;
  
  /// Resizing
  void resize(uint rs, uint as, ushort d, ushort s);

  //! \brief The domain of validity of the Taylor model. 
  BoxType domain() const;
  //! \brief The centre of the derivative expansion. 
  Vector<R> centre() const;
  //! \brief The range of values the Taylor model can take. 
  BoxType range() const;
  
  //! \brief The range of values the Taylor model can take. 
  BoxType bounding_box() const;
  
  //! \brief Evaluate the Taylor model at the point \a x. 
  BoxType evaluate(const BoxType& x) const;
  PointType evaluate(const PointType& x) const;
  
  //! \brief Compute the derivate of the map at a point. 
  Matrix<I> jacobian(const BoxType& s) const;
  Matrix<A> jacobian(const PointType& s) const;
  
  //! \brief Truncate to a model of lower order and/or smoothness, possibly on a different domain. 
  ApproximateTaylorModel truncate(const BoxType& domain, const Vector<R>& centre, 
                                  ushort order, ushort smoothness) const;
  
  //!
  static ApproximateTaylorModel zero(uint rs, uint as);
  //!
  static ApproximateTaylorModel one(uint as);
  //!
  static ApproximateTaylorModel constant(uint as, const R& c);
  //!
  static ApproximateTaylorModel constant(const BoxType& d, const Vector<R>& c, const Vector<R>& x, ushort o=4, ushort s=1);

  //!
  static ApproximateTaylorModel identity(const BoxType& d);
 
  //!
  static ApproximateTaylorModel affine(const I&, const R&, const R&, const R&, ushort, ushort);

  //! \brief Write to an output stream. 
  std::ostream& write(std::ostream& os) const;
  
#ifdef DOXYGEN
  //! \brief Embed in a larger space.
  friend ApproximateTaylorModel embed(const ApproximateTaylorModel&, const BoxType&, const Vector<R>&, uint);
  //! \brief Composition \f$p\circ q(x)=p(q(x))\f$. 
  friend ApproximateTaylorModel compose(const ApproximateTaylorModel&, const ApproximateTaylorModel&);
  //! \brief Inverse function model \f$p^{-1}\f$. 
  friend ApproximateTaylorModel inverse(const ApproximateTaylorModel&);
  //! \brief Implicit function defined by... 
  friend ApproximateTaylorModel implicit(const ApproximateTaylorModel&);
  //! \brief Derivative with respect to variable \a k. 
  friend ApproximateTaylorModel derivative(const ApproximateTaylorModel&, uint k);
  //! \brief Truncate within \a r to a Taylor model of order at most \a d, putting the error into terms of order \a s. 
  friend ApproximateTaylorModel truncate(const ApproximateTaylorModel& p, const Rectangle& bb, uint d, uint );
#endif
 private:
  template<class X> array< array<X> > _powers(const Vector<X>&) const;
  void _compute_jacobian() const;
  void _set_argument_size(uint n);
  uint _compute_maximum_component_size() const;
 private:
  // Data type;
  class Data;
  Data* _data;
 private:
  //BOOST_CONCEPT_ASSERT((FunctionModelConcept<ApproximateTaylorModel>));
};


} //namespace Ariadne

#endif /* ARIADNE_APPROXIMATE_TAYLOR_MODEL_H */
