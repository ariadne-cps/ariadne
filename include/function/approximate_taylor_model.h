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
 
#ifndef ARIADNE_APPROXIMATE_TAYLOR_MODEL_H
#define ARIADNE_APPROXIMATE_TAYLOR_MODEL_H

#include <iosfwd>
#include "function/function_model_concept.h"



namespace Ariadne {

class Slice;
template<class X> class Vector;
template<class X> class Matrix;

class MultiIndex;
template<class X> class SparseDifferential;
template<class X> class SparseDifferentialVector;

template<class R> class ApproximateTaylorModel;
template<class R> class FunctionInterface;



template<class R> ApproximateTaylorModel<R> operator+(const ApproximateTaylorModel<R>&, const ApproximateTaylorModel<R>&);
template<class R> ApproximateTaylorModel<R> operator-(const ApproximateTaylorModel<R>&, const ApproximateTaylorModel<R>&);


template<class R> ApproximateTaylorModel<R> project(const ApproximateTaylorModel<R>&, const Slice& slc);
template<class R> ApproximateTaylorModel<R> recentre(const ApproximateTaylorModel<R>&, const Vector< Interval<R> >& bx, const Vector<R>& pt);
template<class R> ApproximateTaylorModel<R> restrict(const ApproximateTaylorModel<R>&, const Vector< Interval<R> >& bx);
template<class R> ApproximateTaylorModel<R> truncate(const ApproximateTaylorModel<R>&, const Vector<R>&, uint, uint);
template<class R> std::pair<Vector<Interval<R> >, Matrix<R> > affine(const ApproximateTaylorModel<R>&);

template<class R> ApproximateTaylorModel<R> join(const ApproximateTaylorModel<R>& f, const ApproximateTaylorModel<R>& g);

template<class R> ApproximateTaylorModel<R> compose(const ApproximateTaylorModel<R>&, const ApproximateTaylorModel<R>&);
template<class R> ApproximateTaylorModel<R> inverse(const ApproximateTaylorModel<R>&);
template<class R> ApproximateTaylorModel<R> implicit(const ApproximateTaylorModel<R>&);
template<class R> ApproximateTaylorModel<R> derivative(const ApproximateTaylorModel<R>&, uint);
template<class R> ApproximateTaylorModel<R> antiderivative(const ApproximateTaylorModel<R>&, uint);
template<class R> ApproximateTaylorModel<R> flow(const ApproximateTaylorModel<R>&);
template<class R> ApproximateTaylorModel<R> integrate(const ApproximateTaylorModel<R>&, const R& time);
template<class R> ApproximateTaylorModel<R> hitting(const ApproximateTaylorModel<R>& vf, const ApproximateTaylorModel<R>& g);
template<class R> Vector< Interval<R> > solve(const ApproximateTaylorModel<R>&, const Vector<R>&);

template<class R> ApproximateTaylorModel<R> implicit1(const ApproximateTaylorModel<R>&, const Vector<R>&);
template<class R> ApproximateTaylorModel<R> implicit2(const ApproximateTaylorModel<R>&, const Vector<R>&, const Vector<R>&);

template<class R> std::ostream& operator<<(std::ostream&, const ApproximateTaylorModel<R>&);

/* \brief A taylor_model with multivalued output, using a den.
 *  \ingroup FunctionModel<R>
 */
template<class R>
class ApproximateTaylorModel {
  typedef Interval<R> I;
  typedef typename traits<R>::approximate_arithmetic_type A;
 public:
  //! \brief The type used to represent real numbers.
  typedef R real_type;
  
  //! \brief Destructor.
  ~ApproximateTaylorModel<R>();
  //! \brief Default constructor constructs a Taylor model of order zero with no arguments and no result variables. 
  ApproximateTaylorModel<R>();
  //! \brief The zero Taylor model in \a as variables with size \a rs image, order \a o and smoothness \a s, defined on the whole space with centre at the origin. 
  ApproximateTaylorModel<R>(uint rs, uint as, ushort o, ushort s);
  
  //! \brief Construct from a domain, centre, and two derivative expansions, one for the centre and one over the entire domain. 
  ApproximateTaylorModel<R>(const Vector<I>& domain, const Vector<A>& centre, 
                            const SparseDifferentialVector<A>& expansion);
  
  //! \brief Construct from a domain, centre, and two derivative expansions, one for the centre and one over the entire domain. Included for compatibility with TaylorModel<R> class.
  ApproximateTaylorModel<R>(const Vector<I>& domain, const Vector<R>& centre, 
                            const SparseDifferentialVector<I>& centre_expansion,
                            const SparseDifferentialVector<I>& domain_expansion);
  
  //! \brief Construct an approximate model of \a function over \a domain, using a Taylor expansion of the given \a order. The \a smoothness parameter is not used.
  ApproximateTaylorModel<R>(const Vector<I>& domain, const FunctionInterface<R>& function,
                            ushort order, ushort smoothness);
  
  //! \brief Construct from a domain, centre, a function, a maximum order of the polynomial expansion and a dummy \a smoothness parameter. 
  ApproximateTaylorModel<R>(const Vector<I>& domain, const Vector<A>& centre,
                            const FunctionInterface<R>& function,
                            ushort order, ushort smoothness);
  
  //! \brief Copy constructor.
  ApproximateTaylorModel<R>(const ApproximateTaylorModel<R>& atm);
  
  //! \brief Copy assignment operator.
  ApproximateTaylorModel<R>& operator=(const ApproximateTaylorModel<R>& atm);
  
  
  // Data access
  //! \brief The data used to define the centre of the Taylor model. 
  const SparseDifferentialVector<A>& expansion() const;
  
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
  Vector<I> domain() const;
  //! \brief The centre of the derivative expansion. 
  Vector<R> centre() const;
  //! \brief The range of values the Taylor model can take. 
  Vector<I> range() const;
  
  //! \brief Evaluate the Taylor model at the point \a x. 
  Vector<I> evaluate(const Vector<I>& x) const;
  Vector<A> evaluate(const Vector<A>& x) const;
  
  //! \brief Compute the derivate of the map at a point. 
  Matrix<I> jacobian(const Vector<I>& s) const;
  Matrix<A> jacobian(const Vector<A>& s) const;
  
  //! \brief Truncate to a model of lower order and/or smoothness, possibly on a different domain. 
  ApproximateTaylorModel<R> truncate(const Vector<I>& domain, const Vector<R>& centre, 
                                     ushort order, ushort smoothness) const;
  
  //!
  static ApproximateTaylorModel<R> zero(uint rs, uint as);
  //!
  static ApproximateTaylorModel<R> one(uint as);
  //!
  static ApproximateTaylorModel<R> constant(uint as, const R& c);
  //!
  static ApproximateTaylorModel<R> identity(const Vector<I>& d, uint o=1u);
  //!
  static ApproximateTaylorModel<R> identity(uint s, uint o=1u);
 
  //! \brief Write to an output stream. 
  std::ostream& write(std::ostream& os) const;
  
#ifdef DOXYGEN
  //! \brief Composition \f$p\circ q(x)=p(q(x))\f$. 
  friend ApproximateTaylorModel<R> compose(const ApproximateTaylorModel<R>&, const ApproximateTaylorModel<R>&);
  //! \brief Inverse function model \f$p^{-1}\f$. 
  friend ApproximateTaylorModel<R> inverse(const ApproximateTaylorModel<R>&);
  //! \brief Implicit function defined by... 
  friend ApproximateTaylorModel<R> implicit(const ApproximateTaylorModel<R>&);
  //! \brief Derivative with respect to variable \a k. 
  friend ApproximateTaylorModel<R> derivative(const ApproximateTaylorModel<R>&, uint k);
  //! \brief Truncate within \a r to a Taylor model of order at most \a d, putting the error into terms of order \a s. 
  friend ApproximateTaylorModel<R> truncate(const ApproximateTaylorModel<R>& p, const Rectangle& bb, uint d, uint s);
#endif
 private:
  static void _instantiate();
  template<class X> array< array<X> > _powers(const Vector<X>&) const;
  void _compute_jacobian() const;
  void _set_argument_size(uint n);
  uint _compute_maximum_component_size() const;
 private:
  class Data;
  Data* _data;
 private:
  BOOST_CONCEPT_ASSERT((FunctionModelConcept< ApproximateTaylorModel<R> >));
};


} //namespace Ariadne

#endif /* ARIADNE_APPROXIMATE_TAYLOR_MODEL_H */
