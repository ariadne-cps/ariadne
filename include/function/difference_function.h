/***************************************************************************
 *            function/difference_function.h
 *
 *  Copyright  2007  Pieter Collins
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
 
/*! \file function/difference_function.h
 */

#ifndef ARIADNE_DIFFERENCE_FUNCTION_H
#define ARIADNE_DIFFERENCE_FUNCTION_H

#include "base/types.h"

#include "numeric/declarations.h"
#include "numeric/traits.h"

#include "linear_algebra/declarations.h"
#include "linear_algebra/vector.h"

#include "differentiation/taylor_derivative.h"
#include "differentiation/sparse_differential.h"
#include "function/function_interface.h"


namespace Ariadne {


/*!
 * \brief A class representing the difference of a function and the identity.
 *
 * Useful for computing fixed points.
 */
template<class R> class DifferenceFunction
  : public FunctionInterface<R>
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  typedef typename traits<R>::arithmetic_type F;
 public:
  /*!\brief Construct from a map \a f, which must have the same argument dimension as result dimension. */
  DifferenceFunction(const FunctionInterface<R>& f);
  /*! \brief Make a copy (clone) of the vector field. */
  DifferenceFunction<R>* clone() const;
  /*!\brief The dimension of the space the map acts on. */
  virtual smoothness_type smoothness() const;
  /*!\brief The dimension of the space the map acts on. */
  virtual size_type result_size() const;
  /*!\brief The dimension of the space the map acts on. */
  virtual size_type argument_size() const;
  /*!\brief Evaluate the function \f$f(x)-x\f$, where \f$f\f$ is the map used to construct the difference map. */
  virtual Vector<F> evaluate(const Vector<F>& p) const;
  /*!\brief Evaluate the derivative of function \f$f(x)-x\f$, which is \f$Df(x)-I\f$. */
  virtual Matrix<F> jacobian(const Vector<F>& p) const;
  /*!\brief Evaluate the derivative of function \f$f(x)-x\f$, which is \f$Df(x)-I\f$. */
  virtual TaylorDerivative<F> derivative(const Vector<F>& p, const smoothness_type& s) const;
  /*!\brief Evaluate the derivative of function \f$f(x)-x\f$, which is \f$Df(x)-I\f$. */
  virtual SparseDifferentialVector<A> expansion(const Vector<A>& p, const smoothness_type& s) const;
  /*!\brief The name of the class. */
  virtual std::string name() const;
  /*!\brief Write to an output stream. */
  virtual std::ostream& write(std::ostream&) const;
 private:
  const FunctionInterface<R>& _base;
};


template<class R> 
inline
DifferenceFunction<R>::DifferenceFunction(const FunctionInterface<R>& f)
  : _base(f)
{ 
  if(f.argument_size()!=f.result_size()) { 
    throw IncompatibleSizes("DifferenceDifferenceFunction(Map f): The argument dimension must equal the result dimension"); 
  } 
}


template<class R> 
inline
DifferenceFunction<R>* 
DifferenceFunction<R>::clone() const 
{ 
  return new DifferenceFunction<R>(this->_base); 
}


template<class R> 
inline
smoothness_type 
DifferenceFunction<R>::smoothness() const 
{ 
  return _base.smoothness(); 
}


template<class R> 
inline
size_type 
DifferenceFunction<R>::argument_size() const 
{ 
  return _base.argument_size(); 
}

template<class R> 
inline
size_type 
DifferenceFunction<R>::result_size() const 
{ 
  return _base.result_size(); 
}


template<class R> 
inline
Vector<typename DifferenceFunction<R>::F> 
DifferenceFunction<R>::evaluate(const Vector<F>& p) const 
{
  return _base.evaluate(p)-p; 
}


template<class R> 
inline
Matrix<typename DifferenceFunction<R>::F>
DifferenceFunction<R>::jacobian(const Vector<F>& p) const 
{
  Matrix<F> d=_base.jacobian(p);
  Matrix<F> i=Matrix< Interval<R> >::identity(this->result_size());
  return d-i; 
}

template<class R> 
inline
TaylorDerivative<typename DifferenceFunction<R>::F>
DifferenceFunction<R>::derivative(const Vector<F>& p, const smoothness_type& s) const 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R> 
inline
SparseDifferentialVector<typename DifferenceFunction<R>::A>
DifferenceFunction<R>::expansion(const Vector<A>& p, const smoothness_type& s) const 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R> 
inline
std::string 
DifferenceFunction<R>::name() const 
{ 
  return "DifferenceFunction"; 
}

template<class R> 
inline
std::ostream&
DifferenceFunction<R>::write(std::ostream& os) const 
{ 
  return os << "DifferenceFunction"; 
}



} // namespace Ariadne


#endif /* ARIADNE_DIFFERENCE_FUNCTION_H */
