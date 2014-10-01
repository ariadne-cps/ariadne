/***************************************************************************
 *            c1_taylor_function.h
 *
 *  Copyright 2013  Pieter Collins
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

/*! \file c1_taylor_function.h
 *  \brief Over-approximations of continuously-differentiable functions based on Taylor expansions.
 */

#ifndef ARIADNE_C1_TAYLOR_FUNCTION_H
#define ARIADNE_C1_TAYLOR_FUNCTION_H

#include <iosfwd>
#include "declarations.h"
#include "container.h"
#include "numeric.h"
#include "float.h"
#include "expansion.h"
#include "box.h"

namespace Ariadne {


class MultiIndex;
class C1TaylorSeries;
class C1TaylorFunction;


/*! \ingroup FunctionModelSubModule
 *  \brief A C1TaylorSeries is a univariate function \f$f:\R\rightarrow\R\f$ on an interval \f$[a,b]\f$ is approximated by polynomial \f$p\f$ with error bounds \f$e_0 \geq |f(c)-p(c)|\f$ and \f$e_1\geq sup_{\xi\in D} |f'(\xi)-p'(\xi)|\f$.
 */
class C1TaylorSeries
{
  public:
    std::vector<Float> _coefficients;
    Float _zero_error;
    Float _uniform_error;
    Float _derivative_error;
  private:
    C1TaylorSeries(Nat d);
  public:
    C1TaylorSeries();
    static C1TaylorSeries constant(Float);
    static C1TaylorSeries coordinate();
  public:
    Interval domain() const;
    Nat degree() const;
    Void sweep(Float threshold);
  public:
    friend C1TaylorSeries& operator+=(C1TaylorSeries&, ValidatedNumberType);
    friend C1TaylorSeries& operator*=(C1TaylorSeries&, ValidatedNumberType);
    friend C1TaylorSeries operator+(C1TaylorSeries, C1TaylorSeries);
    friend C1TaylorSeries operator*(C1TaylorSeries, C1TaylorSeries);
    friend ValidatedNumberType evaluate(C1TaylorSeries, ValidatedNumberType);
    friend C1TaylorSeries compose(C1TaylorSeries, C1TaylorSeries);
    friend OutputStream& operator<< (OutputStream& os, const C1TaylorSeries& f);
};

class C1TaylorFunction
{
  public:
    typedef Interval NumericType;
  public:
    Expansion<Float> _expansion;
    Float _zero_error;
    Float _uniform_error;
    Array<Float> _derivative_errors;
  private:
  public:
    C1TaylorFunction() { };
    C1TaylorFunction(Nat as);
  public:
    static C1TaylorFunction constant(Nat as, Float c);
    static C1TaylorFunction coordinate(Nat as, Nat ind);
  public:
    Box domain() const;
    Nat argument_size() const;
    Void sweep(Float threshold);
    C1TaylorFunction& operator=(Interval c);
    Void clear();
  public:
    friend C1TaylorFunction& operator+=(C1TaylorFunction& f, Float c);
    friend C1TaylorFunction& operator*=(C1TaylorFunction& f, Float c);
    friend C1TaylorFunction operator+(C1TaylorFunction f1, C1TaylorFunction f2);
    friend C1TaylorFunction operator*(C1TaylorFunction f1, C1TaylorFunction f2);
    friend ValidatedNumberType evaluate(C1TaylorFunction f, Vector<ValidatedNumberType> x);
    friend C1TaylorFunction compose(C1TaylorSeries f, C1TaylorFunction g);
    friend C1TaylorFunction compose(C1TaylorFunction f, Vector<C1TaylorFunction> g);
    friend OutputStream& operator<< (OutputStream& os, const C1TaylorFunction& f);
};

} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_H
