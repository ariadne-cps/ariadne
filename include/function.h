/***************************************************************************
 *            function.h
 *
 *  Copyright 2008-9  Pieter Collins
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

/*! \file function.h
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_H
#define ARIADNE_FUNCTION_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function_interface.h"

#include "macros.h"
#include "pointer.h"
#include "container.h"

#include "vector.h"
#include "matrix.h"
#include "polynomial.h"
#include "affine.h"
#include "taylor_model.h"
#include "differential.h"
#include "real.h"

namespace Ariadne {

class Real;
template<class X> class Vector;
template<class X> class Polynomial;
template<class LHS, class RHS> class Assignment;
template<class T> class Variable;
template<class T> class Expression;
template<class T> class Space;

typedef Variable<Real> RealVariable;
typedef Expression<Real> RealExpression;
typedef Space<Real> RealSpace;
typedef Assignment<RealVariable,RealExpression> RealAssignment;

class ScalarFunction;
class VectorFunction;

//! A scalar function \f$f:\R^n\rightarrow\R\f$.
class ScalarFunction
{
    typedef uint Nat;
    typedef std::ostream OStream;
  public:
    static ScalarFunction constant(Nat n, double c);
    static ScalarFunction constant(Nat n, Real c);
    static ScalarFunction variable(Nat n, uint i);

    explicit ScalarFunction(Nat n=0u);
    ScalarFunction(const Expression<Real>& e, const Space<Real>& s);
    ScalarFunction(const Polynomial<Real>& p);

    ScalarFunction(ScalarFunctionInterface* fptr) : _ptr(fptr) { }
    const ScalarFunctionInterface* pointer() const { return this->_ptr.operator->(); }
    operator const ScalarFunctionInterface& () const { return *this->_ptr; }

    Nat argument_size() const { return this->_ptr->argument_size(); }
    Float evaluate(const Vector<Float>& x) const { return this->_ptr->evaluate(x); }
    Interval evaluate(const Vector<Interval>& x) const { return this->_ptr->evaluate(x); }
    TaylorModel evaluate(const Vector<TaylorModel>& x) const { return this->_ptr->evaluate(x); }

    Differential<Float> evaluate(const Vector< Differential<Float> >& x) const { return this->_ptr->evaluate(x); }
    Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const { return this->_ptr->evaluate(x); }

    Float operator() (const Vector<Float>& x) const { return this->_ptr->evaluate(x); }
    Interval operator() (const Vector<Interval>& x) const { return this->_ptr->evaluate(x); }

    Vector<Float> gradient(const Vector<Float>& x) const;
    Vector<Interval> gradient(const Vector<Interval>& x) const;

    ScalarFunction derivative(Nat j) const;
    Polynomial<Real> polynomial() const;

    std::ostream& write(std::ostream& os) const;
  public:
    friend ScalarFunction compose(const ScalarFunction&, const VectorFunction&);
    friend Real evaluate(const ScalarFunction&, const Vector<Real>&);

    friend ScalarFunction operator+(const ScalarFunction&);
    friend ScalarFunction operator-(const ScalarFunction&);
    friend ScalarFunction operator+(const ScalarFunction&, const ScalarFunction&);
    friend ScalarFunction operator-(const ScalarFunction&, const ScalarFunction&);
    friend ScalarFunction operator*(const ScalarFunction&, const ScalarFunction&);
    friend ScalarFunction operator/(const ScalarFunction&, const ScalarFunction&);
  public:
    const ScalarFunctionInterface* _raw_pointer() const { return this->_ptr.operator->(); }
  private:
    shared_ptr<ScalarFunctionInterface> _ptr;
};

inline Float evaluate_approx(const ScalarFunction& f, const Vector<Float>& x) { return f(x); }
inline Interval evaluate(const ScalarFunction& f, const Vector<Interval>& x) { return f(x); }
inline TaylorModel compose(const ScalarFunction& f, const Vector<TaylorModel>& x) { return f.evaluate(x); }
inline Vector<Float> gradient_approx(const ScalarFunction& f, const Vector<Float>& x) { return f.gradient(x); }
inline Vector<Interval> gradient(const ScalarFunction& f, const Vector<Interval>& x) { return f.gradient(x); }
inline std::ostream& operator<<(std::ostream& os, const ScalarFunction& f) { return f.write(os); }

ScalarFunction embed(const ScalarFunction&, uint i);

ScalarFunction operator+(const ScalarFunction&);
ScalarFunction operator-(const ScalarFunction&);
ScalarFunction operator+(const ScalarFunction&, const ScalarFunction&);
ScalarFunction operator-(const ScalarFunction&, const ScalarFunction&);
ScalarFunction operator*(const ScalarFunction&, const ScalarFunction&);
ScalarFunction operator/(const ScalarFunction&, const ScalarFunction&);
ScalarFunction operator+(const ScalarFunction&, const Real&);
ScalarFunction operator-(const ScalarFunction&, const Real&);
ScalarFunction operator*(const ScalarFunction&, const Real&);
ScalarFunction operator/(const ScalarFunction&, const Real&);
ScalarFunction operator+(const Real&, const ScalarFunction&);
ScalarFunction operator-(const Real&, const ScalarFunction&);
ScalarFunction operator*(const Real&, const ScalarFunction&);
ScalarFunction operator/(const Real&, const ScalarFunction&);

ScalarFunction pow(const ScalarFunction&, int);
ScalarFunction rec(const ScalarFunction&);
ScalarFunction sqr(const ScalarFunction&);
ScalarFunction sqrt(const ScalarFunction&);
ScalarFunction exp(const ScalarFunction&);
ScalarFunction log(const ScalarFunction&);
ScalarFunction sin(const ScalarFunction&);
ScalarFunction cos(const ScalarFunction&);
ScalarFunction tan(const ScalarFunction&);

//! A vector function \f$f:\R^n\rightarrow\R^m\f$.
class VectorFunction
{
    typedef uint Nat;
    typedef uint SizeType;
    typedef ushort SmoothnessType;
    typedef std::ostream OStream;
  public:
    static VectorFunction constant(const Vector<Real>& c, Nat as);
    static VectorFunction identity(Nat n);

    VectorFunction();
    VectorFunction(Nat rs, Nat as);
    VectorFunction(Nat rs, const ScalarFunction& sf);
    VectorFunction(VectorFunctionInterface*);
    const VectorFunctionInterface* pointer() const { return this->_ptr.operator->(); }

    VectorFunction(const List< Expression<Real> >& e, const Space<Real>& s);
    VectorFunction(const Space<Real>& rs, const Map<RealVariable,RealExpression>& e, const Space<Real>& as);

    ScalarFunction operator[](Nat i) const;
    ScalarFunction& operator[](Nat i);

    ScalarFunction get(Nat) const;
    void set(Nat,ScalarFunction);

    VectorFunctionInterface* clone() const { return this->_ptr->clone(); }

    Nat result_size() const { return this->_ptr->result_size(); }
    Nat argument_size() const { return this->_ptr->argument_size(); }
    SmoothnessType smoothness() const { return this->_ptr->smoothness(); }

    Vector<Float> evaluate(const Vector<Float>& x) const { return this->_ptr->evaluate(x); }
    Vector<Interval> evaluate(const Vector<Interval>& x) const { return this->_ptr->evaluate(x); }
    Vector<TaylorModel> evaluate(const Vector<TaylorModel>& x) const { return this->_ptr->evaluate(x); }

    Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const { return this->_ptr->evaluate(x); }
    Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const { return this->_ptr->evaluate(x); }

    Matrix<Float> jacobian(const Vector<Float>& x) const;
    Matrix<Interval> jacobian(const Vector<Interval>& x) const;

    Vector<Float> operator()(const Vector<Float>& x) const;
    Vector<Interval> operator()(const Vector<Interval>& x) const;

    std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }
  public:
    friend VectorFunction join(const ScalarFunction&, const ScalarFunction&);
    friend VectorFunction join(const ScalarFunction&, const VectorFunction&);
    friend VectorFunction join(const VectorFunction&, const ScalarFunction&);
    friend VectorFunction join(const VectorFunction&, const VectorFunction&);

    friend VectorFunction compose(const VectorFunction&, const VectorFunction&);
    friend Vector<Real> evaluate(const VectorFunction&, const Vector<Real>&);

    friend VectorFunction operator+(const VectorFunction&, const VectorFunction&);
    friend VectorFunction operator-(const VectorFunction&, const VectorFunction&);
    friend VectorFunction operator*(const VectorFunction&, const ScalarFunction&);
    friend VectorFunction operator*(const ScalarFunction&, const VectorFunction&);
  public:
    const VectorFunctionInterface* _raw_pointer() const { return this->_ptr.operator->(); }
  private:
    shared_ptr<VectorFunctionInterface> _ptr;
};

inline Vector<Float> evaluate_approx(const VectorFunction& f, const Vector<Float>& x) { return f(x); }
inline Vector<Interval> evaluate(const VectorFunction& f, const Vector<Interval>& x) { return f(x); }
inline Vector<TaylorModel> compose(const VectorFunction& f, const Vector<TaylorModel>& x) { return f.evaluate(x); }
inline Matrix<Float> jacobian_approx(const VectorFunction& f, const Vector<Float>& x);
inline Matrix<Interval> jacobian(const VectorFunction& f, const Vector<Interval>& x);

VectorFunction operator*(const ScalarFunction& sf, const Vector<Real>& e);
VectorFunction operator+(const VectorFunction& f1, const VectorFunction& f2);
VectorFunction operator-(const VectorFunction& f1, const VectorFunction& f2);
VectorFunction operator*(const VectorFunction& vf, const ScalarFunction& sf);
VectorFunction operator*(const ScalarFunction& sf, const VectorFunction& vf);

VectorFunction join(const ScalarFunction&, const ScalarFunction&);
VectorFunction join(const ScalarFunction&, const VectorFunction&);
VectorFunction join(const VectorFunction&, const ScalarFunction&);
VectorFunction join(const VectorFunction&, const VectorFunction&);

ScalarFunction compose(const ScalarFunction& f, const VectorFunction& g);
VectorFunction compose(const VectorFunction& f, const VectorFunction& g);
ScalarFunction lie_derivative(const ScalarFunction& g, const VectorFunction& f);

inline std::ostream& operator<<(std::ostream& os, const VectorFunction& f) { return f.write(os); }


class IdentityFunction
    : public VectorFunction
{
  public:
    //! \brief Construct the identity function in dimension \a n.
    IdentityFunction(uint n);
};


class ProjectionFunction
    : public VectorFunction
{
  public:
    //! \brief Construct the identity function in dimension \a n.
    ProjectionFunction(uint n);
    //! \brief Construct the projection functions \f$f_i(x)=x_{i+k}\f$ for \f$i=0,\ldots,m-1\f$. Precondition: \f$m+k\leq n\f$.
    ProjectionFunction(uint m, uint n, uint k);
    //! \brief Construct the projection function  with \f$f_i(x)=x_{p_i}\f$ for \f$i=0,\ldots,m-1\f$.
    ProjectionFunction(uint m, uint n, const array<uint>& p);
    //! \brief Construct the projection function with \f$f_i(x)=x_{p_i}\f$ for \f$i=0,\ldots,|p|-1\f$.
    ProjectionFunction(const array<uint>& p, uint n);

    const array<uint>& p() const;
    const uint p(uint i) const;
};


class VectorConstantFunction
    : public VectorFunction
{
  public:
    VectorConstantFunction(const Vector<Real>& c, uint as);
};


class VectorAffineFunction
    : public VectorFunction
{
  public:
    VectorAffineFunction(const Matrix<Real>& A, const Vector<Real>& b);
    const Matrix<Real> A() const;
    const Vector<Real> b() const;
};


} // namespace Ariadne

#endif
