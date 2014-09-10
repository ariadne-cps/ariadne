/***************************************************************************
 *            function_interface.h
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

/*! \file function_interface.h
 *  \brief Interface for functions for which derivatives can be computed.
 */
#ifndef ARIADNE_FUNCTION_INTERFACE_H
#define ARIADNE_FUNCTION_INTERFACE_H

#include <iosfwd>

namespace Ariadne {

typedef std::ostream OutputStream;

static const int SMOOTH=255;

typedef void Void;
typedef unsigned int Nat;
typedef int Int;

class Float;
class Interval;
class Real;

class Interval;
class Box;

template<class X> class Vector;
template<class X> class Matrix;
template<class X> class Differential;
template<class X> class TaylorModel;
template<class X> class Formula;
template<class X> class Algebra;

typedef Vector<Float> FloatVector;
typedef Vector<Interval> IntervalVector;

template<class X> class ScalarFunctionInterface;
template<class X> class VectorFunctionInterface;

typedef ScalarFunctionInterface<Float> FloatScalarFunctionInterface;
typedef ScalarFunctionInterface<Interval> IntervalScalarFunctionInterface;
typedef ScalarFunctionInterface<Real> RealScalarFunctionInterface;

typedef VectorFunctionInterface<Float> FloatVectorFunctionInterface;
typedef VectorFunctionInterface<Interval> IntervalVectorFunctionInterface;
typedef VectorFunctionInterface<Real> RealVectorFunctionInterface;

template<class X> class ScalarFunction;
typedef ScalarFunction<Float> FloatScalarFunction;
typedef ScalarFunction<Interval> IntervalScalarFunction;
typedef ScalarFunction<Real> RealScalarFunction;

template<class X> class VectorFunction;
typedef VectorFunction<Float> FloatVectorFunction;
typedef VectorFunction<Interval> IntervalVectorFunction;
typedef VectorFunction<Real> RealVectorFunction;

template<>
class ScalarFunctionInterface<Void>
{
  public:
    //! \brief The type used to describe the number of argument variables.
    typedef Nat SizeType;

    //! \brief Virtual destructor.
    virtual ~ScalarFunctionInterface() { };
    //! \brief The number of arguments to the expression.
    virtual SizeType argument_size() const = 0;

    //! \brief Write a full version to an output stream.
    virtual std::ostream& repr(std::ostream& os) const = 0;
    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const = 0;
  public:
    virtual ScalarFunctionInterface<Void>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\mathbb{F}^n\rightarrow\mathbb{F}\f$ which can only be evaluated approximately.
//! \sa \ref VectorFunctionInterface.
template<>
class ScalarFunctionInterface<Float>
    : public ScalarFunctionInterface<Void>
{
  public:
    //! \brief Return a copy of the function.
    inline ScalarFunction<Float> clone() const;
    //! \brief Compute an approximation to the value of the function at the point \a x.
    virtual Float evaluate(const Vector<Float>& x) const = 0;
    inline Float operator() (const Vector<Float>& x) const;
    //! \brief Evaluate the function over a vector of differentials.
    virtual Differential<Float> evaluate(const Vector< Differential<Float> >& x) const = 0;
    //! \brief Evaluate the function over a vector of approximate Taylor models.
    virtual TaylorModel<Float> evaluate(const Vector< TaylorModel<Float> >& x) const = 0;
    //! \brief Evaluate the function over a vector of formulae to create the composed function.
    virtual Formula<Float> evaluate(const Vector< Formula<Float> >& x) const = 0;
    //! \brief Evaluate the function over a vector of elements of an algebra.
    virtual Algebra<Float> evaluate(const Vector< Algebra<Float> >& x) const = 0;

    Vector<Float> gradient(const Vector<Float>& x) const;
    Differential<Float> differential(const Vector<Float>& x, Nat d) const;

    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    inline ScalarFunction<Float> derivative(Nat i) const;
  private:
    virtual ScalarFunctionInterface<Float>* _derivative(Nat i) const = 0;
  public:
    virtual ScalarFunctionInterface<Float>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\mathbb{I}^n\rightarrow\mathbb{I}\f$ which can be evaluated over intervals.
//! \sa \ref VectorFunctionInterface.
template<>
class ScalarFunctionInterface<Interval>
    : public ScalarFunctionInterface<Float>
{
  public:

    inline ScalarFunction<Interval> clone() const;

    using ScalarFunctionInterface<Float>::evaluate;

    //! \brief Compute an over-approximation to the values of the function over the domain \a x. This method provides an <em>interval extension</em> of the function.
    virtual Interval evaluate(const Vector<Interval>& x) const = 0;
    inline Interval operator() (const Vector<Interval>& x) const;
    //! \brief Evaluate the function over a vector of interval differentials.
    virtual Differential<Interval> evaluate(const Vector< Differential<Interval> >& x) const = 0;
    //! \brief Evaluate the function over a vector of Taylor models with interval error.
    virtual TaylorModel<Interval> evaluate(const Vector< TaylorModel<Interval> >& x) const = 0;

    //! \brief Apply the function to a formula. Can be used to obtain a tree structure from the function.
    virtual Formula<Interval> evaluate(const Vector< Formula<Interval> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Algebra<Interval> evaluate(const Vector< Algebra<Interval> >& x) const = 0;

    using ScalarFunctionInterface<Float>::gradient;
    Vector<Interval> gradient(const Vector<Interval>& x) const;
    using ScalarFunctionInterface<Float>::differential;
    Differential<Interval> differential(const Vector<Interval>& x, Nat d) const;

    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    inline ScalarFunction<Interval> derivative(Nat i) const;
  private:
    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    virtual ScalarFunctionInterface<Interval>* _derivative(Nat i) const = 0;
  public:
    virtual ScalarFunctionInterface<Interval>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\R^n\rightarrow\R\f$ which can be evaluated exactly.
//! \sa \ref VectorFunctionInterface.
template<>
class ScalarFunctionInterface<Real>
    : public ScalarFunctionInterface<Interval>
{
  public:
    inline ScalarFunction<Real> clone() const;

    using ScalarFunctionInterface<Interval>::evaluate;

    //! \brief Evaluate over computable reals.
    virtual Real evaluate(const Vector<Real>& x) const = 0;
    inline Real operator() (const Vector<Real>& x) const;
    //! \brief Apply the function to a formula. Can be used to obtain a tree structure from the function.
    virtual Formula<Real> evaluate(const Vector< Formula<Real> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Algebra<Real> evaluate(const Vector< Algebra<Real> >& x) const = 0;

    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    inline ScalarFunction<Real> derivative(Nat i) const;
  private:
    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    virtual ScalarFunctionInterface<Real>* _derivative(Nat i) const = 0;
  public:
    virtual ScalarFunctionInterface<Real>* _clone() const = 0;
};

//! \relates ScalarFunctionInterface
//! \brief Write to an output stream. Calls the write(std::ostream&) method to perform dynamic dispatching.
inline std::ostream& operator<<(std::ostream& os, const ScalarFunctionInterface<Void>& f) {
    return f.write(os);
}


//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\F^n\rightarrow\F^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<Void>
{
  public:
    //! \brief The type used to describe the number of argument variables.
    typedef Nat SizeType;

    //! \brief Virtual destructor.
    virtual ~VectorFunctionInterface() { };

    //! \brief The number of arguments to the function.
    virtual SizeType argument_size() const = 0;
    //! \brief The number of result variables of the function.
    virtual SizeType result_size() const = 0;

    //! \brief Write a full version to an output stream.
    virtual std::ostream& repr(std::ostream& os) const = 0;
    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const = 0;
  public:
    virtual VectorFunctionInterface<Void>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\F^n\rightarrow\F^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<Float>
    : public VectorFunctionInterface<Void>
{
  public:
    //! \brief Compute an approximation to the value of the function at the point \a x.
    virtual Vector<Float> evaluate(const Vector<Float>& x) const = 0;
    //! \brief Evaluate the function over a vector of differentials.
    virtual Vector< Differential<Float> > evaluate(const Vector< Differential<Float> >& x) const = 0;
    //! \brief Evaluate the function over a vector of approximate Taylor models.
    virtual Vector< TaylorModel<Float> > evaluate(const Vector< TaylorModel<Float> >& x) const = 0;
    //! \brief Evaluate the function over a vector of formulae to create the composed function.
    virtual Vector< Formula<Float> > evaluate(const Vector< Formula<Float> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Vector< Algebra<Float> > evaluate(const Vector< Algebra<Float> >& x) const = 0;

    Matrix<Float> jacobian(const Vector<Float>& x) const;
    Vector< Differential<Float> > differentials(const Vector<Float>& x, Nat d) const;

    //! \brief Get the \a i<sup>th</sup> component function.
    inline ScalarFunction<Float> operator[](Nat i) const;
  public:
    virtual ScalarFunctionInterface<Float>* _get(Nat i) const = 0;
    virtual VectorFunctionInterface<Float>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\I^n\rightarrow\I^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<Interval>
    : public VectorFunctionInterface<Float>
{
  public:
    using VectorFunctionInterface<Float>::evaluate;

    //! \brief Compute an over-approximation to the values of the function over the domain \a x. This method provides an <em>interval extension</em> of the function.
    virtual Vector<Interval> evaluate(const Vector<Interval>& x) const = 0;
    //! \brief Evaluate the function over a vector of interval differentials.
    virtual Vector< Differential<Interval> > evaluate(const Vector< Differential<Interval> >& x) const = 0;
    //! \brief Evaluate the function over a vector of Taylor models with interval error.
    virtual Vector< TaylorModel<Interval> > evaluate(const Vector< TaylorModel<Interval> >& x) const = 0;

    //! \brief Evaluate the function over a vector of formulae.
    virtual Vector< Formula<Interval> > evaluate(const Vector< Formula<Interval> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Vector< Algebra<Interval> > evaluate(const Vector< Algebra<Interval> >& x) const = 0;

    using VectorFunctionInterface<Float>::jacobian;
    Matrix<Interval> jacobian(const Vector<Interval>& x) const;
    using VectorFunctionInterface<Float>::differentials;
    Vector< Differential<Interval> > differentials(const Vector<Interval>& x, Nat d) const;

    //! \brief Get the \a i<sup>th</sup> component function.
    inline ScalarFunction<Interval> operator[](Nat i) const;
  public:
    virtual ScalarFunctionInterface<Interval>* _get(Nat i) const = 0;
    virtual VectorFunctionInterface<Interval>* _clone() const = 0;

};

//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\R^n\rightarrow\R^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<Real>
    : public VectorFunctionInterface<Interval>
{
  public:
    using VectorFunctionInterface<Interval>::evaluate;

    //! \brief Evaluate over computable reals.
    virtual Vector<Real> evaluate(const Vector<Real>& x) const = 0;
    //! \brief Evaluate the function over a vector of formulae.
    virtual Vector< Formula<Real> > evaluate(const Vector< Formula<Real> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Vector< Algebra<Real> > evaluate(const Vector< Algebra<Real> >& x) const = 0;

    //! \brief Get the \a i<sup>th</sup> component function.
    inline ScalarFunction<Real> operator[](Nat i) const;
  public:
    virtual ScalarFunctionInterface<Real>* _get(Nat i) const = 0;
    virtual VectorFunctionInterface<Real>* _clone() const = 0;
};

//! \relates VectorFunctionInterface
//! \brief Write to an output stream. Calls the write(std::ostream&) method to perform dynamic dispatching.
inline std::ostream& operator<<(std::ostream& os, const VectorFunctionInterface<Void>& f) {
    return f.write(os);
}



//! \brief An interface for scalar function models on a restricted domain.
class ScalarModelInterface {
    virtual Box domain() const = 0;
    virtual Interval evaluate(const Vector<Interval>&) const = 0;
};

template<class X> class FunctionFactoryInterface;

template<> class FunctionFactoryInterface<Interval>
{
    typedef IntervalVector DomainType;
  public:
    virtual FunctionFactoryInterface<Interval>* clone() const = 0;
    virtual Void write(OutputStream& os) const = 0;
    inline ScalarFunction<Interval> create(const Box& domain, const ScalarFunctionInterface<Interval>& function) const;
    inline VectorFunction<Interval> create(const Box& domain, const VectorFunctionInterface<Interval>& function) const;
    inline ScalarFunction<Interval> create_zero(const Box& domain) const;
    inline VectorFunction<Interval> create_identity(const Box& domain) const;
  private:
    virtual ScalarFunctionInterface<Interval>* _create(const Box& domain, const ScalarFunctionInterface<Interval>& function) const = 0;
    virtual VectorFunctionInterface<Interval>* _create(const Box& domain, const VectorFunctionInterface<Interval>& function) const = 0;
};

template<class X> inline OutputStream& operator<<(OutputStream& os, const FunctionFactoryInterface<Interval>& factory) {
    factory.write(os); return os;
}

} // namespace Ariadne

#endif
