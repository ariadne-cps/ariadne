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
#include "declarations.h"

namespace Ariadne {

static const int SMOOTH=255;

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
class ScalarFunctionInterface<ApproximateTag>
    : public ScalarFunctionInterface<Void>
{
  public:
    //! \brief Return a copy of the function.
    inline ScalarFunction<ApproximateTag> clone() const;
    //! \brief Compute an approximation to the value of the function at the point \a x.
    virtual ApproximateNumber evaluate(const Vector<ApproximateNumber>& x) const = 0;
    inline ApproximateNumber operator() (const Vector<ApproximateNumber>& x) const;
    //! \brief Evaluate the function over a vector of differentials.
    virtual Differential<ApproximateNumber> evaluate(const Vector< Differential<ApproximateNumber> >& x) const = 0;
    //! \brief Evaluate the function over a vector of approximate Taylor models.
    virtual TaylorModel<ApproximateNumber> evaluate(const Vector< TaylorModel<ApproximateNumber> >& x) const = 0;
    //! \brief Evaluate the function over a vector of formulae to create the composed function.
    virtual Formula<ApproximateNumber> evaluate(const Vector< Formula<ApproximateNumber> >& x) const = 0;
    //! \brief Evaluate the function over a vector of elements of an algebra.
    virtual Algebra<ApproximateNumber> evaluate(const Vector< Algebra<ApproximateNumber> >& x) const = 0;

    Vector<ApproximateNumber> gradient(const Vector<ApproximateNumber>& x) const;
    Differential<ApproximateNumber> differential(const Vector<ApproximateNumber>& x, Nat d) const;

    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    inline ScalarFunction<ApproximateTag> derivative(Nat i) const;
  private:
    virtual ScalarFunctionInterface<ApproximateTag>* _derivative(Nat i) const = 0;
  public:
    virtual ScalarFunctionInterface<ApproximateTag>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\mathbb{I}^n\rightarrow\mathbb{I}\f$ which can be evaluated over intervals.
//! \sa \ref VectorFunctionInterface.
template<>
class ScalarFunctionInterface<ValidatedTag>
    : public ScalarFunctionInterface<ApproximateTag>
{
  public:

    inline ScalarFunction<ValidatedTag> clone() const;

    using ScalarFunctionInterface<ApproximateTag>::evaluate;

    //! \brief Compute an over-approximation to the values of the function over the domain \a x. This method provides an <em>interval extension</em> of the function.
    virtual ValidatedNumber evaluate(const Vector<ValidatedNumber>& x) const = 0;
    inline ValidatedNumber evaluate(const Vector<ExactNumber>& x) const;
    inline ValidatedNumber operator() (const Vector<ValidatedNumber>& x) const;
    //! \brief Evaluate the function over a vector of interval differentials.
    virtual Differential<ValidatedNumber> evaluate(const Vector< Differential<ValidatedNumber> >& x) const = 0;
    //! \brief Evaluate the function over a vector of Taylor models with interval error.
    virtual TaylorModel<ValidatedNumber> evaluate(const Vector< TaylorModel<ValidatedNumber> >& x) const = 0;

    //! \brief Apply the function to a formula. Can be used to obtain a tree structure from the function.
    virtual Formula<ValidatedNumber> evaluate(const Vector< Formula<ValidatedNumber> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Algebra<ValidatedNumber> evaluate(const Vector< Algebra<ValidatedNumber> >& x) const = 0;

    using ScalarFunctionInterface<ApproximateTag>::gradient;
    Vector<ValidatedNumber> gradient(const Vector<ValidatedNumber>& x) const;
    using ScalarFunctionInterface<ApproximateTag>::differential;
    Differential<ValidatedNumber> differential(const Vector<ValidatedNumber>& x, Nat d) const;

    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    inline ScalarFunction<ValidatedTag> derivative(Nat i) const;
  private:
    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    virtual ScalarFunctionInterface<ValidatedTag>* _derivative(Nat i) const = 0;
  public:
    virtual ScalarFunctionInterface<ValidatedTag>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\R^n\rightarrow\R\f$ which can be evaluated exactly.
//! \sa \ref VectorFunctionInterface.
template<>
class ScalarFunctionInterface<EffectiveTag>
    : public ScalarFunctionInterface<ValidatedTag>
{
  public:
    inline ScalarFunction<EffectiveTag> clone() const;

    using ScalarFunctionInterface<ValidatedTag>::evaluate;

    //! \brief Evaluate over computable reals.
    virtual EffectiveNumber evaluate(const Vector<EffectiveNumber>& x) const = 0;
    inline EffectiveNumber operator() (const Vector<EffectiveNumber>& x) const;
    //! \brief Apply the function to a formula. Can be used to obtain a tree structure from the function.
    virtual Formula<EffectiveNumber> evaluate(const Vector< Formula<EffectiveNumber> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Algebra<EffectiveNumber> evaluate(const Vector< Algebra<EffectiveNumber> >& x) const = 0;

    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    inline ScalarFunction<EffectiveTag> derivative(Nat i) const;
  private:
    //! \brief The derivative with respect to the \a j<sup>th</sup> coordinate.
    virtual ScalarFunctionInterface<EffectiveTag>* _derivative(Nat i) const = 0;
  public:
    virtual ScalarFunctionInterface<EffectiveTag>* _clone() const = 0;
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
class VectorFunctionInterface<ApproximateTag>
    : public VectorFunctionInterface<Void>
{
  public:
    //! \brief Compute an approximation to the value of the function at the point \a x.
    virtual Vector<ApproximateNumber> evaluate(const Vector<ApproximateNumber>& x) const = 0;
    //! \brief Evaluate the function over a vector of differentials.
    virtual Vector< Differential<ApproximateNumber> > evaluate(const Vector< Differential<ApproximateNumber> >& x) const = 0;
    //! \brief Evaluate the function over a vector of approximate Taylor models.
    virtual Vector< TaylorModel<ApproximateNumber> > evaluate(const Vector< TaylorModel<ApproximateNumber> >& x) const = 0;
    //! \brief Evaluate the function over a vector of formulae to create the composed function.
    virtual Vector< Formula<ApproximateNumber> > evaluate(const Vector< Formula<ApproximateNumber> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Vector< Algebra<ApproximateNumber> > evaluate(const Vector< Algebra<ApproximateNumber> >& x) const = 0;

    Matrix<ApproximateNumber> jacobian(const Vector<ApproximateNumber>& x) const;
    Vector< Differential<ApproximateNumber> > differentials(const Vector<ApproximateNumber>& x, Nat d) const;

    //! \brief Get the \a i<sup>th</sup> component function.
    inline ScalarFunction<ApproximateTag> operator[](Nat i) const;
  public:
    virtual ScalarFunctionInterface<ApproximateTag>* _get(Nat i) const = 0;
    virtual VectorFunctionInterface<ApproximateTag>* _clone() const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\I^n\rightarrow\I^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<ValidatedTag>
    : public VectorFunctionInterface<ApproximateTag>
{
  public:
    using VectorFunctionInterface<ApproximateTag>::evaluate;

    //! \brief Compute an over-approximation to the values of the function over the domain \a x. This method provides an <em>interval extension</em> of the function.
    virtual Vector<ValidatedNumber> evaluate(const Vector<ValidatedNumber>& x) const = 0;
    inline Vector<ValidatedNumber> evaluate(const Vector<ExactNumber>& x) const;
    //! \brief Evaluate the function over a vector of interval differentials.
    virtual Vector< Differential<ValidatedNumber> > evaluate(const Vector< Differential<ValidatedNumber> >& x) const = 0;
    //! \brief Evaluate the function over a vector of Taylor models with interval error.
    virtual Vector< TaylorModel<ValidatedNumber> > evaluate(const Vector< TaylorModel<ValidatedNumber> >& x) const = 0;

    //! \brief Evaluate the function over a vector of formulae.
    virtual Vector< Formula<ValidatedNumber> > evaluate(const Vector< Formula<ValidatedNumber> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Vector< Algebra<ValidatedNumber> > evaluate(const Vector< Algebra<ValidatedNumber> >& x) const = 0;

    using VectorFunctionInterface<ApproximateTag>::jacobian;
    Matrix<ValidatedNumber> jacobian(const Vector<ValidatedNumber>& x) const;
    Matrix<ValidatedNumber> jacobian(const Vector<ExactNumber>& x) const;
    using VectorFunctionInterface<ApproximateTag>::differentials;
    Vector< Differential<ValidatedNumber> > differentials(const Vector<ValidatedNumber>& x, Nat d) const;

    //! \brief Get the \a i<sup>th</sup> component function.
    inline ScalarFunction<ValidatedTag> operator[](Nat i) const;
  public:
    virtual ScalarFunctionInterface<ValidatedTag>* _get(Nat i) const = 0;
    virtual VectorFunctionInterface<ValidatedTag>* _clone() const = 0;

};

//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\R^n\rightarrow\R^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<>
class VectorFunctionInterface<EffectiveTag>
    : public VectorFunctionInterface<ValidatedTag>
{
  public:
    using VectorFunctionInterface<ValidatedTag>::evaluate;

    //! \brief Evaluate over computable reals.
    virtual Vector<EffectiveNumber> evaluate(const Vector<EffectiveNumber>& x) const = 0;
    //! \brief Evaluate the function over a vector of formulae.
    virtual Vector< Formula<EffectiveNumber> > evaluate(const Vector< Formula<EffectiveNumber> >& x) const = 0;
    //! \brief Apply the function to an algebra.
    virtual Vector< Algebra<EffectiveNumber> > evaluate(const Vector< Algebra<EffectiveNumber> >& x) const = 0;

    //! \brief Get the \a i<sup>th</sup> component function.
    inline ScalarFunction<EffectiveTag> operator[](Nat i) const;
  public:
    virtual ScalarFunctionInterface<EffectiveTag>* _get(Nat i) const = 0;
    virtual VectorFunctionInterface<EffectiveTag>* _clone() const = 0;
};

//! \relates VectorFunctionInterface
//! \brief Write to an output stream. Calls the write(std::ostream&) method to perform dynamic dispatching.
inline std::ostream& operator<<(std::ostream& os, const VectorFunctionInterface<Void>& f) {
    return f.write(os);
}



//! \brief An interface for scalar function models on a restricted domain.
class ScalarModelInterface {
    virtual ExactBox domain() const = 0;
    virtual ValidatedNumber evaluate(const Vector<ValidatedNumber>&) const = 0;
};

template<class X> class FunctionFactoryInterface;

template<> class FunctionFactoryInterface<ValidatedTag>
{
    typedef ExactBox DomainType;
  public:
    virtual FunctionFactoryInterface<ValidatedTag>* clone() const = 0;
    virtual Void write(OutputStream& os) const = 0;
    inline ScalarFunction<ValidatedTag> create(const ExactBox& domain, const ScalarFunctionInterface<ValidatedTag>& function) const;
    inline VectorFunction<ValidatedTag> create(const ExactBox& domain, const VectorFunctionInterface<ValidatedTag>& function) const;
    inline ScalarFunction<ValidatedTag> create_zero(const ExactBox& domain) const;
    inline VectorFunction<ValidatedTag> create_identity(const ExactBox& domain) const;
  private:
    virtual ScalarFunctionInterface<ValidatedTag>* _create(const ExactBox& domain, const ScalarFunctionInterface<ValidatedTag>& function) const = 0;
    virtual VectorFunctionInterface<ValidatedTag>* _create(const ExactBox& domain, const VectorFunctionInterface<ValidatedTag>& function) const = 0;
};

template<class X> inline OutputStream& operator<<(OutputStream& os, const FunctionFactoryInterface<ValidatedTag>& factory) {
    factory.write(os); return os;
}

} // namespace Ariadne

#endif
