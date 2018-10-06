/***************************************************************************
 *            function_interface.hpp
 *
 *  Copyright 2008-17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file function_interface.hpp
 *  \brief Interface for functions for which derivatives can be computed.
 */

#ifndef ARIADNE_FUNCTION_INTERFACE_HPP
#define ARIADNE_FUNCTION_INTERFACE_HPP

#include <iosfwd>

#include "../utility/declarations.hpp"
#include "../function/function.decl.hpp"

namespace Ariadne {

static const Int SMOOTH=255;

template<class S> struct ElementTraits;
template<class S, class X> using ElementType = typename ElementTraits<S>::template Type<X>;
template<class S> using ElementSizeType = decltype(declval<S>().dimension());
template<class S> using ElementIndexType = decltype(declval<S>().dimension());

template<class P, class D, class C> class FunctionInterface;


//! \ingroup FunctionModule
//! \brief Interface for vector functions \f$\F^n\rightarrow\F^m\f$ whose derivatives can be computed.
//! \sa \ref ScalarFunctionInterface
template<class P,class D,class C=BoxDomainType>
class VectorOfFunctionInterface
{
};

template<class D> class VectorOfFunctionInterface<ApproximateTag,D,BoxDomainType> {
  public:
    virtual ~VectorOfFunctionInterface<ApproximateTag,D>() = default;
    virtual ScalarFunctionInterface<ApproximateTag,D>* _get(SizeType i) const = 0;
};

template<class D> class VectorOfFunctionInterface<ValidatedTag,D,BoxDomainType>
    : public virtual VectorOfFunctionInterface<ApproximateTag,D>
{
  public:
    virtual ScalarFunctionInterface<ValidatedTag,D>* _get(SizeType i) const override = 0;
};

template<class D> class VectorOfFunctionInterface<EffectiveTag,D,BoxDomainType>
    : public virtual VectorOfFunctionInterface<ValidatedTag,D>
{
  public:
    virtual ScalarFunctionInterface<EffectiveTag,D>* _get(SizeType i) const override = 0;
};


template<class D, class C>
class FunctionInterface<Void,D,C>
{
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    typedef D DomainType;
    typedef C CodomainType;
    typedef ElementSizeType<DomainType> ArgumentSizeType;
    typedef ElementSizeType<CodomainType> ResultSizeType;
    typedef ElementIndexType<DomainType> ArgumentIndexType;

    virtual ~FunctionInterface() = default;
    virtual ArgumentSizeType argument_size() const = 0;
    virtual ResultSizeType result_size() const = 0;
    virtual DomainType const domain() const = 0;
    virtual CodomainType const codomain() const = 0;

    virtual OutputStream& repr(OutputStream& os) const = 0;
    virtual OutputStream& write(OutputStream& os) const = 0;
  public:
    virtual FunctionInterface<Void,D,C>* _clone() const = 0;
  public:
    friend inline OutputStream& operator<<(OutputStream& os, const FunctionInterface<Void,D,C>& f) { return f.write(os); }
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\mathbb{F}^n\rightarrow\mathbb{F}\f$ which can only be evaluated approximately.
//! \sa \ref VectorFunctionInterface.
template<class D, class C>
class FunctionInterface<ApproximateTag,D,C>
    : public virtual FunctionInterface<Void,D,C>
{
    typedef ApproximateTag P;
  public:
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    virtual Result<FloatDPApproximation> _evaluate(const Argument<FloatDPApproximation>& x) const = 0;
    virtual Result<FloatMPApproximation> _evaluate(const Argument<FloatMPApproximation>& x) const = 0;
    virtual Result<Differential<FloatDPApproximation>> _evaluate(const Argument< Differential<FloatDPApproximation> >& x) const = 0;
    virtual Result<Differential<FloatMPApproximation>> _evaluate(const Argument< Differential<FloatMPApproximation> >& x) const = 0;
    virtual Result<TaylorModel<ApproximateTag,FloatDP>> _evaluate(const Argument< TaylorModel<ApproximateTag,FloatDP> >& x) const = 0;
    virtual Result<TaylorModel<ApproximateTag,FloatMP>> _evaluate(const Argument< TaylorModel<ApproximateTag,FloatMP> >& x) const = 0;
    virtual Result<Formula<ApproximateNumber>> _evaluate(const Argument< Formula<ApproximateNumber> >& x) const = 0;
    virtual Result<Algebra<ApproximateNumber>> _evaluate(const Argument< Algebra<ApproximateNumber> >& x) const = 0;

    virtual FunctionInterface<P,D,C>* _clone() const = 0;
    virtual FunctionInterface<P,D,C>* _derivative(ElementIndexType<D> i) const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\mathbb{I}^n\rightarrow\mathbb{I}\f$ which can be evaluated over intervals.
//! \sa \ref VectorFunctionInterface.
template<class D, class C>
class FunctionInterface<ValidatedTag,D,C>
    : public virtual FunctionInterface<ApproximateTag,D,C>
{
    typedef ValidatedTag P;
    typedef ApproximateTag AP;
  public:
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    using FunctionInterface<AP,D,C>::_evaluate;
    virtual Result<FloatDPBounds> _evaluate(const Argument<FloatDPBounds>& x) const = 0;
    virtual Result<FloatMPBounds> _evaluate(const Argument<FloatMPBounds>& x) const = 0;
    virtual Result<Differential<FloatDPBounds>> _evaluate(const Argument< Differential<FloatDPBounds> >& x) const = 0;
    virtual Result<Differential<FloatMPBounds>> _evaluate(const Argument< Differential<FloatMPBounds> >& x) const = 0;
    virtual Result<TaylorModel<ValidatedTag,FloatDP>> _evaluate(const Argument< TaylorModel<ValidatedTag,FloatDP> >& x) const = 0;
    virtual Result<TaylorModel<ValidatedTag,FloatMP>> _evaluate(const Argument< TaylorModel<ValidatedTag,FloatMP> >& x) const = 0;
    virtual Result<Formula<ValidatedNumber>> _evaluate(const Argument< Formula<ValidatedNumber> >& x) const = 0;
    virtual Result<Algebra<ValidatedNumber>> _evaluate(const Argument< Algebra<ValidatedNumber> >& x) const = 0;

    virtual Result<ScalarFunction<ValidatedTag>> _evaluate(const Argument< ScalarFunction<ValidatedTag> >& x) const = 0;

    inline Result<FloatDPBounds> _evaluate(const Argument<FloatDPValue>& x) const {
        return this->_evaluate(Argument<FloatDPBounds>(x)); }
    inline Result<FloatMPBounds> _evaluate(const Argument<FloatMPValue>& x) const {
        return this->_evaluate(Argument<FloatMPBounds>(x)); }

    virtual FunctionInterface<P,D,C>* _clone() const = 0;
    virtual FunctionInterface<P,D,C>* _derivative(ElementIndexType<D> i) const = 0;
};

//! \ingroup FunctionModule
//! \brief Interface for scalar functions \f$\R^n\rightarrow\R\f$ which can be evaluated exactly.
//! \sa \ref VectorFunctionInterface.
template<class D, class C>
class FunctionInterface<EffectiveTag,D,C>
    : public virtual FunctionInterface<ValidatedTag,D,C>
{
    typedef ValidatedTag WP;
    typedef EffectiveTag P;
  public:
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    using FunctionInterface<WP,D,C>::_evaluate;
    virtual Result<Real> _evaluate(const Argument<Real>& x) const = 0;
    virtual Result<Algebra<Real>> _evaluate(const Argument<Algebra<Real>>& x) const = 0;
    virtual Result<Formula<Real>> _evaluate(const Argument<Formula<Real>>& x) const = 0;
    virtual Result<Algebra<EffectiveNumber>> _evaluate(const Argument<Algebra<EffectiveNumber>>& x) const = 0;
    virtual Result<Formula<EffectiveNumber>> _evaluate(const Argument<Formula<EffectiveNumber>>& x) const = 0;

    virtual FunctionInterface<P,D,C>* _clone() const = 0;
    virtual FunctionInterface<P,D,C>* _derivative(ElementIndexType<D> i) const = 0;
};


template<class I> class VectorInterface {
    ~VectorInterface() = default;
    virtual I* _get(SizeType i) const = 0;
};

template<class X> class FunctionFactoryInterface;

template<> class FunctionFactoryInterface<ValidatedTag>
{
    typedef BoxDomainType DomainType;
  public:
    virtual FunctionFactoryInterface<ValidatedTag>* clone() const = 0;
    virtual OutputStream& write(OutputStream& os) const = 0;
    inline ValidatedScalarFunction create(const BoxDomainType& domain, const ScalarFunctionInterface<ValidatedTag>& function) const;
    inline ValidatedVectorFunction create(const BoxDomainType& domain, const VectorFunctionInterface<ValidatedTag>& function) const;
    inline ValidatedScalarFunction create_zero(const BoxDomainType& domain) const;
    inline ValidatedVectorFunction create_identity(const BoxDomainType& domain) const;
  private:
    virtual ScalarFunctionInterface<ValidatedTag>* _create(const BoxDomainType& domain, const ScalarFunctionInterface<ValidatedTag>& function) const = 0;
    virtual VectorFunctionInterface<ValidatedTag>* _create(const BoxDomainType& domain, const VectorFunctionInterface<ValidatedTag>& function) const = 0;
  public:
    friend inline OutputStream& operator<<(OutputStream& os, const FunctionFactoryInterface<ValidatedTag>& factory) {
        return factory.write(os); }
  public:
    virtual ~FunctionFactoryInterface<ValidatedTag>() = default;
};


} // namespace Ariadne

#endif
