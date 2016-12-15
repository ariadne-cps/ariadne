/***************************************************************************
 *            function_mixin.h
 *
 *  Copyright 2008-10  Pieter Collins
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

#ifndef ARIADNE_FUNCTION_MIXIN_H
#define ARIADNE_FUNCTION_MIXIN_H

#include "function/function_interface.h"

// Adaptors for classes to conform to the Function interface.

namespace Ariadne {

typedef ApproximateNumericType ApproximateNumericType;
typedef ValidatedNumericType ValidatedNumericType;
typedef EffectiveNumericType EffectiveNumericType;
typedef Differential<ApproximateNumericType> ApproximateDifferential;
typedef Differential<ValidatedNumericType> ValidatedDifferential;
typedef UnivariateDifferential<ApproximateNumericType> ApproximateUnivariateDifferential;
typedef UnivariateDifferential<ValidatedNumericType> ValidatedUnivariateDifferential;
typedef TaylorModel<ApproximateTag,Float64> ApproximateTaylorModel;
typedef TaylorModel<ValidatedTag,Float64> ValidatedTaylorModel;
typedef Formula<ApproximateNumericType> ApproximateFormula;
typedef Formula<ValidatedNumericType> ValidatedFormula;
typedef Formula<EffectiveNumericType> EffectiveFormula;
typedef Algebra<ApproximateNumericType> ApproximateAlgebra;
typedef Algebra<ValidatedNumericType> ValidatedAlgebra;
typedef Algebra<EffectiveNumericType> EffectiveAlgebra;

template<class F, class P, class D, class C> class FunctionMixin { };
template<class F, class P, class D=BoxDomain> class ScalarFunctionMixin;
template<class F, class P, class D=BoxDomain> class VectorFunctionMixin;

template<class T> T* heap_copy(const T& t) { return new T(t); }

template<class D> D make_domain(SizeType d);
template<> inline IntervalDomain make_domain(SizeType d) { assert(d==1u); return IntervalDomain(-inf,+inf); }
template<> inline BoxDomain make_domain(SizeType d) { return BoxDomain(d,IntervalDomain(-inf,+inf)); }

template<class F, class D, class C>
class FunctionMixin<F,Void,D,C>
    : public virtual FunctionInterface<Void,D,C>
{
    typedef typename FunctionInterface<Void,D,C>::DomainType DomainType;
    typedef typename FunctionInterface<Void,D,C>::CodomainType CodomainType;
  protected:
    FunctionMixin() { }
    template<class X> ElementType<C,X> _base_evaluate(const ElementType<D,X>& x) const;
  public:
    virtual DomainType const domain() const override { return make_domain<D>(this->argument_size()); }
    virtual CodomainType const codomain() const override { return make_domain<C>(this->result_size()); }
    virtual SizeType argument_size() const override { return this->domain().dimension(); }
    virtual SizeType result_size() const override { return this->codomain().dimension(); }

    virtual OutputStream& write(OutputStream& os) const override = 0;
    virtual OutputStream& repr(OutputStream& os) const override { return this->write(os); }
};

template<class F, class D>
class FunctionMixin<F,Void,D,IntervalDomain>
    : public virtual FunctionInterface<Void,D,IntervalDomain>
{
    typedef IntervalDomain C;
    typedef typename FunctionInterface<Void,D,C>::DomainType DomainType;
    typedef typename FunctionInterface<Void,D,C>::CodomainType CodomainType;
  protected:
    FunctionMixin() { }
    template<class X> X _base_evaluate(const ElementType<D,X>& x) const;
  public:
    virtual DomainType const domain() const override { return make_domain<D>(this->argument_size()); }
    virtual CodomainType const codomain() const override { return make_domain<C>(this->result_size()); }
    virtual SizeType argument_size() const override { return this->domain().dimension(); }
    virtual SizeType result_size() const override { return 1u; }

    virtual OutputStream& write(OutputStream& os) const override = 0;
    virtual OutputStream& repr(OutputStream& os) const override { return this->write(os); }
};


template<class F, class D, class C>
class FunctionMixin<F,ApproximateTag,D,C>
    : public virtual FunctionInterface<ApproximateTag,D,C>
    , public FunctionMixin<F,Void,D,C>
{
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    virtual FunctionInterface<ApproximateTag,D,C>* _clone() const override;
    virtual Result<Float64Approximation> _evaluate(const Argument<Float64Approximation>& x) const override;
    virtual Result<FloatMPApproximation> _evaluate(const Argument<FloatMPApproximation>& x) const override;
    virtual Result<Differential<Float64Approximation>> _evaluate(const Argument<Differential<Float64Approximation>>& x) const override;
    virtual Result<Differential<FloatMPApproximation>> _evaluate(const Argument<Differential<FloatMPApproximation>>& x) const override;
    virtual Result<ApproximateTaylorModel> _evaluate(const Argument<ApproximateTaylorModel>& x) const override;
    virtual Result<ApproximateFormula> _evaluate(const Argument<ApproximateFormula>& x) const override;
    virtual Result<ApproximateAlgebra> _evaluate(const Argument<ApproximateAlgebra>& x) const override;
};

// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F, class D, class C>
class FunctionMixin<F,ValidatedTag,D,C>
    : public virtual FunctionInterface<ValidatedTag,D,C>
    , public FunctionMixin<F,ApproximateTag,D,C>
{
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    using FunctionMixin<F,ApproximateTag,D,C>::_evaluate;
    virtual FunctionInterface<ValidatedTag,D,C>* _clone() const override;
    virtual Result<Float64Bounds> _evaluate(const Argument<Float64Bounds>& x) const override;
    virtual Result<FloatMPBounds> _evaluate(const Argument<FloatMPBounds>& x) const override;
    virtual Result<Differential<Float64Bounds>> _evaluate(const Argument<Differential<Float64Bounds>>& x) const override;
    virtual Result<Differential<FloatMPBounds>> _evaluate(const Argument<Differential<FloatMPBounds>>& x) const override;
    virtual Result<ValidatedFormula> _evaluate(const Argument<ValidatedFormula>& x) const override;
    virtual Result<ValidatedTaylorModel> _evaluate(const Argument<ValidatedTaylorModel>& x) const override;
    virtual Result<ValidatedAlgebra> _evaluate(const Argument<ValidatedAlgebra>& x) const override;

    virtual Result<ValidatedScalarFunction> _evaluate(const Argument<ValidatedScalarFunction>& x) const override;

};

// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F, class D, class C>
class FunctionMixin<F,EffectiveTag,D,C>
    : public virtual FunctionInterface<EffectiveTag,D,C>
    , public FunctionMixin<F,ValidatedTag,D,C>
{
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    using FunctionMixin<F,ValidatedTag,D,C>::_evaluate;
    virtual FunctionInterface<EffectiveTag,D,C>* _clone() const override;
    virtual Result<EffectiveNumericType> _evaluate(const Argument<EffectiveNumericType>& x) const override;
    virtual Result<EffectiveFormula> _evaluate(const Argument<EffectiveFormula>& x) const override;
    virtual Result<EffectiveAlgebra> _evaluate(const Argument<EffectiveAlgebra>& x) const override;
};

template<class F, class P, class D> class ScalarFunctionMixin
    : public FunctionMixin<F,P,D,IntervalDomain> { };

template<class F, class P, class D> class VectorFunctionMixin
    : public FunctionMixin<F,P,D,BoxDomain>
    , public virtual VectorOfFunctionInterface<P,D>
{
    virtual ScalarFunctionInterface<P,D>* _get(SizeType i) const override = 0;
};


} // namespace Ariadne

#endif // ARIADNE_FUNCTION_TEMPLATE_H
