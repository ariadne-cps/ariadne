/***************************************************************************
 *            function/function_mixin.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#ifndef ARIADNE_FUNCTION_MIXIN_HPP
#define ARIADNE_FUNCTION_MIXIN_HPP

#include "../function/function_interface.hpp"

// Adaptors for classes to conform to the Function interface.

namespace Ariadne {

typedef Differential<FloatDPApproximation> ApproximateDifferentialDP;
typedef Differential<FloatDPBounds> ValidatedDifferentialDP;
typedef UnivariateDifferential<FloatDPApproximation> ApproximateUnivariateDifferentialDP;
typedef UnivariateDifferential<FloatDPBounds> ValidatedUnivariateDifferentialDP;
typedef TaylorModel<ApproximateTag,FloatDP> ApproximateTaylorModelDP;
typedef TaylorModel<ValidatedTag,FloatDP> ValidatedTaylorModelDP;
typedef TaylorModel<ApproximateTag,FloatMP> ApproximateTaylorModelMP;
typedef TaylorModel<ValidatedTag,FloatMP> ValidatedTaylorModelMP;

typedef TaylorModel<ValidatedTag,FloatDPUpperInterval> ValidatedIntervalTaylorModelDP;

typedef Formula<ApproximateNumber> ApproximateFormula;
typedef Formula<ValidatedNumber> ValidatedFormula;
typedef Formula<EffectiveNumber> EffectiveFormula;
typedef Algebra<ApproximateNumber> ApproximateAlgebra;
typedef Algebra<ValidatedNumber> ValidatedAlgebra;
typedef Algebra<EffectiveNumber> EffectiveAlgebra;
typedef ElementaryAlgebra<ApproximateNumber> ApproximateElementaryAlgebra;
typedef ElementaryAlgebra<ValidatedNumber> ValidatedElementaryAlgebra;
typedef ElementaryAlgebra<EffectiveNumber> EffectiveElementaryAlgebra;

template<class F, class P, class D, class C> class FunctionMixin { };
template<class F, class P, class D> class ScalarFunctionMixin;
template<class F, class P, class D> class VectorFunctionMixin;
template<class F, class P> using ScalarUnivariateFunctionMixin = ScalarFunctionMixin<F,P,IntervalDomainType>;
template<class F, class P> using VectorUnivariateFunctionMixin = VectorFunctionMixin<F,P,IntervalDomainType>;
template<class F, class P> using ScalarMultivariateFunctionMixin = ScalarFunctionMixin<F,P,BoxDomainType>;
template<class F, class P> using VectorMultivariateFunctionMixin = VectorFunctionMixin<F,P,BoxDomainType>;

template<class T> T* heap_copy(const T& t) { return new T(t); }
template<class T> T* heap_move(T&& t) { return new T(std::move(t)); }

template<class P, class D> ScalarFunctionInterface<P,D>* heap_copy(ScalarFunction<P,D> const& f) { return f.raw_pointer()->_clone(); }


template<class D> D make_domain(SizeType d);
template<> inline IntervalDomainType make_domain(SizeType d) { assert(d==1u); return IntervalDomainType(-inf,+inf); }
template<> inline BoxDomainType make_domain(SizeType d) { return BoxDomainType(d,IntervalDomainType(-inf,+inf)); }

template<class F, class D, class C>
class FunctionMixin<F,Void,D,C>
    : public virtual FunctionInterface<Void,D,C>
{
  public:
    typedef typename FunctionInterface<Void,D,C>::DomainType DomainType;
    typedef typename FunctionInterface<Void,D,C>::CodomainType CodomainType;
    typedef typename FunctionInterface<Void,D,C>::ArgumentSizeType ArgumentSizeType;
    typedef typename FunctionInterface<Void,D,C>::ResultSizeType ResultSizeType;
  protected:
    FunctionMixin() { }
    template<class X> ElementType<C,X> _base_evaluate(const ElementType<D,X>& x) const;
  public:
//    virtual DomainType const domain() const override = 0;
//    virtual CodomainType const codomain() const override = 0;
    virtual DomainType const domain() const override { return make_domain<D>(this->argument_size()); }
    virtual CodomainType const codomain() const override { return make_domain<C>(this->result_size()); }
    virtual ArgumentSizeType argument_size() const override { return this->domain().dimension(); }
    virtual ResultSizeType result_size() const override { return this->codomain().dimension(); }

    virtual OutputStream& _write(OutputStream& os) const override = 0;
    virtual OutputStream& repr(OutputStream& os) const override { return this->_write(os); }
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
    virtual Result<FloatDPApproximation> _evaluate(const Argument<FloatDPApproximation>& x) const override;
    virtual Result<FloatMPApproximation> _evaluate(const Argument<FloatMPApproximation>& x) const override;
    virtual Result<Differential<FloatDPApproximation>> _evaluate(const Argument<Differential<FloatDPApproximation>>& x) const override;
    virtual Result<Differential<FloatMPApproximation>> _evaluate(const Argument<Differential<FloatMPApproximation>>& x) const override;
    virtual Result<TaylorModel<ApproximateTag,FloatDP>> _evaluate(const Argument<TaylorModel<ApproximateTag,FloatDP>>& x) const override;
    virtual Result<TaylorModel<ApproximateTag,FloatMP>> _evaluate(const Argument<TaylorModel<ApproximateTag,FloatMP>>& x) const override;
    virtual Result<Formula<ApproximateNumber>> _evaluate(const Argument<Formula<ApproximateNumber>>& x) const override;
    virtual Result<ElementaryAlgebra<ApproximateNumber>> _evaluate(const Argument<ElementaryAlgebra<ApproximateNumber>>& x) const override;
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
    virtual Result<FloatDPBounds> _evaluate(const Argument<FloatDPBounds>& x) const override;
    virtual Result<FloatMPBounds> _evaluate(const Argument<FloatMPBounds>& x) const override;
    virtual Result<Differential<FloatDPBounds>> _evaluate(const Argument<Differential<FloatDPBounds>>& x) const override;
    virtual Result<Differential<FloatMPBounds>> _evaluate(const Argument<Differential<FloatMPBounds>>& x) const override;
    virtual Result<TaylorModel<ValidatedTag,FloatDP>> _evaluate(const Argument<TaylorModel<ValidatedTag,FloatDP>>& x) const override;
    virtual Result<TaylorModel<ValidatedTag,FloatMP>> _evaluate(const Argument<TaylorModel<ValidatedTag,FloatMP>>& x) const override;
    virtual Result<TaylorModel<ValidatedTag,FloatDPUpperInterval>> _evaluate(const Argument<TaylorModel<ValidatedTag,FloatDPUpperInterval>>& x) const override;
    virtual Result<TaylorModel<ValidatedTag,FloatMPUpperInterval>> _evaluate(const Argument<TaylorModel<ValidatedTag,FloatMPUpperInterval>>& x) const override;
    virtual Result<Formula<ValidatedNumber>> _evaluate(const Argument<Formula<ValidatedNumber>>& x) const override;
    virtual Result<ElementaryAlgebra<ValidatedNumber>> _evaluate(const Argument<ElementaryAlgebra<ValidatedNumber>>& x) const override;

    virtual Result<ValidatedScalarMultivariateFunction> _evaluate(const Argument<ValidatedScalarMultivariateFunction>& x) const override;

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
    virtual Result<Real> _evaluate(const Argument<Real>& x) const override;
    virtual Result<ElementaryAlgebra<Real>> _evaluate(const Argument<ElementaryAlgebra<Real>>& x) const override;
    virtual Result<Formula<Real>> _evaluate(const Argument<Formula<Real>>& x) const override;
    virtual Result<ElementaryAlgebra<EffectiveNumber>> _evaluate(const Argument<ElementaryAlgebra<EffectiveNumber>>& x) const override;
    virtual Result<Formula<EffectiveNumber>> _evaluate(const Argument<Formula<EffectiveNumber>>& x) const override;
};

template<class F, class P, class D> class ScalarFunctionMixin
    : public FunctionMixin<F,P,D,IntervalDomainType> { };

template<class F, class P, class D> class VectorFunctionMixin
    : public FunctionMixin<F,P,D,BoxDomainType>
    , public virtual VectorOfFunctionInterface<P,D>
{
    virtual ScalarFunctionInterface<P,D>* _get(SizeType i) const override {
        auto fi=static_cast<F const&>(*this)[i]; return heap_copy(fi); }
};


template<class F,class D, class C> FunctionInterface<ApproximateTag,D,C>* FunctionMixin<F,ApproximateTag,D,C>::_clone() const {
    return new F(static_cast<const F&>(*this)); }
template<class F,class D, class C> FunctionInterface<ValidatedTag,D,C>* FunctionMixin<F,ValidatedTag,D,C>::_clone() const {
    return new F(static_cast<const F&>(*this)); }
template<class F,class D, class C> FunctionInterface<EffectiveTag,D,C>* FunctionMixin<F,EffectiveTag,D,C>::_clone() const {
    return new F(static_cast<const F&>(*this)); }



} // namespace Ariadne

#endif // ARIADNE_FUNCTION_TEMPLATE_HPP
