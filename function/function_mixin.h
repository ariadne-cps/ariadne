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

typedef ApproximateNumber ApproximateNumber;
typedef ValidatedNumber ValidatedNumber;
typedef EffectiveNumber EffectiveNumber;
typedef Differential<ApproximateNumber> ApproximateDifferential;
typedef Differential<ValidatedNumber> ValidatedDifferential;
typedef UnivariateDifferential<ApproximateNumber> ApproximateUnivariateDifferential;
typedef UnivariateDifferential<ValidatedNumber> ValidatedUnivariateDifferential;
typedef TaylorModel<Approximate,Float64> ApproximateTaylorModel;
typedef TaylorModel<Validated,Float64> ValidatedTaylorModel;
typedef Formula<ApproximateNumber> ApproximateFormula;
typedef Formula<ValidatedNumber> ValidatedFormula;
typedef Formula<EffectiveNumber> EffectiveFormula;
typedef Algebra<ApproximateNumber> ApproximateAlgebra;
typedef Algebra<ValidatedNumber> ValidatedAlgebra;
typedef Algebra<EffectiveNumber> EffectiveAlgebra;

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
    virtual DomainType const domain() const { return make_domain<D>(this->argument_size()); }
    virtual CodomainType const codomain() const override { return make_domain<C>(this->result_size()); }
    virtual SizeType argument_size() const override { return this->domain().dimension(); }
    virtual SizeType result_size() const override { return this->codomain().dimension(); }

    virtual OutputStream& write(OutputStream& os) const = 0;
    virtual OutputStream& repr(OutputStream& os) const { return this->write(os); }
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
    virtual DomainType const domain() const { return make_domain<D>(this->argument_size()); }
    virtual CodomainType const codomain() const override { return make_domain<C>(this->result_size()); }
    virtual SizeType argument_size() const override { return this->domain().dimension(); }
    virtual SizeType result_size() const override { return 1u; }

    virtual OutputStream& write(OutputStream& os) const = 0;
    virtual OutputStream& repr(OutputStream& os) const { return this->write(os); }
};


template<class F, class D, class C>
class FunctionMixin<F,ApproximateTag,D,C>
    : public virtual FunctionInterface<ApproximateTag,D,C>
    , public FunctionMixin<F,Void,D,C>
{
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    virtual FunctionInterface<ApproximateTag,D,C>* _clone() const;
    virtual Result<ApproximateNumber> _evaluate(const Argument<ApproximateNumber>& x) const;
    virtual Result<ApproximateDifferential> _evaluate(const Argument<ApproximateDifferential>& x) const;
    virtual Result<ApproximateFormula> _evaluate(const Argument<ApproximateFormula>& x) const;
    virtual Result<ApproximateTaylorModel> _evaluate(const Argument<ApproximateTaylorModel>& x) const;
    virtual Result<ApproximateAlgebra> _evaluate(const Argument<ApproximateAlgebra>& x) const;
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
    virtual FunctionInterface<ValidatedTag,D,C>* _clone() const;
    virtual Result<ValidatedNumber> _evaluate(const Argument<ValidatedNumber>& x) const;
    virtual Result<ValidatedDifferential> _evaluate(const Argument<ValidatedDifferential>& x) const;
    virtual Result<ValidatedFormula> _evaluate(const Argument<ValidatedFormula>& x) const;
    virtual Result<ValidatedTaylorModel> _evaluate(const Argument<ValidatedTaylorModel>& x) const;
    virtual Result<ValidatedAlgebra> _evaluate(const Argument<ValidatedAlgebra>& x) const;

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
    virtual FunctionInterface<EffectiveTag,D,C>* _clone() const;
    virtual Result<EffectiveNumber> _evaluate(const Argument<EffectiveNumber>& x) const;
    virtual Result<EffectiveFormula> _evaluate(const Argument<EffectiveFormula>& x) const;
    virtual Result<EffectiveAlgebra> _evaluate(const Argument<EffectiveAlgebra>& x) const;
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
