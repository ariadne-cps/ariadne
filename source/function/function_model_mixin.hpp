/***************************************************************************
 *            function_model_mixin.hpp
 *
 *  Copyright 2011-17  Pieter Collins
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

/*! \file function_model_mixin.hpp
 *  \brief Mixin for concrete functions on bounded domains.
 */

#ifndef ARIADNE_FUNCTION_MODEL_MIXIN_HPP
#define ARIADNE_FUNCTION_MODEL_MIXIN_HPP

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "../function/function_model_interface.hpp"
#include "../function/function_model.hpp"

#include "../numeric/operators.hpp"
#include "../numeric/numeric.hpp"
#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/operations.hpp"
#include "../function/domain.hpp"

#include "../function/function_interface.hpp"
#include "../function/function_mixin.hpp"
#include "../function/function.hpp"

namespace Ariadne {

template<class FM, class P, class D, class C, class PR=DoublePrecision, class PRE=PR> class FunctionModelMixin;
template<class FM, class P, class D, class PR=DoublePrecision, class PRE=PR> using ScalarMultivariateFunctionModelMixin = FunctionModelMixin<FM,P,D,IntervalDomainType,PR,PRE>;
template<class FM, class P, class D, class PR=DoublePrecision, class PRE=PR> using VectorMultivariateFunctionModelMixin = FunctionModelMixin<FM,P,D,BoxDomainType,PR,PRE>;

template<class FM, class P, class D, class PR, class PRE> class FunctionModelMixin<FM,P,D,IntervalDomainType,PR,PRE>
    : public virtual ScalarFunctionModelInterface<P,D,PR,PRE>
    , public ScalarFunctionMixin<FM,P,D>
{
    typedef FloatError<PR> NormType;
  public:
    FM apply(OperatorCode op) const;
  public:
    ScalarFunctionModelInterface<P,D,PR,PRE>* _clone() const override {
        return new FM(static_cast<const FM&>(*this)); }
    NormType const _norm() const override {
        return norm(static_cast<const FM&>(*this)); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _antiderivative(SizeType j) const override {
        return new FM(antiderivative(static_cast<const FM&>(*this),j)); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _antiderivative(SizeType j, CanonicalNumericType<P,PR,PRE> c) const override {
        return new FM(antiderivative(static_cast<const FM&>(*this),j,c)); }
     ScalarFunctionModelInterface<P,D,PR,PRE>* _restriction(const BoxDomainType& d) const override {
        return new FM(restriction(static_cast<const FM&>(*this),d)); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _apply(OperatorCode op) const override {
        return new FM(this->apply(op)); }
    CanonicalNumericType<P,PR,PRE> _unchecked_evaluate(const Vector<CanonicalNumericType<P,PR,PRE>>& x) const override {
        return unchecked_evaluate(static_cast<const FM&>(*this),x); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _partial_evaluate(SizeType j, const CanonicalNumericType<P,PR,PRE>& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FM&>(*this),j,c)); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _embed(const BoxDomainType& d1, const BoxDomainType& d2) const override {
        return new FM(embed(d1,static_cast<const FM&>(*this),d2)); }
    Boolean _refines(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return refines(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f)); }
    Boolean _inconsistent(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return inconsistent(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f)); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _refinement(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) const override {
        ARIADNE_ASSERT(dynamic_cast<const FM*>(&f)); return new FM(refinement(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    Void _iadd(const CanonicalNumericType<P,PR,PRE>& c) override {
        static_cast<FM&>(*this)+=c; }
    Void _imul(const CanonicalNumericType<P,PR,PRE>& c) override {
        static_cast<FM&>(*this)*=c; }
    Void _isma(const CanonicalNumericType<P,PR,PRE>& c, const ScalarFunctionModelInterface<P,D,PR,PRE>& f) override {
        static_cast<FM&>(*this)+=c*dynamic_cast<const FM&>(f); }
    Void _ifma(const ScalarFunctionModelInterface<P,D,PR,PRE>& f1, const ScalarFunctionModelInterface<P,D,PR,PRE>& f2) override {
        static_cast<FM&>(*this)+=dynamic_cast<const FM&>(f1)*dynamic_cast<const FM&>(f2); }
};

template<class FM, class P, class D, class PR, class PRE> FM ScalarMultivariateFunctionModelMixin<FM,P,D,PR,PRE>::apply(OperatorCode op) const {
    const FM& f=static_cast<const FM&>(*this);
    switch(op) {
        case OperatorCode::NEG: return neg(f);
        case OperatorCode::REC: return rec(f);
        case OperatorCode::EXP: return exp(f);
        default: ARIADNE_FAIL_MSG("ScalarFunctionModel<P,D,PR,PRE>::apply(OperatorCode op): Operator op="<<op<<" not implemented\n");
    }
}


template<class FM, class P, class D, class PR, class PRE> class FunctionModelMixin<FM,P,D,BoxDomainType,PR,PRE>
    : public virtual VectorFunctionModelInterface<P,D,PR,PRE>
    , public VectorFunctionMixin<FM,P,D>
{
    typedef typename Element<FM>::Type ScalarMultivariateFunctionType;
    typedef FloatError<PR> NormType;
  public:
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _clone() const override { return new FM(static_cast<const FM&>(*this)); }
    virtual Void _set(SizeType i, const ScalarFunctionModelInterface<P,D,PR,PRE>& sf) override {
        if(!dynamic_cast<const typename FM::ScalarMultivariateFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorMultivariateFunctionModel "<<*this<<" to "<<sf<<"\n"); }
        static_cast<FM&>(*this).FM::set(i,dynamic_cast<const ScalarMultivariateFunctionType&>(sf)); }
    virtual VectorFunctionModelInterface<P,D,PR,PRE>* _derivative(SizeType j) const override {
        ARIADNE_NOT_IMPLEMENTED; }
    NormType const _norm() const override {
         return norm(static_cast<const FM&>(*this)); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _embed(const BoxDomainType& d1, const BoxDomainType& d2) const override {
        return heap_copy(embed(d1,static_cast<const FM&>(*this),d2)); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _restriction(const BoxDomainType& d) const override {
        return new FM(restriction(static_cast<const FM&>(*this),d)); }
    Void _adjoin(const ScalarFunctionModelInterface<P,D,PR,PRE>& f) override {
        static_cast<FM&>(*this).FM::adjoin(dynamic_cast<const ScalarMultivariateFunctionType&>(f)); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _join(const VectorFunctionModelInterface<P,D,PR,PRE>& f) const override {
        return heap_copy(join(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _combine(const VectorFunctionModelInterface<P,D,PR,PRE>& f) const override {
        return heap_copy(combine(static_cast<const FM&>(*this),dynamic_cast<const FM&>(f))); }
    Vector<CanonicalNumericType<P,PR,PRE>> _unchecked_evaluate(const Vector<CanonicalNumericType<P,PR,PRE>>& x) const override {
        return unchecked_evaluate(static_cast<const FM&>(*this),x); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _compose(const ScalarMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(compose(f,static_cast<const FM&>(*this))); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _compose(const VectorMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(compose(f,static_cast<const FM&>(*this))); }
    ScalarFunctionModelInterface<P,D,PR,PRE>* _unchecked_compose(const ScalarMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<const ScalarMultivariateFunctionType&>(f),static_cast<const FM&>(*this))); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _unchecked_compose(const VectorMultivariateFunctionInterface<P>& f) const override {
        return heap_copy(unchecked_compose(dynamic_cast<const FM&>(f),static_cast<const FM&>(*this))); }
    VectorFunctionModelInterface<P,D,PR,PRE>* _partial_evaluate(SizeType j, const CanonicalNumericType<P,PR,PRE>& c) const override {
        return heap_copy(partial_evaluate(static_cast<const FM&>(*this),j,c)); }
};


template<class FCTRY, class P, class PR, class PRE> class FunctionModelFactoryMixin
    : public FunctionModelFactoryInterface<P,PR,PRE>
{
    typedef BoxDomainType VD;
    typedef IntervalDomainType SD;
    friend class FunctionModelFactory<P,PR,PRE>;
  public:
    typedef VD VectorDomainType;
    typedef SD ScalarDomainType;

    virtual FunctionModelFactoryInterface<P,PR,PRE>* clone() const override { return new FCTRY(this->upcast()); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << this->upcast(); }
/*
    CanonicalNumericType<P,PR,PRE> create(const Number<P>& number) const;
    ScalarFunctionModel<P,D,PR,PRE> create(const BoxDomainType& domain, const ScalarMultivariateFunctionInterface<P>& function) const;
    VectorFunctionModel<P,D,PR,PRE> create(const BoxDomainType& domain, const VectorMultivariateFunctionInterface<P>& function) const;
    ScalarFunctionModel<P,D,PR,PRE> create_zero(const BoxDomainType& domain) const;
    VectorFunctionModel<P,D,PR,PRE> create_zeros(SizeType result_size, const BoxDomainType& domain) const;
    ScalarFunctionModel<P,D,PR,PRE> create_constant(const BoxDomainType& domain, const Number<P>& value) const;
    ScalarFunctionModel<P,D,PR,PRE> create_constant(const BoxDomainType& domain, const CanonicalNumericType<P,PR,PRE>& value) const;
    VectorFunctionModel<P,D,PR,PRE> create_constants(const BoxDomainType& domain, const Vector<Number<P>>& values) const;
    VectorFunctionModel<P,D,PR,PRE> create_constants(const BoxDomainType& domain, const Vector<CanonicalNumericType<P,PR,PRE>>& values) const;
    ScalarFunctionModel<P,D,PR,PRE> create_coordinate(const BoxDomainType& domain, SizeType index) const;
    ScalarFunctionModel<P,D,PR,PRE> create_identity(const IntervalDomainType& domain) const;
    VectorFunctionModel<P,D,PR,PRE> create_identity(const BoxDomainType& domain) const;
    CanonicalNumericType<P,PR,PRE> create_number(const Number<P>& number) const;
*/
  private:
    template<class T> static inline T* heap_move(T&& t) { return new T(std::forward<T>(t)); }
    inline FCTRY const& upcast() const { return static_cast<FCTRY const&>(*this); }
    virtual CanonicalNumericType<P,PR,PRE> _create(const Number<P>& number) const override {
        return this->upcast().create(number); }
    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create(const VectorDomainType& domain, const ScalarFunctionInterface<P,VD>& function) const override {
        return heap_move(this->upcast().create(domain,function)); };
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create(const VectorDomainType& domain, const VectorFunctionInterface<P,VD>& function) const override {
        return heap_move(this->upcast().create(domain,function)); };

    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create_zero(const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_zero(domain)); };
    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create_constant(const VectorDomainType& domain, const Number<P>& value) const override {
        return heap_move(this->upcast().create_constant(domain,value)); };
    virtual ScalarFunctionModelInterface<P,VD,PR,PRE>* _create_coordinate(const VectorDomainType& domain, SizeType j) const override {
        return heap_move(this->upcast().create_coordinate(domain,j)); };
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create_zeros(SizeType n, const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_zeros(n,domain)); };
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create_constants(const VectorDomainType& domain, const Vector<Number<P>>& values) const override {
        return heap_move(this->upcast().create_constants(domain,values)); };
    virtual VectorFunctionModelInterface<P,VD,PR,PRE>* _create_identity(const VectorDomainType& domain) const override {
        return heap_move(this->upcast().create_identity(domain)); };
};


} // namespace Ariadne

#endif
