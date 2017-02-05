/***************************************************************************
 *            function_model.h
 *
 *  Copyright 2011-17  Pieter Collins
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

/*! \file function_model.h
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_MODEL_H
#define ARIADNE_FUNCTION_MODEL_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function/function.decl.h"
#include "function/function_model_interface.h"

#include "numeric/operators.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/operations.h"
#include "function/domain.h"

#include "function/function_interface.h"
#include "function/function_mixin.h"
#include "function/function.h"

namespace Ariadne {

template<class P, class PR, class PRE> class FunctionModelFactory {
    SharedPointer<const FunctionModelFactoryInterface<P,PR,PRE>> _ptr;
  public:
    typedef P Paradigm;
    typedef PR PrecisionType;
    typedef PRE ErrorPrecisionType;
    typedef BoxDomainType DomainType;

    explicit FunctionModelFactory(const FunctionModelFactoryInterface<P,PR,PRE>* p) : _ptr(p) { }
    explicit FunctionModelFactory(SharedPointer<const FunctionModelFactoryInterface<P,PR,PRE>> p) : _ptr(p) { }

    CanonicalNumericType<P,PR,PRE> create(Number<P> const& c) const { return this->_ptr->_create(c); }
    ScalarFunctionModel<P,PR,PRE> create(DomainType const& dom, ScalarFunction<P> const& f) { return this->_ptr->_create(dom,f); }
    VectorFunctionModel<P,PR,PRE> create(DomainType const& dom, VectorFunction<P> const& f) { return this->_ptr->_create(dom,f); }
    ScalarFunctionModel<P,PR,PRE> create_zero(DomainType const& dom) const { return this->_ptr->_create_zero(dom); }
    ScalarFunctionModel<P,PR,PRE> create_constant(DomainType const& dom, Number<P> const& c) const { return this->_ptr->_create_constant(dom,c); }
    ScalarFunctionModel<P,PR,PRE> create_coordinate(DomainType const& dom, SizeType j) const { return this->_ptr->_create_cordinate(dom,j); }
    VectorFunctionModel<P,PR,PRE> create_zeros(SizeType n, DomainType const& dom) const { return this->_ptr->_create_zeros(n,dom); }
    VectorFunctionModel<P,PR,PRE> create_identity(DomainType const& dom) const { return this->_ptr->_create_identity(dom); }
    ScalarFunctionModel<P,PR,PRE> create_identity(IntervalDomainType const& dom) const { return this->_ptr->_create_coordinate(BoxDomainType(1u,dom),0u); }

    friend OutputStream& operator<<(OutputStream& os, FunctionModelFactory<P,PR,PRE> const& factory) { return factory._ptr->_write(os); }
};

template<class FCTRY> class FunctionModelCreator {
    typedef typename FCTRY::Paradigm P;
    typedef typename FCTRY::PrecisionType PR;
    typedef typename FCTRY::ErrorPrecisionType PRE;
  public:
    typedef FCTRY FactoryType;
    typedef BoxDomainType DomainType;
    typedef P Paradigm;

    explicit FunctionModelCreator(DomainType domain, FactoryType factory) : _factory(factory), _domain(domain) { }

    decltype(auto) create(Number<P> const& c) const { return this->_factory.create(c); }
    decltype(auto) create(ScalarFunction<P> const& f) { return this->_factory.create(this->_domain,f); }
    decltype(auto) create(VectorFunction<P> const& f) { return this->_factory.create(this->_domain,f); }
    decltype(auto) create_zero() { return this->_factory.create_zero(this->_domain); }
    decltype(auto) create_zeros(SizeType n) { return this->_factory.create_zeros(n,this->_domain); }
    decltype(auto) create_identity() { return this->_factory.create_identity(this->_domain); }

//    ScalarFunctionModel<P,PR,PRE> const& create(ScalarFunctionModel<P,PR,PRE> const& f) const { return f; }
    CanonicalNumericType<P,PR,PRE> const& create(CanonicalNumericType<P,PR,PRE> const& c) const { return c; }
  protected:
    FactoryType _factory;
    DomainType _domain;
};

//! \ingroup FunctionModelSubModule
//! \brief Generic scalar functions on singleton domains.
template<class P, class PR, class PRE> class ScalarFunctionModel
    : public DispatchAlgebraOperators<ScalarFunctionModel<P,PR,PRE>, CanonicalNumericType<P,PR,PRE>>
    , public ProvideConcreteGenericArithmeticOperators<ScalarFunctionModel<P,PR,PRE>>
    , public ProvideConcreteGenericArithmeticOperators<ScalarFunctionModel<P,PR,PRE>,Number<P>>
{
  public:
    typedef ScalarFunction<P> GenericType;
    typedef BoxDomainType DomainType;
    typedef IntervalDomainType CodomainType;
    typedef CanonicalCoefficientType<P,PR> CoefficientType;
    typedef CanonicalErrorType<P,PRE> ErrorType;
    typedef CanonicalNumericType<P,PR,PRE> NumericType;
    typedef FloatError<PR> NormType;
    typedef Interval<FloatUpperBound<PR>> RangeType;
  public:
    clone_on_copy_ptr< ScalarFunctionModelInterface<P,PR,PRE> > _ptr;
  public:
    ScalarFunctionModel() : _ptr() { }
    ScalarFunctionModel(ScalarFunctionModelInterface<P,PR,PRE>* p) : _ptr(p) { }
    ScalarFunctionModel(const SharedPointer<const ScalarFunctionModelInterface<P,PR,PRE>> p) : _ptr(p->_clone()) { }
    ScalarFunctionModel(const ScalarFunctionModel<P,PR,PRE>& f) : _ptr(f._ptr) { }
    ScalarFunctionModel(const ScalarFunctionModelInterface<P,PR,PRE>& f) : _ptr(f._clone()) { }
    ScalarFunctionModel(const ScalarFunction<P>& f) : _ptr(dynamic_cast<ScalarFunctionModelInterface<P,PR,PRE>*>(f.raw_pointer()->_clone())) { }
    operator ScalarFunction<P>() const { return ScalarFunction<P>(this->_ptr->_clone()); }
    operator ScalarFunctionModelInterface<P,PR,PRE>& () { return *_ptr; }
    operator const ScalarFunctionModelInterface<P,PR,PRE>& () const { return *_ptr; }
    const ScalarFunctionModelInterface<P,PR,PRE>* raw_pointer() const { return _ptr.operator->(); }
    ScalarFunctionModelInterface<P,PR,PRE>& reference() { return *_ptr; }
    const ScalarFunctionModelInterface<P,PR,PRE>& reference() const { return *_ptr; }
    ScalarFunctionModel<P,PR,PRE>& operator=(const Number<P>& c);
    ScalarFunctionModel<P,PR,PRE>& operator=(const CanonicalNumericType<P,PR,PRE>& c);
    ScalarFunctionModel<P,PR,PRE>& operator=(const ScalarFunction<P>& f);
    ScalarFunctionModel<P,PR,PRE>& operator=(const ScalarFunctionModelInterface<P,PR,PRE>& f);
//    ScalarFunctionModel<P,PR,PRE>& operator=(const ValidatedScalarTaylorFunctionModel64& f);
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    template<class X> X operator() (const Vector<X>& x) const {
        return this->_ptr->_evaluate(x); }
    template<class X> X evaluate(const Vector<X>& x) const {
        return this->_ptr->_evaluate(x); }
    inline DomainType const domain() const { return this->_ptr->domain(); }
    inline CodomainType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const { return this->_ptr->range(); }

    inline CoefficientType value() const { return this->_ptr->value(); }
    inline CoefficientType gradient_value(SizeType j) const { return this->_ptr->gradient_value(j); }
    inline ErrorType error() const { return this->_ptr->error(); }

    inline Void set_error(const ErrorType& e) { return this->_ptr->set_error(e); }
    inline Void set_error(Nat e) { return this->_ptr->set_error(ErrorType(e,this->error().precision())); }
    inline Void clobber() { return this->_ptr->clobber(); }

    inline ScalarFunctionModel<P,PR,PRE> apply(Operator op) const { return this->_ptr->_apply(op); }
    inline Void restrict(const DomainType& d) { *this=restriction(*this,d); }
  public:
    friend FunctionModelCreator<FunctionModelFactory<P,PR,PRE>> factory(ScalarFunctionModel<P,PR,PRE> const& f) {
        FunctionModelFactory<P,PR,PRE> factory(f._ptr->_factory()); return FunctionModelCreator<FunctionModelFactory<P,PR,PRE>>(f.domain(),factory); }
  public:
    friend inline ScalarFunctionModel<P,PR,PRE> operator/(const ScalarFunctionModel<P,PR,PRE>& f1, const ScalarFunctionModel<P,PR,PRE>& f2) {
        return mul(f1,rec(f2)); }
    friend inline ScalarFunctionModel<P,PR,PRE> operator/(const CanonicalNumericType<P,PR,PRE>& c1, const ScalarFunctionModel<P,PR,PRE>& f2) {
        return mul(c1,rec(f2)); }
    friend inline ScalarFunctionModel<P,PR,PRE> neg(ScalarFunctionModel<P,PR,PRE> f) {
        f._ptr->_imul(CanonicalNumericType<P,PR,PRE>(-1)); return std::move(f); }
    friend inline ScalarFunctionModel<P,PR,PRE> add(ScalarFunctionModel<P,PR,PRE> f1, const ScalarFunctionModel<P,PR,PRE>& f2) {
        f1._ptr->_isma(CanonicalNumericType<P,PR,PRE>(+1),f2); return std::move(f1); }
    friend inline ScalarFunctionModel<P,PR,PRE> sub(ScalarFunctionModel<P,PR,PRE> f1, const ScalarFunctionModel<P,PR,PRE>& f2) {
        f1._ptr->_isma(CanonicalNumericType<P,PR,PRE>(-1),f2); return std::move(f1); }
    friend inline ScalarFunctionModel<P,PR,PRE> mul(const ScalarFunctionModel<P,PR,PRE>& f1, const ScalarFunctionModel<P,PR,PRE>& f2) {
        ScalarFunctionModel<P,PR,PRE> r=factory(f1).create_zero(); r._ptr->_ifma(f1,f2); return r; }
    friend inline ScalarFunctionModel<P,PR,PRE> div(const ScalarFunctionModel<P,PR,PRE>& f1, const ScalarFunctionModel<P,PR,PRE>& f2) {
        return mul(f1,rec(f2)); }
    friend inline ScalarFunctionModel<P,PR,PRE> add(ScalarFunctionModel<P,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
        f1._ptr->_iadd(c2); return std::move(f1); }
    friend inline ScalarFunctionModel<P,PR,PRE> mul(ScalarFunctionModel<P,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
        f1._ptr->_imul(c2); return std::move(f1); }
    friend inline ScalarFunctionModel<P,PR,PRE> add(const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,PR,PRE> f2) {
        f2._ptr->_iadd(c1); return std::move(f2); }
    friend inline ScalarFunctionModel<P,PR,PRE> mul(const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,PR,PRE> f2) {
        f2._ptr->_imul(c1); return std::move(f2); }

    friend inline ScalarFunctionModel<P,PR,PRE> sub(ScalarFunctionModel<P,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
        return add(std::move(f1),neg(c2)); }
    friend inline ScalarFunctionModel<P,PR,PRE> div(ScalarFunctionModel<P,PR,PRE> f1, const CanonicalNumericType<P,PR,PRE>& c2) {
        return mul(std::move(f1),rec(c2)); }
    friend inline ScalarFunctionModel<P,PR,PRE> sub(const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,PR,PRE> f2) {
        return add(neg(std::move(f2)),c1); }
    friend inline ScalarFunctionModel<P,PR,PRE> div(const CanonicalNumericType<P,PR,PRE>& c1, ScalarFunctionModel<P,PR,PRE> f2) {
        return mul(rec(std::move(f2)),c1); }

    friend inline ScalarFunctionModel<P,PR,PRE> add(ScalarFunctionModel<P,PR,PRE> f1, const Number<P>& c2) {
        CanonicalNumericType<P,PR,PRE> s2=factory(f1).create(c2); return add(f1,s2); }
    friend inline ScalarFunctionModel<P,PR,PRE> sub(ScalarFunctionModel<P,PR,PRE> f1, const Number<P>& c2) {
        return add(f1,neg(c2)); }
    friend inline ScalarFunctionModel<P,PR,PRE> mul(ScalarFunctionModel<P,PR,PRE> f1, const Number<P>& c2) {
        CanonicalNumericType<P,PR,PRE> s2=factory(f1).create(c2); return mul(f1,s2); }
    friend inline ScalarFunctionModel<P,PR,PRE> div(ScalarFunctionModel<P,PR,PRE> f1, const Number<P>& c2) {
        return mul(f1,rec(c2)); }
    friend inline ScalarFunctionModel<P,PR,PRE> add(const Number<P>& c1, ScalarFunctionModel<P,PR,PRE> f2) {
        return add(f2,c1); }
    friend inline ScalarFunctionModel<P,PR,PRE> sub(const Number<P>& c1, ScalarFunctionModel<P,PR,PRE> f2) {
        return add(neg(f2),c1); }
    friend inline ScalarFunctionModel<P,PR,PRE> mul(const Number<P>& c1, ScalarFunctionModel<P,PR,PRE> f2) {
        return mul(f2,c1); }
    friend inline ScalarFunctionModel<P,PR,PRE> div(const Number<P>& c1, ScalarFunctionModel<P,PR,PRE> f2) {
        return mul(rec(f2),c1); }

    friend inline ScalarFunctionModel<P,PR,PRE> pow(const ScalarFunctionModel<P,PR,PRE>& f1, Int n2) {
        return generic_pow(f1,n2); }

    friend inline ScalarFunctionModel<P,PR,PRE> rec(const ScalarFunctionModel<P,PR,PRE>& f) {
        return f.apply(Rec()); }

  public:
    friend ScalarFunctionModel<P,PR,PRE> partial_evaluate(const ScalarFunctionModel<P,PR,PRE>& f, SizeType j, const CanonicalNumericType<P,PR,PRE>& c) {
        return f._ptr->_partial_evaluate(j,c); }

    friend NormType norm(const ScalarFunctionModel<P,PR,PRE>& f) { return f._ptr->_norm(); }
    friend ScalarFunctionModel<P,PR,PRE> derivative(const ScalarFunctionModel<P,PR,PRE>& f, SizeType j) { return f._ptr->_derivative(j); }
    friend ScalarFunctionModel<P,PR,PRE> antiderivative(const ScalarFunctionModel<P,PR,PRE>& f, SizeType j) { return f._ptr->_antiderivative(j); }
    friend ScalarFunctionModel<P,PR,PRE> antiderivative(const ScalarFunctionModel<P,PR,PRE>& f, SizeType j, CanonicalNumericType<P,PR,PRE> c) { return f._ptr->_antiderivative(j,c); }

    friend ScalarFunctionModel<P,PR,PRE> embed(const DomainType& d1, const ScalarFunctionModel<P,PR,PRE>& f, const DomainType& d2) {
        return f._ptr->_embed(d1,d2); }
    friend ScalarFunctionModel<P,PR,PRE> embed(const DomainType& d, const ScalarFunctionModel<P,PR,PRE>& f) {
        return embed(d,f,DomainType()); }
    friend ScalarFunctionModel<P,PR,PRE> embed(const ScalarFunctionModel<P,PR,PRE>& f, const DomainType& d) {
        return embed(DomainType(),f,d); }
    friend ScalarFunctionModel<P,PR,PRE> embed(const ScalarFunctionModel<P,PR,PRE>& f, const IntervalDomainType& d) {
        return embed(f,DomainType(1,d)); }
    friend ScalarFunctionModel<P,PR,PRE> restrict(const ScalarFunctionModel<P,PR,PRE>& f, const DomainType& d) {
        return f._ptr->_restriction(d); }
    friend ScalarFunctionModel<P,PR,PRE> restriction(const ScalarFunctionModel<P,PR,PRE>& f, const DomainType& d) {
        return f._ptr->_restriction(d); }

    friend VectorFunctionModel<P,PR,PRE> join(const ScalarFunctionModel<P,PR,PRE>& f1, const ScalarFunctionModel<P,PR,PRE>& f2) {
        return join(VectorFunctionModel<P,PR,PRE>(1,f1),f2); }
    friend VectorFunctionModel<P,PR,PRE> combine(const ScalarFunctionModel<P,PR,PRE>& f1, const ScalarFunctionModel<P,PR,PRE>& f2);
  public:
    friend OutputStream& operator<<(OutputStream& os, const ScalarFunctionModel<P,PR,PRE>& f) {
        return os <<  f.operator ScalarFunction<P>(); }
};

inline ScalarFunctionModel64<ValidatedTag> refinement(const ScalarFunctionModel64<ValidatedTag>& f1, const ScalarFunctionModel64<ValidatedTag>& f2) {
    return f1._ptr->_refinement(f2); }
inline Boolean inconsistent(const ScalarFunctionModel64<ValidatedTag>& f1, const ScalarFunctionModel64<ValidatedTag>& f2) {
    return f1._ptr->_inconsistent(f2); }
inline Boolean refines(const ScalarFunctionModel64<ValidatedTag>& f1, const ScalarFunctionModel64<ValidatedTag>& f2) {
    return f1._ptr->_refines(f2); }

// FIXME: Should not be needed since ScalarFunctionModel has a representation
template<class P> inline ScalarFunctionModel64<P> embed(const ScalarFunction<P>& f, const IntervalDomainType& d) {
    return embed(ScalarFunctionModel64<P>(f),d); }

template<class P, class PR, class PRE> inline ScalarFunctionModel<P,PR,PRE>& ScalarFunctionModel<P,PR,PRE>::operator=(const CanonicalNumericType<P,PR,PRE>& c) {
        (*this)*=CanonicalNumericType<P,PR,PRE>(0); (*this)+=c; return *this; }
template<class P, class PR, class PRE> inline ScalarFunctionModel<P,PR,PRE>& ScalarFunctionModel<P,PR,PRE>::operator=(const Number<P>& c) {
        return (*this)=factory(*this).create(c); }
template<class P, class PR, class PRE> inline ScalarFunctionModel<P,PR,PRE>& ScalarFunctionModel<P,PR,PRE>::operator=(const ScalarFunction<P>& f) {
        return (*this)=factory(*this).create(f); }
template<class P, class PR, class PRE> inline ScalarFunctionModel<P,PR,PRE>& ScalarFunctionModel<P,PR,PRE>::operator=(const ScalarFunctionModelInterface<P,PR,PRE>& f) {
        return (*this)=f._clone(); }



template<class V> struct Element;

template<class M> class ScaledFunctionPatch;
template<class M> class VectorScaledFunctionPatch;
template<class M> struct Element<VectorScaledFunctionPatch<M>> { typedef ScaledFunctionPatch<M> Type; };

typedef ScaledFunctionPatch<ValidatedTaylorModel64> ValidatedScalarTaylorFunctionModel64;

template<class P, class PR, class PRE> class VectorFunctionModelElement
    : public DispatchAlgebraOperators<ScalarFunctionModel<P,PR,PRE>, CanonicalNumericType<P,PR,PRE>>
    , public ProvideConcreteGenericArithmeticOperators<ScalarFunctionModel<P,PR,PRE>>
{
    VectorFunctionModel<P,PR,PRE>* _p; SizeType _i;
  public:
    typedef typename ScalarFunctionModel<P,PR,PRE>::GenericType GenericType;
    operator const ScalarFunctionModel<P,PR,PRE> () const;
    VectorFunctionModelElement(VectorFunctionModel<P,PR,PRE>* p, SizeType i) : _p(p), _i(i) { }
    VectorFunctionModelElement<P,PR,PRE>& operator=(const ScalarFunctionModel<P,PR,PRE>& sf) {
        _p->set(_i,sf); return *this; }
    VectorFunctionModelElement<P,PR,PRE>& operator=(const VectorFunctionModelElement<P,PR,PRE>& sf) {
        return this->operator=(static_cast<ScalarFunctionModel<P,PR,PRE>const>(sf)); }
    Void clobber() { ScalarFunctionModel<P,PR,PRE> sf=_p->get(_i); sf.clobber(); _p->set(_i,sf); }
    Void set_error(CanonicalErrorType<P,PRE> e) const { ScalarFunctionModel<P,PR,PRE> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
    Void set_error(Nat e) const { ScalarFunctionModel<P,PR,PRE> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
    friend ScalarFunctionModel<P,PR,PRE> antiderivative(VectorFunctionModelElement<P,PR,PRE> const& f, SizeType k) {
        return antiderivative(ScalarFunctionModel<P,PR,PRE>(f),k); }
};
template<class P, class PR, class PRE> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionModelElement<P,PR,PRE>& function) {
    return os << static_cast< const ScalarFunctionModel<P,PR,PRE> >(function);
}

//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on singleton domains.
template<class P, class PR, class PRE> class VectorFunctionModel
{
  public:
    clone_on_copy_ptr< VectorFunctionModelInterface<P,PR,PRE> > _ptr;
  public:
    typedef BoxDomainType DomainType;
    typedef BoxDomainType CodomainType;
    typedef CanonicalCoefficientType<P,PR> CoefficientType;
    typedef CanonicalErrorType<P,PRE> ErrorType;
    typedef CanonicalNumericType<P,PR,PRE> NumericType;
    typedef Box<Interval<FloatUpperBound<PR>>> RangeType;
  public:
    inline VectorFunctionModel() : _ptr() { }
    inline VectorFunctionModel(SharedPointer<const VectorFunctionModelInterface<P,PR,PRE>> vfp)
        : _ptr(vfp->_clone()) { }
    inline VectorFunctionModel(SizeType n, const ScalarFunctionModelInterface<P,PR,PRE>& sf) {
        FunctionModelFactory<P,PR,PRE> factory(sf._factory()); *this=factory.create_zeros(n,sf.domain());
        for(SizeType i=0; i!=n; ++i) { (*this)[i]=sf; } }
    inline VectorFunctionModel(Array<ScalarFunctionModel<P,PR,PRE>> const& asf)
        : VectorFunctionModel(asf.size(),asf[0]) { for(SizeType i=0; i!=asf.size(); ++i) { (*this)[i]=asf[i]; } }
    inline VectorFunctionModel(List<ScalarFunctionModel<P,PR,PRE>> const& lsf)
        : VectorFunctionModel(lsf.size(),lsf[0]) { for(SizeType i=0; i!=lsf.size(); ++i) { (*this)[i]=lsf[i]; } }
    inline VectorFunctionModel(VectorFunctionModelInterface<P,PR,PRE>* p) : _ptr(p) { }
    inline VectorFunctionModel(const VectorFunctionModelInterface<P,PR,PRE>& f) : _ptr(f._clone()) { }
    inline VectorFunctionModel(const VectorFunctionModel<P,PR,PRE>& f) : _ptr(f._ptr) { }
    inline operator const VectorFunctionModelInterface<P,PR,PRE>& () const { return *_ptr; }
    inline operator ValidatedVectorFunction () const { return ValidatedVectorFunction(*_ptr); }
    inline const VectorFunctionModelInterface<P,PR,PRE>* raw_pointer() const { return _ptr.operator->(); }
    inline const VectorFunctionModelInterface<P,PR,PRE>& reference() const { return *_ptr; }
    inline VectorFunctionModelInterface<P,PR,PRE>& reference() { return *_ptr; }
    inline SizeType result_size() const { return this->_ptr->result_size(); }
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    inline SizeType size() const { return this->_ptr->result_size(); }
    template<class XX> inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    template<class XX> inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    inline ScalarFunctionModel<P,PR,PRE> const get(SizeType i) const { return this->_ptr->_get(i); }
    inline Void set(SizeType i, ScalarFunctionModel<P,PR,PRE> const& sf) { this->_ptr->_set(i,sf); }
    inline ScalarFunctionModel<P,PR,PRE> const operator[](SizeType i) const { return this->get(i); }
    inline VectorFunctionModelElement<P,PR,PRE> operator[](SizeType i) { return VectorFunctionModelElement<P,PR,PRE>(this,i); }
    inline VectorFunctionModel<P,PR,PRE> operator[](Range rng) { VectorFunctionModel<P,PR,PRE> r=factory(*this).create_zeros(rng.size());
        for(SizeType i=0; i!=rng.size(); ++i) { r[i]=this->operator[](rng[i]); } return r; }
    inline DomainType const domain() const { return this->_ptr->domain(); }
    inline CodomainType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const { return this->_ptr->range(); }
    inline Vector<ErrorType> const errors() const { return this->_ptr->errors(); }
    inline ErrorType const error() const { return this->_ptr->error(); }
    inline Void clobber() { this->_ptr->clobber(); }
    inline Matrix<NumericType> const jacobian(const Vector<NumericType>& x) const {
        Vector<Differential<NumericType>> dfx=this->_ptr->_evaluate(Differential<NumericType>::variables(1u,x));
        return dfx.jacobian(); }

    inline Void restrict(const DomainType& d) { this->_ptr->restrict(d); }
  public:
    friend FunctionModelCreator<FunctionModelFactory<P,PR,PRE>> factory(VectorFunctionModel<P,PR,PRE> const& f) {
        FunctionModelFactory<P,PR,PRE> factory(f._ptr->_factory()); return FunctionModelCreator<FunctionModelFactory<P,PR,PRE>>(f.domain(),factory); }
  public:
    friend inline ScalarFunctionModel<P,PR,PRE> compose(const ScalarFunction<P>& f, const VectorFunctionModel<P,PR,PRE>& g) {
        return g._ptr->_compose(f); }
    friend inline ScalarFunctionModel<P,PR,PRE> compose(const ScalarFunctionModel<P,PR,PRE>& f, const VectorFunctionModel<P,PR,PRE>& g) {
        return g._ptr->_compose(f); }
    friend inline VectorFunctionModel<P,PR,PRE> compose(const VectorFunction<P>& f, const VectorFunctionModel<P,PR,PRE>& g) {
        return g._ptr->_compose(f); }
    friend inline VectorFunctionModel<P,PR,PRE> compose(const VectorFunctionModel<P,PR,PRE>& f, const VectorFunctionModel<P,PR,PRE>& g) {
        return g._ptr->_compose(f); }

    friend inline ScalarFunctionModel<P,PR,PRE> unchecked_compose(const ScalarFunctionModel<P,PR,PRE>& f, const VectorFunctionModel<P,PR,PRE>& g) {
        return g._ptr->_unchecked_compose(f); }
    friend inline VectorFunctionModel<P,PR,PRE> unchecked_compose(const VectorFunctionModel<P,PR,PRE>& f, const VectorFunctionModel<P,PR,PRE>& g) {
        return g._ptr->_unchecked_compose(f); }

    friend inline VectorFunctionModel<P,PR,PRE> operator+(const VectorFunctionModel<P,PR,PRE>& f) {
        return f._ptr->_clone(); }
    friend inline VectorFunctionModel<P,PR,PRE> operator-(const VectorFunctionModel<P,PR,PRE>& f) {
        VectorFunctionModel<P,PR,PRE> r=f; for(SizeType i=0; i!=r.size(); ++i) { r[i]=-f[i]; } return r; }
    friend inline VectorFunctionModel<P,PR,PRE> operator+(const VectorFunctionModel<P,PR,PRE>& f1, const VectorFunctionModel<P,PR,PRE>& f2) {
        VectorFunctionModel<P,PR,PRE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+f2[i]; } return r; }
    friend inline VectorFunctionModel<P,PR,PRE> operator-(const VectorFunctionModel<P,PR,PRE>& f1, const VectorFunctionModel<P,PR,PRE>& f2) {
        VectorFunctionModel<P,PR,PRE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
    friend inline VectorFunctionModel<P,PR,PRE> operator+(const VectorFunctionModel<P,PR,PRE>& f1, const Vector<CanonicalNumericType<P,PR,PRE>>& c2) {
        VectorFunctionModel<P,PR,PRE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+c2[i]; } return r; }
    friend inline VectorFunctionModel<P,PR,PRE> operator-(const VectorFunctionModel<P,PR,PRE>& f1, const Vector<CanonicalNumericType<P,PR,PRE>>& c2) {
        VectorFunctionModel<P,PR,PRE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-c2[i]; } return r; }
    friend inline VectorFunctionModel<P,PR,PRE> operator+(const Vector<CanonicalNumericType<P,PR,PRE>>& c1, const VectorFunctionModel<P,PR,PRE>& f2);
    friend inline VectorFunctionModel<P,PR,PRE> operator-(const Vector<CanonicalNumericType<P,PR,PRE>>& c1, const VectorFunctionModel<P,PR,PRE>& f2);
    friend inline VectorFunctionModel<P,PR,PRE> operator*(const VectorFunctionModel<P,PR,PRE>& f1, const CanonicalNumericType<P,PR,PRE>& c2) {
        VectorFunctionModel<P,PR,PRE> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]*c2; } return r; }
    friend inline VectorFunctionModel<P,PR,PRE> operator*(const CanonicalNumericType<P,PR,PRE>& c1, const VectorFunctionModel<P,PR,PRE>& f2) {
        VectorFunctionModel<P,PR,PRE> r=f2; for(SizeType i=0; i!=r.size(); ++i) { r[i]=c1*f2[i]; } return r; }

    friend inline VectorFunctionModel<P,PR,PRE> operator+(const VectorFunctionModel<P,PR,PRE>& f1, const VectorFunction<P>& f2) {
        return f1+factory(f1).create(f2); }
    friend inline VectorFunctionModel<P,PR,PRE> operator-(const VectorFunctionModel<P,PR,PRE>& f1, const VectorFunction<P>& f2) {
        return f1-factory(f1).create(f2); }
    friend inline VectorFunctionModel<P,PR,PRE> operator+(const VectorFunction<P>& f1, const VectorFunctionModel<P,PR,PRE>& f2) {
        return factory(f2).create(f1)+f2; }
    friend inline VectorFunctionModel<P,PR,PRE> operator-(const VectorFunction<P>& f1, const VectorFunctionModel<P,PR,PRE>& f2) {
        return factory(f2).create(f1)-f2; }


  public:
    friend NormType norm(const VectorFunctionModel<P,PR,PRE>& f) {
        return f._ptr->_norm(); }
    friend VectorFunctionModel<P,PR,PRE> embed(const DomainType& d1, const VectorFunctionModel<P,PR,PRE>& f, const DomainType& d2) {
        return f._ptr->_embed(d1,d2); }
    friend VectorFunctionModel<P,PR,PRE> embed(const VectorFunctionModel<P,PR,PRE>& f, const DomainType& d) {
        return embed(DomainType(),f,d); }
    friend VectorFunctionModel<P,PR,PRE> embed(const VectorFunctionModel<P,PR,PRE>& f, const IntervalDomainType& d) {
        return embed(f,DomainType(1,d)); }
    friend VectorFunctionModel<P,PR,PRE> restrict(const VectorFunctionModel<P,PR,PRE>& f, const DomainType& d) {
        VectorFunctionModelInterface<P,PR,PRE>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }
    friend VectorFunctionModel<P,PR,PRE> restriction(const VectorFunctionModel<P,PR,PRE>& f, const DomainType& d) {
        VectorFunctionModelInterface<P,PR,PRE>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }

    friend Vector<CanonicalNumericType<P,PR,PRE>> unchecked_evaluate(const VectorFunctionModel<P,PR,PRE>& f, const Vector<CanonicalNumericType<P,PR,PRE>>& x) {
        return f._ptr->_unchecked_evaluate(x); }

    friend ScalarFunctionModel<P,PR,PRE> unchecked_compose(const ScalarFunction<P>& f, const VectorFunctionModel<P,PR,PRE>& g) {
        ScalarFunctionModelInterface<P,PR,PRE> const* fptr = dynamic_cast<ScalarFunctionModelInterface<P,PR,PRE> const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(ScalarFunctionModel<P,PR,PRE>(*fptr),g); } else { return compose(f,g); } }
    friend VectorFunctionModel<P,PR,PRE> unchecked_compose(const VectorFunction<P>& f, const VectorFunctionModel<P,PR,PRE>& g) {
        VectorFunctionModelInterface<P,PR,PRE> const* fptr = dynamic_cast<VectorFunctionModelInterface<P,PR,PRE> const*>(f.raw_pointer());
        if(fptr) { return unchecked_compose(VectorFunctionModel<P,PR,PRE>(*fptr),g); } else { return compose(f,g); } }

    //friend VectorFunctionModel<P,PR,PRE> join(const ScalarFunctionModel<P,PR,PRE>& f1, const ScalarFunctionModel<P,PR,PRE>& f2);
    friend VectorFunctionModel<P,PR,PRE> join(const ScalarFunctionModel<P,PR,PRE>& f1, const VectorFunctionModel<P,PR,PRE>& f2) {
        return join(VectorFunctionModel<P,PR,PRE>(1u,f1),f2); }
    friend VectorFunctionModel<P,PR,PRE> join(const VectorFunctionModel<P,PR,PRE>& f1, const ScalarFunctionModel<P,PR,PRE>& f2) {
        VectorFunctionModel<P,PR,PRE> r=f1._ptr->_clone(); r._ptr->_adjoin(f2); return r; }
    friend VectorFunctionModel<P,PR,PRE> join(const VectorFunctionModel<P,PR,PRE>& f1, const VectorFunctionModel<P,PR,PRE>& f2) {
        return f1._ptr->_join(f2); }

    friend VectorFunctionModel<P,PR,PRE> combine(const ScalarFunctionModel<P,PR,PRE>& f1, const ScalarFunctionModel<P,PR,PRE>& f2) {
        return VectorFunctionModel<P,PR,PRE>(1,f1)._ptr->_combine(VectorFunctionModel<P,PR,PRE>(1,f2)); };;
    friend VectorFunctionModel<P,PR,PRE> combine(const ScalarFunctionModel<P,PR,PRE>& f1, const VectorFunctionModel<P,PR,PRE>& f2) {
        return VectorFunctionModel<P,PR,PRE>(1,f1)._ptr->_combine(f2); };;
    friend VectorFunctionModel<P,PR,PRE> combine(const VectorFunctionModel<P,PR,PRE>& f1, const ScalarFunctionModel<P,PR,PRE>& f2) {
        return f1._ptr->_combine(VectorFunctionModel<P,PR,PRE>(1,f2)); };
    friend VectorFunctionModel<P,PR,PRE> combine(const VectorFunctionModel<P,PR,PRE>& f1, const VectorFunctionModel<P,PR,PRE>& f2) {
        return f1._ptr->_combine(f2); }

    friend inline VectorFunctionModel<ValidatedTag,PR,PRE> refinement(const VectorFunctionModel<ValidatedTag,PR,PRE>& f1, const VectorFunctionModel<ValidatedTag,PR,PRE>& f2) {
        ARIADNE_ASSERT_MSG(f1.size()==f2.size(),"refinement(f1,f2): f1="<<f1<<", f2="<<f2<<")");
        VectorFunctionModel<ValidatedTag,PR,PRE> r=+f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=refinement(f1[i],f2[i]); } return r; }

    friend VectorFunctionModel<P,PR,PRE> antiderivative(const VectorFunctionModel<P,PR,PRE>& f, SizeType j) {
        VectorFunctionModel<P,PR,PRE> r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j); }
        return r;
    }

    friend VectorFunctionModel<P,PR,PRE> antiderivative(const VectorFunctionModel<P,PR,PRE>& f, SizeType j, const CanonicalNumericType<P,PR,PRE>& c) {
        VectorFunctionModel<P,PR,PRE> r(f);
        for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j,c); }
        return r;
    }

    friend VectorFunctionModel<P,PR,PRE> partial_evaluate(const VectorFunctionModel<P,PR,PRE>& f, SizeType j, const CanonicalNumericType<P,PR,PRE>& c) {
        return f._ptr->_partial_evaluate(j,c); }

    friend OutputStream& operator<<(OutputStream& os, const VectorFunctionModel<P,PR,PRE>& f) {
        return os <<  f.operator VectorFunction<P>(); }

};

template<class P, class PR, class PRE> inline CanonicalNumericType<P,PR,PRE> evaluate(const ScalarFunctionModel<P,PR,PRE>& f, const Vector<CanonicalNumericType<P,PR,PRE>>& x) {
    return f(x); }
template<class P, class PR, class PRE> inline CanonicalNumericType<P,PR,PRE> unchecked_evaluate(const ScalarFunctionModel<P,PR,PRE>& f, const Vector<CanonicalNumericType<P,PR,PRE>>& x) {
    return f._ptr->_unchecked_evaluate(x); }

template<class P> inline CanonicalNumeric64Type<P> unchecked_evaluate(const ScalarFunction<P>& f, const Vector<CanonicalNumeric64Type<P>>& x) {
    ScalarFunctionModel64Interface<P> const* fptr = dynamic_cast<ScalarFunctionModel64Interface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(ScalarFunctionModel64<P>(*fptr),x); } else { return evaluate(f,x); } }

template<class P> inline Vector<CanonicalNumeric64Type<P>> unchecked_evaluate(const VectorFunction<P>& f, const Vector<CanonicalNumeric64Type<P>>& x) {
    VectorFunctionModel64Interface<P> const* fptr = dynamic_cast<VectorFunctionModel64Interface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(VectorFunctionModel64<P>(*fptr),x); } else { return evaluate(f,x); } }

template<class P, class PR, class PRE> inline ScalarFunctionModel<P,PR,PRE> unchecked_compose(const ScalarFunction<P>& f, const VectorFunctionModel<P,PR,PRE>& g) {
    ScalarFunctionModelInterface<P,PR,PRE> const* fptr = dynamic_cast<ScalarFunctionModelInterface<P,PR,PRE> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(ScalarFunctionModel<P,PR,PRE>(*fptr),g); } else { return compose(f,g); } }
template<class P, class PR, class PRE> inline VectorFunctionModel<P,PR,PRE> unchecked_compose(const VectorFunction<P>& f, const VectorFunctionModel<P,PR,PRE>& g) {
    VectorFunctionModelInterface<P,PR,PRE> const* fptr = dynamic_cast<VectorFunctionModelInterface<P,PR,PRE> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(VectorFunctionModel<P,PR,PRE>(*fptr),g); } else { return compose(f,g); } }

// FIXME: Should be unneeded
template<class PR, class PRE> ScalarFunctionModel<ValidatedTag,PR,PRE> unchecked_compose(const ScalarFunctionModel<ValidatedTag,PR,PRE>& f, const VectorFunctionModel<ValidatedTag,PR,PRE>& g) {
    return g._ptr->_unchecked_compose(f); }
template<class PR, class PRE> VectorFunctionModel<ValidatedTag,PR,PRE> unchecked_compose(const VectorFunctionModel<ValidatedTag,PR,PRE>& f, const VectorFunctionModel<ValidatedTag,PR,PRE>& g) {
    return g._ptr->_unchecked_compose(f); }
template<class PR, class PRE> ScalarFunctionModel<ValidatedTag,PR,PRE> restrict(const ScalarFunctionModel<ValidatedTag,PR,PRE>& f, const BoxDomainType& dom) {
    return f._ptr->_restriction(dom); }
template<class PR, class PRE> VectorFunctionModel<ValidatedTag,PR,PRE> restrict(const VectorFunctionModel<ValidatedTag,PR,PRE>& f, const BoxDomainType& dom) {
    return f._ptr->_restriction(dom); }


// Not in function_model_interface.h since DomainType (ExactFloat64Box) is undefined.
template<class P, class PR, class PRE> ScalarFunctionModel<P,PR,PRE> FunctionModelFactoryInterface<P,PR,PRE>::create_identity(const IntervalDomainType& domain) const {
    return _create_coordinate(DomainType(1u,domain),0u); };

template<class P, class PR, class PRE> VectorFunctionModelElement<P,PR,PRE>::operator const ScalarFunctionModel<P,PR,PRE> () const {
    return _p->get(_i); }





// Full output
template<class T> struct Representation { const T* pointer; Representation(const T& t) : pointer(&t) { } const T& reference() const { return *pointer; } };
template<class T> inline Representation<T> representation(const T& t) { return Representation<T>(t); }

template<class T> class HasRepr {
    template<class TT, class = decltype(declval<TT>().repr(declval<OutputStream&>()))> static True test(int);
    template<class TT> static False test(...);
  public:
    static const bool value = decltype(test<T>(1))::value;
};

template<class T, EnableIf<HasRepr<T>> =dummy> inline OutputStream& operator<<(OutputStream& os, const Representation<T>& obj) { obj.reference().repr(os); return os; }
template<class T, DisableIf<HasRepr<T>> =dummy> inline OutputStream& operator<<(OutputStream& os, const Representation<T>& obj) { return os << obj.reference(); }

} // namespace Ariadne

#endif
