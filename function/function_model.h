/***************************************************************************
 *            function_model.h
 *
 *  Copyright 2011  Pieter Collins
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
#include "geometry/box.h"

#include "function/function_interface.h"
#include "function/function_mixin.h"
#include "function/function.h"

namespace Ariadne {

template<class P, class PR=Precision64, class PRE=PR> class FunctionModelBuilder {
    ScalarFunctionModel<P,PR,PRE> const& _prototype;
  public:
    FunctionModelBuilder(ScalarFunctionModel<P,PR,PRE> const& f);
    FunctionModelBuilder(VectorFunctionModel<P,PR,PRE> const& f);

    ScalarFunctionModel<P,PR,PRE> create(ScalarFunction<P> const& f);
    VectorFunctionModel<P,PR,PRE> create(VectorFunction<P> const& f);
    CanonicalNumericType<P,PR,PRE> create(Number<P> const& c);

    ScalarFunctionModel<P,PR,PRE> const& create(ScalarFunctionModel<P,PR,PRE> const& f) { return f; }
    CanonicalNumericType<P,PR,PRE> const& create(CanonicalNumericType<P,PR,PRE> const& c) { return c; }
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
    typedef ExactBoxType DomainType;
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
    ScalarFunctionModel<P,PR,PRE> create_zero() const;
    ScalarFunctionModel<P,PR,PRE> create_constant(const Number<P>& c) const;
    ScalarFunctionModel<P,PR,PRE> create_coordinate(SizeType j) const;
    VectorFunctionModel<P,PR,PRE> create_identity() const;
    ScalarFunctionModel<P,PR,PRE> create(const ScalarFunction<P>& f) const;
    VectorFunctionModel<P,PR,PRE> create(const VectorFunction<P>& f) const;
    CanonicalNumericType<P,PR,PRE> create(const Number<P>& f) const;
    Vector<ScalarFunctionModel<P,PR,PRE>> create_coordinates(DomainType const&) const;
    ScalarFunctionModel<P,PR,PRE>& operator=(const Number<P>& c);
    ScalarFunctionModel<P,PR,PRE>& operator=(const CanonicalNumericType<P,PR,PRE>& c);
    ScalarFunctionModel<P,PR,PRE>& operator=(const ScalarFunction<P>& f);
    ScalarFunctionModel<P,PR,PRE>& operator=(const ScalarFunctionModelInterface<P,PR,PRE>& f);
//    ScalarFunctionModel<P,PR,PRE>& operator=(const ScalarTaylorFunction& f);
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    template<class X> X operator() (const Vector<X>& x) const {
        return this->_ptr->_evaluate(x); }
    template<class X> X evaluate(const Vector<X>& x) const {
        return this->_ptr->_evaluate(x); }
    inline ExactBoxType const domain() const { return this->_ptr->domain(); }
    inline ExactIntervalType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const { return this->_ptr->range(); }

    inline CoefficientType value() const { return this->_ptr->value(); }
    inline CoefficientType gradient_value(SizeType j) const { return this->_ptr->gradient_value(j); }
    inline ErrorType error() const { return this->_ptr->error(); }

    inline Void set_error(const ErrorType& e) { return this->_ptr->set_error(e); }
    inline Void set_error(Nat e) { return this->_ptr->set_error(ErrorType(e,this->error().precision())); }
    inline Void clobber() { return this->_ptr->clobber(); }

    inline ScalarFunctionModel<P,PR,PRE> apply(Operator op) const { return this->_ptr->_apply(op); }
    inline Void restrict(const ExactBoxType& d) { *this=restriction(*this,d); }
  public:
    friend FunctionModelBuilder<P,PR,PRE> factory(ScalarFunctionModel<P,PR,PRE> const& f) {
        return FunctionModelBuilder<P,PR,PRE>(f); }
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
        ScalarFunctionModel<P,PR,PRE> r=f1.create_zero(); r._ptr->_ifma(f1,f2); return r; }
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

    friend ScalarFunctionModel<P,PR,PRE> embed(const ExactBoxType& d1, const ScalarFunctionModel<P,PR,PRE>& f, const ExactBoxType& d2) {
        return f._ptr->_embed(d1,d2); }
    friend ScalarFunctionModel<P,PR,PRE> embed(const ExactBoxType& d, const ScalarFunctionModel<P,PR,PRE>& f) {
        return embed(d,f,ExactBoxType()); }
    friend ScalarFunctionModel<P,PR,PRE> embed(const ScalarFunctionModel<P,PR,PRE>& f, const ExactBoxType& d) {
        return embed(ExactBoxType(),f,d); }
    friend ScalarFunctionModel<P,PR,PRE> embed(const ScalarFunctionModel<P,PR,PRE>& f, const ExactIntervalType& d) {
        return embed(f,ExactBoxType(1,d)); }
    friend ScalarFunctionModel<P,PR,PRE> restrict(const ScalarFunctionModel<P,PR,PRE>& f, const ExactBoxType& d) {
        return f._ptr->_restriction(d); }
    friend ScalarFunctionModel<P,PR,PRE> restriction(const ScalarFunctionModel<P,PR,PRE>& f, const ExactBoxType& d) {
        return f._ptr->_restriction(d); }

    friend VectorFunctionModel<P,PR,PRE> join(const ScalarFunctionModel<P,PR,PRE>& f1, const ScalarFunctionModel<P,PR,PRE>& f2) {
        return join(VectorFunctionModel<P,PR,PRE>(1,f1),f2); }
    friend VectorFunctionModel<P,PR,PRE> combine(const ScalarFunctionModel<P,PR,PRE>& f1, const ScalarFunctionModel<P,PR,PRE>& f2);
  public:
    friend OutputStream& operator<<(OutputStream& os, const ScalarFunctionModel<P,PR,PRE>& f) {
        return os <<  f.operator ScalarFunction<P>(); }
};

inline ScalarFunctionModel<ValidatedTag> refinement(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_refinement(f2); }
inline Boolean inconsistent(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_inconsistent(f2); }
inline Boolean refines(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_refines(f2); }

// FIXME: Should not be needed since ScalarFunctionModel has a representation
template<class P> inline ScalarFunctionModel<P> embed(const ScalarFunction<P>& f, const ExactIntervalType& d) {
    return embed(ScalarFunctionModel<P>(f),d); }

template<class P, class PR, class PRE> inline ScalarFunctionModel<P,PR,PRE>& ScalarFunctionModel<P,PR,PRE>::operator=(const CanonicalNumericType<P,PR,PRE>& c) {
        (*this)*=CanonicalNumericType<P,PR,PRE>(0); (*this)+=c; return *this; }
template<class P, class PR, class PRE> inline ScalarFunctionModel<P,PR,PRE>& ScalarFunctionModel<P,PR,PRE>::operator=(const Number<P>& c) {
        return (*this)=factory(*this).create(c); }
template<class P, class PR, class PRE> inline ScalarFunctionModel<P,PR,PRE>& ScalarFunctionModel<P,PR,PRE>::operator=(const ScalarFunction<P>& f) {
        return (*this)=factory(*this).create(f); }
template<class P, class PR, class PRE> inline ScalarFunctionModel<P,PR,PRE>& ScalarFunctionModel<P,PR,PRE>::operator=(const ScalarFunctionModelInterface<P,PR,PRE>& f) {
        ARIADNE_NOT_IMPLEMENTED; }



template<class V> struct Element;

template<class M> class FunctionPatch;
template<class M> class VectorFunctionPatch;
template<class M> struct Element<VectorFunctionPatch<M>> { typedef FunctionPatch<M> Type; };

typedef FunctionPatch<ValidatedTaylorModel64> ScalarTaylorFunction;

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
    typedef CanonicalCoefficientType<P,PR> CoefficientType;
    typedef CanonicalErrorType<P,PRE> ErrorType;
    typedef CanonicalNumericType<P,PR,PRE> NumericType;
    typedef Box<Interval<FloatUpperBound<PR>>> RangeType;
  public:
    inline VectorFunctionModel() : _ptr() { }
    inline VectorFunctionModel(SharedPointer<const VectorFunctionModelInterface<P,PR,PRE>> vfp)
        : _ptr(vfp->_clone()) { }
    inline VectorFunctionModel(SizeType n, const ScalarFunctionModelInterface<P,PR,PRE>& sf)
        : _ptr(sf._create_vector(n)) { for(SizeType i=0; i!=n; ++i) { (*this)[i]=sf; } }
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
    inline ScalarFunctionModel<P,PR,PRE> create_zero() const { return this->_ptr->_create_zero(); }
    inline VectorFunctionModel<P,PR,PRE> create_identity() const { return this->_ptr->_create_identity(); }
    inline VectorFunctionModel<P,PR,PRE> create(VectorFunction<P> const& vf) const { return compose(vf,this->create_identity()); }
    inline SizeType result_size() const { return this->_ptr->result_size(); }
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    inline SizeType size() const { return this->_ptr->result_size(); }
    template<class XX> inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    template<class XX> inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    inline ScalarFunctionModel<P,PR,PRE> const get(SizeType i) const { return this->_ptr->_get(i); }
    inline Void set(SizeType i, ScalarFunctionModel<P,PR,PRE> const& sf) { this->_ptr->_set(i,sf); }
    inline ScalarFunctionModel<P,PR,PRE> const operator[](SizeType i) const { return this->get(i); }
    inline VectorFunctionModelElement<P,PR,PRE> operator[](SizeType i) { return VectorFunctionModelElement<P,PR,PRE>(this,i); }
    inline VectorFunctionModel<P,PR,PRE> operator[](Range rng) { VectorFunctionModel<P,PR,PRE> r(rng.size(),this->create_zero());
        for(SizeType i=0; i!=rng.size(); ++i) { r[i]=this->operator[](rng[i]); } return r; }
    inline ExactBoxType const domain() const { return this->_ptr->domain(); }
    inline ExactBoxType const codomain() const { return this->_ptr->codomain(); }
    inline RangeType const range() const { return this->_ptr->range(); }
    inline Vector<ErrorType> const errors() const { return this->_ptr->errors(); }
    inline ErrorType const error() const { return this->_ptr->error(); }
    inline Void clobber() { this->_ptr->clobber(); }
    inline Matrix<NumericType> const jacobian(const Vector<NumericType>& x) const {
        Vector<Differential<NumericType>> dfx=this->_ptr->_evaluate(Differential<NumericType>::variables(1u,x));
        return dfx.jacobian(); }

    inline Void restrict(const ExactBoxType& d) { this->_ptr->restrict(d); }
  public:
    friend FunctionModelBuilder<P,PR,PRE> factory(VectorFunctionModel<P,PR,PRE> const& f) {
        return FunctionModelBuilder<P,PR,PRE>(f); }
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
        return f1+f1.create(f2); }
    friend inline VectorFunctionModel<P,PR,PRE> operator-(const VectorFunctionModel<P,PR,PRE>& f1, const VectorFunction<P>& f2) {
        return f1-f1.create(f2); }
    friend inline VectorFunctionModel<P,PR,PRE> operator+(const VectorFunction<P>& f1, const VectorFunctionModel<P,PR,PRE>& f2) {
        return f2.create(f1)+f2; }
    friend inline VectorFunctionModel<P,PR,PRE> operator-(const VectorFunction<P>& f1, const VectorFunctionModel<P,PR,PRE>& f2) {
        return f2.create(f1)-f2; }


  public:
    friend NormType norm(const VectorFunctionModel<P,PR,PRE>& f) {
        return f._ptr->_norm(); }
    friend VectorFunctionModel<P,PR,PRE> embed(const ExactBoxType& d1, const VectorFunctionModel<P,PR,PRE>& f, const ExactBoxType& d2) {
        return f._ptr->_embed(d1,d2); }
    friend VectorFunctionModel<P,PR,PRE> embed(const VectorFunctionModel<P,PR,PRE>& f, const ExactBoxType& d) {
        return embed(ExactBoxType(),f,d); }
    friend VectorFunctionModel<P,PR,PRE> embed(const VectorFunctionModel<P,PR,PRE>& f, const ExactIntervalType& d) {
        return embed(f,ExactBoxType(1,d)); }
    friend VectorFunctionModel<P,PR,PRE> restrict(const VectorFunctionModel<P,PR,PRE>& f, const ExactBoxType& d) {
        VectorFunctionModelInterface<P,PR,PRE>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }
    friend VectorFunctionModel<P,PR,PRE> restriction(const VectorFunctionModel<P,PR,PRE>& f, const ExactBoxType& d) {
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


template<class P> inline CanonicalNumericType<P> unchecked_evaluate(const ScalarFunction<P>& f, const Vector<CanonicalNumericType<P>>& x) {
    ScalarFunctionModelInterface<P> const* fptr = dynamic_cast<ScalarFunctionModelInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(ScalarFunctionModel<P>(*fptr),x); } else { return evaluate(f,x); } }

template<class P> inline Vector<CanonicalNumericType<P>> unchecked_evaluate(const VectorFunction<P>& f, const Vector<CanonicalNumericType<P>>& x) {
    VectorFunctionModelInterface<P> const* fptr = dynamic_cast<VectorFunctionModelInterface<P> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(VectorFunctionModel<P>(*fptr),x); } else { return evaluate(f,x); } }

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
template<class PR, class PRE> ScalarFunctionModel<ValidatedTag,PR,PRE> restrict(const ScalarFunctionModel<ValidatedTag,PR,PRE>& f, const ExactBoxType& dom) {
    return f._ptr->_restriction(dom); }
template<class PR, class PRE> VectorFunctionModel<ValidatedTag,PR,PRE> restrict(const VectorFunctionModel<ValidatedTag,PR,PRE>& f, const ExactBoxType& dom) {
    return f._ptr->_restriction(dom); }



template<class P, class PR, class PRE> VectorFunctionModelElement<P,PR,PRE>::operator const ScalarFunctionModel<P,PR,PRE> () const {
    return _p->get(_i); }



template<class P, class PR, class PRE> inline FunctionModelBuilder<P,PR,PRE>::FunctionModelBuilder(ScalarFunctionModel<P,PR,PRE> const& f)
    : _prototype(f) { }
template<class P, class PR, class PRE> inline FunctionModelBuilder<P,PR,PRE>::FunctionModelBuilder(VectorFunctionModel<P,PR,PRE> const& f)
    : _prototype(f.create_zero()) { }
template<class P, class PR, class PRE> inline ScalarFunctionModel<P,PR,PRE> FunctionModelBuilder<P,PR,PRE>::create(ScalarFunction<P> const& f) {
    return _prototype.create(f); }
template<class P, class PR, class PRE> inline VectorFunctionModel<P,PR,PRE> FunctionModelBuilder<P,PR,PRE>::create(VectorFunction<P> const& f) {
    return _prototype.create(f); }
template<class P, class PR, class PRE> inline CanonicalNumericType<P,PR,PRE> FunctionModelBuilder<P,PR,PRE>::create(Number<P> const& c) {
    return CanonicalNumericType<P,PR,PRE>(c,_prototype.value().precision()); }


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
