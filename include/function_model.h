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

#include "function_interface.h"
#include "function_mixin.h"
#include "function.h"

#include "operators.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"

namespace Ariadne {


template<class X> class ScalarFunctionModelInterface;
template<class F,class X> class ScalarFunctionModelMixin;
template<class X> class ScalarFunctionModel;

template<class X> class VectorFunctionModelInterface;
template<class F,class X> class VectorFunctionModelMixin;
template<class X> class VectorFunctionModel;

typedef ScalarFunctionModelInterface<ValidatedTag> ValidatedScalarFunctionModelInterface;
typedef VectorFunctionModelInterface<ValidatedTag> ValidatedVectorFunctionModelInterface;

typedef ScalarFunctionModel<ValidatedTag> ValidatedScalarFunctionModel;
typedef VectorFunctionModel<ValidatedTag> ValidatedVectorFunctionModel;

template<class X> class FunctionModelFactory;
typedef FunctionModelFactory<ValidatedTag> ValidatedFunctionModelFactory;

class ScalarTaylorFunction;
class VectorTaylorFunction;

template<> class ScalarFunctionModelInterface<ValidatedTag>
    : public virtual ScalarFunctionInterface<ValidatedTag>
{
  public:
    virtual ExactBox const& domain() const = 0;
    virtual ExactInterval const codomain() const = 0;
    virtual UpperInterval range() const = 0;

    virtual CoefficientType const& value() const = 0;
    virtual CoefficientType const gradient_value(Nat i) const = 0;
    virtual ErrorType const& error() const = 0;

    virtual Void set_error(const ErrorType& e) = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual ScalarFunctionModelInterface<ValidatedTag>* _apply(OperatorCode op) const = 0;
    virtual ValidatedNumber _unchecked_evaluate(const Vector<ValidatedNumber>& x) const = 0;

    virtual ScalarFunctionModelInterface<ValidatedTag>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create() const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create_vector(Nat i) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _embed(const ExactBox& d1, const ExactBox& d2) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create_identity() const = 0;
    virtual Void restrict(const ExactBox& d) = 0;

    virtual ScalarFunctionModelInterface<ValidatedTag>* _derivative(Nat j) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _antiderivative(Nat j) const = 0;

    virtual Tribool _refines(const ScalarFunctionModelInterface<ValidatedTag>& f) const = 0;
    virtual Tribool _disjoint(const ScalarFunctionModelInterface<ValidatedTag>& f) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _intersection(const ScalarFunctionModelInterface<ValidatedTag>& f) const = 0;

    virtual Void _iadd(const ValidatedNumber& c) = 0;
    virtual Void _imul(const ValidatedNumber& c) = 0;
    virtual Void _isma(const ValidatedNumber& c, const ScalarFunctionModelInterface<ValidatedTag>& f) = 0;
    virtual Void _ifma(const ScalarFunctionModelInterface<ValidatedTag>& f1, const ScalarFunctionModelInterface<ValidatedTag>& f2) = 0;
};


template<class F> class ScalarFunctionModelMixin<F,ValidatedTag>
    : public virtual ScalarFunctionModelInterface<ValidatedTag>
    , public ScalarFunctionMixin<F,ValidatedTag>
{
  public:
    F apply(OperatorCode op) const;
  public:
    ScalarFunctionModelInterface<ValidatedTag>* _clone() const {
        return new F(static_cast<const F&>(*this)); }
    NormType const _norm() const {
        return norm(static_cast<const F&>(*this)); }
    ScalarFunctionModelInterface<ValidatedTag>* _antiderivative(Nat j) const {
        return new F(antiderivative(static_cast<const F&>(*this),j)); }
    ScalarFunctionModelInterface<ValidatedTag>* _apply(OperatorCode op) const {
        return new F(this->apply(op)); }
    ValidatedNumber _unchecked_evaluate(const Vector<ValidatedNumber>& x) const {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<ValidatedTag>* _embed(const ExactBox& d1, const ExactBox& d2) const {
        return new F(embed(d1,static_cast<const F&>(*this),d2)); }
    Tribool _refines(const ScalarFunctionModelInterface<ValidatedTag>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return refines(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    Tribool _disjoint(const ScalarFunctionModelInterface<ValidatedTag>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return disjoint(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    ScalarFunctionModelInterface<ValidatedTag>* _intersection(const ScalarFunctionModelInterface<ValidatedTag>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return new F(intersection(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Void _iadd(const ValidatedNumber& c) {
        static_cast<F&>(*this)+=c; }
    Void _imul(const ValidatedNumber& c) {
        static_cast<F&>(*this)*=c; }
    Void _isma(const ValidatedNumber& c, const ScalarFunctionModelInterface<ValidatedTag>& f) {
        static_cast<F&>(*this)+=c*dynamic_cast<const F&>(f); }
    Void _ifma(const ScalarFunctionModelInterface<ValidatedTag>& f1, const ScalarFunctionModelInterface<ValidatedTag>& f2) {
        static_cast<F&>(*this)+=dynamic_cast<const F&>(f1)*dynamic_cast<const F&>(f2); }
};

//! \ingroup FunctionModelSubModule
//! \brief Generic scalar functions on bounded domains.
template<> class ScalarFunctionModel<ValidatedTag>
{
  public:
    clone_on_copy_ptr< ScalarFunctionModelInterface<ValidatedTag> > _ptr;
  public:
    ScalarFunctionModel() : _ptr() { }
    ScalarFunctionModel(ScalarFunctionModelInterface<ValidatedTag>* p) : _ptr(p) { }
    ScalarFunctionModel(const ScalarFunctionModel<ValidatedTag>& f) : _ptr(f._ptr) { }
    ScalarFunctionModel(const ScalarFunctionModelInterface<ValidatedTag>& f) : _ptr(f._clone()) { }
    ScalarFunctionModel(const ValidatedScalarFunction& f) : _ptr(dynamic_cast<ScalarFunctionModelInterface<ValidatedTag>*>(f.raw_pointer()->_clone())) { }
    operator ValidatedScalarFunction() const { return ValidatedScalarFunction(this->_ptr->_clone()); }
    operator ScalarFunctionModelInterface<ValidatedTag>& () { return *_ptr; }
    operator const ScalarFunctionModelInterface<ValidatedTag>& () const { return *_ptr; }
    const ScalarFunctionModelInterface<ValidatedTag>* raw_pointer() const { return _ptr.operator->(); }
    ScalarFunctionModelInterface<ValidatedTag>& reference() { return *_ptr; }
    const ScalarFunctionModelInterface<ValidatedTag>& reference() const { return *_ptr; }
    ScalarFunctionModel<ValidatedTag> create_zero() const { return this->_ptr->_create(); }
    VectorFunctionModel<ValidatedTag> create_identity() const;
    ScalarFunctionModel<ValidatedTag>& operator=(const ValidatedNumber& c);
    ScalarFunctionModel<ValidatedTag>& operator=(const ValidatedScalarFunction& f);
    ScalarFunctionModel<ValidatedTag>& operator=(const ScalarTaylorFunction& f);
    inline Nat argument_size() const { return this->_ptr->argument_size(); }
    template<class XX> inline XX operator()(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    template<class XX> inline XX evaluate(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    inline ExactBox const domain() const { return this->_ptr->domain(); }
    inline ExactInterval const codomain() const { return this->_ptr->codomain(); }
    inline UpperInterval const range() const { return this->_ptr->range(); }

    inline CoefficientType value() const { return this->_ptr->value(); }
    inline CoefficientType gradient_value(Nat j) const { return this->_ptr->gradient_value(j); }
    inline ErrorType error() const { return this->_ptr->error(); }

    inline Void set_error(const ErrorType& e) { return this->_ptr->set_error(e); }
    inline Void clobber() { return this->_ptr->clobber(); }

    inline ScalarFunctionModel<ValidatedTag> apply(Operator op) const { return this->_ptr->_apply(op); }
    inline Void restrict(const ExactBox& d) { this->_ptr->restrict(d); }
};

// inline ScalarFunctionModel<ValidatedTag>& ScalarFunctionModel<ValidatedTag>::operator=(const ValidatedScalarFunction& f) { (*this)=this->_ptr->_create()+f; return *this; }

inline NormType norm(const ScalarFunctionModel<ValidatedTag>& f) { return f._ptr->_norm(); }
inline ScalarFunctionModel<ValidatedTag> derivative(const ScalarFunctionModel<ValidatedTag>& f, Nat j) { return f._ptr->_derivative(j); }
inline ScalarFunctionModel<ValidatedTag> antiderivative(const ScalarFunctionModel<ValidatedTag>& f, Nat j) { return f._ptr->_antiderivative(j); }

inline ScalarFunctionModel<ValidatedTag> embed(const ExactBox& d1, const ScalarFunctionModel<ValidatedTag>& f, const ExactBox& d2) {
    return f._ptr->_embed(d1,d2); }
inline ScalarFunctionModel<ValidatedTag> embed(const ExactBox& d, const ScalarFunctionModel<ValidatedTag>& f) {
    return embed(d,f,ExactBox()); }
inline ScalarFunctionModel<ValidatedTag> embed(const ScalarFunctionModel<ValidatedTag>& f, const ExactBox& d) {
    return embed(ExactBox(),f,d); }
inline ScalarFunctionModel<ValidatedTag> embed(const ScalarFunctionModel<ValidatedTag>& f, const ExactInterval& d) {
    return embed(f,ExactBox(1,d)); }
inline ScalarFunctionModel<ValidatedTag> restrict(const ScalarFunctionModel<ValidatedTag>& f, const ExactBox& d) {
    ScalarFunctionModelInterface<ValidatedTag>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }
inline ScalarFunctionModel<ValidatedTag> intersection(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_intersection(f2); }

inline Tribool disjoint(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_disjoint(f2); }
inline Tribool refines(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_refines(f2); }

inline ScalarFunctionModel<ValidatedTag>& operator+=(ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) { f1._ptr->_isma(ValidatedTag(+1),f2); return  f1; }
inline ScalarFunctionModel<ValidatedTag>& operator-=(ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) { f1._ptr->_isma(ValidatedTag(-1),f2); return  f1; }
inline ScalarFunctionModel<ValidatedTag>& operator+=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) { f1._ptr->_iadd(c2); return f1; }
inline ScalarFunctionModel<ValidatedTag>& operator-=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) { f1._ptr->_iadd(neg(c2)); return f1; }
inline ScalarFunctionModel<ValidatedTag>& operator*=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) { f1._ptr->_imul(c2); return f1; }
inline ScalarFunctionModel<ValidatedTag>& operator/=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) { f1._ptr->_imul(rec(c2)); return f1; }

inline ScalarFunctionModel<ValidatedTag> neg(const ScalarFunctionModel<ValidatedTag>& f) { return f.apply(Neg()); }
inline ScalarFunctionModel<ValidatedTag> rec(const ScalarFunctionModel<ValidatedTag>& f) { return f.apply(Rec()); }

inline ScalarFunctionModel<ValidatedTag> operator+(const ScalarFunctionModel<ValidatedTag>& f) {
    return f._ptr->_clone(); }
inline ScalarFunctionModel<ValidatedTag> operator-(const ScalarFunctionModel<ValidatedTag>& f) {
    ScalarFunctionModel<ValidatedTag> r=f; r*=-1; return r; }
inline ScalarFunctionModel<ValidatedTag> operator+(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    ScalarFunctionModel<ValidatedTag> r=f1; r+=f2; return r; }
inline ScalarFunctionModel<ValidatedTag> operator-(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    ScalarFunctionModel<ValidatedTag> r=f1; r-=f2; return r; }
inline ScalarFunctionModel<ValidatedTag> operator*(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    ScalarFunctionModel<ValidatedTag> r=f1.create_zero(); r._ptr->_ifma(f1,f2); return r; }
inline ScalarFunctionModel<ValidatedTag> operator/(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1*rec(f2); }

inline ScalarFunctionModel<ValidatedTag> operator+(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) {
    ScalarFunctionModel<ValidatedTag> r=f1; r+=c2; return r; }
inline ScalarFunctionModel<ValidatedTag> operator-(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) {
    ScalarFunctionModel<ValidatedTag> r=f1; r-=c2; return r; }
inline ScalarFunctionModel<ValidatedTag> operator*(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) {
    ScalarFunctionModel<ValidatedTag> r=f1; r*=c2; return r; }
inline ScalarFunctionModel<ValidatedTag> operator/(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) {
    ScalarFunctionModel<ValidatedTag> r=f1; r/=c2; return r; }
inline ScalarFunctionModel<ValidatedTag> operator+(const ValidatedNumber& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
    ScalarFunctionModel<ValidatedTag> r=f2; r+=c1; return r; }
inline ScalarFunctionModel<ValidatedTag> operator-(const ValidatedNumber& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
    ScalarFunctionModel<ValidatedTag> r=neg(f2); r+=c1; return r; }
inline ScalarFunctionModel<ValidatedTag> operator*(const ValidatedNumber& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
    ScalarFunctionModel<ValidatedTag> r=f2; r*=c1; return r; }
inline ScalarFunctionModel<ValidatedTag> operator/(const ValidatedNumber& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
    ScalarFunctionModel<ValidatedTag> r=rec(f2); r*=c1; return r; }

inline ScalarFunctionModel<ValidatedTag>& ScalarFunctionModel<ValidatedTag>::operator=(const ValidatedNumber& c) { (*this)*=0.0; (*this)+=c; return *this; }

template<class F> F ScalarFunctionModelMixin<F,ValidatedTag>::apply(OperatorCode op) const {
    const F& f=static_cast<const F&>(*this);
    switch(op) {
        case NEG: return neg(f);
        case REC: return rec(f);
        default: ARIADNE_FAIL_MSG("ScalarFunctionModel<ValidatedTag>::apply(OperatorCode op): Operator op="<<op<<" not implemented\n");
    }
}

//inline ScalarFunctionModel<ValidatedTag>::ScalarFunctionModel(const ValidatedScalarFunction& f)
//    : _ptr(dynamic_cast<const ScalarFunctionModelInterface<ValidatedTag>&>(*f.raw_pointer())._clone()) { }




//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on bounded domains.
template<> class VectorFunctionModelInterface<ValidatedTag>
    : public virtual VectorFunctionInterface<ValidatedTag>
{
  public:
    virtual ExactBox const& domain() const = 0;
    virtual ExactBox const codomain() const = 0;
    virtual UpperBox const range() const = 0;
    virtual Vector<ErrorType> const errors() const = 0;
    virtual ErrorType const error() const = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual VectorFunctionModelInterface<ValidatedTag>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create_zero() const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create_identity() const = 0;
    virtual Void _set(Nat, ScalarFunctionModelInterface<ValidatedTag> const&) = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _get(Nat) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _embed(const ExactBox& d1, const ExactBox& d2) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _join(const VectorFunctionModelInterface<ValidatedTag>& f2) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _combine(const VectorFunctionModelInterface<ValidatedTag>& f2) const = 0;
    virtual Void _adjoin(const ScalarFunctionModelInterface<ValidatedTag>& f2) = 0;
    virtual Vector<ValidatedNumber> _unchecked_evaluate(const Vector<ValidatedNumber>& x) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _compose(const ScalarFunctionInterface<ValidatedTag>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _compose(const VectorFunctionInterface<ValidatedTag>& f) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _unchecked_compose(const ScalarFunctionInterface<ValidatedTag>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _unchecked_compose(const VectorFunctionInterface<ValidatedTag>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _partial_evaluate(Nat j, const ValidatedNumber& c) const = 0;
    virtual Void restrict(const ExactBox& d) = 0;
};

class ScalarTaylorFunction;
class VectorTaylorFunction;
template<class V> struct Element;
template<> struct Element<VectorTaylorFunction> { typedef ScalarTaylorFunction Type; };

template<class F> class VectorFunctionModelMixin<F,ValidatedTag>
    : public virtual VectorFunctionModelInterface<ValidatedTag>
    , public  VectorFunctionMixin<F,ValidatedTag>
{
    typedef typename Element<F>::Type ScalarFunctionType;
  public:
    virtual VectorFunctionModelInterface<ValidatedTag>* _clone() const { return new F(static_cast<const F&>(*this)); }
    virtual Void _set(Nat i, const ScalarFunctionModelInterface<ValidatedTag>& sf) {
        if(!dynamic_cast<const typename F::ScalarFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorFunctionModel "<<*this<<" to "<<sf<<"\n"); }
        static_cast<F&>(*this).F::set(i,dynamic_cast<const ScalarFunctionType&>(sf)); }
    NormType const _norm() const {
         return norm(static_cast<const F&>(*this)); }
    VectorFunctionModelInterface<ValidatedTag>* _embed(const ExactBox& d1, const ExactBox& d2) const {
        return heap_copy(embed(d1,static_cast<const F&>(*this),d2)); }
    Void _adjoin(const ScalarFunctionModelInterface<ValidatedTag>& f) {
        static_cast<F&>(*this).F::adjoin(dynamic_cast<const ScalarFunctionType&>(f)); }
    VectorFunctionModelInterface<ValidatedTag>* _join(const VectorFunctionModelInterface<ValidatedTag>& f) const {
        return heap_copy(join(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    VectorFunctionModelInterface<ValidatedTag>* _combine(const VectorFunctionModelInterface<ValidatedTag>& f) const {
        return heap_copy(combine(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Vector<ValidatedNumber> _unchecked_evaluate(const Vector<ValidatedNumber>& x) const {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<ValidatedTag>* _compose(const ScalarFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedTag>* _compose(const VectorFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    ScalarFunctionModelInterface<ValidatedTag>* _unchecked_compose(const ScalarFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(unchecked_compose(dynamic_cast<const ScalarFunctionType&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedTag>* _unchecked_compose(const VectorFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(unchecked_compose(dynamic_cast<const F&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedTag>* _partial_evaluate(Nat j, const ValidatedNumber& c) const {
        return heap_copy(partial_evaluate(static_cast<const F&>(*this),j,c)); }
};

template<class X> class VectorFunctionModelElement {
    VectorFunctionModel<X>* _p; Nat _i;
  public:
    operator const ScalarFunctionModel<X> () const;
    VectorFunctionModelElement(VectorFunctionModel<X>* p, Nat i) : _p(p), _i(i) { }
    VectorFunctionModelElement<X>& operator=(const ScalarFunctionModel<X>& sf) { _p->set(_i,sf); return *this; }
    VectorFunctionModelElement<X>& operator=(const VectorFunctionModelElement<X>& sf) { return this->operator=(static_cast<ScalarFunctionModel<X>const>(sf)); }
    Void clobber() { ScalarFunctionModel<ValidatedTag> sf=_p->get(_i); sf.clobber(); _p->set(_i,sf); }
    Void set_error(ErrorType e) const { ScalarFunctionModel<ValidatedTag> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
};
template<class X> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionModelElement<X>& function) {
    return os << static_cast< const ScalarFunctionModel<X> >(function);
}

//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on bounded domains.
template<> class VectorFunctionModel<ValidatedTag>
{
  public:
    clone_on_copy_ptr< VectorFunctionModelInterface<ValidatedTag> > _ptr;
  public:
    inline VectorFunctionModel() : _ptr() { }
    inline VectorFunctionModel(Nat n, const ScalarFunctionModelInterface<ValidatedTag>& sf)
        : _ptr(sf._create_vector(n)) { for(uint i=0; i!=n; ++i) { (*this)[i]=sf; } }
    inline VectorFunctionModel(VectorFunctionModelInterface<ValidatedTag>* p) : _ptr(p) { }
    inline VectorFunctionModel(const VectorFunctionModelInterface<ValidatedTag>& f) : _ptr(f._clone()) { }
    inline VectorFunctionModel(const VectorFunctionModel<ValidatedTag>& f) : _ptr(f._ptr) { }
    inline operator const VectorFunctionModelInterface<ValidatedTag>& () const { return *_ptr; }
    inline operator ValidatedVectorFunction () const { return ValidatedVectorFunction(*_ptr); }
    inline const VectorFunctionModelInterface<ValidatedTag>* raw_pointer() const { return _ptr.operator->(); }
    inline const VectorFunctionModelInterface<ValidatedTag>& reference() const { return *_ptr; }
    inline VectorFunctionModelInterface<ValidatedTag>& reference() { return *_ptr; }
    inline ScalarFunctionModel<ValidatedTag> create_zero() const { return this->_ptr->_create_zero(); }
    inline VectorFunctionModel<ValidatedTag> create_identity() const { return this->_ptr->_create_identity(); }
    inline Nat result_size() const { return this->_ptr->result_size(); }
    inline Nat argument_size() const { return this->_ptr->argument_size(); }
    inline Nat size() const { return this->_ptr->result_size(); }
    template<class XX> inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    template<class XX> inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    inline ScalarFunctionModel<ValidatedTag> const get(Nat i) const { return this->_ptr->_get(i); }
    inline Void set(Nat i, ScalarFunctionModel<ValidatedTag> const& sf) { this->_ptr->_set(i,sf); }
    inline ScalarFunctionModel<ValidatedTag> const operator[](Nat i) const { return this->get(i); }
    inline VectorFunctionModelElement<ValidatedTag> operator[](Nat i) { return VectorFunctionModelElement<ValidatedTag>(this,i); }
    inline ExactBox const domain() const { return this->_ptr->domain(); }
    inline ExactBox const codomain() const { return this->_ptr->codomain(); }
    inline UpperBox const range() const { return this->_ptr->range(); }
    inline Vector<ErrorType> const errors() const { return this->_ptr->errors(); }
    inline ErrorType const error() const { return this->_ptr->error(); }
    inline Void clobber() { this->_ptr->clobber(); }
    inline Matrix<ValidatedNumber> const jacobian(const Vector<ValidatedNumber>& x) const { return this->_ptr->jacobian(x); }

    inline Void restrict(const ExactBox& d) { this->_ptr->restrict(d); }

};


inline NormType norm(const VectorFunctionModel<ValidatedTag>& f) {
    return f._ptr->_norm(); }
inline VectorFunctionModel<ValidatedTag> embed(const ExactBox& d1, const VectorFunctionModel<ValidatedTag>& f, const ExactBox& d2) {
    return f._ptr->_embed(d1,d2); }
inline VectorFunctionModel<ValidatedTag> embed(const VectorFunctionModel<ValidatedTag>& f, const ExactBox& d) {
    return embed(ExactBox(),f,d); }
inline VectorFunctionModel<ValidatedTag> embed(const VectorFunctionModel<ValidatedTag>& f, const ExactInterval& d) {
    return embed(f,ExactBox(1,d)); }
inline VectorFunctionModel<ValidatedTag> restrict(const VectorFunctionModel<ValidatedTag>& f, const ExactBox& d) {
    VectorFunctionModelInterface<ValidatedTag>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }

inline VectorFunctionModel<ValidatedTag> operator+(const VectorFunctionModel<ValidatedTag>& f) {
    return f._ptr->_clone(); }
inline VectorFunctionModel<ValidatedTag> operator-(const VectorFunctionModel<ValidatedTag>& f) {
    VectorFunctionModel<ValidatedTag> r=f; for(uint i=0; i!=r.size(); ++i) { r[i]=-f[i]; } return r; }
inline VectorFunctionModel<ValidatedTag> operator+(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
    VectorFunctionModel<ValidatedTag> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
inline VectorFunctionModel<ValidatedTag> operator-(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
    VectorFunctionModel<ValidatedTag> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
inline VectorFunctionModel<ValidatedTag> operator+(const VectorFunctionModel<ValidatedTag>& f1, const Vector<ValidatedNumber>& c2) {
    VectorFunctionModel<ValidatedTag> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]+c2[i]; } return r; }
inline VectorFunctionModel<ValidatedTag> operator-(const VectorFunctionModel<ValidatedTag>& f1, const Vector<ValidatedNumber>& c2) {
    VectorFunctionModel<ValidatedTag> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]-c2[i]; } return r; }
inline VectorFunctionModel<ValidatedTag> operator+(const Vector<ValidatedNumber>& c1, const VectorFunctionModel<ValidatedTag>& f2);
inline VectorFunctionModel<ValidatedTag> operator-(const Vector<ValidatedNumber>& c1, const VectorFunctionModel<ValidatedTag>& f2);
inline VectorFunctionModel<ValidatedTag> operator*(const VectorFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) {
    VectorFunctionModel<ValidatedTag> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]*c2; } return r; }
inline VectorFunctionModel<ValidatedTag> operator*(const ValidatedNumber& c1, const VectorFunctionModel<ValidatedTag>& f2) {
    VectorFunctionModel<ValidatedTag> r=f2; for(uint i=0; i!=r.size(); ++i) { r[i]=c1*f2[i]; } return r; }

inline ValidatedNumber evaluate(const ScalarFunctionModel<ValidatedTag>& f, const Vector<ValidatedNumber>& x) { return f._ptr->evaluate(x); }
inline Vector<ValidatedNumber> evaluate(const VectorFunctionModel<ValidatedTag>& f, const Vector<ValidatedNumber>& x) { return f._ptr->evaluate(x); }

inline ValidatedNumber unchecked_evaluate(const ScalarFunctionModel<ValidatedTag>& f, const Vector<ValidatedNumber>& x) { return f._ptr->_unchecked_evaluate(x); }
inline Vector<ValidatedNumber> unchecked_evaluate(const VectorFunctionModel<ValidatedTag>& f, const Vector<ValidatedNumber>& x) { return f._ptr->_unchecked_evaluate(x); }

inline ScalarFunctionModel<ValidatedTag> compose(const ValidatedScalarFunction& f, const VectorFunctionModel<ValidatedTag>& g) { return g._ptr->_compose(f); }
inline ScalarFunctionModel<ValidatedTag> compose(const ScalarFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag>& g) { return g._ptr->_compose(f); }
inline VectorFunctionModel<ValidatedTag> compose(const ValidatedVectorFunction& f, const VectorFunctionModel<ValidatedTag>& g) { return g._ptr->_compose(f); }
inline VectorFunctionModel<ValidatedTag> compose(const VectorFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag>& g) { return g._ptr->_compose(f); }

inline ScalarFunctionModel<ValidatedTag> unchecked_compose(const ScalarFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag>& g) { return g._ptr->_unchecked_compose(f); }
inline VectorFunctionModel<ValidatedTag> unchecked_compose(const VectorFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag>& g) { return g._ptr->_unchecked_compose(f); }

inline ValidatedNumber unchecked_evaluate(const ValidatedScalarFunction& f, const Vector<ValidatedNumber>& x) {
    ScalarFunctionModelInterface<ValidatedTag> const* fptr = dynamic_cast<ScalarFunctionModelInterface<ValidatedTag> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(ScalarFunctionModel<ValidatedTag>(*fptr),x); } else { return evaluate(f,x); } }
inline Vector<ValidatedNumber> unchecked_evaluate(const ValidatedVectorFunction& f, const Vector<ValidatedNumber>& x) {
    VectorFunctionModelInterface<ValidatedTag> const* fptr = dynamic_cast<VectorFunctionModelInterface<ValidatedTag> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(VectorFunctionModel<ValidatedTag>(*fptr),x); } else { return evaluate(f,x); } }
inline ScalarFunctionModel<ValidatedTag> unchecked_compose(const ValidatedScalarFunction& f, const VectorFunctionModel<ValidatedTag>& g) {
    ScalarFunctionModelInterface<ValidatedTag> const* fptr = dynamic_cast<ScalarFunctionModelInterface<ValidatedTag> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(ScalarFunctionModel<ValidatedTag>(*fptr),g); } else { return compose(f,g); } }
inline VectorFunctionModel<ValidatedTag> unchecked_compose(const ValidatedVectorFunction& f, const VectorFunctionModel<ValidatedTag>& g) {
    VectorFunctionModelInterface<ValidatedTag> const* fptr = dynamic_cast<VectorFunctionModelInterface<ValidatedTag> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(VectorFunctionModel<ValidatedTag>(*fptr),g); } else { return compose(f,g); } }

//inline VectorFunctionModel<ValidatedTag> compose(const VectorFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag>& g) { return g._ptr->_compose(f); }

inline VectorFunctionModel<ValidatedTag> operator+(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunction<ValidatedTag>& f2) {
    return f1+compose(f2,f1.create_identity()); }
inline VectorFunctionModel<ValidatedTag> operator-(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunction<ValidatedTag>& f2) {
    return f1-compose(f2,f1.create_identity()); }
inline VectorFunctionModel<ValidatedTag> operator+(const VectorFunction<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
    return compose(f1,f2.create_identity())+f2; }
inline VectorFunctionModel<ValidatedTag> operator-(const VectorFunction<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
    return compose(f1,f2.create_identity())-f2; }

inline VectorFunctionModel<ValidatedTag> join(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2);
inline VectorFunctionModel<ValidatedTag> join(const ScalarFunctionModel<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2);
inline VectorFunctionModel<ValidatedTag> join(const VectorFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    VectorFunctionModel<ValidatedTag> r=f1._ptr->_clone(); r._ptr->_adjoin(f2); return r; }
inline VectorFunctionModel<ValidatedTag> join(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_join(f2); }

inline VectorFunctionModel<ValidatedTag> combine(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return VectorFunctionModel<ValidatedTag>(1,f1)._ptr->_combine(VectorFunctionModel<ValidatedTag>(1,f2)); };;
inline VectorFunctionModel<ValidatedTag> combine(const ScalarFunctionModel<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
    return VectorFunctionModel<ValidatedTag>(1,f1)._ptr->_combine(f2); };;
inline VectorFunctionModel<ValidatedTag> combine(const VectorFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_combine(VectorFunctionModel<ValidatedTag>(1,f2)); };
inline VectorFunctionModel<ValidatedTag> combine(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_combine(f2); }

inline VectorFunctionModel<ValidatedTag> intersection(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
    ARIADNE_ASSERT_MSG(f1.size()==f2.size(),"intersection(f1,f2): f1="<<f1<<", f2="<<f2<<")");
    VectorFunctionModel<ValidatedTag> r=+f1; for(uint i=0; i!=r.size(); ++i) { r[i]=intersection(f1[i],f2[i]); } return r; }

inline VectorFunctionModel<ValidatedTag> antiderivative(const VectorFunctionModel<ValidatedTag>& f, Nat j) {
    VectorFunctionModel<ValidatedTag> r(f);
    for(uint i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j); }
    return r;
}

inline VectorFunctionModel<ValidatedTag> partial_evaluate(const VectorFunctionModel<ValidatedTag>& f, Nat j, const ValidatedNumber& c) {
    return f._ptr->_partial_evaluate(j,c); }


template<class X> VectorFunctionModelElement<X>::operator const ScalarFunctionModel<X> () const {
    return _p->get(_i); }


inline VectorFunctionModel<ValidatedTag> ScalarFunctionModel<ValidatedTag>::create_identity() const { return this->_ptr->_create_identity(); }
inline ScalarFunctionModel<ValidatedTag> operator+(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedScalarFunction& f2) {
    return f1+compose(f2,f1.create_identity()); }
inline ScalarFunctionModel<ValidatedTag> operator-(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedScalarFunction& f2) {
    return f1-compose(f2,f1.create_identity()); }




// Exact output
template<class T> struct Representation { const T* pointer; Representation(const T& t) : pointer(&t) { } const T& reference() const { return *pointer; } };
template<class T> inline Representation<T> representation(const T& t) { return Representation<T>(t); }
template<class T> inline OutputStream& operator<<(OutputStream& os, const Representation<T>& obj) { obj.reference().repr(os); return os; }

template<class X> class FunctionModelFactoryInterface;

template<> class FunctionModelFactoryInterface<ValidatedTag>
{
    typedef ExactBox DomainType;
  public:
    virtual FunctionModelFactoryInterface<ValidatedTag>* clone() const = 0;
    virtual Void write(OutputStream& os) const = 0;
    inline ScalarFunctionModel<ValidatedTag> create(const ExactBox& domain, const ScalarFunctionInterface<ValidatedTag>& function) const;
    inline VectorFunctionModel<ValidatedTag> create(const ExactBox& domain, const VectorFunctionInterface<ValidatedTag>& function) const;
    inline ScalarFunctionModel<ValidatedTag> create_zero(const ExactBox& domain) const;
    inline VectorFunctionModel<ValidatedTag> create_zeros(Nat result_size, const ExactBox& domain) const;
    inline ScalarFunctionModel<ValidatedTag> create_constant(const ExactBox& domain, const ValidatedNumber& value) const;
    inline VectorFunctionModel<ValidatedTag> create_constants(const ExactBox& domain, const Vector<ValidatedNumber>& values) const;
    inline ScalarFunctionModel<ValidatedTag> create_coordinate(const ExactBox& domain, Nat index) const;
    inline ScalarFunctionModel<ValidatedTag> create_identity(const ExactInterval& domain) const;
    inline VectorFunctionModel<ValidatedTag> create_identity(const ExactBox& domain) const;
  private:
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create(const ExactBox& domain, const ScalarFunctionInterface<ValidatedTag>& function) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create(const ExactBox& domain, const VectorFunctionInterface<ValidatedTag>& function) const = 0;
};

inline OutputStream& operator<<(OutputStream& os, const FunctionModelFactoryInterface<ValidatedTag>& factory) {
    factory.write(os); return os;
}

inline ScalarFunctionModel<ValidatedTag>
FunctionModelFactoryInterface<ValidatedTag>::create(const ExactBox& domain,
                                                const ScalarFunctionInterface<ValidatedTag>& function) const
{
    return this->_create(domain,function);
}

VectorFunctionModel<ValidatedTag>
FunctionModelFactoryInterface<ValidatedTag>::create(const ExactBox& domain,
                                                const VectorFunctionInterface<ValidatedTag>& function) const
{
    return this->_create(domain,function);
}


} // namespace Ariadne

#include "function.h"

namespace Ariadne {

ScalarFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_zero(const ExactBox& domain) const {
    return this->_create(domain,EffectiveScalarFunction::zero(domain.size())); }
VectorFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_zeros(Nat result_size, const ExactBox& domain) const {
    return this->_create(domain,EffectiveVectorFunction::zeros(result_size,domain.size())); }
ScalarFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_constant(const ExactBox& domain, const ValidatedNumber& value) const {
    return ScalarFunctionModel<ValidatedTag>(this->_create(domain,EffectiveScalarFunction::zero(domain.size())))+value; };
VectorFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_constants(const ExactBox& domain, const Vector<ValidatedNumber>& values) const {
    return VectorFunctionModel<ValidatedTag>(this->_create(domain,EffectiveVectorFunction::zeros(values.size(),domain.size())))+values; };
ScalarFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_coordinate(const ExactBox& domain, Nat index) const {
    return ScalarFunctionModel<ValidatedTag>(this->_create(domain,EffectiveScalarFunction::coordinate(domain.size(),index))); };
ScalarFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_identity(const ExactInterval& domain) const {
    return this->_create(ExactBox(1,domain),EffectiveScalarFunction::coordinate(1,0)); };
VectorFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_identity(const ExactBox& domain) const {
    return this->_create(domain,EffectiveVectorFunction::identity(domain.size())); };




} // namespace Ariadne

#endif
