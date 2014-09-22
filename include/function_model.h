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

typedef ScalarFunctionModelInterface<ValidatedNumberType> ValidatedScalarFunctionModelInterface;
typedef VectorFunctionModelInterface<ValidatedNumberType> ValidatedVectorFunctionModelInterface;

typedef ScalarFunctionModel<ValidatedNumberType> ValidatedScalarFunctionModel;
typedef VectorFunctionModel<ValidatedNumberType> ValidatedVectorFunctionModel;

template<class X> class FunctionModelFactory;
typedef FunctionModelFactory<ValidatedNumberType> ValidatedFunctionModelFactory;

class ScalarTaylorFunction;
class VectorTaylorFunction;

template<> class ScalarFunctionModelInterface<ValidatedNumberType>
    : public virtual ScalarFunctionInterface<ValidatedNumberType>
{
  public:
    virtual Box const& domain() const = 0;
    virtual Interval range() const = 0;
    virtual Interval const codomain() const = 0;

    virtual CoefficientType const& value() const = 0;
    virtual CoefficientType const gradient_value(Nat i) const = 0;
    virtual ErrorType const& error() const = 0;

    virtual Void set_error(const ErrorType& e) = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _apply(OperatorCode op) const = 0;
    virtual ValidatedNumberType _unchecked_evaluate(const Vector<ValidatedNumberType>& x) const = 0;

    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _create() const = 0;
    virtual VectorFunctionModelInterface<ValidatedNumberType>* _create_vector(Nat i) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _embed(const Box& d1, const Box& d2) const = 0;
    virtual Void restrict(const Box& d) = 0;

    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _derivative(Nat j) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _antiderivative(Nat j) const = 0;

    virtual Tribool _refines(const ScalarFunctionModelInterface<ValidatedNumberType>& f) const = 0;
    virtual Tribool _disjoint(const ScalarFunctionModelInterface<ValidatedNumberType>& f) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _intersection(const ScalarFunctionModelInterface<ValidatedNumberType>& f) const = 0;

    virtual Void _iadd(const ValidatedNumberType& c) = 0;
    virtual Void _imul(const ValidatedNumberType& c) = 0;
    virtual Void _isma(const ValidatedNumberType& c, const ScalarFunctionModelInterface<ValidatedNumberType>& f) = 0;
    virtual Void _ifma(const ScalarFunctionModelInterface<ValidatedNumberType>& f1, const ScalarFunctionModelInterface<ValidatedNumberType>& f2) = 0;
};


template<class F> class ScalarFunctionModelMixin<F,ValidatedNumberType>
    : public virtual ScalarFunctionModelInterface<ValidatedNumberType>
    , public ScalarFunctionMixin<F,ValidatedNumberType>
{
  public:
    F apply(OperatorCode op) const;
  public:
    ScalarFunctionModelInterface<ValidatedNumberType>* _clone() const {
        return new F(static_cast<const F&>(*this)); }
    NormType const _norm() const {
        return norm(static_cast<const F&>(*this)); }
    ScalarFunctionModelInterface<ValidatedNumberType>* _antiderivative(Nat j) const {
        return new F(antiderivative(static_cast<const F&>(*this),j)); }
    ScalarFunctionModelInterface<ValidatedNumberType>* _apply(OperatorCode op) const {
        return new F(this->apply(op)); }
    ValidatedNumberType _unchecked_evaluate(const Vector<ValidatedNumberType>& x) const {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<ValidatedNumberType>* _embed(const Box& d1, const Box& d2) const {
        return new F(embed(d1,static_cast<const F&>(*this),d2)); }
    Tribool _refines(const ScalarFunctionModelInterface<ValidatedNumberType>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return refines(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    Tribool _disjoint(const ScalarFunctionModelInterface<ValidatedNumberType>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return disjoint(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    ScalarFunctionModelInterface<ValidatedNumberType>* _intersection(const ScalarFunctionModelInterface<ValidatedNumberType>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return new F(intersection(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Void _iadd(const ValidatedNumberType& c) {
        static_cast<F&>(*this)+=c; }
    Void _imul(const ValidatedNumberType& c) {
        static_cast<F&>(*this)*=c; }
    Void _isma(const ValidatedNumberType& c, const ScalarFunctionModelInterface<ValidatedNumberType>& f) {
        static_cast<F&>(*this)+=c*dynamic_cast<const F&>(f); }
    Void _ifma(const ScalarFunctionModelInterface<ValidatedNumberType>& f1, const ScalarFunctionModelInterface<ValidatedNumberType>& f2) {
        static_cast<F&>(*this)+=dynamic_cast<const F&>(f1)*dynamic_cast<const F&>(f2); }
};

//! \ingroup FunctionModelSubModule
//! \brief Generic scalar functions on bounded domains.
template<> class ScalarFunctionModel<ValidatedNumberType>
{
  public:
    clone_on_copy_ptr< ScalarFunctionModelInterface<ValidatedNumberType> > _ptr;
  public:
    ScalarFunctionModel() : _ptr() { }
    ScalarFunctionModel(ScalarFunctionModelInterface<ValidatedNumberType>* p) : _ptr(p) { }
    ScalarFunctionModel(const ScalarFunctionModel<ValidatedNumberType>& f) : _ptr(f._ptr) { }
    ScalarFunctionModel(const ScalarFunctionModelInterface<ValidatedNumberType>& f) : _ptr(f._clone()) { }
    ScalarFunctionModel(const ValidatedScalarFunction& f) : _ptr(dynamic_cast<ScalarFunctionModelInterface<ValidatedNumberType>*>(f.raw_pointer()->_clone())) { }
    operator ValidatedScalarFunction() const { return ValidatedScalarFunction(this->_ptr->_clone()); }
    operator ScalarFunctionModelInterface<ValidatedNumberType>& () { return *_ptr; }
    operator const ScalarFunctionModelInterface<ValidatedNumberType>& () const { return *_ptr; }
    const ScalarFunctionModelInterface<ValidatedNumberType>* raw_pointer() const { return _ptr.operator->(); }
    ScalarFunctionModelInterface<ValidatedNumberType>& reference() { return *_ptr; }
    const ScalarFunctionModelInterface<ValidatedNumberType>& reference() const { return *_ptr; }
    ScalarFunctionModel<ValidatedNumberType> create_zero() const { return this->_ptr->_create(); }
    ScalarFunctionModel<ValidatedNumberType>& operator=(const ValidatedNumberType& c);
    ScalarFunctionModel<ValidatedNumberType>& operator=(const ValidatedScalarFunction& f);
    ScalarFunctionModel<ValidatedNumberType>& operator=(const ScalarTaylorFunction& f);
    inline Nat argument_size() const { return this->_ptr->argument_size(); }
    template<class XX> inline XX operator()(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    template<class XX> inline XX evaluate(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    inline Box const domain() const { return this->_ptr->domain(); }
    inline Interval const range() const { return this->_ptr->range(); }
    inline Interval const codomain() const { return this->_ptr->codomain(); }

    inline CoefficientType value() const { return this->_ptr->value(); }
    inline CoefficientType gradient_value(Nat j) const { return this->_ptr->gradient_value(j); }
    inline ErrorType error() const { return this->_ptr->error(); }

    inline Void set_error(const ErrorType& e) { return this->_ptr->set_error(e); }
    inline Void clobber() { return this->_ptr->clobber(); }

    inline ScalarFunctionModel<ValidatedNumberType> apply(Operator op) const { return this->_ptr->_apply(op); }
    inline Void restrict(const Box& d) { this->_ptr->restrict(d); }
};

// inline ScalarFunctionModel<ValidatedNumberType>& ScalarFunctionModel<ValidatedNumberType>::operator=(const ValidatedScalarFunction& f) { (*this)=this->_ptr->_create()+f; return *this; }

inline NormType norm(const ScalarFunctionModel<ValidatedNumberType>& f) { return f._ptr->_norm(); }
inline ScalarFunctionModel<ValidatedNumberType> derivative(const ScalarFunctionModel<ValidatedNumberType>& f, Nat j) { return f._ptr->_derivative(j); }
inline ScalarFunctionModel<ValidatedNumberType> antiderivative(const ScalarFunctionModel<ValidatedNumberType>& f, Nat j) { return f._ptr->_antiderivative(j); }

inline ScalarFunctionModel<ValidatedNumberType> embed(const Box& d1, const ScalarFunctionModel<ValidatedNumberType>& f, const Box& d2) {
    return f._ptr->_embed(d1,d2); }
inline ScalarFunctionModel<ValidatedNumberType> embed(const Box& d, const ScalarFunctionModel<ValidatedNumberType>& f) {
    return embed(d,f,Box()); }
inline ScalarFunctionModel<ValidatedNumberType> embed(const ScalarFunctionModel<ValidatedNumberType>& f, const Box& d) {
    return embed(Box(),f,d); }
inline ScalarFunctionModel<ValidatedNumberType> embed(const ScalarFunctionModel<ValidatedNumberType>& f, const Interval& d) {
    return embed(f,Box(1,d)); }
inline ScalarFunctionModel<ValidatedNumberType> restrict(const ScalarFunctionModel<ValidatedNumberType>& f, const Box& d) {
    ScalarFunctionModelInterface<ValidatedNumberType>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }
inline ScalarFunctionModel<ValidatedNumberType> intersection(const ScalarFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) {
    return f1._ptr->_intersection(f2); }

inline Tribool disjoint(const ScalarFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) {
    return f1._ptr->_disjoint(f2); }
inline Tribool refines(const ScalarFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) {
    return f1._ptr->_refines(f2); }

inline ScalarFunctionModel<ValidatedNumberType>& operator+=(ScalarFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) { f1._ptr->_isma(ValidatedNumberType(+1),f2); return  f1; }
inline ScalarFunctionModel<ValidatedNumberType>& operator-=(ScalarFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) { f1._ptr->_isma(ValidatedNumberType(-1),f2); return  f1; }
inline ScalarFunctionModel<ValidatedNumberType>& operator+=(ScalarFunctionModel<ValidatedNumberType>& f1, const ValidatedNumberType& c2) { f1._ptr->_iadd(c2); return f1; }
inline ScalarFunctionModel<ValidatedNumberType>& operator-=(ScalarFunctionModel<ValidatedNumberType>& f1, const ValidatedNumberType& c2) { f1._ptr->_iadd(neg(c2)); return f1; }
inline ScalarFunctionModel<ValidatedNumberType>& operator*=(ScalarFunctionModel<ValidatedNumberType>& f1, const ValidatedNumberType& c2) { f1._ptr->_imul(c2); return f1; }
inline ScalarFunctionModel<ValidatedNumberType>& operator/=(ScalarFunctionModel<ValidatedNumberType>& f1, const ValidatedNumberType& c2) { f1._ptr->_imul(rec(c2)); return f1; }

inline ScalarFunctionModel<ValidatedNumberType> neg(const ScalarFunctionModel<ValidatedNumberType>& f) { return f.apply(Neg()); }
inline ScalarFunctionModel<ValidatedNumberType> rec(const ScalarFunctionModel<ValidatedNumberType>& f) { return f.apply(Rec()); }

inline ScalarFunctionModel<ValidatedNumberType> operator+(const ScalarFunctionModel<ValidatedNumberType>& f) {
    return f._ptr->_clone(); }
inline ScalarFunctionModel<ValidatedNumberType> operator-(const ScalarFunctionModel<ValidatedNumberType>& f) {
    ScalarFunctionModel<ValidatedNumberType> r=f; r*=-1; return r; }
inline ScalarFunctionModel<ValidatedNumberType> operator+(const ScalarFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) {
    ScalarFunctionModel<ValidatedNumberType> r=f1; r+=f2; return r; }
inline ScalarFunctionModel<ValidatedNumberType> operator-(const ScalarFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) {
    ScalarFunctionModel<ValidatedNumberType> r=f1; r-=f2; return r; }
inline ScalarFunctionModel<ValidatedNumberType> operator*(const ScalarFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) {
    ScalarFunctionModel<ValidatedNumberType> r=f1.create_zero(); r._ptr->_ifma(f1,f2); return r; }
inline ScalarFunctionModel<ValidatedNumberType> operator/(const ScalarFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) {
    return f1*rec(f2); }

inline ScalarFunctionModel<ValidatedNumberType> operator+(const ScalarFunctionModel<ValidatedNumberType>& f1, const ValidatedScalarFunction& f2) {
    ScalarFunctionModel<ValidatedNumberType> r=f1; r+=f2; return r; }

inline ScalarFunctionModel<ValidatedNumberType> operator+(const ScalarFunctionModel<ValidatedNumberType>& f1, const ValidatedNumberType& c2) { ScalarFunctionModel<ValidatedNumberType> r=f1; r+=c2; return r; }
inline ScalarFunctionModel<ValidatedNumberType> operator-(const ScalarFunctionModel<ValidatedNumberType>& f1, const ValidatedNumberType& c2) { ScalarFunctionModel<ValidatedNumberType> r=f1; r-=c2; return r; }
inline ScalarFunctionModel<ValidatedNumberType> operator*(const ScalarFunctionModel<ValidatedNumberType>& f1, const ValidatedNumberType& c2) { ScalarFunctionModel<ValidatedNumberType> r=f1; r*=c2; return r; }
inline ScalarFunctionModel<ValidatedNumberType> operator/(const ScalarFunctionModel<ValidatedNumberType>& f1, const ValidatedNumberType& c2) { ScalarFunctionModel<ValidatedNumberType> r=f1; r/=c2; return r; }
inline ScalarFunctionModel<ValidatedNumberType> operator+(const ValidatedNumberType& c1, const ScalarFunctionModel<ValidatedNumberType>& f2) { ScalarFunctionModel<ValidatedNumberType> r=f2; r+=c1; return r; }
inline ScalarFunctionModel<ValidatedNumberType> operator-(const ValidatedNumberType& c1, const ScalarFunctionModel<ValidatedNumberType>& f2) { ScalarFunctionModel<ValidatedNumberType> r=neg(f2); r+=c1; return r; }
inline ScalarFunctionModel<ValidatedNumberType> operator*(const ValidatedNumberType& c1, const ScalarFunctionModel<ValidatedNumberType>& f2) { ScalarFunctionModel<ValidatedNumberType> r=f2; r*=c1; return r; }
inline ScalarFunctionModel<ValidatedNumberType> operator/(const ValidatedNumberType& c1, const ScalarFunctionModel<ValidatedNumberType>& f2) { ScalarFunctionModel<ValidatedNumberType> r=rec(f2); r*=c1; return r; }

inline ScalarFunctionModel<ValidatedNumberType>& ScalarFunctionModel<ValidatedNumberType>::operator=(const ValidatedNumberType& c) { (*this)*=0.0; (*this)+=c; return *this; }

template<class F> F ScalarFunctionModelMixin<F,ValidatedNumberType>::apply(OperatorCode op) const {
    const F& f=static_cast<const F&>(*this);
    switch(op) {
        case NEG: return neg(f);
        case REC: return rec(f);
        default: ARIADNE_FAIL_MSG("ScalarFunctionModel<ValidatedNumberType>::apply(OperatorCode op): Operator op="<<op<<" not implemented\n");
    }
}

//inline ScalarFunctionModel<ValidatedNumberType>::ScalarFunctionModel(const ValidatedScalarFunction& f)
//    : _ptr(dynamic_cast<const ScalarFunctionModelInterface<ValidatedNumberType>&>(*f.raw_pointer())._clone()) { }




//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on bounded domains.
template<> class VectorFunctionModelInterface<ValidatedNumberType>
    : public virtual VectorFunctionInterface<ValidatedNumberType>
{
  public:
    virtual Box const& domain() const = 0;
    virtual Box const range() const = 0;
    virtual Box const codomain() const = 0;
    virtual Vector<ErrorType> const errors() const = 0;
    virtual ErrorType const error() const = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual VectorFunctionModelInterface<ValidatedNumberType>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _create_zero() const = 0;
    virtual VectorFunctionModelInterface<ValidatedNumberType>* _create_identity() const = 0;
    virtual Void _set(Nat, ScalarFunctionModelInterface<ValidatedNumberType> const&) = 0;
    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _get(Nat) const = 0;
    virtual VectorFunctionModelInterface<ValidatedNumberType>* _embed(const Box& d1, const Box& d2) const = 0;
    virtual VectorFunctionModelInterface<ValidatedNumberType>* _join(const VectorFunctionModelInterface<ValidatedNumberType>& f2) const = 0;
    virtual VectorFunctionModelInterface<ValidatedNumberType>* _combine(const VectorFunctionModelInterface<ValidatedNumberType>& f2) const = 0;
    virtual Void _adjoin(const ScalarFunctionModelInterface<ValidatedNumberType>& f2) = 0;
    virtual Vector<ValidatedNumberType> _unchecked_evaluate(const Vector<ValidatedNumberType>& x) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _compose(const ScalarFunctionInterface<ValidatedNumberType>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedNumberType>* _compose(const VectorFunctionInterface<ValidatedNumberType>& f) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _unchecked_compose(const ScalarFunctionInterface<ValidatedNumberType>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedNumberType>* _unchecked_compose(const VectorFunctionInterface<ValidatedNumberType>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedNumberType>* _partial_evaluate(Nat j, const ValidatedNumberType& c) const = 0;
    virtual Void restrict(const Box& d) = 0;
};

class ScalarTaylorFunction;
class VectorTaylorFunction;
template<class V> struct Element;
template<> struct Element<VectorTaylorFunction> { typedef ScalarTaylorFunction Type; };

template<class F> class VectorFunctionModelMixin<F,ValidatedNumberType>
    : public virtual VectorFunctionModelInterface<ValidatedNumberType>
    , public  VectorFunctionMixin<F,ValidatedNumberType>
{
    typedef typename Element<F>::Type ScalarFunctionType;
  public:
    virtual VectorFunctionModelInterface<ValidatedNumberType>* _clone() const { return new F(static_cast<const F&>(*this)); }
    virtual Void _set(Nat i, const ScalarFunctionModelInterface<ValidatedNumberType>& sf) {
        if(!dynamic_cast<const typename F::ScalarFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorFunctionModel "<<*this<<" to "<<sf<<"\n"); }
        static_cast<F&>(*this).F::set(i,dynamic_cast<const ScalarFunctionType&>(sf)); }
    NormType const _norm() const {
         return norm(static_cast<const F&>(*this)); }
    VectorFunctionModelInterface<ValidatedNumberType>* _embed(const Box& d1, const Box& d2) const {
        return heap_copy(embed(d1,static_cast<const F&>(*this),d2)); }
    Void _adjoin(const ScalarFunctionModelInterface<ValidatedNumberType>& f) {
        static_cast<F&>(*this).F::adjoin(dynamic_cast<const ScalarFunctionType&>(f)); }
    VectorFunctionModelInterface<ValidatedNumberType>* _join(const VectorFunctionModelInterface<ValidatedNumberType>& f) const {
        return heap_copy(join(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    VectorFunctionModelInterface<ValidatedNumberType>* _combine(const VectorFunctionModelInterface<ValidatedNumberType>& f) const {
        return heap_copy(combine(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Vector<ValidatedNumberType> _unchecked_evaluate(const Vector<ValidatedNumberType>& x) const {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<ValidatedNumberType>* _compose(const ScalarFunctionInterface<ValidatedNumberType>& f) const {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedNumberType>* _compose(const VectorFunctionInterface<ValidatedNumberType>& f) const {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    ScalarFunctionModelInterface<ValidatedNumberType>* _unchecked_compose(const ScalarFunctionInterface<ValidatedNumberType>& f) const {
        return heap_copy(unchecked_compose(dynamic_cast<const ScalarFunctionType&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedNumberType>* _unchecked_compose(const VectorFunctionInterface<ValidatedNumberType>& f) const {
        return heap_copy(unchecked_compose(dynamic_cast<const F&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedNumberType>* _partial_evaluate(Nat j, const ValidatedNumberType& c) const {
        return heap_copy(partial_evaluate(static_cast<const F&>(*this),j,c)); }
};

template<class X> class VectorFunctionModelElement {
    VectorFunctionModel<X>* _p; Nat _i;
  public:
    operator const ScalarFunctionModel<X> () const;
    VectorFunctionModelElement(VectorFunctionModel<X>* p, Nat i) : _p(p), _i(i) { }
    VectorFunctionModelElement<X>& operator=(const ScalarFunctionModel<X>& sf) { _p->set(_i,sf); return *this; }
    VectorFunctionModelElement<X>& operator=(const VectorFunctionModelElement<X>& sf) { return this->operator=(static_cast<ScalarFunctionModel<X>const>(sf)); }
    Void clobber() { ScalarFunctionModel<ValidatedNumberType> sf=_p->get(_i); sf.clobber(); _p->set(_i,sf); }
    Void set_error(ErrorType e) const { ScalarFunctionModel<ValidatedNumberType> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
};
template<class X> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionModelElement<X>& function) {
    return os << static_cast< const ScalarFunctionModel<X> >(function);
}

//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on bounded domains.
template<> class VectorFunctionModel<ValidatedNumberType>
{
  public:
    clone_on_copy_ptr< VectorFunctionModelInterface<ValidatedNumberType> > _ptr;
  public:
    inline VectorFunctionModel() : _ptr() { }
    inline VectorFunctionModel(Nat n, const ScalarFunctionModelInterface<ValidatedNumberType>& sf)
        : _ptr(sf._create_vector(n)) { for(uint i=0; i!=n; ++i) { (*this)[i]=sf; } }
    inline VectorFunctionModel(VectorFunctionModelInterface<ValidatedNumberType>* p) : _ptr(p) { }
    inline VectorFunctionModel(const VectorFunctionModelInterface<ValidatedNumberType>& f) : _ptr(f._clone()) { }
    inline VectorFunctionModel(const VectorFunctionModel<ValidatedNumberType>& f) : _ptr(f._ptr) { }
    inline operator const VectorFunctionModelInterface<ValidatedNumberType>& () const { return *_ptr; }
    inline operator ValidatedVectorFunction () const { return ValidatedVectorFunction(*_ptr); }
    inline const VectorFunctionModelInterface<ValidatedNumberType>* raw_pointer() const { return _ptr.operator->(); }
    inline const VectorFunctionModelInterface<ValidatedNumberType>& reference() const { return *_ptr; }
    inline VectorFunctionModelInterface<ValidatedNumberType>& reference() { return *_ptr; }
    inline ScalarFunctionModel<ValidatedNumberType> create_zero() const { return this->_ptr->_create_zero(); }
    inline VectorFunctionModel<ValidatedNumberType> create_identity() const { return this->_ptr->_create_identity(); }
    inline Nat result_size() const { return this->_ptr->result_size(); }
    inline Nat argument_size() const { return this->_ptr->argument_size(); }
    inline Nat size() const { return this->_ptr->result_size(); }
    template<class XX> inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    template<class XX> inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    inline ScalarFunctionModel<ValidatedNumberType> const get(Nat i) const { return this->_ptr->_get(i); }
    inline Void set(Nat i, ScalarFunctionModel<ValidatedNumberType> const& sf) { this->_ptr->_set(i,sf); }
    inline ScalarFunctionModel<ValidatedNumberType> const operator[](Nat i) const { return this->get(i); }
    inline VectorFunctionModelElement<ValidatedNumberType> operator[](Nat i) { return VectorFunctionModelElement<ValidatedNumberType>(this,i); }
    inline Box const domain() const { return this->_ptr->domain(); }
    inline Box const range() const { return this->_ptr->range(); }
    inline Box const codomain() const { return this->_ptr->codomain(); }
    inline Vector<ErrorType> const errors() const { return this->_ptr->errors(); }
    inline ErrorType const error() const { return this->_ptr->error(); }
    inline Void clobber() { this->_ptr->clobber(); }
    inline Matrix<ValidatedNumberType> const jacobian(const Vector<ValidatedNumberType>& x) const { return this->_ptr->jacobian(x); }

    inline Void restrict(const Box& d) { this->_ptr->restrict(d); }

};

inline NormType norm(const VectorFunctionModel<ValidatedNumberType>& f) {
    return f._ptr->_norm(); }
inline VectorFunctionModel<ValidatedNumberType> embed(const Box& d1, const VectorFunctionModel<ValidatedNumberType>& f, const Box& d2) {
    return f._ptr->_embed(d1,d2); }
inline VectorFunctionModel<ValidatedNumberType> embed(const VectorFunctionModel<ValidatedNumberType>& f, const Box& d) {
    return embed(Box(),f,d); }
inline VectorFunctionModel<ValidatedNumberType> embed(const VectorFunctionModel<ValidatedNumberType>& f, const Interval& d) {
    return embed(f,Box(1,d)); }
inline VectorFunctionModel<ValidatedNumberType> restrict(const VectorFunctionModel<ValidatedNumberType>& f, const Box& d) {
    VectorFunctionModelInterface<ValidatedNumberType>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }

inline VectorFunctionModel<ValidatedNumberType> operator+(const VectorFunctionModel<ValidatedNumberType>& f) {
    return f._ptr->_clone(); }
inline VectorFunctionModel<ValidatedNumberType> operator-(const VectorFunctionModel<ValidatedNumberType>& f) {
    VectorFunctionModel<ValidatedNumberType> r=f; for(uint i=0; i!=r.size(); ++i) { r[i]=-f[i]; } return r; }
inline VectorFunctionModel<ValidatedNumberType> operator+(const VectorFunctionModel<ValidatedNumberType>& f1, const VectorFunctionModel<ValidatedNumberType>& f2) {
    VectorFunctionModel<ValidatedNumberType> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
inline VectorFunctionModel<ValidatedNumberType> operator-(const VectorFunctionModel<ValidatedNumberType>& f1, const VectorFunctionModel<ValidatedNumberType>& f2) {
    VectorFunctionModel<ValidatedNumberType> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
inline VectorFunctionModel<ValidatedNumberType> operator+(const VectorFunctionModel<ValidatedNumberType>& f1, const Vector<ValidatedNumberType>& c2) {
    VectorFunctionModel<ValidatedNumberType> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]+c2[i]; } return r; }
inline VectorFunctionModel<ValidatedNumberType> operator-(const VectorFunctionModel<ValidatedNumberType>& f1, const Vector<ValidatedNumberType>& c2) {
    VectorFunctionModel<ValidatedNumberType> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]-c2[i]; } return r; }
inline VectorFunctionModel<ValidatedNumberType> operator+(const Vector<ValidatedNumberType>& c1, const VectorFunctionModel<ValidatedNumberType>& f2);
inline VectorFunctionModel<ValidatedNumberType> operator-(const Vector<ValidatedNumberType>& c1, const VectorFunctionModel<ValidatedNumberType>& f2);
inline VectorFunctionModel<ValidatedNumberType> operator*(const VectorFunctionModel<ValidatedNumberType>& f1, const ValidatedNumberType& c2) {
    VectorFunctionModel<ValidatedNumberType> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]*c2; } return r; }
inline VectorFunctionModel<ValidatedNumberType> operator*(const ValidatedNumberType& c1, const VectorFunctionModel<ValidatedNumberType>& f2) {
    VectorFunctionModel<ValidatedNumberType> r=f2; for(uint i=0; i!=r.size(); ++i) { r[i]=c1*f2[i]; } return r; }

inline ValidatedNumberType evaluate(const ScalarFunctionModel<ValidatedNumberType>& f, const Vector<ValidatedNumberType>& x) { return f._ptr->evaluate(x); }
inline Vector<ValidatedNumberType> evaluate(const VectorFunctionModel<ValidatedNumberType>& f, const Vector<ValidatedNumberType>& x) { return f._ptr->evaluate(x); }

inline ValidatedNumberType unchecked_evaluate(const ScalarFunctionModel<ValidatedNumberType>& f, const Vector<ValidatedNumberType>& x) { return f._ptr->_unchecked_evaluate(x); }
inline Vector<ValidatedNumberType> unchecked_evaluate(const VectorFunctionModel<ValidatedNumberType>& f, const Vector<ValidatedNumberType>& x) { return f._ptr->_unchecked_evaluate(x); }

inline ScalarFunctionModel<ValidatedNumberType> compose(const ValidatedScalarFunction& f, const VectorFunctionModel<ValidatedNumberType>& g) { return g._ptr->_compose(f); }
inline ScalarFunctionModel<ValidatedNumberType> compose(const ScalarFunctionModel<ValidatedNumberType>& f, const VectorFunctionModel<ValidatedNumberType>& g) { return g._ptr->_compose(f); }
inline VectorFunctionModel<ValidatedNumberType> compose(const ValidatedVectorFunction& f, const VectorFunctionModel<ValidatedNumberType>& g) { return g._ptr->_compose(f); }
inline VectorFunctionModel<ValidatedNumberType> compose(const VectorFunctionModel<ValidatedNumberType>& f, const VectorFunctionModel<ValidatedNumberType>& g) { return g._ptr->_compose(f); }

inline ScalarFunctionModel<ValidatedNumberType> unchecked_compose(const ScalarFunctionModel<ValidatedNumberType>& f, const VectorFunctionModel<ValidatedNumberType>& g) { return g._ptr->_unchecked_compose(f); }
inline VectorFunctionModel<ValidatedNumberType> unchecked_compose(const VectorFunctionModel<ValidatedNumberType>& f, const VectorFunctionModel<ValidatedNumberType>& g) { return g._ptr->_unchecked_compose(f); }

inline ValidatedNumberType unchecked_evaluate(const ValidatedScalarFunction& f, const Vector<ValidatedNumberType>& x) {
    ScalarFunctionModelInterface<ValidatedNumberType> const* fptr = dynamic_cast<ScalarFunctionModelInterface<ValidatedNumberType> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(ScalarFunctionModel<ValidatedNumberType>(*fptr),x); } else { return evaluate(f,x); } }
inline Vector<ValidatedNumberType> unchecked_evaluate(const ValidatedVectorFunction& f, const Vector<ValidatedNumberType>& x) {
    VectorFunctionModelInterface<ValidatedNumberType> const* fptr = dynamic_cast<VectorFunctionModelInterface<ValidatedNumberType> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(VectorFunctionModel<ValidatedNumberType>(*fptr),x); } else { return evaluate(f,x); } }
inline ScalarFunctionModel<ValidatedNumberType> unchecked_compose(const ValidatedScalarFunction& f, const VectorFunctionModel<ValidatedNumberType>& g) {
    ScalarFunctionModelInterface<ValidatedNumberType> const* fptr = dynamic_cast<ScalarFunctionModelInterface<ValidatedNumberType> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(ScalarFunctionModel<ValidatedNumberType>(*fptr),g); } else { return compose(f,g); } }
inline VectorFunctionModel<ValidatedNumberType> unchecked_compose(const ValidatedVectorFunction& f, const VectorFunctionModel<ValidatedNumberType>& g) {
    VectorFunctionModelInterface<ValidatedNumberType> const* fptr = dynamic_cast<VectorFunctionModelInterface<ValidatedNumberType> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(VectorFunctionModel<ValidatedNumberType>(*fptr),g); } else { return compose(f,g); } }

//inline VectorFunctionModel<ValidatedNumberType> compose(const VectorFunctionModel<ValidatedNumberType>& f, const VectorFunctionModel<ValidatedNumberType>& g) { return g._ptr->_compose(f); }

inline VectorFunctionModel<ValidatedNumberType> operator-(const VectorFunctionModel<ValidatedNumberType>& f1, const VectorFunctionInterface<ValidatedNumberType>& f2) {
    return f1-compose(f2,f1.create_identity()); }

inline VectorFunctionModel<ValidatedNumberType> join(const ScalarFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2);
inline VectorFunctionModel<ValidatedNumberType> join(const ScalarFunctionModel<ValidatedNumberType>& f1, const VectorFunctionModel<ValidatedNumberType>& f2);
inline VectorFunctionModel<ValidatedNumberType> join(const VectorFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) {
    VectorFunctionModel<ValidatedNumberType> r=f1._ptr->_clone(); r._ptr->_adjoin(f2); return r; }
inline VectorFunctionModel<ValidatedNumberType> join(const VectorFunctionModel<ValidatedNumberType>& f1, const VectorFunctionModel<ValidatedNumberType>& f2) {
    return f1._ptr->_join(f2); }

inline VectorFunctionModel<ValidatedNumberType> combine(const ScalarFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) {
    return VectorFunctionModel<ValidatedNumberType>(1,f1)._ptr->_combine(VectorFunctionModel<ValidatedNumberType>(1,f2)); };;
inline VectorFunctionModel<ValidatedNumberType> combine(const ScalarFunctionModel<ValidatedNumberType>& f1, const VectorFunctionModel<ValidatedNumberType>& f2) {
    return VectorFunctionModel<ValidatedNumberType>(1,f1)._ptr->_combine(f2); };;
inline VectorFunctionModel<ValidatedNumberType> combine(const VectorFunctionModel<ValidatedNumberType>& f1, const ScalarFunctionModel<ValidatedNumberType>& f2) {
    return f1._ptr->_combine(VectorFunctionModel<ValidatedNumberType>(1,f2)); };
inline VectorFunctionModel<ValidatedNumberType> combine(const VectorFunctionModel<ValidatedNumberType>& f1, const VectorFunctionModel<ValidatedNumberType>& f2) {
    return f1._ptr->_combine(f2); }

inline VectorFunctionModel<ValidatedNumberType> intersection(const VectorFunctionModel<ValidatedNumberType>& f1, const VectorFunctionModel<ValidatedNumberType>& f2) {
    ARIADNE_ASSERT_MSG(f1.size()==f2.size(),"intersection(f1,f2): f1="<<f1<<", f2="<<f2<<")");
    VectorFunctionModel<ValidatedNumberType> r=+f1; for(uint i=0; i!=r.size(); ++i) { r[i]=intersection(f1[i],f2[i]); } return r; }

inline VectorFunctionModel<ValidatedNumberType> antiderivative(const VectorFunctionModel<ValidatedNumberType>& f, Nat j) {
    VectorFunctionModel<ValidatedNumberType> r(f);
    for(uint i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j); }
    return r;
}

inline VectorFunctionModel<ValidatedNumberType> partial_evaluate(const VectorFunctionModel<ValidatedNumberType>& f, Nat j, const ValidatedNumberType& c) {
    return f._ptr->_partial_evaluate(j,c); }


template<class X> VectorFunctionModelElement<X>::operator const ScalarFunctionModel<X> () const {
    return _p->get(_i); }




// Exact output
template<class T> struct Representation { const T* pointer; Representation(const T& t) : pointer(&t) { } const T& reference() const { return *pointer; } };
template<class T> inline Representation<T> representation(const T& t) { return Representation<T>(t); }
template<class T> inline OutputStream& operator<<(OutputStream& os, const Representation<T>& obj) { obj.reference().repr(os); return os; }

template<class X> class FunctionModelFactoryInterface;

template<> class FunctionModelFactoryInterface<ValidatedNumberType>
{
    typedef Box DomainType;
  public:
    virtual FunctionModelFactoryInterface<ValidatedNumberType>* clone() const = 0;
    virtual Void write(OutputStream& os) const = 0;
    inline ScalarFunctionModel<ValidatedNumberType> create(const Box& domain, const ScalarFunctionInterface<ValidatedNumberType>& function) const;
    inline VectorFunctionModel<ValidatedNumberType> create(const Box& domain, const VectorFunctionInterface<ValidatedNumberType>& function) const;
    inline ScalarFunctionModel<ValidatedNumberType> create_zero(const Box& domain) const;
    inline VectorFunctionModel<ValidatedNumberType> create_zeros(Nat result_size, const Box& domain) const;
    inline ScalarFunctionModel<ValidatedNumberType> create_constant(const Box& domain, const ValidatedNumberType& value) const;
    inline VectorFunctionModel<ValidatedNumberType> create_constants(const Box& domain, const Vector<ValidatedNumberType>& values) const;
    inline ScalarFunctionModel<ValidatedNumberType> create_coordinate(const Box& domain, Nat index) const;
    inline ScalarFunctionModel<ValidatedNumberType> create_identity(const Interval& domain) const;
    inline VectorFunctionModel<ValidatedNumberType> create_identity(const Box& domain) const;
  private:
    virtual ScalarFunctionModelInterface<ValidatedNumberType>* _create(const Box& domain, const ScalarFunctionInterface<ValidatedNumberType>& function) const = 0;
    virtual VectorFunctionModelInterface<ValidatedNumberType>* _create(const Box& domain, const VectorFunctionInterface<ValidatedNumberType>& function) const = 0;
};

inline OutputStream& operator<<(OutputStream& os, const FunctionModelFactoryInterface<ValidatedNumberType>& factory) {
    factory.write(os); return os;
}

inline ScalarFunctionModel<ValidatedNumberType>
FunctionModelFactoryInterface<ValidatedNumberType>::create(const Box& domain,
                                                const ScalarFunctionInterface<ValidatedNumberType>& function) const
{
    return this->_create(domain,function);
}

VectorFunctionModel<ValidatedNumberType>
FunctionModelFactoryInterface<ValidatedNumberType>::create(const Box& domain,
                                                const VectorFunctionInterface<ValidatedNumberType>& function) const
{
    return this->_create(domain,function);
}


} // namespace Ariadne

#include "function.h"

namespace Ariadne {

ScalarFunctionModel<ValidatedNumberType> FunctionModelFactoryInterface<ValidatedNumberType>::create_zero(const Box& domain) const {
    return this->_create(domain,EffectiveScalarFunction::zero(domain.size())); }
VectorFunctionModel<ValidatedNumberType> FunctionModelFactoryInterface<ValidatedNumberType>::create_zeros(Nat result_size, const Box& domain) const {
    return this->_create(domain,EffectiveVectorFunction::zeros(result_size,domain.size())); }
ScalarFunctionModel<ValidatedNumberType> FunctionModelFactoryInterface<ValidatedNumberType>::create_constant(const Box& domain, const ValidatedNumberType& value) const {
    return ScalarFunctionModel<ValidatedNumberType>(this->_create(domain,EffectiveScalarFunction::zero(domain.size())))+value; };
VectorFunctionModel<ValidatedNumberType> FunctionModelFactoryInterface<ValidatedNumberType>::create_constants(const Box& domain, const Vector<ValidatedNumberType>& values) const {
    return VectorFunctionModel<ValidatedNumberType>(this->_create(domain,EffectiveVectorFunction::zeros(values.size(),domain.size())))+values; };
ScalarFunctionModel<ValidatedNumberType> FunctionModelFactoryInterface<ValidatedNumberType>::create_coordinate(const Box& domain, Nat index) const {
    return ScalarFunctionModel<ValidatedNumberType>(this->_create(domain,EffectiveScalarFunction::coordinate(domain.size(),index))); };
ScalarFunctionModel<ValidatedNumberType> FunctionModelFactoryInterface<ValidatedNumberType>::create_identity(const Interval& domain) const {
    return this->_create(Box(1,domain),EffectiveScalarFunction::coordinate(1,0)); };
VectorFunctionModel<ValidatedNumberType> FunctionModelFactoryInterface<ValidatedNumberType>::create_identity(const Box& domain) const {
    return this->_create(domain,EffectiveVectorFunction::identity(domain.size())); };




} // namespace Ariadne

#endif
