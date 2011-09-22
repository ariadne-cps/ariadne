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

namespace Ariadne {

class Float;
class Interval;
class Real;

template<class X> class ScalarFunctionInterface;
typedef ScalarFunctionInterface<Interval> IntervalScalarFunctionInterface;
template<class X> class VectorFunctionInterface;
typedef VectorFunctionInterface<Interval> IntervalVectorFunctionInterface;

template<class X> class ScalarFunctionModelInterface;
template<class F,class X> class ScalarFunctionModelMixin;
template<class X> class ScalarFunctionModel;

template<class X> class VectorFunctionModelInterface;
template<class F,class X> class VectorFunctionModelMixin;
template<class X> class VectorFunctionModel;

typedef ScalarFunctionModelInterface<Interval> IntervalScalarFunctionModelInterface;
typedef VectorFunctionModelInterface<Interval> IntervalVectorFunctionModelInterface;

typedef ScalarFunctionModel<Interval> IntervalScalarFunctionModel;
typedef VectorFunctionModel<Interval> IntervalVectorFunctionModel;

template<class X> class FunctionModelFactory;
typedef FunctionModelFactory<Interval> IntervalFunctionModelFactory;

class ScalarTaylorFunction;

template<> class ScalarFunctionModelInterface<Interval>
    : public virtual ScalarFunctionInterface<Interval>
{
  public:
    virtual IntervalVector const& domain() const = 0;
    virtual Interval range() const = 0;
    virtual Interval const codomain() const = 0;

    virtual Float const& value() const = 0;
    virtual Float const gradient_value(Nat i) const = 0;
    virtual Float const& error() const = 0;

    virtual Void set_error(const Float& e) = 0;
    virtual Void clobber() = 0;

    virtual Float const _norm() const = 0;

    virtual ScalarFunctionModelInterface<Interval>* _apply(OperatorCode op) const = 0;
    virtual Interval _unchecked_evaluate(const Vector<Interval>& x) const = 0;

    virtual ScalarFunctionModelInterface<Interval>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<Interval>* _create() const = 0;
    virtual VectorFunctionModelInterface<Interval>* _create_vector(Nat i) const = 0;
    virtual ScalarFunctionModelInterface<Interval>* _embed(const IntervalVector& d1, const IntervalVector& d2) const = 0;
    virtual Void restrict(const IntervalVector& d) = 0;

    virtual ScalarFunctionModelInterface<Interval>* _derivative(Nat j) const = 0;
    virtual ScalarFunctionModelInterface<Interval>* _antiderivative(Nat j) const = 0;

    virtual Tribool _refines(const ScalarFunctionModelInterface<Interval>& f) const = 0;
    virtual Tribool _disjoint(const ScalarFunctionModelInterface<Interval>& f) const = 0;
    virtual ScalarFunctionModelInterface<Interval>* _intersection(const ScalarFunctionModelInterface<Interval>& f) const = 0;

    virtual Void _iadd(const Interval& c) = 0;
    virtual Void _imul(const Interval& c) = 0;
    virtual Void _isma(const Interval& c, const ScalarFunctionModelInterface<Interval>& f) = 0;
    virtual Void _ifma(const ScalarFunctionModelInterface<Interval>& f1, const ScalarFunctionModelInterface<Interval>& f2) = 0;
};


template<class F> class ScalarFunctionModelMixin<F,Interval>
    : public virtual ScalarFunctionModelInterface<Interval>
    , public ScalarFunctionMixin<F,Interval>
{
  public:
    F apply(OperatorCode op) const;
  public:
    ScalarFunctionModelInterface<Interval>* _clone() const {
        return new F(static_cast<const F&>(*this)); }
    Float const _norm() const {
        return norm(static_cast<const F&>(*this)); }
    ScalarFunctionModelInterface<Interval>* _antiderivative(Nat j) const {
        return new F(antiderivative(static_cast<const F&>(*this),j)); }
    ScalarFunctionModelInterface<Interval>* _apply(OperatorCode op) const {
        return new F(this->apply(op)); }
    Interval _unchecked_evaluate(const Vector<Interval>& x) const {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<Interval>* _embed(const IntervalVector& d1, const IntervalVector& d2) const {
        return new F(embed(d1,static_cast<const F&>(*this),d2)); }
    Tribool _refines(const ScalarFunctionModelInterface<Interval>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return refines(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    Tribool _disjoint(const ScalarFunctionModelInterface<Interval>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return disjoint(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    ScalarFunctionModelInterface<Interval>* _intersection(const ScalarFunctionModelInterface<Interval>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return new F(intersection(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Void _iadd(const Interval& c) {
        static_cast<F&>(*this)+=c; }
    Void _imul(const Interval& c) {
        static_cast<F&>(*this)*=c; }
    Void _isma(const Interval& c, const ScalarFunctionModelInterface<Interval>& f) {
        static_cast<F&>(*this)+=c*dynamic_cast<const F&>(f); }
    Void _ifma(const ScalarFunctionModelInterface<Interval>& f1, const ScalarFunctionModelInterface<Interval>& f2) {
        static_cast<F&>(*this)+=dynamic_cast<const F&>(f1)*dynamic_cast<const F&>(f2); }
};

template<> class ScalarFunctionModel<Interval>
{
  public:
    clone_on_copy_ptr< ScalarFunctionModelInterface<Interval> > _ptr;
  public:
    ScalarFunctionModel() : _ptr() { }
    ScalarFunctionModel(ScalarFunctionModelInterface<Interval>* p) : _ptr(p) { }
    ScalarFunctionModel(const ScalarFunctionModel<Interval>& f) : _ptr(f._ptr) { }
    ScalarFunctionModel(const ScalarFunctionModelInterface<Interval>& f) : _ptr(f._clone()) { }
    ScalarFunctionModel(const ScalarFunction<Interval>& f) : _ptr(dynamic_cast<ScalarFunctionModelInterface<Interval>*>(f.raw_pointer()->_clone())) { }
    operator ScalarFunction<Interval>() const { return ScalarFunction<Interval>(this->_ptr->_clone()); }
    operator ScalarFunctionModelInterface<Interval>& () { return *_ptr; }
    operator const ScalarFunctionModelInterface<Interval>& () const { return *_ptr; }
    const ScalarFunctionModelInterface<Interval>* raw_pointer() const { return _ptr.operator->(); }
    ScalarFunctionModelInterface<Interval>& reference() { return *_ptr; }
    const ScalarFunctionModelInterface<Interval>& reference() const { return *_ptr; }
    ScalarFunctionModel<Interval> create_zero() const { return this->_ptr->_create(); }
    ScalarFunctionModel<Interval>& operator=(const Float& c) { return this->operator=(c); }
    ScalarFunctionModel<Interval>& operator=(const Interval& c);
    ScalarFunctionModel<Interval>& operator=(const ScalarFunction<Interval>& f) { (*this)=this->_ptr->_create()+f; return *this; }
    ScalarFunctionModel<Interval>& operator=(const ScalarTaylorFunction& f);
    inline Nat argument_size() const { return this->_ptr->argument_size(); }
    template<class XX> inline XX operator()(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    template<class XX> inline XX evaluate(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    inline IntervalVector const& domain() const { return this->_ptr->domain(); }
    inline Interval const range() const { return this->_ptr->range(); }
    inline Interval const codomain() const { return this->_ptr->codomain(); }

    inline Float value() const { return this->_ptr->value(); }
    inline Float gradient_value(Nat j) const { return this->_ptr->gradient_value(j); }
    inline Float error() const { return this->_ptr->error(); }

    inline Void set_error(const Float& e) { return this->_ptr->set_error(e); }
    inline Void clobber() { return this->_ptr->clobber(); }

    inline ScalarFunctionModel<Interval> apply(Operator op) const { return this->_ptr->_apply(op); }
    inline Void restrict(const IntervalVector& d) { this->_ptr->restrict(d); }
};

inline Float norm(const ScalarFunctionModel<Interval>& f) { return f._ptr->_norm(); }
inline ScalarFunctionModel<Interval> derivative(const ScalarFunctionModel<Interval>& f, Nat j) { return f._ptr->_derivative(j); }
inline ScalarFunctionModel<Interval> antiderivative(const ScalarFunctionModel<Interval>& f, Nat j) { return f._ptr->_antiderivative(j); }

inline ScalarFunctionModel<Interval> embed(const IntervalVector& d1, const ScalarFunctionModel<Interval>& f, const IntervalVector& d2) {
    return f._ptr->_embed(d1,d2); }
inline ScalarFunctionModel<Interval> embed(const IntervalVector& d, const ScalarFunctionModel<Interval>& f) {
    return embed(d,f,IntervalVector()); }
inline ScalarFunctionModel<Interval> embed(const ScalarFunctionModel<Interval>& f, const IntervalVector& d) {
    return embed(IntervalVector(),f,d); }
inline ScalarFunctionModel<Interval> embed(const ScalarFunctionModel<Interval>& f, const Interval& d) {
    return embed(f,IntervalVector(1,d)); }
inline ScalarFunctionModel<Interval> restrict(const ScalarFunctionModel<Interval>& f, const IntervalVector& d) {
    ScalarFunctionModelInterface<Interval>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }
inline ScalarFunctionModel<Interval> intersection(const ScalarFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) {
    return f1._ptr->_intersection(f2); }

inline Tribool disjoint(const ScalarFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) {
    return f1._ptr->_disjoint(f2); }
inline Tribool refines(const ScalarFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) {
    return f1._ptr->_refines(f2); }

inline ScalarFunctionModel<Interval>& operator+=(ScalarFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) { f1._ptr->_isma(Interval(+1),f2); return  f1; }
inline ScalarFunctionModel<Interval>& operator-=(ScalarFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) { f1._ptr->_isma(Interval(-1),f2); return  f1; }
inline ScalarFunctionModel<Interval>& operator+=(ScalarFunctionModel<Interval>& f1, const Interval& c2) { f1._ptr->_iadd(c2); return f1; }
inline ScalarFunctionModel<Interval>& operator-=(ScalarFunctionModel<Interval>& f1, const Interval& c2) { f1._ptr->_iadd(neg(c2)); return f1; }
inline ScalarFunctionModel<Interval>& operator*=(ScalarFunctionModel<Interval>& f1, const Interval& c2) { f1._ptr->_imul(c2); return f1; }
inline ScalarFunctionModel<Interval>& operator/=(ScalarFunctionModel<Interval>& f1, const Interval& c2) { f1._ptr->_imul(rec(c2)); return f1; }

inline ScalarFunctionModel<Interval> neg(const ScalarFunctionModel<Interval>& f) { return f.apply(Neg()); }
inline ScalarFunctionModel<Interval> rec(const ScalarFunctionModel<Interval>& f) { return f.apply(Rec()); }

inline ScalarFunctionModel<Interval> operator+(const ScalarFunctionModel<Interval>& f) {
    return f._ptr->_clone(); }
inline ScalarFunctionModel<Interval> operator-(const ScalarFunctionModel<Interval>& f) {
    ScalarFunctionModel<Interval> r=f; r*=-1; return r; }
inline ScalarFunctionModel<Interval> operator+(const ScalarFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) {
    ScalarFunctionModel<Interval> r=f1; r+=f2; return r; }
inline ScalarFunctionModel<Interval> operator-(const ScalarFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) {
    ScalarFunctionModel<Interval> r=f1; r-=f2; return r; }
inline ScalarFunctionModel<Interval> operator*(const ScalarFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) {
    ScalarFunctionModel<Interval> r=f1.create_zero(); r._ptr->_ifma(f1,f2); return r; }
inline ScalarFunctionModel<Interval> operator/(const ScalarFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) {
    return f1*rec(f2); }

inline ScalarFunctionModel<Interval> operator+(const ScalarFunctionModel<Interval>& f1, const ScalarFunction<Interval>& f2) {
    ScalarFunctionModel<Interval> r=f1; r+=f2; return r; }

inline ScalarFunctionModel<Interval> operator+(const ScalarFunctionModel<Interval>& f1, const Interval& c2) { ScalarFunctionModel<Interval> r=f1; r+=c2; return r; }
inline ScalarFunctionModel<Interval> operator-(const ScalarFunctionModel<Interval>& f1, const Interval& c2) { ScalarFunctionModel<Interval> r=f1; r-=c2; return r; }
inline ScalarFunctionModel<Interval> operator*(const ScalarFunctionModel<Interval>& f1, const Interval& c2) { ScalarFunctionModel<Interval> r=f1; r*=c2; return r; }
inline ScalarFunctionModel<Interval> operator/(const ScalarFunctionModel<Interval>& f1, const Interval& c2) { ScalarFunctionModel<Interval> r=f1; r/=c2; return r; }
inline ScalarFunctionModel<Interval> operator+(const Interval& c1, const ScalarFunctionModel<Interval>& f2) { ScalarFunctionModel<Interval> r=f2; r+=c1; return r; }
inline ScalarFunctionModel<Interval> operator-(const Interval& c1, const ScalarFunctionModel<Interval>& f2) { ScalarFunctionModel<Interval> r=neg(f2); r+=c1; return r; }
inline ScalarFunctionModel<Interval> operator*(const Interval& c1, const ScalarFunctionModel<Interval>& f2) { ScalarFunctionModel<Interval> r=f2; r*=c1; return r; }
inline ScalarFunctionModel<Interval> operator/(const Interval& c1, const ScalarFunctionModel<Interval>& f2) { ScalarFunctionModel<Interval> r=rec(f2); r*=c1; return r; }

inline ScalarFunctionModel<Interval>& ScalarFunctionModel<Interval>::operator=(const Interval& c) { (*this)*=0.0; (*this)+=c; return *this; }

template<class F> F ScalarFunctionModelMixin<F,Interval>::apply(OperatorCode op) const {
    const F& f=static_cast<const F&>(*this);
    switch(op) {
        case NEG: return neg(f);
        case REC: return rec(f);
        default: ARIADNE_FAIL_MSG("ScalarFunctionModel<Interval>::apply(OperatorCode op): Operator op="<<op<<" not implemented\n");
    }
}

//inline ScalarFunctionModel<Interval>::ScalarFunctionModel(const ScalarFunction<Interval>& f)
//    : _ptr(dynamic_cast<const ScalarFunctionModelInterface<Interval>&>(*f.raw_pointer())._clone()) { }




template<> class VectorFunctionModelInterface<Interval>
    : public virtual VectorFunctionInterface<Interval>
{
  public:
    virtual IntervalVector const& domain() const = 0;
    virtual IntervalVector const range() const = 0;
    virtual IntervalVector const codomain() const = 0;
    virtual Vector<Float> const errors() const = 0;
    virtual Float const error() const = 0;

    virtual Float const _norm() const = 0;

    virtual VectorFunctionModelInterface<Interval>* _clone() const = 0;
    virtual VectorFunctionModelInterface<Interval>* _create_identity() const = 0;
    virtual Void _set(Nat, ScalarFunctionModelInterface<Interval> const&) = 0;
    virtual ScalarFunctionModelInterface<Interval>* _get(Nat) const = 0;
    virtual VectorFunctionModelInterface<Interval>* _embed(const IntervalVector& d1, const IntervalVector& d2) const = 0;
    virtual VectorFunctionModelInterface<Interval>* _join(const VectorFunctionModelInterface<Interval>& f2) const = 0;
    virtual VectorFunctionModelInterface<Interval>* _combine(const VectorFunctionModelInterface<Interval>& f2) const = 0;
    virtual Void _adjoin(const ScalarFunctionModelInterface<Interval>& f2) = 0;
    virtual Vector<Interval> _unchecked_evaluate(const Vector<Interval>& x) const = 0;
    virtual ScalarFunctionModelInterface<Interval>* _compose(const ScalarFunctionInterface<Interval>& f) const = 0;
    virtual VectorFunctionModelInterface<Interval>* _compose(const VectorFunctionInterface<Interval>& f) const = 0;
    virtual ScalarFunctionModelInterface<Interval>* _unchecked_compose(const ScalarFunctionInterface<Interval>& f) const = 0;
    virtual VectorFunctionModelInterface<Interval>* _unchecked_compose(const VectorFunctionInterface<Interval>& f) const = 0;
    virtual VectorFunctionModelInterface<Interval>* _partial_evaluate(Nat j, const Interval& c) const = 0;
    virtual Void restrict(const IntervalVector& d) = 0;
};

class ScalarTaylorFunction;
class VectorTaylorFunction;
template<class V> struct Element;
template<> struct Element<VectorTaylorFunction> { typedef ScalarTaylorFunction Type; };

template<class F> class VectorFunctionModelMixin<F,Interval>
    : public virtual VectorFunctionModelInterface<Interval>
    , public  VectorFunctionMixin<F,Interval>
{
    typedef typename Element<F>::Type ScalarFunctionType;
  public:
    virtual VectorFunctionModelInterface<Interval>* _clone() const { return new F(static_cast<const F&>(*this)); }
    virtual Void _set(Nat i, const ScalarFunctionModelInterface<Interval>& sf) {
        if(!dynamic_cast<const typename F::ScalarFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorFunctionModel "<<*this<<" to "<<sf<<"\n"); }
        static_cast<F&>(*this).F::set(i,dynamic_cast<const ScalarFunctionType&>(sf)); }
    Float const _norm() const {
         return norm(static_cast<const F&>(*this)); }
    VectorFunctionModelInterface<Interval>* _embed(const IntervalVector& d1, const IntervalVector& d2) const {
        return heap_copy(embed(d1,static_cast<const F&>(*this),d2)); }
    Void _adjoin(const ScalarFunctionModelInterface<Interval>& f) {
        static_cast<F&>(*this).F::adjoin(dynamic_cast<const ScalarFunctionType&>(f)); }
    VectorFunctionModelInterface<Interval>* _join(const VectorFunctionModelInterface<Interval>& f) const {
        return heap_copy(join(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    VectorFunctionModelInterface<Interval>* _combine(const VectorFunctionModelInterface<Interval>& f) const {
        return heap_copy(combine(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Vector<Interval> _unchecked_evaluate(const Vector<Interval>& x) const {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<Interval>* _compose(const ScalarFunctionInterface<Interval>& f) const {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<Interval>* _compose(const VectorFunctionInterface<Interval>& f) const {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    ScalarFunctionModelInterface<Interval>* _unchecked_compose(const ScalarFunctionInterface<Interval>& f) const {
        return heap_copy(unchecked_compose(dynamic_cast<const ScalarFunctionType&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<Interval>* _unchecked_compose(const VectorFunctionInterface<Interval>& f) const {
        return heap_copy(unchecked_compose(dynamic_cast<const F&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<Interval>* _partial_evaluate(Nat j, const Interval& c) const {
        return heap_copy(partial_evaluate(static_cast<const F&>(*this),j,c)); }
};

template<class X> class VectorFunctionModelElement {
    VectorFunctionModel<X>* _p; Nat _i;
  public:
    operator const ScalarFunctionModel<X> () const;
    VectorFunctionModelElement(VectorFunctionModel<X>* p, Nat i) : _p(p), _i(i) { }
    VectorFunctionModelElement<X>& operator=(const ScalarFunctionModel<X>& sf) { _p->set(_i,sf); return *this; }
    VectorFunctionModelElement<X>& operator=(const VectorFunctionModelElement<X>& sf) { return this->operator=(static_cast<ScalarFunctionModel<X>const>(sf)); }
    Void clobber() { ScalarFunctionModel<Interval> sf=_p->get(_i); sf.clobber(); _p->set(_i,sf); }
    Void set_error(Float e) const { ScalarFunctionModel<Interval> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
};
template<class X> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionModelElement<X>& function) {
    return os << static_cast< const ScalarFunctionModel<X> >(function);
}

template<> class VectorFunctionModel<Interval>
{
  public:
    clone_on_copy_ptr< VectorFunctionModelInterface<Interval> > _ptr;
  public:
    inline VectorFunctionModel() : _ptr() { }
    inline VectorFunctionModel(Nat n, const ScalarFunctionModelInterface<Interval>& sf)
        : _ptr(sf._create_vector(n)) { for(uint i=0; i!=n; ++i) { (*this)[i]=sf; } }
    inline VectorFunctionModel(VectorFunctionModelInterface<Interval>* p) : _ptr(p) { }
    inline VectorFunctionModel(const VectorFunctionModelInterface<Interval>& f) : _ptr(f._clone()) { }
    inline VectorFunctionModel(const VectorFunctionModel<Interval>& f) : _ptr(f._ptr) { }
    inline operator const VectorFunctionModelInterface<Interval>& () const { return *_ptr; }
    inline operator VectorFunction<Interval> () const { return VectorFunction<Interval>(*_ptr); }
    inline const VectorFunctionModelInterface<Interval>* raw_pointer() const { return _ptr.operator->(); }
    inline const VectorFunctionModelInterface<Interval>& reference() const { return *_ptr; }
    inline VectorFunctionModelInterface<Interval>& reference() { return *_ptr; }
    inline VectorFunctionModel<Interval> create_identity() const { return this->_ptr->_create_identity(); }
    inline Nat result_size() const { return this->_ptr->result_size(); }
    inline Nat argument_size() const { return this->_ptr->argument_size(); }
    inline Nat size() const { return this->_ptr->result_size(); }
    template<class XX> inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    template<class XX> inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->evaluate(v); }
    inline ScalarFunctionModel<Interval> const get(Nat i) const { return this->_ptr->_get(i); }
    inline Void set(Nat i, ScalarFunctionModel<Interval> const& sf) { this->_ptr->_set(i,sf); }
    inline ScalarFunctionModel<Interval> const operator[](Nat i) const { return this->get(i); }
    inline VectorFunctionModelElement<Interval> operator[](Nat i) { return VectorFunctionModelElement<Interval>(this,i); }
    inline IntervalVector const& domain() const { return this->_ptr->domain(); }
    inline IntervalVector const range() const { return this->_ptr->range(); }
    inline IntervalVector const codomain() const { return this->_ptr->codomain(); }
    inline Vector<Float> const errors() const { return this->_ptr->errors(); }
    inline Float const error() const { return this->_ptr->error(); }
    inline Matrix<Interval> const jacobian(const Vector<Interval>& x) const {
        return this->_ptr->evaluate(Differential<Interval>::variables(1u,x)).jacobian(); }

    inline Void restrict(const IntervalVector& d) { this->_ptr->restrict(d); }

};

inline Float norm(const VectorFunctionModel<Interval>& f) {
    return f._ptr->_norm(); }
inline VectorFunctionModel<Interval> embed(const IntervalVector& d1, const VectorFunctionModel<Interval>& f, const IntervalVector& d2) {
    return f._ptr->_embed(d1,d2); }
inline VectorFunctionModel<Interval> embed(const VectorFunctionModel<Interval>& f, const IntervalVector& d) {
    return embed(IntervalVector(),f,d); }
inline VectorFunctionModel<Interval> embed(const VectorFunctionModel<Interval>& f, const Interval& d) {
    return embed(f,IntervalVector(1,d)); }
inline VectorFunctionModel<Interval> restrict(const VectorFunctionModel<Interval>& f, const IntervalVector& d) {
    VectorFunctionModelInterface<Interval>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }

inline VectorFunctionModel<Interval> operator+(const VectorFunctionModel<Interval>& f) {
    return f._ptr->_clone(); }
inline VectorFunctionModel<Interval> operator-(const VectorFunctionModel<Interval>& f) {
    VectorFunctionModel<Interval> r=f; for(uint i=0; i!=r.size(); ++i) { r[i]=-f[i]; } return r; }
inline VectorFunctionModel<Interval> operator+(const VectorFunctionModel<Interval>& f1, const VectorFunctionModel<Interval>& f2) {
    VectorFunctionModel<Interval> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
inline VectorFunctionModel<Interval> operator-(const VectorFunctionModel<Interval>& f1, const VectorFunctionModel<Interval>& f2) {
    VectorFunctionModel<Interval> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
inline VectorFunctionModel<Interval> operator+(const VectorFunctionModel<Interval>& f1, const Vector<Interval>& c2) {
    VectorFunctionModel<Interval> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]+c2[i]; } return r; }
inline VectorFunctionModel<Interval> operator-(const VectorFunctionModel<Interval>& f1, const Vector<Interval>& c2) {
    VectorFunctionModel<Interval> r=f1; for(uint i=0; i!=r.size(); ++i) { r[i]=f1[i]-c2[i]; } return r; }
inline VectorFunctionModel<Interval> operator+(const Vector<Interval>& c1, const VectorFunctionModel<Interval>& f2);
inline VectorFunctionModel<Interval> operator-(const Vector<Interval>& c1, const VectorFunctionModel<Interval>& f2);

inline Interval evaluate(const ScalarFunction<Interval>& f, const Vector<Interval>& x);
inline Vector<Interval> evaluate(const VectorFunction<Interval>& f, const Vector<Interval>& x);

inline Interval unchecked_evaluate(const ScalarFunctionModel<Interval>& f, const Vector<Interval>& x) { return f._ptr->_unchecked_evaluate(x); }
inline Vector<Interval> unchecked_evaluate(const VectorFunctionModel<Interval>& f, const Vector<Interval>& x) { return f._ptr->_unchecked_evaluate(x); }

inline ScalarFunctionModel<Interval> compose(const ScalarFunction<Interval>& f, const VectorFunctionModel<Interval>& g) { return g._ptr->_compose(f); }
inline VectorFunctionModel<Interval> compose(const VectorFunction<Interval>& f, const VectorFunctionModel<Interval>& g) { return g._ptr->_compose(f); }

inline ScalarFunctionModel<Interval> unchecked_compose(const ScalarFunctionModel<Interval>& f, const VectorFunctionModel<Interval>& g) { return g._ptr->_unchecked_compose(f); }
inline VectorFunctionModel<Interval> unchecked_compose(const VectorFunctionModel<Interval>& f, const VectorFunctionModel<Interval>& g) { return g._ptr->_unchecked_compose(f); }

//inline VectorFunctionModel<Interval> compose(const VectorFunctionModel<Interval>& f, const VectorFunctionModel<Interval>& g) { return g._ptr->_compose(f); }

inline VectorFunctionModel<Interval> operator-(const VectorFunctionModel<Interval>& f1, const VectorFunctionInterface<Interval>& f2) {
    return f1-compose(f2,f1.create_identity()); }

inline VectorFunctionModel<Interval> join(const ScalarFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2);
inline VectorFunctionModel<Interval> join(const ScalarFunctionModel<Interval>& f1, const VectorFunctionModel<Interval>& f2);
inline VectorFunctionModel<Interval> join(const VectorFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) {
    VectorFunctionModel<Interval> r=f1._ptr->_clone(); r._ptr->_adjoin(f2); return r; }
inline VectorFunctionModel<Interval> join(const VectorFunctionModel<Interval>& f1, const VectorFunctionModel<Interval>& f2) {
    return f1._ptr->_join(f2); }

inline VectorFunctionModel<Interval> combine(const ScalarFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) {
    return VectorFunctionModel<Interval>(1,f1)._ptr->_combine(VectorFunctionModel<Interval>(1,f2)); };;
inline VectorFunctionModel<Interval> combine(const ScalarFunctionModel<Interval>& f1, const VectorFunctionModel<Interval>& f2) {
    return VectorFunctionModel<Interval>(1,f1)._ptr->_combine(f2); };;
inline VectorFunctionModel<Interval> combine(const VectorFunctionModel<Interval>& f1, const ScalarFunctionModel<Interval>& f2) {
    return f1._ptr->_combine(VectorFunctionModel<Interval>(1,f2)); };
inline VectorFunctionModel<Interval> combine(const VectorFunctionModel<Interval>& f1, const VectorFunctionModel<Interval>& f2) {
    return f1._ptr->_combine(f2); }

inline VectorFunctionModel<Interval> intersection(const VectorFunctionModel<Interval>& f1, const VectorFunctionModel<Interval>& f2) {
    ARIADNE_ASSERT_MSG(f1.size()==f2.size(),"intersection(f1,f2): f1="<<f1<<", f2="<<f2<<")");
    VectorFunctionModel<Interval> r=+f1; for(uint i=0; i!=r.size(); ++i) { r[i]=intersection(f1[i],f2[i]); } return r; }

inline VectorFunctionModel<Interval> partial_evaluate(const VectorFunctionModel<Interval>& f, Nat j, const Interval& c) {
    return f._ptr->_partial_evaluate(j,c); }


template<class X> VectorFunctionModelElement<X>::operator const ScalarFunctionModel<X> () const {
    return _p->get(_i); }




// Sanitised output
template<class T, class D=Void> struct Representation;
template<class T, class D> struct Representation { const T* pointer; D data; };
template<class T> struct Representation<T> { const T* pointer; };
template<class T> Representation<T> repr(const T& t) { Representation<T> r={&t}; return r; }
template<class T, class D> Representation<T,D> repr(const T& t, const D& d) { Representation<T,D> r={&t,d}; return r; }

template<class X> class FunctionModelFactoryInterface;

template<> class FunctionModelFactoryInterface<Interval>
{
    typedef IntervalVector DomainType;
  public:
    virtual FunctionModelFactoryInterface<Interval>* clone() const = 0;
    virtual Void write(OutputStream& os) const = 0;
    inline ScalarFunctionModel<Interval> create(const IntervalVector& domain, const ScalarFunctionInterface<Interval>& function) const;
    inline VectorFunctionModel<Interval> create(const IntervalVector& domain, const VectorFunctionInterface<Interval>& function) const;
    inline ScalarFunctionModel<Interval> create_zero(const IntervalVector& domain) const;
    inline VectorFunctionModel<Interval> create_zeros(Nat result_size, const IntervalVector& domain) const;
    inline ScalarFunctionModel<Interval> create_constant(const IntervalVector& domain, const Interval& value) const;
    inline ScalarFunctionModel<Interval> create_constant(const IntervalVector& domain, const Float& value) const;
    inline ScalarFunctionModel<Interval> create_constant(const IntervalVector& domain, double value) const;
    inline ScalarFunctionModel<Interval> create_coordinate(const IntervalVector& domain, Nat index) const;
    inline ScalarFunctionModel<Interval> create_identity(const Interval& domain) const;
    inline VectorFunctionModel<Interval> create_identity(const IntervalVector& domain) const;
  private:
    virtual ScalarFunctionModelInterface<Interval>* _create(const IntervalVector& domain, const ScalarFunctionInterface<Interval>& function) const = 0;
    virtual VectorFunctionModelInterface<Interval>* _create(const IntervalVector& domain, const VectorFunctionInterface<Interval>& function) const = 0;
};

inline OutputStream& operator<<(OutputStream& os, const FunctionModelFactoryInterface<Interval>& factory) {
    factory.write(os); return os;
}

inline ScalarFunctionModel<Interval>
FunctionModelFactoryInterface<Interval>::create(const IntervalVector& domain,
                                                const ScalarFunctionInterface<Interval>& function) const
{
    return this->_create(domain,function);
}

VectorFunctionModel<Interval>
FunctionModelFactoryInterface<Interval>::create(const IntervalVector& domain,
                                                const VectorFunctionInterface<Interval>& function) const
{
    return this->_create(domain,function);
}


} // namespace Ariadne

#include "function.h"

namespace Ariadne {

ScalarFunctionModel<Interval> FunctionModelFactoryInterface<Interval>::create_zero(const IntervalVector& domain) const {
    return this->_create(domain,RealScalarFunction::zero(domain.size())); }
VectorFunctionModel<Interval> FunctionModelFactoryInterface<Interval>::create_zeros(Nat result_size, const IntervalVector& domain) const {
    return this->_create(domain,RealVectorFunction(result_size,domain.size())); }
ScalarFunctionModel<Interval> FunctionModelFactoryInterface<Interval>::create_constant(const IntervalVector& domain, double value) const {
    return this->create_constant(domain,numeric_cast<Interval>(value)); }
ScalarFunctionModel<Interval> FunctionModelFactoryInterface<Interval>::create_constant(const IntervalVector& domain, const Float& value) const {
    return this->create_constant(domain,numeric_cast<Interval>(value)); }
ScalarFunctionModel<Interval> FunctionModelFactoryInterface<Interval>::create_constant(const IntervalVector& domain, const Interval& value) const {
    return ScalarFunctionModel<Interval>(this->_create(domain,RealScalarFunction::zero(domain.size())))+value; };
ScalarFunctionModel<Interval> FunctionModelFactoryInterface<Interval>::create_coordinate(const IntervalVector& domain, Nat index) const {
    return ScalarFunctionModel<Interval>(this->_create(domain,RealScalarFunction::coordinate(domain.size(),index))); };
ScalarFunctionModel<Interval> FunctionModelFactoryInterface<Interval>::create_identity(const Interval& domain) const {
    return this->_create(IntervalVector(1,domain),RealScalarFunction::coordinate(1,0)); };
VectorFunctionModel<Interval> FunctionModelFactoryInterface<Interval>::create_identity(const IntervalVector& domain) const {
    return this->_create(domain,RealVectorFunction::identity(domain.size())); };




} // namespace Ariadne

#endif
