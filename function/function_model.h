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

#include "function/function_interface.h"
#include "function/function_mixin.h"
#include "function/function.h"

#include "numeric/operators.h"
#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/operations.h"
#include "geometry/box.h"

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

typedef Float64Value ValidatedCoefficientType;
typedef Float64Error ValidatedErrorType;

template<class P, class F> class TaylorModel;
template<class M> class FunctionPatch;
typedef FunctionPatch<TaylorModel<ValidatedTag,Float64>> ScalarTaylorFunction;

template<> class ScalarFunctionModelInterface<ValidatedTag>
    : public virtual ScalarFunctionInterface<ValidatedTag>
{
  public:
    typedef Ariadne::ValidatedCoefficientType CoefficientType;
    typedef Ariadne::ValidatedErrorType ErrorType;
  public:
    virtual UpperIntervalType range() const = 0;

    virtual CoefficientType const& value() const = 0;
    virtual CoefficientType const gradient_value(SizeType i) const = 0;
    virtual ErrorType const& error() const = 0;

    virtual Void set_error(const ErrorType& e) = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual ScalarFunctionModelInterface<ValidatedTag>* _apply(OperatorCode op) const = 0;
    virtual ValidatedNumericType _unchecked_evaluate(const Vector<ValidatedNumericType>& x) const = 0;

    virtual ScalarFunctionModelInterface<ValidatedTag>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create() const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create_vector(SizeType i) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create_identity() const = 0;
    virtual Void restrict(const ExactBoxType& d) = 0;

    virtual ScalarFunctionModelInterface<ValidatedTag>* _create_zero(DomainType const& dom) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create_constant(DomainType const& dom, ValidatedNumericType const& c) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create_coordinate(DomainType const& dom, SizeType j) const = 0;

    virtual ScalarFunctionModelInterface<ValidatedTag>* _derivative(SizeType j) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _antiderivative(SizeType j) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _antiderivative(SizeType j, ValidatedNumericType c) const = 0;

    virtual Boolean _refines(const ScalarFunctionModelInterface<ValidatedTag>& f) const = 0;
    virtual Boolean _disjoint(const ScalarFunctionModelInterface<ValidatedTag>& f) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _intersection(const ScalarFunctionModelInterface<ValidatedTag>& f) const = 0;

    virtual Void _iadd(const ValidatedNumericType& c) = 0;
    virtual Void _imul(const ValidatedNumericType& c) = 0;
    virtual Void _isma(const ValidatedNumericType& c, const ScalarFunctionModelInterface<ValidatedTag>& f) = 0;
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
    ScalarFunctionModelInterface<ValidatedTag>* _antiderivative(SizeType j) const {
        return new F(antiderivative(static_cast<const F&>(*this),j)); }
    ScalarFunctionModelInterface<ValidatedTag>* _antiderivative(SizeType j, ValidatedNumericType c) const {
        return new F(antiderivative(static_cast<const F&>(*this),j,c)); }
    ScalarFunctionModelInterface<ValidatedTag>* _apply(OperatorCode op) const {
        return new F(this->apply(op)); }
    ValidatedNumericType _unchecked_evaluate(const Vector<ValidatedNumericType>& x) const {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<ValidatedTag>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const {
        return new F(embed(d1,static_cast<const F&>(*this),d2)); }
    Boolean _refines(const ScalarFunctionModelInterface<ValidatedTag>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return refines(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    Boolean _disjoint(const ScalarFunctionModelInterface<ValidatedTag>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return disjoint(static_cast<const F&>(*this),dynamic_cast<const F&>(f)); }
    ScalarFunctionModelInterface<ValidatedTag>* _intersection(const ScalarFunctionModelInterface<ValidatedTag>& f) const {
        ARIADNE_ASSERT(dynamic_cast<const F*>(&f)); return new F(intersection(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Void _iadd(const ValidatedNumericType& c) {
        static_cast<F&>(*this)+=c; }
    Void _imul(const ValidatedNumericType& c) {
        static_cast<F&>(*this)*=c; }
    Void _isma(const ValidatedNumericType& c, const ScalarFunctionModelInterface<ValidatedTag>& f) {
        static_cast<F&>(*this)+=c*dynamic_cast<const F&>(f); }
    Void _ifma(const ScalarFunctionModelInterface<ValidatedTag>& f1, const ScalarFunctionModelInterface<ValidatedTag>& f2) {
        static_cast<F&>(*this)+=dynamic_cast<const F&>(f1)*dynamic_cast<const F&>(f2); }
};

//! \ingroup FunctionModelSubModule
//! \brief Generic scalar functions on singleton domains.
template<> class ScalarFunctionModel<ValidatedTag>
{
  public:
    typedef ExactBoxType DomainType;
    typedef Ariadne::ValidatedCoefficientType CoefficientType;
    typedef Ariadne::ValidatedErrorType ErrorType;
  public:
    clone_on_copy_ptr< ScalarFunctionModelInterface<ValidatedTag> > _ptr;
  public:
    ScalarFunctionModel() : _ptr() { }
    ScalarFunctionModel(ScalarFunctionModelInterface<ValidatedTag>* p) : _ptr(p) { }
    ScalarFunctionModel(const shared_ptr<const ScalarFunctionModelInterface<ValidatedTag>> p) : _ptr(p->_clone()) { }
    ScalarFunctionModel(const ScalarFunctionModel<ValidatedTag>& f) : _ptr(f._ptr) { }
    ScalarFunctionModel(const ScalarFunctionModelInterface<ValidatedTag>& f) : _ptr(f._clone()) { }
    ScalarFunctionModel(const ValidatedScalarFunction& f) : _ptr(dynamic_cast<ScalarFunctionModelInterface<ValidatedTag>*>(f.raw_pointer()->_clone())) { }
    operator ValidatedScalarFunction() const { return ValidatedScalarFunction(this->_ptr->_clone()); }
    operator ScalarFunctionModelInterface<ValidatedTag>& () { return *_ptr; }
    operator const ScalarFunctionModelInterface<ValidatedTag>& () const { return *_ptr; }
    const ScalarFunctionModelInterface<ValidatedTag>* raw_pointer() const { return _ptr.operator->(); }
    ScalarFunctionModelInterface<ValidatedTag>& reference() { return *_ptr; }
    const ScalarFunctionModelInterface<ValidatedTag>& reference() const { return *_ptr; }
    ScalarFunctionModel<ValidatedTag> create_zero() const { return ScalarFunctionModel<ValidatedTag>(this->_ptr->_create_zero(this->domain())); }
    ScalarFunctionModel<ValidatedTag> create_constant(const ValidatedNumericType& c) const;
    VectorFunctionModel<ValidatedTag> create_identity() const;
    ScalarFunctionModel<ValidatedTag> create(const ValidatedScalarFunction& f) const;
    Vector<ScalarFunctionModel<ValidatedTag>> create_coordinates(DomainType const&) const;
    ScalarFunctionModel<ValidatedTag>& operator=(const ValidatedNumber& c);
    ScalarFunctionModel<ValidatedTag>& operator=(const ValidatedNumericType& c);
    ScalarFunctionModel<ValidatedTag>& operator=(const ValidatedScalarFunction& f);
    ScalarFunctionModel<ValidatedTag>& operator=(const ScalarTaylorFunction& f);
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    template<class XX> XX operator() (const Vector<XX>& x) const {
        return this->_ptr->_evaluate(x); }
    template<class XX> XX evaluate(const Vector<XX>& x) const {
        return this->_ptr->_evaluate(x); }
    inline ExactBoxType const domain() const { return this->_ptr->domain(); }
    inline ExactIntervalType const codomain() const { return this->_ptr->codomain(); }
    inline UpperIntervalType const range() const { return this->_ptr->range(); }

    inline CoefficientType value() const { return this->_ptr->value(); }
    inline CoefficientType gradient_value(SizeType j) const { return this->_ptr->gradient_value(j); }
    inline ErrorType error() const { return this->_ptr->error(); }

    inline Void set_error(const ErrorType& e) { return this->_ptr->set_error(e); }
    inline Void set_error(Nat e) { return this->_ptr->set_error(ErrorType(e,this->error().precision())); }
    inline Void clobber() { return this->_ptr->clobber(); }

    inline ScalarFunctionModel<ValidatedTag> apply(Operator op) const { return this->_ptr->_apply(op); }
    inline Void restrict(const ExactBoxType& d) { this->_ptr->restrict(d); }
  public:
    friend inline ScalarFunctionModel<ValidatedTag>& operator+=(ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) { f1._ptr->_isma(ValidatedNumericType(+1),f2); return  f1; }
    friend inline ScalarFunctionModel<ValidatedTag>& operator-=(ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) { f1._ptr->_isma(ValidatedNumericType(-1),f2); return  f1; }
    friend inline ScalarFunctionModel<ValidatedTag>& operator+=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) { f1._ptr->_iadd(c2); return f1; }
    friend inline ScalarFunctionModel<ValidatedTag>& operator-=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) { f1._ptr->_iadd(neg(c2)); return f1; }
    friend inline ScalarFunctionModel<ValidatedTag>& operator*=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) { f1._ptr->_imul(c2); return f1; }
    friend inline ScalarFunctionModel<ValidatedTag>& operator/=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) { f1._ptr->_imul(rec(c2)); return f1; }
    friend inline ScalarFunctionModel<ValidatedTag>& operator+=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) { return f1+=ValidatedNumericType(c2,f1.value().precision()); }
    friend inline ScalarFunctionModel<ValidatedTag>& operator-=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) { return f1-=ValidatedNumericType(c2,f1.value().precision()); }
    friend inline ScalarFunctionModel<ValidatedTag>& operator*=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) { return f1*=ValidatedNumericType(c2,f1.value().precision()); }
    friend inline ScalarFunctionModel<ValidatedTag>& operator/=(ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) { return f1/=ValidatedNumericType(c2,f1.value().precision()); }

    friend inline ScalarFunctionModel<ValidatedTag> pos(const ScalarFunctionModel<ValidatedTag>& f) { return f.apply(Pos()); }
    friend inline ScalarFunctionModel<ValidatedTag> neg(const ScalarFunctionModel<ValidatedTag>& f) { return f.apply(Neg()); }
    friend inline ScalarFunctionModel<ValidatedTag> sqr(const ScalarFunctionModel<ValidatedTag>& f) { return f.apply(Sqr()); }
    friend inline ScalarFunctionModel<ValidatedTag> rec(const ScalarFunctionModel<ValidatedTag>& f) { return f.apply(Rec()); }
    friend inline ScalarFunctionModel<ValidatedTag> add(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r+=f2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> sub(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r-=f2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> mul(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f1.create_zero(); r._ptr->_ifma(f1,f2); return r; }
    friend inline ScalarFunctionModel<ValidatedTag> div(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
        return mul(f1,rec(f2)); }

    friend inline ScalarFunctionModel<ValidatedTag> operator+(const ScalarFunctionModel<ValidatedTag>& f) {
        return f._ptr->_clone(); }
    friend inline ScalarFunctionModel<ValidatedTag> operator-(const ScalarFunctionModel<ValidatedTag>& f) {
        ScalarFunctionModel<ValidatedTag> r=f; r*=-1; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator+(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r+=f2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator-(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r-=f2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator*(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f1.create_zero(); r._ptr->_ifma(f1,f2); return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator/(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
        return f1*rec(f2); }

    friend inline ScalarFunctionModel<ValidatedTag> add(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r+=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> sub(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r-=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> mul(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r*=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> div(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r/=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> add(const ValidatedNumericType& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f2; r+=c1; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> sub(const ValidatedNumericType& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=neg(f2); r+=c1; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> mul(const ValidatedNumericType& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f2; r*=c1; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> div(const ValidatedNumericType& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=rec(f2); r*=c1; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> pow(const ScalarFunctionModel<ValidatedTag>& f, Nat m) {
        ScalarFunctionModel<ValidatedTag> r=f; r*=0; r+=1; ScalarFunctionModel<ValidatedTag> p=f; while(m!=0) { if(m%2==1) { r=r*p; } p=p*p; m=m/2; } return r; }
    friend inline ScalarFunctionModel<ValidatedTag> pow(const ScalarFunctionModel<ValidatedTag>& f, Int n) {
        return n>=0 ? pow(f,uint(n)) : rec(pow(f,uint(-n))); }

    friend inline ScalarFunctionModel<ValidatedTag> operator+(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r+=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator-(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r-=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator*(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r*=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator/(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r/=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator+(const ValidatedNumericType& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f2; r+=c1; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator-(const ValidatedNumericType& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=neg(f2); r+=c1; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator*(const ValidatedNumericType& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f2; r*=c1; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator/(const ValidatedNumericType& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=rec(f2); r*=c1; return r; }

    friend inline ScalarFunctionModel<ValidatedTag> operator+(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r+=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator-(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r-=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator*(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r*=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator/(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedNumber& c2) {
        ScalarFunctionModel<ValidatedTag> r=f1; r/=c2; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator+(const ValidatedNumber& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f2; r+=c1; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator-(const ValidatedNumber& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=neg(f2); r+=c1; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator*(const ValidatedNumber& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=f2; r*=c1; return r; }
    friend inline ScalarFunctionModel<ValidatedTag> operator/(const ValidatedNumber& c1, const ScalarFunctionModel<ValidatedTag>& f2) {
        ScalarFunctionModel<ValidatedTag> r=rec(f2); r*=c1; return r; }

    friend inline ScalarFunctionModel<ValidatedTag> operator+(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedScalarFunction& g2) {
        return f1+f1.create(g2); }
    friend inline ScalarFunctionModel<ValidatedTag> operator-(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedScalarFunction& g2) {
        return f1-f1.create(g2); }
    friend inline ScalarFunctionModel<ValidatedTag> operator*(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedScalarFunction& g2) {
        return f1*f1.create(g2); }
    friend inline ScalarFunctionModel<ValidatedTag> operator/(const ScalarFunctionModel<ValidatedTag>& f1, const ValidatedScalarFunction& g2) {
        return f1/f1.create(g2); }
    friend inline ScalarFunctionModel<ValidatedTag> operator+(const ValidatedScalarFunction& g1, const ScalarFunctionModel<ValidatedTag>& f2) {
        return f2.create(g1)+f2; }
    friend inline ScalarFunctionModel<ValidatedTag> operator-(const ValidatedScalarFunction& g1, const ScalarFunctionModel<ValidatedTag>& f2) {
        return f2.create(g1)-f2; }
    friend inline ScalarFunctionModel<ValidatedTag> operator*(const ValidatedScalarFunction& g1, const ScalarFunctionModel<ValidatedTag>& f2) {
        return f2.create(g1)*f2; }
    friend inline ScalarFunctionModel<ValidatedTag> operator/(const ValidatedScalarFunction& g1, const ScalarFunctionModel<ValidatedTag>& f2) {
        return f2.create(g1)/f2; }

};

// inline ScalarFunctionModel<ValidatedTag>& ScalarFunctionModel<ValidatedTag>::operator=(const ValidatedScalarFunction& f) { (*this)=this->_ptr->_create()+f; return *this; }

inline NormType norm(const ScalarFunctionModel<ValidatedTag>& f) { return f._ptr->_norm(); }
inline ScalarFunctionModel<ValidatedTag> derivative(const ScalarFunctionModel<ValidatedTag>& f, SizeType j) { return f._ptr->_derivative(j); }
inline ScalarFunctionModel<ValidatedTag> antiderivative(const ScalarFunctionModel<ValidatedTag>& f, SizeType j) { return f._ptr->_antiderivative(j); }
inline ScalarFunctionModel<ValidatedTag> antiderivative(const ScalarFunctionModel<ValidatedTag>& f, SizeType j, ValidatedNumericType c) { return f._ptr->_antiderivative(j,c); }
inline ScalarFunctionModel<ValidatedTag> restriction(const ScalarFunctionModel<ValidatedTag>& f, const ExactBoxType& d) {
    ScalarFunctionModel<ValidatedTag> r(f); r.restrict(d); return r; }

inline ScalarFunctionModel<ValidatedTag> embed(const ExactBoxType& d1, const ScalarFunctionModel<ValidatedTag>& f, const ExactBoxType& d2) {
    return f._ptr->_embed(d1,d2); }
inline ScalarFunctionModel<ValidatedTag> embed(const ExactBoxType& d, const ScalarFunctionModel<ValidatedTag>& f) {
    return embed(d,f,ExactBoxType()); }
inline ScalarFunctionModel<ValidatedTag> embed(const ScalarFunctionModel<ValidatedTag>& f, const ExactBoxType& d) {
    return embed(ExactBoxType(),f,d); }
inline ScalarFunctionModel<ValidatedTag> embed(const ScalarFunctionModel<ValidatedTag>& f, const ExactIntervalType& d) {
    return embed(f,ExactBoxType(1,d)); }
inline ScalarFunctionModel<ValidatedTag> restrict(const ScalarFunctionModel<ValidatedTag>& f, const ExactBoxType& d) {
    ScalarFunctionModelInterface<ValidatedTag>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }
inline ScalarFunctionModel<ValidatedTag> intersection(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_intersection(f2); }

inline Boolean disjoint(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_disjoint(f2); }
inline Boolean refines(const ScalarFunctionModel<ValidatedTag>& f1, const ScalarFunctionModel<ValidatedTag>& f2) {
    return f1._ptr->_refines(f2); }


inline ScalarFunctionModel<ValidatedTag>& ScalarFunctionModel<ValidatedTag>::operator=(const ValidatedNumericType& c) { (*this)*=ValidatedNumericType(0); (*this)+=c; return *this; }

template<class F> F ScalarFunctionModelMixin<F,ValidatedTag>::apply(OperatorCode op) const {
    const F& f=static_cast<const F&>(*this);
    switch(op) {
        case OperatorCode::NEG: return neg(f);
        case OperatorCode::REC: return rec(f);
        case OperatorCode::EXP: return exp(f);
        default: ARIADNE_FAIL_MSG("ScalarFunctionModel<ValidatedTag>::apply(OperatorCode op): Operator op="<<op<<" not implemented\n");
    }
}

template<class P> inline OutputStream& operator<<(OutputStream& os, const ScalarFunctionModel<P>& f) {
    return os <<  f.operator ScalarFunction<P>(); }
//inline ScalarFunctionModel<ValidatedTag>::ScalarFunctionModel(const ValidatedScalarFunction& f)
//    : _ptr(dynamic_cast<const ScalarFunctionModelInterface<ValidatedTag>&>(*f.raw_pointer())._clone()) { }




//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on singleton domains.
template<> class VectorFunctionModelInterface<ValidatedTag>
    : public virtual VectorFunctionInterface<ValidatedTag>
{
  public:
    virtual UpperBoxType const range() const = 0;
    virtual Vector<ErrorType> const errors() const = 0;
    virtual ErrorType const error() const = 0;
    virtual Void clobber() = 0;

    virtual NormType const _norm() const = 0;

    virtual VectorFunctionModelInterface<ValidatedTag>* _clone() const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create_zero() const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create_identity() const = 0;
    virtual Void _set(SizeType, ScalarFunctionModelInterface<ValidatedTag> const&) = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _get(SizeType) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _join(const VectorFunctionModelInterface<ValidatedTag>& f2) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _combine(const VectorFunctionModelInterface<ValidatedTag>& f2) const = 0;
    virtual Void _adjoin(const ScalarFunctionModelInterface<ValidatedTag>& f2) = 0;
    virtual Vector<ValidatedNumericType> _unchecked_evaluate(const Vector<ValidatedNumericType>& x) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _compose(const ScalarFunctionInterface<ValidatedTag>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _compose(const VectorFunctionInterface<ValidatedTag>& f) const = 0;
    virtual ScalarFunctionModelInterface<ValidatedTag>* _unchecked_compose(const ScalarFunctionInterface<ValidatedTag>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _unchecked_compose(const VectorFunctionInterface<ValidatedTag>& f) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _partial_evaluate(SizeType j, const ValidatedNumericType& c) const = 0;
    virtual Void restrict(const ExactBoxType& d) = 0;
};

template<class V> struct Element;

template<class M> class FunctionPatch;
template<class M> class VectorFunctionPatch;
template<class M> struct Element<VectorFunctionPatch<M>> { typedef FunctionPatch<M> Type; };

typedef FunctionPatch<ValidatedTaylorModel> ScalarTaylorFunction;
//
template<class F> class VectorFunctionModelMixin<F,ValidatedTag>
    : public virtual VectorFunctionModelInterface<ValidatedTag>
    , public  VectorFunctionMixin<F,ValidatedTag>
{
    typedef typename Element<F>::Type ScalarFunctionType;
  public:
    virtual VectorFunctionModelInterface<ValidatedTag>* _clone() const { return new F(static_cast<const F&>(*this)); }
    virtual Void _set(SizeType i, const ScalarFunctionModelInterface<ValidatedTag>& sf) {
        if(!dynamic_cast<const typename F::ScalarFunctionType*>(&sf)) {
            ARIADNE_FAIL_MSG("Cannot set element of VectorFunctionModel "<<*this<<" to "<<sf<<"\n"); }
        static_cast<F&>(*this).F::set(i,dynamic_cast<const ScalarFunctionType&>(sf)); }
    virtual VectorFunctionModelInterface<ValidatedTag>* _derivative(SizeType j) const {
        ARIADNE_NOT_IMPLEMENTED; }
    NormType const _norm() const {
         return norm(static_cast<const F&>(*this)); }
    VectorFunctionModelInterface<ValidatedTag>* _embed(const ExactBoxType& d1, const ExactBoxType& d2) const {
        return heap_copy(embed(d1,static_cast<const F&>(*this),d2)); }
    Void _adjoin(const ScalarFunctionModelInterface<ValidatedTag>& f) {
        static_cast<F&>(*this).F::adjoin(dynamic_cast<const ScalarFunctionType&>(f)); }
    VectorFunctionModelInterface<ValidatedTag>* _join(const VectorFunctionModelInterface<ValidatedTag>& f) const {
        return heap_copy(join(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    VectorFunctionModelInterface<ValidatedTag>* _combine(const VectorFunctionModelInterface<ValidatedTag>& f) const {
        return heap_copy(combine(static_cast<const F&>(*this),dynamic_cast<const F&>(f))); }
    Vector<ValidatedNumericType> _unchecked_evaluate(const Vector<ValidatedNumericType>& x) const {
        return unchecked_evaluate(static_cast<const F&>(*this),x); }
    ScalarFunctionModelInterface<ValidatedTag>* _compose(const ScalarFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedTag>* _compose(const VectorFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(compose(f,static_cast<const F&>(*this))); }
    ScalarFunctionModelInterface<ValidatedTag>* _unchecked_compose(const ScalarFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(unchecked_compose(dynamic_cast<const ScalarFunctionType&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedTag>* _unchecked_compose(const VectorFunctionInterface<ValidatedTag>& f) const {
        return heap_copy(unchecked_compose(dynamic_cast<const F&>(f),static_cast<const F&>(*this))); }
    VectorFunctionModelInterface<ValidatedTag>* _partial_evaluate(SizeType j, const ValidatedNumericType& c) const {
        return heap_copy(partial_evaluate(static_cast<const F&>(*this),j,c)); }
};

template<class X> class VectorFunctionModelElement
    : public DeclareArithmeticOperators<ScalarFunctionModel<ValidatedTag>>
    , public DeclareMixedArithmeticOperators<ScalarFunctionModel<ValidatedTag>,ValidatedNumericType>
    , public DeclareMixedArithmeticOperators<ScalarFunctionModel<ValidatedTag>,ValidatedScalarFunction>
{
    VectorFunctionModel<X>* _p; SizeType _i;
  public:
    operator const ScalarFunctionModel<X> () const;
    VectorFunctionModelElement(VectorFunctionModel<X>* p, SizeType i) : _p(p), _i(i) { }
    VectorFunctionModelElement<X>& operator=(const ScalarFunctionModel<X>& sf) { _p->set(_i,sf); return *this; }
    VectorFunctionModelElement<X>& operator=(const VectorFunctionModelElement<X>& sf) { return this->operator=(static_cast<ScalarFunctionModel<X>const>(sf)); }
    Void clobber() { ScalarFunctionModel<ValidatedTag> sf=_p->get(_i); sf.clobber(); _p->set(_i,sf); }
    Void set_error(ErrorType e) const { ScalarFunctionModel<ValidatedTag> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
    Void set_error(Nat e) const { ScalarFunctionModel<ValidatedTag> sf=_p->get(_i); sf.set_error(e); _p->set(_i,sf); }
};
template<class X> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionModelElement<X>& function) {
    return os << static_cast< const ScalarFunctionModel<X> >(function);
}

//! \ingroup FunctionModelSubModule
//! \brief Generic vector functions on singleton domains.
template<> class VectorFunctionModel<ValidatedTag>
{
  public:
    clone_on_copy_ptr< VectorFunctionModelInterface<ValidatedTag> > _ptr;
  public:
    inline VectorFunctionModel() : _ptr() { }
    inline VectorFunctionModel(SizeType n, const ScalarFunctionModelInterface<ValidatedTag>& sf)
        : _ptr(sf._create_vector(n)) { for(SizeType i=0; i!=n; ++i) { (*this)[i]=sf; } }
    inline VectorFunctionModel(Array<ScalarFunctionModel<ValidatedTag>> const& asf)
        : VectorFunctionModel(asf.size(),asf[0]) { for(SizeType i=0; i!=asf.size(); ++i) { (*this)[i]=asf[i]; } }
    inline VectorFunctionModel(List<ScalarFunctionModel<ValidatedTag>> const& lsf)
        : VectorFunctionModel(lsf.size(),lsf[0]) { for(SizeType i=0; i!=lsf.size(); ++i) { (*this)[i]=lsf[i]; } }
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
    inline VectorFunctionModel<ValidatedTag> create(VectorFunction<ValidatedTag> const& vf) const { return compose(vf,this->create_identity()); }
    inline SizeType result_size() const { return this->_ptr->result_size(); }
    inline SizeType argument_size() const { return this->_ptr->argument_size(); }
    inline SizeType size() const { return this->_ptr->result_size(); }
    template<class XX> inline Vector<XX> operator()(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    template<class XX> inline Vector<XX> evaluate(const Vector<XX>& v) const { return this->_ptr->_evaluate(v); }
    inline ScalarFunctionModel<ValidatedTag> const get(SizeType i) const { return this->_ptr->_get(i); }
    inline Void set(SizeType i, ScalarFunctionModel<ValidatedTag> const& sf) { this->_ptr->_set(i,sf); }
    inline ScalarFunctionModel<ValidatedTag> const operator[](SizeType i) const { return this->get(i); }
    inline VectorFunctionModelElement<ValidatedTag> operator[](SizeType i) { return VectorFunctionModelElement<ValidatedTag>(this,i); }
    inline ExactBoxType const domain() const { return this->_ptr->domain(); }
    inline ExactBoxType const codomain() const { return this->_ptr->codomain(); }
    inline UpperBoxType const range() const { return this->_ptr->range(); }
    inline Vector<ErrorType> const errors() const { return this->_ptr->errors(); }
    inline ErrorType const error() const { return this->_ptr->error(); }
    inline Void clobber() { this->_ptr->clobber(); }
    inline Matrix<ValidatedNumericType> const jacobian(const Vector<ValidatedNumericType>& x) const {
        Vector<Differential<ValidatedNumericType>> dfx=this->_ptr->_evaluate(Differential<ValidatedNumericType>::variables(1u,x));
        return dfx.jacobian(); }

    inline Void restrict(const ExactBoxType& d) { this->_ptr->restrict(d); }
  public:
    friend inline ScalarFunctionModel<ValidatedTag> compose(const ValidatedScalarFunction& f, const VectorFunctionModel<ValidatedTag>& g) {
        return g._ptr->_compose(f); }
    friend inline ScalarFunctionModel<ValidatedTag> compose(const ScalarFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag>& g) {
        return g._ptr->_compose(f); }
    friend inline VectorFunctionModel<ValidatedTag> compose(const ValidatedVectorFunction& f, const VectorFunctionModel<ValidatedTag>& g) {
        return g._ptr->_compose(f); }
    friend inline VectorFunctionModel<ValidatedTag> compose(const VectorFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag>& g) {
        return g._ptr->_compose(f); }

    friend inline ScalarFunctionModel<ValidatedTag> unchecked_compose(const ScalarFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag>& g) {
        return g._ptr->_unchecked_compose(f); }
    friend inline VectorFunctionModel<ValidatedTag> unchecked_compose(const VectorFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag>& g) {
        return g._ptr->_unchecked_compose(f); }

    friend inline VectorFunctionModel<ValidatedTag> operator+(const VectorFunctionModel<ValidatedTag>& f) {
        return f._ptr->_clone(); }
    friend inline VectorFunctionModel<ValidatedTag> operator-(const VectorFunctionModel<ValidatedTag>& f) {
        VectorFunctionModel<ValidatedTag> r=f; for(SizeType i=0; i!=r.size(); ++i) { r[i]=-f[i]; } return r; }
    friend inline VectorFunctionModel<ValidatedTag> operator+(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
        VectorFunctionModel<ValidatedTag> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+f2[i]; } return r; }
    friend inline VectorFunctionModel<ValidatedTag> operator-(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
        VectorFunctionModel<ValidatedTag> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-f2[i]; } return r; }
    friend inline VectorFunctionModel<ValidatedTag> operator+(const VectorFunctionModel<ValidatedTag>& f1, const Vector<ValidatedNumericType>& c2) {
        VectorFunctionModel<ValidatedTag> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]+c2[i]; } return r; }
    friend inline VectorFunctionModel<ValidatedTag> operator-(const VectorFunctionModel<ValidatedTag>& f1, const Vector<ValidatedNumericType>& c2) {
        VectorFunctionModel<ValidatedTag> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]-c2[i]; } return r; }
    friend inline VectorFunctionModel<ValidatedTag> operator+(const Vector<ValidatedNumericType>& c1, const VectorFunctionModel<ValidatedTag>& f2);
    friend inline VectorFunctionModel<ValidatedTag> operator-(const Vector<ValidatedNumericType>& c1, const VectorFunctionModel<ValidatedTag>& f2);
    friend inline VectorFunctionModel<ValidatedTag> operator*(const VectorFunctionModel<ValidatedTag>& f1, const ValidatedNumericType& c2) {
        VectorFunctionModel<ValidatedTag> r=f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=f1[i]*c2; } return r; }
    friend inline VectorFunctionModel<ValidatedTag> operator*(const ValidatedNumericType& c1, const VectorFunctionModel<ValidatedTag>& f2) {
        VectorFunctionModel<ValidatedTag> r=f2; for(SizeType i=0; i!=r.size(); ++i) { r[i]=c1*f2[i]; } return r; }

    friend inline VectorFunctionModel<ValidatedTag> operator+(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunction<ValidatedTag>& f2) {
        return f1+f1.create(f2); }
    friend inline VectorFunctionModel<ValidatedTag> operator-(const VectorFunctionModel<ValidatedTag>& f1, const VectorFunction<ValidatedTag>& f2) {
        return f1-f1.create(f2); }
    friend inline VectorFunctionModel<ValidatedTag> operator+(const VectorFunction<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
        return f2.create(f1)+f2; }
    friend inline VectorFunctionModel<ValidatedTag> operator-(const VectorFunction<ValidatedTag>& f1, const VectorFunctionModel<ValidatedTag>& f2) {
        return f2.create(f1)-f2; }

};

template<class P> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionModel<P>& f) {
    return os <<  f.operator VectorFunction<P>(); }

inline NormType norm(const VectorFunctionModel<ValidatedTag>& f) {
    return f._ptr->_norm(); }
inline VectorFunctionModel<ValidatedTag> embed(const ExactBoxType& d1, const VectorFunctionModel<ValidatedTag>& f, const ExactBoxType& d2) {
    return f._ptr->_embed(d1,d2); }
inline VectorFunctionModel<ValidatedTag> embed(const VectorFunctionModel<ValidatedTag>& f, const ExactBoxType& d) {
    return embed(ExactBoxType(),f,d); }
inline VectorFunctionModel<ValidatedTag> embed(const VectorFunctionModel<ValidatedTag>& f, const ExactIntervalType& d) {
    return embed(f,ExactBoxType(1,d)); }
inline VectorFunctionModel<ValidatedTag> restrict(const VectorFunctionModel<ValidatedTag>& f, const ExactBoxType& d) {
    VectorFunctionModelInterface<ValidatedTag>* rptr=f._ptr->_clone(); rptr->restrict(d); return rptr; }

inline ValidatedNumericType evaluate(const ScalarFunctionModel<ValidatedTag>& f, const Vector<ValidatedNumericType>& x) { return f(x); }
//inline Vector<ValidatedNumericType> evaluate(const VectorFunctionModel<ValidatedTag>& f, const Vector<ValidatedNumericType>& x) {
//    std::cerr<<"evaluate(const VectorFunctionModel<ValidatedTag>& f, const Vector<ValidatedNumericType>& x)\n"; return f._ptr->evaluate(x); }

inline ValidatedNumericType unchecked_evaluate(const ScalarFunctionModel<ValidatedTag>& f, const Vector<ValidatedNumericType>& x) { return f._ptr->_unchecked_evaluate(x); }
inline Vector<ValidatedNumericType> unchecked_evaluate(const VectorFunctionModel<ValidatedTag>& f, const Vector<ValidatedNumericType>& x) { return f._ptr->_unchecked_evaluate(x); }

inline ValidatedNumericType unchecked_evaluate(const ValidatedScalarFunction& f, const Vector<ValidatedNumericType>& x) {
    ScalarFunctionModelInterface<ValidatedTag> const* fptr = dynamic_cast<ScalarFunctionModelInterface<ValidatedTag> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(ScalarFunctionModel<ValidatedTag>(*fptr),x); } else { return evaluate(f,x); } }
inline Vector<ValidatedNumericType> unchecked_evaluate(const ValidatedVectorFunction& f, const Vector<ValidatedNumericType>& x) {
    VectorFunctionModelInterface<ValidatedTag> const* fptr = dynamic_cast<VectorFunctionModelInterface<ValidatedTag> const*>(f.raw_pointer());
    if(fptr) { return unchecked_evaluate(VectorFunctionModel<ValidatedTag>(*fptr),x); } else { return evaluate(f,x); } }
inline ScalarFunctionModel<ValidatedTag> unchecked_compose(const ValidatedScalarFunction& f, const VectorFunctionModel<ValidatedTag>& g) {
    ScalarFunctionModelInterface<ValidatedTag> const* fptr = dynamic_cast<ScalarFunctionModelInterface<ValidatedTag> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(ScalarFunctionModel<ValidatedTag>(*fptr),g); } else { return compose(f,g); } }
inline VectorFunctionModel<ValidatedTag> unchecked_compose(const ValidatedVectorFunction& f, const VectorFunctionModel<ValidatedTag>& g) {
    VectorFunctionModelInterface<ValidatedTag> const* fptr = dynamic_cast<VectorFunctionModelInterface<ValidatedTag> const*>(f.raw_pointer());
    if(fptr) { return unchecked_compose(VectorFunctionModel<ValidatedTag>(*fptr),g); } else { return compose(f,g); } }

//inline VectorFunctionModel<ValidatedTag> compose(const VectorFunctionModel<ValidatedTag>& f, const VectorFunctionModel<ValidatedTag>& g) { return g._ptr->_compose(f); }

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
    VectorFunctionModel<ValidatedTag> r=+f1; for(SizeType i=0; i!=r.size(); ++i) { r[i]=intersection(f1[i],f2[i]); } return r; }

inline VectorFunctionModel<ValidatedTag> antiderivative(const VectorFunctionModel<ValidatedTag>& f, SizeType j) {
    VectorFunctionModel<ValidatedTag> r(f);
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j); }
    return r;
}

inline VectorFunctionModel<ValidatedTag> antiderivative(const VectorFunctionModel<ValidatedTag>& f, SizeType j, const ValidatedNumericType& c) {
    VectorFunctionModel<ValidatedTag> r(f);
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=antiderivative(f[i],j,c); }
    return r;
}

inline VectorFunctionModel<ValidatedTag> partial_evaluate(const VectorFunctionModel<ValidatedTag>& f, SizeType j, const ValidatedNumericType& c) {
    return f._ptr->_partial_evaluate(j,c); }


template<class X> VectorFunctionModelElement<X>::operator const ScalarFunctionModel<X> () const {
    return _p->get(_i); }


inline Vector<ScalarFunctionModel<ValidatedTag>> ScalarFunctionModel<ValidatedTag>::create_coordinates(DomainType const& dom) const {
    Vector<ScalarFunctionModel<ValidatedTag>> res(dom.dimension(),this->_ptr->_create_zero(dom));
    for(SizeType i=0; i!=dom.dimension(); ++i) { res[i]=this->_ptr->_create_coordinate(dom,i); }
    return res;
}


inline VectorFunctionModel<ValidatedTag> ScalarFunctionModel<ValidatedTag>::create_identity() const { return this->_ptr->_create_identity(); }
inline ScalarFunctionModel<ValidatedTag> ScalarFunctionModel<ValidatedTag>::create(const ValidatedScalarFunction& g) const {
    ScalarFunctionModel<ValidatedTag> const& f=*this; return compose(g,f.create_identity()); }




// ExactTag output
template<class T> struct Representation { const T* pointer; Representation(const T& t) : pointer(&t) { } const T& reference() const { return *pointer; } };
template<class T> inline Representation<T> representation(const T& t) { return Representation<T>(t); }
template<class T> inline OutputStream& operator<<(OutputStream& os, const Representation<T>& obj) { obj.reference().repr(os); return os; }

template<class X> class FunctionModelFactoryInterface;

template<> class FunctionModelFactoryInterface<ValidatedTag>
{
    typedef ExactBoxType DomainType;
  public:
    virtual FunctionModelFactoryInterface<ValidatedTag>* clone() const = 0;
    virtual Void write(OutputStream& os) const = 0;
    inline ScalarFunctionModel<ValidatedTag> create(const ExactBoxType& domain, const ScalarFunctionInterface<ValidatedTag>& function) const;
    inline VectorFunctionModel<ValidatedTag> create(const ExactBoxType& domain, const VectorFunctionInterface<ValidatedTag>& function) const;
    inline ScalarFunctionModel<ValidatedTag> create_zero(const ExactBoxType& domain) const;
    inline VectorFunctionModel<ValidatedTag> create_zeros(SizeType result_size, const ExactBoxType& domain) const;
    inline ScalarFunctionModel<ValidatedTag> create_constant(const ExactBoxType& domain, const ValidatedNumber& value) const;
    inline ScalarFunctionModel<ValidatedTag> create_constant(const ExactBoxType& domain, const ValidatedNumericType& value) const;
    inline VectorFunctionModel<ValidatedTag> create_constants(const ExactBoxType& domain, const Vector<ValidatedNumericType>& values) const;
    inline ScalarFunctionModel<ValidatedTag> create_coordinate(const ExactBoxType& domain, SizeType index) const;
    inline ScalarFunctionModel<ValidatedTag> create_identity(const ExactIntervalType& domain) const;
    inline VectorFunctionModel<ValidatedTag> create_identity(const ExactBoxType& domain) const;
    inline CanonicalNumericType<ValidatedTag> create_number(const ValidatedNumber& number) const;
  private:
    virtual ScalarFunctionModelInterface<ValidatedTag>* _create(const ExactBoxType& domain, const ScalarFunctionInterface<ValidatedTag>& function) const = 0;
    virtual VectorFunctionModelInterface<ValidatedTag>* _create(const ExactBoxType& domain, const VectorFunctionInterface<ValidatedTag>& function) const = 0;
};

inline OutputStream& operator<<(OutputStream& os, const FunctionModelFactoryInterface<ValidatedTag>& factory) {
    factory.write(os); return os;
}

inline ScalarFunctionModel<ValidatedTag>
FunctionModelFactoryInterface<ValidatedTag>::create(const ExactBoxType& domain,
                                                const ScalarFunctionInterface<ValidatedTag>& function) const
{
    return this->_create(domain,function);
}

VectorFunctionModel<ValidatedTag>
FunctionModelFactoryInterface<ValidatedTag>::create(const ExactBoxType& domain,
                                                const VectorFunctionInterface<ValidatedTag>& function) const
{
    return this->_create(domain,function);
}


} // namespace Ariadne

#include "function/function.h"

namespace Ariadne {

ScalarFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_zero(const ExactBoxType& domain) const {
    return this->_create(domain,EffectiveScalarFunction::zero(domain.size())); }
VectorFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_zeros(SizeType result_size, const ExactBoxType& domain) const {
    return this->_create(domain,EffectiveVectorFunction::zeros(result_size,domain.size())); }
ScalarFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_constant(const ExactBoxType& domain, const ValidatedNumber& value) const {
        ValidatedNumericType concrete_value(value,Precision64()); return this->create_constant(domain,concrete_value); }
ScalarFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_constant(const ExactBoxType& domain, const ValidatedNumericType& value) const {
    return ScalarFunctionModel<ValidatedTag>(this->_create(domain,EffectiveScalarFunction::zero(domain.size())))+value; };
VectorFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_constants(const ExactBoxType& domain, const Vector<ValidatedNumericType>& values) const {
    return VectorFunctionModel<ValidatedTag>(this->_create(domain,EffectiveVectorFunction::zeros(values.size(),domain.size())))+values; };
ScalarFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_coordinate(const ExactBoxType& domain, SizeType index) const {
    return ScalarFunctionModel<ValidatedTag>(this->_create(domain,EffectiveScalarFunction::coordinate(domain.size(),index))); };
ScalarFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_identity(const ExactIntervalType& domain) const {
    return this->_create(ExactBoxType(1,domain),EffectiveScalarFunction::coordinate(1,0)); };
VectorFunctionModel<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_identity(const ExactBoxType& domain) const {
    return this->_create(domain,EffectiveVectorFunction::identity(domain.size())); };

inline CanonicalNumericType<ValidatedTag> FunctionModelFactoryInterface<ValidatedTag>::create_number(const ValidatedNumber& number) const {
    return CanonicalNumericType<ValidatedTag>(number,typename CanonicalNumericType<ValidatedTag>::PrecisionType()); }


} // namespace Ariadne

#endif
