/***************************************************************************
 *            function.h
 *
 *  Copyright 2008-12  Pieter Collins
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

/*! \file function.h
 *  \brief Built-in and user functions and expressions
 */

#ifndef ARIADNE_FUNCTION_H
#define ARIADNE_FUNCTION_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>

#include "function_interface.h"

#include "macros.h"
#include "pointer.h"
#include "container.h"
#include "metaprogramming.h"

#include "numeric.h"
#include "vector.h"

#include "box.h"

namespace Ariadne {

template<class X> class VectorFunctionElementReference;

//! \ingroup FunctionModule
//! \brief A generic scalar function which can be evaluated over the number type \a X,  \f$f:\X^n\rightarrow\X\f$.
template<class P>
class ScalarFunction
{
    typedef P X;
  private:
    std::shared_ptr< const ScalarFunctionInterface<P> > _ptr;
  public:
    typedef P InformationTag;
    typedef X NumericType;

    static ScalarFunction<P> zero(Nat m);
    static ScalarFunction<P> constant(Nat m, NumericType c);
    static ScalarFunction<P> coordinate(Nat m, Nat j);
    static List< ScalarFunction<P> > coordinates(Nat n);

    explicit ScalarFunction(Nat as);
    explicit ScalarFunction(Nat as, Formula<NumericType> f);

    ScalarFunction();
    ScalarFunction(ScalarFunctionInterface<P>* p) : _ptr(p) { }
    ScalarFunction(const ScalarFunctionInterface<P>& t) : _ptr(t._clone()) { }
    ScalarFunction(std::shared_ptr< const ScalarFunctionInterface<P> > p) : _ptr(p) { }
    ScalarFunction<P>& operator=(const ScalarFunctionInterface<P>& f) {
        _ptr=std::shared_ptr< const ScalarFunctionInterface<P> >(f._clone()); return *this; }

    template<class PP> ScalarFunction(const ScalarFunction<PP>& f, typename EnableIf<IsSafelyConvertible<PP,P>,Void>::Type* = 0)
        : _ptr(std::dynamic_pointer_cast< const ScalarFunctionInterface<P> >(f.managed_pointer())) { }
    template<class PP> inline ScalarFunction(const VectorFunctionElementReference<PP>& vfe,
                                             typename EnableIf<IsSafelyConvertible<PP,P>,Void>::Type* = 0);

    std::shared_ptr< const ScalarFunctionInterface<P> > managed_pointer() const  { return _ptr; }
    const ScalarFunctionInterface<P>* raw_pointer() const  { return _ptr.operator->(); }
    const ScalarFunctionInterface<P>& reference() const  { return _ptr.operator*(); }
    operator const ScalarFunctionInterface<P>& () const { return _ptr.operator*(); }

    Nat argument_size() const { return this->reference().argument_size(); }
    template<class XX> XX evaluate(const Vector<XX>& x) const { return this->reference().evaluate(x); }
    template<class XX> XX operator() (const Vector<XX>& x) const { return this->reference().evaluate(x); }

    ScalarFunction<P> derivative(Nat j) const { return this->reference().derivative(j); }

    template<class XX> Vector<XX> gradient(const Vector<XX>& x) const { return this->reference().gradient(x); }
    template<class XX> Differential<XX> differential(const Vector<XX>& x, Nat d) const { return this->_ptr->differential(x,d); }

    OutputStream& write(OutputStream& os) const { return this->_ptr->write(os); }
};

template<class P> template<class PP> inline
ScalarFunction<P>::ScalarFunction(const VectorFunctionElementReference<PP>& vfe,
                                  typename EnableIf<IsSafelyConvertible<PP,P>,Void>::Type*)
    : _ptr(vfe._vf.raw_pointer()->_get(vfe._i))
{
}

template<class P> inline OutputStream& operator<<(OutputStream& os, const ScalarFunction<P>& f) { return f.write(os); }
template<class P, class XX> inline XX evaluate(const ScalarFunction<P>& f, const Vector<XX>& x) { return f(x); }
template<class P, class XX> inline Vector<XX> gradient(const ScalarFunction<P>& f, const Vector<XX>& x) { return f.gradient(x); }

ApproximateScalarFunction ScalarFunctionInterface<ApproximateNumberType>::clone() const { return this->_clone(); }
ValidatedScalarFunction ScalarFunctionInterface<ValidatedNumberType>::clone() const { return this->_clone(); }
EffectiveScalarFunction ScalarFunctionInterface<EffectiveNumberType>::clone() const { return this->_clone(); }
inline ValidatedNumberType ScalarFunctionInterface<ValidatedNumberType>::evaluate(const Vector<ExactNumberType>& x) const {
    return this->evaluate(Vector<ValidatedNumberType>(x)); }
inline ApproximateNumberType ScalarFunctionInterface<ApproximateNumberType>::operator() (const Vector<ApproximateNumberType>& x) const {
    return this->evaluate(x); }
inline ValidatedNumberType ScalarFunctionInterface<ValidatedNumberType>::operator() (const Vector<ValidatedNumberType>& x) const {
    return this->evaluate(x); }
inline Interval ScalarFunctionInterface<ValidatedNumberType>::operator() (const Vector<Interval>& x) const {
    return this->evaluate(x); }
inline EffectiveNumberType ScalarFunctionInterface<EffectiveNumberType>::operator() (const Vector<EffectiveNumberType>& x) const {
    return this->evaluate(x); }
inline ApproximateScalarFunction ScalarFunctionInterface<ApproximateNumberType>::derivative(Nat j) const {
    return this->_derivative(j); }
inline ValidatedScalarFunction ScalarFunctionInterface<ValidatedNumberType>::derivative(Nat j) const {
    return this->_derivative(j); }
inline EffectiveScalarFunction ScalarFunctionInterface<EffectiveNumberType>::derivative(Nat j) const {
    return this->_derivative(j); }


/*
template<class X> ScalarFunction<X> operator+(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator-(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator+(const ScalarFunction<X>&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator-(const ScalarFunction<X>&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator*(const ScalarFunction<X>&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator/(const ScalarFunction<X>&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator+(const ScalarFunction<X>&, const X&);
template<class X> ScalarFunction<X> operator-(const ScalarFunction<X>&, const X&);
template<class X> ScalarFunction<X> operator*(const ScalarFunction<X>&, const X&);
template<class X> ScalarFunction<X> operator/(const ScalarFunction<X>&, const X&);
template<class X> ScalarFunction<X> operator+(const X&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator-(const X&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator*(const X&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator/(const X&, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator+(const ScalarFunction<X>&, double);
template<class X> ScalarFunction<X> operator-(const ScalarFunction<X>&, double);
template<class X> ScalarFunction<X> operator*(const ScalarFunction<X>&, double);
template<class X> ScalarFunction<X> operator/(const ScalarFunction<X>&, double);
template<class X> ScalarFunction<X> operator+(double, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator-(double, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator*(double, const ScalarFunction<X>&);
template<class X> ScalarFunction<X> operator/(double, const ScalarFunction<X>&);

template<class X> ScalarFunction<X> pow(const ScalarFunction<X>&, Int);
template<class X> ScalarFunction<X> neg(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> rec(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> sqr(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> sqrt(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> exp(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> log(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> sin(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> cos(const ScalarFunction<X>&);
template<class X> ScalarFunction<X> tan(const ScalarFunction<X>&);
*/

EffectiveScalarFunction operator+(const EffectiveScalarFunction&);
EffectiveScalarFunction operator-(const EffectiveScalarFunction&);
EffectiveScalarFunction operator+(const EffectiveScalarFunction&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator-(const EffectiveScalarFunction&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator*(const EffectiveScalarFunction&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator/(const EffectiveScalarFunction&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator+(const EffectiveScalarFunction&, const EffectiveNumberType&);
EffectiveScalarFunction operator-(const EffectiveScalarFunction&, const EffectiveNumberType&);
EffectiveScalarFunction operator*(const EffectiveScalarFunction&, const EffectiveNumberType&);
EffectiveScalarFunction operator/(const EffectiveScalarFunction&, const EffectiveNumberType&);
EffectiveScalarFunction operator+(const EffectiveNumberType&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator-(const EffectiveNumberType&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator*(const EffectiveNumberType&, const EffectiveScalarFunction&);
EffectiveScalarFunction operator/(const EffectiveNumberType&, const EffectiveScalarFunction&);

EffectiveScalarFunction pow(const EffectiveScalarFunction&, Int);
EffectiveScalarFunction neg(const EffectiveScalarFunction&);
EffectiveScalarFunction rec(const EffectiveScalarFunction&);
EffectiveScalarFunction sqr(const EffectiveScalarFunction&);
EffectiveScalarFunction sqrt(const EffectiveScalarFunction&);
EffectiveScalarFunction exp(const EffectiveScalarFunction&);
EffectiveScalarFunction log(const EffectiveScalarFunction&);
EffectiveScalarFunction sin(const EffectiveScalarFunction&);
EffectiveScalarFunction cos(const EffectiveScalarFunction&);
EffectiveScalarFunction tan(const EffectiveScalarFunction&);




//! \ingroup FunctionModule
//! \brief A generic vector function which can be evaluated over the number type \a X,  \f$f:\X^n\rightarrow\X^m\f$.
template<class P>
class VectorFunction
{
    typedef P X;
  private:
    std::shared_ptr< const VectorFunctionInterface<X> > _ptr;
  public:
    typedef X NumericType;

    static VectorFunction<P> zeros(Nat rs, Nat as);
    static VectorFunction<P> identity(Nat n);

    VectorFunction();
    VectorFunction(Nat rs, Nat as);
    VectorFunction(Nat as, const List< Formula<X> >& flst);
    VectorFunction(Nat rs, ScalarFunction<P> sf);

    VectorFunction(VectorFunctionInterface<P>* fptr) : _ptr(fptr) { }
    VectorFunction(std::shared_ptr< VectorFunctionInterface<P> > fptr) : _ptr(fptr) { }
    VectorFunction(const VectorFunctionInterface<P>& fref) : _ptr(fref._clone()) { }
    std::shared_ptr< const VectorFunctionInterface<P> > managed_pointer() const { return this->_ptr; }
    const VectorFunctionInterface<P>* raw_pointer() const { return this->_ptr.operator->(); }
    const VectorFunctionInterface<P>& reference() const { return this->_ptr.operator*(); }
    operator const VectorFunctionInterface<P>& () const { return *this->_ptr; }

    VectorFunction(const List< ScalarFunction<P> >& lsf);
    VectorFunction(std::initializer_list< ScalarFunction<P> > lsf);
    template<class PP> VectorFunction(const List< ScalarFunction<PP> >& lsf, typename EnableIf< IsSafelyConvertible<PP,P>, Void >::Type* = 0) {
        *this=VectorFunction<P>(List< ScalarFunction<P> >(lsf)); }
    template<class PP> VectorFunction(const VectorFunction<PP>& vf, typename EnableIf< IsSafelyConvertible<PP,P>, Void >::Type* = 0)
        : _ptr(std::dynamic_pointer_cast< const VectorFunctionInterface<P> >(vf.managed_pointer())) { }

    ScalarFunction<P> get(Nat i) const { return this->_ptr->_get(i); }
    //Void set(Nat i, ScalarFunction<P> f) { this->_ptr->_set(i,f); };
    Void set(Nat i, ScalarFunction<P> f);

    ScalarFunction<P> operator[](Nat i) const { return this->get(i); }
    VectorFunctionElementReference<P> operator[](Nat i);

    Nat result_size() const { return this->_ptr->result_size(); }
    Nat argument_size() const { return this->_ptr->argument_size(); }

    template<class XX> auto evaluate(const Vector<XX>& x) const -> decltype(this->_ptr->evaluate(x)) { return this->_ptr->evaluate(x); }
    template<class XX> auto operator() (const Vector<XX>& x) const -> decltype(this->_ptr->evaluate(x)) { return this->_ptr->evaluate(x); }

    template<class XX> auto jacobian(const Vector<XX>& x) const -> decltype(this->_ptr->jacobian(x)) {
        return this->_ptr->jacobian(x); }
    template<class XX> auto differentials(const Vector<XX>& x, Nat d) const -> decltype(this->_ptr->differentials(x,d)) {
        return this->_ptr->differentials(x,d); }

    OutputStream& write(OutputStream& os) const { return this->_ptr->write(os); }
};

template<class P> inline OutputStream& operator<<(OutputStream& os, const VectorFunction<P>& f) { return f.write(os); }

template<class P, class XX> inline Vector<XX> evaluate(const VectorFunction<P>& f, const Vector<XX>& x) { return f(x); }
template<class P, class XX> inline Matrix<XX> jacobian(const VectorFunction<P>& f, const Vector<XX>& x) { return f.jacobian(x); }

/*
template<class X> VectorFunction<P> operator*(const ScalarFunction<X>& sf, const Vector<X>& e);
template<class X> VectorFunction<X> operator+(const VectorFunction<X>& f1, const VectorFunction<X>& f2);
template<class X> VectorFunction<X> operator-(const VectorFunction<X>& f1, const VectorFunction<X>& f2);
template<class X> VectorFunction<X> operator*(const VectorFunction<X>& vf, const ScalarFunction<X>& sf);
template<class X> VectorFunction<X> operator*(const ScalarFunction<X>& sf, const VectorFunction<X>& vf);

template<class X> ScalarFunction<X> embed(Nat as1, const ScalarFunction<X>& f2, Nat as3);
template<class X> VectorFunction<X> embed(Nat as1, const VectorFunction<X>& f2, Nat as3);

template<class X> VectorFunction<X> join(const ScalarFunction<X>& f1, const ScalarFunction<X>& f2);
template<class X> VectorFunction<X> join(const ScalarFunction<X>& f1, const VectorFunction<X>& f2);
template<class X> VectorFunction<X> join(const VectorFunction<X>& f1, const ScalarFunction<X>& f2);
template<class X> VectorFunction<X> join(const VectorFunction<X>& f1, const VectorFunction<X>& f2);

template<class X> ScalarFunction<X> compose(const ScalarFunction<X>& f, const VectorFunction<X>& g);
template<class X> VectorFunction<X> compose(const VectorFunction<X>& f, const VectorFunction<X>& g);

template<class X> ScalarFunction<X> lie_derivative(const ScalarFunction<X>& g, const VectorFunction<X>& f);
*/

EffectiveVectorFunction operator*(const EffectiveScalarFunction& sf, const Vector<EffectiveNumberType>& e);
EffectiveVectorFunction operator+(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2);
EffectiveVectorFunction operator-(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2);
EffectiveVectorFunction operator*(const EffectiveVectorFunction& vf, const EffectiveScalarFunction& sf);
EffectiveVectorFunction operator*(const EffectiveScalarFunction& sf, const EffectiveVectorFunction& vf);
EffectiveVectorFunction operator*(const EffectiveNumberType& c, const EffectiveVectorFunction& vf);

EffectiveScalarFunction embed(Nat as1, const EffectiveScalarFunction& f2, Nat as3);
EffectiveVectorFunction embed(Nat as1, const EffectiveVectorFunction& f2, Nat as3);

EffectiveVectorFunction join(const EffectiveScalarFunction& f1, const EffectiveScalarFunction& f2);
EffectiveVectorFunction join(const EffectiveScalarFunction& f1, const EffectiveVectorFunction& f2);
EffectiveVectorFunction join(const EffectiveVectorFunction& f1, const EffectiveScalarFunction& f2);
EffectiveVectorFunction join(const EffectiveVectorFunction& f1, const EffectiveVectorFunction& f2);

EffectiveScalarFunction compose(const EffectiveScalarFunction& f, const EffectiveVectorFunction& g);
EffectiveVectorFunction compose(const EffectiveVectorFunction& f, const EffectiveVectorFunction& g);

EffectiveScalarFunction lie_derivative(const EffectiveScalarFunction& g, const EffectiveVectorFunction& f);

Formula<EffectiveNumberType> formula(const EffectiveScalarFunction& f);
Vector< Formula<EffectiveNumberType> > formula(const EffectiveVectorFunction& f);

inline Vector<ValidatedNumberType> VectorFunctionInterface<ValidatedTag>::evaluate(const Vector<ExactNumberType> & x) const {
    return this->evaluate(Vector<ValidatedNumberType>(x)); }


ValidatedScalarFunction operator-(const ValidatedScalarFunction&, const ValidatedScalarFunction&);
ValidatedScalarFunction operator-(const ValidatedScalarFunction&, const ValidatedNumberType&);
ValidatedScalarFunction operator-(const ValidatedNumberType&, const ValidatedScalarFunction&);
ValidatedVectorFunction operator-(const ValidatedVectorFunction&, const ValidatedVectorFunction&);
ValidatedVectorFunction join(const ValidatedVectorFunction& f1, const ValidatedVectorFunction& f2);
ValidatedScalarFunction compose(const ValidatedScalarFunction& f, const ValidatedVectorFunction& g);
ValidatedVectorFunction compose(const ValidatedVectorFunction& f, const ValidatedVectorFunction& g);

template<class P>
struct VectorFunctionElementReference {
    VectorFunction<P>& _vf; Nat _i;
    VectorFunctionElementReference<P>(VectorFunction<P>& vf, Nat i) : _vf(vf), _i(i) { }
    void operator=(const ScalarFunction<P>& sf);
    VectorFunctionElementReference<P>& operator=(const VectorFunctionElementReference<P>& sfr);
    template<class XX> XX evaluate(const Vector<XX> & x) const;
    template<class XX> XX operator()(const Vector<XX> & x) const;
};

template<class P> inline VectorFunctionElementReference<P> VectorFunction<P>::operator[](Nat i) { return VectorFunctionElementReference<P>(*this,i); }
template<class P> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionElementReference<P>& vfe) {
    return  os << static_cast< ScalarFunction<P> >(vfe); }

template<class P> inline void VectorFunctionElementReference<P>::operator=(const ScalarFunction<P>& sf) { _vf.set(_i,sf); }
template<class P> inline VectorFunctionElementReference<P>& VectorFunctionElementReference<P>::operator=(const VectorFunctionElementReference<P>& sfr) {
    _vf.set(_i,static_cast< ScalarFunction<P> >(sfr)); return *this; }
template<class P> template<class XX> inline XX VectorFunctionElementReference<P>::evaluate(const Vector<XX> & x) const {
    return static_cast< ScalarFunction<P> >(*this).evaluate(x); }
template<class P> template<class XX> inline XX VectorFunctionElementReference<P>::operator()(const Vector<XX> & x) const {
    return static_cast< ScalarFunction<P> >(*this).evaluate(x); }

inline ApproximateScalarFunction VectorFunctionInterface<ApproximateNumberType>::operator[](Nat i) const { return this->_get(i); }
inline ValidatedScalarFunction VectorFunctionInterface<ValidatedNumberType>::operator[](Nat i) const { return this->_get(i); }
inline EffectiveScalarFunction VectorFunctionInterface<EffectiveNumberType>::operator[](Nat i) const { return this->_get(i); }


inline List< EffectiveScalarFunction > operator,(const EffectiveNumberType& c1, const EffectiveScalarFunction& sf2) {
    return (EffectiveScalarFunction::constant(sf2.argument_size(),c1),sf2); }
inline List< EffectiveScalarFunction > operator,(const EffectiveScalarFunction& sf1, const EffectiveNumberType& c2) {
    return (sf1,EffectiveScalarFunction::constant(sf1.argument_size(),c2)); }
inline List< EffectiveScalarFunction > operator,(const List< EffectiveScalarFunction >& vf1, const EffectiveNumberType& c2) {
    return (vf1,EffectiveScalarFunction::constant(vf1.back().argument_size(),c2)); }




template<class P> class FunctionFactory;
typedef FunctionFactory<ValidatedTag> ValidatedFunctionFactory;

template<>
class FunctionFactory<ValidatedTag>
{
    std::shared_ptr< const FunctionFactoryInterface<ValidatedTag> > _ptr;
  public:
    FunctionFactory(const FunctionFactoryInterface<ValidatedTag>& ref) : _ptr(ref.clone()) { }
    FunctionFactory(const FunctionFactoryInterface<ValidatedTag>* ptr) : _ptr(ptr) { }
    FunctionFactory(std::shared_ptr< const FunctionFactoryInterface<ValidatedTag> > ptr) : _ptr(ptr) { }
    inline ValidatedScalarFunction create(const Box& d, const ValidatedScalarFunctionInterface& f) const;
    inline ValidatedVectorFunction create(const Box& d, const ValidatedVectorFunctionInterface& f) const;
    inline ValidatedScalarFunction create_zero(const Box& d) const;
    inline ValidatedVectorFunction create_identity(const Box& d) const;
    friend OutputStream& operator<<(OutputStream& os, const FunctionFactory<ValidatedTag>& factory);
};

inline ValidatedScalarFunction FunctionFactoryInterface<ValidatedTag>::create(const Box& domain, const ValidatedScalarFunctionInterface& function) const {
    return this->_create(domain,function); }
inline ValidatedVectorFunction FunctionFactoryInterface<ValidatedTag>::create(const Box& domain, const ValidatedVectorFunctionInterface& function) const {
    return this->_create(domain,function); }
inline ValidatedScalarFunction FunctionFactoryInterface<ValidatedTag>::create_zero(const Box& domain) const {
    return this->_create(domain,EffectiveScalarFunction::zero(domain.size())); }
inline ValidatedVectorFunction FunctionFactoryInterface<ValidatedTag>::create_identity(const Box& domain) const {
    return this->_create(domain,EffectiveVectorFunction::identity(domain.size())); }

inline ValidatedScalarFunction FunctionFactory<ValidatedTag>::create(const Box& domain, const ValidatedScalarFunctionInterface& function) const {
    return this->_ptr->create(domain,function); }
inline ValidatedVectorFunction FunctionFactory<ValidatedTag>::create(const Box& domain, const ValidatedVectorFunctionInterface& function) const {
    return this->_ptr->create(domain,function); }
inline ValidatedScalarFunction FunctionFactory<ValidatedTag>::create_zero(const Box& domain) const {
    return this->_ptr->create_zero(domain); }
inline ValidatedVectorFunction FunctionFactory<ValidatedTag>::create_identity(const Box& domain) const {
    return this->_ptr->create_identity(domain); }

inline OutputStream& operator<<(OutputStream& os, const ValidatedFunctionFactory& factory) {
    factory._ptr->write(os); return os; }

} // namespace Ariadne

#endif
