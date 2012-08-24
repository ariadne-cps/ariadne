/***************************************************************************
 *            function.h
 *
 *  Copyright 2008-9  Pieter Collins
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

namespace Ariadne {

typedef uint Nat;
typedef int Int;

template<class X> class Vector;
template<class X> class Matrix;

template<class X> class Differential;
template<class X> class Vector< Differential<X> >;

template<class X> class Formula;

template<class X> class ScalarFunction;
typedef ScalarFunction<Float> FloatScalarFunction;
typedef ScalarFunction<Interval> IntervalScalarFunction;
typedef ScalarFunction<Real> RealScalarFunction;

template<class X> class VectorFunction;
typedef VectorFunction<Float> FloatVectorFunction;
typedef VectorFunction<Interval> IntervalVectorFunction;
typedef VectorFunction<Real> RealVectorFunction;

template<class X> class VectorFunctionElementReference;

//! \ingroup FunctionModule
//! \brief A generic scalar function which can be evaluated over the number type \a X,  \f$f:\X^n\rightarrow\X\f$.
template<class X>
class ScalarFunction
{
  private:
    std::shared_ptr< const ScalarFunctionInterface<X> > _ptr;
  public:
    static ScalarFunction<X> zero(Nat m);
    static ScalarFunction<X> constant(Nat m, X c);
    static ScalarFunction<X> coordinate(Nat m, Nat j);
    static List< ScalarFunction<X> > coordinates(Nat n);

    explicit ScalarFunction(Nat as);
    explicit ScalarFunction(Nat as, Formula<X> f);

    ScalarFunction();
    ScalarFunction(ScalarFunctionInterface<X>* p) : _ptr(p) { }
    ScalarFunction(const ScalarFunctionInterface<X>& t) : _ptr(t._clone()) { }
    ScalarFunction(std::shared_ptr< const ScalarFunctionInterface<X> > p) : _ptr(p) { }
    ScalarFunction<X>& operator=(const ScalarFunctionInterface<X>& f) { _ptr=std::shared_ptr< const ScalarFunctionInterface<X> >(f._clone()); return *this; }

    template<class XX> ScalarFunction(const ScalarFunction<XX>& f, typename EnableIf<IsSafelyConvertible<XX,X>,Void>::Type* = 0)
        : _ptr(std::dynamic_pointer_cast< const ScalarFunctionInterface<X> >(f.managed_pointer())) { }
    template<class XX> inline ScalarFunction(const VectorFunctionElementReference<XX>& vfe,
                                             typename EnableIf<IsSafelyConvertible<XX,X>,Void>::Type* = 0);

    std::shared_ptr< const ScalarFunctionInterface<X> > managed_pointer() const  { return _ptr; }
    const ScalarFunctionInterface<X>* raw_pointer() const  { return _ptr.operator->(); }
    const ScalarFunctionInterface<X>& reference() const  { return _ptr.operator*(); }
    operator const ScalarFunctionInterface<X>& () const { return _ptr.operator*(); }

    Nat argument_size() const { return this->reference().argument_size(); }
    template<class XX> XX evaluate(const Vector<XX>& x) const { return this->reference().evaluate(x); }
    template<class XX> XX operator() (const Vector<XX>& x) const { return this->reference().evaluate(x); }

    ScalarFunction<X> derivative(Nat j) const { return this->reference().derivative(j); }

    template<class XX> Vector<XX> gradient(const Vector<XX>& x) const { return this->reference().gradient(x); }
    template<class XX> Differential<XX> differential(const Vector<XX>& x, Nat d) const { return this->_ptr->differential(x,d); }

    OutputStream& write(OutputStream& os) const { return this->_ptr->write(os); }
};

template<class X> template<class XX> inline
ScalarFunction<X>::ScalarFunction(const VectorFunctionElementReference<XX>& vfe,
                                  typename EnableIf<IsSafelyConvertible<XX,X>,Void>::Type*)
    : _ptr(vfe._vf.raw_pointer()->_get(vfe._i))
{
}

template<class X> inline OutputStream& operator<<(OutputStream& os, const ScalarFunction<X>& f) { return f.write(os); }
template<class X, class XX> inline XX evaluate(const ScalarFunction<X>& f, const Vector<XX>& x) { return f(x); }
template<class X, class XX> inline Vector<XX> gradient(const ScalarFunction<X>& f, const Vector<XX>& x) { return f.gradient(x); }

//template<class X> inline ScalarFunction<X> ScalarFunctionInterface<X>::clone() const { return this->_clone(); }
//template<class X> inline X ScalarFunctionInterface<X>::operator() (const Vector<X>& x) const { return this->evaluate(x); }
//template<class X> inline X ScalarFunctionInterface<X>::derivative(Nat j) const { return this->_derivative(j); }
inline ScalarFunction<Float> ScalarFunctionInterface<Float>::clone() const { return this->_clone(); }
inline ScalarFunction<Interval> ScalarFunctionInterface<Interval>::clone() const { return this->_clone(); }
inline ScalarFunction<Real> ScalarFunctionInterface<Real>::clone() const { return this->_clone(); }
inline Float ScalarFunctionInterface<Float>::operator() (const Vector<Float>& x) const { return this->evaluate(x); }
inline Interval ScalarFunctionInterface<Interval>::operator() (const Vector<Interval>& x) const { return this->evaluate(x); }
inline Real ScalarFunctionInterface<Real>::operator() (const Vector<Real>& x) const { return this->evaluate(x); }
inline ScalarFunction<Float> ScalarFunctionInterface<Float>::derivative(Nat j) const { return this->_derivative(j); }
inline ScalarFunction<Interval> ScalarFunctionInterface<Interval>::derivative(Nat j) const { return this->_derivative(j); }
inline ScalarFunction<Real> ScalarFunctionInterface<Real>::derivative(Nat j) const { return this->_derivative(j); }

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

ScalarFunction<Real> operator+(const ScalarFunction<Real>&);
ScalarFunction<Real> operator-(const ScalarFunction<Real>&);
ScalarFunction<Real> operator+(const ScalarFunction<Real>&, const ScalarFunction<Real>&);
ScalarFunction<Real> operator-(const ScalarFunction<Real>&, const ScalarFunction<Real>&);
ScalarFunction<Real> operator*(const ScalarFunction<Real>&, const ScalarFunction<Real>&);
ScalarFunction<Real> operator/(const ScalarFunction<Real>&, const ScalarFunction<Real>&);
ScalarFunction<Real> operator+(const ScalarFunction<Real>&, const Real&);
ScalarFunction<Real> operator-(const ScalarFunction<Real>&, const Real&);
ScalarFunction<Real> operator*(const ScalarFunction<Real>&, const Real&);
ScalarFunction<Real> operator/(const ScalarFunction<Real>&, const Real&);
ScalarFunction<Real> operator+(const Real&, const ScalarFunction<Real>&);
ScalarFunction<Real> operator-(const Real&, const ScalarFunction<Real>&);
ScalarFunction<Real> operator*(const Real&, const ScalarFunction<Real>&);
ScalarFunction<Real> operator/(const Real&, const ScalarFunction<Real>&);

ScalarFunction<Real> pow(const ScalarFunction<Real>&, Int);
ScalarFunction<Real> neg(const ScalarFunction<Real>&);
ScalarFunction<Real> rec(const ScalarFunction<Real>&);
ScalarFunction<Real> sqr(const ScalarFunction<Real>&);
ScalarFunction<Real> sqrt(const ScalarFunction<Real>&);
ScalarFunction<Real> exp(const ScalarFunction<Real>&);
ScalarFunction<Real> log(const ScalarFunction<Real>&);
ScalarFunction<Real> sin(const ScalarFunction<Real>&);
ScalarFunction<Real> cos(const ScalarFunction<Real>&);
ScalarFunction<Real> tan(const ScalarFunction<Real>&);




//! \ingroup FunctionModule
//! \brief A generic vector function which can be evaluated over the number type \a X,  \f$f:\X^n\rightarrow\X^m\f$.
template<class X>
class VectorFunction
{
  public:
    static VectorFunction<X> zeros(Nat rs, Nat as);
    static VectorFunction<X> identity(Nat n);

    VectorFunction();
    VectorFunction(Nat rs, Nat as);
    VectorFunction(Nat as, const List< Formula<X> >& flst);
    VectorFunction(Nat rs, ScalarFunction<X> sf);

    VectorFunction(VectorFunctionInterface<X>* fptr) : _ptr(fptr) { }
    VectorFunction(std::shared_ptr< VectorFunctionInterface<X> > fptr) : _ptr(fptr) { }
    VectorFunction(const VectorFunctionInterface<X>& fref) : _ptr(fref._clone()) { }
    std::shared_ptr< const VectorFunctionInterface<X> > managed_pointer() const { return this->_ptr; }
    const VectorFunctionInterface<X>* raw_pointer() const { return this->_ptr.operator->(); }
    const VectorFunctionInterface<X>& reference() const { return this->_ptr.operator*(); }
    operator const VectorFunctionInterface<X>& () const { return *this->_ptr; }

    VectorFunction(const List< ScalarFunction<X> >& lsf);
    VectorFunction(std::initializer_list< ScalarFunction<X> > lsf);
    template<class XX> VectorFunction(const List< ScalarFunction<XX> >& lsf, typename EnableIf< IsSafelyConvertible<XX,X>, Void >::Type* = 0) {
        *this=VectorFunction<X>(List< ScalarFunction<X> >(lsf)); }
    template<class XX> VectorFunction(const VectorFunction<XX>& vf, typename EnableIf< IsSafelyConvertible<XX,X>, Void >::Type* = 0)
        : _ptr(std::dynamic_pointer_cast< const VectorFunctionInterface<X> >(vf.managed_pointer())) { }

    ScalarFunction<X> get(Nat i) const { return this->_ptr->_get(i); }
    //Void set(Nat i, ScalarFunction<X> f) { this->_ptr->_set(i,f); };
    Void set(Nat i, ScalarFunction<X> f);

    ScalarFunction<X> operator[](Nat i) const { return this->get(i); }
    VectorFunctionElementReference<X> operator[](Nat i);

    Nat result_size() const { return this->_ptr->result_size(); }
    Nat argument_size() const { return this->_ptr->argument_size(); }

    template<class XX> Vector<XX> evaluate(const Vector<XX>& x) const { return this->_ptr->evaluate(x); }
    template<class XX> Vector<XX> operator() (const Vector<XX>& x) const { return this->_ptr->evaluate(x); }

    template<class XX> Matrix<XX> jacobian(const Vector<XX>& x) const { return this->_ptr->jacobian(x); }
    template<class XX> Vector< Differential<XX> > differentials(const Vector<XX>& x, Nat d) const { return this->_ptr->differentials(x,d); }

    OutputStream& write(OutputStream& os) const { return this->_ptr->write(os); }
  private:
    std::shared_ptr< const VectorFunctionInterface<X> > _ptr;
};

template<class X> inline OutputStream& operator<<(OutputStream& os, const VectorFunction<X>& f) { return f.write(os); }

template<class X, class XX> inline Vector<XX> evaluate(const VectorFunction<X>& f, const Vector<XX>& x) { return f(x); }
template<class X, class XX> inline Matrix<XX> jacobian(const VectorFunction<X>& f, const Vector<XX>& x) { return f.jacobian(x); }

/*
template<class X> VectorFunction<X> operator*(const ScalarFunction<X>& sf, const Vector<X>& e);
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

VectorFunction<Real> operator*(const ScalarFunction<Real>& sf, const Vector<Real>& e);
VectorFunction<Real> operator+(const VectorFunction<Real>& f1, const VectorFunction<Real>& f2);
VectorFunction<Real> operator-(const VectorFunction<Real>& f1, const VectorFunction<Real>& f2);
VectorFunction<Real> operator*(const VectorFunction<Real>& vf, const ScalarFunction<Real>& sf);
VectorFunction<Real> operator*(const ScalarFunction<Real>& sf, const VectorFunction<Real>& vf);

ScalarFunction<Real> embed(Nat as1, const ScalarFunction<Real>& f2, Nat as3);
VectorFunction<Real> embed(Nat as1, const VectorFunction<Real>& f2, Nat as3);

VectorFunction<Real> join(const ScalarFunction<Real>& f1, const ScalarFunction<Real>& f2);
VectorFunction<Real> join(const ScalarFunction<Real>& f1, const VectorFunction<Real>& f2);
VectorFunction<Real> join(const VectorFunction<Real>& f1, const ScalarFunction<Real>& f2);
VectorFunction<Real> join(const VectorFunction<Real>& f1, const VectorFunction<Real>& f2);

ScalarFunction<Real> compose(const ScalarFunction<Real>& f, const VectorFunction<Real>& g);
VectorFunction<Real> compose(const VectorFunction<Real>& f, const VectorFunction<Real>& g);

ScalarFunction<Real> lie_derivative(const ScalarFunction<Real>& g, const VectorFunction<Real>& f);

Formula<Real> formula(const ScalarFunction<Real>& f);
Vector< Formula<Real> > formula(const VectorFunction<Real>& f);


ScalarFunction<Interval> operator-(const ScalarFunction<Interval>&, const ScalarFunction<Interval>&);
ScalarFunction<Interval> operator-(const ScalarFunction<Interval>&, const Interval&);
ScalarFunction<Interval> operator-(const Interval&, const ScalarFunction<Interval>&);
VectorFunction<Interval> operator-(const VectorFunction<Interval>&, const VectorFunction<Interval>&);
VectorFunction<Interval> join(const VectorFunction<Interval>& f1, const VectorFunction<Interval>& f2);
ScalarFunction<Interval> compose(const ScalarFunction<Interval>& f, const VectorFunction<Interval>& g);
VectorFunction<Interval> compose(const VectorFunction<Interval>& f, const VectorFunction<Interval>& g);


template<class X>
struct VectorFunctionElementReference {
    VectorFunction<X>& _vf; Nat _i;
    VectorFunctionElementReference<X>(VectorFunction<X>& vf, Nat i) : _vf(vf), _i(i) { }
    void operator=(const ScalarFunction<X>& sf);
    VectorFunctionElementReference<X>& operator=(const VectorFunctionElementReference<X>& sfr);
    template<class XX> XX evaluate(const Vector<XX> & x) const;
    template<class XX> XX operator()(const Vector<XX> & x) const;
};

template<class X> inline VectorFunctionElementReference<X> VectorFunction<X>::operator[](Nat i) { return VectorFunctionElementReference<X>(*this,i); }
template<class X> inline OutputStream& operator<<(OutputStream& os, const VectorFunctionElementReference<X>& vfe) {
    return  os << static_cast< ScalarFunction<X> >(vfe); }

template<class X> inline void VectorFunctionElementReference<X>::operator=(const ScalarFunction<X>& sf) { _vf.set(_i,sf); }
template<class X> inline VectorFunctionElementReference<X>& VectorFunctionElementReference<X>::operator=(const VectorFunctionElementReference<X>& sfr) {
    _vf.set(_i,static_cast< ScalarFunction<X> >(sfr)); return *this; }
template<class X> template<class XX> inline XX VectorFunctionElementReference<X>::evaluate(const Vector<XX> & x) const {
    return static_cast< ScalarFunction<X> >(*this).evaluate(x); }
template<class X> template<class XX> inline XX VectorFunctionElementReference<X>::operator()(const Vector<XX> & x) const {
    return static_cast< ScalarFunction<X> >(*this).evaluate(x); }

inline ScalarFunction<Float> VectorFunctionInterface<Float>::operator[](Nat i) const { return this->_get(i); }
inline ScalarFunction<Interval> VectorFunctionInterface<Interval>::operator[](Nat i) const { return this->_get(i); }
inline ScalarFunction<Real> VectorFunctionInterface<Real>::operator[](Nat i) const { return this->_get(i); }


inline List< ScalarFunction<Real> > operator,(const Real& c1, const ScalarFunction<Real>& sf2) {
    return (ScalarFunction<Real>::constant(sf2.argument_size(),c1),sf2); }
inline List< ScalarFunction<Real> > operator,(const ScalarFunction<Real>& sf1, const Real& c2) {
    return (sf1,ScalarFunction<Real>::constant(sf1.argument_size(),c2)); }
inline List< ScalarFunction<Real> > operator,(const List< ScalarFunction<Real> >& vf1, const Real& c2) {
    return (vf1,ScalarFunction<Real>::constant(vf1.back().argument_size(),c2)); }




template<class X> class FunctionFactory;
typedef FunctionFactory<Interval> IntervalFunctionFactory;

template<>
class FunctionFactory<Interval>
{
    std::shared_ptr< const FunctionFactoryInterface<Interval> > _ptr;
  public:
    FunctionFactory(const FunctionFactoryInterface<Interval>& ref) : _ptr(ref.clone()) { }
    FunctionFactory(const FunctionFactoryInterface<Interval>* ptr) : _ptr(ptr) { }
    FunctionFactory(std::shared_ptr< const FunctionFactoryInterface<Interval> > ptr) : _ptr(ptr) { }
    inline ScalarFunction<Interval> create(const IntervalVector& d, const ScalarFunctionInterface<Interval>& f) const;
    inline VectorFunction<Interval> create(const IntervalVector& d, const VectorFunctionInterface<Interval>& f) const;
    inline ScalarFunction<Interval> create_zero(const IntervalVector& d) const;
    inline VectorFunction<Interval> create_identity(const IntervalVector& d) const;
    friend OutputStream& operator<<(OutputStream& os, const FunctionFactory<Interval>& factory);
};

inline ScalarFunction<Interval> FunctionFactoryInterface<Interval>::create(const IntervalVector& domain, const ScalarFunctionInterface<Interval>& function) const {
    return this->_create(domain,function); }
inline VectorFunction<Interval> FunctionFactoryInterface<Interval>::create(const IntervalVector& domain, const VectorFunctionInterface<Interval>& function) const {
    return this->_create(domain,function); }
inline ScalarFunction<Interval> FunctionFactoryInterface<Interval>::create_zero(const IntervalVector& domain) const {
    return this->_create(domain,RealScalarFunction::zero(domain.size())); }
inline VectorFunction<Interval> FunctionFactoryInterface<Interval>::create_identity(const IntervalVector& domain) const {
    return this->_create(domain,RealVectorFunction::identity(domain.size())); }

inline ScalarFunction<Interval> FunctionFactory<Interval>::create(const IntervalVector& domain, const ScalarFunctionInterface<Interval>& function) const {
    return this->_ptr->create(domain,function); }
inline VectorFunction<Interval> FunctionFactory<Interval>::create(const IntervalVector& domain, const VectorFunctionInterface<Interval>& function) const {
    return this->_ptr->create(domain,function); }
inline ScalarFunction<Interval> FunctionFactory<Interval>::create_zero(const IntervalVector& domain) const {
    return this->_ptr->create_zero(domain); }
inline VectorFunction<Interval> FunctionFactory<Interval>::create_identity(const IntervalVector& domain) const {
    return this->_ptr->create_identity(domain); }

inline OutputStream& operator<<(OutputStream& os, const FunctionFactory<Interval>& factory) {
    factory._ptr->write(os); return os; }

} // namespace Ariadne

#endif
