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

#include "numeric.h"
#include "vector.h"
#include "matrix.h"

class D;
namespace Ariadne {

typedef uint Nat;
typedef int Int;

template<class X> class Vector;
template<class X> class Matrix;

template<class X> class Polynomial;
template<class X> class Formula;
template<class X> class Algebra;

template<class X> class ScalarFunction;
typedef ScalarFunction<Float> FloatScalarFunction;
typedef ScalarFunction<Interval> IntervalScalarFunction;
typedef ScalarFunction<Real> RealScalarFunction;

template<class X> class VectorFunction;
typedef VectorFunction<Float> FloatVectorFunction;
typedef VectorFunction<Interval> IntervalVectorFunction;
typedef VectorFunction<Real> RealVectorFunction;

//! \ingroup FunctionModule
//! \brief A scalar function which can only be evaluated approximately, \f$f:\F^n\rightarrow\F\f$.
template<>
class ScalarFunction<Float>
{
  public:
    explicit ScalarFunction<Float>(Nat n=0u);
    ScalarFunction<Float>(ScalarFunctionInterface<Float>* fptr) : _ptr(fptr) { }
    ScalarFunction<Float>(shared_ptr< ScalarFunctionInterface<Float> > fptr) : _ptr(fptr) { }
    ScalarFunction<Float>(const ScalarFunctionInterface<Float>& fref) : _ptr(fref._clone()) { }

    const ScalarFunctionInterface<Float>* raw_pointer() const { return this->_ptr.operator->(); }
    operator const ScalarFunctionInterface<Float>& () const { return *this->_ptr; }

    uint argument_size() const { return this->_ptr->argument_size(); }

    template<class X> X evaluate(const Vector<X>& x) const { return this->_ptr->evaluate(x); }
    Float operator() (const Vector<Float>& x) const { return this->_ptr->evaluate(x); }

    template<class X> Vector<X> gradient(const Vector<X>& x) const {
        return this->evaluate(Differential<X>::variables(1u,x)).gradient(); }

    std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }

  private:
    shared_ptr< ScalarFunctionInterface<Float> > _ptr;
};

//! \ingroup FunctionModule
//! \brief A scalar function which can only be evaluated approximately or on intervals, \f$f:\I^n\rightarrow\I\f$.
template<>
class ScalarFunction<Interval>
{
  public:
    explicit ScalarFunction<Interval>(Nat n=0u);
    ScalarFunction<Interval>(ScalarFunctionInterface<Interval>* fptr) : _ptr(fptr) { }
    ScalarFunction<Interval>(shared_ptr< ScalarFunctionInterface<Interval> > fptr) : _ptr(fptr) { }
    ScalarFunction<Interval>(const ScalarFunctionInterface<Interval>& fref) : _ptr(fref._clone()) { }
    ScalarFunction<Interval>(const ScalarFunction<Real>& f);

    const ScalarFunctionInterface<Interval>* raw_pointer() const { return this->_ptr.operator->(); }
    operator const ScalarFunctionInterface<Interval>& () const { return *this->_ptr; }
    operator ScalarFunction<Float> () const { return ScalarFunction<Float>(this->_ptr); }

    uint argument_size() const { return this->_ptr->argument_size(); }
    template<class X> X evaluate(const Vector<X>& x) const { return this->_ptr->evaluate(x); }
    Float operator() (const Vector<Float>& x) const { return this->_ptr->evaluate(x); }
    Interval operator() (const Vector<Interval>& x) const { return this->_ptr->evaluate(x); }

    template<class X> Vector<X> gradient(const Vector<X>& x) const {
        return this->evaluate(Differential<X>::variables(1u,x)).gradient(); }

    ScalarFunction<Interval> derivative(uint j) const;

    std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }

  private:
    shared_ptr< ScalarFunctionInterface<Interval> > _ptr;
};


//! \ingroup FunctionModule
//! \brief A scalar function \f$f:\R^n\rightarrow\R\f$.
template<>
class ScalarFunction<Real>
{
    typedef uint Nat;
    typedef uint SizeType;
    typedef std::ostream OStream;
  public:
    static ScalarFunction<Real> zero(Nat n);
    static ScalarFunction<Real> constant(Nat n, Real c);
    static ScalarFunction<Real> coordinate(Nat n, uint i);
    static List< ScalarFunction<Real> > coordinates(Nat n);

    explicit ScalarFunction<Real>(Nat n=0u);
    ScalarFunction<Real>(Nat as, const Formula<Real>& e);
    ScalarFunction<Real>(const Polynomial<Real>& p);

    ScalarFunction<Real>(ScalarFunctionInterface<Real>* fptr) : _ptr(fptr) { }
    const ScalarFunctionInterface<Real>* raw_pointer() const { return this->_ptr.operator->(); }
    const shared_ptr< ScalarFunctionInterface<Real> > shared_pointer() const { return this->_ptr; }
    operator const ScalarFunctionInterface<Real>& () const { return *this->_ptr; }

    operator ScalarFunction<Float> () const { return ScalarFunction<Float>(this->_ptr); }
    //operator ScalarFunction<Interval> () const { return ScalarFunction<Interval>(this->_ptr); }

    Nat argument_size() const { return this->_ptr->argument_size(); }
    template<class X> X evaluate(const Vector<X>& x) const { return this->_ptr->evaluate(x); }

    Float operator() (const Vector<Float>& x) const { return this->_ptr->evaluate(x); }
    Interval operator() (const Vector<Interval>& x) const { return this->_ptr->evaluate(x); }
    Real operator() (const Vector<Real>& x) const { return this->_ptr->evaluate(x); }

    template<class X> Vector<X> gradient(const Vector<X>& x) const {
        return this->evaluate(Differential<X>::variables(1u,x)).gradient(); }

    std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }

  public:
    ScalarFunction<Real> derivative(Nat j) const;
    Polynomial<Real> polynomial() const;
  public:
    friend ScalarFunction<Real> compose(const ScalarFunction<Real>&, const VectorFunction<Real>&);
    friend Real evaluate(const ScalarFunction<Real>&, const Vector<Real>&);
    friend ScalarFunction<Real> embed(const ScalarFunction<Real>&, uint i);
  public:
    friend ScalarFunction<Real> operator+(const ScalarFunction<Real>&);
    friend ScalarFunction<Real> operator-(const ScalarFunction<Real>&);
    friend ScalarFunction<Real> operator+(const ScalarFunction<Real>&, const ScalarFunction<Real>&);
    friend ScalarFunction<Real> operator-(const ScalarFunction<Real>&, const ScalarFunction<Real>&);
    friend ScalarFunction<Real> operator*(const ScalarFunction<Real>&, const ScalarFunction<Real>&);
    friend ScalarFunction<Real> operator/(const ScalarFunction<Real>&, const ScalarFunction<Real>&);
  protected:
    virtual void _check_type(const ScalarFunctionInterface<Real>*) { }
  private:
    shared_ptr< ScalarFunctionInterface<Real> > _ptr;
};


inline ScalarFunction<Interval>::ScalarFunction(const ScalarFunction<Real>& f)
    : _ptr(boost::dynamic_pointer_cast< ScalarFunctionInterface<Real>  >(f.shared_pointer())) { }

// FIXME: This conversion to a RealScalarFunction is a hack. The kind of derivative available should
// depend on the function type, and not require casting.
inline ScalarFunction<Interval> ScalarFunction<Interval>::derivative(uint j) const {
    return ScalarFunction<Interval>(static_cast<const ScalarFunction<Interval>&>(boost::dynamic_pointer_cast< ScalarFunctionInterface<Real> >(_ptr)->derivative(j))); }

ScalarFunction<Interval> restrict(ScalarFunction<Interval> const& f, const IntervalVector& dom);

inline Float evaluate_approx(const ScalarFunction<Real>& f, const Vector<Float>& x) { return f(x); }
inline Interval evaluate(const ScalarFunction<Real>& f, const Vector<Interval>& x) { return f(x); }
inline Real evaluate(const ScalarFunction<Real>& f, const Vector<Real>& x) { return f(x); }
inline std::ostream& operator<<(std::ostream& os, const ScalarFunction<Real>& f) { return f.write(os); }

inline FloatScalarFunction ScalarFunctionInterface<Float>::clone() const { return this->_clone(); }
inline IntervalScalarFunction ScalarFunctionInterface<Interval>::clone() const { return this->_clone(); }
inline RealScalarFunction ScalarFunctionInterface<Real>::clone() const { return this->_clone(); }

inline Float ScalarFunctionInterface<Float>::operator() (const Vector<Float>& x) const { return this->evaluate(x); }
inline Interval ScalarFunctionInterface<Interval>::operator() (const Vector<Interval>& x) const { return this->evaluate(x); }
inline Real ScalarFunctionInterface<Real>::operator() (const Vector<Real>& x) const { return this->evaluate(x); }

inline FloatScalarFunction ScalarFunctionInterface<Float>::derivative(uint j) const { return this->_derivative(j); }
inline IntervalScalarFunction ScalarFunctionInterface<Interval>::derivative(uint j) const { return this->_derivative(j); }
inline RealScalarFunction ScalarFunctionInterface<Real>::derivative(uint j) const { return this->_derivative(j); }

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

ScalarFunction<Real> pow(const ScalarFunction<Real>&, int);
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
//! \brief A vector function which can only be evaluated approximately, \f$f:\F^n\rightarrow\F^m\f$.
template<>
class VectorFunction<Float>
{
  public:
    VectorFunction<Float>(VectorFunctionInterface<Float>* fptr) : _ptr(fptr) { }
    VectorFunction<Float>(shared_ptr< VectorFunctionInterface<Float> > fptr) : _ptr(fptr) { }
    VectorFunction<Float>(const VectorFunctionInterface<Float>& fref) : _ptr(fref._clone()) { }
    const VectorFunctionInterface<Float>* raw_pointer() const { return this->_ptr.operator->(); }
    operator const VectorFunctionInterface<Float>& () const { return *this->_ptr; }

    uint result_size() const { return this->_ptr->result_size(); }
    uint argument_size() const { return this->_ptr->argument_size(); }

    template<class X> Vector<X> evaluate(const Vector<X>& x) const { return this->_ptr->evaluate(x); }
    Vector<Float> operator() (const Vector<Float>& x) const { return this->_ptr->evaluate(x); }

    template<class X> Matrix<X> jacobian(const Vector<X>& x) const {
        return this->_ptr->evaluate(Differential<X>::variables(1u,x)).jacobian(); }

    std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }

  private:
    shared_ptr< VectorFunctionInterface<Float> > _ptr;
};

//! \ingroup FunctionModule
//! \brief A scalar function which can only be evaluated approximately or on intervals, \f$f:\I^n\rightarrow\I\f$.
template<>
class VectorFunction<Interval>
{
    friend class VectorFunction<Float>;
  public:
    static VectorFunction<Interval> identity(Nat n);

    VectorFunction<Interval>();
    VectorFunction<Interval>(const List<ScalarFunction<Interval> >& lf);
    VectorFunction<Interval>(const List<ScalarFunction<Real> >& lf);
    VectorFunction<Interval>(Nat n, const ScalarFunction<Interval>& f);

    VectorFunction<Interval>(VectorFunctionInterface<Interval>* fptr) : _ptr(fptr) { }
    VectorFunction<Interval>(shared_ptr< VectorFunctionInterface<Interval> > fptr) : _ptr(fptr) { }
    VectorFunction<Interval>(const VectorFunctionInterface<Interval>& fref) : _ptr(fref._clone()) { }
    const VectorFunctionInterface<Interval>* raw_pointer() const { return this->_ptr.operator->(); }
    operator const VectorFunctionInterface<Interval>& () const { return *this->_ptr; }
    operator VectorFunction<Float> () const { return VectorFunction<Float>(this->_ptr); }

    uint result_size() const { return this->_ptr->result_size(); }
    uint argument_size() const { return this->_ptr->argument_size(); }

    template<class X> Vector<X> evaluate(const Vector<X>& x) const { return this->_ptr->evaluate(x); }
    Vector<Float> operator() (const Vector<Float>& x) const { return this->_ptr->evaluate(x); }
    Vector<Interval> operator() (const Vector<Interval>& x) const { return this->_ptr->evaluate(x); }

    template<class X> Matrix<X> jacobian(const Vector<X>& x) const {
        return this->_ptr->evaluate(Differential<X>::variables(1u,x)).jacobian(); }

    ScalarFunction<Interval> operator[](uint i) const { return _ptr->operator[](i); }

    std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }

  private:
    shared_ptr< VectorFunctionInterface<Interval> > _ptr;
};

ScalarFunction<Interval> operator+(const ScalarFunction<Interval>& f);
ScalarFunction<Interval> operator-(const ScalarFunction<Interval>& f);
ScalarFunction<Interval> operator+(const ScalarFunction<Interval>& f1, const ScalarFunction<Interval>& f2);
ScalarFunction<Interval> operator-(const ScalarFunction<Interval>& f1, const ScalarFunction<Interval>& f2);
ScalarFunction<Interval> operator*(const ScalarFunction<Interval>& f1, const ScalarFunction<Interval>& f2);
ScalarFunction<Interval> operator/(const ScalarFunction<Interval>& f1, const ScalarFunction<Interval>& f2);
ScalarFunction<Interval> operator+(const ScalarFunction<Interval>& f1, const Interval& c2);
ScalarFunction<Interval> operator+(const Interval& c1, const ScalarFunction<Interval>& f2);
ScalarFunction<Interval> operator-(const ScalarFunction<Interval>& f1, const Interval& c2);
ScalarFunction<Interval> operator-(const Interval& c1, const ScalarFunction<Interval>& f2);
ScalarFunction<Interval> operator*(const ScalarFunction<Interval>& f1, const Interval& c2);
ScalarFunction<Interval> operator*(const Interval& c1, const ScalarFunction<Interval>& f2);
ScalarFunction<Interval> operator/(const ScalarFunction<Interval>& f1, const Interval& c2);
ScalarFunction<Interval> operator/(const Interval& c1, const ScalarFunction<Interval>& f2);
ScalarFunction<Interval> pow(const ScalarFunction<Interval>& f, Int n);
ScalarFunction<Interval> neg(const ScalarFunction<Interval>& f);
ScalarFunction<Interval> rec(const ScalarFunction<Interval>& f);
ScalarFunction<Interval> sqr(const ScalarFunction<Interval>& f);
ScalarFunction<Interval> sqrt(const ScalarFunction<Interval>& f);
ScalarFunction<Interval> exp(const ScalarFunction<Interval>& f);
ScalarFunction<Interval> log(const ScalarFunction<Interval>& f);
ScalarFunction<Interval> sin(const ScalarFunction<Interval>& f);
ScalarFunction<Interval> cos(const ScalarFunction<Interval>& f);
ScalarFunction<Interval> tan(const ScalarFunction<Interval>& f);

VectorFunction<Interval> operator*(const ScalarFunction<Interval>& sf, const Vector<Interval>& e);
VectorFunction<Interval> operator+(const VectorFunction<Interval>& f1, const VectorFunction<Interval>& f2);
VectorFunction<Interval> operator-(const VectorFunction<Interval>& f1, const VectorFunction<Interval>& f2);
VectorFunction<Interval> operator*(const VectorFunction<Interval>& vf, const ScalarFunction<Interval>& sf);
VectorFunction<Interval> operator*(const ScalarFunction<Interval>& sf, const VectorFunction<Interval>& vf);

VectorFunction<Interval> join(const ScalarFunction<Interval>& f1, const ScalarFunction<Interval>& f2);
VectorFunction<Interval> join(const ScalarFunction<Interval>& f1, const VectorFunction<Interval>& f2);
VectorFunction<Interval> join(const VectorFunction<Interval>& f1, const ScalarFunction<Interval>& f2);
VectorFunction<Interval> join(const VectorFunction<Interval>& f1, const VectorFunction<Interval>& f2);

ScalarFunction<Interval> compose(const ScalarFunction<Interval>& f, const VectorFunction<Interval>& g);
VectorFunction<Interval> compose(const VectorFunction<Interval>& f, const VectorFunction<Interval>& g);
ScalarFunction<Interval> lie_derivative(const ScalarFunction<Interval>& g, const VectorFunction<Interval>& f);

inline Interval evaluate(const ScalarFunction<Interval>& f, const Vector<Interval>& x) { return f(x); }
inline Vector<Interval> evaluate(const VectorFunction<Interval>& f, const Vector<Interval>& x) { return f(x); }


template<class X>
struct VectorFunctionElementReference {
    VectorFunction<X>& _vf; uint _i;
    VectorFunctionElementReference<X>(VectorFunction<X>& vf, uint i) : _vf(vf), _i(i) { }
    operator ScalarFunction<X>() const;
    void operator=(const ScalarFunction<X>& sf);
    VectorFunctionElementReference<X>& operator=(const VectorFunctionElementReference<X>& sfr);
    template<class XX> XX evaluate(const Vector<XX> & x) const;
    template<class XX> XX operator()(const Vector<XX> & x) const;
};


//! \ingroup FunctionModule
//! \brief A vector function \f$f:\R^n\rightarrow\R^m\f$.
template<>
class VectorFunction<Real>
{
    friend class VectorFunction<Float>;
    friend class VectorFunction<Interval>;
    typedef uint SizeType;
    typedef ushort SmoothnessType;
    typedef std::ostream OStream;
  public:
    static VectorFunction<Real> constant(const Vector<Real>& c, Nat as);
    static VectorFunction<Real> identity(Nat n);

    VectorFunction<Real>();
    VectorFunction<Real>(Nat rs, Nat as);
    VectorFunction<Real>(Nat rs, const ScalarFunction<Real>& sf);
    VectorFunction<Real>(const List< ScalarFunction<Real> >& lsf);
    VectorFunction<Real>(Nat as, const List< Formula<Real> >& e);
    VectorFunction<Real>(const Vector< Polynomial<Real> >& p);

    VectorFunction<Real>(VectorFunctionInterface<Real>*);
    VectorFunction<Real>(shared_ptr< VectorFunctionInterface<Real> > fptr) : _ptr(fptr) { }

    const VectorFunctionInterface<Real>* raw_pointer() const { return this->_ptr.operator->(); }
    operator const VectorFunctionInterface<Real>& () const { return *this->_ptr; }
    operator const VectorFunction<Float> () const { return VectorFunction<Float>(this->_ptr); }
    operator const VectorFunction<Interval> () const { return VectorFunction<Interval>(this->_ptr); }

    VectorFunction<Real>& operator=(const VectorFunction<Real>& f) {
        this->_check_type(f._raw_pointer()); this->_ptr=f._ptr; return *this; }

    Nat result_size() const { return this->_ptr->result_size(); }
    Nat argument_size() const { return this->_ptr->argument_size(); }

    template<class X> Vector<X> evaluate(const Vector<X>& x) const { return this->_ptr->evaluate(x); }
    Vector<Float> operator()(const Vector<Float>& x) const { return this->evaluate(x); }
    Vector<Interval> operator()(const Vector<Interval>& x) const { return this->evaluate(x); };
    Vector<Real> operator()(const Vector<Real>& x) const { return this->evaluate(x); };

    template<class X> Matrix<X> jacobian(const Vector<X>& x) const {
        return this->_ptr->evaluate(Differential<X>::variables(1u,x)).jacobian(); }

    std::ostream& write(std::ostream& os) const { return this->_ptr->write(os); }

  public:
    ScalarFunction<Real> operator[](Nat i) const;
    VectorFunctionElementReference<Real> operator[](Nat i) { return VectorFunctionElementReference<Real>(*this,i); }

    ScalarFunction<Real> get(Nat) const;
    void set(Nat,ScalarFunction<Real>);

    VectorFunction<Real> polynomial() const;
  public:
    friend VectorFunction<Real> join(const ScalarFunction<Real>&, const ScalarFunction<Real>&);
    friend VectorFunction<Real> join(const ScalarFunction<Real>&, const VectorFunction<Real>&);
    friend VectorFunction<Real> join(const VectorFunction<Real>&, const ScalarFunction<Real>&);
    friend VectorFunction<Real> join(const VectorFunction<Real>&, const VectorFunction<Real>&);

    friend VectorFunction<Real> embed(const VectorFunction<Real>& f1, Nat n2);
    friend VectorFunction<Real> compose(const VectorFunction<Real>&, const VectorFunction<Real>&);
    friend Vector<Real> evaluate(const VectorFunction<Real>&, const Vector<Real>&);

    friend VectorFunction<Real> operator+(const VectorFunction<Real>&, const VectorFunction<Real>&);
    friend VectorFunction<Real> operator-(const VectorFunction<Real>&, const VectorFunction<Real>&);
    friend VectorFunction<Real> operator*(const VectorFunction<Real>&, const ScalarFunction<Real>&);
    friend VectorFunction<Real> operator*(const ScalarFunction<Real>&, const VectorFunction<Real>&);
  public:
    const VectorFunctionInterface<Real>* _raw_pointer() const { return this->_ptr.operator->(); }
  protected:
    virtual void _check_type(const VectorFunctionInterface<Real>*) { }
  private:
    shared_ptr< VectorFunctionInterface<Real> > _ptr;
};

template<> inline VectorFunctionElementReference<Real>::operator ScalarFunction<Real>() const { return _vf.get(_i); }
template<> inline void VectorFunctionElementReference<Real>::operator=(const ScalarFunction<Real>& sf) { _vf.set(_i,sf); }
template<> inline VectorFunctionElementReference<Real>& VectorFunctionElementReference<Real>::operator=(const VectorFunctionElementReference<Real>& sfr) {
    _vf.set(_i,static_cast< ScalarFunction<Real> >(sfr)); return *this; }
template<class X> template<class XX> inline XX VectorFunctionElementReference<X>::evaluate(const Vector<XX> & x) const {
    return static_cast< ScalarFunction<X> >(*this).evaluate(x); }
template<class X> template<class XX> inline XX VectorFunctionElementReference<X>::operator()(const Vector<XX> & x) const {
    return static_cast< ScalarFunction<X> >(*this).evaluate(x); }

inline Vector<Float> evaluate_approx(const VectorFunction<Real>& f, const Vector<Float>& x) { return f(x); }
inline Vector<Interval> evaluate(const VectorFunction<Real>& f, const Vector<Interval>& x) { return f(x); }
inline Vector<Real> evaluate(const VectorFunction<Real>& f, const Vector<Real>& x) { return f(x); }
inline Matrix<Float> jacobian_approx(const VectorFunction<Real>& f, const Vector<Float>& x);
inline Matrix<Interval> jacobian(const VectorFunction<Real>& f, const Vector<Interval>& x);


inline ScalarFunction<Float> VectorFunctionInterface<Float>::operator[](Nat i) const { return this->_get(i); }
inline ScalarFunction<Interval> VectorFunctionInterface<Interval>::operator[](Nat i) const { return this->_get(i); }
inline ScalarFunction<Real> VectorFunctionInterface<Real>::operator[](Nat i) const { return this->_get(i); }

VectorFunction<Real> operator*(const ScalarFunction<Real>& sf, const Vector<Real>& e);
VectorFunction<Real> operator+(const VectorFunction<Real>& f1, const VectorFunction<Real>& f2);
VectorFunction<Real> operator-(const VectorFunction<Real>& f1, const VectorFunction<Real>& f2);
VectorFunction<Real> operator*(const VectorFunction<Real>& vf, const ScalarFunction<Real>& sf);
VectorFunction<Real> operator*(const ScalarFunction<Real>& sf, const VectorFunction<Real>& vf);

VectorFunction<Real> join(const ScalarFunction<Real>&, const ScalarFunction<Real>&);
VectorFunction<Real> join(const ScalarFunction<Real>&, const VectorFunction<Real>&);
VectorFunction<Real> join(const VectorFunction<Real>&, const ScalarFunction<Real>&);
VectorFunction<Real> join(const VectorFunction<Real>&, const VectorFunction<Real>&);

ScalarFunction<Real> compose(const ScalarFunction<Real>& f, const VectorFunction<Real>& g);
VectorFunction<Real> compose(const VectorFunction<Real>& f, const VectorFunction<Real>& g);
ScalarFunction<Real> lie_derivative(const ScalarFunction<Real>& g, const VectorFunction<Real>& f);

inline List< ScalarFunction<Real> > operator,(const Real& c1, const ScalarFunction<Real>& sf2) {
    return (ScalarFunction<Real>::constant(sf2.argument_size(),c1),sf2); }
inline List< ScalarFunction<Real> > operator,(const ScalarFunction<Real>& sf1, const Real& c2) {
    return (sf1,ScalarFunction<Real>::constant(sf1.argument_size(),c2)); }
inline List< ScalarFunction<Real> > operator,(const ScalarFunction<Real>& sf1, double c2) {
    return (sf1,ScalarFunction<Real>::constant(sf1.argument_size(),c2)); }
inline List< ScalarFunction<Real> > operator,(const List< ScalarFunction<Real> >& vf1, const Real& c2) {
    return (vf1,ScalarFunction<Real>::constant(vf1.back().argument_size(),c2)); }

inline std::ostream& operator<<(std::ostream& os, const VectorFunction<Real>& f) { return f.write(os); }

Vector< Formula<Real> > formula(const VectorFunction<Real>& f);

//! \ingroup FunctionModule
//! \brief A scalar function of the form \f$f(x)=c\f$.
class ScalarConstantFunction
    : public ScalarFunction<Real>
{
  public:
    //! \brief Construct the constant function \f$f(x)=\sum a_ix_i+b\f$.
    ScalarConstantFunction(uint n, const Real& c);
  protected:
    virtual void _check_type(const ScalarFunctionInterface<Real>* ptr) const;
};

//! \ingroup FunctionModule
//! \brief A scalar function of the form \f$f(x)=c\f$.
class CoordinateFunction
    : public ScalarFunction<Real>
{
  public:
    //! \brief Construct the constant function \f$f(x)=\sum a_ix_i+b\f$.
    CoordinateFunction(uint n, uint j);
  protected:
    virtual void _check_type(const ScalarFunctionInterface<Real>* ptr) const;
};

//! \ingroup FunctionModule
//! \brief A scalar function of the form \f$f(x)=\sum_{i=1}^{n}a_ix_i+b\f$.
class ScalarAffineFunction
    : public ScalarFunction<Real>
{
  public:
    //! \brief Construct the affine function \f$f(x)=\sum a_ix_i+b\f$.
    ScalarAffineFunction(const Vector<Real>& a, const Real& b);
  protected:
    virtual void _check_type(const ScalarFunctionInterface<Real>* ptr) const;
};

//! \ingroup FunctionModule
//! \brief A scalar function of the form \f$f(x)=\sum_{\alpha} c_\alpha x_1^{\alpha_1}\cdots x_n^{\alpha_n}\f$.
class ScalarPolynomialFunction
    : public ScalarFunction<Real>
{
  public:
    //! \brief Construct the affine function \f$f(x)=\sum a_ix_i+b\f$.
    ScalarPolynomialFunction(const Polynomial<Real>& p);
  protected:
    virtual void _check_type(const ScalarFunctionInterface<Real>* ptr) const;
};



class VectorConstantFunction
    : public VectorFunction<Real>
{
  public:
    VectorConstantFunction(const Vector<Real>& c, uint as);
  protected:
    virtual void _check_type(const VectorFunctionInterface<Real>* ptr) const;
};


class VectorAffineFunction
    : public VectorFunction<Real>
{
  public:
    VectorAffineFunction(const Matrix<Real>& A, const Vector<Real>& b);
    const Matrix<Real> A() const;
    const Vector<Real> b() const;
  protected:
    virtual void _check_type(const VectorFunctionInterface<Real>* ptr) const;
};



//! \ingroup FunctionModule
//! \brief The identiy function on \f$\R^n\f$.
class IdentityFunction
    : public VectorFunction<Real>
{
  public:
    //! \brief Construct the identity function in dimension \a n.
    IdentityFunction(uint n);
  protected:
    virtual void _check_type(const VectorFunctionInterface<Real>* ptr) const;
};


class ProjectionFunction
    : public VectorFunction<Real>
{
  public:
    //! \brief Construct the identity function in dimension \a n.
    ProjectionFunction(uint n);
    //! \brief Construct the projection functions \f$f_i(x)=x_{i+k}\f$ for \f$i=0,\ldots,m-1\f$. Precondition: \f$m+k\leq n\f$.
    ProjectionFunction(uint m, uint n, uint k);
    //! \brief Construct the projection function  with \f$f_i(x)=x_{p_i}\f$ for \f$i=0,\ldots,m-1\f$.
    ProjectionFunction(uint m, uint n, const Array<uint>& p);
    //! \brief Construct the projection function with \f$f_i(x)=x_{p_i}\f$ Ffor \f$i=0,\ldots,|p|-1\f$.
    ProjectionFunction(const Array<uint>& p, uint n);

    const Array<uint>& p() const;
    const uint p(uint i) const;
  protected:
    virtual void _check_type(const VectorFunctionInterface<Real>* ptr) const;
};

template<class X> class FunctionFactory;
typedef FunctionFactory<Interval> IntervalFunctionFactory;

template<>
class FunctionFactory<Interval>
{
    shared_ptr< const FunctionFactoryInterface<Interval> > _ptr;
  public:
    FunctionFactory(const FunctionFactoryInterface<Interval>& ref) : _ptr(ref.clone()) { }
    FunctionFactory(const FunctionFactoryInterface<Interval>* ptr) : _ptr(ptr) { }
    FunctionFactory(shared_ptr< const FunctionFactoryInterface<Interval> > ptr) : _ptr(ptr) { }
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
