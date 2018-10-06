/***************************************************************************
 *            algebra.hpp
 *
 *  Copyright 2010-17  Pieter Collins
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

/*! \file algebra.hpp
 *  \brief Algebras and normed algebras.
 */
#ifndef ARIADNE_ALGEBRA_HPP
#define ARIADNE_ALGEBRA_HPP

#include <iosfwd>
#include <iostream>
#include "../numeric/numeric.hpp"
#include "../utility/pointer.hpp"
#include "../numeric/operators.hpp"
#include "../algebra/algebra_interface.hpp"

namespace Ariadne {

template<class R, class A> R dynamic_handle_cast(A a) {
    typedef decltype(*declval<R>()._ptr) RP;
    return R(std::dynamic_pointer_cast<RP>(a.pointer()));
}

//! \brief Generic class for elements of unital algebras.
template<class X> class Algebra
    : public DeclareMixedArithmeticOperators<Algebra<X>,Int>
{
  private:
  public:
    std::shared_ptr< AlgebraInterface<X> > _ptr;
    std::shared_ptr< AlgebraInterface<X> > pointer() const { return _ptr; }
  public:
    typedef X ScalarType;
    typedef typename X::Paradigm Paradigm;
    typedef typename X::NumericType NumericType;
  public:
    Algebra() : _ptr() { }
    explicit Algebra(AlgebraInterface<X>* p) : _ptr(p) { }
    explicit Algebra(std::shared_ptr< AlgebraInterface<X> > p) : _ptr(p) { }
    Algebra(const AlgebraInterface<X>& a) : _ptr(a._create_copy()) { }
    Algebra(const Algebra<X>& a) : _ptr(a._ptr->_create_copy()) { }
    template<class A> A extract() const { A const* ap=dynamic_cast<A const*>(this->_ptr.operator->()); assert(ap); return *ap; }
    Algebra<X>& operator=(const X& c) { return *this = this->create_constant(c); }
    Algebra<X>& operator=(const Algebra<X>& a) { this->_ptr=std::shared_ptr< AlgebraInterface<X> >(a._ptr->_create_copy()); return *this; }
    operator const AlgebraInterface<X>& () const { return *this->_ptr; }
    Algebra<X> create() const { return Algebra<X>(this->_ptr->_create_zero()); }
    Algebra<X> clone() const { return Algebra<X>(this->_ptr->_create_copy()); }
    Algebra<X> create_zero() const { return Algebra<X>(this->_ptr->_create_zero()); }
    Algebra<X> create_constant(X const& c) const { return Algebra<X>(this->_ptr->_create_constant(c)); }
    OutputStream& write(OutputStream& os) const { return _ptr->write(os); }
  public:
    Void iadd(const X& c) { _ptr->_iadd(c); }
    Void imul(const X& c) { _ptr->_imul(c); }
    Void isma(const X& c, const Algebra<X>& a) { _ptr->_isma(c,*a._ptr); }
    Void ifma(const Algebra<X>& a1, const Algebra<X>& a2) { _ptr->_ifma(*a1._ptr,*a2._ptr); }

    friend Algebra<X> pos(Algebra<X> const& a) { return Algebra<X>(a); }
    friend Algebra<X> neg(Algebra<X> const& a) { return Algebra<X>(a._ptr->_neg()); }
    friend Algebra<X> add(Algebra<X> const& a1, Algebra<X> const& a2) { return Algebra<X>(a1._ptr->_add(*a2._ptr)); }
    friend Algebra<X> sub(Algebra<X> const& a1, Algebra<X> const& a2) { return Algebra<X>(a1._ptr->_sub(*a2._ptr)); }
    friend Algebra<X> mul(Algebra<X> const& a1, Algebra<X> const& a2) { return Algebra<X>(a1._ptr->_mul(*a2._ptr)); }
    friend Algebra<X> add(Algebra<X> const& a1, X const& c2) { return Algebra<X>(a1._ptr->_add(c2)); }
    friend Algebra<X> sub(Algebra<X> const& a1, X const& c2) { return Algebra<X>(a1._ptr->_sub(c2)); }
    friend Algebra<X> mul(Algebra<X> const& a1, X const& c2) { return Algebra<X>(a1._ptr->_mul(c2)); }
    friend Algebra<X> div(Algebra<X> const& a1, X const& c2) { return Algebra<X>(a1._ptr->_div(c2)); }
    friend Algebra<X> add(X const& c1, Algebra<X> const& a2) { return Algebra<X>(a2._ptr->_radd(c1)); }
    friend Algebra<X> sub(X const& c1, Algebra<X> const& a2) { return Algebra<X>(a2._ptr->_rsub(c1)); }
    friend Algebra<X> mul(X const& c1, Algebra<X> const& a2) { return Algebra<X>(a2._ptr->_rmul(c1)); }

    friend Algebra<X> operator+(Algebra<X> const& a) { return pos(a); }
    friend Algebra<X> operator-(Algebra<X> const& a) { return neg(a); }
    friend Algebra<X> operator+(Algebra<X> const& a1, Algebra<X> const& a2) { return add(a1,a2); }
    friend Algebra<X> operator-(Algebra<X> const& a1, Algebra<X> const& a2) { return sub(a1,a2); }
    friend Algebra<X> operator*(Algebra<X> const& a1, Algebra<X> const& a2) { return mul(a1,a2); }
    friend Algebra<X> operator+(Algebra<X> const& a1, X const& c2) { return add(a1,c2); }
    friend Algebra<X> operator-(Algebra<X> const& a1, X const& c2) { return sub(a1,c2); }
    friend Algebra<X> operator*(Algebra<X> const& a1, X const& c2) { return mul(a1,c2); }
    friend Algebra<X> operator/(Algebra<X> const& a1, X const& c2) { return div(a1,c2); }
    friend Algebra<X> operator+(X const& c1, Algebra<X> const& a2) { return add(c1,a2); }
    friend Algebra<X> operator-(X const& c1, Algebra<X> const& a2) { return sub(c1,a2); }
    friend Algebra<X> operator*(X const& c1, Algebra<X> const& a2) { return mul(c1,a2); }
    friend Algebra<X>& operator+=(Algebra<X>& a1, X const& c2) { return a1=add(a1,c2); }
    friend Algebra<X>& operator-=(Algebra<X>& a1, X const& c2) { return a1=sub(a1,c2); }
    friend Algebra<X>& operator*=(Algebra<X>& a1, X const& c2) { return a1=mul(a1,c2); }
    friend Algebra<X>& operator/=(Algebra<X>& a1, X const& c2) { return a1=div(a1,c2); }
    friend Algebra<X>& operator+=(Algebra<X>& a1, Algebra<X> const& a2) { return a1=add(a1,a2); }
    friend Algebra<X>& operator-=(Algebra<X>& a1, Algebra<X> const& a2) { return a1=sub(a1,a2); }
    friend Algebra<X>& operator*=(Algebra<X>& a1, Algebra<X> const& a2) { return a1=mul(a1,a2); }
    friend Algebra<X> sqr(Algebra<X> const& a) { return mul(a,a); }
    friend Algebra<X> pow(Algebra<X> const& a, Nat m) { return Algebra<X>(a._ptr->_pow(m)); }

    // DEPRECATED
    friend Algebra<X> operator+(Algebra<X> const& a1, Int const& c2) { return add(a1,X(c2)); }

    friend Algebra<X> operator/(Algebra<X> const& a1, Algebra<X> const& a2) { return mul(a1,rec(a2)); }
    friend Algebra<X> operator/(X const& c1, Algebra<X> const& a2) { return mul(c1,rec(a2)); }
    friend Algebra<X>& operator/=(Algebra<X>& a1, Algebra<X> const& a2) { return a1=div(a1,a2); }
    friend Algebra<X> pow(Algebra<X> const& a, Int n) {
        return n>=0 ? pow(a,Nat(n)) : rec(pow(a,Nat(-n))); }

    template<class OP> friend Algebra<X> apply(OP op, Algebra<X> a1, Algebra<X> a2) { return Algebra<X>(a1._ptr->_apply(op,*a2._ptr)); }
    template<class OP> friend Algebra<X> apply(OP op, Algebra<X> a1, X c2) { return Algebra<X>(a1._ptr->_apply(op,c2)); }
    template<class OP> friend Algebra<X> apply(OP op, Algebra<X> a) {
        TranscendentalAlgebraInterface<X> const* tap=dynamic_cast<TranscendentalAlgebraInterface<X>const*>(a._ptr.operator->());
        if(!tap) { ARIADNE_THROW(std::runtime_error,"apply(OP,Algebra<X>)","a="<<a<<" does not support transcendental operations."); }
        return Algebra<X>(tap->_apply(op)); }
    friend Algebra<X> apply(Neg op, Algebra<X> a) { return Algebra<X>(a._ptr->_apply(op)); }
};

template<class X> Algebra<X> rec(const Algebra<X>& a) { return apply(Rec(),a); }
template<class X> Algebra<X> sqrt(const Algebra<X>& a) { return apply(Sqrt(),a); }
template<class X> Algebra<X> exp(const Algebra<X>& a) { return apply(Exp(),a); }
template<class X> Algebra<X> log(const Algebra<X>& a) { return apply(Log(),a); }
template<class X> Algebra<X> sin(const Algebra<X>& a) { return apply(Sin(),a); }
template<class X> Algebra<X> cos(const Algebra<X>& a) { return apply(Cos(),a); }
template<class X> Algebra<X> tan(const Algebra<X>& a) { return apply(Tan(),a); }
template<class X> Algebra<X> atan(const Algebra<X>& a) { return apply(Atan(),a); }

// FIXME: Eliminate use of clone with std::dynamic_pointer_cast
template<class X, class OP> Algebra<X> apply(OP op, const Algebra<X>& a) {
    std::shared_ptr<TranscendentalAlgebraInterface<X> const> tap
        = std::dynamic_pointer_cast<TranscendentalAlgebraInterface<X>const>(a._ptr);
    if(!tap) { ARIADNE_FAIL_MSG("Cannot apply operator "<<op<<" to "<<a<<"\n"); }
    return Algebra<X>(tap->_apply(op));
}



template<class X> class NormedAlgebra
{
    static_assert(not IsSame<X,Real>::value,"");
  private:
    std::shared_ptr< NormedAlgebraInterface<X> > _ptr;
  public:
    typedef X ScalarType;
    typedef typename X::NumericType NumericType;
    NormedAlgebra(const Algebra<X>& a) : _ptr(std::dynamic_pointer_cast(a._ptr)) { }
    NormedAlgebra(NormedAlgebraInterface<X>* p) : _ptr(p) { }
    NormedAlgebra(std::shared_ptr< NormedAlgebraInterface<X> > p) : _ptr(p) { }
    NormedAlgebra(const NormedAlgebraInterface<X>& a) : _ptr(a.clone()) { }
    operator Algebra<X> () const { return _ptr->_create_copy(); }
    operator const NormedAlgebraInterface<X>& () const { return *_ptr; }
    NormedAlgebra<X>& operator=(Int c) { *this = this->create(); this->iadd(c); return *this; }
    NormedAlgebra<X>& operator=(const X& c) { *this = this->create(); this->iadd(c); return *this; }
    NormedAlgebra<X>& operator=(const NormedAlgebra<X>& a) { this->_ptr=std::shared_ptr< NormedAlgebraInterface<X> >(a._ptr->_create_copy()); return *this; }
    //! \brief The norm of the element.
    typename AlgebraTraits<X>::NormType norm() const { return _ptr->norm(); }
    //! \brief A value \a c minimising |a-c|.
    typename AlgebraTraits<X>::ValueType average() const { return _ptr->average(); }
    //! \brief The tolerance used to determine the truncation error when applying an analytic function.
    FloatDP tolerance() const { return _ptr->tolerance(); }
    //! \brief A value \c r such that \c |a-c1|<=r.
    typename AlgebraTraits<X>::NormType radius() const { return _ptr->radius(); }
    NormedAlgebra<X> create_zero() const { return NormedAlgebra<X>(_ptr->_create_zero()); }
    NormedAlgebra<X> create_constant(X c) const { return NormedAlgebra<X>(_ptr->_create_constant(c)); }
    NormedAlgebra<X> create_ball(FloatDPError r) const { return NormedAlgebra<X>(_ptr->_create_ball(r)); }
    NormedAlgebra<X> create_ball(FloatDPApproximation r) const {
        return create_ball(FloatDPError(r.raw())); }
    NormedAlgebra<X> create() const { return NormedAlgebra<X>(_ptr->_create_zero()); }
    NormedAlgebra<X> clone() const { return NormedAlgebra<X>(_ptr->_create_copy()); }
    Void clear() { this->imul(0); }
    Void iadd(const X& c) { _ptr->_iadd(c); }
    Void imul(const X& c) { _ptr->_imul(c); }
    Void isma(const X& c, const NormedAlgebra<X>& x) { _ptr->_isma(c,*x._ptr); }
    Void ifma(const NormedAlgebra<X>& x1, const NormedAlgebra<X>& x2) { _ptr->_ifma(*x1._ptr,*x2._ptr); }
    OutputStream& write(OutputStream& os) const { return _ptr->write(os); }
};

//! \brief Generic class for elements of unital algebras.
template<class X> class GradedAlgebra
{
  private:
  public:
    std::shared_ptr< GradedAlgebraInterface<X> > _ptr;
  public:
    typedef X ScalarType;
    typedef typename X::NumericType NumericType;
    GradedAlgebra() : _ptr() { }
    GradedAlgebra(GradedAlgebraInterface<X>* p) : _ptr(p) { }
    GradedAlgebra(const GradedAlgebraInterface<X>& a) : _ptr(a._create_copy()) { }
    GradedAlgebra(const GradedAlgebra<X>& a) : _ptr(a._ptr->_create_copy()) { }
    GradedAlgebra<X>& operator=(Int c) { *this = this->create(); this->iadd(c); return *this; }
    GradedAlgebra<X>& operator=(const X& c) { *this = this->create(); this->iadd(c); return *this; }
    GradedAlgebra<X>& operator=(const Algebra<X>& a) { this->_ptr=std::shared_ptr< AlgebraInterface<X> >(a._ptr->_create_copy()); return *this; }
    operator Algebra<X> () const { return _ptr->_create_copy(); }
    operator const GradedAlgebraInterface<X>& () const { return *_ptr; }
    GradedAlgebra<X> create() const { return GradedAlgebra<X>(_ptr->_create_zero()); }
    GradedAlgebra<X> clone() const { return GradedAlgebra<X>(_ptr->_create_copy()); }
    //! \brief The degree of the algebra.
    Nat degree() const { return _ptr->degree(); }
    //! \brief The value in the null grade.
    const X& value() const { return _ptr->value(); }
    OutputStream& write(OutputStream& os) const { return _ptr->write(os); }
  public:
    Void iadd(const X& c) { _ptr->_iadd(c); }
    Void imul(const X& c) { _ptr->_imul(c); }
    Void isma(const X& c, const GradedAlgebra<X>& x) { _ptr->_isma(c,*x._ptr); }
    Void ifma(const GradedAlgebra<X>& x1, const GradedAlgebra<X>& x2) { _ptr->_ifma(*x1._ptr,*x2._ptr); }
};

//! \brief Generic class for elements of unital algebras.
template<class X> class SymbolicAlgebra
{
  private:
  public:
    std::shared_ptr< SymbolicAlgebraInterface<X> > _ptr;
  public:
    typedef X ScalarType;
    typedef typename X::NumericType NumericType;
    SymbolicAlgebra() : _ptr() { }
    explicit SymbolicAlgebra(SymbolicAlgebraInterface<X>* p) : _ptr(p) { }
    SymbolicAlgebra(const SymbolicAlgebraInterface<X>& a) : _ptr(a._create_copy()) { }
    SymbolicAlgebra(const SymbolicAlgebra<X>& a) : _ptr(a._ptr->_create_copy()) { }
    //! \brief Create the representation of the operator \a op applied to \a a.
    SymbolicAlgebra(OperatorCode op, const SymbolicAlgebra<X>& a) : _ptr(a._ptr->_apply(op)) { }
    SymbolicAlgebra<X>& operator=(Int c) { *this = this->create(); this->iadd(c); return *this; }
    SymbolicAlgebra<X>& operator=(const X& c) { *this = this->create(); this->iadd(c); return *this; }
    SymbolicAlgebra<X>& operator=(const Algebra<X>& a) { this->_ptr=std::shared_ptr< AlgebraInterface<X> >(a._ptr->_create_copy()); return *this; }
    operator Algebra<X> () const { return _ptr->_create_copy(); }
    operator const SymbolicAlgebraInterface<X>& () const { return *_ptr; }
    SymbolicAlgebra<X> create() const { return SymbolicAlgebra<X>(_ptr->_create_zero()); }
    SymbolicAlgebra<X> create_zero() const { return SymbolicAlgebra<X>(_ptr->_create_zero()); }
    SymbolicAlgebra<X> clone() const { return SymbolicAlgebra<X>(_ptr->_create_copy()); }
    OutputStream& write(OutputStream& os) const { return _ptr->write(os); }
  public:
    Void iadd(const X& c) { _ptr->_iadd(c); }
    Void imul(const X& c) { _ptr->_imul(c); }
    Void isma(const X& c, const SymbolicAlgebra<X>& x) { _ptr->_isma(c,*x._ptr); }
    Void ifma(const SymbolicAlgebra<X>& x1, const SymbolicAlgebra<X>& x2) { _ptr->_ifma(*x1._ptr,*x2._ptr); }
};

template<class X> OutputStream& operator<<(OutputStream& os, const Algebra<X>& x) { return x.write(os); }
template<class X> OutputStream& operator<<(OutputStream& os, const NormedAlgebra<X>& x) { return x.write(os); }
} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_HPP */
