/***************************************************************************
 *            algebra.h
 *
 *  Copyright 2010-15  Pieter Collins
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

/*! \file algebra.h
 *  \brief Algebras and normed algebras.
 */
#ifndef ARIADNE_ALGEBRA_H
#define ARIADNE_ALGEBRA_H

#include <iosfwd>
#include <iostream>
#include "numeric/numeric.h"
#include "utility/pointer.h"
#include "expression/operators.h"
#include "algebra/algebra_interface.h"

namespace Ariadne {

//! \brief Generic class for elements of unital algebras.
template<class X> class Algebra
{
  private:
  public:
    std::shared_ptr< AlgebraInterface<X> > _ptr;
  public:
    typedef X ScalarType;
    typedef typename X::Paradigm Paradigm;
    typedef typename X::NumericType NumericType;
    Algebra() : _ptr() { }
    Algebra(AlgebraInterface<X>* p) : _ptr(p) { }
    Algebra(const AlgebraInterface<X>& a) : _ptr(a._clone()) { }
    Algebra(const Algebra<X>& a) : _ptr(a._ptr->_clone()) { }
    Algebra<X>& operator=(Int c) { *this = this->create(); this->iadd(c); return *this; }
    Algebra<X>& operator=(const X& c) { *this = this->create(); this->iadd(c); return *this; }
    Algebra<X>& operator=(const Algebra<X>& a) { this->_ptr=std::shared_ptr< AlgebraInterface<X> >(a._ptr->_clone()); return *this; }
    operator const AlgebraInterface<X>& () const { return *_ptr; }
    Algebra<X> create() const { return Algebra<X>(_ptr->_create()); }
    Algebra<X> clone() const { return Algebra<X>(_ptr->_clone()); }
    Algebra<X> create_zero() const { return Algebra<X>(_ptr->_create()); }
    OutputStream& write(OutputStream& os) const { return _ptr->write(os); }
  public:
    Void iadd(const X& c) { _ptr->_iadd(c); }
    Void imul(const X& c) { _ptr->_imul(c); }
    Void isma(const X& c, const Algebra<X>& x) { _ptr->_isma(c,*x._ptr); }
    Void ifma(const Algebra<X>& x1, const Algebra<X>& x2) { _ptr->_ifma(*x1._ptr,*x2._ptr); }
};

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
    operator Algebra<X> () const { return _ptr->_clone(); }
    operator const NormedAlgebraInterface<X>& () const { return *_ptr; }
    NormedAlgebra<X>& operator=(Int c) { *this = this->create(); this->iadd(c); return *this; }
    NormedAlgebra<X>& operator=(const X& c) { *this = this->create(); this->iadd(c); return *this; }
    NormedAlgebra<X>& operator=(const NormedAlgebra<X>& a) { this->_ptr=std::shared_ptr< NormedAlgebraInterface<X> >(a._ptr->_clone()); return *this; }
    //! \brief The norm of the element.
    typename AlgebraTraits<X>::NormType norm() const { return _ptr->norm(); }
    //! \brief A value \a c minimising |a-c|.
    typename AlgebraTraits<X>::ValueType average() const { return _ptr->average(); }
    //! \brief The tolerance used to determine the truncation error when applying an analytic function.
    Float64 tolerance() const { return _ptr->tolerance(); }
    //! \brief A value \c r such that \c |a-c1|<=r.
    typename AlgebraTraits<X>::NormType radius() const { return _ptr->radius(); }
    NormedAlgebra<X> create_zero() const { return NormedAlgebra<X>(_ptr->_create()); }
    NormedAlgebra<X> create_constant(X c) const { return NormedAlgebra<X>(_ptr->_create_constant(c)); }
    NormedAlgebra<X> create_ball(ErrorFloat64 r) const { return NormedAlgebra<X>(_ptr->_create_ball(r)); }
    NormedAlgebra<X> create_ball(ApproximateFloat64 r) const {
        return create_ball(ErrorFloat64(r.raw())); }
    NormedAlgebra<X> create() const { return NormedAlgebra<X>(_ptr->_create()); }
    NormedAlgebra<X> clone() const { return NormedAlgebra<X>(_ptr->_clone()); }
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
    GradedAlgebra(const GradedAlgebraInterface<X>& a) : _ptr(a._clone()) { }
    GradedAlgebra(const GradedAlgebra<X>& a) : _ptr(a._ptr->_clone()) { }
    GradedAlgebra<X>& operator=(Int c) { *this = this->create(); this->iadd(c); return *this; }
    GradedAlgebra<X>& operator=(const X& c) { *this = this->create(); this->iadd(c); return *this; }
    GradedAlgebra<X>& operator=(const Algebra<X>& a) { this->_ptr=std::shared_ptr< AlgebraInterface<X> >(a._ptr->_clone()); return *this; }
    operator Algebra<X> () const { return _ptr->_clone(); }
    operator const GradedAlgebraInterface<X>& () const { return *_ptr; }
    GradedAlgebra<X> create() const { return GradedAlgebra<X>(_ptr->_create()); }
    GradedAlgebra<X> clone() const { return GradedAlgebra<X>(_ptr->_clone()); }
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
    SymbolicAlgebra(SymbolicAlgebraInterface<X>* p) : _ptr(p) { }
    SymbolicAlgebra(const SymbolicAlgebraInterface<X>& a) : _ptr(a._clone()) { }
    SymbolicAlgebra(const SymbolicAlgebra<X>& a) : _ptr(a._ptr->_clone()) { }
    //! \brief Create the representation of the operator \a op applied to \a a.
    SymbolicAlgebra(OperatorCode op, const SymbolicAlgebra<X>& a) : _ptr(a._ptr->_apply(op)) { }
    SymbolicAlgebra<X>& operator=(Int c) { *this = this->create(); this->iadd(c); return *this; }
    SymbolicAlgebra<X>& operator=(const X& c) { *this = this->create(); this->iadd(c); return *this; }
    SymbolicAlgebra<X>& operator=(const Algebra<X>& a) { this->_ptr=std::shared_ptr< AlgebraInterface<X> >(a._ptr->_clone()); return *this; }
    operator Algebra<X> () const { return _ptr->_clone(); }
    operator const SymbolicAlgebraInterface<X>& () const { return *_ptr; }
    SymbolicAlgebra<X> create() const { return SymbolicAlgebra<X>(_ptr->_create()); }
    SymbolicAlgebra<X> clone() const { return SymbolicAlgebra<X>(_ptr->_clone()); }
    OutputStream& write(OutputStream& os) const { return _ptr->write(os); }
  public:
    Void iadd(const X& c) { _ptr->_iadd(c); }
    Void imul(const X& c) { _ptr->_imul(c); }
    Void isma(const X& c, const SymbolicAlgebra<X>& x) { _ptr->_isma(c,*x._ptr); }
    Void ifma(const SymbolicAlgebra<X>& x1, const SymbolicAlgebra<X>& x2) { _ptr->_ifma(*x1._ptr,*x2._ptr); }
};

template<class X> inline Algebra<X> AlgebraInterface<X>::create() const { return this->_create(); }
template<class X> inline Algebra<X> AlgebraInterface<X>::clone() const { return this->_clone(); }
template<class X> inline NormedAlgebra<X> NormedAlgebraInterface<X>::create() const { return this->_create(); }
template<class X> inline NormedAlgebra<X> NormedAlgebraInterface<X>::clone() const { return this->_clone(); }

template<class X> OutputStream& operator<<(OutputStream& os, const Algebra<X>& x) { return x.write(os); }
template<class X> OutputStream& operator<<(OutputStream& os, const NormedAlgebra<X>& x) { return x.write(os); }

template<class X> Algebra<X> create(const Algebra<X>& a) { return a.create(); }
template<class X> Algebra<X> copy(const Algebra<X>& a) { return a.clone(); }

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
    ARIADNE_FAIL_MSG("Cannot apply operator "<<op<<" to "<<a<<"\n");
}

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_H */
