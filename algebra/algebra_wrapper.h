/***************************************************************************
 *            algebra_mixin.h
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

/*! \file algebra_wrapper.h
 *  \brief Mixin class providing operations for (Banach) algebras.
 */

#ifndef ARIADNE_ALGEBRA_WRAPPER_H
#define ARIADNE_ALGEBRA_WRAPPER_H

#include "algebra/algebra_interface.h"
#include "algebra/algebra_mixin.h"
#include "numeric/operators.h"

namespace Ariadne {

template<class A, class X=NumericType<A>> class AlgebraWrapper
    : public AlgebraMixin<AlgebraWrapper<A,X>,X>
    , public A
{
  public:
    using A::A;
    AlgebraWrapper(A const& a) : A(a) { }
    friend OutputStream& operator<<(OutputStream& os, AlgebraWrapper<A,X> const& a) { return os << static_cast<A const&>(a); }
};
/*
template<class A, class X=NumericType<A>> class AlgebraWrapper
    : public virtual AlgebraInterface<X>
    , public A
{
    static AlgebraInterface<X>* move_heap(A&& a) { return new AlgebraWrapper<A,X>(std::forward<A>(a)); }
  public:
    typedef X NumericType;

    AlgebraWrapper(A&& a) : A(std::move(a)) { }
    AlgebraWrapper(A const& a) : A(a) { }
    virtual AlgebraInterface<X>* _create_zero() const override { return move_heap(this->A::create()); }
    virtual AlgebraInterface<X>* _create_constant(X const& c) const override { return move_heap(this->A::create_constant(c)); }
    virtual AlgebraInterface<X>* _create_copy() const { return new A(*this); }
    virtual Void _iadd(const X& c) { static_cast<A*>(this)->A::iadd(c); }
    virtual Void _imul(const X& c) { static_cast<A*>(this)->A::imul(c); }
    virtual Void _isma(const X& c, const AlgebraInterface<X>& x) {
        static_cast<A*>(this)->A::isma(c,dynamic_cast<const A&>(x)); }
    virtual Void _ifma(const AlgebraInterface<X>& x1, const AlgebraInterface<X>& x2)  {
        static_cast<A*>(this)->A::ifma(dynamic_cast<const A&>(x1),dynamic_cast<const A&>(x2)); }
    virtual OutputStream& write(OutputStream& os) const { os << static_cast<const A&>(*this); return os; }
};
*/

template<class A, class X=NumericType<A>> class NormedAlgebraWrapper
    : public virtual NormedAlgebraInterface<X>
    , public AlgebraWrapper<A,X>
{
    using AlgebraWrapper<A,X>::AlgebraWrapper;
    virtual NormedAlgebraInterface<X>* _create_ball(ErrorType r) const { return new A(static_cast<const A&>(*this).A::create_ball(r)); }
    virtual NormedAlgebraInterface<X>* _create_constant(X c) const { return new A(static_cast<const A&>(*this).A::create_constant(c)); }
    virtual NormedAlgebraInterface<X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual NormedAlgebraInterface<X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
};

template<class A, class X=NumericType<A>> class GradedAlgebraWrapper
    : public virtual GradedAlgebraInterface<X>
    , public AlgebraWrapper<A,X>
{
    using AlgebraWrapper<A,X>::AlgebraWrapper;
    virtual GradedAlgebraInterface<X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual GradedAlgebraInterface<X>* _create_copy() const { return new A(static_cast<const A&>(*this)); }
    virtual GradedAlgebraInterface<X>* _apply(const Series<X>& f) const { return new A(compose(f,static_cast<const A&>(*this))); }
};

template<class A, class X=NumericType<A>> class SymbolicAlgebraWrapper
    : public virtual SymbolicAlgebraInterface<X>
    , public AlgebraWrapper<A,X>
{
    using AlgebraWrapper<A,X>::AlgebraWrapper;
    virtual SymbolicAlgebraInterface<X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual SymbolicAlgebraInterface<X>* _create_copy() const { return new SymbolicAlgebraWrapper<A>(*this); }
    virtual SymbolicAlgebraInterface<X>* _apply(OperatorCode op) { return new A(op,static_cast<const A&>(*this)); }
};

} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_WRAPPER_H */
