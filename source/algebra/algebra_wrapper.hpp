/***************************************************************************
 *            algebra_mixin.hpp
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

/*! \file algebra_wrapper.hpp
 *  \brief Mixin class providing operations for (Banach) algebras.
 */

#ifndef ARIADNE_ALGEBRA_WRAPPER_HPP
#define ARIADNE_ALGEBRA_WRAPPER_HPP

#include "../algebra/algebra_interface.hpp"
#include "../algebra/algebra_mixin.hpp"
#include "../numeric/operators.hpp"

namespace Ariadne {

template<class A, class X=NumericType<A>> class AlgebraWrapper
    : public AlgebraMixin<AlgebraWrapper<A,X>,X>
    , public A
{
  public:
    using A::A;
    using AlgebraMixin<AlgebraWrapper<A,X>,X>::_apply;
    using AlgebraMixin<AlgebraWrapper<A,X>,X>::_rapply;
    AlgebraWrapper(A const& a) : A(a) { }
    friend OutputStream& operator<<(OutputStream& os, AlgebraWrapper<A,X> const& a) { return os << static_cast<A const&>(a); }
};

template<class X> template<class A> A ElementaryAlgebra<X>::extract() const {
    AlgebraInterface<X> const* p = this->_ptr.operator->();
    if constexpr(IsBaseOf<AlgebraInterface<X>,A>::value) {
        A const* ap=dynamic_cast<A const*>(p);
        if (!ap) { std::cerr << "*p=" << *p << "; " << typeid(p).name() << ", " << typeid(*p).name() << "\n"; }
        assert(ap); return *ap; }
    else {
        AlgebraWrapper<A,X> const* ap=dynamic_cast<AlgebraWrapper<A,X> const*>(this->_ptr.operator->());
        if (!ap) { std::cerr << "*p=" << *p << "; " << typeid(p).name() << ", " << typeid(*p).name() << "\n"; }
        assert(ap); return *ap; }
}

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
    using SymbolicAlgebraInterface<X>::_apply;
    using AlgebraWrapper<A,X>::_apply;
    virtual SymbolicAlgebraInterface<X>* _create_zero() const { return new A(static_cast<const A&>(*this).A::create()); }
    virtual SymbolicAlgebraInterface<X>* _create_copy() const { return new SymbolicAlgebraWrapper<A>(*this); }
    virtual SymbolicAlgebraInterface<X>* _apply(UnaryElementaryOperator op) { return new A(op,static_cast<const A&>(*this)); }
};


template<class A, class X> class ElementaryAlgebraWrapper
    : public virtual ElementaryAlgebraInterface<X>
    , public AlgebraWrapper<A,X>
{
  private:
    static A const& _cast(AlgebraInterface<X> const& a) { return static_cast<A const&>(dynamic_cast<ElementaryAlgebraWrapper<A,X>const&>(a)); }
    static A const& _cast(ElementaryAlgebraWrapper<A,X> const& a) { return static_cast<A const&>(a); }
    static ElementaryAlgebraInterface<X>* _make(A&& a) { return new ElementaryAlgebraWrapper<A,X>(std::move(a)); }
    template<class OP> static ElementaryAlgebraInterface<X>* _eval(OP op, ElementaryAlgebraWrapper<A,X> const& aw1, AlgebraInterface<X> const& ai2) {
        return _make(op(_cast(aw1),_cast(ai2))); }
    template<class OP> static ElementaryAlgebraInterface<X>* _eval(OP op, ElementaryAlgebraWrapper<A,X> const& aw1, X const& c2) {
        return _make(op(_cast(aw1),c2)); }
    template<class OP> static ElementaryAlgebraInterface<X>* _eval(OP op, X const& c1, ElementaryAlgebraWrapper<A,X> const& aw2) {
        return _make(op(c1,_cast(aw2))); }
    template<class OP> static ElementaryAlgebraInterface<X>* _eval(OP op, ElementaryAlgebraWrapper<A,X> const& aw) {
        return _make(op(_cast(aw))); }
  public:
    ElementaryAlgebraWrapper(A const& a) : AlgebraWrapper<A,X>(a) { }
    virtual ElementaryAlgebraInterface<X>* _create_zero() const { return new ElementaryAlgebraWrapper<A,X>(A()); }
    virtual ElementaryAlgebraInterface<X>* _create_constant(X const& c) const { return new ElementaryAlgebraWrapper<A,X>(A::constant(c)); }
    virtual ElementaryAlgebraInterface<X>* _create_copy() const { return new ElementaryAlgebraWrapper<A,X>(*this); }

    using AlgebraWrapper<A,X>::_apply;
    virtual ElementaryAlgebraInterface<X>* _apply(BinaryElementaryOperator op_, AlgebraInterface<X>const& other) const final {
        return _eval(op_,*this,other); }
    virtual ElementaryAlgebraInterface<X>* _apply(BinaryElementaryOperator op_, X const& cnst) const { return _eval(op_,*this,cnst); }
    virtual ElementaryAlgebraInterface<X>* _apply(BinaryFieldOperator op_, X const& cnst) const { return _eval(op_,*this,cnst); }
    virtual ElementaryAlgebraInterface<X>* _rapply(BinaryRingOperator op_, X const& cnst) const { return _eval(op_,cnst,*this); }
    virtual OutputStream& _write(OutputStream& os) const { return os << _cast(*this); }
    virtual ElementaryAlgebraInterface<X>* _apply(UnaryElementaryOperator op_) const {
        return new ElementaryAlgebraWrapper<A,X>(op_(static_cast<A const&>(*this))); }
    virtual ElementaryAlgebraInterface<X>* _apply(GradedElementaryOperator op_, Int n) const {
        return new ElementaryAlgebraWrapper<A,X>(op_(static_cast<A const&>(*this),n)); }
};



} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_WRAPPER_HPP */
