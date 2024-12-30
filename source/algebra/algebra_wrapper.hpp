/***************************************************************************
 *            algebra/algebra_wrapper.hpp
 *
 *  Copyright  2010-20  Pieter Collins
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

/*! \file algebra/algebra_wrapper.hpp
 *  \brief Mixin class providing operations for (Banach) algebras.
 */

#ifndef ARIADNE_ALGEBRA_WRAPPER_HPP
#define ARIADNE_ALGEBRA_WRAPPER_HPP

#include "utility/macros.hpp"

#include "algebra/algebra_interface.hpp"
#include "algebra/algebra_concepts.hpp"
#include "numeric/operators.hpp"

namespace Ariadne {

template<class A, class X=NumericType<A>> class AlgebraWrapper;

template<class A, class X> requires AnAlgebraOver<A,X> class AlgebraWrapper<A,X>
    : public virtual AlgebraInterface<X>
    , public A
{
    static AlgebraInterface<X>* heap_move(A&& a) { return new AlgebraWrapper<A,X>(std::forward<A>(a)); }
    static A& static_algebra_cast(AlgebraWrapper<A,X>& a) { return static_cast<A&>(a); }
    static A const& static_algebra_cast(AlgebraWrapper<A,X> const& a) { return static_cast<A const&>(a); }
    static A const& dynamic_algebra_cast(AlgebraInterface<X> const& a) {
        A const* ap = dynamic_cast<AlgebraWrapper<A,X>const*>(&a);
        if (ap == nullptr) { ARIADNE_THROW(std::runtime_error,"bad_cast","Cannot cast AlgebraInterface "<<a<<" to "<<class_name<A>();); }
        return *ap; }
  public:
    typedef X NumericType;

    AlgebraWrapper(A&& a) : A(std::move(a)) { }
    AlgebraWrapper(A const& a) : A(a) { }
    virtual AlgebraInterface<X>* _copy() const override { return new AlgebraWrapper<A>(*this); }
    virtual AlgebraInterface<X>* _create_zero() const override { return heap_move(nul(*this)); }
    virtual AlgebraInterface<X>* _create_constant(X const& c) const override { A a=nul(*this); a=c; return heap_move(std::move(a)); }
    virtual AlgebraInterface<X>* _create_copy() const override { return new AlgebraWrapper<A>(*this); }
    virtual Void _iadd(const X& c) override{
        static_algebra_cast(*this) += c; }
    virtual Void _imul(const X& c) override{
        static_algebra_cast(*this) *= c; }
    virtual Void _isma(const X& c, const AlgebraInterface<X>& a) override {
        static_algebra_cast(*this) += c * dynamic_algebra_cast(a); }
    virtual Void _ifma(const AlgebraInterface<X>& a1, const AlgebraInterface<X>& a2) override {
        static_algebra_cast(*this) += dynamic_algebra_cast(a1) * dynamic_algebra_cast(a2); }
    virtual AlgebraInterface<X>* _apply(UnaryRingOperator op) const override {
        return heap_move(op(static_algebra_cast(*this))); }
    virtual AlgebraInterface<X>* _apply(BinaryRingOperator op, AlgebraInterface<X> const& a) const override {
        return heap_move(op(static_algebra_cast(*this),dynamic_algebra_cast(a))); }
    virtual AlgebraInterface<X>* _apply(BinaryFieldOperator op, X const& c) const override {
        return heap_move(op(static_algebra_cast(*this),c)); }
    virtual AlgebraInterface<X>* _rapply(BinaryRingOperator op, X const& c) const override {
        return heap_move(op(c,static_algebra_cast(*this))); };
    virtual AlgebraInterface<X>* _apply(GradedRingOperator op, Nat m) const override {
        return heap_move(op(static_algebra_cast(*this),m)); };
    virtual OutputStream& _write(OutputStream& os) const override { os << static_algebra_cast(*this); return os; }
};

template<class A, class X> requires AnInplaceAlgebraOver<A,X> class AlgebraWrapper<A,X>
    : public virtual AlgebraInterface<X>
    , public A
{
    static AlgebraInterface<X>* heap_move(A&& a) { return new AlgebraWrapper<A,X>(std::forward<A>(a)); }
  public:
    typedef X NumericType;

    AlgebraWrapper(A&& a) : A(std::move(a)) { }
    AlgebraWrapper(A const& a) : A(a) { }
    virtual AlgebraInterface<X>* _copy() const override { return new AlgebraWrapper<A>(*this); }
    virtual AlgebraInterface<X>* _create_zero() const override { return heap_move(this->A::create()); }
    virtual AlgebraInterface<X>* _create_constant(X const& c) const override { return heap_move(this->A::create_constant(c)); }
    virtual AlgebraInterface<X>* _create_copy() const override { return new AlgebraWrapper<A>(*this); }
    virtual Void _iadd(const X& c) override{ static_cast<A*>(this)->A::iadd(c); }
    virtual Void _imul(const X& c) override{ static_cast<A*>(this)->A::imul(c); }
    virtual Void _isma(const X& c, const AlgebraInterface<X>& x) override {
        static_cast<A*>(this)->A::isma(c,dynamic_cast<const A&>(x)); }
    virtual Void _ifma(const AlgebraInterface<X>& x1, const AlgebraInterface<X>& x2) override {
        static_cast<A*>(this)->A::ifma(dynamic_cast<const A&>(x1),dynamic_cast<const A&>(x2)); }
    virtual AlgebraInterface<X>* _apply(Neg op) const override {
        ARIADNE_NOT_IMPLEMENTED; }
    virtual AlgebraInterface<X>* _apply(BinaryRingOperator op, AlgebraInterface<X> const& c) const override {
        ARIADNE_NOT_IMPLEMENTED; }
    virtual AlgebraInterface<X>* _apply(BinaryFieldOperator op, X const& c) const override {
        ARIADNE_NOT_IMPLEMENTED; }
    virtual AlgebraInterface<X>* _rapply(BinaryRingOperator op, X const& c) const override {
        ARIADNE_NOT_IMPLEMENTED; }
    virtual AlgebraInterface<X>* _apply(GradedRingOperator op, Nat m) const override {
        ARIADNE_NOT_IMPLEMENTED; }
    virtual OutputStream& _write(OutputStream& os) const override { os << static_cast<const A&>(*this); return os; }
};


template<class A, class X=NumericType<A>> class NormedAlgebraWrapper
    : public virtual NormedAlgebraInterface<X>
    , public AlgebraWrapper<A,X>
{
    using typename AlgebraInterface<X>::ErrorType;
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


template<class A, class X> requires ATranscendentalAlgebraOver<A,X> class TranscendentalAlgebraWrapper
    : public virtual TranscendentalAlgebraInterface<X>
    , public A
{
    static TranscendentalAlgebraInterface<X>* heap_move(A&& a) { return new TranscendentalAlgebraWrapper<A,X>(std::forward<A>(a)); }
    static A const& static_algebra_cast(TranscendentalAlgebraWrapper<A,X> const& aw) { return static_cast<A const&>(aw); }
    static A const& dynamic_algebra_cast(TranscendentalAlgebraInterface<X> const& a) {
        A const* ap = dynamic_cast<TranscendentalAlgebraWrapper<A,X>const*>(&a);
        if (ap == nullptr) { ARIADNE_THROW(std::runtime_error,"bad_cast","Cannot cast TranscendentalAlgebraInterface "<<a<<" to "<<class_name<A>();); }
        return *ap; }
  public:
    typedef X NumericType;

    TranscendentalAlgebraWrapper(A&& a) : A(std::move(a)) { }
    TranscendentalAlgebraWrapper(A const& a) : A(a) { }
    virtual TranscendentalAlgebraInterface<X>* _copy() const override { return new TranscendentalAlgebraWrapper<A,X>(*this); }
    virtual TranscendentalAlgebraInterface<X>* _create_zero() const override { return heap_move(nul(*this)); }
    virtual TranscendentalAlgebraInterface<X>* _create_constant(X const& c) const override { A a=nul(*this); a=c; return heap_move(std::move(a)); }
    virtual TranscendentalAlgebraInterface<X>* _create_copy() const override { return new TranscendentalAlgebraWrapper<A,X>(*this); }
    virtual TranscendentalAlgebraInterface<X>* _apply(UnaryTranscendentalOperator op) const override {
        return heap_move(op(static_algebra_cast(*this))); }
    virtual TranscendentalAlgebraInterface<X>* _apply(BinaryFieldOperator op, TranscendentalAlgebraInterface<X> const& a) const override {
        return heap_move(op(static_algebra_cast(*this),dynamic_algebra_cast(a))); }
    virtual TranscendentalAlgebraInterface<X>* _apply(BinaryFieldOperator op, X const& c) const override {
        return heap_move(op(static_algebra_cast(*this),c)); }
    virtual TranscendentalAlgebraInterface<X>* _rapply(BinaryFieldOperator op, X const& c) const override {
        return heap_move(op(c,static_algebra_cast(*this))); };
    virtual TranscendentalAlgebraInterface<X>* _apply(GradedFieldOperator op, Int n) const override {
        return heap_move(op(static_algebra_cast(*this),n)); };
    virtual OutputStream& _write(OutputStream& os) const override { os << static_algebra_cast(*this); return os; }
};

template<class A, class X> requires AnElementaryAlgebraOver<A,X> class ElementaryAlgebraWrapper
    : public virtual ElementaryAlgebraInterface<X>
    , public A
{
    static ElementaryAlgebraInterface<X>* heap_move(A&& a) { return new ElementaryAlgebraWrapper<A,X>(std::forward<A>(a)); }
    static A const& static_algebra_cast(ElementaryAlgebraWrapper<A,X> const& aw) { return static_cast<A const&>(aw); }
    static A const& dynamic_algebra_cast(ElementaryAlgebraInterface<X> const& a) {
        A const* ap = dynamic_cast<ElementaryAlgebraWrapper<A,X>const*>(&a);
        if (ap == nullptr) { ARIADNE_THROW(std::runtime_error,"bad_cast","Cannot cast ElementaryAlgebraInterface "<<a<<" to "<<class_name<A>();); }
        return *ap; }
  public:
    typedef X NumericType;

    ElementaryAlgebraWrapper(A&& a) : A(std::move(a)) { }
    ElementaryAlgebraWrapper(A const& a) : A(a) { }
    virtual ElementaryAlgebraInterface<X>* _copy() const override { return new ElementaryAlgebraWrapper<A,X>(*this); }
    virtual ElementaryAlgebraInterface<X>* _create_zero() const override { return heap_move(nul(*this)); }
    virtual ElementaryAlgebraInterface<X>* _create_constant(X const& c) const override { A a=nul(*this); a=c; return heap_move(std::move(a)); }
    virtual ElementaryAlgebraInterface<X>* _create_copy() const override { return new ElementaryAlgebraWrapper<A,X>(*this); }
    virtual ElementaryAlgebraInterface<X>* _apply(UnaryElementaryOperator op) const override {
        return heap_move(op(static_algebra_cast(*this))); }
    virtual ElementaryAlgebraInterface<X>* _apply(BinaryElementaryOperator op, ElementaryAlgebraInterface<X> const& a) const override {
        return heap_move(op(static_algebra_cast(*this),dynamic_algebra_cast(a))); }
    virtual ElementaryAlgebraInterface<X>* _apply(BinaryElementaryOperator op, X const& c) const override {
        return heap_move(op(static_algebra_cast(*this),c)); }
    virtual ElementaryAlgebraInterface<X>* _rapply(BinaryElementaryOperator op, X const& c) const override {
        return heap_move(op(c,static_algebra_cast(*this))); };
    virtual ElementaryAlgebraInterface<X>* _apply(GradedElementaryOperator op, Int n) const override {
        return heap_move(op(static_algebra_cast(*this),n)); };
    virtual OutputStream& _write(OutputStream& os) const override { os << static_algebra_cast(*this); return os; }
};


} // namespace Ariadne

#endif /* ARIADNE_ALGEBRA_WRAPPER_HPP */
