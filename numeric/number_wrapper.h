/***************************************************************************
 *            numeric/number_wrapper.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file numeric/number_wrapper.h
 *  \brief
 */



#ifndef ARIADNE_NUMBER_WRAPPER_H
#define ARIADNE_NUMBER_WRAPPER_H

#include "utility/module.h"
#include "numeric/paradigm.h"

#include "number_interface.h"

#include "number.h"
#include "logical.h"
#include "float64.h"
#include "floatmp.h"
#include "float-user.h"

#include "expression/templates.h"
#include "numeric/operators.h"

namespace Ariadne {

/************ Number *********************************************************/

template<class... AWS> class Aware;
template<class X> class NumberMixin;
template<class X> class NumberWrapper;


template<class I, class OP, class Y> struct OperableInterface {
    virtual ~OperableInterface() = default;
    virtual I* _apply_left(OP op, Y const& y) const = 0;
    virtual I* _apply_right(OP op, Y const& y) const = 0;
};
template<class X, class I, class OP, class Y> struct OperableMixin : virtual OperableInterface<I,OP,Y> {
    static X const& _cast(OperableMixin<X,I,OP,Y> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    template<class R> static I* _make_wrapper(R&& r) { return new NumberWrapper<R>(r); }
    virtual I* _apply_left(OP op, Y const& other) const { return _make_wrapper(op(other,_cast(*this))); }
    virtual I* _apply_right(OP op, Y const& other) const { return _make_wrapper(op(_cast(*this),other)); }
};
template<class X, class I, class OP, class AW> struct Operable;
template<class X, class I, class OP, class Y, class... YS> struct Operable<X,I,OP,Aware<Y,YS...>>
    : OperableMixin<X,I,OP,Y>, Operable<X,I,OP,Aware<YS...>> { };
template<class X, class I, class OP> struct Operable<X,I,OP,Aware<>> { };


template<class OP> inline NumberInterface* make_symbolic(OP op, NumberInterface const* yp1, NumberInterface const* yp2) {
    Handle<NumberInterface> y1(const_cast<NumberInterface*>(yp1)->shared_from_this());
    Handle<NumberInterface> y2(const_cast<NumberInterface*>(yp2)->shared_from_this());
    return nullptr;
};

template<class I, class X, class OP> inline I* _apply(X const& self, OP op, I const* self_ptr, I const* other_ptr) {
    auto aware_other_ptr=dynamic_cast<OperableInterface<I,OP,X>const*>(other_ptr);
    if(aware_other_ptr) { return aware_other_ptr->_apply_right(op,self); }
    else { other_ptr->_rapply(op,self_ptr); }
}
template<class I, class X, class OP> inline I* _rapply(X const& self, OP op, I const* self_ptr, I const* other_ptr) {
    auto aware_other_ptr=dynamic_cast<OperableInterface<I,OP,X>const*>(other_ptr);
    if(aware_other_ptr) { return aware_other_ptr->_apply_left(op,self); }
    else { return make_symbolic(op,self_ptr,other_ptr); }
}



template<class X, class I, class OP> struct UnaryOperationMixin : public virtual I {
    static X const& _cast(UnaryOperationMixin<X,I,OP> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    template<class R> static I* _make_wrapper(R&& r) { return new NumberWrapper<R>(r); }
    virtual I* _apply(OP op) const final { return _make_wrapper(pos(_cast(*this))); }
};

template<class I, class OP> struct BinaryOperationInterface {
    virtual ~BinaryOperationInterface() = default;
    virtual I* _apply(OP op, I const* other) const = 0;
    virtual I* _rapply(OP op, I const* other) const = 0;
};

template<class X, class I, class OP, class J=I> struct BinaryOperationMixin : public virtual J {
    static X const& _cast(BinaryOperationMixin<X,I,OP,J> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    virtual I* _apply(OP op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(OP op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
};


template<class X, class I, class AW> struct FieldAware
    : Operable<X,I,Add,AW>, Operable<X,I,Sub,AW>, Operable<X,I,Mul,AW>, Operable<X,I,Div,AW> {
};
template<class X, class I, class AW> struct LatticeAware
    : Operable<X,I,Max,AW>, Operable<X,I,Min,AW> {
};
template<class X, class I, class AW> struct LatticeFieldAware
    : FieldAware<X,I,AW>, LatticeAware<X,I,AW> {
};

template<class X, class I> struct UnaryOperationsMixin : public virtual I {
    static X const& _cast(UnaryOperationsMixin<X,I> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    template<class R> static I* _make_wrapper(R&& r) { return new NumberWrapper<R>(r); }
    virtual I* _apply(Pos op) const final { return _make_wrapper(pos(_cast(*this))); }
    virtual I* _apply(Neg op) const final { return _make_wrapper(neg(_cast(*this))); }
    virtual I* _apply(Sqr op) const final { return _make_wrapper(sqr(_cast(*this))); }
    virtual I* _apply(Rec op) const final { return _make_wrapper(rec(_cast(*this))); }
    virtual I* _apply(Pow op, Int n) const final { return _make_wrapper(pow(_cast(*this),n)); }
    virtual I* _apply(Sqrt op) const final { return _make_wrapper(sqrt(_cast(*this))); }
    virtual I* _apply(Exp op) const final { return _make_wrapper(exp(_cast(*this))); }
    virtual I* _apply(Log op) const final { return _make_wrapper(log(_cast(*this))); }
    virtual I* _apply(Sin op) const final { return _make_wrapper(sin(_cast(*this))); }
    virtual I* _apply(Cos op) const final { return _make_wrapper(cos(_cast(*this))); }
    virtual I* _apply(Tan op) const final { return _make_wrapper(tan(_cast(*this))); }
    virtual I* _apply(Atan op) const final { return _make_wrapper(atan(_cast(*this))); }
    virtual I* _apply(Abs op) const final { return _make_wrapper(abs(_cast(*this))); }
};


template<class X, class I, class J=I> struct AwareFieldMixin : public virtual J {
    static X const& _cast(AwareFieldMixin<X,I,J> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    virtual I* _apply(Add op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _apply(Sub op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _apply(Mul op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _apply(Div op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(Add op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(Sub op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(Mul op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(Div op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
};

template<class X, class I, class J=I> struct AwareLatticeMixin : public virtual J {
    static X const& _cast(AwareLatticeMixin<X,I,J> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    virtual I* _apply(Max op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _apply(Min op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(Max op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(Min op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
};


template<class X, class I, class W> struct SameArithmeticMixin : public virtual I {
    X const& _cast(X const& self) { return self; }
    X const& _cast(I const& other) { return dynamic_cast<Wrapper<X,I>const&>(other); }
    I* _heap_move(X&& x) { return new W(x); }
    virtual I* _add(I const* other) const final { return _heap_move(add(_cast(*this),_cast(*other))); }
    virtual I* _sub(I const* other) const final { return _heap_move(add(_cast(*this),_cast(*other))); }
    virtual I* _mul(I const* other) const final { return _heap_move(add(_cast(*this),_cast(*other))); }
    virtual I* _div(I const* other) const final { return _heap_move(add(_cast(*this),_cast(*other))); }
    virtual I* _radd(I const* other) const final { return _heap_move(add(_cast(*other),_cast(*this))); }
    virtual I* _rsub(I const* other) const final { return _heap_move(add(_cast(*other),_cast(*this))); }
    virtual I* _rmul(I const* other) const final { return _heap_move(add(_cast(*other),_cast(*this))); }
    virtual I* _rdiv(I const* other) const final { return _heap_move(add(_cast(*other),_cast(*this))); }
};


template<class X> struct NumberGetterMixin : public virtual NumberInterface {
  //  operator X const& () const { return static_cast<NumberMixin<X>const&>(*this); }
    static X const& _cast(NumberGetterMixin<X> const& self) { return static_cast<NumberMixin<X>const&>(self); }
    static X& _cast(NumberGetterMixin<X>& self) { return static_cast<NumberMixin<X>&>(self); }

    typedef Paradigm<X> P;
    friend class Number<P>;

    virtual NumberInterface* _copy() const { return new NumberWrapper<X>(_cast(*this)); }
    virtual NumberInterface* _move() { return new NumberWrapper<X>(std::move(_cast(*this))); }

    virtual LogicalValue _equals(NumberInterface const& y) const {
        if (this->_paradigm() == ParadigmCode::VALIDATED && y._paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(BoundedTag(),Precision64()) == y._get(BoundedTag(),Precision64())); }
        else { return LogicalValue(this->_get(ApproximateTag(),Precision64()) == y._get(ApproximateTag(),Precision64())); } }
    virtual LogicalValue _less(NumberInterface const& y) const {
        if (this->_paradigm() == ParadigmCode::VALIDATED && y._paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(BoundedTag(),Precision64()) < y._get(BoundedTag(),Precision64())); }
        else { return LogicalValue(this->_get(ApproximateTag(),Precision64()) < y._get(ApproximateTag(),Precision64())); } }

    virtual Float64Ball _get(MetricTag,Precision64 pr) const {
        return this->_get_as<Float64Ball>(pr); }
    virtual Float64Bounds _get(BoundedTag,Precision64 pr) const {
        return this->_get_as<Float64Bounds>(pr); }
    virtual Float64UpperBound _get(UpperTag,Precision64 pr) const {
        return this->_get_as<Float64UpperBound>(pr); }
    virtual Float64LowerBound _get(LowerTag,Precision64 pr) const {
        return this->_get_as<Float64LowerBound>(pr); }
    virtual Float64Approximation _get(ApproximateTag,Precision64 pr) const {
        return this->_get_as<Float64Approximation>(pr); }
    virtual FloatMPBall _get(MetricTag, PrecisionMP pr) const {
        return this->_get_as<FloatMPBall>(pr); }
    virtual FloatMPBounds _get(BoundedTag, PrecisionMP pr) const {
        return this->_get_as<FloatMPBounds>(pr); }
    virtual FloatMPUpperBound _get(UpperTag, PrecisionMP pr) const {
        return this->_get_as<FloatMPUpperBound>(pr); }
    virtual FloatMPLowerBound _get(LowerTag, PrecisionMP pr) const {
        return this->_get_as<FloatMPLowerBound>(pr); }
    virtual FloatMPApproximation _get(ApproximateTag, PrecisionMP pr) const {
        return this->_get_as<FloatMPApproximation>(pr); }

    virtual ParadigmCode _paradigm() const { return P::code(); }
    virtual String _class_name() const { return class_name<X>(); }
    virtual OutputStream& _write(OutputStream& os) const { return os << _cast(*this); }
    virtual OutputStream& write(OutputStream& os) const { return os << _cast(*this); }

  private:
    template<class R, EnableIf<IsConvertible<X,R>> = dummy>
        inline R _get_as() const { return static_cast<R>(_cast(*this)); }
    template<class R, DisableIf<IsConvertible<X,R>> = dummy>
        inline R _get_as() const { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << "\n"; throw ParadigmError(); }
    template<class R, class PR, EnableIf<IsConstructible<R,X,PR>> = dummy>
        inline R _get_as(PR pr) const { return R(_cast(*this),pr); }
    template<class R, class PR, DisableIf<IsConstructible<R,X,PR>> = dummy>
        inline R _get_as(PR pr) const { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << " with precision" << pr << "\n"; throw ParadigmError(); }
};

template<class X> struct DispatchingTraits { typedef Aware<X> AwareOfTypes; };
template<class X> using Awares = typename DispatchingTraits<X>::AwareOfTypes;

template<class X> struct NumberMixin
    : public AwareFieldMixin<X,NumberInterface>
    , public AwareLatticeMixin<X,NumberInterface>
    , public UnaryOperationsMixin<X,NumberInterface>
    , public FieldAware<X,NumberInterface,Awares<X>>
    , public LatticeAware<X,NumberInterface,Awares<X>>
    , public NumberGetterMixin<X>
{
    operator X const& () const { return static_cast<NumberWrapper<X>const&>(*this); }
    operator X& () { return static_cast<NumberWrapper<X>&>(*this); }
};


template<class X> class NumberWrapper
    : public X, public NumberMixin<X>
{
    inline static const X& _cast(const NumberWrapper<X>& x) { return static_cast<const X&>(x); }
    static_assert(Not<IsSame<X,Handle<NumberInterface>>>::value,"X must be a concrete number, not a handle");
    static_assert(Not<IsSame<X,Number<Paradigm<X>>>>::value,"X must be a concrete number, not a generic number");
  public:
    NumberWrapper(const X& a) : X(a) { }
    NumberWrapper(X&& a) : X(std::forward<X>(a)) { }
};


} // namespace Ariadne

#endif /* ARIADNE_NUMBER_WRAPPER_H */
