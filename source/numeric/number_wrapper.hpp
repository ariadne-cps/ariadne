/***************************************************************************
 *            numeric/number_wrapper.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/number_wrapper.hpp
 *  \brief
 */



#ifndef ARIADNE_NUMBER_WRAPPER_HPP
#define ARIADNE_NUMBER_WRAPPER_HPP

#include "../utility/module.hpp"
#include "../numeric/paradigm.hpp"

#include "number_interface.hpp"

#include "number.hpp"

#include "logical.hpp"
#include "builtin.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"
#include "float_value.hpp"
#include "float_ball.hpp"
#include "float_bounds.hpp"
#include "float_upper_bound.hpp"
#include "float_lower_bound.hpp"
#include "float_approximation.hpp"
#include "float_error.hpp"

#include "../symbolic/templates.hpp"
#include "../numeric/operators.hpp"

namespace Ariadne {

/************ Number *********************************************************/

template<class... AWS> struct Aware;
template<class X> class NumberMixin;
template<class X> class NumberWrapper;

template<class X> inline X const* extract(NumberInterface const* y) {
     return dynamic_cast<NumberWrapper<X>const*>(y);
}

inline OutputStream& operator<<(OutputStream& os, NumberInterface const& y) { return y._write(os); }


// FIXME: Should test for other potential infinities
inline Comparison cmp(NumberInterface const& y1, NumberInterface const& y2) {
    Comparison res;
    FloatDPValue const* x1=extract<FloatDPValue>(&y1);
    FloatDPValue const* x2=extract<FloatDPValue>(&y2);
    if(x1) {
        if(x2) { res= cmp(ExactDouble(x1->raw().get_d()),ExactDouble(x2->raw().get_d())); }
        else { res= cmp(ExactDouble(x1->raw().get_d()),y2._get_q()); }
    } else {
        if(x2) { res= cmp(y1._get_q(),ExactDouble(x2->raw().get_d())); }
        else { res= cmp(y1._get_q(),y2._get_q()); }
    }
    return res;
}

template<class I, class OP, class Y> struct OperableInterface {
    virtual ~OperableInterface() = default;
    virtual I* _apply_left(OP op, Y const& y) const = 0;
    virtual I* _apply_right(OP op, Y const& y) const = 0;
};
template<class X, class I, class OP, class Y> struct OperableMixin : virtual OperableInterface<I,OP,Y> {
    static X const& _cast(OperableMixin<X,I,OP,Y> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    template<class R> static I* _make_wrapper(R&& r) { return new NumberWrapper<R>(r); }
    virtual I* _apply_left(OP op, Y const& other) const { return _make_wrapper(op(_cast(*this),other)); }
    virtual I* _apply_right(OP op, Y const& other) const { return _make_wrapper(op(other,_cast(*this))); }
};
template<class X, class I, class OP, class AW> struct Operable;
template<class X, class I, class OP, class Y, class... YS> struct Operable<X,I,OP,Aware<Y,YS...>>
    : OperableMixin<X,I,OP,Y>, Operable<X,I,OP,Aware<YS...>> { };
template<class X, class I, class OP> struct Operable<X,I,OP,Aware<>> { };


template<class OP> inline NumberInterface* make_symbolic(OP op, NumberInterface const* yp1, NumberInterface const* yp2) {
    Handle<NumberInterface> y1(const_cast<NumberInterface*>(yp1)->shared_from_this());
    Handle<NumberInterface> y2(const_cast<NumberInterface*>(yp2)->shared_from_this());
    String yc1=yp1->_class_name(); String yc2=yp2->_class_name();
    ARIADNE_THROW(DispatchException,op<<"(Number y1, Number y2) with y1="<<*yp1<<", y2="<<*yp2,"No dispatch for "<<op<<"("<<yc1<<", "<<yc2<<")");
}


template<class I, class X, class OP> inline I* _apply(X const& self, OP op, I const* self_ptr, I const* other_ptr) {
    auto aware_other_ptr=dynamic_cast<OperableInterface<I,OP,X>const*>(other_ptr);
    if(aware_other_ptr) { return aware_other_ptr->_apply_right(op,self); }
    else { return other_ptr->_rapply(op,self_ptr); }
}
template<class I, class X, class OP> inline I* _rapply(X const& self, OP op, I const* self_ptr, I const* other_ptr) {
    auto aware_other_ptr=dynamic_cast<OperableInterface<I,OP,X>const*>(other_ptr);
    if(aware_other_ptr) { return aware_other_ptr->_apply_left(op,self); }
    else { return make_symbolic(op,other_ptr,self_ptr); }
}



template<class X, class I, class OP> struct UnaryOperationMixin : public virtual I {
    static X const& _cast(UnaryOperationMixin<X,I,OP> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    template<class R> static I* _make_wrapper(R&& r) { return new NumberWrapper<R>(r); }
    virtual I* _apply(OP op) const final { return _make_wrapper(op(_cast(*this))); }
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


template<class X, class I, class AW> struct FieldLatticeAware
    : Operable<X,I,BinaryElementaryOperator,AW> {
};

template<class X, class I, class N> struct ElementaryGradedOperationsMixin : public virtual I {
    using I::_apply;
    static X const& _cast(ElementaryGradedOperationsMixin<X,I,N> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    template<class R> static I* _make_wrapper(R&& r) { return new NumberWrapper<R>(r); }
    virtual I* _apply(GradedElementaryOperator op, N n) const final { return _make_wrapper(op(std::forward<const X>(_cast(*this)),n)); }
};

template<class X, class I> struct ElementaryUnaryOperationsMixin : public virtual I {
    using I::_apply;
    static X const& _cast(ElementaryUnaryOperationsMixin<X,I> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    template<class R> static I* _make_wrapper(R&& r) { return new NumberWrapper<R>(r); }
    virtual I* _apply(UnaryElementaryOperator op) const final { return _make_wrapper(op(_cast(*this))); }
};


template<class X, class I, class J=I> struct ElementaryBinaryOperationsMixin : public virtual J {
    using J::_rapply; using J::_apply;
    static X const& _cast(ElementaryBinaryOperationsMixin<X,I,J> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    virtual I* _apply(BinaryElementaryOperator op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(BinaryElementaryOperator op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
};

/*
template<class X, class I, class J=I> struct AwareFieldMixin : public virtual J {
    using J::_rapply; using J::_apply;
    static X const& _cast(AwareFieldMixin<X,I,J> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    virtual I* _apply(BinaryElementaryOperator op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(BinaryElementaryOperator op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
};

template<class X, class I, class J=I> struct AwareLatticeMixin : public virtual J {
    using J::_rapply; using J::_apply;
    static X const& _cast(AwareLatticeMixin<X,I,J> const& self) { return static_cast<NumberMixin<X> const&>(self); }
    virtual I* _apply(Max op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _apply(Min op, I const* other) const final { return Ariadne::_apply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(Max op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
    virtual I* _rapply(Min op, I const* other) const final { return Ariadne::_rapply<I,X>(_cast(*this),op,this,other); }
};
*/

template<class X, class I, class W> struct SameArithmeticMixin : public virtual I {
    X const& _cast(X const& self) { return self; }
    X const& _cast(I const& other) { return dynamic_cast<Wrapper<X,I>const&>(other); }
    I* _heap_move(X&& x) { return new W(x); }
    virtual I* _apply(BinaryArithmeticOperator op, I const* other) const final { return _heap_move(op(_cast(*this),_cast(*other))); }
    virtual I* _rapply(BinaryArithmeticOperator op, I const* other) const final { return _heap_move(op(_cast(*other),_cast(*this))); }
};

template<class X> class NumberGetterMixin : public virtual NumberInterface {
  public:
  //  operator X const& () const { return static_cast<NumberMixin<X>const&>(*this); }
    static X const& _cast(NumberGetterMixin<X> const& self) { return static_cast<NumberWrapper<X>const&>(static_cast<NumberMixin<X>const&>(self)); }
    static X& _cast(NumberGetterMixin<X>& self) { return static_cast<NumberWrapper<X>&>(static_cast<NumberMixin<X>&>(self)); }

    typedef Paradigm<X> P;
    friend class Number<P>;

    virtual NumberInterface* _copy() const override { return new NumberWrapper<X>(_cast(*this)); }
    virtual NumberInterface* _move() override { return new NumberWrapper<X>(std::move(_cast(*this))); }

    // FIXME: Proper comparisons for ExactNumber.
    virtual LogicalValue _equals(NumberInterface const& y) const override {
        if (this->_paradigm() == ParadigmCode::EXACT && y._paradigm() == ParadigmCode::EXACT) {
            return LogicalValue( cmp(*this,y)==Comparison::EQUAL ? LogicalValue::TRUE : LogicalValue::FALSE ); }
        if (this->_paradigm() == ParadigmCode::VALIDATED && y._paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(BoundedTag(),dp) == y._get(BoundedTag(),dp)); }
        else {
            return LogicalValue(this->_get(ApproximateTag(),dp) == y._get(ApproximateTag(),dp)); }
    }
    virtual LogicalValue _less(NumberInterface const& y) const override {
        if (this->_paradigm() == ParadigmCode::EXACT && y._paradigm() == ParadigmCode::EXACT) {
            return LogicalValue( cmp(*this,y)==Comparison::LESS ? LogicalValue::TRUE : LogicalValue::FALSE ); }
        else if (this->_paradigm() == ParadigmCode::VALIDATED && y._paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(BoundedTag(),dp) < y._get(BoundedTag(),dp)); }
        else {
            return LogicalValue(this->_get(ApproximateTag(),dp) < y._get(ApproximateTag(),dp));
        }
    }

    virtual Rational _get_q() const override {
        return this->_get_as<Rational>(); }

    virtual FloatDPBall _get(MetricTag,DoublePrecision pr,DoublePrecision pre) const override {
        return this->_get_as<FloatDPBall>(pr); }
    virtual FloatDPBounds _get(BoundedTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPBounds>(pr); }
    virtual FloatDPUpperBound _get(UpperTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPUpperBound>(pr); }
    virtual FloatDPLowerBound _get(LowerTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPLowerBound>(pr); }
    virtual FloatDPApproximation _get(ApproximateTag,DoublePrecision pr) const override {
        return this->_get_as<FloatDPApproximation>(pr); }
    virtual FloatMPDPBall _get(MetricTag,MultiplePrecision pr, DoublePrecision pre) const override {
        return this->_get_as<FloatMPDPBall>(pr,pre); }
    virtual FloatMPBall _get(MetricTag, MultiplePrecision pr, MultiplePrecision pre) const override {
        return this->_get_as<FloatMPBall>(pr,pre); }
    virtual FloatMPBounds _get(BoundedTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPBounds>(pr); }
    virtual FloatMPUpperBound _get(UpperTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPUpperBound>(pr); }
    virtual FloatMPLowerBound _get(LowerTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPLowerBound>(pr); }
    virtual FloatMPApproximation _get(ApproximateTag, MultiplePrecision pr) const override {
        return this->_get_as<FloatMPApproximation>(pr); }

    virtual ParadigmCode _paradigm() const override { return P::code(); }
    virtual String _class_name() const override { return class_name<X>(); }
    virtual OutputStream& _write(OutputStream& os) const override { return os << _cast(*this); }

  private:
    template<class R, EnableIf<IsConstructible<R,X>> = dummy>
        inline R _get_as() const { return static_cast<R>(_cast(*this)); }
    template<class R, DisableIf<IsConstructible<R,X>> = dummy>
        inline R _get_as() const { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << "\n"; throw ParadigmError(); }
    template<class R, class PR, EnableIf<IsConstructible<R,X,PR>> = dummy>
        inline R _get_as(PR pr) const { return R(_cast(*this),pr); }
    template<class R, class PR, DisableIf<IsConstructible<R,X,PR>> = dummy>
        inline R _get_as(PR pr) const { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << " with precision " << pr << "\n"; throw ParadigmError(); }
    template<class R, class PR, class PRE, EnableIf<IsConstructible<R,X,PR,PRE>> = dummy>
        inline R _get_as(PR pr, PRE pre) const { return R(_cast(*this),pr,pre); }
    template<class R, class PR, class PRE, DisableIf<IsConstructible<R,X,PR,PRE>> = dummy>
        inline R _get_as(PR pr, PRE pre) const { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << " with precision " << pr << " and error precision " << pre << "\n"; throw ParadigmError(); }
};

template<class X> struct DispatchingTraits { typedef Aware<X> AwareOfTypes; };
template<class X> using Awares = typename DispatchingTraits<X>::AwareOfTypes;

template<class X> class NumberMixin
    : public ElementaryBinaryOperationsMixin<X,NumberInterface>
    , public ElementaryUnaryOperationsMixin<X,NumberInterface>
    , public ElementaryGradedOperationsMixin<X,NumberInterface,Int>
    , public FieldLatticeAware<X,NumberInterface,Awares<X>>
    , public NumberGetterMixin<X>
{
  public:
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

#endif /* ARIADNE_NUMBER_WRAPPER_HPP */
