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

#include "numeric/operators.h"

namespace Ariadne {

/************ Number *********************************************************/


class UnaryOperatorInterface {
  public:
    virtual NumberInterface* _compute(Real) const = 0;
    virtual NumberInterface* _compute(Float64Value) const = 0;
    virtual NumberInterface* _compute(Float64Ball) const = 0;
    virtual NumberInterface* _compute(Float64Bounds) const = 0;
    virtual NumberInterface* _compute(Float64Approximation) const = 0;
    virtual NumberInterface* _compute(FloatMPValue) const = 0;
    virtual NumberInterface* _compute(FloatMPBall) const = 0;
    virtual NumberInterface* _compute(FloatMPBounds) const = 0;
    virtual NumberInterface* _compute(FloatMPApproximation) const = 0;
};

template<class X> class NumberWrapper;

template<class R> inline NumberInterface* _heap_move_number(R&& r) {
    return new NumberWrapper<R>(r); }

template<class R> inline NumberInterface* _heap_move_number(R const& r) {
    return new NumberWrapper<R>(r); }

template<class O, class A1> class UnaryOperator : public UnaryOperatorInterface {
    O _op; A1 _arg1;
  public:
    UnaryOperator(A1 a1) : _arg1(a1) { }
    UnaryOperator(O o, A1 a1) : _op(o), _arg1(a1) { }
    virtual NumberInterface* _compute(Real a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(Float64Value a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(Float64Ball a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(Float64Bounds a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(Float64Approximation a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(FloatMPValue a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(FloatMPBall a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(FloatMPBounds a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(FloatMPApproximation a2) const { return _heap_move_number(_op(_arg1,a2)); }
};

template<class T1, class T2> inline auto add(T1&& t1, T2&& t2) -> decltype(t1+t2) { return std::forward<T1>(t1) + std::forward<T2>(t2); }
template<class T1, class T2> inline auto sub(T1&& t1, T2&& t2) -> decltype(t1-t2) { return std::forward<T1>(t1) - std::forward<T2>(t2); }
template<class T1, class T2> inline auto mul(T1&& t1, T2&& t2) -> decltype(t1*t2) { return std::forward<T1>(t1) * std::forward<T2>(t2); }
template<class T1, class T2> inline auto div(T1&& t1, T2&& t2) -> decltype(t1/t2) { return std::forward<T1>(t1) / std::forward<T2>(t2); }

template<class O, class P> class UnaryOperator<O,Float<P,Precision64>> : public UnaryOperatorInterface {
    typedef Float<P,Precision64> A1;
    O _op; A1 _arg1;
  public:
    UnaryOperator(A1 a1) : _arg1(a1) { }
    UnaryOperator(O o, A1 a1) : _op(o), _arg1(a1) { }
    virtual NumberInterface* _compute(Real a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(Float64Value a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(Float64Ball a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(Float64Bounds a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(Float64Approximation a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(FloatMPValue a2) const { assert(false); }
    virtual NumberInterface* _compute(FloatMPBall a2) const { assert(false); }
    virtual NumberInterface* _compute(FloatMPBounds a2) const { assert(false); }
    virtual NumberInterface* _compute(FloatMPApproximation a2) const { assert(false); }
};

template<class O, class P> class UnaryOperator<O,Float<P,PrecisionMP>> : public UnaryOperatorInterface {
    typedef Float<P,PrecisionMP> A1;
    O _op; A1 _arg1;
  public:
    UnaryOperator(A1 a1) : _arg1(a1) { }
    UnaryOperator(O o, A1 a1) : _op(o), _arg1(a1) { }
    virtual NumberInterface* _compute(Real a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(Float64Value a2) const { assert(false); }
    virtual NumberInterface* _compute(Float64Ball a2) const { assert(false); }
    virtual NumberInterface* _compute(Float64Bounds a2) const { assert(false); }
    virtual NumberInterface* _compute(Float64Approximation a2) const { assert(false); }
    virtual NumberInterface* _compute(FloatMPValue a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(FloatMPBall a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(FloatMPBounds a2) const { return _heap_move_number(_op(_arg1,a2)); }
    virtual NumberInterface* _compute(FloatMPApproximation a2) const { return _heap_move_number(_op(_arg1,a2)); }
};

template<class X> class NumberWrapper
    : public virtual NumberInterface
    , public X
{
    static_assert(Not<IsSame<X,Handle<NumberInterface>>>::value,"X must be a concrete number, not a handle");
    typedef Paradigm<X> P;
    static_assert(Not<IsSame<X,Number<P>>>::value,"X must be a concrete number, not a generic number");
    friend class Number<P>;
  public:
    NumberWrapper(const X& a) : X(a) { }
    NumberWrapper(X&& a) : X(std::forward<X>(a)) { }
  private:
    inline static const X& _cast(const NumberWrapper<X>& x) {
        return static_cast<const X&>(x); }
    inline static const X& _cast(const NumberInterface& x) {
        const NumberWrapper<X>* p=dynamic_cast<const NumberWrapper<X>*>(&x);
        if(p==nullptr) { std::cerr << "Number "<<x<<" of type " << x._class_name() << " does not have the same type as "<<class_name<X>()<<"\n"; throw std::bad_cast(); }
        return static_cast<const X&>(*p); }
    virtual NumberInterface* _copy() const { return new NumberWrapper<X>(static_cast<const X&>(*this)); }
    virtual NumberInterface* _move() { return new NumberWrapper<X>(std::move(static_cast<X&>(*this))); }
    virtual NumberInterface* _nul() const {
        return _heap_move_number(0*_cast(*this)); }
    virtual NumberInterface* _pos() const {
        return _heap_move_number(pos(_cast(*this))); }
    virtual NumberInterface* _sqr() const {
        return _heap_move_number(sqr(_cast(*this))); }
    virtual NumberInterface* _neg() const {
        return _heap_move_number(neg(_cast(*this))); }
    virtual NumberInterface* _rec() const {
        return _heap_move_number(rec(_cast(*this))); }
    virtual NumberInterface* _pow(Int n) const {
        return _heap_move_number(pow(_cast(*this),n)); }
    virtual NumberInterface* _add(NumberInterface const& y) const {
        const NumberWrapper<X>* p=dynamic_cast<const NumberWrapper<X>*>(&y);
        if(p==nullptr) { return y._apply(UnaryOperator<Add,X>(_cast(*this))); }
        return _heap_move_number(add(_cast(*this),_cast(*p))); }
    virtual NumberInterface* _sub(NumberInterface const& y) const {
        const NumberWrapper<X>* p=dynamic_cast<const NumberWrapper<X>*>(&y);
        if(p==nullptr) { return y._apply(UnaryOperator<Sub,X>(_cast(*this))); }
        return _heap_move_number(sub(_cast(*this),_cast(y))); }
    virtual NumberInterface* _mul(NumberInterface const& y) const {
        const NumberWrapper<X>* p=dynamic_cast<const NumberWrapper<X>*>(&y);
        if(p==nullptr) { return y._apply(UnaryOperator<Mul,X>(_cast(*this))); }
        return _heap_move_number(mul(_cast(*this),_cast(y))); }
    virtual NumberInterface* _div(NumberInterface const& y) const {
        const NumberWrapper<X>* p=dynamic_cast<const NumberWrapper<X>*>(&y);
        if(p==nullptr) { return y._apply(UnaryOperator<Div,X>(_cast(*this))); }
        return _heap_move_number(div(_cast(*this),_cast(y))); }

    virtual NumberInterface* _sqrt() const {
        return _heap_move_number(sqrt(_cast(*this))); }
    virtual NumberInterface* _exp() const {
        return _heap_move_number(exp(_cast(*this))); }
    virtual NumberInterface* _log() const {
        return _heap_move_number(log(_cast(*this))); }
    virtual NumberInterface* _sin() const {
        return _heap_move_number(sin(_cast(*this))); }
    virtual NumberInterface* _cos() const {
        return _heap_move_number(cos(_cast(*this))); }
    virtual NumberInterface* _tan() const {
        return _heap_move_number(tan(_cast(*this))); }
    virtual NumberInterface* _atan() const {
        return _heap_move_number(atan(_cast(*this))); }
    virtual NumberInterface* _abs() const {
        return _heap_move_number(abs(_cast(*this))); }
    virtual NumberInterface* _max(NumberInterface const& y) const {
        const NumberWrapper<X>* p=&dynamic_cast<const NumberWrapper<X>&>(y);
        return _heap_move_number(max(_cast(*this),_cast(*p))); }
    virtual NumberInterface* _min(NumberInterface const& y) const {
        const NumberWrapper<X>* p=&dynamic_cast<const NumberWrapper<X>&>(y);
        return _heap_move_number(min(_cast(*this),_cast(*p))); }
    virtual LogicalValue _equals(NumberInterface const& y) const {
        if (this->_paradigm() == ParadigmCode::VALIDATED && y._paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(BoundedTag(),Precision64()) == y._get(BoundedTag(),Precision64())); }
        else { return LogicalValue(this->_get(ApproximateTag(),Precision64()) == y._get(ApproximateTag(),Precision64())); } }
    virtual LogicalValue _less(NumberInterface const& y) const {
        if (this->_paradigm() == ParadigmCode::VALIDATED && y._paradigm() == ParadigmCode::VALIDATED) {
            return LogicalValue(this->_get(BoundedTag(),Precision64()) < y._get(BoundedTag(),Precision64())); }
        else { return LogicalValue(this->_get(ApproximateTag(),Precision64()) < y._get(ApproximateTag(),Precision64())); } }
    virtual Float64Ball _get(MetricTag,Precision64) const {
        return this->_get_as<Float64Ball>(); }
    virtual Float64Bounds _get(BoundedTag,Precision64) const {
        return this->_get_as<Float64Bounds>(); }
    virtual Float64UpperBound _get(UpperTag,Precision64) const {
        return this->_get_as<Float64UpperBound>(); }
    virtual Float64LowerBound _get(LowerTag,Precision64) const {
        return this->_get_as<Float64LowerBound>(); }
    virtual Float64Approximation _get(ApproximateTag,Precision64) const {
        return this->_get_as<Float64Approximation>(); }
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
    virtual OutputStream& _write(OutputStream& os) const { return os << static_cast<const X&>(*this); }
    virtual OutputStream& write(OutputStream& os) const { return os << static_cast<const X&>(*this); }
    virtual NumberInterface* _apply(const UnaryOperatorInterface& o) const { return o._compute(static_cast<const X&>(*this)); }

  private:
    template<class R, typename std::enable_if<std::is_constructible<R,X>::value,Int>::type = 0>
        inline R _get_as() const { return static_cast<R>(static_cast<const X&>(*this)); }
    template<class R, typename std::enable_if<!std::is_constructible<R,X>::value,Int>::type = 0>
        inline R _get_as() const { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << "\n"; throw ParadigmError(); }
    template<class R, typename std::enable_if<std::is_constructible<R,X,PrecisionMP>::value,Int>::type = 0>
        inline R _get_as(PrecisionMP pr) const { return R(static_cast<const X&>(*this),pr); }
    template<class R, typename std::enable_if<!std::is_constructible<R,X,PrecisionMP>::value,Int>::type = 0>
        inline R _get_as(PrecisionMP) const { std::cerr<<"Warning: Cannot convert " << _cast(*this) << " of type " << this->_class_name() << " to " << class_name<R>() << " with given precision\n"; throw ParadigmError(); }
};



} // namespace Ariadne

#endif /* ARIADNE_NUMBER_WRAPPER_H */
