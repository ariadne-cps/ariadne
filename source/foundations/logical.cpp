/***************************************************************************
 *            foundations/logical.cpp
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

/*! \file foundations/logical.cpp
 *  \brief
 */

#include "utility/stdlib.hpp"
#include "utility/string.hpp"
#include "utility/macros.hpp"
#include "logical.hpp"

namespace Ariadne {


namespace Detail {

class LogicalConstant : public LogicalInterface {
    LogicalValue _value;
  public:
    explicit LogicalConstant(LogicalValue value) : _value(value) { }
    operator LogicalValue() const { return _value; }
  private:
    LogicalInterface* _copy() const override { return new LogicalConstant(*this); }
    LogicalValue _check(Effort) const override { return _value; }
    OutputStream& _write(OutputStream& os) const override { return os << _value; }
};

enum class LogicalOperation { NOT, AND, OR, XOR, EQUAL };

inline char const* operation_name(LogicalOperation op) {
    switch(op) {
        case LogicalOperation::NOT: return "not";
        case LogicalOperation::AND: return "and";
        case LogicalOperation::OR: return "or";
        case LogicalOperation::XOR: return "xor";
        case LogicalOperation::EQUAL: return "equal";
    }
    return "logical";
}
inline LogicalValue apply(LogicalOperation op, LogicalValue v) {
    return op==LogicalOperation::NOT ? !v : LogicalValue::INDETERMINATE;
}
inline LogicalValue apply(LogicalOperation op, LogicalValue l, LogicalValue r) {
    switch(op) {
        case LogicalOperation::AND: return l&&r;
        case LogicalOperation::OR: return l||r;
        case LogicalOperation::XOR: return l^r;
        case LogicalOperation::EQUAL: return l==r;
        case LogicalOperation::NOT: break;
    }
    return LogicalValue::INDETERMINATE;
}
class UnaryLogicalExpression : public LogicalInterface {
    LogicalOperation _op; LogicalHandle _arg;
  public:
    UnaryLogicalExpression(LogicalOperation op, LogicalHandle arg):_op(op),_arg(arg){}
  private:
    LogicalInterface* _copy() const override { return new UnaryLogicalExpression(*this); }
    LogicalValue _check(Effort e) const override { return apply(_op,_arg.check(e)); }
    OutputStream& _write(OutputStream& os) const override { return os<<operation_name(_op)<<"("<<_arg<<")"; }
};
class BinaryLogicalExpression : public LogicalInterface {
    LogicalOperation _op; LogicalHandle _lhs; LogicalHandle _rhs;
  public:
    BinaryLogicalExpression(LogicalOperation op, LogicalHandle lhs, LogicalHandle rhs):_op(op),_lhs(lhs),_rhs(rhs){}
  private:
    LogicalInterface* _copy() const override { return new BinaryLogicalExpression(*this); }
    LogicalValue _check(Effort e) const override { return apply(_op,_lhs.check(e),_rhs.check(e)); }
    OutputStream& _write(OutputStream& os) const override { return os<<operation_name(_op)<<"("<<_lhs<<","<<_rhs<<")"; }
};


LogicalInterface* new_logical_pointer_from_value(LogicalValue v) {
    return new LogicalConstant(v);
}

LogicalValue logical_value_from_pointer(LogicalInterface* ptr) {
    auto vlptr=dynamic_cast<LogicalConstant*>(ptr);
    if (!vlptr) { throw std::runtime_error("logical_type_from_pointer: No conversion from abstract to concrete logical value"); }
    return *vlptr;
}

LogicalHandle LogicalHandle::constant(LogicalValue l) {
    return LogicalHandle(make_handle<const LogicalConstant>(l));
}

LogicalHandle operator&&(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(make_handle<const BinaryLogicalExpression>(LogicalOperation::AND,l1,l2));
}

LogicalHandle operator||(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(make_handle<const BinaryLogicalExpression>(LogicalOperation::OR,l1,l2));
}

LogicalHandle operator==(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(make_handle<const BinaryLogicalExpression>(LogicalOperation::EQUAL,l1,l2));
}

LogicalHandle operator^(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(make_handle<const BinaryLogicalExpression>(LogicalOperation::XOR,l1,l2));
}

LogicalHandle operator!(LogicalHandle l) {
    return LogicalHandle(make_handle<const UnaryLogicalExpression>(LogicalOperation::NOT,l));
}

LogicalHandle conjunction(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(make_handle<const BinaryLogicalExpression>(LogicalOperation::AND,l1,l2));
}

LogicalHandle disjunction(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(make_handle<const BinaryLogicalExpression>(LogicalOperation::OR,l1,l2));
}

LogicalHandle negation(LogicalHandle l) {
    return LogicalHandle(make_handle<const UnaryLogicalExpression>(LogicalOperation::NOT,l));
}

LogicalHandle equality(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(make_handle<const BinaryLogicalExpression>(LogicalOperation::EQUAL,l1,l2));
}


LogicalHandle exclusive(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(make_handle<const BinaryLogicalExpression>(LogicalOperation::XOR,l1,l2));
}


LogicalValue operator==(LogicalValue l1, LogicalValue l2) {
    switch (l1) {
        case LogicalValue::TRUE:
            return l2;
        case LogicalValue::LIKELY:
            switch (l2) { case LogicalValue::TRUE: return LogicalValue::LIKELY; case LogicalValue::FALSE: return LogicalValue::UNLIKELY; default: return l2; }
        case LogicalValue::INDETERMINATE:
            return LogicalValue::INDETERMINATE;
        case LogicalValue::UNLIKELY:
            switch (l2) { case LogicalValue::TRUE: return LogicalValue::UNLIKELY; case LogicalValue::FALSE: return LogicalValue::LIKELY; default: return not l2; }
        case LogicalValue::FALSE:
            return not l2;
        default:
            return LogicalValue::INDETERMINATE;
    }
}

OutputStream& operator<<(OutputStream& os, LogicalValue l) {
    switch(l) {
        case LogicalValue::TRUE: os << "true"; break;
        case LogicalValue::LIKELY: os << "likely";  break;
        case LogicalValue::INDETERMINATE: os << "indeterminate";  break;
        case LogicalValue::UNLIKELY: os << "unlikely"; break;
        case LogicalValue::FALSE: os << "false"; break;
        default: ARIADNE_FAIL_MSG("Unhandled LogicalValue for output streaming.");
    }
    return os;
}

} // namespace Detail

Nat Effort::_default = 0u;

const Indeterminate indeterminate = Indeterminate();

Bool NondeterministicBoolean::_choose(LowerKleenean p1, LowerKleenean p2) {
    Effort eff(0u);
    while(true) {
        if(definitely(p1.check(eff))) { return true; }
        if(definitely(p2.check(eff))) { return false; }
        ++eff;
    }
}


template<> String class_name<ExactTag>() { return "Exact"; }
template<> String class_name<EffectiveTag>() { return "Effective"; }
template<> String class_name<ValidatedTag>() { return "Validated"; }
//template<> String class_name<UpperTag>() { return "Upper"; }
//template<> String class_name<LowerTag>() { return "Lower"; }
template<> String class_name<ApproximateTag>() { return "Approximate"; }

template<> String class_name<Bool>() { return "Bool"; }
template<> String class_name<Boolean>() { return "Boolean"; }
template<> String class_name<Sierpinskian>() { return "Sierpinskian"; }
template<> String class_name<NegatedSierpinskian>() { return "NegatedSierpinskian"; }
template<> String class_name<Kleenean>() { return "Kleenean"; }
template<> String class_name<LowerKleenean>() { return "LowerKleenean"; }
template<> String class_name<UpperKleenean>() { return "UpperKleenean"; }
template<> String class_name<ValidatedSierpinskian>() { return "ValidatedSierpinskian"; }
template<> String class_name<ValidatedNegatedSierpinskian>() { return "ValidatedNegatedSierpinskian"; }
template<> String class_name<ValidatedKleenean>() { return "ValidatedKleenean"; }
template<> String class_name<ValidatedLowerKleenean>() { return "ValidatedLowerKleenean"; }
template<> String class_name<ValidatedUpperKleenean>() { return "ValidatedUpperKleenean"; }
template<> String class_name<ApproximateKleenean>() { return "ApproximateKleenean"; }

} // namespace Ariadne

#include "utility/array.hpp"

namespace Ariadne {

SizeType nondeterministic_choose_index(Array<LowerKleenean> const& p) {
    Effort eff(0u);
    while(true) {
        for (SizeType i=0; i!=p.size(); ++i) {
            if(definitely(p[i].check(eff))) { return i; }
        }
        ++eff;
    }
}

} // namespace Ariadne
