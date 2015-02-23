/***************************************************************************
 *            logical.cc
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

/*! \file logical.cc
 *  \brief
 */

#include "utility/stdlib.h"
#include "utility/string.h"
#include "expression/templates.h"

#include "logical.h"

namespace Ariadne {

namespace{
inline LogicalValue operator!(LogicalValue lv) { return negation(lv); }
inline LogicalValue operator&&(LogicalValue lv1, LogicalValue lv2) { return disjunction(lv1,lv2); }
inline LogicalValue operator||(LogicalValue lv1, LogicalValue lv2) { return conjunction(lv1,lv2); }
} // namespace

class LogicalConstant : public LogicalInterface {
    LogicalValue _v;
  public:
    LogicalConstant(LogicalValue v) : _v(v) { };
    virtual LogicalValue _check(Effort e) const { return _v; }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_v; }
};

template<class OP, class... ARGS> struct LogicalExpression;

template<class OP, class ARG> struct LogicalExpression<OP,ARG>
    : virtual LogicalInterface, ExpressionTemplate<OP,ARG>
{
    using ExpressionTemplate<OP,ARG>::ExpressionTemplate;
    virtual LogicalValue _check(Effort e) const { return this->_op(check(this->_arg,e)); }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << static_cast<ExpressionTemplate<OP,ARG>const&>(*this); }
};

template<class OP, class ARG1, class ARG2> struct LogicalExpression<OP,ARG1,ARG2>
    : virtual LogicalInterface, ExpressionTemplate<OP,ARG1,ARG2>
{
    using ExpressionTemplate<OP,ARG1,ARG2>::ExpressionTemplate;
    virtual LogicalValue _check(Effort e) const { return this->_op(check(this->_arg1,e),check(this->_arg2,e)); }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << static_cast<ExpressionTemplate<OP,ARG1,ARG2>const&>(*this); }
};

LogicalHandle::LogicalHandle(LogicalValue l)
    : _ptr(std::make_shared<LogicalConstant>(l)) {
}

LogicalHandle disjunction(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(std::make_shared<LogicalExpression<AndOp,LogicalHandle,LogicalHandle>>(AndOp(),l1,l2));
};

LogicalHandle conjunction(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(std::make_shared<LogicalExpression<OrOp,LogicalHandle,LogicalHandle>>(OrOp(),l1,l2));
};

LogicalHandle negation(LogicalHandle l) {
    return LogicalHandle(std::make_shared<LogicalExpression<NotOp,LogicalHandle>>(NotOp(),l));
};

LogicalValue equal(LogicalValue l1, LogicalValue l2) {
    switch (l1) {
        case LogicalValue::TRUE:
            return l2;
        case LogicalValue::LIKELY:
            switch (l2) { case LogicalValue::TRUE: return LogicalValue::LIKELY; case LogicalValue::FALSE: return LogicalValue::UNLIKELY; default: return l2; }
        case LogicalValue::INDETERMINATE:
            return LogicalValue::INDETERMINATE;
        case LogicalValue::UNLIKELY:
            switch (l2) { case LogicalValue::TRUE: return LogicalValue::UNLIKELY; case LogicalValue::FALSE: return LogicalValue::LIKELY; default: return negation(l2); }
        case LogicalValue::FALSE:
            return negation(l2);
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
    }
    return os;
}

template<> String class_name<Exact>() { return "Exact"; }
template<> String class_name<Effective>() { return "Effective"; }
template<> String class_name<Validated>() { return "Validated"; }
template<> String class_name<Bounded>() { return "Bounded"; }
template<> String class_name<Upper>() { return "Upper"; }
template<> String class_name<Lower>() { return "Lower"; }
template<> String class_name<Approximate>() { return "Approximate"; }

template<> String class_name<Bool>() { return "Bool"; }
template<> String class_name<Boolean>() { return "Boolean"; }
template<> String class_name<Tribool>() { return "Tribool"; }
template<> String class_name<Sierpinski>() { return "Sierpinski"; }
template<> String class_name<NegSierpinski>() { return "NegSierpinski"; }
template<> String class_name<Fuzzy>() { return "Fuzzy"; }

} // namespace Ariadne
