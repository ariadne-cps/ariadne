/***************************************************************************
 *            numeric/logical.cpp
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

/*! \file numeric/logical.cpp
 *  \brief
 */

#include "../utility/stdlib.hpp"
#include "../utility/string.hpp"
#include "../utility/macros.hpp"
#include "../numeric/sequence.hpp"
#include "../symbolic/templates.hpp"

#include "logical.hpp"
#include "integer.hpp"

namespace Ariadne {

namespace Detail {

class LogicalConstant : public LogicalInterface {
    LogicalValue _v;
  public:
    LogicalConstant(LogicalValue v) : _v(v) { };
    virtual LogicalValue _check(Effort e) const { return _v; }
    virtual OutputStream& _write(OutputStream& os) const { return os << this->_v; }
};

template<class OP, class... ARGS> struct LogicalExpression;

template<class OP, class ARG> struct LogicalExpression<OP,ARG>
    : virtual LogicalInterface, Symbolic<OP,ARG>
{
    using Symbolic<OP,ARG>::Symbolic;
    virtual LogicalValue _check(Effort e) const { return this->_op(check(this->_arg,e)); }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << static_cast<Symbolic<OP,ARG>const&>(*this); }
};

template<class OP, class ARG1, class ARG2> struct LogicalExpression<OP,ARG1,ARG2>
    : virtual LogicalInterface, Symbolic<OP,ARG1,ARG2>
{
    using Symbolic<OP,ARG1,ARG2>::Symbolic;
    virtual LogicalValue _check(Effort e) const { return this->_op(check(this->_arg1,e),check(this->_arg2,e)); }
    virtual OutputStream& _write(OutputStream& os) const {
        return os << static_cast<Symbolic<OP,ARG1,ARG2>const&>(*this); }
};

LogicalHandle::LogicalHandle(LogicalValue l)
    : _ptr(std::make_shared<LogicalConstant>(l)) {
}

LogicalHandle operator&&(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(std::make_shared<LogicalExpression<AndOp,LogicalHandle,LogicalHandle>>(AndOp(),l1,l2));
}

LogicalHandle operator||(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(std::make_shared<LogicalExpression<OrOp,LogicalHandle,LogicalHandle>>(OrOp(),l1,l2));
}

LogicalHandle operator==(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(std::make_shared<LogicalExpression<Equal,LogicalHandle,LogicalHandle>>(Equal(),l1,l2));
}

LogicalHandle operator^(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(std::make_shared<LogicalExpression<XOrOp,LogicalHandle,LogicalHandle>>(XOrOp(),l1,l2));
}

LogicalHandle operator!(LogicalHandle l) {
    return LogicalHandle(std::make_shared<LogicalExpression<NotOp,LogicalHandle>>(NotOp(),l));
}

LogicalHandle conjunction(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(std::make_shared<LogicalExpression<AndOp,LogicalHandle,LogicalHandle>>(AndOp(),l1,l2));
}

LogicalHandle disjunction(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(std::make_shared<LogicalExpression<OrOp,LogicalHandle,LogicalHandle>>(OrOp(),l1,l2));
}

LogicalHandle negation(LogicalHandle l) {
    return LogicalHandle(std::make_shared<LogicalExpression<NotOp,LogicalHandle>>(NotOp(),l));
}

LogicalHandle equality(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(std::make_shared<LogicalExpression<Equal,LogicalHandle,LogicalHandle>>(Equal(),l1,l2));
}


LogicalHandle exclusive(LogicalHandle l1, LogicalHandle l2) {
    return LogicalHandle(std::make_shared<LogicalExpression<XOrOp,LogicalHandle,LogicalHandle>>(XOrOp(),l1,l2));
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
        default: ARIADNE_FAIL_MSG("Unhandled LogicalValue for output streaming.\n");
    }
    return os;
}

template<> struct LogicalExpression<OrOp,Sequence<LowerKleenean>> : public LogicalInterface {
    Sequence<LowerKleenean> _seq;
  public:
    LogicalExpression(OrOp, Sequence<LowerKleenean> seq) : _seq(seq) { }
    LogicalValue _check(Effort eff) const {
        for(Natural k=0u; k!=eff.work(); ++k) {
            if ( definitely(_seq[k].check(eff)) ) { return LogicalValue::TRUE; }
        }
        return LogicalValue::INDETERMINATE;
    }
    OutputStream& _write(OutputStream& os) const {
        return os << "disjunction(" << _seq[0u] << "," << _seq[1u] << "," << _seq[2u] << ",...)";
    }
};


template<> struct LogicalExpression<AndOp,Sequence<UpperKleenean>> : public LogicalInterface {
    Sequence<UpperKleenean> _seq;
  public:
    LogicalExpression(AndOp, Sequence<UpperKleenean> seq) : _seq(seq) { }
    LogicalValue _check(Effort eff) const {
        for(Natural k=0u; k!=eff.work(); ++k) {
            if ( definitely(not _seq[k].check(eff)) ) { return LogicalValue::FALSE; }
        }
        return LogicalValue::INDETERMINATE;
    }
    OutputStream& _write(OutputStream& os) const {
        return os << "conjunction(" << _seq[0u] << "," << _seq[1u] << "," << _seq[2u] << ",...)";
    }
};

} // namespace Detail


LowerKleenean disjunction(Sequence<LowerKleenean> const& l) {
    return LowerKleenean(LogicalHandle(std::make_shared<Detail::LogicalExpression<OrOp,Sequence<LowerKleenean>>>(OrOp(),l))) ;
}
UpperKleenean conjunction(Sequence<UpperKleenean> const& l) {
    return UpperKleenean(LogicalHandle(std::make_shared<Detail::LogicalExpression<AndOp,Sequence<UpperKleenean>>>(AndOp(),l)));
}

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
template<> String class_name<BoundedTag>() { return "Bounded"; }
template<> String class_name<UpperTag>() { return "Upper"; }
template<> String class_name<LowerTag>() { return "Lower"; }
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
