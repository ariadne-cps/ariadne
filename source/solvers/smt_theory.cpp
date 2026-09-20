/***************************************************************************
 *            solvers/smt_theory.cpp
 *
 *  Copyright  2026  Ariadne developers
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 */

#include "solvers/smt_theory.hpp"

#include "utility/exceptions.hpp"

namespace Ariadne {

SmtTheoryLiteral SmtTheoryLiteral::negated() const
{
    switch(_relation) {
        case SmtTheoryRelation::EQ: return SmtTheoryLiteral(_lhs,SmtTheoryRelation::NEQ,_rhs);
        case SmtTheoryRelation::NEQ: return SmtTheoryLiteral(_lhs,SmtTheoryRelation::EQ,_rhs);
        case SmtTheoryRelation::LEQ: return SmtTheoryLiteral(_lhs,SmtTheoryRelation::GT,_rhs);
        case SmtTheoryRelation::GEQ: return SmtTheoryLiteral(_lhs,SmtTheoryRelation::LT,_rhs);
        case SmtTheoryRelation::LT: return SmtTheoryLiteral(_lhs,SmtTheoryRelation::GEQ,_rhs);
        case SmtTheoryRelation::GT: return SmtTheoryLiteral(_lhs,SmtTheoryRelation::LEQ,_rhs);
        default: ARIADNE_FAIL_MSG("Unknown SMT theory relation");
    }
}

SmtTheoryLiteral make_smt_theory_literal(ContinuousPredicate const& predicate)
{
    SmtTheoryRelation relation;
    switch(predicate.code()) {
        case OperatorCode::EQ: relation=SmtTheoryRelation::EQ; break;
        case OperatorCode::NEQ: relation=SmtTheoryRelation::NEQ; break;
        case OperatorCode::LEQ: relation=SmtTheoryRelation::LEQ; break;
        case OperatorCode::GEQ: relation=SmtTheoryRelation::GEQ; break;
        case OperatorCode::LT: relation=SmtTheoryRelation::LT; break;
        case OperatorCode::GT: relation=SmtTheoryRelation::GT; break;
        default:
            ARIADNE_FAIL_MSG("Expected real comparison atom, got operator "<<predicate.code());
    }

    return SmtTheoryLiteral(predicate.cmp1<Real>(),relation,predicate.cmp2<Real>());
}

SmtTheoryAlternatives normalize_smt_theory_literal(SmtTheoryLiteral const& literal)
{
    RealExpression difference=literal.lhs()-literal.rhs();
    RealExpression opposite=literal.rhs()-literal.lhs();

    switch(literal.relation()) {
        case SmtTheoryRelation::EQ:
            return {{{difference,SmtTheoryPrimitiveRelation::EQ_ZERO}}};
        case SmtTheoryRelation::GEQ:
            return {{{difference,SmtTheoryPrimitiveRelation::GEQ_ZERO}}};
        case SmtTheoryRelation::GT:
            return {{{difference,SmtTheoryPrimitiveRelation::GT_ZERO}}};
        case SmtTheoryRelation::LEQ:
            return {{{opposite,SmtTheoryPrimitiveRelation::GEQ_ZERO}}};
        case SmtTheoryRelation::LT:
            return {{{opposite,SmtTheoryPrimitiveRelation::GT_ZERO}}};
        case SmtTheoryRelation::NEQ:
            return {
                {{difference,SmtTheoryPrimitiveRelation::GT_ZERO}},
                {{opposite,SmtTheoryPrimitiveRelation::GT_ZERO}}
            };
        default:
            ARIADNE_FAIL_MSG("Unknown SMT theory relation");
    }
}

SmtTheoryWeakLiteral weaken_smt_theory_literal(SmtTheoryPrimitiveLiteral const& literal,
                                               ExactDouble epsilon)
{
    ARIADNE_PRECONDITION(epsilon>0.0_x);

    switch(literal.relation()) {
        case SmtTheoryPrimitiveRelation::EQ_ZERO:
            return SmtTheoryWeakLiteral(
                literal.expression(),SmtTheoryWeakRelation::ABS_LEQ_EPSILON,epsilon);
        case SmtTheoryPrimitiveRelation::GEQ_ZERO:
            return SmtTheoryWeakLiteral(
                literal.expression(),SmtTheoryWeakRelation::GEQ_MINUS_EPSILON,epsilon);
        case SmtTheoryPrimitiveRelation::GT_ZERO:
            return SmtTheoryWeakLiteral(
                literal.expression(),SmtTheoryWeakRelation::GT_MINUS_EPSILON,epsilon);
        default:
            ARIADNE_FAIL_MSG("Unknown SMT primitive theory relation");
    }
}

OutputStream& operator<<(OutputStream& os, SmtTheoryRelation relation)
{
    switch(relation) {
        case SmtTheoryRelation::EQ: return os << "=";
        case SmtTheoryRelation::NEQ: return os << "!=";
        case SmtTheoryRelation::LEQ: return os << "<=";
        case SmtTheoryRelation::GEQ: return os << ">=";
        case SmtTheoryRelation::LT: return os << "<";
        case SmtTheoryRelation::GT: return os << ">";
        default: ARIADNE_FAIL_MSG("Unknown SMT theory relation");
    }
}

OutputStream& operator<<(OutputStream& os, SmtTheoryPrimitiveRelation relation)
{
    switch(relation) {
        case SmtTheoryPrimitiveRelation::EQ_ZERO: return os << "=0";
        case SmtTheoryPrimitiveRelation::GEQ_ZERO: return os << ">=0";
        case SmtTheoryPrimitiveRelation::GT_ZERO: return os << ">0";
        default: ARIADNE_FAIL_MSG("Unknown SMT primitive theory relation");
    }
}

OutputStream& operator<<(OutputStream& os, SmtTheoryWeakRelation relation)
{
    switch(relation) {
        case SmtTheoryWeakRelation::ABS_LEQ_EPSILON: return os << "abs<=epsilon";
        case SmtTheoryWeakRelation::GEQ_MINUS_EPSILON: return os << ">=-epsilon";
        case SmtTheoryWeakRelation::GT_MINUS_EPSILON: return os << ">-epsilon";
        default: ARIADNE_FAIL_MSG("Unknown SMT weak theory relation");
    }
}

} // namespace Ariadne
