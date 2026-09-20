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

} // namespace Ariadne
