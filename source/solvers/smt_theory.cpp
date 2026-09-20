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

namespace {

SmtTheoryRelation make_smt_theory_relation(OperatorCode code)
{
    switch(code) {
        case OperatorCode::EQ: return SmtTheoryRelation::EQ;
        case OperatorCode::NEQ: return SmtTheoryRelation::NEQ;
        case OperatorCode::LEQ: return SmtTheoryRelation::LEQ;
        case OperatorCode::GEQ: return SmtTheoryRelation::GEQ;
        case OperatorCode::LT: return SmtTheoryRelation::LT;
        case OperatorCode::GT: return SmtTheoryRelation::GT;
        default:
            ARIADNE_FAIL_MSG("Expected real comparison operator, got "<<code);
    }
}

SmtTheoryLiteral make_smt_theory_literal_node(BinaryExpressionNode<Kleenean,Real,Real> const& node)
{
    return SmtTheoryLiteral(
        node.arg1(),
        make_smt_theory_relation(node.op().code()),
        node.arg2());
}

template<class E>
SmtTheoryLiteral make_smt_theory_literal_node(E const&)
{
    ARIADNE_FAIL_MSG("Expected binary real comparison node");
}

} // namespace

SmtTheoryLiteral make_smt_theory_literal(ContinuousPredicate const& predicate)
{
    return predicate.node_ref().accept(
        [](auto const& node) { return make_smt_theory_literal_node(node); });
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
