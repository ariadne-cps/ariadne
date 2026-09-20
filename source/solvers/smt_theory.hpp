/***************************************************************************
 *            solvers/smt_theory.hpp
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

#ifndef ARIADNE_SMT_THEORY_HPP
#define ARIADNE_SMT_THEORY_HPP

#include "numeric/numeric.hpp"
#include "symbolic/expression.hpp"

namespace Ariadne {

//! \ingroup Solvers
//! \brief Relation of an atomic real SMT literal.
enum class SmtTheoryRelation {
    EQ,
    NEQ,
    LEQ,
    GEQ,
    LT,
    GT
};

//! \ingroup Solvers
//! \brief Atomic real theory literal preserving strictness and disequality.
class SmtTheoryLiteral {
  public:
    SmtTheoryLiteral(RealExpression lhs, SmtTheoryRelation relation, RealExpression rhs)
        : _lhs(std::move(lhs)), _relation(relation), _rhs(std::move(rhs)) { }

    RealExpression const& lhs() const { return _lhs; }
    RealExpression const& rhs() const { return _rhs; }
    SmtTheoryRelation relation() const { return _relation; }

    SmtTheoryLiteral negated() const;

  private:
    RealExpression _lhs;
    SmtTheoryRelation _relation;
    RealExpression _rhs;
};

//! \brief Convert a comparison node of a continuous predicate into a theory literal.
//! \pre The predicate root is one of EQ, NEQ, LEQ, GEQ, LT or GT.
SmtTheoryLiteral make_smt_theory_literal(ContinuousPredicate const& predicate);

OutputStream& operator<<(OutputStream& os, SmtTheoryRelation relation);

} // namespace Ariadne

#endif // ARIADNE_SMT_THEORY_HPP
