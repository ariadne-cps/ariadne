/***************************************************************************
 *            solvers/smt_theory.hpp
 *
 *  Copyright  2026  Luca Geretti
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

#include <vector>

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

//! \brief Primitive relation used by the delta-SMT theory layer.
enum class SmtTheoryPrimitiveRelation {
    EQ_ZERO,
    GEQ_ZERO,
    GT_ZERO
};

//! \brief Primitive theory literal f(x) op 0.
class SmtTheoryPrimitiveLiteral {
  public:
    SmtTheoryPrimitiveLiteral(RealExpression expression, SmtTheoryPrimitiveRelation relation)
        : _expression(std::move(expression)), _relation(relation) { }

    RealExpression const& expression() const { return _expression; }
    SmtTheoryPrimitiveRelation relation() const { return _relation; }

  private:
    RealExpression _expression;
    SmtTheoryPrimitiveRelation _relation;
};

//! \brief One disjunctive alternative obtained when normalizing a theory literal.
using SmtTheoryAlternative = std::vector<SmtTheoryPrimitiveLiteral>;

//! \brief Disjunction of normalized alternatives.
//!
//! All relations except NEQ yield one alternative containing one primitive.
//! NEQ yields two alternatives, corresponding to f<0 or f>0.
using SmtTheoryAlternatives = std::vector<SmtTheoryAlternative>;

//! \brief Relation after epsilon weakening of a primitive theory literal.
enum class SmtTheoryWeakRelation {
    ABS_LEQ_EPSILON,
    GEQ_MINUS_EPSILON,
    GT_MINUS_EPSILON
};

//! \brief Epsilon-weakened primitive preserving strict inequalities.
class SmtTheoryWeakLiteral {
  public:
    SmtTheoryWeakLiteral(RealExpression expression, SmtTheoryWeakRelation relation, ExactDouble epsilon)
        : _expression(std::move(expression)), _relation(relation), _epsilon(epsilon) { }

    RealExpression const& expression() const { return _expression; }
    SmtTheoryWeakRelation relation() const { return _relation; }
    ExactDouble epsilon() const { return _epsilon; }

  private:
    RealExpression _expression;
    SmtTheoryWeakRelation _relation;
    ExactDouble _epsilon;
};

//! \brief Convert a comparison node of a continuous predicate into a theory literal.
//! \pre The predicate root is one of EQ, NEQ, LEQ, GEQ, LT or GT.
SmtTheoryLiteral make_smt_theory_literal(ContinuousPredicate const& predicate);

//! \brief Normalize a literal to primitives f=0, f>=0 and f>0.
//!
//! Disequality is expanded into two alternatives, as required by the standard
//! delta-SMT normalization.
SmtTheoryAlternatives normalize_smt_theory_literal(SmtTheoryLiteral const& literal);

//! \brief Apply epsilon weakening to a normalized primitive.
//! \pre epsilon is strictly positive.
SmtTheoryWeakLiteral weaken_smt_theory_literal(SmtTheoryPrimitiveLiteral const& literal,
                                               ExactDouble epsilon);

OutputStream& operator<<(OutputStream& os, SmtTheoryRelation relation);
OutputStream& operator<<(OutputStream& os, SmtTheoryPrimitiveRelation relation);
OutputStream& operator<<(OutputStream& os, SmtTheoryWeakRelation relation);

} // namespace Ariadne

#endif // ARIADNE_SMT_THEORY_HPP
