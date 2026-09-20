/***************************************************************************
 *            solvers/smt_boolean.hpp
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
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file solvers/smt_boolean.hpp
 *  \brief Tseitin encoding of continuous predicates for SMT solving.
 */

#ifndef ARIADNE_SMT_BOOLEAN_HPP
#define ARIADNE_SMT_BOOLEAN_HPP

#include <vector>

#include "numeric/numeric.hpp"
#include "symbolic/expression.hpp"

namespace Ariadne {

//! \ingroup Solvers
//! \brief A CNF encoding of a continuous predicate.
//!
//! SAT variables are numbered from one. Positive integers denote positive
//! literals and negative integers denote negated literals.
class SmtBooleanEncoding {
  public:
    using Literal = Int;
    using Clause = std::vector<Literal>;

    std::vector<Clause> const& clauses() const { return _clauses; }
    SizeType variable_count() const { return _variable_count; }
    SizeType atom_count() const { return _atoms.size(); }

    ContinuousPredicate const& atom(SizeType i) const { return _atoms[i]; }
    SizeType atom_variable(SizeType i) const { return _atom_variables[i]; }

  private:
    friend class SmtBooleanEncoder;

    Void _add_clause(Clause clause) { _clauses.push_back(std::move(clause)); }
    Int _new_variable() {
        ++_variable_count;
        return static_cast<Int>(_variable_count);
    }
    Void _add_atom(ContinuousPredicate const& atom, SizeType variable) {
        _atoms.push_back(atom);
        _atom_variables.push_back(variable);
    }

    std::vector<Clause> _clauses;
    std::vector<ContinuousPredicate> _atoms;
    std::vector<SizeType> _atom_variables;
    SizeType _variable_count = 0u;
};

//! \ingroup Solvers
//! \brief Encode AND/OR/NOT combinations of real comparisons into CNF.
class SmtBooleanEncoder {
  public:
    SmtBooleanEncoding encode(ContinuousPredicate const& predicate) const;
};

} // namespace Ariadne

#endif // ARIADNE_SMT_BOOLEAN_HPP
