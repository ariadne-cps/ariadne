/***************************************************************************
 *            solvers/smt_boolean.cpp
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
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "solvers/smt_boolean.hpp"

#include <map>

#include "utility/exceptions.hpp"
#include "solvers/smt_theory.hpp"

namespace Ariadne {

namespace {

struct CanonicalTheoryAtom {
    RealExpression lhs;
    SmtTheoryRelation relation;
    RealExpression rhs;
};

CanonicalTheoryAtom canonical_theory_atom(ContinuousPredicate const& predicate)
{
    SmtTheoryLiteral literal=make_smt_theory_literal(predicate);
    RealExpression lhs=literal.lhs();
    RealExpression rhs=literal.rhs();
    SmtTheoryRelation relation=literal.relation();

    switch(relation) {
        case SmtTheoryRelation::LEQ:
            return {rhs,SmtTheoryRelation::GEQ,lhs};
        case SmtTheoryRelation::LT:
            return {rhs,SmtTheoryRelation::GT,lhs};
        case SmtTheoryRelation::EQ:
        case SmtTheoryRelation::NEQ:
            if(before(rhs,lhs)) {
                std::swap(lhs,rhs);
            }
            return {lhs,relation,rhs};
        case SmtTheoryRelation::GEQ:
        case SmtTheoryRelation::GT:
            return {lhs,relation,rhs};
        default:
            ARIADNE_FAIL_MSG("Unknown SMT theory relation");
    }
}

Bool equivalent_theory_atoms(ContinuousPredicate const& lhs,
                             ContinuousPredicate const& rhs)
{
    CanonicalTheoryAtom left=canonical_theory_atom(lhs);
    CanonicalTheoryAtom right=canonical_theory_atom(rhs);
    return left.relation==right.relation
        && identical(left.lhs,right.lhs)
        && identical(left.rhs,right.rhs);
}

} // namespace

class TseitinBuilder {
  public:
    SmtBooleanEncoding build(ContinuousPredicate const& predicate) {
        Int root=this->_encode(predicate);
        _encoding._add_clause({root});
        return _encoding;
    }

  private:
    Int _new_variable() {
        return _encoding._new_variable();
    }

    std::optional<Bool> _constant_value(ContinuousPredicate const& predicate) const {
        if(predicate.code()!=OperatorCode::CNST) {
            return std::nullopt;
        }
        Kleenean const& value=predicate.val();
        if(definitely(value)) {
            return true;
        }
        if(definitely(!value)) {
            return false;
        }
        ARIADNE_FAIL_MSG("Indeterminate constant in SMT Boolean encoding");
    }

    Int _encode_constant(ContinuousPredicate const& predicate) {
        Int variable=this->_new_variable();
        Kleenean const& value=predicate.val();
        if(definitely(value)) {
            _encoding._add_clause({variable});
        } else if(definitely(!value)) {
            _encoding._add_clause({-variable});
        } else {
            ARIADNE_FAIL_MSG("Indeterminate constant in SMT Boolean encoding");
        }
        return variable;
    }

    Int _encode_atom(ContinuousPredicate const& predicate) {
        const Void* key=static_cast<const Void*>(predicate.node_raw_ptr());
        auto iter=_atom_variables.find(key);
        if(iter!=_atom_variables.end()) {
            return iter->second;
        }

        for(SizeType i=0u; i!=_encoding.atom_count(); ++i) {
            if(identical(_encoding.atom(i),predicate)
                    || equivalent_theory_atoms(_encoding.atom(i),predicate)) {
                Int variable=static_cast<Int>(_encoding.atom_variable(i));
                _atom_variables.emplace(key,variable);
                return variable;
            }
        }

        Int variable=this->_new_variable();
        _atom_variables.emplace(key,variable);
        _encoding._add_atom(predicate,static_cast<SizeType>(variable));
        return variable;
    }

    Int _encode(ContinuousPredicate const& predicate) {
        switch(predicate.code()) {
            case OperatorCode::CNST:
                return this->_encode_constant(predicate);

            case OperatorCode::SGN:
            case OperatorCode::EQ:
            case OperatorCode::NEQ:
            case OperatorCode::GEQ:
            case OperatorCode::LEQ:
            case OperatorCode::GT:
            case OperatorCode::LT:
                return this->_encode_atom(predicate);

            case OperatorCode::NOT: {
                auto value=this->_constant_value(predicate.arg());
                if(value.has_value()) {
                    return this->_encode_constant(ContinuousPredicate(!*value));
                }
                return -this->_encode(predicate.arg());
            }

            case OperatorCode::AND: {
                auto lhs_value=this->_constant_value(predicate.arg1());
                auto rhs_value=this->_constant_value(predicate.arg2());
                if(lhs_value.has_value()) {
                    return *lhs_value ? this->_encode(predicate.arg2())
                                      : this->_encode_constant(ContinuousPredicate(false));
                }
                if(rhs_value.has_value()) {
                    return *rhs_value ? this->_encode(predicate.arg1())
                                      : this->_encode_constant(ContinuousPredicate(false));
                }

                Int lhs=this->_encode(predicate.arg1());
                Int rhs=this->_encode(predicate.arg2());
                Int variable=this->_new_variable();

                _encoding._add_clause({-variable,lhs});
                _encoding._add_clause({-variable,rhs});
                _encoding._add_clause({variable,-lhs,-rhs});
                return variable;
            }

            case OperatorCode::OR: {
                auto lhs_value=this->_constant_value(predicate.arg1());
                auto rhs_value=this->_constant_value(predicate.arg2());
                if(lhs_value.has_value()) {
                    return *lhs_value ? this->_encode_constant(ContinuousPredicate(true))
                                      : this->_encode(predicate.arg2());
                }
                if(rhs_value.has_value()) {
                    return *rhs_value ? this->_encode_constant(ContinuousPredicate(true))
                                      : this->_encode(predicate.arg1());
                }

                Int lhs=this->_encode(predicate.arg1());
                Int rhs=this->_encode(predicate.arg2());
                Int variable=this->_new_variable();

                _encoding._add_clause({variable,-lhs});
                _encoding._add_clause({variable,-rhs});
                _encoding._add_clause({-variable,lhs,rhs});
                return variable;
            }

            default:
                ARIADNE_FAIL_MSG("Unsupported operator in SMT Boolean encoding: "<<predicate.code());
        }
    }

    SmtBooleanEncoding _encoding;
    std::map<const Void*,Int> _atom_variables;
};

SmtBooleanEncoding SmtBooleanEncoder::encode(ContinuousPredicate const& predicate) const
{
    return TseitinBuilder().build(predicate);
}

} // namespace Ariadne
