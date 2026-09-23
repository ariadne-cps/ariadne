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

[[noreturn]] Void fail_indeterminate_smt_constant()
{
    throw std::runtime_error("Indeterminate constant in SMT Boolean encoding");
}

} // namespace

namespace SmtBooleanTestSupport {

SmtTheoryRelation canonical_relation(SmtTheoryRelation relation)
{
    switch(relation) {
        case SmtTheoryRelation::LEQ: return SmtTheoryRelation::GEQ;
        case SmtTheoryRelation::LT: return SmtTheoryRelation::GT;
        case SmtTheoryRelation::EQ:
        case SmtTheoryRelation::NEQ:
        case SmtTheoryRelation::GEQ:
        case SmtTheoryRelation::GT:
            return relation;
        default:
            throw std::runtime_error("Unknown SMT theory relation");
    }
}

Bool supported_operator(OperatorCode code)
{
    switch(code) {
        case OperatorCode::CNST:
        case OperatorCode::SGN:
        case OperatorCode::EQ:
        case OperatorCode::NEQ:
        case OperatorCode::GEQ:
        case OperatorCode::LEQ:
        case OperatorCode::GT:
        case OperatorCode::LT:
        case OperatorCode::NOT:
        case OperatorCode::AND:
        case OperatorCode::OR:
            return true;
        default:
            throw std::runtime_error("Unsupported operator in SMT Boolean encoding");
    }
}

} // namespace SmtBooleanTestSupport

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

    SmtTheoryRelation canonical=
        SmtBooleanTestSupport::canonical_relation(relation);
    if(relation==SmtTheoryRelation::LEQ || relation==SmtTheoryRelation::LT) {
        return {rhs,canonical,lhs};
    }
    if(relation==SmtTheoryRelation::EQ || relation==SmtTheoryRelation::NEQ) {
        if(before(rhs,lhs)) {
            std::swap(lhs,rhs);
        }
    }
    return {lhs,canonical,rhs};
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
        std::optional<Bool> result;
        if(definitely(value)) {
            result=true;
        } else if(definitely(!value)) {
            result=false;
        } else {
            fail_indeterminate_smt_constant();
        }
        return result;
    }

    Int _encode_constant(ContinuousPredicate const& predicate) {
        Int variable=this->_new_variable();
        Kleenean const& value=predicate.val();
        if(not definitely(value) && not definitely(!value)) {
            fail_indeterminate_smt_constant();
        }
        _encoding._add_clause({definitely(value) ? variable : -variable});
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
        OperatorCode code=predicate.code();
        SmtBooleanTestSupport::supported_operator(code);

        if(code==OperatorCode::CNST) {
            return this->_encode_constant(predicate);
        }

        if(code==OperatorCode::SGN
           || code==OperatorCode::EQ
           || code==OperatorCode::NEQ
           || code==OperatorCode::GEQ
           || code==OperatorCode::LEQ
           || code==OperatorCode::GT
           || code==OperatorCode::LT) {
            return this->_encode_atom(predicate);
        }

        if(code==OperatorCode::NOT) {
            auto value=this->_constant_value(predicate.arg());
            if(value.has_value()) {
                return this->_encode_constant(ContinuousPredicate(!*value));
            }
            return -this->_encode(predicate.arg());
        }

        if(code==OperatorCode::AND) {
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

    SmtBooleanEncoding _encoding;
    std::map<const Void*,Int> _atom_variables;
};

SmtBooleanEncoding SmtBooleanEncoder::encode(ContinuousPredicate const& predicate) const
{
    return TseitinBuilder().build(predicate);
}

} // namespace Ariadne
