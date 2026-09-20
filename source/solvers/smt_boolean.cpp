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

namespace Ariadne {

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

    Int _encode_atom(ContinuousPredicate const& predicate) {
        const Void* key=static_cast<const Void*>(predicate.node_raw_ptr());
        auto iter=_atom_variables.find(key);
        if(iter!=_atom_variables.end()) {
            return iter->second;
        }

        Int variable=this->_new_variable();
        _atom_variables.emplace(key,variable);
        _encoding._add_atom(predicate,static_cast<SizeType>(variable));
        return variable;
    }

    Int _encode(ContinuousPredicate const& predicate) {
        switch(predicate.code()) {
            case OperatorCode::EQ:
            case OperatorCode::NEQ:
            case OperatorCode::GEQ:
            case OperatorCode::LEQ:
            case OperatorCode::GT:
            case OperatorCode::LT:
                return this->_encode_atom(predicate);

            case OperatorCode::NOT:
                return -this->_encode(predicate.arg());

            case OperatorCode::AND: {
                Int lhs=this->_encode(predicate.arg1());
                Int rhs=this->_encode(predicate.arg2());
                Int variable=this->_new_variable();

                _encoding._add_clause({-variable,lhs});
                _encoding._add_clause({-variable,rhs});
                _encoding._add_clause({variable,-lhs,-rhs});
                return variable;
            }

            case OperatorCode::OR: {
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
