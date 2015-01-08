/***************************************************************************
 *            predicate.h
 *
 *  Copyright 2008  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file predicate.h
 *  \brief Predicates.
 */

#ifndef ARIADNE_PREDICATE_H
#define ARIADNE_PREDICATE_H

#include <cstdarg>
#include <iosfwd>
#include <iostream>
#include "function/function_interface.h"

#include "utility/macros.h"
#include "utility/pointer.h"

#include "algebra/vector.h"
#include "algebra/matrix.h"
#include "algebra/differential.h"
#include "function/taylor_model.h"

namespace Ariadne {


class PredicateInterface {
  public:
    typedef std::size_t SizeType;

    virtual SizeType argument_size() const = 0;
    virtual Tribool evaluate(const Vector<Float>& x) const = 0;
    virtual Tribool evaluate(const Vector<ValidatedNumber>& x) const = 0;
};

class ExpressionPredicate
    : public PredicateInterface
{
    friend ExpressionPredicate operator!(const ExpressionPredicate&);
  public:
    ExpressionPredicate(const EffectiveScalarFunction& expression)
        : _expression(expression), _sign(+1) { }
    const EffectiveScalarFunction& expression() const { return _expression; }
    Int sign() const { return _sign; }

    Bool same(const ExpressionPredicate& ep2) const {
        const ExpressionPredicate& ep1=*this;
        return ep1._expression.pointer()==ep2._expression.pointer() && ep1._sign==ep2._sign; }
    Bool opposite(const ExpressionPredicate& ep2) const {
        const ExpressionPredicate& ep1=*this;
        return ep1._expression.pointer()==ep2._expression.pointer() && ep1._sign!=ep2._sign; }
    Bool operator<(const ExpressionPredicate& p) const {
        return (_expression.pointer()) < (const Void*)(p._expression.pointer()); }
    SizeType argument_size() const { return _expression.argument_size(); }
    Tribool evaluate(const Vector<Float>& x) const {
        Float value=_expression.evaluate(x)*_sign;
        if(value<0) { return true; }
        else if(value>0) { return false; }
        else { return indeterminate; } }
    Tribool evaluate(const Vector<ValidatedNumber>& x) const {
        ExactInterval range=_expression.evaluate(x)*_sign;
        if(range.upper()<0) { return true; }
        else if(range.lower()>0) { return false; }
        else { return indeterminate; } }
  private:
    EffectiveScalarFunction _expression;
    Int _sign;
};

//! \brief A predicate which is a disjunction of simple predicates.
//! An example of a disjunctive predicate is \f$p_1\vee p_2\vee p_3\f$.
class DisjunctivePredicate
    : public PredicateInterface
{
  public:
    SizeType size() const { return _predicates.size(); }
    const ExpressionPredicate& operator[](SizeType i) const { return _predicates[i]; }

    Bool vacuous() const { return _predicates.empty(); }
    Bool tautologous() const { return _tautology; }

    DisjunctivePredicate& operator|=(ExpressionPredicate p) {
        if(_tautology) { return *this; }
        for(Nat i=0; i!=_predicates.size(); ++i) {
            if(p.same(_predicates[i])) { return *this; }
            if(p.opposite(_predicates[i])) { _predicates.clear(); _tautology=true; return *this; }
        }
        _predicates.push_back(p);
        return *this;
    }

    DisjunctivePredicate& operator|=(const DisjunctivePredicate& p) {
        for(Nat i=0; i!=p.size(); ++i) { (*this) |= p[i]; } return *this; }

    virtual SizeType argument_size() const;
    virtual Tribool evaluate(const Vector<Float>& x) const;
    virtual Tribool evaluate(const Vector<ValidatedNumber>& x) const;

  private:
    std::vector<ExpressionPredicate> _predicates;
    Bool _tautology;
};

//! \brief A predicate in <em>conjunctive normal form</em>.
//!
//! An example of conjunctive normal form is \f$(p_1\vee p_2)\wedge(p_2\vee p_3\vee p_4)\f$.
class ConjunctiveNormalFormPredicate
    : public PredicateInterface
{
    friend ConjunctiveNormalFormPredicate
        operator&&(const ConjunctiveNormalFormPredicate& p1, const ConjunctiveNormalFormPredicate& p2);
    friend ConjunctiveNormalFormPredicate
        operator||(const ConjunctiveNormalFormPredicate& p1, const ConjunctiveNormalFormPredicate& p2);
  public:
    ConjunctiveNormalFormPredicate& operator&=(DisjunctivePredicate p) {
        _cnf.push_back(p); return *this; }
    ConjunctiveNormalFormPredicate& operator|=(DisjunctivePredicate p) {
        for(Nat i=0; i!=_cnf.size(); ++i) { _cnf[i] |= p; } return *this; }

    virtual SizeType argument_size() const;
    virtual Tribool evaluate(const Vector<Float>& x) const;
    virtual Tribool evaluate(const Vector<ValidatedNumber>& x) const;
  private:
    std::vector<DisjunctivePredicate> _cnf;
};

ExpressionPredicate operator!(const ExpressionPredicate& p) {
    ExpressionPredicate res(p); res._sign=-p._sign; return res;
}



} // namespace Ariadne

#endif
