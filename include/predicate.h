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
#include "function_interface.h"

#include "macros.h"
#include "pointer.h"

#include "vector.h"
#include "matrix.h"
#include "differential.h"
#include "taylor_model.h"

namespace Ariadne {


class PredicateInterface {
  public:
    typedef unsigned int size_type;

    virtual size_type argument_size() const = 0;
    virtual tribool evaluate(const Vector<Float>& x) const = 0;
    virtual tribool evaluate(const Vector<ValidatedNumberType>& x) const = 0;
};

class ExpressionPredicate
    : public PredicateInterface
{
    friend ExpressionPredicate operator!(const ExpressionPredicate&);
  public:
    ExpressionPredicate(const EffectiveScalarFunction& expression)
        : _expression(expression), _sign(+1) { }
    const EffectiveScalarFunction& expression() const { return _expression; }
    int sign() const { return _sign; }

    bool same(const ExpressionPredicate& ep2) const {
        const ExpressionPredicate& ep1=*this;
        return ep1._expression.pointer()==ep2._expression.pointer() && ep1._sign==ep2._sign; }
    bool opposite(const ExpressionPredicate& ep2) const {
        const ExpressionPredicate& ep1=*this;
        return ep1._expression.pointer()==ep2._expression.pointer() && ep1._sign!=ep2._sign; }
    bool operator<(const ExpressionPredicate& p) const {
        return (_expression.pointer()) < (const void*)(p._expression.pointer()); }
    size_type argument_size() const { return _expression.argument_size(); }
    tribool evaluate(const Vector<Float>& x) const {
        Float value=_expression.evaluate(x)*_sign;
        if(value<0) { return true; }
        else if(value>0) { return false; }
        else { return indeterminate; } }
    tribool evaluate(const Vector<ValidatedNumberType>& x) const {
        Interval range=_expression.evaluate(x)*_sign;
        if(range.upper()<0) { return true; }
        else if(range.lower()>0) { return false; }
        else { return indeterminate; } }
  private:
    EffectiveScalarFunction _expression;
    int _sign;
};

//! \brief A predicate which is a disjunction of simple predicates.
//! An example of a disjunctive predicate is \f$p_1\vee p_2\vee p_3\f$.
class DisjunctivePredicate
    : public PredicateInterface
{
  public:
    size_type size() const { return _predicates.size(); }
    const ExpressionPredicate& operator[](size_type i) const { return _predicates[i]; }

    bool vacuous() const { return _predicates.empty(); }
    bool tautologous() const { return _tautology; }

    DisjunctivePredicate& operator|=(ExpressionPredicate p) {
        if(_tautology) { return *this; }
        for(uint i=0; i!=_predicates.size(); ++i) {
            if(p.same(_predicates[i])) { return *this; }
            if(p.opposite(_predicates[i])) { _predicates.clear(); _tautology=true; return *this; }
        }
        _predicates.push_back(p);
        return *this;
    }

    DisjunctivePredicate& operator|=(const DisjunctivePredicate& p) {
        for(uint i=0; i!=p.size(); ++i) { (*this) |= p[i]; } return *this; }

    virtual size_type argument_size() const;
    virtual tribool evaluate(const Vector<Float>& x) const;
    virtual tribool evaluate(const Vector<ValidatedNumberType>& x) const;

  private:
    std::vector<ExpressionPredicate> _predicates;
    bool _tautology;
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
        for(uint i=0; i!=_cnf.size(); ++i) { _cnf[i] |= p; } return *this; }

    virtual size_type argument_size() const;
    virtual tribool evaluate(const Vector<Float>& x) const;
    virtual tribool evaluate(const Vector<ValidatedNumberType>& x) const;
  private:
    std::vector<DisjunctivePredicate> _cnf;
};

ExpressionPredicate operator!(const ExpressionPredicate& p) {
    ExpressionPredicate res(p); res._sign=-p._sign; return res;
}



} // namespace Ariadne

#endif
