/***************************************************************************
 *            assignment.hpp
 *
 *  Copyright 2008-17 Pieter Collins
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


/*! \file assignment.hpp
 *  \brief Assignment expressions
 */

#ifndef ARIADNE_ASSIGNMENT_HPP
#define ARIADNE_ASSIGNMENT_HPP

#include <cstdarg>
#include <iostream>
#include <string>


#include "utility/macros.hpp"
#include "utility/pointer.hpp"
#include "utility/container.hpp"
#include "utility/stlio.hpp"
#include "utility/string.hpp"

#include "numeric/numeric.hpp"

#include "expression/variables.hpp"
#include "expression/expression.hpp"

namespace Ariadne {

class Integer;
class Real;
class String;

class StateSpace;

template<class T> class Constant;
template<class T> class Variable;
template<class T> class Expression;
template<class T> class Space;

typedef Space<Real> RealSpace;

typedef Expression<String> StringExpression;
typedef Expression<Integer> IntegerExpression;
typedef Expression<Real> RealExpression;

class Identifier;
typedef Set<UntypedVariable> VariableSet;
typedef Set<Variable<Real>> RealVariableSet;

// Simplifying typedefs
typedef Assignment<StringVariable,StringExpression> StringAssignment;
typedef Assignment<PrimedStringVariable,StringExpression> PrimedStringAssignment;
typedef Assignment<IntegerVariable,IntegerExpression> IntegerAssignment;
typedef Assignment<PrimedIntegerVariable,IntegerExpression> PrimedIntegerAssignment;
typedef Assignment<RealVariable,RealExpression> RealAssignment;
typedef Assignment<PrimedRealVariable,RealExpression> PrimedRealAssignment;
typedef Assignment<DottedRealVariable,RealExpression> DottedRealAssignment;

typedef List<PrimedStringAssignment> PrimedStringAssignments;
typedef List<RealAssignment> RealAssignments;
typedef List<PrimedRealAssignment> PrimedRealAssignments;
typedef List<DottedRealAssignment> DottedRealAssignments;

typedef Assignment<RealVariable,Real> RealConstantAssignment;


template<class F, class T> List<ResultOf<F(T)>> zip(F const& f, List<T> const& l) {
    typedef ResultOf<F(T)> R;
    List<R> r; r.reserve(l.size());
    for(SizeType i=0; i!=l.size(); ++i) { r.append(f(l[i])); }
    return r;
}

template<class F, class T1, class T2> List<ResultOf<F(T1,T2)>> zip(F const& f, List<T1> const& l1, List<T2> const& l2) {
    typedef ResultOf<F(T1,T2)> R;
    assert(l1.size()==l2.size());
    List<R> r; r.reserve(l1.size());
    for(SizeType i=0; i!=l1.size(); ++i) { r.append(f(l1[i],l2[i])); }
    return r;
}








//! \ingroup ExpressionModule
//! \brief An assignment statement.
template<class LHS, class RHS>
class Assignment
{
  public:
    Assignment(const LHS& l, const RHS& r) : lhs(l), rhs(r) { }
    const LHS& variable() const { return this->lhs; }
    const RHS& expression() const { return this->rhs; }
    const LHS& left_hand_side() const { return this->lhs; }
    const RHS& right_hand_side() const { return this->rhs; }
    LHS lhs; RHS rhs;
};
template<class LHS, class RHS> Bool operator<(const Assignment<LHS,RHS>& a1, const Assignment<LHS,RHS>& a2) {
    return a1.lhs < a2.lhs;
}
template<class LHS, class RHS> inline OutputStream& operator<<(OutputStream& os, const Assignment<LHS,RHS>& a) {
    return os<<a.lhs<<"="<<a.rhs;
}

template<class LHS, class RHS> typename LHS::BaseType left_hand_side(const Assignment<LHS,RHS>& assignment) {
    return assignment.lhs.base();
}
template<class LHS, class RHS> List<typename LHS::BaseType> left_hand_sides(const List<Assignment<LHS,RHS>>& assignments) {
    return zip([](Assignment<LHS,RHS>const&a){return a.lhs.base();},assignments);
}


template<class T>
class Assignment<Variable<T>,T>
{
  public:
    Assignment(const Variable<T>& l, const T& r) : lhs(l), rhs(r) { }
    const Variable<T>& variable() const { return this->lhs; }
    const T& value() const { return this->rhs; }
    const Expression<T>& expression() const { return Expression<T>(this->rhs); }
    const Variable<T>& left_hand_side() const { return this->lhs; }
    const T& right_hand_side() const { return this->rhs; }
    inline operator Valuation<T,T> () const;
    Variable<T> lhs; T rhs;
};
template<class T> inline Assignment<Variable<T>,T>
Variable<T>::operator=(const T& val) const {
    return Assignment<Variable<T>,T>(*this,val);
}



template<class T> inline Assignment<Variable<T>,Expression<T>>
LetVariable<T>::operator=(const Expression<T>& expr) const {
    return Assignment<Variable<T>,Expression<T>>(this->base(),Expression<T>(expr)); }
template<class T> inline Assignment<Variable<T>,Expression<T>>
LetVariable<T>::operator=(const Variable<T>& var) const {
    return this->operator=(Expression<T>(var)); }
template<class T> inline Assignment<Variable<T>,Expression<T>>
LetVariable<T>::operator=(const T& cnst) const {
    return this->operator=(Expression<T>(cnst)); }

template<class T> inline Assignment<PrimedVariable<T>,Expression<T>>
PrimedVariable<T>::operator=(const Expression<T>& expr) const {
    return Assignment<PrimedVariable<T>,Expression<T>>(*this,Expression<T>(expr)); }
template<class T> inline Assignment<PrimedVariable<T>,Expression<T>>
PrimedVariable<T>::operator=(const T& cnst) const {
    return this->operator=(Expression<T>(cnst)); }

template<class T> inline Assignment<DottedVariable<T>,Expression<T>>
DottedVariable<T>::operator=(const Expression<T>& expr) const {
    return Assignment<DottedVariable<T>,Expression<T>>(*this,expr); }
template<class T> inline Assignment<DottedVariable<T>,Expression<T>>
DottedVariable<T>::operator=(const T& cnst) const {
    return this->operator=(Expression<T>(cnst)); }


template<class T> inline List<Assignment<Variable<T>,Expression<T>>> LetVariables<T>::operator=(const List<Expression<T>>& rhs) {
    return zip([](Variable<T>const&l,Expression<T>const&r){return let(l)=r;},this->_lhs,rhs);
}

template<class T> inline List<Assignment<PrimedVariable<T>,Expression<T>>> PrimedVariables<T>::operator=(const List<Expression<T>>& rhs) {
    return zip([](Variable<T>const&l,Expression<T>const&r){return prime(l)=r;},this->_lhs,rhs);
}

template<class T> inline List<Assignment<DottedVariable<T>,Expression<T>>> DottedVariables<T>::operator=(const List<Expression<T>>& rhs) {
    return zip([](Variable<T>const&l,Expression<T>const&r){return dot(l)=r;},this->_lhs,rhs);
}


} // namespace Ariadne

#include "expression/valuation.hpp"
namespace Ariadne {
template<class T> inline Assignment<Variable<T>,T>::operator Valuation<T> () const { Valuation<T> r; r.insert(this->lhs,this->rhs); return r; }
} // namespace Ariadne


#endif // ARIADNE_ASSIGNMENT_HPP
