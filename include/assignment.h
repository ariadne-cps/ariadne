/***************************************************************************
 *            assignment.h
 *
 *  Copyright 2008-9  Pieter Collins
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


/*! \file assignment.h
 *  \brief Assignment expressions
 */

#ifndef ARIADNE_ASSIGNMENT_H
#define ARIADNE_ASSIGNMENT_H

#include <cstdarg>
#include <iostream>
#include <string>


#include "macros.h"
#include "pointer.h"
#include "container.h"
#include "stlio.h"

#include "numeric.h"

#include "variables.h"
#include "expression.h"

namespace Ariadne {

class Integer;
class Real;

//! \brief An ASCII string.
typedef std::string String;

class StateSpace;

template<class T> class Constant;
template<class T> class Variable;
template<class T> class Expression;
template<class T> class Space;

typedef Constant<Real> RealConstant;

typedef Variable<String> StringVariable;
typedef Variable<Integer> IntegerVariable;
typedef Variable<Real> RealVariable;

typedef Space<Real> RealSpace;

typedef Expression<String> StringExpression;
typedef Expression<Integer> IntegerExpression;
typedef Expression<Real> RealExpression;

typedef Expression<bool> DiscretePredicate;
typedef Expression<tribool> ContinuousPredicate;

class Identifier;
typedef Set<UntypedVariable> VariableSet;
typedef Set< Variable<Real> > RealVariableSet;


template<class T> class Expression;







//! \brief An assignment statement.
template<class LHS, class RHS>
struct Assignment
{
    template<class XLHS> explicit Assignment(const Assignment<XLHS,RHS>& a)
        : lhs(a.lhs), rhs(a.rhs) { }
    Assignment(const LHS& l, const RHS& r) : lhs(l), rhs(r) { }
    operator List< Assignment<LHS,RHS> >() const { return List< Assignment<LHS,RHS> >(1u,*this); }
    const LHS& variable() const { return this->lhs; }
    const RHS& expression() const { return this->rhs; }
    const LHS& left_hand_size() const { return this->lhs; }
    const RHS& right_hand_size() const { return this->rhs; }
    LHS lhs; RHS rhs;
};

template<class LHS, class RHS> bool operator<(const Assignment<LHS,RHS>& a1, const Assignment<LHS,RHS>& a2) {
    return a1.lhs < a2.lhs;
}

template<class LHS, class RHS> inline std::ostream& operator<<(std::ostream& os, const Assignment<LHS,RHS>& a) {
    return os<<a.lhs<<"="<<a.rhs;
}

template<class T> inline Assignment< Variable<T>, Expression<T> >
Variable<T>::operator=(const T& val) const {
    return Assignment< Variable<T>, Expression<T> >(*this,Expression<T>(val)); }

template<class T> inline Assignment< Variable<T>, Expression<T> >
Variable<T>::operator=(const Constant<T>& cnst) const {
    return Assignment< Variable<T>, Expression<T> >(*this,Expression<T>(cnst)); }

template<class T> inline Assignment< Variable<T>, Expression<T> >
Variable<T>::operator=(const Variable<T>& var) const {
    return Assignment< Variable<T>, Expression<T> >(*this,Expression<T>(var)); }

template<class T> inline Assignment< Variable<T>, Expression<T> >
Variable<T>::operator=(const Expression<T>& expr) const {
    return Assignment< Variable<T>, Expression<T> >(*this,Expression<T>(expr)); }

template<class T> template<class D> inline typename EnableIfRealDouble<T,D,Assignment< Variable<T>, Expression<T> > >::Type
Variable<T>::operator=(D c) const {
    return this->operator=(Real(c)); }



template<class T> inline Assignment< PrimedVariable<T>, Expression<T> >
PrimedVariable<T>::operator=(const T& val) const {
    return Assignment< PrimedVariable<T>, Expression<T> >(*this,Expression<T>(val)); }

template<class T> inline Assignment< PrimedVariable<T>, Expression<T> >
PrimedVariable<T>::operator=(const Constant<T>& cnst) const {
    return Assignment< PrimedVariable<T>, Expression<T> >(*this,Expression<T>(cnst)); }

template<class T> inline Assignment< PrimedVariable<T>, Expression<T> >
PrimedVariable<T>::operator=(const Variable<T>& var) const {
    return Assignment< PrimedVariable<T>, Expression<T> >(*this,Expression<T>(var)); }

template<class T> inline Assignment< PrimedVariable<T>, Expression<T> >
PrimedVariable<T>::operator=(const Expression<T>& expr) const {
    return Assignment< PrimedVariable<T>, Expression<T> >(*this,Expression<T>(expr)); }

template<class T> template<class D> inline typename EnableIfRealDouble<T,D,Assignment< PrimedVariable<T>, Expression<T> > >::Type
PrimedVariable<T>::operator=(D c) const {
    return this->operator=(Real(c)); }


template<class T> inline Assignment< DottedVariable<T>, Expression<T> >
DottedVariable<T>::operator=(const T& val) const {
    return Assignment< DottedVariable<T>, Expression<T> >(*this,Expression<T>(val)); }

template<class T> inline Assignment< DottedVariable<T>, Expression<T> >
DottedVariable<T>::operator=(const Constant<T>& cnst) const {
    return Assignment< DottedVariable<T>, Expression<T> >(*this,Expression<T>(cnst)); }

template<class T> inline Assignment< DottedVariable<T>, Expression<T> >
DottedVariable<T>::operator=(const Variable<T>& var) const {
    return Assignment< DottedVariable<T>, Expression<T> >(*this,Expression<T>(var)); }

template<class T> inline Assignment< DottedVariable<T>, Expression<T> >
DottedVariable<T>::operator=(const Expression<T>& expr) const {
    return Assignment< DottedVariable<T>, Expression<T> >(*this,Expression<T>(expr)); }

template<class T> template<class D> inline typename EnableIfRealDouble<T,D,Assignment< DottedVariable<T>, Expression<T> > >::Type
DottedVariable<T>::operator=(D c) const {
    return this->operator=(Real(c)); }

template<class LHS, class RHS> List<Identifier> left_hand_sides(const List<Assignment<LHS,RHS> >& assignments) {
    List<Identifier> result;
    result.reserve(assignments.size());
    for(typename List< Assignment<LHS,RHS> >::const_iterator assignment_iter=assignments.begin();
        assignment_iter!=assignments.end(); ++assignment_iter)
    {
        result.append(assignment_iter->lhs.base().name());
    }
    return result;
}

// Simplifying typedefs
typedef Assignment<IntegerVariable,IntegerExpression> IntegerAssignment;
typedef Assignment<PrimedIntegerVariable,IntegerExpression> PrimedIntegerAssignment;
typedef Assignment<StringVariable,StringExpression> StringAssignment;
typedef Assignment<PrimedStringVariable,StringExpression> PrimedStringAssignment;
typedef Assignment<RealVariable,RealExpression> RealAssignment;
typedef Assignment<PrimedRealVariable,RealExpression> PrimedRealAssignment;
typedef Assignment<DottedRealVariable,RealExpression> DottedRealAssignment;

// Deprecated
typedef Assignment<PrimedIntegerVariable,IntegerExpression> IntegerUpdate;
typedef Assignment<PrimedStringVariable,StringExpression> StringUpdate;







} // namespace Ariadne

#endif // ARIADNE_FORMULA_H
