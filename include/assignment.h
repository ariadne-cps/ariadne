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

typedef List<RealExpression> RealExpressions;

typedef Expression<bool> DiscretePredicate;
typedef Expression<tribool> ContinuousPredicate;

class Identifier;
typedef Set<UntypedVariable> VariableSet;
typedef Set< Variable<Real> > RealVariableSet;









//! \ingroup ExpressionModule
//! \brief An assignment statement.
template<class LHS, class RHS>
struct Assignment
{
//    template<class XLHS> explicit Assignment(const Assignment<XLHS,RHS>& a)
//        : lhs(a.lhs), rhs(a.rhs) { }
    Assignment(const LHS& l, const RHS& r) : lhs(l), rhs(r) { }
//    operator List< Assignment<LHS,RHS> >() const { return List< Assignment<LHS,RHS> >(1u,*this); }
    const LHS& variable() const { return this->lhs; }
    const RHS& expression() const { return this->rhs; }
    const LHS& left_hand_side() const { return this->lhs; }
    const RHS& right_hand_side() const { return this->rhs; }
    LHS lhs; RHS rhs;
};

template<class T>
struct Assignment< Variable<T>, T>
{
    Assignment(const Variable<T>& l, const T& r) : lhs(l), rhs(r) { }
    operator Assignment< Variable<T>,Expression<T> > () const { return Assignment<Variable<T>,Expression<T> >(this->lhs,this->rhs); }
    operator List< Assignment<Variable<T>,T > >() const { return List< Assignment<Variable<T>,T > >(1u,*this); }
    operator List< Assignment<Variable<T>, Expression<T> > >() const { return List< Assignment<Variable<T>,Expression<T> > >(1u,*this); }
    const Variable<T>& variable() const { return this->lhs; }
    const T& value() const { return this->rhs; }
    const Expression<T>& expression() const { return Expression<T>(this->rhs); }
    const Variable<T>& left_hand_side() const { return this->lhs; }
    const T& right_hand_side() const { return this->rhs; }
    inline operator Valuation<T,T> () const;
    Variable<T> lhs; T rhs;
};

template<class LHS, class RHS> bool operator<(const Assignment<LHS,RHS>& a1, const Assignment<LHS,RHS>& a2) {
    return a1.lhs < a2.lhs;
}

template<class LHS, class RHS> inline std::ostream& operator<<(std::ostream& os, const Assignment<LHS,RHS>& a) {
    return os<<a.lhs<<"="<<a.rhs;
}

//template<class T> inline Assignment< Variable<T>, Expression<T> >
//Variable<T>::operator=(const T& val) const {
//    return Assignment< Variable<T>, Expression<T> >(*this,Expression<T>(val)); }

template<class T> inline Assignment< Variable<T>, T >
Variable<T>::operator=(const T& val) const {
    return Assignment< Variable<T>, T >(*this,val); }

template<class T> inline Assignment< Variable<T>, Expression<T> >
Variable<T>::operator=(const Constant<T>& cnst) const {
    return Assignment< Variable<T>, Expression<T> >(*this,Expression<T>(cnst)); }

template<class T> inline Assignment< Variable<T>, Expression<T> >
Variable<T>::operator=(const Variable<T>& var) const {
    return Assignment< Variable<T>, Expression<T> >(*this,Expression<T>(var)); }

template<class T> inline Assignment< Variable<T>, Expression<T> >
Variable<T>::operator=(const Expression<T>& expr) const {
    return Assignment< Variable<T>, Expression<T> >(*this,Expression<T>(expr)); }

template<class T> template<class D> inline typename EnableIfRealBuiltin<T,D,Assignment< Variable<T>, T > >::Type
Variable<T>::operator=(D c) const {
    return this->operator=(Real(c)); }


template<class T> inline Assignment< Variable<T>, Expression<T> >
LetVariable<T>::operator=(const T& c) const {
    return Assignment< Variable<T>, Expression<T> >(this->base(),Expression<T>(c)); }

template<class T> template<class D> inline typename EnableIfRealBuiltin<T,D,Assignment< Variable<T>, Expression<T> > >::Type
LetVariable<T>::operator=(D c) const {
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

template<class T> template<class D> inline typename EnableIfRealBuiltin<T,D,Assignment< PrimedVariable<T>, Expression<T> > >::Type
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

template<class T> template<class D> inline typename EnableIfRealBuiltin<T,D,Assignment< DottedVariable<T>, Expression<T> > >::Type
DottedVariable<T>::operator=(D c) const {
    return this->operator=(Real(c)); }



template<class T> struct Let {
    const List< Variable<T> >& _lhs;
    Let(const List< Variable<T> >& lhs) : _lhs(lhs) { }
    List< Assignment<Variable<T>, Expression<T> > > operator=(const List<Expression<T> >&);
};
template<class T> inline Let<T> let(const List<Variable<T> >& lhs) { return Let<T>(lhs); }
template<class T> inline List< Assignment<Variable<T>, Expression<T> > > Let<T>::operator=(const List<Expression<T> >& rhs) {
    assert(this->_lhs.size()==rhs.size());
    List< Assignment<Variable<T>,Expression<T> > > result;
    for(uint i=0; i!=rhs.size(); ++i) { result.append(let(this->_lhs[i])=rhs[i]); }
    return result;
}

template<class T> struct Dot {
    const List< Variable<T> >& _lhs;
    Dot(const List< Variable<T> >& lhs) : _lhs(lhs) { }
    List< Assignment<DottedVariable<T>, Expression<T> > > operator=(const List<Expression<T> >&);
};
template<class T> inline Dot<T> dot(const List<Variable<T> >& lhs) { return Dot<T>(lhs); }
template<class T> inline List< Assignment<DottedVariable<T>, Expression<T> > > Dot<T>::operator=(const List<Expression<T> >& rhs) {
    assert(this->_lhs.size()==rhs.size());
    List< Assignment<DottedVariable<T>,Expression<T> > > result;
    for(uint i=0; i!=rhs.size(); ++i) { result.append(dot(this->_lhs[i])=rhs[i]); }
    return result;
}

template<> struct Dot<void> {
    const List<Identifier>& _lhs;
    Dot(const List<Identifier>& lhs) : _lhs(lhs) { }
    template<class T> List< Assignment<DottedVariable<T>, Expression<T> > > operator=(const List<Expression<T> >&);
};
inline Dot<void> dot(const List<String>& lhs) { return Dot<void>(lhs); }
template<class T> inline List< Assignment<DottedVariable<T>, Expression<T> > > Dot<void>::operator=(const List<Expression<T> >& rhs) {
    List< Assignment<DottedVariable<T>,Expression<T> > > result;
    for(uint i=0; i!=rhs.size(); ++i) { result.append(DottedVariable<T>(this->_lhs[i])=rhs[i]); }
    return result;
}

template<class T> struct Primed {
    const List< Variable<T> >& _lhs;
    Primed(const List< Variable<T> >& lhs) : _lhs(lhs) { }
    List< Assignment<PrimedVariable<T>, Expression<T> > > operator=(const List<Expression<T> >&);
    List< Assignment<PrimedVariable<T>, Expression<T> > > operator=(const List<Variable<T> >&);
};
template<class T> inline Primed<T> next(const List<Variable<T> >& lhs) { return Primed<T>(lhs); }
template<class T> inline List< Assignment<PrimedVariable<T>, Expression<T> > > Primed<T>::operator=(const List<Expression<T> >& rhs) {
    assert(this->_lhs.size()==rhs.size());
    List< Assignment<PrimedVariable<T>,Expression<T> > > result;
    for(uint i=0; i!=rhs.size(); ++i) { result.append(next(this->_lhs[i])=rhs[i]); }
    return result;
}
template<class T> inline List< Assignment<PrimedVariable<T>, Expression<T> > > Primed<T>::operator=(const List<Variable<T> >& rhs) {
    return this->operator=(List<Expression<T> >(rhs)); }

template<class LHS, class RHS> List<typename LHS::BaseType> left_hand_sides(const List<Assignment<LHS,RHS> >& assignments) {
    List<typename LHS::BaseType> result;
    result.reserve(assignments.size());
    for(typename List< Assignment<LHS,RHS> >::const_iterator assignment_iter=assignments.begin();
        assignment_iter!=assignments.end(); ++assignment_iter)
    {
        result.append(assignment_iter->lhs.base());
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

typedef List<PrimedStringAssignment> PrimedStringAssignments;
typedef List<RealAssignment> RealAssignments;
typedef List<DottedRealAssignment> DottedRealAssignments;
typedef List<PrimedRealAssignment> PrimedRealAssignments;

typedef Assignment<RealVariable,Real> RealConstantAssignment;

inline List<RealAssignment> operator,(const List<RealAssignment>& la, const RealConstantAssignment& ac) { return operator,(la,RealAssignment(ac)); }
inline List<RealAssignment> operator,(const List<RealConstantAssignment>& lac, const RealAssignment& a) { return operator,(List<RealAssignment>(lac),a); }
inline List<RealAssignment> operator,(const RealConstantAssignment& ac, const RealAssignment& a) { return operator,(RealAssignment(ac),a); }
inline List<RealAssignment> operator,(const RealAssignment& a, const RealConstantAssignment& ac) { return operator,(a,RealAssignment(ac)); }

inline List<RealConstantAssignment> operator,(const List<RealConstantAssignment>& lac, const RealConstantAssignment& ac) {
    List<RealConstantAssignment> r(lac); r.append(ac); return r; }
inline List<RealConstantAssignment> operator,(const RealConstantAssignment& ac1, const RealConstantAssignment& ac2) {
    List<RealConstantAssignment> r; r.append(ac1); r.append(ac2); return r; }

} // namespace Ariadne

#include "valuation.h"
namespace Ariadne {
template<class T> inline Assignment< Variable<T>, T>::operator Valuation<T> () const { Valuation<T> r; r.insert(this->lhs,this->rhs); return r; }
} // namespace Ariadne


#endif // ARIADNE_ASSIGNMENT_H
