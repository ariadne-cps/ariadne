/***************************************************************************
 *            expression.decl.hpp
 *
 *  Copyright 2008-17 Pieter Collins
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

/*! \file expression.decl.hpp
 *  \brief Declarations of expression classes and type aliases.
 */

#ifndef ARIADNE_EXPRESSION_DECL_HPP
#define ARIADNE_EXPRESSION_DECL_HPP

namespace Ariadne {

class String;
class Integer;
class Real;

class Identifier;

template<class T> class Constant;
template<class T> class Variable;
template<class T> class Variables;
template<class T> class Space;
template<class T> class Expression;

template<class T> class LetVariable;
template<class T> class DottedVariable;
template<class T> class PrimedVariable;
template<class LHS,class RHS> class Assignment;

using BooleanVariable = Variable<Boolean>;
using KleeneanVariable = Variable<Kleenean>;
using StringVariable = Variable<String>;
using IntegerVariable = Variable<Integer>;
using RealVariable = Variable<Real>;
using RealVariables = Variables<Real>;

using PrimedStringVariable = PrimedVariable<String>;
using LetIntegerVariable = LetVariable<Integer>;
using PrimedIntegerVariable = PrimedVariable<Integer>;
using LetRealVariable = LetVariable<Real>;
using PrimedRealVariable = PrimedVariable<Real>;
using DottedRealVariable = DottedVariable<Real>;

using RealSpace = Space<Real>;

template<class T, class X=T> class Valuation;
using StringValuation = Valuation<String>;
using IntegerValuation = Valuation<Integer>;
class DiscreteValuation;
template<class X> class ContinuousValuation;

using DiscretePredicate = Expression<Boolean>;
using ContinuousPredicate = Expression<Kleenean>;
using BooleanExpression = Expression<Boolean>;
using KleeneanExpression = Expression<Kleenean>;
using StringExpression = Expression<String>;
using IntegerExpression = Expression<Integer>;
using RealExpression = Expression<Real>;
using RealExpressions = List<Expression<Real>>;

using StringAssignment = Assignment<StringVariable,StringExpression>;
using PrimedStringAssignment = Assignment<PrimedStringVariable,StringExpression>;
using IntegerAssignment = Assignment<IntegerVariable,IntegerExpression>;
using PrimedIntegerAssignment = Assignment<PrimedIntegerVariable,IntegerExpression>;
using RealAssignment = Assignment<RealVariable,RealExpression>;
using PrimedRealAssignment = Assignment<PrimedRealVariable,RealExpression>;
using DottedRealAssignment = Assignment<DottedRealVariable,RealExpression>;

using PrimedStringAssignments = List<PrimedStringAssignment>;
using RealAssignments = List<RealAssignment>;
using PrimedRealAssignments = List<PrimedRealAssignment>;
using DottedRealAssignments = List<DottedRealAssignment>;

using RealConstantAssignment = Assignment<RealVariable,Real>;



template<class UB> class Interval;
template<class UB> class VariableInterval;
template<class IVL> class VariablesBox;
using RealInterval = Interval<Real>;
using RealVariableInterval = VariableInterval<Real>;
using RealVariableIntervals = List<RealVariableInterval>;
using RealVariablesBox = VariablesBox<RealInterval>;


} // namespace Ariadne

#endif /* ARIADNE_EXPRESSION_DECL_HPP */
