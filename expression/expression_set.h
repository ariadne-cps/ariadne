/***************************************************************************
 *            expression_set.h
 *
 *  Copyright 2011  Pieter Collins
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

/*! \file expression_set.h
 *  \brief Sets defined using expressions over real variables.
 */

#ifndef ARIADNE_EXPRESSION_SET_H
#define ARIADNE_EXPRESSION_SET_H

#include <iostream>

#include "utility/container.h"
#include "utility/declarations.h"
#include "expression/variables.h"
#include "expression/expression.h"
#include "expression/assignment.h"
#include "expression/space.h"
#include "numeric/float.h"
#include "numeric/real.h"
#include "geometry/box.h"

namespace Ariadne {

class ConstraintSet;
class BoundedConstraintSet;

class ValidatedConstrainedImageSet;

template<class X> class Variable;
typedef Variable<Real> RealVariable;
template<class X> class Expression;
typedef Expression<Real> RealExpression;
template<class X> class Space;
typedef Space<Real> RealSpace;

class RealVariableInterval;
class RealVariablesBox;
class ExactFloat64VariablesBox;
class ExpressionConstraintSet;

typedef ExactFloat64VariablesBox ExactVariablesBoxType;

Set<Identifier> arguments(const List<ContinuousPredicate>& c);

struct RealVariableLowerInterval {
    Real _lower; RealVariable _variable;
    RealVariableLowerInterval(const Real& l, const RealVariable& v) : _lower(l), _variable(v) { }
    Real lower() const { return _lower; }
    const RealVariable& variable() const { return _variable; }
    operator Expression<Kleenean>() const { return ( Real(this->_lower) <= RealExpression(this->_variable) ); };
};

struct RealVariableUpperInterval {
    RealVariable _variable; Real _upper;
    RealVariableUpperInterval(const RealVariable& v, const Real& u) : _variable(v), _upper(u)  { }
    const RealVariable& variable() const { return _variable; }
    Real upper() const { return _upper; }
    operator Expression<Kleenean>() const { return ( RealExpression(this->_variable) <= Real(this->_upper) ); }
};


//! \ingroup ExpressionSetSubModule
//! \brief An interval range for a real variable.
class RealVariableInterval {
  private:
    Real _lower;
    Variable<Real> _variable;
    Real _upper;
  public:
    RealVariableInterval(const Real& l, const Variable<Real>& v, const Real& u)
        : _lower(l), _variable(v), _upper(u) { ARIADNE_ASSERT_MSG(definitely(l<=u),"RealInterval("<<l<<","<<u<<") not provably nonempty"); }
    RealVariableInterval(const RealVariableLowerInterval& lv)
        : _lower(lv._lower), _variable(lv._variable), _upper(+infty) { }
    RealVariableInterval(const RealVariableUpperInterval& vu)
        : _lower(-infty), _variable(vu._variable), _upper(vu._upper) { }
    Variable<Real> const& variable() const { return this->_variable; }
    const RealInterval interval() const;
    const Real lower() const { return this->_lower; }
    const Real upper() const { return this->_upper; }
};

OutputStream& operator<<(OutputStream& os, const RealVariableInterval& eivl);

template<class T> template<class XL, class XU> inline RealVariableInterval Variable<T>::in(const XL& l, const XU& u) {
    //ARIADNE_FAIL_MESSAGE("Can't create interval in variable "<<*this<<" of type "<<name<T>()<<"\n");
    ARIADNE_FAIL_MESSAGE("Can't create interval in variable "<<*this<<"\n");
    return RealVariableInterval(l,*this,u);
}

template<> template<class XL, class XU> inline RealVariableInterval Variable<Real>::in(const XL& l, const XU& u) {
    return RealVariableInterval(l,*this,u);
}

inline RealVariableInterval operator<=(const RealVariableLowerInterval& lv, const Real& u) {
    return RealVariableInterval(lv.lower(),lv.variable(),u); }
inline RealVariableInterval operator>=(const Real& u, const RealVariableLowerInterval& lv) {
    return RealVariableInterval(lv.lower(),lv.variable(),u); }
inline RealVariableInterval operator<=(const Real& l, const RealVariableUpperInterval& vu) {
    return RealVariableInterval(l,vu.variable(),vu.upper()); }
inline RealVariableInterval operator>=(const RealVariableUpperInterval& vu, const Real& l) {
    return RealVariableInterval(l,vu.variable(),vu.upper()); }

inline RealVariableLowerInterval operator<=(const Real& l, const RealVariable& v) {
    return RealVariableLowerInterval(l,v); }
inline RealVariableLowerInterval operator>=(const RealVariable& v, const Real& l) {
    return RealVariableLowerInterval(l,v); }
inline RealVariableUpperInterval operator<=(const RealVariable& v, const Real& u) {
    return RealVariableUpperInterval(v,u); }
inline RealVariableUpperInterval operator>=(const Real& u, const RealVariable& v) {
    return RealVariableUpperInterval(v,u); }

inline RealVariableInterval operator==(const RealVariable& v, const Real& x) {
    return RealVariableInterval(x,v,x); }
inline RealVariableInterval operator==(const Real& x, const RealVariable& v) {
    return RealVariableInterval(x,v,x); }

inline RealVariableLowerInterval operator<=(Int l, RealVariable const& v) {
    return RealVariableLowerInterval(l,v); }

//! \ingroup ExpressionSetSubModule
//! \brief An box defining ranges for a collection of real variables.
class RealVariablesBox {
    Map<RealVariable,RealInterval> _bounds;
  public:
    RealVariablesBox() : _bounds() { };
    RealVariablesBox(const InitializerList<RealVariableInterval>& lst);
    RealVariablesBox(const List<RealVariableInterval>& lst);
    RealVariablesBox(const Map<RealVariable,RealInterval>& bnds) : _bounds(bnds) { }
    RealVariablesBox(const RealSpace& spc, const RealBox& bx);
    Set<RealVariable> variables() const { return _bounds.keys(); }
    RealSpace canonical_space() const { Set<RealVariable> vars=this->variables(); return RealSpace(List<RealVariable>(vars.begin(),vars.end())); }
    Map<RealVariable,RealInterval> bounds() const { return this->_bounds; }
    const RealInterval& operator[](const RealVariable& v) const { return this->_bounds[v]; }
    RealBox box(const RealSpace& spc) const;
    RealBox euclidean_set(const RealSpace& spc) const;
    friend ExactFloat64VariablesBox approximation(const RealVariablesBox& set);
    friend ExactFloat64VariablesBox over_approximation(const RealVariablesBox& vbx);
    friend ExactFloat64VariablesBox under_approximation(const RealVariablesBox& set);
    friend OutputStream& operator<<(OutputStream& os, const RealVariablesBox& ebx);
};


template<class T> template<class IVL> inline RealVariablesBox Variables<T>::in(const List<IVL>& bx) const {
    static_assert(IsSame<T,Real>::value,"Can only make box in Real variables.");
    assert(false);
}

template<> template<class IVL> inline RealVariablesBox Variables<Real>::in(const List<IVL>& bx) const {
    ARIADNE_PRECONDITION(this->size()==bx.size());
    Map<RealVariable,RealInterval> bnds;
    for(SizeType i=0; i!=this->size(); ++i) { bnds.insert((*this)[i],bx[i]); }
    return std::move(bnds);
}

//! \ingroup ExpressionSetSubModule
//! \brief An interval range for a real variable.
class ExactFloat64VariableInterval {
  private:
    Float64Value _lower;
    Variable<Real> _variable;
    Float64Value _upper;
  public:
    ExactFloat64VariableInterval(const Variable<Real>& v, const Float64ExactInterval& ivl)
        : _lower(ivl.lower()), _variable(v), _upper(ivl.upper()) { }
    ExactFloat64VariableInterval(const Float64Value& l, const Variable<Real>& v, const Float64Value& u)
        : _lower(l), _variable(v), _upper(u) { ARIADNE_ASSERT_MSG(definitely(l<=u),"Float64Value("<<l<<","<<u<<") not provably nonempty"); }
    Variable<Real> const& variable() const { return this->_variable; }
    const Float64ExactInterval interval() const { return Float64ExactInterval(this->_lower,this->_upper); }
    const Float64Value lower() const { return this->_lower; }
    const Float64Value upper() const { return this->_upper; }
    friend OutputStream& operator<<(OutputStream& os, const ExactFloat64VariableInterval& eivl);
};

//! \brief An box defining ranges for a collection of real variables.
class ExactFloat64VariablesBox {
    RealSpace _spc;
    Float64ExactBox _bx;
  public:
    ExactFloat64VariablesBox() : _spc(), _bx(0) { }
    ExactFloat64VariablesBox(const Map<RealVariable,Float64ExactInterval>& bnds) : _spc(List<RealVariable>(bnds.keys())), _bx(_spc.dimension()) {
        for(Nat i=0; i!=this->_bx.dimension(); ++i) { this->_bx[i] = bnds[this->_spc[i]]; } }
    ExactFloat64VariablesBox(const List<RealVariableInterval>& bnds) : _spc(), _bx(bnds.size()) {
        for(Nat i=0; i!=bnds.size(); ++i) {
            this->_spc.append(bnds[i].variable()); this->_bx[i]=cast_exact_interval(ApproximateIntervalType(bnds[i].interval())); } }
    ExactFloat64VariablesBox(const RealSpace& spc, const Float64ExactBox& bx) : _spc(spc), _bx(bx) { ARIADNE_ASSERT(spc.dimension()==bx.dimension()); }
    Map<RealVariable,Float64ExactInterval> bounds() const { Map<RealVariable,Float64ExactInterval> bnds;
        for(SizeType i=0; i!=_spc.size(); ++i) { bnds.insert(_spc[i],_bx[i]); } return bnds; }
    Set<RealVariable> variables() const { return Set<RealVariable>(_spc.variables()); }
    RealSpace const& space() const { return this->_spc; }
    Float64ExactBox const& continuous_set() const { return this->_bx; }
    Float64ExactBox const& box() const { return this->_bx; }
    const Float64ExactInterval& operator[](const RealVariable& v) const { return this->_bx[this->_spc.index(v)]; }
    Float64ExactBox euclidean_set(const RealSpace& spc) const {
        Float64ExactBox res(spc.dimension()); for(Nat i=0; i!=res.dimension(); ++i) { res[i]=this->_bx[this->_spc.index(spc[i])]; } return res; }
    friend OutputStream& operator<<(OutputStream& os, const ExactFloat64VariablesBox& ebx) {
        return os << "ExactFloat64VariablesBox( space=" << ebx.space() << ", box=" << ebx.box() << " )"; }
};



//! \ingroup ExpressionSetSubModule
//! \brief A set defined as the preimage of a box under a continuous function.
//! The set is described as \f$S=g^{-1}(C)\f$ where \f$g\f$ the constraint function and \f$C\f$ the codomain.
class RealExpressionConstraintSet
{
    List<ContinuousPredicate> _constraints;
  public:
    RealExpressionConstraintSet();
    RealExpressionConstraintSet(const List<ContinuousPredicate>& constraints);
    Set<RealVariable> variables() const { return Set<RealVariable>(arguments(this->_constraints)); }
    List<ContinuousPredicate> const& constraints() const { return this->_constraints; }
    ConstraintSet euclidean_set(const RealSpace& space) const;
    friend OutputStream& operator<<(OutputStream& os, const RealExpressionConstraintSet& eset);
};

//! \ingroup ExpressionSetSubModule
//! \brief A set defined as the intersection of a box and the preimage of a box under a continuous function.
//! The set is described as \f$S=D\cap g^{-1}(C)\f$ where \f$D\f$ is the domain, \f$g\f$ the constraint function and \f$C\f$ the codomain.
class RealExpressionBoundedConstraintSet
{
    Map<RealVariable,RealInterval> _bounds;
    List<ContinuousPredicate> _constraints;
  public:
    RealExpressionBoundedConstraintSet(const InitializerList<RealVariableInterval>& domain);
    RealExpressionBoundedConstraintSet(const List<RealVariableInterval>& domain);
    RealExpressionBoundedConstraintSet(const List<RealVariableInterval>& domain, const List<ContinuousPredicate>& constraints);
    RealExpressionBoundedConstraintSet(const RealVariablesBox& box) : _bounds(box.bounds()) { }
    Set<RealVariable> variables() const { return this->_bounds.keys(); }
    Map<RealVariable,RealInterval> bounds() const { return this->_bounds; }
    List<ContinuousPredicate> const& constraints() const { return this->_constraints; }
    BoundedConstraintSet euclidean_set(const RealSpace& space) const;
    friend OutputStream& operator<<(OutputStream& os, const RealExpressionBoundedConstraintSet& eset);
};


}



#endif
