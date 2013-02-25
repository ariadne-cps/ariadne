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

#include "container.h"
#include "variables.h"
#include "expression.h"
#include "assignment.h"
#include "space.h"
#include "float.h"
#include "real.h"
#include "box.h"

namespace Ariadne {

class Real;
class IntervalSet;

class BoxSet;
class ConstraintSet;
class BoundedConstraintSet;

class Box;
class IntervalConstrainedImageSet;

template<class X> class Variable;
typedef Variable<Real> RealVariable;
template<class X> class Expression;
typedef Expression<Real> RealExpression;
template<class X> class Space;
typedef Space<Real> RealSpace;

class RealVariableInterval;
class RealVariablesBox;
class VariablesBox;
class ExpressionConstraintSet;

Set<Identifier> arguments(const List<ContinuousPredicate>& c);

struct RealVariableLowerInterval {
    Real _lower; RealVariable _variable;
    RealVariableLowerInterval(const Real& l, const RealVariable& v) : _lower(l), _variable(v) { }
    Real lower() const { return _lower; }
    const RealVariable& variable() const { return _variable; }
    operator Expression<Tribool>() const { return ( Real(this->_lower) <= RealExpression(this->_variable) ); };
};

struct RealVariableUpperInterval {
    RealVariable _variable; Real _upper;
    RealVariableUpperInterval(const RealVariable& v, const Real& u) : _variable(v), _upper(u)  { }
    const RealVariable& variable() const { return _variable; }
    Real upper() const { return _upper; }
    operator Expression<Tribool>() const { return ( RealExpression(this->_variable) <= Real(this->_upper) ); }
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
        : _lower(l), _variable(v), _upper(u) { ARIADNE_ASSERT_MSG(l<=u,"Interval("<<l<<","<<u<<") not provably nonempty"); }
    RealVariableInterval(const RealVariableLowerInterval& lv)
        : _lower(lv._lower), _variable(lv._variable), _upper(+inf) { }
    RealVariableInterval(const RealVariableUpperInterval& vu)
        : _lower(-inf), _variable(vu._variable), _upper(vu._upper) { }
    Variable<Real> const& variable() const { return this->_variable; }
    const Interval approximate_interval() const;
    const IntervalSet interval() const;
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

inline RealVariableLowerInterval operator<=(double l, const RealVariable& v) { return Real(l)<=v; }
inline RealVariableLowerInterval operator>=(const RealVariable& v, double l) { return v>=Real(l); }
inline RealVariableUpperInterval operator<=(const RealVariable& v, double u) { return v<=Real(u); }
inline RealVariableUpperInterval operator>=(double u, const RealVariable& v) { return Real(u)>=v; }

inline RealVariableInterval operator==(const RealVariable& v, double x) { return v==Real(x); }
inline RealVariableInterval operator==(double x, const RealVariable& v) { return Real(x)==v; }

//! \ingroup ExpressionSetSubModule
//! \brief An box defining ranges for a collection of real variables.
class RealVariablesBox {
    Map<RealVariable,IntervalSet> _bounds;
  public:
    RealVariablesBox() : _bounds() { };
    RealVariablesBox(const List<RealVariableInterval>& lst);
    RealVariablesBox(const Map<RealVariable,IntervalSet>& bnds) : _bounds(bnds) { }
    RealVariablesBox(const RealSpace& spc, const BoxSet& bx);
    Set<RealVariable> variables() const { return _bounds.keys(); }
    RealSpace canonical_space() const { Set<RealVariable> vars=this->variables(); return RealSpace(List<RealVariable>(vars.begin(),vars.end())); }
    Map<RealVariable,IntervalSet> bounds() const { return this->_bounds; }
    const IntervalSet& operator[](const RealVariable& v) const { return this->_bounds[v]; }
    BoxSet box(const RealSpace& spc) const;
    BoxSet euclidean_set(const RealSpace& spc) const;
    friend VariablesBox approximation(const RealVariablesBox& set);
    friend VariablesBox over_approximation(const RealVariablesBox& vbx);
    friend VariablesBox under_approximation(const RealVariablesBox& set);
    friend OutputStream& operator<<(OutputStream& os, const RealVariablesBox& ebx);
};

//! \brief An box defining ranges for a collection of real variables.
class VariablesBox {
    RealSpace _spc;
    Box _bx;
  public:
    VariablesBox() : _spc(), _bx(0) { }
    VariablesBox(const Map<RealVariable,Interval>& bnds) : _spc(make_list(bnds.keys())), _bx(_spc.dimension()) {
        for(uint i=0; i!=this->_bx.dimension(); ++i) { this->_bx[i] = bnds[this->_spc[i]]; } }
    VariablesBox(const List<RealVariableInterval>& bnds) : _spc(), _bx(bnds.size()) {
        for(uint i=0; i!=bnds.size(); ++i) { this->_spc.append(bnds[i].variable()); this->_bx[i]=bnds[i].approximate_interval(); } }
    VariablesBox(const RealSpace& spc, const Box& bx) : _spc(spc), _bx(bx) { ARIADNE_ASSERT(spc.dimension()==bx.dimension()); }
    Set<RealVariable> variables() const { return make_set(_spc.variables()); }
    RealSpace const& space() const { return this->_spc; }
    Box const& continuous_set() const { return this->_bx; }
    Box const& box() const { return this->_bx; }
    const Interval& operator[](const RealVariable& v) const { return this->_bx[this->_spc.index(v)]; }
    Box euclidean_set(const RealSpace& spc) const {
        Box res(spc.dimension()); for(uint i=0; i!=res.dimension(); ++i) { res[i]=this->_bx[this->_spc.index(spc[i])]; } return res; }
    friend OutputStream& operator<<(OutputStream& os, const VariablesBox& ebx) {
        return os << "VariablesBox( space=" << ebx.space() << ", box=" << ebx.box() << " )"; }
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
    friend std::ostream& operator<<(std::ostream& os, const RealExpressionConstraintSet& eset);
};

//! \ingroup ExpressionSetSubModule
//! \brief A set defined as the intersection of a box and the preimage of a box under a continuous function.
//! The set is described as \f$S=D\cap g^{-1}(C)\f$ where \f$D\f$ is the domain, \f$g\f$ the constraint function and \f$C\f$ the codomain.
class RealExpressionBoundedConstraintSet
{
    Map<RealVariable,IntervalSet> _bounds;
    List<ContinuousPredicate> _constraints;
  public:
    RealExpressionBoundedConstraintSet(const List<RealVariableInterval>& domain);
    RealExpressionBoundedConstraintSet(const List<RealVariableInterval>& domain, const List<ContinuousPredicate>& constraints);
    RealExpressionBoundedConstraintSet(const RealVariablesBox& box) : _bounds(box.bounds()) { }
    Set<RealVariable> variables() const { return this->_bounds.keys(); }
    Map<RealVariable,IntervalSet> bounds() const { return this->_bounds; }
    List<ContinuousPredicate> const& constraints() const { return this->_constraints; }
    BoundedConstraintSet euclidean_set(const RealSpace& space) const;
    friend std::ostream& operator<<(std::ostream& os, const RealExpressionBoundedConstraintSet& eset);
};


}



#endif
