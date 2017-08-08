/***************************************************************************
 *            expression_set.hpp
 *
 *  Copyright 2011-17  Pieter Collins
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

/*! \file expression_set.hpp
 *  \brief Sets defined using expressions over real variables.
 */

#ifndef ARIADNE_EXPRESSION_SET_HPP
#define ARIADNE_EXPRESSION_SET_HPP

#include <iostream>

#include "utility/container.hpp"
#include "utility/declarations.hpp"
#include "expression/variables.hpp"
#include "expression/expression.hpp"
#include "expression/assignment.hpp"
#include "expression/space.hpp"
#include "numeric/float.hpp"
#include "numeric/real.hpp"
#include "function/projection.hpp"
#include "geometry/box.hpp"

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

template<class UB> class VariableInterval;
template<class IVL> class VariablesBox;
template<class S> class ExpressionSet;

typedef VariableInterval<Real> RealVariableInterval;
typedef VariablesBox<RealInterval> RealVariablesBox;
typedef VariablesBox<ExactIntervalType> ExactVariablesBoxType;

class RealExpressionConstraintSet;
class RealExpressionBoundedConstraintSet;

Set<Identifier> arguments(const List<ContinuousPredicate>& c);

template<class C, class P> decltype(auto) any(C const& c, P const& p) {
    typedef decltype(p(*c.begin())) R;
    R r=false; for(auto e:c) { r=r or p(e); if(definitely(r)) { return r; } } return r;
}
template<class C, class P> decltype(auto) all(C const& c, P const& p) {
    typedef decltype(p(*c.begin())) R;
    R r=true; for(auto e:c) { r=r and p(e); if(definitely(not r)) { return r; } } return r;
}

//! \ingroup ExpressionSetSubModule
//! \brief An interval range for a real variable.
template<class UB> class VariableInterval {
    typedef decltype(-declval<UB>()) LB;
  public:
    typedef LB LowerBoundType;
    typedef UB UpperBoundType;
    typedef Interval<UB> IntervalType;
  public:
    template<class U,EnableIf<IsConstructible<UB,U>> =dummy> VariableInterval(VariableInterval<U> const& eivl)
        : VariableInterval(eivl.variable(), Interval<UB>(eivl.interval())) { }
    VariableInterval(const RealVariable& v, const IntervalType& ivl) : _variable(v), _ivl(ivl) { }
    VariableInterval(const LowerBoundType& l, const Variable<Real>& v, const UpperBoundType& u) : _variable(v), _ivl(l,u) { }
    RealVariable const& variable() const { return this->_variable; }
    const IntervalType interval() const { return this->_ivl; }
    const LowerBoundType lower() const { return this->_ivl.lower(); }
    const UpperBoundType upper() const { return this->_ivl.upper(); }
    friend OutputStream& operator<<(OutputStream& os, const VariableInterval<UB>& eivl) {
        return os << eivl._variable << ".in(" << eivl._ivl << "\n"; }
  private:
    RealVariable _variable;
    IntervalType _ivl;
};

//! \ingroup ExpressionSetSubModule
//! \brief An interval range for a real variable.
template<class T> template<class XL, class XU> inline VariableInterval<XU> Variable<T>::in(const XL& l, const XU& u) {
    //static_assert(IsSame<XL,Real>::value,"Can only make box in Real variables.");
    ARIADNE_FAIL_MESSAGE("Can't create interval in non-real variable "<<*this<<"\n");
    assert(false);
}

template<> template<class XL, class XU> inline VariableInterval<XU> Variable<Real>::in(const XL& l, const XU& u) {
    return VariableInterval<XU>(l,*this,u);
}

template<class UB> struct VariableLowerInterval {
    typedef NegationType<UB> LB;
    RealVariable _variable; LB _lower;
    template<class U,EnableIf<IsConstructible<UB,U>> =dummy> VariableLowerInterval(VariableLowerInterval<U> const& lv)
        : VariableLowerInterval<UB>(UB(lv.lower()),lv.variable()) { }
    VariableLowerInterval(const LB& l, const RealVariable& v) : _variable(v),  _lower(l) { }
    operator VariableInterval<UB>() const { return VariableInterval<UB>(_lower,_variable,+infinity); }
    LB lower() const { return _lower; }
    const RealVariable& variable() const { return _variable; }
    operator Expression<Kleenean>() const {
        static_assert(IsSame<UB,Real>::value,""); return ( Real(this->_lower) <= RealExpression(this->_variable) ); };
};

template<class UB> struct VariableUpperInterval {
    typedef NegationType<UB> LB;
    RealVariable _variable; UB _upper;
    VariableUpperInterval<UB>(const RealVariable& v, const UB& u) : _variable(v), _upper(u)  { }
    operator VariableInterval<UB>() const { return VariableInterval<UB>(-infinity,_variable,_upper); }
    const RealVariable& variable() const { return _variable; }
    UB upper() const { return _upper; }
    operator Expression<Kleenean>() const {
        static_assert(IsSame<UB,Real>::value,""); return ( RealExpression(this->_variable) <= Real(this->_upper) ); }
};

typedef VariableLowerInterval<Real> RealVariableLowerInterval;
typedef VariableUpperInterval<Real> RealVariableUpperInterval;

template<class X> using PosType = decltype(+declval<X>());
template<class X> using NegType = decltype(-declval<X>());

template<class UB> inline VariableInterval<UB> operator<=(const VariableLowerInterval<UB>& lv, const PosType<UB>& u) {
    return VariableInterval<UB>(lv.lower(),lv.variable(),u); }
template<class UB> inline VariableInterval<UB> operator>=(const PosType<UB>& u, const VariableLowerInterval<UB>& lv) {
    return VariableInterval<UB>(lv.lower(),lv.variable(),u); }
template<class UB> inline VariableInterval<UB> operator<=(const NegType<UB>& l, const VariableUpperInterval<UB>& vu) {
    return VariableInterval<UB>(l,vu.variable(),vu.upper()); }
template<class UB> inline VariableInterval<UB> operator>=(const VariableUpperInterval<UB>& vu, const NegType<UB>& l) {
    return VariableInterval<UB>(l,vu.variable(),vu.upper()); }

inline VariableLowerInterval<Real> operator<=(const Real& l, const RealVariable& v) {
    return VariableLowerInterval<Real>(l,v); }
inline VariableLowerInterval<Real> operator>=(const RealVariable& v, const Real& l) {
    return VariableLowerInterval<Real>(l,v); }
inline VariableUpperInterval<Real> operator<=(const RealVariable& v, const Real& u) {
    return VariableUpperInterval<Real>(v,u); }
inline VariableUpperInterval<Real> operator>=(const Real& u, const RealVariable& v) {
    return VariableUpperInterval<Real>(v,u); }
inline VariableInterval<Real> operator==(const RealVariable& v, const Real& x) {
    return VariableInterval<Real>(x,v,x); }
inline VariableInterval<Real> operator==(const Real& x, const RealVariable& v) {
    return VariableInterval<Real>(x,v,x); }
inline VariableInterval<Real> operator<=(const VariableLowerInterval<Real>& lv, const Real& u) {
    return VariableInterval<Real>(lv.lower(),lv.variable(),u); }

inline VariableLowerInterval<FloatDPValue> operator<=(const FloatDPValue& l, const RealVariable& v) { return VariableLowerInterval<FloatDPValue>(l,v); }
inline VariableLowerInterval<Dyadic> operator<=(const Dyadic& l, const RealVariable& v) { return VariableLowerInterval<Dyadic>(l,v); }
inline VariableLowerInterval<Dyadic> operator<=(const int& l, const RealVariable& v) { return Dyadic(l)<=v; }

inline VariableLowerInterval<FloatDP> operator<=(const FloatDP& l, const RealVariable& v) { return VariableLowerInterval<FloatDP>(l,v); }
inline VariableLowerInterval<FloatDP> operator<=(const double& l, const RealVariable& v) { return VariableLowerInterval<FloatDP>(l,v); }

//! \brief An box defining ranges for a collection of real variables.
template<class IVL> class VariablesBox {
    typedef typename IVL::UpperBoundType UB;
    Map<RealVariable,IVL> _bnds;
  public:
    typedef IVL IntervalType;
    typedef VariableInterval<UB> VariableIntervalType;
    typedef Box<IVL> BoxType;
    typedef VariablesBox<IVL> VariablesBoxType;

    VariablesBox() : _bnds() { }
    VariablesBox(const RealSpace& spc, const Box<IVL>& bx);
    VariablesBox(const Map<RealVariable,IntervalType>& bnds) : _bnds(bnds) { }
    VariablesBox(const List<VariableIntervalType>& bnds) :  _bnds() {
        for(auto bnd : bnds) { this->_bnds.insert(bnd.variable(),bnd.interval()); } }
    VariablesBox(const InitializerList<VariableIntervalType>& lst) : VariablesBox(List<VariableIntervalType>(lst)) { }
    template<class I,EnableIf<IsConstructible<IVL,I>> =dummy> VariablesBox(VariablesBox<I> const& ebx) : _bnds(ebx.bounds()) { }
    Map<RealVariable,IntervalType> bounds() const { return this->_bnds; }
    Set<RealVariable> variables() const { return this->_bnds.keys(); }
    const IntervalType& operator[](const RealVariable& v) const { return this->_bnds[v]; }
    BoxType euclidean_set(const RealSpace& spc) const {
        BoxType res(spc.dimension()); for(SizeType i=0; i!=res.dimension(); ++i) { res[i]=this->_bnds[spc[i]]; } return res; }
    decltype(auto) is_empty() const { return any(_bnds,[](auto e){return e.second.is_empty();}); }
    friend OutputStream& operator<<(OutputStream& os, const VariablesBoxType& ebx) {
        return os << "VariablesBox"<<class_name<IVL>()<<"( bounds=" << ebx.bounds() << " )"; }
    explicit operator ExpressionSet<Box<IVL>> () const;
};

template<class IVL> VariablesBox<IVL>::VariablesBox(const RealSpace& spc, const Box<IVL>& bx) {
    ARIADNE_PRECONDITION(spc.dimension()==bx.dimension());
    for(SizeType i=0; i!=spc.dimension(); ++i) {
        this->_bnds.insert(spc.variable(i),bx[i]);
    }
}

template<class T> template<class IVL> inline VariablesBox<IVL> Variables<T>::in(const List<IVL>& bx) const {
    static_assert(IsSame<T,Real>::value,"Can only make box in Real variables.");
    assert(false);
}

template<> template<class IVL> inline VariablesBox<IVL> Variables<Real>::in(const List<IVL>& bx) const {
    ARIADNE_PRECONDITION(this->size()==bx.size());
    Map<RealVariable,IVL> bnds;
    for(SizeType i=0; i!=this->size(); ++i) { bnds.insert((*this)[i],bx[i]); }
    return std::move(bnds);
}

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

    friend RealExpressionBoundedConstraintSet intersection(RealVariablesBox const&, RealExpressionConstraintSet const&);

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
    RealExpressionBoundedConstraintSet(const Map<RealVariable,RealInterval>& domain, const List<ContinuousPredicate>& constraints);
    RealExpressionBoundedConstraintSet(const RealVariablesBox& box) : _bounds(box.bounds()) { }
    Set<RealVariable> variables() const { return this->_bounds.keys(); }
    Map<RealVariable,RealInterval> bounds() const { return this->_bounds; }
    List<ContinuousPredicate> const& constraints() const { return this->_constraints; }
    BoundedConstraintSet euclidean_set(const RealSpace& space) const;
    friend OutputStream& operator<<(OutputStream& os, const RealExpressionBoundedConstraintSet& eset);
};


//! \brief An box defining ranges for a collection of real variables.
template<class S> class ExpressionSet {
    RealSpace _spc;
    S _set;
  public:
    typedef S EuclideanSetType;

    ExpressionSet(const RealSpace& spc, const EuclideanSetType& set) : _spc(spc), _set(set) { ARIADNE_ASSERT(spc.dimension()==set.dimension()); }
    Set<RealVariable> variables() const { return _spc.variables(); }
    RealSpace const& space() const { return this->_spc; }
    RealSpace canonical_space() const; [[deprecated]]
    EuclideanSetType const& continuous_set() const { return this->_set; }
    EuclideanSetType const& euclidean_set() const { return this->_set; }
    EuclideanSetType& euclidean_set() { return this->_set; }
    EuclideanSetType euclidean_set(const RealSpace& vars) const {
        Array<SizeType> prj(vars.dimension()); for(SizeType i=0; i!=vars.dimension(); ++i) { prj[i]=this->_spc.index(vars.variable(i)); } return project(this->_set,prj); }
    friend OutputStream& operator<<(OutputStream& os, const ExpressionSet<S>& eset) {
        return os << "ExpressionSet( space=" << eset.space() << ", set=" << eset.continuous_set() << " )"; }
};

template<class S> EqualsType<S> operator==(ExpressionSet<S> const& eset1, ExpressionSet<S> const& eset2) {
    if(eset1.variables()!=eset2.variables()) { return false; }
    ARIADNE_ASSERT(eset1.space()==eset2.space());
    return eset1.euclidean_set()==eset2.euclidean_set();
}

template<class IVL> VariablesBox<IVL>::operator ExpressionSet<Box<IVL>>() const {
    RealSpace spc(List<RealVariable>(this->variables()));
    return ExpressionSet<Box<IVL>>(spc,this->euclidean_set(spc));
}

} // namespace Ariadne



#endif
