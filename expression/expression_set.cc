/***************************************************************************
 *            expression_set.cc
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

#include "utility/standard.h"
#include "config.h"

#include "expression/expression_set.h"
#include "expression/expression.h"
#include "expression/space.h"
#include "geometry/box.h"
#include "geometry/function_set.h"



namespace Ariadne {

template<class T> struct KeyValue { };

template<> struct KeyValue<RealVariableInterval> {
    typedef RealVariable KeyType; typedef RealInterval ValueType;
    KeyType const& key(const RealVariableInterval& ivl) { return ivl.variable(); }
    ValueType value(const RealVariableInterval& ivl) { return ivl.interval(); }
};

template<class T> Bool unique_keys(const List<T>& lst) {
    typedef typename KeyValue<T>::KeyType K;
    Set<K> keys;
    for(typename List<T>::ConstIterator lst_iter=lst.begin(); lst_iter!=lst.end(); ++lst_iter) {
        const K& k=KeyValue<T>().key(*lst_iter);
        if(keys.contains(k)) {
            return false;
        }
        keys.insert(k);
    }
    return true;
}

template<class T> Map<typename KeyValue<T>::KeyType,typename KeyValue<T>::ValueType>
make_key_value_map(const List<T>& lst) {
    typedef typename KeyValue<T>::KeyType K;
    typedef typename KeyValue<T>::ValueType V;
    Map<K,V> map;
    for(typename List<T>::ConstIterator lst_iter=lst.begin(); lst_iter!=lst.end(); ++lst_iter) {
        map.insert(KeyValue<T>().key(*lst_iter),KeyValue<T>().value(*lst_iter));
    }
    return map;
}

template<class K, class V> List<V> order(const Map<K,V>& map, const List<K>& ordering) {
    List<V> result;
    result.reserve(ordering.size());
    for(Nat i=0; i!=ordering.size(); ++i) {
        result.append(map[ordering[i]]);
    }
    return result;
}

Set<Identifier> variables(const List<RealVariableInterval>& b) {
    Set<Identifier> r;
    for(Nat i=0; i!=b.size(); ++i) { r.insert(b[i].variable().name()); }
    return r;
}

Set<Identifier> arguments(const List<ContinuousPredicate>& c) {
    Set<Identifier> r;
    for(Nat i=0; i!=c.size(); ++i) { r.adjoin(arguments(c[i])); }
    return r;
}


const RealInterval RealVariableInterval::interval() const {
    return RealInterval(this->_lower,this->_upper);
}


OutputStream& operator<<(OutputStream& os, const RealVariableInterval& eivl) {
    return os << eivl.variable() << ".in(" << eivl.lower() << "," << eivl.upper() << ")";
}


RealVariablesBox::RealVariablesBox(const InitializerList<RealVariableInterval>& lst)
    : RealVariablesBox(List<RealVariableInterval>(lst))
{
}

RealVariablesBox::RealVariablesBox(const List<RealVariableInterval>& lst)
{
    for(Nat i=0; i!=lst.size(); ++i) {
        _bounds.insert(lst[i].variable(),lst[i].interval());
    }
}

RealVariablesBox::RealVariablesBox(const RealSpace& spc, const RealBox& bx)
{
    ARIADNE_ASSERT(spc.size()==bx.size());
    for(Nat i=0; i!=spc.size(); ++i) {
        _bounds.insert(spc[i],bx[i]);
    }
}

RealBox RealVariablesBox::euclidean_set(const RealSpace& spc) const {
    RealBox bx(spc.size());
    for(Nat i=0; i!=spc.size(); ++i) {
        bx[i]=(*this)[spc[i]];
    }
    return bx;
}

RealBox RealVariablesBox::box(const RealSpace& spc) const {
    return this->euclidean_set(spc);
}

OutputStream& operator<<(OutputStream& os, const RealVariablesBox& ebx) {
    return os << ebx._bounds;
}

ExactInterval over_approximation(RealInterval ivl) {
    return cast_exact_interval(UpperInterval(ivl));
}

ExactInterval under_approximation(RealInterval ivl) {
    return cast_exact_interval(LowerInterval(ivl));
}

ExactInterval approximation(RealInterval ivl) {
    return cast_exact_interval(ApproximateInterval(ivl));
}


VariablesBox over_approximation(const RealVariablesBox& ebx) {
    Map<RealVariable,ExactInterval> result;
    for(Map<RealVariable,RealInterval>::ConstIterator iter=ebx.bounds().begin();
        iter!=ebx.bounds().end(); ++iter)
    {
        result[iter->first]=over_approximation(iter->second);
    }
    return result;
}

VariablesBox approximation(const RealVariablesBox& ebx) {
    Map<RealVariable,ExactInterval> result;
    for(Map<RealVariable,RealInterval>::ConstIterator iter=ebx.bounds().begin();
        iter!=ebx.bounds().end(); ++iter)
    {
        result[iter->first]=approximation(iter->second);
    }
    return result;
}

VariablesBox under_approximation(const RealVariablesBox& ebx) {
    Map<RealVariable,ExactInterval> result;
    for(Map<RealVariable,RealInterval>::ConstIterator iter=ebx.bounds().begin();
        iter!=ebx.bounds().end(); ++iter)
    {
        result[iter->first]=under_approximation(iter->second);
    }
    return result;
}



RealExpressionConstraintSet::RealExpressionConstraintSet(const List<ContinuousPredicate>& constraints)
    : _constraints(constraints)
{
}

ConstraintSet RealExpressionConstraintSet::euclidean_set(const RealSpace& space) const {
    const RealExpressionConstraintSet& set = *this;
    List<RealConstraint> constraints;
    for(Nat i=0; i!=set.constraints().size(); ++i) {
        RealExpression constraint_expression=indicator(set.constraints()[i],NEGATIVE);
        EffectiveScalarFunction constraint_function( Ariadne::make_function(constraint_expression,space) );
        constraints.append( constraint_function <= Real(0) );
    }
    return ConstraintSet(constraints);
}

OutputStream& operator<<(OutputStream& os, const RealExpressionConstraintSet& eset) {
    os << "[";
    for(List<ContinuousPredicate>::ConstIterator iter=eset._constraints.begin(); iter!=eset._constraints.end(); ++iter) {
        os << (iter==eset._constraints.begin()?"":",") << *iter; }
    return os << "]";
}



RealExpressionBoundedConstraintSet::RealExpressionBoundedConstraintSet(const InitializerList<RealVariableInterval>& bounds)
    : RealExpressionBoundedConstraintSet(List<RealVariableInterval>(bounds))
{
}

RealExpressionBoundedConstraintSet::RealExpressionBoundedConstraintSet(const List<RealVariableInterval>& bounds)
    : _bounds(make_key_value_map(bounds)), _constraints()
{
    ARIADNE_ASSERT(unique_keys(bounds));
}

RealExpressionBoundedConstraintSet::RealExpressionBoundedConstraintSet(const List<RealVariableInterval>& bounds, const List<ContinuousPredicate>& constraints)
    : _bounds(make_key_value_map(bounds)), _constraints(constraints)
{
    ARIADNE_ASSERT(unique_keys(bounds));
    ARIADNE_ASSERT( subset(arguments(constraints),Ariadne::variables(bounds)) );
}

BoundedConstraintSet RealExpressionBoundedConstraintSet::euclidean_set(const RealSpace& space) const {
    const RealExpressionBoundedConstraintSet& set = *this;
    RealBox domain=RealVariablesBox(set.bounds()).euclidean_set(space);
    List<RealConstraint> constraints;
    for(Nat i=0; i!=set.constraints().size(); ++i) {
        RealExpression constraint_expression=indicator(set.constraints()[i],NEGATIVE);
        EffectiveScalarFunction constraint_function( Ariadne::make_function(constraint_expression,space) );
        constraints.append( constraint_function <= Real(0) );
    }
    return BoundedConstraintSet(domain,constraints);}

OutputStream& operator<<(OutputStream& os, const RealExpressionBoundedConstraintSet& eset) {
    os << "[";
    for(Map<RealVariable,RealInterval>::ConstIterator iter=eset._bounds.begin(); iter!=eset._bounds.end(); ++iter) {
        os << (iter==eset._bounds.begin()?"":",") << *iter; }
    os << ";";
    for(List<ContinuousPredicate>::ConstIterator iter=eset._constraints.begin(); iter!=eset._constraints.end(); ++iter) {
        os << (iter==eset._constraints.begin()?"":",") << *iter; }
    return os << "]";
    return os << eset._bounds << eset._constraints;
}

ValidatedConstrainedImageSet approximate_euclidean_set(const RealExpressionBoundedConstraintSet& set, const RealSpace& space) {
    ExactIntervalVector domain=cast_exact_box(ApproximateBox(RealVariablesBox(set.bounds()).euclidean_set(space)));
    ValidatedVectorFunction identity=ValidatedVectorFunction::identity(domain.size());

    ValidatedConstrainedImageSet result(domain,identity);
    //List<ValidatedConstraint> constraints;
    for(Nat i=0; i!=set.constraints().size(); ++i) {
        RealExpression constraint_expression=indicator(set.constraints()[i],NEGATIVE);
        ValidatedScalarFunction constraint_function( Ariadne::make_function(constraint_expression,space) );
        result.new_parameter_constraint(constraint_function <= ExactNumber(0) );
        //constraints.append( constraint_function <= 0.0 );
    }
    return result;
    //return ValidatedConstrainedImageSet(domain,identity,constraints);
}





} // namespace Ariadne
