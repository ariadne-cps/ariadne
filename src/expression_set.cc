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

#include <iostream>

#include "config.h"

#include "expression_set.h"
#include "expression.h"
#include "space.h"
#include "box.h"
#include "function_set.h"



namespace Ariadne {

template<class T> struct KeyValue { };

template<> struct KeyValue<RealVariableInterval> {
    typedef RealVariable KeyType; typedef RealIntervalSet ValueType;
    KeyType const& key(const RealVariableInterval& ivl) { return ivl.variable(); }
    ValueType value(const RealVariableInterval& ivl) { return ivl.interval(); }
};

template<class T> bool unique_keys(const List<T>& lst) {
    typedef typename KeyValue<T>::KeyType K;
    Set<K> keys;
    for(typename List<T>::const_iterator lst_iter=lst.begin(); lst_iter!=lst.end(); ++lst_iter) {
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
    for(typename List<T>::const_iterator lst_iter=lst.begin(); lst_iter!=lst.end(); ++lst_iter) {
        map.insert(KeyValue<T>().key(*lst_iter),KeyValue<T>().value(*lst_iter));
    }
    return map;
}

template<class K, class V> List<V> order(const Map<K,V>& map, const List<K>& ordering) {
    List<V> result;
    result.reserve(ordering.size());
    for(uint i=0; i!=ordering.size(); ++i) {
        result.append(map[ordering[i]]);
    }
    return result;
}

Set<Identifier> variables(const List<RealVariableInterval>& b) {
    Set<Identifier> r;
    for(uint i=0; i!=b.size(); ++i) { r.insert(b[i].variable().name()); }
    return r;
}

Set<Identifier> arguments(const List<ContinuousPredicate>& c) {
    Set<Identifier> r;
    for(uint i=0; i!=c.size(); ++i) { r.adjoin(arguments(c[i])); }
    return r;
}


const RealIntervalSet RealVariableInterval::interval() const {
    return RealIntervalSet(this->_lower,this->_upper);
}

const Interval RealVariableInterval::approximate_interval() const {
    return approximation(this->interval());
}

OutputStream& operator<<(OutputStream& os, const RealVariableInterval& eivl) {
    return os << eivl.variable() << ".in(" << eivl.lower() << "," << eivl.upper() << ")";
}


RealVariablesBox::RealVariablesBox(const List<RealVariableInterval>& lst)
{
    for(uint i=0; i!=lst.size(); ++i) {
        _bounds.insert(lst[i].variable(),lst[i].interval());
    }
}

RealVariablesBox::RealVariablesBox(const RealSpace& spc, const RealBoxSet& bx)
{
    ARIADNE_ASSERT(spc.size()==bx.size());
    for(uint i=0; i!=spc.size(); ++i) {
        _bounds.insert(spc[i],bx[i]);
    }
}

RealBoxSet RealVariablesBox::euclidean_set(const RealSpace& spc) const {
    RealBoxSet bx(spc.size());
    for(uint i=0; i!=spc.size(); ++i) {
        bx[i]=(*this)[spc[i]];
    }
    return bx;
}

RealBoxSet RealVariablesBox::box(const RealSpace& spc) const {
    return this->euclidean_set(spc);
}

OutputStream& operator<<(OutputStream& os, const RealVariablesBox& ebx) {
    return os << ebx._bounds;
}



VariablesBox over_approximation(const RealVariablesBox& ebx) {
    Map<RealVariable,Interval> result;
    for(Map<RealVariable,RealIntervalSet>::const_iterator iter=ebx.bounds().begin();
        iter!=ebx.bounds().end(); ++iter)
    {
        result[iter->first]=over_approximation(iter->second);
    }
    return result;
}

VariablesBox approximation(const RealVariablesBox& ebx) {
    Map<RealVariable,Interval> result;
    for(Map<RealVariable,RealIntervalSet>::const_iterator iter=ebx.bounds().begin();
        iter!=ebx.bounds().end(); ++iter)
    {
        result[iter->first]=approximation(iter->second);
    }
    return result;
}

VariablesBox under_approximation(const RealVariablesBox& ebx) {
    Map<RealVariable,Interval> result;
    for(Map<RealVariable,RealIntervalSet>::const_iterator iter=ebx.bounds().begin();
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

RealConstraintSet RealExpressionConstraintSet::euclidean_set(const RealSpace& space) const {
    const RealExpressionConstraintSet& set = *this;
    List<RealConstraint> constraints;
    for(uint i=0; i!=set.constraints().size(); ++i) {
        RealExpression constraint_expression=indicator(set.constraints()[i],NEGATIVE);
        RealScalarFunction constraint_function( Ariadne::make_function(constraint_expression,space) );
        constraints.append( constraint_function <= Real(0) );
    }
    return RealConstraintSet(constraints);
}

OutputStream& operator<<(OutputStream& os, const RealExpressionConstraintSet& eset) {
    os << "[";
    for(List<ContinuousPredicate>::const_iterator iter=eset._constraints.begin(); iter!=eset._constraints.end(); ++iter) {
        os << (iter==eset._constraints.begin()?"":",") << *iter; }
    return os << "]";
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

RealBoundedConstraintSet RealExpressionBoundedConstraintSet::euclidean_set(const RealSpace& space) const {
    const RealExpressionBoundedConstraintSet& set = *this;
    RealBoxSet domain=RealVariablesBox(set.bounds()).euclidean_set(space);
    List<RealConstraint> constraints;
    for(uint i=0; i!=set.constraints().size(); ++i) {
        RealExpression constraint_expression=indicator(set.constraints()[i],NEGATIVE);
        RealScalarFunction constraint_function( Ariadne::make_function(constraint_expression,space) );
        constraints.append( constraint_function <= Real(0) );
    }
    return RealBoundedConstraintSet(domain,constraints);}

OutputStream& operator<<(OutputStream& os, const RealExpressionBoundedConstraintSet& eset) {
    os << "[";
    for(Map<RealVariable,RealIntervalSet>::const_iterator iter=eset._bounds.begin(); iter!=eset._bounds.end(); ++iter) {
        os << (iter==eset._bounds.begin()?"":",") << *iter; }
    os << ";";
    for(List<ContinuousPredicate>::const_iterator iter=eset._constraints.begin(); iter!=eset._constraints.end(); ++iter) {
        os << (iter==eset._constraints.begin()?"":",") << *iter; }
    return os << "]";
    return os << eset._bounds << eset._constraints;
}

IntervalConstrainedImageSet approximate_euclidean_set(const RealExpressionBoundedConstraintSet& set, const RealSpace& space) {
    IntervalVector domain=approximation(RealVariablesBox(set.bounds()).euclidean_set(space));
    IntervalVectorFunction identity=IntervalVectorFunction::identity(domain.size());

    IntervalConstrainedImageSet result(domain,identity);
    //List<IntervalConstraint> constraints;
    for(uint i=0; i!=set.constraints().size(); ++i) {
        RealExpression constraint_expression=indicator(set.constraints()[i],NEGATIVE);
        IntervalScalarFunction constraint_function( Ariadne::make_function(constraint_expression,space) );
        result.new_parameter_constraint(constraint_function <= Float(0) );
        //constraints.append( constraint_function <= 0.0 );
    }
    return result;
    //return IntervalConstrainedImageSet(domain,identity,constraints);
}





} // namespace Ariadne
