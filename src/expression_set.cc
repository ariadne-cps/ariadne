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
    typedef RealVariable KeyType; typedef IntervalSet ValueType;
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
    for(uint i=0; i!=c.size(); ++i) { r.adjoin(arguments(c)); }
    return r;
}


const IntervalSet RealVariableInterval::interval() const {
    return IntervalSet(this->_lower,this->_upper);
}

const Interval RealVariableInterval::approximate_interval() const {
    return approximation(this->interval());
}

OutputStream& operator<<(OutputStream& os, const RealVariableInterval& eivl) {
    return os << eivl.variable() << ".in(" << eivl.lower() << "," << eivl.upper() << ")";
}


RealVariableBox::RealVariableBox(const List<RealVariableInterval>& lst)
{
    for(uint i=0; i!=lst.size(); ++i) {
        _bounds.insert(lst[i].variable(),lst[i].interval());
    }
}

RealVariableBox::RealVariableBox(const RealSpace& spc, const BoxSet& bx)
{
    ARIADNE_ASSERT(spc.size()==bx.size());
    for(uint i=0; i!=spc.size(); ++i) {
        _bounds.insert(spc[i],bx[i]);
    }
}

BoxSet RealVariableBox::box(const List<RealVariable>& spc) const {
    const RealVariableBox& ebx = *this;
    BoxSet bx(spc.size());
    for(uint i=0; i!=spc.size(); ++i) {
        bx[i]=ebx[spc[i]];
    }
    return bx;
}

BoxSet euclidean_set(const RealVariableBox& ebx, const List<RealVariable>& spc) {
    BoxSet bx(spc.size());
    for(uint i=0; i!=spc.size(); ++i) {
        bx[i]=ebx[spc[i]];
    }
    return bx;
}

Box approximate_euclidean_set(const RealVariableBox& ebx, const List<RealVariable>& spc) {
    Box bx(spc.size());
    for(uint i=0; i!=spc.size(); ++i) {
        bx[i]=approximation(ebx[spc[i]]);
    }
    return bx;
}

OutputStream& operator<<(OutputStream& os, const RealVariableBox& ebx) {
    return os << ebx._bounds;
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
    ARIADNE_ASSERT( subset(arguments(constraints),variables(bounds)) );
}

OutputStream& operator<<(OutputStream& os, const RealExpressionBoundedConstraintSet& eset) {
    os << "[";
    for(Map<RealVariable,IntervalSet>::const_iterator iter=eset._bounds.begin(); iter!=eset._bounds.end(); ++iter) {
        os << (iter==eset._bounds.begin()?"":",") << *iter; }
    os << ";";
    for(List<ContinuousPredicate>::const_iterator iter=eset._constraints.begin(); iter!=eset._constraints.end(); ++iter) {
        os << (iter==eset._constraints.begin()?"":",") << *iter; }
    return os << "]";
    return os << eset._bounds << eset._constraints;
}

BoundedConstraintSet euclidean_set(const RealExpressionBoundedConstraintSet& set, const RealSpace& space) {
    BoxSet domain=euclidean_set(RealVariableBox(set.bounds()),space.variables());
    List<RealConstraint> constraints;
    for(uint i=0; i!=set.constraints().size(); ++i) {
        RealExpression constraint_expression=indicator(set.constraints()[i],NEGATIVE);
        RealScalarFunction constraint_function( Ariadne::make_function(constraint_expression,space) );
        constraints.append( constraint_function <= Real(0) );
    }
    return BoundedConstraintSet(domain,constraints);
}

ValidatedConstrainedImageSet approximate_euclidean_set(const RealExpressionBoundedConstraintSet& set, const RealSpace& space) {
    IntervalVector domain=approximation(euclidean_set(RealVariableBox(set.bounds()),space.variables()));
    IntervalVectorFunction identity=IntervalVectorFunction::identity(domain.size());

    ValidatedConstrainedImageSet result(domain,identity);
    //List<IntervalConstraint> constraints;
    for(uint i=0; i!=set.constraints().size(); ++i) {
        RealExpression constraint_expression=indicator(set.constraints()[i],NEGATIVE);
        IntervalScalarFunction constraint_function( Ariadne::make_function(constraint_expression,space) );
        result.new_parameter_constraint(constraint_function <= Float(0) );
        //constraints.append( constraint_function <= 0.0 );
    }
    return result;
    //return ValidatedConstrainedImageSet(domain,identity,constraints);
}





} // namespace Ariadne
