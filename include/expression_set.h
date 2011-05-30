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

#include <iosfwd>

#include <boost/shared_ptr.hpp>

#include "expression.h"
#include "space.h"
#include "function_set.h"
#include "formula.h"
#include <boost/iterator/iterator_concepts.hpp>

namespace Ariadne {

class Box;

template<class T> struct KeyValue { };

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

template<class K, class V> List<V> order(const Map<K,V>& map, const List<K>& ordering) {
    List<V> result;
    result.reserve(ordering.size());
    for(uint i=0; i!=ordering.size(); ++i) {
        result.append(map[ordering[i]]);
    }
    return result;
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

template<class X> class ExpressionInterval;

template<class X> struct GeometricInterval {
    X _lower;
    X _upper;
    GeometricInterval() : _lower(-infty), _upper(+infty) { }
    GeometricInterval(const X& l, const X& u) : _lower(l), _upper(u) { }
    operator Interval() const { return Interval(Interval(this->_lower).lower(),Interval(this->_upper).upper()); }
};
template<class X> inline std::ostream& operator<<(std::ostream& os, const GeometricInterval<X>& ivl) {
    return os << "[" << ivl._lower << "," << ivl._upper << "]";
}
typedef GeometricInterval<Real> RealInterval;

template<class X> struct ExpressionLowerHalfspace {
    X _lower;
    Variable<X> _variable;
    ExpressionLowerHalfspace(const X& l, const Variable<X>& v) : _lower(l), _variable(v) { }
    operator Expression<tribool>() const { return (this->_lower <= Expression<X>(this->_variable)); };
};

template<class X> struct ExpressionUpperHalfspace {
    Variable<X> _variable;
    X _upper;
    ExpressionUpperHalfspace(const Variable<X>& v, const X& u) : _variable(v), _upper(u)  { }
    operator Expression<tribool>() const { return (Expression<X>(this->_variable) <= this->_upper); }
};

template<class X> class ExpressionInterval {
  public:
    X _lower;
    Variable<X> _variable;
    X _upper;
  public:
    ExpressionInterval<X>(const X& l, const Variable<X>& v, const X& u)
        : _lower(l), _variable(v), _upper(u) { ARIADNE_ASSERT_MSG(l<=u,"Interval("<<l<<","<<u<<") not provably nonempty"); }
    ExpressionInterval<X>(const ExpressionLowerHalfspace<X>& lv)
        : _lower(lv._lower), _variable(lv._variable), _upper(infty) { }
    ExpressionInterval<X>(const ExpressionUpperHalfspace<X>& vu)
        : _lower(-infty), _variable(vu._variable), _upper(vu._upper) { }
    Variable<X> const& variable() const { return this->_variable; }
};
template<class X> inline std::ostream& operator<<(std::ostream& os, const ExpressionInterval<X>& eivl) {
    //return os << eivl._lower << "<=" << eivl._variable << "<=" << eivl._upper;
    return os << eivl._variable << ".in(" << eivl._lower << "," << eivl._upper << ")";
}


template<class X> inline ExpressionInterval<X> Variable<X>::in(const X& l, const X& u) {
    return ExpressionInterval<X>(l,*this,u); }

template<class X> struct KeyValue< ExpressionInterval<X> > {
    typedef Variable<X> KeyType;
    typedef GeometricInterval<X> ValueType;
    const Variable<X>& key(const ExpressionInterval<X>& eivl) { return eivl.variable(); }
    RealInterval value(const ExpressionInterval<X>& eivl) const { return RealInterval(eivl._lower,eivl._upper); }
};

typedef ExpressionLowerHalfspace<Real> RealExpressionLowerHalfspace;
typedef ExpressionUpperHalfspace<Real> RealExpressionUpperHalfspace;
typedef ExpressionInterval<Real> RealExpressionInterval;

inline RealExpressionInterval operator<=(const RealExpressionLowerHalfspace& lv, const Real& u) {
    return RealExpressionInterval(lv._lower,lv._variable,u);
}

inline RealExpressionLowerHalfspace operator<=(const Real& l, const RealVariable& v) {
    return RealExpressionLowerHalfspace(l,v);
}

inline RealExpressionLowerHalfspace operator>=(const RealVariable& v, const Real& l) {
    return RealExpressionLowerHalfspace(l,v);
}

inline RealExpressionUpperHalfspace operator<=(const RealVariable& v, const Real& u) {
    return RealExpressionUpperHalfspace(v,u);
}

inline RealExpressionUpperHalfspace operator>=(const Real& u, const RealVariable& v) {
    return RealExpressionUpperHalfspace(v,u);
}

class ExpressionBox {
    Map<RealVariable,RealInterval> _intervals;
  public:
    ExpressionBox(const List<RealExpressionInterval>& lst) {
        for(uint i=0; i!=lst.size(); ++i) { _intervals[lst[i]._variable]=RealInterval(lst[i]._lower,lst[i]._upper); } 
    }
    RealSpace space() const { return RealSpace(make_list(_intervals.keys())); }
    Box euclidean_box() const;
    const RealInterval& operator[](const RealVariable& v) const { return this->_intervals[v]; }
    friend std::ostream& operator<<(std::ostream& os, const ExpressionBox& ebx) {
        return os << ebx._intervals; }
};

inline Box euclidean_box(const ExpressionBox& ebx, const List<RealVariable>& spc) {
    Box bx(spc.size());
    for(uint i=0; i!=spc.size(); ++i) {
        bx[i]=Interval(ebx[spc[i]]);
    }
    return bx;
}

inline Box ExpressionBox::euclidean_box() const {
    return Ariadne::euclidean_box(*this,this->space().variables());
}

//! \ingroup GeometryModule ExactSetSubModule
//! \brief A set defined as the image of a box under a continuous function.
//! The set is described as \f$S=h(D) = \{ h(s) \mid s \in D\}\f$ where \f$D\f$ is the domain and \f$h\f$ the function.
class ExpressionSet
{
    List<RealExpressionInterval> _domain;
    List<ContinuousPredicate> _constraints;
  public:
    ExpressionSet(const List<RealExpressionInterval>& domain, const List<ContinuousPredicate>& constraints)
        : _domain(domain), _constraints(constraints)
    {
        ARIADNE_ASSERT(unique_keys(domain));
    }
    ExpressionSet(const List<RealExpressionInterval>& domain)
        : _domain(domain), _constraints()
    {
        ARIADNE_ASSERT(unique_keys(domain));
    }
    Map<RealVariable,RealInterval> domain() const { return make_key_value_map(this->_domain); }
    List<ContinuousPredicate> const& constraints() const { return this->_constraints; }
    friend std::ostream& operator<<(std::ostream& os, const ExpressionSet& eset) {
        os << "[";
        for(uint i=0; i!=eset._domain.size(); ++i) {
            os << (i==0?"":",") << eset._domain[i]; }
        os << ";";
        for(uint i=0; i!=eset._constraints.size(); ++i) {
            os << (i==0?"":",") << eset._constraints[i]; }
        return os << "]";
        return os << eset._domain << eset._constraints; }
};

inline BoundedConstraintSet euclidean_set(const ExpressionSet& set, const RealSpace& space) {
    IntervalVector domain(order(set.domain(),space.variables()));
    List<RealNonlinearConstraint> constraints;
    for(uint i=0; i!=set.constraints().size(); ++i) {
        RealExpression constraint_expression=indicator(set.constraints()[i],NEGATIVE);
        RealScalarFunction constraint_function( Ariadne::dimension(space),Ariadne::formula(constraint_expression,space) );
        constraints.append( constraint_function <= Real(0) );
    }
    return BoundedConstraintSet(domain,constraints);
}



} //namespace Ariadne



#endif
