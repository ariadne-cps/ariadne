/***************************************************************************
 *            valuation.h
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


/*! \file valuation.h
 *  \brief Valuations over named variables.
 */

#ifndef ARIADNE_VALUATION_H
#define ARIADNE_VALUATION_H

#include <cstdarg>
#include <iostream>
#include <string>

#include "macros.h"
#include "container.h"
#include "tribool.h"

#include "expression.h"

namespace Ariadne {

typedef std::string String;

// Sequencing operators to make Valuation objects or objects convertible to Valuations.
template<class X> inline Map<Identifier,X> operator|(const Variable<X>& v, const X& c) { Map<Identifier,X> r; r.insert(v.name(),c); return r; }
template<class X> inline Map<Identifier,X> operator|(const Variable<X>& v, const Constant<X>& c) { Map<Identifier,X> r; r.insert(v.name(),c); return r; }
template<class K, class V> inline  Map<K,V> operator,(const Map<K,V>& m1, const Map<K,V>& m2) {
    Map<K,V> r=m1; for(typename Map<K,V>::const_iterator iter=m2.begin(); iter!=m2.end(); ++iter) { r.insert(*iter); } return r; }
Map<Identifier,String> inline operator|(const Variable<String>& v, const char* c) { return v|String(c); }

template<class V, class X=V> class Valuation;
typedef Valuation<String> StringValuation;
typedef Valuation<Integer> IntegerValuation;

template<class V, class X>
class Valuation
{
  public:
    typedef typename Map<Identifier,X>::const_iterator const_iterator;
    typedef V Type;
    typedef X ValueType;
  public:
    Valuation() { }
    Valuation(const Map<Identifier,ValueType>& m) : _values(m) { }
    void insert(const Variable<Type>& v, const ValueType& s) { this->_values.insert(v.name(),s); }
    void set(const Variable<Type>& v, const ValueType& s) { this->_values[v.name()]=s; }
    const ValueType& get(const Variable<Type>& v) const { return _values[v.name()]; }
    const ValueType& operator[](const Variable<Type>& v) const { return _values[v.name()]; }
    ValueType& operator[](const Variable<Type>& v) { return _values[v.name()]; }
    const Map<Identifier,ValueType>& values() const { return _values; }
    Set<Identifier> defined() const { return _values.keys(); }
    const_iterator begin() const { return _values.begin(); }
    const_iterator end() const { return _values.end(); }
  public:
    Map<Identifier,ValueType> _values;
};

template<class V, class X> bool operator==(const Valuation<V,X>& v1, const Valuation<V,X>& v2) {
    bool identical = true;
    const Map<Identifier,X>& v1sm=v1.values();
    const Map<Identifier,X>& v2sm=v2.values();
    typename Map<Identifier,X>::const_iterator v1iter=v1sm.begin();
    typename Map<Identifier,X>::const_iterator v2iter=v2sm.begin();

    while(v1iter!=v1sm.end() && v2iter!=v2sm.end()) {
        if(v1iter->first==v2iter->first) {
            if(v1iter->second != v2iter->second) { return false; }
            ++v1iter; ++v2iter;
        } else if(v1iter->first<v2iter->first) {
            identical=false; ++v1iter;
        } else {
            identical=false; ++v2iter;
        }
    }
    if(v1iter!=v1sm.end() || v2iter!=v2sm.end()) { identical=false; }
    if(!identical) { ARIADNE_THROW(std::runtime_error,"operator==(DiscreteLocation,DiscreteLocation)",
                                   "Valuations "<<v1<<" and "<<v2<<" are not identical, but no values differ."); }
    return true;
}



//! \brief
class DiscreteValuation
    : public Valuation<String>, public Valuation<Integer>
{
    typedef String StringType;
    typedef Integer IntegerType;
  public:
    DiscreteValuation() { }
    DiscreteValuation(const Map<Identifier,StringType>& sm) : Valuation<String>(sm) { }
    DiscreteValuation(const Map<Identifier,IntegerType>& zm) : Valuation<Integer>(zm) { }
    DiscreteValuation(const Map<Identifier,StringType>& sm,const Map<Identifier,IntegerType>& zm) : Valuation<String>(sm), Valuation<Integer>(zm) { }
    using Valuation<String>::insert; using Valuation<Integer>::insert;
    using Valuation<String>::get; using Valuation<Integer>::get;
    using Valuation<String>::set; using Valuation<Integer>::set;
    using Valuation<String>::operator[]; using Valuation<Integer>::operator[];
    const Map<Identifier,StringType>& string_values() const { return Valuation<String>::values(); }
    const Map<Identifier,IntegerType>& integer_values() const { return Valuation<Integer>::values(); }
};

inline bool operator==(const DiscreteValuation& v1, const DiscreteValuation& v2) {
    return static_cast<const StringValuation&>(v1)==static_cast<const StringValuation&>(v2)
        && static_cast<const IntegerValuation&>(v1)==static_cast<const IntegerValuation&>(v2); }

template<class X>
class ContinuousValuation
    : public Valuation<Real,X>
{
  public:
    typedef X RealType;
};


template<class X> class HybridValuation
    : public DiscreteValuation
    , public ContinuousValuation<X>
{
  public:
    using DiscreteValuation::set;
    using ContinuousValuation<X>::set;
    using DiscreteValuation::get;
    using ContinuousValuation<X>::get;
    using DiscreteValuation::operator[];
    using ContinuousValuation<X>::operator[];
};

inline std::ostream& operator<<(std::ostream& os, const DiscreteValuation& val) {
    const char open='('; const char mid='|'; char const close=')';
    //const char open='{'; const char mid=':'; char const close='}';
    char sep=open;
    for(Map<Identifier,DiscreteValuation::StringType>::const_iterator iter=val.string_values().begin(); iter!=val.string_values().end(); ++iter) {
        os << sep << iter->first << mid << iter->second; sep=','; }
    for(Map<Identifier,DiscreteValuation::IntegerType>::const_iterator iter=val.integer_values().begin(); iter!=val.integer_values().end(); ++iter) {
        os << sep << iter->first << mid << iter->second; }
    return os << close;
}

template<class X> inline std::ostream& operator<<(std::ostream& os, const ContinuousValuation<X>& val) {
    return os << val.real_values();
}

template<class X> inline std::ostream& operator<<(std::ostream& os, const Valuation<X>& val) {
    const char open='('; const char mid ='|'; const char close=')'; const char sep=',';
    //const char open='{'; const char mid=':'; char const close='}'; char const sep=",";
    os << open;
    for(typename Map<Identifier,X>::const_iterator iter=val.values().begin(); iter!=val.values().end(); ++iter) {
        if(iter!=val.values().begin()) { os << sep; }
        os << iter->first << mid << iter->second;
    }
    return os << close;
}


Boolean evaluate(const Expression<Boolean>&, const StringValuation&);
String evaluate(const Expression<String>&, const StringValuation&);
Integer evaluate(const Expression<Integer>&, const IntegerValuation&);
Boolean evaluate(const Expression<Boolean>&, const DiscreteValuation&);
template<class X> X evaluate(const Expression<Real>& e, const ContinuousValuation<X>&);
template<class X> Tribool evaluate(const Expression<Tribool>&, const ContinuousValuation<X>&);



} // namespace Ariadne

#endif // ARIADNE_VALUATION_H
