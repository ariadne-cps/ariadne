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
#include "expression.h"

namespace Ariadne {




//! \brief
class DiscreteValuation {
    typedef std::string StringType;
    typedef int IntegerType;
  public:
    void set(const Variable<String>& v, const StringType& s) { this->_string_values[v.name()]=s; }
    void set(const Variable<Integer>& v, const IntegerType& n) { this->_integer_values[v.name()]=n; }
    void set(const Variable<EnumeratedValue>& v, const EnumeratedValue& e) { ARIADNE_NOT_IMPLEMENTED; }
    const StringType& get(const Variable<String>& v) const { return _string_values[v.name()]; }
    const IntegerType& get(const Variable<Integer>& v) const { return _integer_values[v.name()]; }
    const Map<Identifier,StringType>& string_values() const { return _string_values; }
    const Map<Identifier,IntegerType>& integer_values() const { return _integer_values; }
  public:
    Map<Identifier,StringType> _string_values;
    Map<Identifier,IntegerType> _integer_values;
};


template<class X> class ContinuousValuation {
  public:
    typedef X RealType;
    void set(const Variable<Real>& v, const RealType& x) { this->_real_values[v.name()]=x; }
    const RealType& get(const Variable<Real>& v) const { return _real_values[v.name()]; }
    Map<Identifier,RealType>& real_values() { return _real_values; }
    const Map<Identifier,RealType>& real_values() const { return _real_values; }
    RealType& operator[](const Variable<Real>& v) { return _real_values[v.name()]; }
    const RealType& operator[](const Variable<Real>& v) const { return _real_values[v.name()]; }
  private:
    Map<Identifier,RealType> _real_values;
};


template<class X> class Valuation
    : public DiscreteValuation
    , public ContinuousValuation<X>
{
  public:
    using DiscreteValuation::set;
    using ContinuousValuation<X>::set;
    using DiscreteValuation::get;
    using ContinuousValuation<X>::get;
    using ContinuousValuation<X>::operator[];
};

inline std::ostream& operator<<(std::ostream& os, const DiscreteValuation& val) {
    char sep='{';
    for(Map<Identifier,DiscreteValuation::StringType>::const_iterator iter=val.string_values().begin(); iter!=val.string_values().end(); ++iter) {
        os << sep << iter->first << ":" << iter->second; sep=','; }
    for(Map<Identifier,DiscreteValuation::IntegerType>::const_iterator iter=val.integer_values().begin(); iter!=val.integer_values().end(); ++iter) {
        os << sep << iter->first << ":" << iter->second; }
    return os << '}';
}

template<class X> inline std::ostream& operator<<(std::ostream& os, const ContinuousValuation<X>& val) {
    return os << val.real_values();
}

template<class X> inline std::ostream& operator<<(std::ostream& os, const Valuation<X>& val) {
    char sep='{';
    for(Map<Identifier,DiscreteValuation::StringType>::const_iterator iter=val._string_values.begin(); iter!=val._string_values.end(); ++iter) {
        os << sep << iter->first << ":" << iter->second; sep=','; }
    for(Map<Identifier,DiscreteValuation::IntegerType>::const_iterator iter=val._integer_values.begin(); iter!=val._integer_values.end(); ++iter) {
        os << sep << iter->first << ":" << iter->second; }
    for(typename Map<Identifier,X>::const_iterator iter=val._real_values.begin(); iter!=val._real_values.end(); ++iter) {
        os << sep << iter->first << ":" << iter->second; }
    return os << '}';
}


} // namespace Ariadne

#endif // ARIADNE_VALUATION_H
