/***************************************************************************
 *            discrete_location.h
 *
 *  Copyright  2004-9  Alberto Casagrande, Pieter Collins
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

/*! \file discrete_location.h
 *  \brief Class representing a discrete location of a hybrid automaton.
 */

#ifndef ARIADNE_DISCRETE_LOCATION_H
#define ARIADNE_DISCRETE_LOCATION_H

#include "container.h"
#include "variables.h"

namespace Ariadne {

enum Urgency { urgent, permissive };

typedef Variable<String> StringVariable;

class IncompleteLocationError : public std::runtime_error {
  public:
      IncompleteLocationError(const std::string& what) : std::runtime_error(what) { }
};

//! \brief Type of a  discrete location of a hybrid system.
class DiscreteLocation
{
  public:
    DiscreteLocation() { }
    explicit DiscreteLocation(const int& n) {
        this->_valuation.insert(StringVariable("q"),to_string(n)); }
    explicit DiscreteLocation(const std::string& s) {
        this->_valuation.insert(StringVariable("q"),s); }
    void insert(const StringVariable& v, const String& s) {
        if(_valuation.has_key(v)) { ARIADNE_THROW(std::runtime_error,"DiscreteLocation::insert","Location "<<*this<<" already has variable "<<v<<"\n"); }
        this->_valuation.insert(v,s); }
    void adjoin(const DiscreteLocation& q) {
        for(Map<StringVariable,String>::const_iterator iter=q._valuation.begin(); iter!=q._valuation.end(); ++iter) {
            this->insert(iter->first,iter->second); } }
    bool has_key(const StringVariable& v) { return this->_valuation.has_key(v); }
    const String& operator[](const StringVariable& v) const { return this->_valuation[v]; }
    friend bool operator==(const DiscreteLocation& q1, const DiscreteLocation& q2);
    friend bool operator<(const DiscreteLocation& q1, const DiscreteLocation& q2);
    friend DiscreteLocation operator,(const DiscreteLocation& os, const DiscreteLocation& q);
    friend std::ostream& operator<<(std::ostream& os, const DiscreteLocation& q);
  private:
    Map<StringVariable,String> _valuation;
};

inline bool operator==(const DiscreteLocation& q1, const DiscreteLocation& q2) {
    tribool result=true;
    Map<StringVariable,String>::const_iterator q1iter=q1._valuation.begin();
    Map<StringVariable,String>::const_iterator q2iter=q2._valuation.begin();
    while(q1iter!=q1._valuation.end() && q2iter!=q2._valuation.end()) {
        if(q1iter==q1._valuation.end() || q2iter==q2._valuation.end()) {
            result=indeterminate;
            break;
        }
        if(q1iter->first==q2iter->first) {
            if(q1iter->second != q2iter->second) {
                return false;
            }
            ++q1iter;
            ++q2iter;
        } else {
            result=indeterminate;
            if(q1iter->first < q2iter->first) {
                ++q1iter;
            } else {
                ++q2iter;
            }
        }
    }
    if(definitely(result)) {
        return true;
    }
    ARIADNE_THROW(IncompleteLocationError,"DiscreteLocation::operator==","q1="<<q1<<" and q2="<<q2<<" cannot be compared.");
    return result;
}

inline bool operator!=(const DiscreteLocation& q1, const DiscreteLocation& q2) {
    return !(q1==q2);
}

inline bool operator<(const DiscreteLocation& q1, const DiscreteLocation& q2) {
    return q1._valuation < q2._valuation;
}

inline DiscreteLocation operator%=(const StringVariable& v, const String& s) {
    DiscreteLocation loc; loc.insert(v,s); return loc;
}

inline DiscreteLocation operator,(const DiscreteLocation& q1, const DiscreteLocation& q2) {
    DiscreteLocation r(q1);
    for(Map<StringVariable,String>::const_iterator iter=q2._valuation.begin(); iter!=q2._valuation.end(); ++iter) {
        r.insert(iter->first,iter->second);
    }
    return r;
}

inline std::ostream& operator<<(std::ostream& os, const DiscreteLocation& q) {
    return os << q._valuation;
}

template<class A> inline void serialize(A& archive, DiscreteLocation& location, const uint version) {
    archive & location._valuation;
}


} //namespace Ariadne

#endif /* ARIADNE_DISCRETE_LOCATION_H */
