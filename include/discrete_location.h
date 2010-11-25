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

namespace Ariadne {

template<class T> class List;

//! \brief Type of a  discrete location of an atomic hybrid automaton.
class AtomicDiscreteLocation {
  public:
    AtomicDiscreteLocation() : _id("q?") { }
    AtomicDiscreteLocation(int n) : _id(std::string("q"+to_str(n))) { }
    AtomicDiscreteLocation(const std::string& s) : _id(s) { }
    std::string name() const { return this->_id; }
    bool operator==(const AtomicDiscreteLocation& q) const { return this->_id==q._id; }
    bool operator!=(const AtomicDiscreteLocation& q) const { return this->_id!=q._id; }
    bool operator<=(const AtomicDiscreteLocation& q) const { return this->_id<=q._id; }
    bool operator>=(const AtomicDiscreteLocation& q) const { return this->_id>=q._id; }
    bool operator< (const AtomicDiscreteLocation& q) const { return this->_id< q._id; }
    bool operator> (const AtomicDiscreteLocation& q) const { return this->_id> q._id; }
    friend std::ostream& operator<<(std::ostream& os, const AtomicDiscreteLocation& q);
  private:
    std::string _id;
};

inline std::ostream& operator<<(std::ostream& os, const AtomicDiscreteLocation& q) {
    return os << q._id; }

//! \brief Type of a  discrete location of a hybrid system.
//! Currently implemented as a list of String elements, but this may change...
class DiscreteLocation
{
    typedef size_t size_type;
  public:
    DiscreteLocation() : _lst() { }
    explicit DiscreteLocation(int n): _lst(1u,AtomicDiscreteLocation(n)) { }
    explicit DiscreteLocation(const std::string& s) : _lst(1u,AtomicDiscreteLocation(s)) { }
    DiscreteLocation(const AtomicDiscreteLocation& q) : _lst(1u,q) { }
    DiscreteLocation(const List<AtomicDiscreteLocation>& l) : _lst(l) { }
    size_type size() const { return _lst.size(); }
    const AtomicDiscreteLocation& operator[](size_type i) const { assert(i<_lst.size()); return _lst[i]; }
    void append(AtomicDiscreteLocation q) { _lst.append(q); }
    friend std::ostream& operator<<(std::ostream& os, const DiscreteLocation& q) { return os << q._lst; }
  private:
    List<AtomicDiscreteLocation> _lst;
};

inline bool operator==(const DiscreteLocation& q1, const DiscreteLocation& q2) {
    if(q1.size()!=q2.size()) { return false; }
    for(uint i=0; i!=q1.size(); ++i) { if(q1[i]!=q2[i]) { return false; } }
    return true;
}

inline bool operator!=(const DiscreteLocation& q1, const DiscreteLocation& q2) {
    return !(q1==q2);
}

inline bool operator<(const DiscreteLocation& q1, const DiscreteLocation& q2) {
    if(q1.size()!=q2.size()) { return q1.size() < q2.size(); }
    for(uint i=0; i!=q1.size(); ++i) { if(q1[i]!=q2[i]) { return q1[i]<q2[i]; } }
    return false;
}

inline DiscreteLocation operator,(const AtomicDiscreteLocation& q1, const AtomicDiscreteLocation& q2) {
    DiscreteLocation loc; loc.append(q1); loc.append(q2); return loc; }


template<class A> inline void serialize(A& archive, AtomicDiscreteLocation& state, const uint version) {
    std::string& id=reinterpret_cast<std::string&>(state);
    archive & id;
}

template<class A> inline void serialize(A& archive, DiscreteLocation& location, const uint version) {
    std::vector<AtomicDiscreteLocation>& vec=static_cast<std::vector<AtomicDiscreteLocation>&>(location);
    archive & vec;
}

} //namespace Ariadne

#endif /* ARIADNE_DISCRETE_LOCATION_H */
