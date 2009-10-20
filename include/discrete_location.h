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

template<class T> std::string to_str(const T& t) {
    std::stringstream ss; ss<<t; return ss.str(); }
template<class T> std::string to_string(const T& t) {
    std::stringstream ss; ss<<t; return ss.str(); }

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

class DiscreteLocation
    : public List<AtomicDiscreteLocation>
{
  public:
    DiscreteLocation() : List<AtomicDiscreteLocation>() { }
    DiscreteLocation(const AtomicDiscreteLocation& q) : List<AtomicDiscreteLocation>(1u,q) { }
    DiscreteLocation(const List<AtomicDiscreteLocation>& l) : List<AtomicDiscreteLocation>(l) { }
};

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
