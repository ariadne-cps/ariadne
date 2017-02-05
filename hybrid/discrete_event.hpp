/***************************************************************************
 *            discrete_event.h
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

/*! \file discrete_event.h
 *  \brief Class representing a discrete event.
 */

#ifndef ARIADNE_DISCRETE_EVENT_H
#define ARIADNE_DISCRETE_EVENT_H

namespace Ariadne {

//! \ingroup SystemModule
//! \brief Type of a  discrete event of a hybrid system.
class DiscreteEvent {
  public:
    DiscreteEvent() : _id("e?") { }
    DiscreteEvent(Int n) : _id(StringType("e"+to_str(n))) { if(n<0) { _id=StringType("i"+to_str(-n)); } }
    DiscreteEvent(const StringType& s) : _id(s) { }
    StringType name() const { return this->_id; }
    Bool operator==(const DiscreteEvent& e) const { return this->_id==e._id; }
    Bool operator!=(const DiscreteEvent& e) const { return this->_id!=e._id; }
    Bool operator<=(const DiscreteEvent& e) const { return this->_id<=e._id; }
    Bool operator>=(const DiscreteEvent& e) const { return this->_id>=e._id; }
    Bool operator< (const DiscreteEvent& e) const { return this->_id< e._id; }
    Bool operator> (const DiscreteEvent& e) const { return this->_id> e._id; }
    friend OutputStream& operator<<(OutputStream& os, const DiscreteEvent& e) {
        return os << e._id; }
  private:
    StringType _id;
};

#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> inline Void serialize(A& archive, DiscreteEvent& event, const Nat version) {
      StringType& id=reinterpret_cast<StringType&>(event);
      archive & id;
  }
#endif /* ARIADNE_ENABLE_SERIALIZATION */

} //namespace Ariadne


#endif /* ARIADNE_DISCRETE_EVENT_H */
