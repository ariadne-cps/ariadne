/***************************************************************************
 *            discrete_state.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
#ifndef ARIADNE_DISCRETE_STATE_H
#define ARIADNE_DISCRETE_STATE_H

namespace Ariadne {
  namespace Geometry {

    class DiscreteState;
    std::ostream& operator<<(std::ostream& os, const DiscreteState& ds);

    /*! \brief A class representing a discrete state of a hybrid system. */
    class DiscreteState {
     public:
      //! \brief Default constructor. Used for serialization. 
      DiscreteState() : _id() { }
      //! \brief Construct from an identifier. 
      explicit DiscreteState(const id_type& id) : _id(id) { }
      const id_type& id() const { return _id; }
      bool operator==(const DiscreteState& other) const { return this->_id==other._id; }
      bool operator!=(const DiscreteState& other) const { return this->_id!=other._id; }
      bool operator<(const DiscreteState& other) const { return this->_id<other._id; }
     private:
      friend std::ostream& operator<<(std::ostream& os, const DiscreteState& ds);
     private:
      id_type _id;
    };


    inline 
    std::ostream& 
    operator<<(std::ostream& os, const DiscreteState& state) { 
      return os << state._id; 
    }

    template<class A> void serialize(A& a, Geometry::DiscreteState& ds, const uint v) {
      id_type& id = const_cast<id_type&>(ds.id()); a & id;
    }

  }
}

#endif
