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

#include <string>
#include "utility/container.h"
#include "expression/valuation.h"
#include <boost/iterator/iterator_concepts.hpp>

namespace Ariadne {

//! \ingroup SystemModule
//! \brief Type of a  discrete location of a hybrid system.
class DiscreteLocation
    : public StringValuation
{
  public:
    //! \brief Construct a location which does not set any discrete variables.
    DiscreteLocation() : StringValuation() { }
    explicit DiscreteLocation(const std::string& str) : StringValuation() { this->insert(StringVariable("q"),str); }
    //! \brief Construct a location for which the string variable named \a var is given value \a val.
    explicit DiscreteLocation(const Identifier& var, const String& val) : StringValuation() { this->insert(StringVariable(var),val); }
    explicit DiscreteLocation(const std::string& var, const std::string& val) : StringValuation() { this->insert(StringVariable(var),val); }
    explicit DiscreteLocation(const int& num) : StringValuation() { this->insert(StringVariable("q"),to_str(num)); }
    DiscreteLocation(const DiscreteLocation&) = default;
    DiscreteLocation(const Map<Identifier,String>& sm) : StringValuation(sm) { }
    void adjoin(const DiscreteLocation& loc) { this->_values.adjoin(loc._values); }
};

//! \relates DiscreteLocation \brief Combine the values of variables of two locations.
DiscreteLocation operator,(const DiscreteLocation& loc1, const DiscreteLocation& loc2);
bool operator==(const DiscreteLocation& loc1, const DiscreteLocation& loc2);
bool operator!=(const DiscreteLocation& loc1, const DiscreteLocation& loc2);
bool operator<(const DiscreteLocation& loc1, const DiscreteLocation& loc2);

//! \relates DiscreteLocation \brief Test if two locations are distinguishable i.e. specify a common variable which takes different values.
inline bool are_distinguishable(const DiscreteLocation& location1, const DiscreteLocation& location2) {
    for(Map<Identifier,String>::const_iterator iter1=location1._values.begin(); iter1!=location1._values.end(); ++iter1) {
        if(location2._values.has_key(iter1->first) && location2._values[iter1->first] != iter1->second) {
            return true;
        }
    }
    return false;
}

//! \relates DiscreteLocation \brief Test if \a partial_location is defined by a restricted set of variables of \a full_location.
inline bool is_restriction(const DiscreteLocation& partial_location, const DiscreteLocation& full_location) {
    for(Map<Identifier,String>::const_iterator value_iter=partial_location._values.begin(); value_iter!=partial_location._values.end(); ++value_iter) {
        if(!full_location._values.has_key(value_iter->first) || full_location[StringVariable(value_iter->first)]!=value_iter->second) {
            return false;
        }
    }
    return true;
}

inline DiscreteLocation restrict(const DiscreteLocation& location, const Set<Identifier>& variables) {
    ARIADNE_ASSERT_MSG(subset(location._values.keys(),variables)," location "<<location<<" variables not a subset of "<<variables);
    return restrict_keys(location._values,variables);
}

inline DiscreteLocation restrict(const DiscreteLocation& location, const List<Identifier>& variables) {
    return restrict(location._values,Set<Identifier>(variables));
}

} //namespace Ariadne

#endif /* ARIADNE_DISCRETE_LOCATION_H */
