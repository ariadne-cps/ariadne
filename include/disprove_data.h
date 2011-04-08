/***************************************************************************
 *            disprove_data.h
 *
 *  Copyright 2010  Luca Geretti
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

/*! \file disprove_data.h
 *  \brief Data structure for holding and manipulating disprovement results.
 */

#ifndef DISPROVE_DATA_H_
#define DISPROVE_DATA_H_

#include "discrete_state.h"
#include "box.h"
#include "hybrid_set.h"
#include "vector.h"

namespace Ariadne {

/**
 * \brief This is a structure for keeping data obtained during lower chain reach, in order to get info for falsification.
 */
struct DisproveData
{
private:

	bool _isDisproved;
	HybridBoxes _reachBounds;
	HybridFloatVector _epsilon;

public:

	DisproveData(const HybridSpace& space): _isDisproved(false) {
		for (HybridSpace::const_iterator space_it = space.begin(); space_it != space.end(); ++space_it) {
			_reachBounds.insert(std::pair<DiscreteState,Box>(space_it->first,Box::empty_box(space_it->second)));
			_epsilon.insert(std::pair<DiscreteState,Vector<Float> >(space_it->first,Vector<Float>(space_it->second)));
		}
	}

	DisproveData(bool isDisproved, HybridBoxes reachBounds, HybridFloatVector epsilon): _isDisproved(isDisproved),
																							 _reachBounds(reachBounds),
																							 _epsilon(epsilon) { }

	const bool& getIsDisproved() const { return _isDisproved; }
	const HybridBoxes& getReachBounds() const { return _reachBounds; }
	const HybridFloatVector& getEpsilon() const { return _epsilon; }

	/** Updates the content based on another info object */
	void updateWith(const DisproveData& otherInfo) {
		_isDisproved = _isDisproved || otherInfo.getIsDisproved();
		_reachBounds = hull(_reachBounds,otherInfo.getReachBounds());
		_epsilon = max(_epsilon,otherInfo.getEpsilon());
	}

	void updateIsDisproved(const bool& isDisproved) {
		_isDisproved = _isDisproved || isDisproved;
	}

	void updateReachBounds(const DiscreteState& location, const Box& bounds) {
		ARIADNE_ASSERT_MSG(_reachBounds.find(location) != _reachBounds.end(),"The reach bounds location provided is not present into the target falsification info.");
		_reachBounds[location] = hull(_reachBounds[location],bounds);
	}

	void updateEpsilon(const DiscreteState& location, const Vector<Float>& epsilon) {
		ARIADNE_ASSERT_MSG(_epsilon.find(location) != _epsilon.end(),"The epsilon location provided is not present into the target falsification info.");
		_epsilon[location] = max(_epsilon[location],epsilon);
	}

	virtual std::ostream& write(std::ostream& os) const {
		os << "Disproved: " << _isDisproved << "; Reach bounds: " << _reachBounds << "; Epsilon: " << _epsilon;
		return os;
	}
};

inline std::ostream& operator<<(std::ostream& os, const DisproveData& falsInfo) {
    return falsInfo.write(os); }

}

#endif /* DISPROVE_DATA_H_ */
