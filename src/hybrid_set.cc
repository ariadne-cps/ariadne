/***************************************************************************
 *            hybrid_set.cc
 *
 *  Copyright 2008  Pieter Collins
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
 
#include "hybrid_set.h"

namespace Ariadne {

HybridBoxes
hull(const HybridBoxes& box1, const HybridBoxes& box2)
{
	ARIADNE_ASSERT_MSG(box1.size() == box2.size(),"The two hybrid boxes must have the same number of locations ("
												  << box1.size() << " vs " << box2.size() << ").");
	HybridBoxes result;
	for (HybridBoxes::const_iterator box1_it = box1.begin(); box1_it != box1.end(); ++box1_it) {
		HybridBoxes::const_iterator box2_it = box2.find(box1_it->first);
		ARIADNE_ASSERT_MSG(box2_it != box2.end(),"The location " << box1_it->first.name() << " is not present in both hybrid boxes.");
		result.insert(std::pair<DiscreteState,Box>(box1_it->first,hull(box1_it->second,box2_it->second)));
	}

	return result;
}

HybridBoxes
shrink_in(const HybridBoxes& box, const HybridFloatVector& epsilon)
{
	HybridBoxes result;

	for (HybridBoxes::const_iterator loc_it = box.begin(); loc_it != box.end(); ++loc_it) {
		HybridFloatVector::const_iterator epsilon_it = epsilon.find(loc_it->first);
		ARIADNE_ASSERT_MSG(epsilon_it != epsilon.end(),"The location " << loc_it->first.name() << " is not present in the epsilon map.");
		result.insert(std::pair<DiscreteState,Box>(loc_it->first,loc_it->second.shrink_in(epsilon_it->second)));
	}

	return result;
}

HybridBoxes
shrink_out(const HybridBoxes& box, const HybridFloatVector& epsilon)
{
	HybridBoxes result;

	for (HybridBoxes::const_iterator loc_it = box.begin(); loc_it != box.end(); ++loc_it) {
		HybridFloatVector::const_iterator epsilon_it = epsilon.find(loc_it->first);
		ARIADNE_ASSERT_MSG(epsilon_it != epsilon.end(),"The location " << loc_it->first.name() << " is not present in the epsilon map.");
		result.insert(std::pair<DiscreteState,Box>(loc_it->first,loc_it->second.shrink_out(epsilon_it->second)));
	}

	return result;
}

Box
project(const HybridBoxes& box, const std::vector<uint>& dimensions)
{
	Box result = Box::empty_box(dimensions.size());

	for (HybridBoxes::const_iterator loc_it = box.begin(); loc_it != box.end(); ++loc_it)
		result = hull(result,loc_it->second.project(dimensions));

	return result;
}

HybridBoxes
unbounded_hybrid_box(const HybridSpace& hspace)
{
	HybridBoxes result;

	for (HybridSpace::const_iterator space_it = hspace.begin(); space_it != hspace.end(); ++space_it)
		result.insert(std::pair<DiscreteState,Box>(space_it->first,unbounded_box(space_it->second)));

	return result;
}



} // namespace Ariadne
