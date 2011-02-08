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
	ARIADNE_ASSERT_MSG(box1.size() == box2.size(),"The two hybrid boxes must have the same number of locations.");

	HybridBoxes result;
	for (HybridBoxes::const_iterator box1_it = box1.begin(); box1_it != box2.end(); ++box1_it) {
		HybridBoxes::const_iterator box2_it = box2.find(box1_it->first);
		ARIADNE_ASSERT_MSG(box2_it != box2.end(),"The location " << box1_it->first.name() << " is not present in both hybrid boxes.");
		result.insert(std::pair<DiscreteState,Box>(box1_it->first,hull(box1_it->second,box2_it->second)));
	}

	return result;
}


} // namespace Ariadne
