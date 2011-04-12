/***************************************************************************
 *            verification_input.cc
 *
 *  Copyright 2011  Luca Geretti
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

#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <string>

#include "verification_input.h"

namespace Ariadne {


VerificationInput::VerificationInput(
		HybridAutomaton& system,
		HybridImageSet& initial_set,
		HybridBoxes& domain) :
		_system(system),
		_initial_set(initial_set),
		_domain(domain)
{
	_check_fields();
}


void VerificationInput::_check_fields() const
{
	HybridSpace hspace = _system.state_space();
	for (HybridImageSet::const_iterator it = _initial_set.begin(); it != _initial_set.end(); ++it) {
		HybridSpace::const_iterator hspace_it = hspace.find(it->first);
		ARIADNE_ASSERT_MSG(hspace_it != hspace.end(),
						   "The location " << it->first.name() << "is not present into the system.");
		ARIADNE_ASSERT_MSG(hspace_it->second == it->second.dimension(),
						   "The dimension of the continuous space in the initial set for location " << it->first.name() << " does not match the system space");
	}
	for (HybridSpace::const_iterator hspace_it = hspace.begin(); hspace_it != hspace.end(); ++hspace_it) {
		HybridBoxes::const_iterator domain_it = _domain.find(hspace_it->first);
		ARIADNE_ASSERT_MSG(domain_it != _domain.end(),
						   "The location " << hspace_it->first.name() << "is not present into the domain.");
		ARIADNE_ASSERT_MSG(hspace_it->second == domain_it->second.dimension(),
						   "The dimension of the continuous space in the domain for location " << hspace_it->first.name() << " does not match the system space");
	}
}


std::ostream&
VerificationInput::write(std::ostream& os) const
{
	os << "(System: " << _system << "; Initial set: " << _initial_set << "; Domain: " << _domain << ")";
	return os;
}


SafetyVerificationInput::SafetyVerificationInput(
		HybridAutomaton& system,
		HybridImageSet& initial_set,
		HybridBoxes& domain,
		HybridBoxes& safe_region) :
		VerificationInput(system,initial_set,domain),
	    _safe_region(safe_region)
{
	_check_fields();
}


void SafetyVerificationInput::_check_fields() const
{
	VerificationInput::_check_fields();

	HybridSpace hspace = getSystem().state_space();
	for (HybridSpace::const_iterator hspace_it = hspace.begin(); hspace_it != hspace.end(); ++hspace_it) {
		HybridBoxes::const_iterator safe_it = _safe_region.find(hspace_it->first);
		ARIADNE_ASSERT_MSG(safe_it != _safe_region.end(),
						   "The location " << hspace_it->first.name() << "is not present into the safe region.");
		ARIADNE_ASSERT_MSG(hspace_it->second == safe_it->second.dimension(),
						   "The dimension of the continuous space in the safe region for location " << hspace_it->first.name() << " does not match the system space");
	}
}


std::ostream&
SafetyVerificationInput::write(std::ostream& os) const
{
	os << "(System: " << getSystem() << "; Initial set: " << getInitialSet() << "; Domain: " <<
			getDomain() << "; Safe region: " << _safe_region << ")";
	return os;
}


DominanceVerificationInput::DominanceVerificationInput(
		HybridAutomaton& system,
		HybridImageSet& initial_set,
		HybridBoxes& domain,
		std::vector<uint>& projection) :
		VerificationInput(system,initial_set,domain),
		_projection(projection)
{
}


std::ostream&
DominanceVerificationInput::write(std::ostream& os) const
{
	os << "(System: " << getSystem() << "; Initial set: " << getInitialSet() << "; Domain: " <<
			getDomain() << "; Projection: " << _projection << ")";
	return os;
}



}
