/***************************************************************************
 *            trace.h
 *
 *  Thu Aug 31 16:20:01 2004
 *  Copyright  2004  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
#ifndef _TRACE_H
#define _TRACE_H

#include <list>

#include "location.h"

namespace Ariadne {

class LocationTrace {
		Location *location;
		double time;
	public:
		LocationTrace(Location *location, double time);

		~LocationTrace();

		Location *get_location() const;

		double *get_time() const;
};

typedef std::list<LocationTrace> Trace;

}

#endif /* _TRACE_H */
