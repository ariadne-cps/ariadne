/***************************************************************************
 *            basic_maintainer.h
 *
 *  Mon Apr 19 15:29:15 2004
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
 
#ifndef _BASIC_MAINTAINER_H
#define _BASIC_MAINTAINER_H

#include "classes.h"
#include "cluster_list.h"

namespace Ariadne {
	
	
/*! \brief The simplest structure memorizing a region of the space.
 *
 * It's an interface for classes that maintain a region of the 
 * space.
 */
class BasicMaintainer
{
		/*! \brief Stores the changes of the maintain system
		 *
		 * This member stores the changes of the maintain system
		 * from the last call of the method \a clean_changes.
		 */
		_Maintain_Changes *changes;

	public:
		/*! \brief This is a \a BasicMaintainer class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a BasicMaintainer.
		 * \param dim is the dimention of the space that should 
		 * maintained. 
		 */
		BasicMaintainer(const unsigned int dim);
	
		/*! \brief This is a \a BasicMaintainer class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a BasicMaintainer using as template the object \a templ.
		 * \param templ is the template object for the constructor. 
		 */
		BasicMaintainer(const BasicMaintainer &templ);
		
		/*! \brief Returns the grid box of maintained region.
		 *
		 * This method returns the grid box of maintained region.
		 * \return The cluster list of maintained region.
		 */
		virtual ClusterList* get_maintained() = 0;
		 
		 /*! \brief Returns the grid box of a subset of the 
		 * maintained region.
		 *
		 * This method returns the grid box of the maintained 
		 * region included in the space box passed as 
		 * parameter. 
		 * \param sbox indicates the space box that should
		 * be represented in the output.		 
		 * \return The cluster list of maintained 
		 * region contained into the space box \a sbox.
		 */
		virtual ClusterList* get_maintained(const SpaceBox &sbox) = 0;
		 
		 /*! \brief Copies the object.
		 *
		 * This method copies the object and returns
		 * a pointer to the copy.
		 * \return A pointer to the copy of the object.
		 */
		virtual BasicMaintainer* copy() = 0;
};

}

#endif /* _BASIC_MAINTAINER_H */
