/****************************************************************************
 *            arc.h
 *
 *  Fri Apr  2 20:12:20 2004
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
 *  Foundation, Inc., 59 Temple Place - Suite 330,
 */

#ifndef _ARC_H
#define _ARC_H

#include "classes.h"

namespace Ariadne {
	
/*! \brief Represents arcs leaving a location.
 *
 * Its members and a source location identify an arc.
 */
class LeavingArc 
{
	protected:
		/*! \brief The destination of the leaving arc. */
		Location *destination;   
	
		/*! \brief The activation region of the leaving arc. */ 
		ASet *activation; 

		/*! \brief The reset of the leaving arc. */ 
		Map *reset;
	
	public:
		/*! \brief This is a \a LeavingArc class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * leaving_arc.
		 * \param dest is the destination of the current leaving arc.
		 * \param act is the activation region of the current leaving 
		 * arc.
		 * \param reset is the reset of the current leaving arc.
		 */
		LeavingArc(Location *dest, ASet *act, Map *reset);
		
		/*! \brief This is the destructor of the class \a LeavingArc.
		 * 
		 * This destructor deletes in a safe way an object of the class 
		 * \a LeavingArc.
		 */
		~LeavingArc();
	
		/*! \brief Return the destination of the arc.
		 * 
		 * This method return the destination of the
		 * arc.
		 * \return The destination of the arc.
		 */
		Location* get_destination() const; 
	
		/*! \brief Return the activation region of the arc.
		 * 
		 * This method return the activation region
		 * of the leaving arc.
		 * \return The activation region of the arc.
		 */
      		const ASet* get_activation() const;

		/*! \brief Return the reset.
		 * 
		 * This method return the reset function of the
		 * leaving arc.
		 * \return The reset function of the
		 * leaving arc.
		 */
		const Map* get_reset() const;
 };

/*! \brief Represents arcs.
 *
 * It adds to the \a LeavingArc class's members the source location
 * of the arc.
 */
class Arc : public LeavingArc
{
	
		/*! \brief This is the source of the current arc. */ 
		Location *source;   
	
	public:
		/*! \brief This is a arc class constructor.
		 *
		 * This constructor initializes the object of the class arc.
		 * @see LeavingArc()
		 * \param source is the source of the current arc.
		 * \param dest is the destination of the current arc.
		 * \param act is the activation region of the current arc.
		 * \param reset is the reset of the current arc.
		 */
		Arc(Location *source, Location *dest, 
				ASet *act, Map *reset); 

		/*! \brief This is the destructor of the class arc.
		 * 
		 * This destructor deletes in a safe way an object of the class 
		 * arc.
 		 */
		~Arc(); 
	
		/*! \brief Return the source of the arc.
		 *
		 * This method return the source of the
		 * arc.
		 * \return The source of the arc.
		 */
		Location* get_source() const; 

};

}
 
#endif /* _ARC_H */
