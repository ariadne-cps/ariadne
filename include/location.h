/***************************************************************************
 *            location.h
 *
 *  Sun Mar 28 10:21:32 2004
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

#ifndef _LOCATION_H
#define _LOCATION_H

#include <string>	
#include <list>

#include "classes.h"

namespace Ariadne {
	
/*! \typedef LocationID
 *  \brief It's the type of the location's univocal identifier. 
 */	
typedef unsigned int LocationID; 
	
/*! \brief This class represents locations.
 *
 * It associates to each location an univocal identifier, a 
 * vector field and an invariant.
 */
class Location
{
		/*! \brief The location's name.	*/
		std::string name;	
	
		/*! \brief Univocal identifier of the location. */
		LocationID id;   

		/*! \brief The smallest unused identifier. */
		static LocationID smallest_unused_id;
	
		/*! \brief The vector field which rules the flow in 
		 * the location. */
		VectorField *field;
	
		/*! \brief The invariant region of the location. */ 
		ASet *invariant;

		/*! \brief The list of the arcs leaving the 
		 * location.*/
		std::list<LeavingArc *> leaving_arcs;
		
	public:
		/*! \brief This is a location class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a Location.
		 * \param name is the name of the location.
		 * \param field is the flow of the location.
		 * \param inv is the invariant of the location.
		 */
		Location(const std::string &name, VectorField *field, 
				ASet *inv);
		
		/*! \brief This is the destructor of the class location.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * \a Location.
		 */
		~Location();
	
		/*! \brief Returns the location's vector field.
		 *
		 * \return The vector field of the location.
		 */
		const VectorField* get_vectorfield() const; 
		
		/*! \brief Returns the location's name.
		 *
		 * \return The name of the location.
		 */
		const std::string& get_name() const; 
	
		/*! \brief Add a leaving arc to the location.
		 *
		 * \param dest is the destination of the new leaving arc.
		 * \param act is the activation region of the leaving arc.
		 * \param reset is the reset of the leaving arc.
		 * \return A pointer to the new leaving arc.
		 */
		const LeavingArc* add_a_leaving_arc(Location *dest, 
				ASet *act, Map *reset);
		
		/*! \brief Initialize location's class.
		 *
		 * Initialize the static members of location 
		 * class. The objects of the class Location
		 * can NOT be initialized before the call
		 * this method.
		 */
		void Init(); 
		
		/*! \brief Returns the location's invariant.
		 *
		 * \return The invariant of the location.
		 */
		const ASet* get_invariant() const;
	
		/*! \brief Returns the location's leaving 
		 * arcs.
		 *
		 * \return The list of the arcs leaving the 
		 * location.
		 */
		const std::list<LeavingArc *>& get_leaving_arcs() const;
		
		/*! \brief Checks if two locations are the same.
		 *
		 * This methods checks if two locations are the same.
		 * \param a is the first location.
		 * \return \a TRUE if \a a is equal to the current 
		 * object, \a TRUE otherwise.
		 */
		bool operator==(const Location &a) const;
		
		/*! \brief Checks if two locations are different.
		 *
		 * This methods checks if two locations are different.
		 * \param a is the first location.
		 * \return \a FALSE if \a a is equal to the current 
		 * object, \a TRUE otherwise.
		 */
		bool operator!=(const Location &a) const; 
};

}

#endif /* _LOCATION_H */
