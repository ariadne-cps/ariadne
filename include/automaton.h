/***************************************************************************
 *            automaton.h
 *
 *  Tue Mar 23 14:12:31 2004
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
 
#ifndef _AUTOMATON_H
#define _AUTOMATON_H

#include <string>
#include <list>

#include "classes.h"
#include "io_error.h"
#include "variable.h"

namespace Ariadne {	

/*! \typedef AutomatonID
 *  \brief It's the type of the automaton's univocal identifier. 
 */	
typedef unsigned int AutomatonID;

/*! \class Automaton
 *  \brief This is the class of automata.
 * 
 * It represents automata and its methods evaluate both reachability
 * properties and region.
 */
class Automaton
{
		/*! \brief Automaton's name. */
		std::string name;

		/*! \brief Automaton's identificator*/
		AutomatonID id;

		/*! \brief The smallest unused automaton's identificator*/
		static AutomatonID smallest_unused_id;

		/*! \brief The list of the automaton's locations.
		 *
		 * This is the list of the automaton's locations.
		 * Each location is link to the others by arcs
		 */	
		std::list<Location *> *locations;
	
		/*! \brief The initial regions.
		 *
		 * This is the list of the automaton's initial regions. 
		 */
		std::list<HSet> *i_regions;

		/*! \brief The variable vector.
		 *
		 * This vector maintains the name, the id and the row
		 * on the equation system of each automaton's variable
		 * (ex. if the variable "x" as id 5 and its flow 
		 * equations are at the 7th row of the system, it will
		 * be 7th element of v_variable.
		 */
		VecVariable *v_variable;

		/* TODO: write the method set_v_variable */
	public:
		/*! \brief Initialize automaton's class.
		 *
		 * Initialize the static members of automaton's
		 * class. The objects of the class Automaton
		 * can NOT be initialized before the call
		 * this method.
		 */
		void Init(); 
		
		/*! \brief This is a automaton class constructor.
		 *  
		 * This constructor initializes the object of the class 
		 * automaton.
		 * \param name is the name of the automaton.
		 * \param locations is the automaton's location graph.
		 * \param init is the list of the automaton's initial
		 * regions.
		 * \param v_variable maintains the variable's names and
		 * identificators.
		 */
		Automaton(const std::string &name, 
				std::list<Location *> *locations, 
				std::list<HSet> *init,
				VecVariable *v_variable);
		
		/*! \brief This is a automaton class constructor.
		 *  
		 * This constructor initializes the object of the class 
		 * automaton.
		 * \param name is the name of the automaton.
		 * \param init is the list of the automaton's initial
		 * regions.
		 */
		Automaton(const std::string &name, std::list<HSet> *init);
		
		/*! \brief This is a automaton class constructor.
		 *  
		 * This constructor initializes the object of the class 
		 * automaton.
		 * \param name is the name of the automaton.
		 */
		Automaton(const std::string &name); 
	
		/*! \brief  This is the destructor of the class automaton.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * automaton.
		 */
     		~Automaton();

		/*! \brief  Load an automaton definition.
		 *
		 * This method loads an automaton definition from a file.
		 * The name of the file containing the definition is
		 * stored into the parameter.
		 * \param file is the file name.
		 */
     		Read_Error load_automaton(const std::string &file);
		
		/*! \brief  Save an automaton definition.
		 *
		 * This method saves an automaton definition into a file.
		 * The name of the file containing the definition is
		 * stored into the parameter.
		 * \param file is the file name.
		 */
     		Write_Error save_automaton(const std::string &file);

		/*! \brief Returns the automaton's name.
		 *
		 * \return The name of the automaton.
		 */
		const std::string &get_name();

		/*! \brief Checks if two automata are the same.
		 *
		 * This methods checks if two automata are the same.
		 * \param a is the first automaton.
		 * \return \a TRUE if \a a  is equal to the current
		 * object, \a FALSE otherwise.
		 */
		bool operator==(const Automaton &a); 
		
		/*! \brief Checks if two automata are different.
		 *
		 * This methods checks if two automata are different.
		 * \param a is the first automaton.
		 * \return \a FALSE if \a a is equal to the current
		 * object, \a TRUE otherwise.
		 */
		bool operator!=(const Automaton &a); 
		
		friend class Solver;
};

}

#endif /* _AUTOMATON_H */
