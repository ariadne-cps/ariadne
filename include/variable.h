/***************************************************************************
 *            variable.h
 *
 *  Sun Aug 01 18:55:31 2004
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
 
#ifndef _VARIABLE_H
#define _VARIABLE_H

#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
	
namespace Ariadne {	

/*! \typedef VariableID
 *  \brief It's the type of the variable's univocal identifier. 
 *
 *  This type is used to distinquish two different variables
 *  during the "compilation" of the system definition also.
 */	
typedef unsigned int VariableID; 

/* TODO: add comments */
typedef enum{
	REAL_VAR,
	SYMBOLIC_CONST
} VariableType;

	
/*! \brief This class associates to each variable, an id and a
 * name.
 *
 * This class is used to distinquish two different variables
 * during the "compilation" of the system definition also.
 */
class Variable {
		/*! \brief The variable's name */
		std::string name;

		/*! \brief The variable's unique id */
		VariableID id;

		/*! \brief The variable's type 
		 *
		 * This member will be used in the future to 
		 * distinguish between a real variable and 
		 * symbolic constant.
		 */
		VariableType type;
		
		/*! \brief The smallest unused variable id */
		static VariableID smallest_unused_id;

	public:
		/*! \brief This is a variable class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a Variable.
		 * \param name is the name of the variable.
		 */
		Variable(const std::string &name);

		/*! \brief This is a variable class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a Variable.
		 * \param name is the name of the variable.
		 * \param id is the variable's unique identifier.
		 */
		Variable(const std::string &name, VariableID id);
		
		/*! \brief This is a variable class constructor.
		 *
		 * This constructor initializes the object of the class 
		 * \a Variable.
		 * \param name is the name of the variable.
		 * \param type is the variable's type.
		 */
		Variable(const std::string &name, VariableType type);

		/*! \brief This is the destructor of the class variable.
		 *
		 * This destructor deletes in a safe way an object of the class 
		 * \a Variable.
		 */
		~Variable();

		/*! \brief Check if the variable has id \a id.
		 *
		 * \param id is the potential variable's id.
		 * \return FALSE if the variable's id is different
		 * from \a id, TRUE otherwise.
		 */
		bool has_id(const VariableID id) const;
		
		/*! \brief Check if the variable has name \a name.
		 *
		 * \param name is the potential variable's name.
		 * \return FALSE if the variable's name is different
		 * from \a name, TRUE otherwise.
		 */
		bool has_name(const std::string &name) const;

		/*! \brief Set new name and id for a variable.
		 *
		 * \param name is the new name of the variable.
		 * \param id is the new variable's unique identifier.
		 */
		void set(const std::string &name, VariableID id);

};

typedef std::vector<Variable> VecVariable;

}

#endif /* VARIABLE_H */
