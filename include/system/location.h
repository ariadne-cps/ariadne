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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifndef _LOCATION_H
#define _LOCATION_H

#include <string>
#include <list>


namespace Ariadne {

namespace SystemDescription {

/*! \typedef DiscreteLocationID
 *  \brief It's the type of the discrete location's univocal identifier. 
 */  
typedef unsigned int DiscreteLocationID; 
  
/*! \brief This class represents discrete locations.
 *
 * It associates to each discrete location an univocal identifier, a 
 * vector field and an invariant.
 */
class DiscreteLocation
{
    /*! \brief The discrete location's name.  */
    std::string _name;  
  
    /*! \brief Univocal identifier of the discrete location. */
    DiscreteLocationID _id;   

    /*! \brief The smallest unused identifier. */
    static DiscreteLocationID _smallest_unused_id;
  
    /*! \brief The vector field which rules the flow in 
     * the discrete location. */
    VectorField *_field;
  
    /*! \brief The invariant region of the discrete location. */ 
    Ariadne::Geometry::AbstractDenotableSet *_invariant;

    /*! \brief The list of the discrete transitions leaving the 
     * discrete location.*/
    std::list<DiscreteTransition *> _discrete_transitions;
    
  public:
    /*! \brief This is a discrete location class constructor.
     *
     * This constructor initializes the object of the class 
     * \a DiscreteLocation.
     * \param name is the name of the discrete location.
     * \param field is the flow of the discrete location.
     * \param inv is the invariant of the discrete location.
     */
    DiscreteLocation(const std::string &name, VectorField *field, 
        Ariadne::Geometry::AbstractDenotableSet *inv);
    
    /*! \brief This is the destructor of the discrete 
     * locationi class.
     *
     * This destructor deletes in a safe way an object of the class 
     * \a DiscreteLocation.
     */
    ~DiscreteLocation();
  
    /*! \brief Returns the discrete location's vector field.
     *
     * \return The vector field of the discrete location.
     */
    const VectorField* field() const {
      return this->_field;
    }
    
    /*! \brief Returns the discrete location's name.
     *
     * \return The name of the discrete location.
     */
    const std::string& name() const {

      return this->_name;
    }
  
    /*! \brief Add a discrete transition from a discrete location.
     *
     * \param dest is the destination of the new discrete 
     * transition.
     * \param act is the activation region of the discrete 
     * transition.
     * \param reset is the reset of the discrete transition.
     * \return A pointer to the new discrete transition.
     */
    const DiscreteTransition* add_discrete_transition(
        DiscreteLocation *dest, 
        Ariadne::Geometry::ASet *act, Map *reset);
    
    /*! \brief Initialize discrete location's class.
     *
     * Initialize the static members of discrete location 
     * class. The objects of the class DiscreteLocation
     * can NOT be initialized before the call
     * this method.
     */
    static void Init() { 
      this->_smallest_unused_id=0;
    }

    /*! \brief Returns the discrete location's invariant.
     *
     * \return The invariant of the discrete location.
     */
    const Ariadne::Geometry::AbstractDenotableSet* invariant() const {
      return this->_invariant;
    }
  
    /*! \brief Returns the discrete transitions leaving the
     * discrete location.
     *
     * \return The list of the discrete transitions leaving the 
     * discrete location.
     */
    const std::list<DiscreteTransition *>& discrete_transitions() const {
      return this->_discrete_transitions;
    }
    
    /*! \brief Checks if two discrete locations are the same.
     *
     * This methods checks if two discrete locations are the same.
     * \param a is the first discrete location.
     * \return \a TRUE if \a a is equal to the current 
     * object, \a TRUE otherwise.
     */
    bool operator==(const DiscreteLocation &a) const {
      return (this->_id==a._id);  
    }
    
    /*! \brief Checks if two discrete locations are different.
     *
     * This methods checks if two discrete locations are different.
     * \param a is the first discrete location.
     * \return \a FALSE if \a a is equal to the current 
     * object, \a TRUE otherwise.
     */
    bool operator!=(const DiscreteLocation &a) const {
      return (this->_id!=a._id);  
    } 
 
};

}

}

#endif /* _LOCATION_H */
