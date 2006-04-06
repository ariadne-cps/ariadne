/***************************************************************************
 *            discrete_location.h
 *
 *  Thu Jun 17 15:46:37 2004
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
 
#ifndef _DISCRETE_LOCATION_H
#define _DISCRETE_LOCATION_H

#include <string>

namespace Ariadne {

namespace HybridDefinitions {

template < typename LDT > class HybridAutomaton;
  
/*! \typedef DiscreteLocationID
 *  \brief It's the type of the discrete location's univocal identifier. 
 */  
typedef size_t DiscreteLocationID;  

/*! \brief The discrete location type. */
template <typename VF >
class DiscreteLocation{

  public:
    typedef VF vector_typeField;
    typedef typename vector_typeField::vector_typeFieldMap vector_typeFieldMap;
    typedef typename vector_typeField::DenotableSet DenotableSet;
    typedef typename DenotableSet::BasicSet BasicSet;
    typedef typename BasicSet::state_type state_type;
    typedef typename state_type::real_type real_type;
    
  
    template < typename LDT >
    friend class HybridAutomaton;
      
    template < typename LDS >
    friend class HybridDenotableSet;
  private:
    /*! \brief The discrete location's identificator. */
    DiscreteLocationID _id;
  
    /*! \brief The discrete location's name. */
    std::string _name;
  
    /*! \brief The discrete location's vector field. */
    vector_typeField _vfield;
  
    /*! \brief The discrete location's invariant. */
    DenotableSet _invariant;
  
    /*! \brief Indicates if this location has been inserted into an automaton. */
    bool _in_automaton;
  
        
    /*! \brief Sets the discrete location's id. 
     *
     * \param new_id is the new location's id.
     */
    inline void _set_id(const DiscreteLocationID new_id) {
      this->_id=new_id;
      
      this->_in_automaton=true;
    }
  
  public:  
  
    /*! \brief This is a discrete location class constructor.
     *  
     * This constructor initializes the object of the 
     * discrete location class.
     * \param name is the name of the discrete location.
     * \param vfield is the location's vector field.
     */
    DiscreteLocation(const std::string &name, 
                    const vector_typeField &vfield):
        _name(name), _vfield(vfield), 
        _invariant(vfield.dimension()), _in_automaton(false) {}
          
    /*! \brief This is a discrete location class constructor.
     *  
     * This constructor initializes the object of the 
     * discrete location class.
     * \param name is the name of the discrete location.
     * \param vfield is the location's vector field.
     * \param invariant is the location's invariant.
     */
    DiscreteLocation(const std::string &name, 
                    const vector_typeField &vfield, 
                    const DenotableSet &invariant):
          
        _name(name), _vfield(vfield), 
        _invariant(invariant), _in_automaton(false) {
          
          
      if (vfield.dimension()!=invariant.dimension()) {
        throw std::invalid_argument("The invariant and the vector field have different space dimensions.");
      }
    }
      
    /*! \brief This is a discrete location class constructor.
     *  
     * This constructor initializes the object of the 
     * discrete location class.
     * \param name is the name of the discrete location.
     * \param vfield is the location's vector field.
     * \param invariant is the location's invariant.
     */
    DiscreteLocation(const std::string &name, 
                    const vector_typeField &vfield, 
                    const BasicSet &invariant):
          
        _name(name), _vfield(vfield), 
        _in_automaton(false) {
          
          
      if (vfield.dimension()!=invariant.dimension()) {
        throw std::invalid_argument("The invariant and the vector field have different space dimensions.");
      }
      
      DenotableSet DS_inv(invariant);
      
      this->_invariant=DS_inv;
    }
    
    /*! \brief This is a discrete location class constructor.
     *  
     * This constructor initializes the object of the 
     * discrete location class.
     * \param orig is the original object.
     */
    DiscreteLocation(
          const DiscreteLocation<VF> &orig): 

      _id(orig._id), _name(orig._name), _vfield(orig._vfield),
      _invariant(orig._invariant), _in_automaton(orig._in_automaton) {}
    
    /*! \brief Gets the discrete location's id. 
     *
     * \return The location's id.
     */
    inline const DiscreteLocationID &id() const{
      return this->_id;
    }
    
    /*! \brief Returns the discrete location's name. 
     *
     * \return The discrete location's name.
     */
    inline const std::string& name() const {
      return this->_name;  
    }
    
    /*! \brief Tests if this location is part of an automaton. 
     *
     * \return \a true if this location is part of an automaton, \a false
     * otherwise.
     */
    inline bool in_automaton() const {
      return this->_in_automaton;  
    }
    
    /*! \brief Returns the discrete location's vector field. 
     *
     * \return The discrete location's vector field.
     */
    inline const vector_typeField& vector_field() const {
      return this->_vfield;  
    }
    
    /*! \brief Returns the discrete location's invariant. 
     *
     * \return The discrete location's invariant.
     */
    inline const DenotableSet& invariant() const{
      return this->_invariant;  
    }
    
    /*! \brief Copies a discrete location. 
     *
     * \param orig is the original object.
     * \return The discrete location reference.
     */
    inline const DiscreteLocation<VF>& operator=(
        const DiscreteLocation<VF> &orig){
        
      this->_id=orig._id;
      this->_name=orig._name;
      this->_vfield=orig._vfield;
      this->_invariant=orig._invariant;
      this->_in_automaton=orig._in_automaton;
  
      return *this;      
    }
    
};
  
}
}

#endif /* _DISCRETE_LOCATION_H */
