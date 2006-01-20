/****************************************************************************
 *            discrete_trans.h
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#ifndef _DISCRETE_TRANSITION_H
#define _DISCRETE_TRANSITION_H

#include <discrete_location.h>

namespace Ariadne {
namespace HybridDefinitions {
    
  
/*! \brief Represents discrete transition leaving a discrete location.*/
template <typename LOC, typename RESET >
class LeavingDiscreteTransition 
{
  public:
    typedef RESET ResetMap;
    typedef LOC DiscreteLocation;
    typedef typename DiscreteLocation::VectorField VectorField;
    typedef typename ResetMap::DenotableSet DenotableSet;
    typedef typename DenotableSet::BasicSet BasicSet;
    typedef typename BasicSet::State State;
    typedef typename State::Real Real;
  
  private:
  
    /*! \brief The destination of the leaving discrete 
     * transition. */
    DiscreteLocation _destination;   
  
    /*! \brief The activation region of the leaving discrete 
     * transition. */ 
    DenotableSet _activation; 

    /*! \brief The reset of the leaving discrete transition. */ 
    ResetMap _reset;
  
  public:
  
    /*! \brief This is a \a LeavingDiscreteTransition class 
     * constructor.
     *
     * This constructor initializes the object of the  
     * leaving discrete transition class.
     * \param dest is the destination of the current leaving 
     * discrete transition.
     * \param act is the activation region of the current leaving 
     * discrete transition.
     * \param reset is the reset of the current leaving discrete 
     * transition.
     */
    LeavingDiscreteTransition(const DiscreteLocation &dest, 
            const DenotableSet &act, const ResetMap &reset):
        _destination(dest), _activation(act), _reset(reset) {}
    
    /*! \brief This is a \a LeavingDiscreteTransition class 
     * constructor.
     *
     * This constructor initializes the object of the  
     * leaving discrete transition class.
     * \param dest is the destination of the current leaving 
     * discrete transition.
     * \param act is the activation region of the current leaving 
     * discrete transition.
     * \param reset is the reset of the current leaving discrete 
     * transition.
     */
    LeavingDiscreteTransition(const DiscreteLocation &dest, 
            const BasicSet &act, const ResetMap &reset):
        _destination(dest),_reset(reset) {
    
      DenotableSet DS_act(act);
          
          
      this->_activation=DS_act;
          
    }
          
    /*! \brief This is a \a LeavingDiscreteTransition class 
     * constructor.
     *
     * This constructor initializes the object of the  
     * leaving discrete transition class.
     * \param orig is the original copy of the leaving arc.
     */
    LeavingDiscreteTransition(const
        LeavingDiscreteTransition< LOC , RESET > &orig):
        _destination(orig._destination), _activation(orig._activation),
        _reset(orig._reset) {}
  
    /*! \brief Return the destination of the discrete transition.
     * 
     * This method return the destination of the discrete 
     * transition.
     * \return The destination of the discrete transition.
     */
    inline const DiscreteLocation &destination() const { 
      return this->_destination;
    }
  
    /*! \brief Return the activation region of the discrete 
     * transition.
     * 
     * This method return the activation region of the leaving 
     * discrete transition.
     * \return The activation region of the discrete transition.
     */
        inline const DenotableSet &activation() const{ 
      return this->_activation;
    }

    /*! \brief Return the reset.
     * 
     * This method return the reset function of the leaving 
     * discrete transition.
     * \return The reset function of the leaving discrete 
     * transition.
     */
    inline const ResetMap &reset() const { 
      return this->_reset;
    }
    
    /*! \brief Copies a leaving arc.
     *
     * This method copies a leaving arc.
     * \param orig is the original copy of the leaving arc.
     */
    inline const LeavingDiscreteTransition< LOC , RESET >& 
        operator=(const LeavingDiscreteTransition< LOC , RESET > &orig) {
      
      this->_activation=orig._activation;
      this->_destination=orig._destination;
          
      this->_reset=orig._reset;
          
      return *this;
    }
 };

/*! \brief Represents discrete transitions.
 *
 * It adds to the \a LeavingDiscreteTransition class's members the source 
 * discrete location of the discrete transition.
 */
template < typename LOC, typename RESET >
class DiscreteTransition
{
  public:
    typedef RESET ResetMap;
    typedef LOC DiscreteLocation;
    typedef typename DiscreteLocation::VectorField VectorField;
    typedef typename ResetMap::DenotableSet DenotableSet;
    typedef typename DenotableSet::BasicSet BasicSet;
    typedef typename BasicSet::State State;
    typedef typename State::Real Real;
  
  private:
    /*! \brief This is the source of the current discrete 
     * transition. */ 
    DiscreteLocation _source;   
  
    /*! \brief The destination of the leaving discrete 
     * transition. */
    DiscreteLocation _destination;   
  
    /*! \brief The activation region of the leaving discrete 
     * transition. */ 
    DenotableSet _activation; 

    /*! \brief The reset of the leaving discrete transition. */ 
    ResetMap _reset;  
    
  public:
    /*! \brief This is a discrete transition class constructor.
     *
     * This constructor initializes the object of the discrete 
     * transitioni class.
     * @see LeavingDiscreteTransition()
     * \param source is the source of the current discrete 
     * transition.
     * \param dest is the destination of the current discrete 
     * transition.
     * \param act is the activation region of the current 
     * discrete transition.
     * \param reset is the reset of the current discrete 
     * transition.
     */
    DiscreteTransition(const DiscreteLocation &source, 
        const DiscreteLocation &dest, 
        const DenotableSet &act, ResetMap reset):
              _source(source), _destination(dest), 
              _activation(act), _reset(reset) {}

  
    /*! \brief This is a \a LeavingDiscreteTransition class 
     * constructor.
     *
     * This constructor initializes the object of the  
     * leaving discrete transition class.
     * \param orig is the original copy of the leaving arc.
     */
    DiscreteTransition(const
        DiscreteTransition< LOC, RESET > &orig):
        _source(orig._source), _destination(orig._destination), 
        _activation(orig._activation), _reset(orig._reset) {}
                
    /*! \brief Return the source of the discrete transition.
     *
     * This method return the source of the discrete transition.
     * \return The source of the discrete transition.
     */
    inline const DiscreteLocation &source() const {
      this->_source;
    }      

      /*! \brief Return the destination of the discrete transition.
     * 
     * This method return the destination of the discrete 
     * transition.
     * \return The destination of the discrete transition.
     */
    inline const DiscreteLocation &destination() const { 
      return this->_destination;
    }
  
    /*! \brief Return the activation region of the discrete 
     * transition.
     * 
     * This method return the activation region of the leaving 
     * discrete transition.
     * \return The activation region of the discrete transition.
     */
        inline const DenotableSet &activation() const{ 
      return this->_activation;
    }

    /*! \brief Return the reset.
     * 
     * This method return the reset function of the leaving 
     * discrete transition.
     * \return The reset function of the leaving discrete 
     * transition.
     */
    inline const ResetMap &reset() const{ 
      return this->_reset;
    }
    
    /*! \brief Copies a leaving arc.
     *
     * This method copies a leaving arc.
     * \param orig is the original copy of the leaving arc.
     */
    inline const DiscreteTransition< LOC , RESET >& 
        operator=(const DiscreteTransition< LOC , RESET > &orig) {
      
      this->_activation=orig._activation;
      this->_destination=orig._destination;
          
      this->_reset=orig._reset;
      this->_source=orig._source;
            
      return *this;
    }
};

}
}
 
#endif /* _DISCRETE_TRANSITION_H */
