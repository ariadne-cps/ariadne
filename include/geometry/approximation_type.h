/***************************************************************************
 *            approx_type.h
 *
 *  Wed Apr 21 12:37:01 2004
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
 
#ifndef _APPROX_TYPE_H
#define _APPROX_TYPE_H

#include "../base/numerical_type.h"

namespace Ariadne {
namespace Geometry {

enum ApproxKind {
  OVER, /*!< Over-approximation */
  UNDER, /*!< Under-approximation */
  OUTER, /*!< Outer-approximation */
  INNER, /*!< Inner-approximation */
  NONE /*!< No approximation */
};

/*! \brief Defines the approximation type */ 
template < typename R>
class ApproximationType {
  public:
    typedef R Real;
  
  private:

    /*! \brief Indicates the type of approximation */
    ApproxKind _kind;

    /*! \brief Indicates the value of approximation */
    Real _value; 

  public:
    
    /*! \brief A costructor for the ApproximationType class.*/ 
    ApproximationType() {

      /* an approx type has NONE kind by default */
      this->_kind=NONE;
      this->_value=0.0;
    }
    
    /*! \brief A costructor for the ApproximationType class.
     *
     * \param orig is the object used a framework for
     * the new object.*/ 
    ApproximationType(const ApproximationType<Real> &orig) {
      this->_kind=orig._kind;
      this->_value=orig._value;
    }
    
    /*! \brief A costructor for the ApproximationType class. 
     *
     * This method is a costructor of the ApproximationType class.
     * \param kind is the kind of the new ApproximationType object.
     * \param value is the value of the approximation.
     */
    ApproximationType(const ApproxKind kind, const Real &value) {
      this->_kind=kind;
      this->_value=value;
    }

    /*! \brief Returns the object's value.
     *
     * \return The value of the approximation.
     */
    Real &value() {
      return (this->_value);
    }

    /*! \brief Returns the kind of approximation.
     *
     * \return The object's kind.
     */
    ApproxKind &kind() {
      return (this->_kind);
    }

    /*! \brief Assigns an ApproximationType's object
     * to an other.
     *
     * \param orig is the object which has to be 
     * copied.
     * \return A reference to the new object.
     */
    const ApproximationType& operator=(const ApproximationType &orig) {
      this->_kind=orig._kind;
      this->_value=orig._value;
      
      return *this;
    }

    
};

}

}

#endif /* _APPROX_TYPE_H */
