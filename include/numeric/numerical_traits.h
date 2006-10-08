/***************************************************************************
 *            numerical_traits.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file numerical_traits.h
 *  \brief Traits classes to define properties of numerical types.
 */

#ifndef _ARIADNE_NUMERICAL_TRAITS_H
#define _ARIADNE_NUMERICAL_TRAITS_H

#include <string>

namespace Ariadne {
  namespace Numeric {
    
    /* numerical traits */
    /*! \brief Tags a class representing a ring. */
    class ring_tag { };
    /*! \brief Tags a class representing a field. */
    class field_tag { };
      
    /*! \brief Typedef's describing a numerical type. */
    template<typename T> class numerical_traits { };

    template<> class numerical_traits<double> {
     public:
      typedef field_tag algebraic_category;
      typedef double field_extension_type;
    };
  
    //! \name Numerical type description.
    //@{
    /*! \brief The name of class T. */
    template<typename T> inline std::string name();

    //! \name Standard conversion operations. (Deprecated) 
    /*! \brief Approximate \a x by an element of Res. */
    template<typename Res, typename Arg> inline Res convert_to(const Arg& x) { return Res(x); }
    
    /*! \brief Approximate \a x by an element of Res with accuracy \a e. */
    template<typename Res, typename Arg, typename Err> Res approximate(const Arg& x, const Err& e);
    //@}
  }   
}
  

#endif /* _ARIADNE_NUMERICAL_TRAITS_H */
