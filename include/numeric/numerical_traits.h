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

namespace Ariadne {
  namespace Numeric {
    
    /* numerical traits */
    /*! \brief Tags a class representing a ring. */
    class ring_tag { };
    /*! \brief Tags a class representing a field. */
    class field_tag { };
      
    /*! \brief Typedef's describing a numerical type. */
    template<typename T> class numerical_traits;
  }
}
  
namespace Ariadne {
  namespace Base {
    template<typename R> inline std::string name() { 
      throw std::runtime_error("name(): Unknow type");
    }
  
  }
}

#endif /* _ARIADNE_NUMERICAL_TRAITS_H */
