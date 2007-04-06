/***************************************************************************
 *            combinatoric/exceptions.h
 *
 *  Copyright  2005-7  Pieter Collins, Alberto Casagrande
 *  Email  Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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
 
/*! \file combinatoric/exceptions.h
 *  \brief Exceptions, error handling and assertions for the Combinatoric module.
 */

#ifndef ARIADNE_COMBINATORIC_EXCEPTIONS_H
#define ARIADNE_COMBINATORIC_EXCEPTIONS_H

#include <stdexcept>

#include "../base/exceptions.h"
#include "../geometry/exceptions.h"

namespace Ariadne {
  namespace Combinatoric {
      
    //@{ \name Exceptions
    /*! \brief The dimensions of two geometric objects are incompatible. */
    struct IncompatibleDimensions : public std::runtime_error {
      IncompatibleDimensions(const std::string& str) : std::runtime_error(str) { }
    };

    /*! \brief The coordinate in  a geometric object was invalid. */
    struct InvalidCoordinate : public std::runtime_error {
      InvalidCoordinate(const std::string& str) : std::runtime_error(str) { }
    };

    /*! \brief An invalid index to the set of vertices. */
    struct InvalidVertex : public std::runtime_error {
      InvalidVertex(const std::string& str) : std::runtime_error(str) { }
    };

    /*! \brief Attempting to perform an operation on an unbounded set. */
    struct UnboundedSet : public std::runtime_error {
      UnboundedSet(const std::string& str) : std::runtime_error(str) { }
    };

    
    //@}

  }
}

#endif /* ARIADNE_COMBINATORIC_EXCEPTIONS_H */
