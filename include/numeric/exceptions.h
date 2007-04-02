/***************************************************************************
 *            numeric.except.h.h
 *
 *  Copyright  2005-6  Pieter Collins, Alberto Casagrande
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
 
/*! \file numeric/exceptions.h
 *  \brief Exceptions, error handling and assertions for the Numeric module.
 */

#ifndef ARIADNE_NUMERIC_EXCEPTIONS_H
#define ARIADNE_NUMERIC_EXCEPTIONS_H

#include <stdexcept>
#include <iosfwd>

#include "base/types.h"

namespace Ariadne {
  namespace Numeric {
  
    //@{ \name Exceptions
    /*! \brief A division by zero has occurred. */
    class DivideByZeroException : public std::exception { };
    //@}
    
  }
}


#endif /* ARIADNE_NUMERIC_EXCEPTIONS_H */
