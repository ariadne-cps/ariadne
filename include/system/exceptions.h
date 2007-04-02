/***************************************************************************
 *            system/exceptions.h
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
 
/*! \file system/exceptions.h
 *  \brief Exceptions, error handling and assertions for the System module.
 */

#ifndef ARIADNE_SYSTEM_EXCEPTIONS_H
#define ARIADNE_SYSTEM_EXCEPTIONS_H

#include <stdexcept>
#include <iosfwd>

#include "../geometry/exceptions.h"

namespace Ariadne {
  namespace System {
  
    //@{ \name Exceptions
    template<class M, class S> inline
    void check_argument_dimension(const M& m1, const S& s2, const char* where="") {
      if(m1.argument_dimension()!=s2.dimension()) { throw Geometry::IncompatibleDimensions(where); }
    }

    template<class M, class S> inline
    void check_result_dimension(const M& m1, const S& s2, const char* where="") {
      if(m1.result_dimension()!=s2.dimension()) { throw Geometry::IncompatibleDimensions(where); }
    }
    //@}
  
  }
}

#endif /* ARIADNE_SYSTEM_EXCEPTIONS_H */
