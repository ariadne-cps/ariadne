/***************************************************************************
 *            system/exceptions.h
 *
 *  Copyright  2005-7  Pieter Collins, Alberto Casagrande
 *  Email  Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or61
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
    struct InvalidParameters : public std::runtime_error { InvalidParameters(const std::string& what) : std::runtime_error(what) { } };
    using Geometry::IncompatibleDimensions; 
  }
}

#define ARIADNE_CHECK_ARGUMENT_DIMENSION(map,set,func) \
  { if((map).argument_dimension()!=(set).dimension()) { ARIADNE_THROW(IncompatibleDimensions,func,#map".argument_dimension()="<<(map).argument_dimension()<<", "#set"="<<set); } }

#define ARIADNE_CHECK_RESULT_DIMENSION(map,set,func) \
{ if((map).result_dimension()!=(set).dimension()) { ARIADNE_THROW(IncompatibleDimensions,func,#map".result_dimension()="<<(map).result_dimension()<<", "#set"="<<set); } }

#endif /* ARIADNE_SYSTEM_EXCEPTIONS_H */
