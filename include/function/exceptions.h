/***************************************************************************
 *            function/exceptions.h
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
 
/*! \file function/exceptions.h
 *  \brief Exceptions, error handling and assertions for the Function module.
 */

#ifndef ARIADNE_FUNCTION_EXCEPTIONS_H
#define ARIADNE_FUNCTION_EXCEPTIONS_H

#include <stdexcept>
#include <iosfwd>

#include "linear_algebra/exceptions.h"

namespace Ariadne {
  namespace Function {
    using LinearAlgebra::IncompatibleSizes;
    struct InvalidParameters : public std::runtime_error { InvalidParameters(const std::string& what) : std::runtime_error(what) { } };
  }
}

#define ARIADNE_CHECK_ARGUMENT_SIZE(function,argument,func) \
  { if((function).argument_size()!=(argument).size()) { ARIADNE_THROW(IncompatibleSizes,func,#function".argument_size()="<<(function).argument_size()<<", "#argument"="<<argument); } }


#endif /* ARIADNE_FUNCTION_EXCEPTIONS_H */
