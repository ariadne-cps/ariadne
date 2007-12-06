/***************************************************************************
 *            python/export_exceptions.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <typeinfo>

#include "base/exceptions.h"

#include "python/utilities.h"
using namespace Ariadne;
using namespace Ariadne::Base;
using namespace std;

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
using namespace boost::python;


template<class E>
void translator(const E& e) {
  PyErr_SetString(PyExc_Exception, e.what());
}

template<>
void translator<NotImplemented>(const NotImplemented& e) {
  PyErr_SetString(PyExc_NotImplementedError, e.what());
}


/* No need to export all exceptions explicitly; just need to put enough information in C++ exceptions. */

void export_exceptions()
{
  register_exception_translator<NotImplemented>(translator<NotImplemented>);
  /*
  register_exception_translator<IncompatibleDimensions>(translator<IncompatibleDimensions>);
  register_exception_translator<InvalidCoordinate>(translator<InvalidCoordinate>);
  register_exception_translator<InvalidVertex>(translator<InvalidVertex>);
  register_exception_translator<InvalidGenerators>(translator<InvalidGenerators>);
  register_exception_translator<IncompatibleGrids>(translator<IncompatibleGrids>);
  register_exception_translator<UnboundedSet>(translator<UnboundedSet>);
  */
}

