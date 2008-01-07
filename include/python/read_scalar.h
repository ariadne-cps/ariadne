/***************************************************************************
 *            python/read_scalar.h
 *
 *  Copyright  2007   Pieter Collins
 *  Pieter.Collins@cwi.nl
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

/*! \file read_scalar.h
 *  Method to read a scalar value from a Python object
 */
 
#ifndef ARIADNE_PYTHON_READ_SCALAR_H
#define ARIADNE_PYTHON_READ_SCALAR_H

#include "numeric/declarations.h"
#include <boost/python.hpp>

namespace Ariadne {
  namespace Python {

    void read_scalar(bool&, const boost::python::object&);
    void read_scalar(int&, const boost::python::object&);
    void read_scalar(uint&, const boost::python::object&);
    void read_scalar(double&, const boost::python::object&);
    void read_scalar(Numeric::Integer&, const boost::python::object&);
    void read_scalar(Numeric::Rational&, const boost::python::object&);
    template<class T> void read_scalar(Numeric::Float<T>&, const boost::python::object&);
    template<class R> void read_scalar(Numeric::Interval<R>&, const boost::python::object&);
    
  }
}

#endif /* ARIADNE_PYTHON_READ_SCALAR_H */
