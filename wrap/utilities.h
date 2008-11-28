/***************************************************************************
 *            utilities.h
 *
 *  Copyright  2005-8  Alberto Casagrande, Pieter Collins
 *
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

/*! \file utilities.h
 *  Commonly used inline methods for the Python interface.
 */
 
#ifndef ARIADNE_PYTHON_UTILITIES_H
#define ARIADNE_PYTHON_UTILITIES_H

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#include <string>
#include <sstream>
#include <iostream>

#include "array.h"
#include "numeric.h"

namespace Ariadne {

template<class T> std::ostream& repr(std::ostream& os, const T& t) {
    return os << t;
}


template<class T> std::string __repr__(const T& t) {
    std::stringstream ss;
    repr(ss,t);
    return ss.str();
}


void read_scalar(bool&, const boost::python::object&);
void read_scalar(int&, const boost::python::object&);
void read_scalar(long int&, const boost::python::object&);
void read_scalar(unsigned int&, const boost::python::object&);
void read_scalar(unsigned long int&, const boost::python::object&);
void read_scalar(double&, const boost::python::object&);
void read_scalar(Float&, const boost::python::object&);
void read_scalar(Interval&, const boost::python::object&);

#ifdef HAVE_GMPXX_H
//void read_scalar(Integer&, const boost::python::object&);
void read_scalar(Rational&, const boost::python::object&);
#endif


// Read a array variable of type X from a Python object
template<class X> 
void
read_array(array<X>& ary, const boost::python::object& obj)
{
    // See "Extracting C++ objects" in the Boost Python tutorial
    boost::python::list elements=boost::python::extract<boost::python::list>(obj);
    int n=boost::python::len(elements);
    ary.resize(n);
    for(int i=0; i!=n; ++i) {
        read_scalar(ary[i], elements[i]);
    }
}

}

#endif /* ARIADNE_PYTHON_UTILITIES_H */
