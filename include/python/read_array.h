/***************************************************************************
 *            python/read_array.h
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

/*! \file read_array.h
 *  Method to read a array value from a Python object
 */
 
#ifndef ARIADNE_PYTHON_READ_ARRAY_H
#define ARIADNE_PYTHON_READ_ARRAY_H

#include "base/array.h"
#include "numeric/traits.h"
#include "python/read_scalar.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

namespace Ariadne {
namespace Python {


template<class X> 
array<X>*
read_array(const boost::python::object& obj);

template<class X> 
void
read_array(array<X>&, const boost::python::object& obj);


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
    ary[i]=read_scalar<X>(elements[i]);
  }
}

template<class X> 
void
read_tuple_array(array<X>& ary, const boost::python::object& obj)
{
  // See "Extracting C++ objects" in the Boost Python tutorial
  boost::python::tuple elements=boost::python::extract<boost::python::tuple>(obj);
  int n=boost::python::len(elements);
  ary.resize(n);
  for(int i=0; i!=n; ++i) {
    ary[i]=read_scalar<X>(elements[i]);
  }
}

template<class X> inline
array<X>*
read_array(const boost::python::object& obj)
{
  boost::python::list elements=boost::python::extract<boost::python::list>(obj);
  int n=boost::python::len(elements);
  array<X>& ary=*new array<X>(n);
  for(int i=0; i!=n; ++i) {
    ary[i]=read_scalar<X>(elements[i]);
  }
  return &ary;
}


}
}

#endif /* ARIADNE_PYTHON_READ_ARRAY_H */
