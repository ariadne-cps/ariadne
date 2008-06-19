/***************************************************************************
 *            python/read_matrix.h
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

/*! \file read_matrix.h
 *  Method to read a matrix value from a Python object
 */
 
#ifndef ARIADNE_PYTHON_READ_MATRIX_H
#define ARIADNE_PYTHON_READ_MATRIX_H

#include "numeric/traits.h"
#include "linear_algebra/matrix.h"
#include "python/read_scalar.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

namespace Ariadne {
namespace Python {


template<class X>  
void
read_matrix(Matrix<X>& A, const boost::python::object& obj) 
{
  // See "Extracting C++ objects" in the Boost Python tutorial
  typedef typename traits<X>::number_type R;
  boost::python::list elements=boost::python::extract<boost::python::list>(obj);
  int m=boost::python::len(elements);
  boost::python::list row=boost::python::extract<boost::python::list>(elements[0]);
  int n=boost::python::len(row);
  A.resize(m,n);
  for(int i=0; i!=m; ++i) {
    row=boost::python::extract<boost::python::list>(elements[i]);
    if(boost::python::len(row)!=n) {
      throw std::runtime_error("Matrix with rows of different sizes");
    }
    for(int j=0; j!=n; ++j) {
      read_scalar(A(i,j),row[j]);
    }
  }
}


}
}

#endif /* ARIADNE_PYTHON_READ_MATRIX_H */
