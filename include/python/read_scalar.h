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

#include "numeric/numerical_traits.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

namespace Ariadne {
namespace Python {


template<class X> 
X read_scalar(const boost::python::object& obj);


// Read a scalar variable of type X from a Python object
template<class R> inline
void
read_scalar(R& x, const boost::python::object& obj)
{
  boost::python::extract<std::string> sx(obj);
  boost::python::extract<double> dx(obj);
  boost::python::extract<R> rx(obj);
  if(sx.check()) {
    x=static_cast<R>(sx());
  } else if(dx.check()) {
    x=static_cast<R>(dx());
  } else {
    x=rx();
  }
}
 
template<class R> inline
Numeric::Interval<R>
read_interval(const boost::python::list& pair)
{
  if(boost::python::len(pair)!=2) {
    throw std::runtime_error("Interval must be list of pairs representing intervals");
  }
  R l=read_scalar<R>(pair[0]);
  R u=read_scalar<R>(pair[1]);
  return Numeric::Interval<R>(l,u);
}

template<class R> inline
Numeric::Interval<R>
read_interval(const boost::python::tuple& pair)
{
  if(boost::python::len(pair)!=2) {
    throw std::runtime_error("Interval must be list of pairs representing intervals");
  }
  R l=read_scalar<R>(pair[0]);
  R u=read_scalar<R>(pair[1]);
  return Numeric::Interval<R>(l,u);
}

template<class R> inline
void
read_scalar(Numeric::Interval<R>& x, const boost::python::object& obj)
{
  typedef Numeric::Interval<R> I;
  boost::python::extract<std::string> sx(obj);
  boost::python::extract<boost::python::list> lx(obj);
  boost::python::extract<boost::python::tuple> tx(obj);
  boost::python::extract<double> dx(obj);
  boost::python::extract<R> rx(obj);
  boost::python::extract<I> ix(obj);
  if(sx.check()) {
    x=static_cast<I>(sx());
  } else if(lx.check()) {
    x=read_interval<R>(lx());
  } else if(tx.check()) {
    x=read_interval<R>(tx());
  } else if(dx.check()) {
    x=static_cast<I>(dx());
  } else if(rx.check()) {
    x=static_cast<I>(rx());
  } else {
    x=ix();
  }
}

// Read a scalar variable of type X from a Python object
template<class X>
X
read_scalar(const boost::python::object& obj)
{
  X x;
  read_scalar(x,obj);
  return x;
}


}
}

#endif /* ARIADNE_PYTHON_READ_SCALAR_H */
