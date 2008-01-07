/***************************************************************************
 *            python/read_box.h
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

/*! \file read_box.h
 *  Method to read a box from a Python object
 */
 
#ifndef ARIADNE_PYTHON_READ_BOX_H
#define ARIADNE_PYTHON_READ_BOX_H

#include "numeric/traits.h"
#include "geometry/box.h" 
#include "read_array.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>

namespace Ariadne {
namespace Python {

template<class R>  
void
read_box(Geometry::Box<R>& bx, const boost::python::object& obj) 
{
  // See "Extracting C++ objects" in the Boost Python tutorial
  array< Numeric::Interval<R> > ary;
  read_array(ary,obj);
  bx=Geometry::Box<R>(ary);
}

template<class R>  
void
old_read_box(Geometry::Box<R>& bx, const boost::python::object& obj) 
{
  boost::python::list elements=boost::python::extract<const boost::python::list>(obj);
  int d=boost::python::len(elements);
  bx.resize(d);
  R l,u;
  for(int i=0; i!=d; ++i) {
    boost::python::extract<const boost::python::list> extract_list(elements[i]);
    if(extract_list.check()) {
      boost::python::list pair=extract_list();
      if(boost::python::len(pair)!=2) {
        throw std::runtime_error("Box must be list of pairs representing intervals");
      }
      read_scalar(l,pair[0]);
      read_scalar(u,pair[1]);
      bx.set_lower_bound(i,l);
      bx.set_upper_bound(i,u);
    } else {
      boost::python::extract<std::string> extract_string(elements[i]);
      if(extract_string.check()) {
        Numeric::Interval<R> interval(extract_string());
        bx.set_lower_bound(i,interval.lower());
        bx.set_upper_bound(i,interval.upper());
      } else {
        boost::python::extract< Numeric::Interval<R> > extract_interval(elements[i]);
        Numeric::Interval<R> interval=extract_interval();
        bx.set_lower_bound(i,interval.lower());
        bx.set_upper_bound(i,interval.upper());
      } 
    }
  }
}



}
}

#endif /* ARIADNE_PYTHON_READ_BOX_H */
