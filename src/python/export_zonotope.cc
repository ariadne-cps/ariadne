/***************************************************************************
 *            python/export_zonotope.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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


#include "real_typedef.h"

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"

#include "python/python_utilities.h"
using namespace Ariadne;
using namespace Ariadne::Geometry;

#include <boost/python.hpp>
using namespace boost::python;

template<class BS1, class BS2>
inline
BS1
touching_intersection(const BS1 &a, 
                      const BS2 &b) 
{
  if (disjoint(a,b)) {
    return BS1(a.dimension());
  }
  return a;
}

template<class R> inline
Zonotope<R> over_approximation_of_minkowski_sum(const Zonotope<R>& z1, const Zonotope<R>& z2) {
  return over_approximation(minkowski_sum(z1,z2));
}

template<class R> inline
Zonotope<R> over_approximation_of_minkowski_difference(const Zonotope<R>& z1, const Zonotope<R>& z2) {
  return over_approximation(minkowski_difference(z1,z2));
}


template<class R>
void export_zonotope() 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  
  typedef LinearAlgebra::Matrix<R> RMatrix;
  typedef Point<R> RPoint;
  typedef Zonotope<R> RZonotope;
  typedef Parallelotope<R> RParallelotope;
  typedef Rectangle<R> RRectangle;
  
  typedef Zonotope<F> FZonotope;
  
  def("disjoint", (tribool(*)(const RZonotope&,const RZonotope&))(&disjoint));
  def("disjoint", (tribool(*)(const RZonotope&,const RRectangle&))(&disjoint));
  def("disjoint", (tribool(*)(const RRectangle&,const RZonotope&))(&disjoint));
  def("subset", (tribool(*)(const RZonotope&,const RZonotope&))(&subset));
  def("subset", (tribool(*)(const RZonotope&,const RRectangle&))(&subset));
  def("subset", (tribool(*)(const RRectangle&,const RZonotope&))(&subset));

  def("touching_intersection", (RZonotope(*)(const RZonotope&,const RZonotope&))(&touching_intersection));
  def("touching_intersection", (RZonotope(*)(const RZonotope&,const RRectangle&))(&touching_intersection));

  def("minkowski_sum", &over_approximation_of_minkowski_sum<R>);
  def("minkowski_difference", &over_approximation_of_minkowski_difference<R>);

  class_<RZonotope>("Zonotope",init<int>())
    .def(init<RPoint,RMatrix>())
    .def(init<RZonotope>())
    .def(init<RParallelotope>())
    .def(init<RRectangle>())
    .def(init<std::string>())
    .def("centre",&RZonotope::centre)
    .def("generators",&RZonotope::generators)
    .def("dimension", &RZonotope::dimension)
    .def("empty", &RZonotope::empty)
    .def("contains", &RZonotope::contains)
    .def("bounding_box", &RZonotope::bounding_box)
    .def("subdivide", &RZonotope::subdivide)
    .def("__add__", &over_approximation_of_minkowski_sum<R>)
    .def("__sub__", &over_approximation_of_minkowski_difference<R>)
    .def(self_ns::str(self))
  ;
}


template<class R>
void export_interval_zonotope() 
{
  typedef Interval<R> I;
  typedef Rectangle<R> RRectangle;
  typedef Zonotope<R> RZonotope;
  typedef Zonotope<I> IZonotope;
    
  class_<IZonotope>("IntervalZonotope",init<int>())
    .def(init<RZonotope>())
    .def(init<RRectangle>())
    .def(init<std::string>())
    .def("centre",&IZonotope::centre)
    .def("generators",&IZonotope::generators)
    .def("dimension", &IZonotope::dimension)
    .def("empty", &IZonotope::empty)
    .def("contains", &IZonotope::contains)
    .def("bounding_box", &IZonotope::bounding_box)
    .def("subdivide", &IZonotope::subdivide)
    .def(self_ns::str(self))
  ;
    
  def("interval_over_approximation", (IZonotope(*)(const IZonotope&))(&interval_over_approximation));
  def("over_approximation", (RZonotope(*)(const IZonotope&))(&over_approximation));
  def("approximation", (RZonotope(*)(const IZonotope&))(&approximation));

}

template void export_zonotope<Real>();
template void export_interval_zonotope<Real>();
