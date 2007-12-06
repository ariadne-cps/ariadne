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


#include "python/float.h"

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"

#include "python/utilities.h"
using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::Python;

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
Zonotope<R,R> over_approximation_of_minkowski_sum(const Zonotope<R,R>& z1, const Zonotope<R,R>& z2) {
  return over_approximation(minkowski_sum(z1,z2));
}

template<class R> inline
Zonotope<R,R> over_approximation_of_minkowski_difference(const Zonotope<R,R>& z1, const Zonotope<R,R>& z2) {
  return over_approximation(minkowski_difference(z1,z2));
}


template<class R>
void export_zonotope() 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  
  typedef LinearAlgebra::Matrix<R> RMatrix;
  typedef LinearAlgebra::MatrixSlice<R> RMatrixSlice;
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
    .def("centre",(const RPoint&(RZonotope::*)()const)&RZonotope::centre,return_value_policy<copy_const_reference>())
    .def("generators",(const RMatrix&(RZonotope::*)()const)&RZonotope::generators,return_value_policy<copy_const_reference>())
    .def("dimension", &RZonotope::dimension)
    .def("empty", &RZonotope::empty)
    .def("contains", (tribool(RZonotope::*)(const RPoint&)const)&RZonotope::contains)
    .def("bounding_box", &RZonotope::bounding_box)
    .def("subdivide", &RZonotope::subdivide)
    .def("divide", &RZonotope::subdivide)
    .def("__add__", &over_approximation_of_minkowski_sum<R>)
    .def("__sub__", &over_approximation_of_minkowski_difference<R>)
    .def(self_ns::str(self))
  ;
}


template<class R>
void export_interval_zonotope() 
{
  typedef Interval<R> I;
  typedef Matrix<R> RMatrix;
  typedef Matrix<I> IMatrix;
  typedef Point<R> RPoint;
  typedef Point<I> IPoint;
  typedef Rectangle<R> RRectangle;
  typedef Zonotope<R,R> RZonotope;
  typedef Zonotope<I,R> EZonotope;
  typedef Zonotope<I,I> IZonotope;
    
  class_<EZonotope>("ErrorZonotope",init<int>())
    .def(init<EZonotope>())
    .def(init<RZonotope>())
    .def(init<RRectangle>())
    .def("centre",(const IPoint&(EZonotope::*)()const)&EZonotope::centre,return_value_policy<copy_const_reference>())
    .def("generators",(const RMatrix&(EZonotope::*)()const)&EZonotope::generators,return_value_policy<copy_const_reference>())
    .def("dimension", &EZonotope::dimension)
    .def("empty", &EZonotope::empty)
    .def("contains", (tribool(EZonotope::*)(const RPoint&)const)&EZonotope::contains)
    .def("bounding_box", &EZonotope::bounding_box)
    .def("subdivide", &EZonotope::subdivide)
    .def("divide", &EZonotope::divide)
    .def(self_ns::str(self))
  ;


  class_<IZonotope>("FuzzyZonotope",init<int>())
    .def(init<IZonotope>())
    .def(init<EZonotope>())
    .def(init<RZonotope>())
    .def(init<RRectangle>())
    .def("centre",(const IPoint&(IZonotope::*)()const)&IZonotope::centre,return_value_policy<copy_const_reference>())
    .def("generators",(const IMatrix&(IZonotope::*)()const)&IZonotope::generators,return_value_policy<copy_const_reference>())
    .def("dimension", &IZonotope::dimension)
    .def("empty", &IZonotope::empty)
    .def("contains", (tribool(IZonotope::*)(const RPoint&)const)&IZonotope::contains)
    .def("bounding_box", &IZonotope::bounding_box)
    .def("subdivide", &IZonotope::subdivide)
    .def("divide", &IZonotope::divide)
    .def(self_ns::str(self))
  ;

  def("approximation", (Zonotope<R,R>(*)(const Zonotope<I,I>&))(&approximation));
  def("approximation", (Zonotope<R,R>(*)(const Zonotope<I,R>&))(&approximation));
  def("over_approximation", (Zonotope<R,R>(*)(const Zonotope<I,R>&))(&over_approximation));
  def("over_approximation", (Zonotope<I,R>(*)(const Zonotope<I,I>&))(&over_approximation));

/*
  def("interval_over_approximation", (EZonotope(*)(const IZonotope&))(&interval_over_approximation));
  def("orthogonal_over_approximation", (RZonotope(*)(const IZonotope&))(&orthogonal_over_approximation));
*/
}  

template void export_zonotope<FloatPy>();
template void export_interval_zonotope<FloatPy>();
