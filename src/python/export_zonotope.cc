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


template<class R>
void export_zonotope() 
{
  typedef Interval<R> I;
  typedef Matrix<R> RMatrix;
  typedef Matrix<I> IMatrix;
  typedef Point<R> RPoint;
  typedef Point<I> IPoint;
  typedef Rectangle<R> RRectangle;
  typedef Zonotope<R,ExactTag> RZonotope;
  typedef Zonotope<R,UniformErrorTag> EZonotope;
  typedef Zonotope<R,IntervalTag> IZonotope;
    

  class_<RZonotope> rzonotope_class("Zonotope",init<int>());
  rzonotope_class.def(init<RPoint,RMatrix>());
  rzonotope_class.def(init<RZonotope>());
  rzonotope_class.def(init<RRectangle>());
  rzonotope_class.def("centre",(const RPoint&(RZonotope::*)()const)&RZonotope::centre,return_value_policy<copy_const_reference>());
  rzonotope_class.def("generators",(const RMatrix&(RZonotope::*)()const)&RZonotope::generators,return_value_policy<copy_const_reference>());
  rzonotope_class.def("dimension", &RZonotope::dimension);
  rzonotope_class.def("empty", &RZonotope::empty);
  rzonotope_class.def("contains", (tribool(*)(const RZonotope&,const RPoint&))&Geometry::contains);
  rzonotope_class.def("bounding_box", &RZonotope::bounding_box);
  rzonotope_class.def(self_ns::str(self));
 
  def("contains", (tribool(*)(const RZonotope&,const RPoint&))(&contains));
  def("disjoint", (tribool(*)(const RZonotope&,const RRectangle&))(&disjoint));
  def("subset", (tribool(*)(const RZonotope&,const RRectangle&))(&subset));
  def("superset", (tribool(*)(const RZonotope&,const RRectangle&))(&superset));
  def("subset", (tribool(*)(const RRectangle&,const RZonotope&))(&subset));


  class_<EZonotope> ezonotope_class("ErrorZonotope",init<int>());
  ezonotope_class.def(init<EZonotope>());
  ezonotope_class.def(init<RZonotope>());
  ezonotope_class.def(init<RRectangle>());
  ezonotope_class.def("centre",(const IPoint&(EZonotope::*)()const)&EZonotope::centre,return_value_policy<copy_const_reference>());
  ezonotope_class.def("generators",(const RMatrix&(EZonotope::*)()const)&EZonotope::generators,return_value_policy<copy_const_reference>());
  ezonotope_class.def("dimension", &EZonotope::dimension);
  ezonotope_class.def("empty", &EZonotope::empty);
  rzonotope_class.def("contains", (tribool(*)(const EZonotope&,const RPoint&))&Geometry::contains);
  ezonotope_class.def("bounding_box", &EZonotope::bounding_box);
  ezonotope_class.def("subdivide", (ListSet<EZonotope>(*)(const EZonotope&)) &Geometry::subdivide);
  ezonotope_class.def(self_ns::str(self));

  def("contains", (tribool(*)(const EZonotope&,const RPoint&))(&contains));


  class_<IZonotope> izonotope_class("IntervalZonotope",init<int>());
  izonotope_class.def("centre",(const IPoint&(IZonotope::*)()const)&IZonotope::centre,return_value_policy<copy_const_reference>());
  izonotope_class.def("generators",(const IMatrix&(IZonotope::*)()const)&IZonotope::generators,return_value_policy<copy_const_reference>());
  izonotope_class.def("dimension", &IZonotope::dimension);
  izonotope_class.def(self_ns::str(self));


  def("approximation", (RZonotope(*)(const IZonotope&))(&::approximation));
  def("approximation", (RZonotope(*)(const EZonotope&))(&::approximation));
  def("over_approximation", (RZonotope(*)(const EZonotope&))(&::approximation));
  def("over_approximation", (EZonotope(*)(const IZonotope&))(&::over_approximation));


}  


template<>
void export_zonotope<Rational>() 
{
  typedef Rational Q;
  
  typedef Matrix<Q> QMatrix;
  typedef Point<Q> QPoint;
  typedef Rectangle<Q> QRectangle;
  typedef Zonotope<Q> QZonotope;
  
  class_<QZonotope> zonotope_class("Zonotope",init<int>());
  zonotope_class.def(init<QPoint,QMatrix>());
  zonotope_class.def(init<QZonotope>());
  zonotope_class.def(init<QRectangle>());
  zonotope_class.def("centre",&QZonotope::centre,return_value_policy<copy_const_reference>());
  zonotope_class.def("generators",&QZonotope::generators,return_value_policy<copy_const_reference>());
  zonotope_class.def("dimension", &QZonotope::dimension);
  zonotope_class.def("empty", &QZonotope::empty);
  zonotope_class.def("contains", (tribool(QZonotope::*)(const QPoint&)const)&QZonotope::contains);
  zonotope_class.def("bounding_box", &QZonotope::bounding_box);
  zonotope_class.def("subdivide", (ListSet<QZonotope>(*)(const QZonotope&))&Geometry::subdivide);
  zonotope_class.def(self_ns::str(self));

  def("contains", (tribool(*)(const QZonotope&,const QPoint&))(&contains));
  def("disjoint", (tribool(*)(const QZonotope&,const QRectangle&))(&disjoint));
  def("subset", (tribool(*)(const QZonotope&,const QRectangle&))(&subset));
  def("superset", (tribool(*)(const QZonotope&,const QRectangle&))(&superset));
  def("subset", (tribool(*)(const QRectangle&,const QZonotope&))(&subset));

}


template void export_zonotope<FloatPy>();
