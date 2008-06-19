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
#include "python/read_array.h"

#include "function/affine_model.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"

#include "python/utilities.h"
using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
std::string
__str__(const Zonotope<R>& z)
{
  std::stringstream ss;
  ss << z;
  return ss.str();
}

template<class R>
std::string
__repr__(const Zonotope<R>& z)
{
  std::stringstream ss;
  ss << "Zonotope(" << z << ")";
  return ss.str();
}

template<class R>  
Zonotope<R>* 
make_zonotope(const boost::python::object& obj) 
{
  // See "Extracting C++ objects" in the Boost Python tutorial
  boost::python::list elements=extract<list>(obj);
  int ng=boost::python::len(elements)-1;
  array<R> ary;
  Point<R> c;
  read_tuple_array(const_cast<array<R>&>(c.data()),elements[0]);
  int d=c.dimension();
  Matrix<R> G(d,ng);
  Vector<R> g(d);
  for(int j=0; j!=ng; ++j) {
    read_array(g.data(),elements[j+1]);
    G.column(j)=g;
  }
  return new Zonotope<R>(c,G);
}

template<class R>
void export_zonotope() 
{
  typedef Interval<R> I;
  typedef Matrix<R> RMatrix;
  typedef Matrix<I> IMatrix;
  typedef Point<R> RPoint;
  typedef Point<I> IPoint;
  typedef Box<R> RBox;
  typedef Zonotope<R> RZonotope;
  typedef ListSet< Zonotope<R> > RZonotopeListSet;
  typedef AffineModel<R> RAffineModel;

  class_<RZonotope> zonotope_class("Zonotope",init<int>());
  zonotope_class.def("__init__", make_constructor(&make_zonotope<R>) );
  zonotope_class.def(init<RPoint,RMatrix>());
  zonotope_class.def(init<RZonotope>());
  zonotope_class.def(init<RBox>());
  zonotope_class.def("centre",&RZonotope::centre,return_value_policy<copy_const_reference>());
  zonotope_class.def("generators",&RZonotope::generators,return_value_policy<copy_const_reference>());
  zonotope_class.def("error",&RZonotope::error,return_value_policy<copy_const_reference>());
  zonotope_class.def("dimension", &RZonotope::dimension);
  zonotope_class.def("contains", &RZonotope::contains);
  zonotope_class.def("bounding_box", &RZonotope::bounding_box);
  zonotope_class.def("radius", &RZonotope::radius);
  zonotope_class.def("__str__", (std::string(*)(const RZonotope&))&__str__);
  zonotope_class.def("__repr__", (std::string(*)(const RZonotope&))&__repr__);
 
  def("radius", (R(*)(const RZonotope&))(&radius));
  def("bounding_box", (RBox(*)(const RZonotope&))(&bounding_box));

  def("contains", (tribool(*)(const RZonotope&,const RPoint&))(&contains));
  def("disjoint", (tribool(*)(const RZonotope&,const RBox&))(&disjoint));
  def("subset", (tribool(*)(const RZonotope&,const RBox&))(&subset));
  def("superset", (tribool(*)(const RZonotope&,const RBox&))(&superset));
  def("subset", (tribool(*)(const RBox&,const RZonotope&))(&subset));

  def("approximation", (RZonotope(*)(const RZonotope&))(&approximation));
  def("over_approximation", (RZonotope(*)(const RZonotope&))(&over_approximation));
  def("error_free_over_approximation", (RZonotope(*)(const RZonotope&))(&error_free_over_approximation));
  def("orthogonal_over_approximation", (RZonotope(*)(const RZonotope&))(&orthogonal_over_approximation));
  def("nonsingular_over_approximation", (RZonotope(*)(const RZonotope&))(&nonsingular_over_approximation));
  def("cascade_over_approximation",(RZonotope(*)(const RZonotope&,size_type)) &cascade_over_approximation);

  def("split", (RZonotopeListSet(*)(const RZonotope&))(&::split));

  def("apply", (RZonotope(*)(const RAffineModel&, const RZonotope&))(&apply));
}  


template<>
void export_zonotope<Rational>() 
{
  typedef Rational Q;
  
  typedef Matrix<Q> QMatrix;
  typedef Point<Q> QPoint;
  typedef Box<Q> QBox;
  typedef Zonotope<Q> QZonotope;
  
  class_<QZonotope> zonotope_class("Zonotope",init<int>());
  zonotope_class.def(init<QPoint,QMatrix>());
  zonotope_class.def(init<QZonotope>());
  zonotope_class.def(init<QBox>());
  zonotope_class.def("centre",&QZonotope::centre,return_value_policy<copy_const_reference>());
  zonotope_class.def("generators",&QZonotope::generators,return_value_policy<copy_const_reference>());
  zonotope_class.def("dimension", &QZonotope::dimension);
  zonotope_class.def("contains", (tribool(QZonotope::*)(const QPoint&)const)&QZonotope::contains);
  zonotope_class.def("bounding_box", &QZonotope::bounding_box);
  zonotope_class.def("__str__", (std::string(*)(const QZonotope&)) &__str__);
  zonotope_class.def("__repr__", (std::string(*)(const QZonotope&)) &__repr__);

  def("contains", (tribool(*)(const QZonotope&,const QPoint&))(&contains));
  def("disjoint", (tribool(*)(const QZonotope&,const QBox&))(&disjoint));
  def("subset", (tribool(*)(const QZonotope&,const QBox&))(&subset));
  def("superset", (tribool(*)(const QZonotope&,const QBox&))(&superset));
  def("subset", (tribool(*)(const QBox&,const QZonotope&))(&subset));

}


template void export_zonotope<FloatPy>();
