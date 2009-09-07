/***************************************************************************
 *            geometry_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include "config.h"

#include "geometry.h"
#include "point.h"
#include "curve.h"
#include "box.h"
#include "zonotope.h"
#include "polytope.h"
#include "polyhedron.h"
#include "taylor_set.h"
#include "grid_set.h"

#include "taylor_function.h"

#include "utilities.h"

#include <boost/python.hpp>

using namespace boost::python;

using namespace Ariadne;

namespace Ariadne {


template<>
struct from_python<Point> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<Point>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        extract<boost::python::tuple> xtup(obj_ptr);
        extract<boost::python::list> xlst(obj_ptr);
        Point pt;
        if(xtup.check()) {
            boost::python::tuple tup=xtup(); pt=Point(len(tup));
            for(int i=0; i!=len(tup); ++i) { pt[i]=extract<double>(tup[i]); }
        } else if(xlst.check()) {
            boost::python::list lst=xlst(); pt=Point(len(lst));
            for(int i=0; i!=len(lst); ++i) { pt[i]=extract<double>(lst[i]); }
        }
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        new (storage) Point(pt);
        data->convertible = storage;
    }
};

template<>
struct from_python<Box> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<Box>()); }
    static void* convertible(PyObject* obj_ptr) { std::cerr<<"Checking from_python<Box>::convertible\n"; if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        boost::python::list lst=extract<boost::python::list>(obj_ptr);
        Box* bx_ptr = new (storage) Box(len(lst));
        for(int i=0; i!=len(lst); ++i) { (*bx_ptr)[i]=extract<Interval>(lst[i]); }
        data->convertible = storage;
    }
};




class OpenSetWrapper
  : public OpenSetInterface, public wrapper< OpenSetInterface >
{
  public:
    OpenSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool covers(const Box& r) const { return this->get_override("covers")(); }
    tribool overlaps(const Box& r) const { return this->get_override("overlaps")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

class ClosedSetWrapper
  : public ClosedSetInterface, public wrapper< ClosedSetInterface >
{
  public:
    ClosedSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool disjoint(const Box& r) const { return this->get_override("disjoint")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


class OvertSetWrapper
  : public OvertSetInterface, public wrapper< OvertSetInterface >
{
  public:
    OvertSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool overlaps(const Box& r) const { return this->get_override("overlaps")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


class CompactSetWrapper
  : public CompactSetInterface, public wrapper< CompactSetInterface >
{
  public:
    CompactSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool disjoint(const Box& r) const { return this->get_override("disjoint")(); }
    tribool inside(const Box& r) const { return this->get_override("inside")(); }
    tribool bounded() const { return this->get_override("bounded")(); }
    Box bounding_box() const { return this->get_override("bounding_box")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

class LocatedSetWrapper
  : public LocatedSetInterface, public wrapper< LocatedSetInterface >
{
  public:
    LocatedSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool covers(const Box& r) const { return this->get_override("covers")(); }
    tribool overlaps(const Box& r) const { return this->get_override("overlaps")(); }
    tribool disjoint(const Box& r) const { return this->get_override("disjoint")(); }
    tribool inside(const Box& r) const { return this->get_override("inside")(); }
    tribool bounded() const { return this->get_override("bounded")(); }
    Box bounding_box() const { return this->get_override("bounding_box")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

}


void export_set_interface() {
    class_<OpenSetInterface, boost::noncopyable> open_set_wrapper_class("OpenSetInterface", no_init);
    open_set_wrapper_class.def("covers",&OpenSetInterface::covers);
    open_set_wrapper_class.def("overlaps",&OpenSetInterface::overlaps);

    class_<CompactSetInterface, boost::noncopyable> compact_set_wrapper_class("CompactSetInterface", no_init);
    compact_set_wrapper_class.def("disjoint",&CompactSetInterface::disjoint);
    compact_set_wrapper_class.def("inside",&CompactSetInterface::inside);
    compact_set_wrapper_class.def("bounding_box",&CompactSetInterface::bounding_box);

    class_<LocatedSetInterface, boost::noncopyable> located_set_wrapper_class("LocatedSetInterface", no_init);
}


void export_point()
{
    class_<Point> point_class("Point",init<Point>());
    point_class.def(init<uint>());
    point_class.def("__getitem__", &__getitem__<Point,int,double>);
    point_class.def(self_ns::str(self));

    from_python<Point>();
    implicitly_convertible<Vector<Float>,Point>();

}

void export_box()
{
    typedef Vector<Interval> IVector;

    class_<Box,bases<CompactSetInterface,OpenSetInterface,Vector<Interval> > > box_class("Box",init<Box>());
    box_class.def(init<uint>());
    box_class.def(init< Vector<Interval> >());
    box_class.def("dimension", (uint(Box::*)()const) &Box::dimension);
    box_class.def("centre", (Point(Box::*)()const) &Box::centre);
    box_class.def("separated", (tribool(Box::*)(const Box&)const) &Box::disjoint);
    box_class.def("overlaps", (tribool(Box::*)(const Box&)const) &Box::overlaps);
    box_class.def("covers", (tribool(Box::*)(const Box&)const) &Box::covers);
    box_class.def("inside", (tribool(Box::*)(const Box&)const) &Box::inside);
    box_class.def("widen", (void(Box::*)()) &Box::widen);
    box_class.def(self_ns::str(self));

    def("split", (std::pair<Box,Box>(*)(const Box&)) &split);
    def("disjoint", (bool(*)(const IVector&,const IVector&)) &disjoint);
    def("subset", (bool(*)(const IVector&,const IVector&)) &subset);

    from_python<Box>();
    to_python< std::pair<Box,Box> >();
    implicitly_convertible<Vector<Interval>,Box>();

}

void export_zonotope()
{
    class_<Zonotope,bases<CompactSetInterface,OpenSetInterface> > zonotope_class("Zonotope",init<Zonotope>());
    zonotope_class.def(init< Point, Matrix<Float>, Vector<Float> >());
    zonotope_class.def(init< Point, Matrix<Float> >());
    zonotope_class.def(init< Box >());
    zonotope_class.def("centre",&Zonotope::centre,return_value_policy<copy_const_reference>());
    zonotope_class.def("generators",&Zonotope::generators,return_value_policy<copy_const_reference>());
    zonotope_class.def("error",&Zonotope::error,return_value_policy<copy_const_reference>());
    zonotope_class.def("contains",&Zonotope::contains);
    zonotope_class.def("__str__",&__cstr__<Zonotope>);

    def("contains", (tribool(*)(const Zonotope&,const Point&)) &contains);
    def("disjoint", (tribool(*)(const Zonotope&,const Box&)) &disjoint);
    def("overlaps", (tribool(*)(const Zonotope&,const Box&)) &overlaps);
    def("disjoint", (tribool(*)(const Zonotope&,const Zonotope&)) &disjoint);
}

void export_polytope()
{
    class_<Polytope,bases<LocatedSetInterface> > polytope_class("Polytope",init<Polytope>());
    polytope_class.def(init<int>());
    polytope_class.def("new_vertex",&Polytope::new_vertex);
    polytope_class.def("__iter__",boost::python::range(&Polytope::vertices_begin,&Polytope::vertices_end));
    polytope_class.def(self_ns::str(self));

}

void export_curve()
{
    to_python< std::pair<const Float,Point> >();

    class_<InterpolatedCurve> interpolated_curve_class("InterpolatedCurve",init<InterpolatedCurve>());
    interpolated_curve_class.def(init<Float,Point>());
    interpolated_curve_class.def("insert", &InterpolatedCurve::insert);
    interpolated_curve_class.def("__iter__",boost::python::range(&InterpolatedCurve::begin,&InterpolatedCurve::end));
    interpolated_curve_class.def(self_ns::str(self));


}

void export_taylor_set()
{
    class_<TaylorSet,bases<CompactSetInterface> > taylor_set_class("TaylorSet",init<TaylorSet>());
    taylor_set_class.def(init<uint>());
    taylor_set_class.def(init<Box>());
    //taylor_set_class.def(init<Zonotope>());
    taylor_set_class.def("bounding_box", &TaylorSet::bounding_box);
    taylor_set_class.def("range", &TaylorSet::bounding_box);
    taylor_set_class.def(self_ns::str(self));

    def("split", (std::pair<TaylorSet,TaylorSet>(TaylorSet::*)()const) &TaylorSet::split);
    def("outer_approximation", (GridTreeSet(*)(const TaylorSet&,const Grid&,uint)) &outer_approximation);
    def("adjoin_outer_approximation", (void(*)(GridTreeSet&,const TaylorSet&,uint)) &adjoin_outer_approximation);
    def("zonotope", (Zonotope(*)(const TaylorSet&)) &zonotope);

    def("apply",(TaylorModel(*)(const ScalarFunctionInterface&,const TaylorSet&)) &apply);
    def("apply",(TaylorModel(*)(const TaylorExpression&,const TaylorSet&)) &apply);
    def("apply",(TaylorSet(*)(const FunctionInterface&,const TaylorSet&)) &apply);
    def("apply",(TaylorSet(*)(const TaylorFunction&,const TaylorSet&)) &apply);

    implicitly_convertible<Box,TaylorSet>();
}


void geometry_submodule() {
    export_set_interface();
    export_point();
    export_box();
    export_zonotope();
    export_polytope();
    export_taylor_set();
    export_curve();
}

