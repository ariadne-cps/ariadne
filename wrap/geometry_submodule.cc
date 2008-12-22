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

#include "point.h"
#include "box.h"
#include "zonotope.h"
#include "polytope.h"
#include "polyhedron.h"
#include "taylor_set.h"
#include "grid_set.h"

#include "utilities.h"

#include <boost/python.hpp>

using namespace boost::python;

using namespace Ariadne;

namespace Ariadne {

template<class X> void read(Matrix<X>& A, const boost::python::object& obj);

void
read(Point& pt, const boost::python::object& obj) 
{
  array<Float> ary;
  read_tuple_array(ary,obj);
  pt=Point(Vector<Float>(ary.size(),ary.begin()));
}


void
read(Box& bx, const boost::python::object& obj) 
{
  array<Interval> ary;
  read_list_array(ary,obj);
  bx=Box(Vector<Interval>(ary.size(),ary.begin()));
}

void
read(Zonotope& z, const boost::python::object& obj) 
{
    boost::python::tuple tup=boost::python::extract<boost::python::tuple>(obj);
    Point c;
    Matrix<Float> GT;
    assert(len(tup)==2);
    read(c,tup[0]);
    read(GT,tup[1]);
    z=Zonotope(c,transpose(GT));
}


void read(TaylorVariable& tv, const boost::python::object& obj);


void
read(TaylorSet& ts, const boost::python::object& obj) 
{
    boost::python::list lst=extract<boost::python::list>(obj);
    ts=TaylorSet(len(lst));
    for(uint i=0; i!=len(lst); ++i) {
        TaylorVariable tv;
        read(tv,lst[i]);
        ts[i]=tv;
    }
}
    
template<class SET> boost::python::tuple split(const SET& s, uint i) {
    std::pair<SET,SET> res=s.split(i);
    return boost::python::make_tuple(res.first,res.second);
}


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
    class_<OpenSetWrapper, boost::noncopyable> open_set_wrapper_class("OpenSetInterface");
    class_<CompactSetWrapper, boost::noncopyable> compact_set_wrapper_class("CompactSetInterface");
    class_<LocatedSetWrapper, boost::noncopyable> located_set_wrapper_class("LocatedSetInterface");
}


void export_point() 
{
    class_<Point> point_class("Point",no_init);
    point_class.def("__init__", make_constructor(&make<Point>) );
    point_class.def("__str__",&__str__<Point>);
}

void export_box() 
{
    class_<Box,bases<CompactSetInterface,OpenSetInterface> > box_class("Box",init<>());
    box_class.def("__init__", make_constructor(&make<Box>) );
    box_class.def("__str__",&__str__<Box>);
}

void export_zonotope() 
{
    class_<Zonotope,bases<CompactSetInterface,OpenSetInterface> > zonotope_class("Zonotope",no_init);
    zonotope_class.def("__init__", make_constructor(&make<Zonotope>) );
    zonotope_class.def("centre",&Zonotope::centre,return_value_policy<copy_const_reference>());
    zonotope_class.def("generators",&Zonotope::generators,return_value_policy<copy_const_reference>());
    zonotope_class.def("error",&Zonotope::error,return_value_policy<copy_const_reference>());
    zonotope_class.def("__str__",&__str__<Zonotope>);

    def("contains", (tribool(*)(const Zonotope&,const Point&)) &contains);
    def("disjoint", (tribool(*)(const Zonotope&,const Box&)) &disjoint);
    def("overlaps", (tribool(*)(const Zonotope&,const Box&)) &overlaps);
    def("disjoint", (tribool(*)(const Zonotope&,const Zonotope&)) &disjoint);
}

void export_taylor_set() 
{
    class_<TaylorSet,bases<CompactSetInterface> > taylor_set_class("TaylorSet",no_init);
    taylor_set_class.def("__init__", make_constructor(&make<TaylorSet>) );
    taylor_set_class.def("split", &split<TaylorSet>);
    taylor_set_class.def("outer_approximation", (GridTreeSet(*)(const TaylorSet&,uint)) &outer_approximation);
    taylor_set_class.def("__str__",&__str__<TaylorSet>);
    
    def("zonotope", (Zonotope(*)(const TaylorSet&)) &zonotope);
}


void geometry_submodule() {
    export_set_interface();
    export_point();
    export_box();
    export_zonotope();
    export_taylor_set();
}

