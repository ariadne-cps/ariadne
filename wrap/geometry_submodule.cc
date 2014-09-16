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

#include "utilities.h"

#include "config.h"

#include <boost/python.hpp>

#include "geometry.h"
#include "geometry2d.h"
#include "point.h"
#include "curve.h"
#include "interval.h"
#include "box.h"
#include "zonotope.h"
#include "polytope.h"
#include "polyhedron.h"
#include "grid_set.h"
#include "function_set.h"
#include "affine_set.h"

#include "discrete_location.h"
#include "hybrid_set.h"



using namespace boost::python;
using namespace Ariadne;

namespace Ariadne {


template<>
struct from_python_dict<Interval> {
    from_python_dict() { converter::registry::push_back(&convertible,&construct,type_id<Interval>()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyDict_Check(obj_ptr) || len(boost::python::extract<boost::python::dict>(obj_ptr))!=1) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::dict dct = boost::python::extract<boost::python::dict>(obj_ptr);
        boost::python::list lst=dct.items();
        assert(boost::python::len(lst)==1);
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        new (storage) Interval(boost::python::extract<Float>(lst[0][0]),boost::python::extract<Float>(lst[0][1]));
        data->convertible = storage;
    }
};


template<>
struct from_python_list<Interval> {
    from_python_list() { converter::registry::push_back(&convertible,&construct,type_id<Interval>()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) || len(boost::python::extract<boost::python::list>(obj_ptr))!=2) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        assert(boost::python::len(lst)==2);
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        new (storage) Interval(boost::python::extract<Float>(lst[0]),boost::python::extract<Float>(lst[1]));
        data->convertible = storage;
    }
};

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
            for(int i=0; i!=len(tup); ++i) { pt[i]=Float(extract<double>(tup[i])); }
        } else if(xlst.check()) {
            boost::python::list lst=xlst(); pt=Point(len(lst));
            for(int i=0; i!=len(lst); ++i) { pt[i]=Float(extract<double>(lst[i])); }
        }
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        new (storage) Point(pt);
        data->convertible = storage;
    }
};

template<>
struct from_python<Box> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<Box>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        boost::python::list lst=extract<boost::python::list>(obj_ptr);
        Box* bx_ptr = new (storage) Box(len(lst));
        for(int i=0; i!=len(lst); ++i) { (*bx_ptr)[i]=extract<Interval>(lst[i]); }
        data->convertible = storage;
    }
};


template<class ES>
struct to_python< ListSet<ES> > {
    to_python() { boost::python::to_python_converter< ListSet<ES>, to_python< ListSet<ES> > >(); }

    static PyObject* convert(const ListSet<ES>& ls) {
        boost::python::list result;
        for(typename ListSet<ES>::const_iterator iter=ls.begin(); iter!=ls.end(); ++iter) {
            result.append(boost::python::object(*iter));
        }
        return boost::python::incref(boost::python::list(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

template<class ES>
struct to_python< ListSet< HybridBasicSet<ES> > > {
    typedef ListSet< HybridBasicSet<ES> > SetType;
    to_python() { boost::python::to_python_converter< SetType, to_python<SetType> >(); }

    static PyObject* convert(const SetType& hls) {
        boost::python::dict result;
        for(typename SetType::locations_const_iterator iter=hls.locations_begin(); iter!=hls.locations_end(); ++iter) {
            result[iter->first]=iter->second;
        }
        return boost::python::incref(boost::python::dict(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyDict_Type; }
};

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<ValidatedFloat>& x);
std::ostream& operator<<(std::ostream& os, const PythonRepresentation<Interval>& x) {
    return os << PythonRepresentation<ValidatedFloat>(ValidatedFloat(x.reference()));
}



class OpenSetWrapper
  : public virtual OpenSetInterface, public wrapper< OpenSetInterface >
{
  public:
    OpenSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool covers(const Box& r) const { return this->get_override("covers")(); }
    tribool overlaps(const Box& r) const { return this->get_override("overlaps")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

class ClosedSetWrapper
  : public virtual ClosedSetInterface, public wrapper< ClosedSetInterface >
{
  public:
    ClosedSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool separated(const Box& r) const { return this->get_override("separated")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


class OvertSetWrapper
  : public virtual OvertSetInterface, public wrapper< OvertSetInterface >
{
  public:
    OvertSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool overlaps(const Box& r) const { return this->get_override("overlaps")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


class CompactSetWrapper
  : public virtual CompactSetInterface, public wrapper< CompactSetInterface >
{
  public:
    CompactSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool separated(const Box& r) const { return this->get_override("separated")(); }
    tribool inside(const Box& r) const { return this->get_override("inside")(); }
    tribool bounded() const { return this->get_override("bounded")(); }
    Box bounding_box() const { return this->get_override("bounding_box")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

class RegularSetWrapper
  : public virtual LocatedSetInterface, public wrapper< RegularSetWrapper >
{
  public:
    RegularSetWrapper* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool overlaps(const Box& r) const { return this->get_override("overlaps")(); }
    tribool covers(const Box& r) const { return this->get_override("covers")(); }
    tribool separated(const Box& r) const { return this->get_override("separated")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

class LocatedSetWrapper
  : public virtual LocatedSetInterface, public wrapper< LocatedSetInterface >
{
  public:
    LocatedSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool overlaps(const Box& r) const { return this->get_override("overlaps")(); }
    tribool separated(const Box& r) const { return this->get_override("separated")(); }
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

    class_<ClosedSetInterface, boost::noncopyable> closed_set_wrapper_class("ClosedSetInterface", no_init);
    closed_set_wrapper_class.def("separated",&ClosedSetInterface::separated);

    class_<OvertSetInterface, boost::noncopyable> overt_set_wrapper_class("OvertSetInterface", no_init);
    overt_set_wrapper_class.def("overlaps",&OvertSetInterface::overlaps);

    class_<CompactSetInterface, boost::noncopyable> compact_set_wrapper_class("CompactSetInterface", no_init);
    compact_set_wrapper_class.def("separated",&CompactSetInterface::separated);
    compact_set_wrapper_class.def("inside",&CompactSetInterface::inside);
    compact_set_wrapper_class.def("bounding_box",&CompactSetInterface::bounding_box);

    class_<LocatedSetInterface, boost::noncopyable> located_set_wrapper_class("LocatedSetInterface", no_init);
    class_<RegularSetInterface, boost::noncopyable> regular_set_wrapper_class("RegularSetInterface", no_init);

    class_<DrawableInterface,boost::noncopyable>("DrawableInterface",no_init);

}


void export_point()
{
    class_<Point,bases<DrawableInterface>> point_class("Point",init<Point>());
    point_class.def(init<uint>());
    point_class.def("__getitem__", &__getitem__<Point,int,Float>);
    point_class.def(self_ns::str(self));

    from_python<Point>();
    implicitly_convertible<Vector<Float>,Point>();

}


void export_interval()
{
    using boost::python::class_;
    using boost::python::init;
    using boost::python::self;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::def;

    class_< Interval > interval_class("Interval");
    interval_class.def(init<double>());
    interval_class.def(init<double,double>());
    interval_class.def(init<Float,Float>());
    interval_class.def(init<Real,Real>());
    interval_class.def(init<Decimal>());
    interval_class.def(init<Dyadic>());
    interval_class.def(init<ValidatedFloat>());
    interval_class.def(init<Float>());
#ifdef HAVE_GMPXX_H
    interval_class.def(init<Rational>());
    interval_class.def(init<Rational,Rational>());
#endif

    interval_class.def(self == self);
    interval_class.def(self != self);
    interval_class.def("lower", &Interval::lower, return_value_policy<copy_const_reference>());
    interval_class.def("upper", &Interval::upper, return_value_policy<copy_const_reference>());
    interval_class.def("midpoint", &Interval::midpoint);
    interval_class.def("radius", &Interval::radius);
    interval_class.def("width", &Interval::width);
    interval_class.def("contains", (bool(*)(Interval,Float)) &contains);
    interval_class.def("empty", (bool(Interval::*)()const) &Interval::empty);
    interval_class.def(boost::python::self_ns::str(self));

    from_python_dict<Interval>();
    from_python_list<Interval>();
    //from_python_str<Interval>();

    def("midpoint", &Interval::midpoint);
    def("radius", &Interval::radius);
    def("width", &Interval::width);

    def("disjoint", (bool(*)(Interval,Interval)) &disjoint);
    def("subset", (bool(*)(Interval,Interval)) &subset);
    def("intersection", (Interval(*)(Interval,Interval)) &intersection);

    def("hull", (Interval(*)(Interval,Interval)) &hull);
}

void export_box()
{
    typedef Vector<Interval> IVector;

    class_<Box,bases<CompactSetInterface,OpenSetInterface,Vector<Interval>,DrawableInterface > > box_class("Box",init<Box>());
    box_class.def(init<uint>());
    box_class.def(init< Vector<Interval> >());
    box_class.def("__eq__", (bool(*)(const Vector<Interval>&,const Vector<Interval>&)) &operator==);
    box_class.def("dimension", (uint(Box::*)()const) &Box::dimension);
    box_class.def("centre", (Point(Box::*)()const) &Box::centre);
    box_class.def("radius", (Float(Box::*)()const) &Box::radius);
    box_class.def("separated", (tribool(Box::*)(const Box&)const) &Box::separated);
    box_class.def("overlaps", (tribool(Box::*)(const Box&)const) &Box::overlaps);
    box_class.def("covers", (tribool(Box::*)(const Box&)const) &Box::covers);
    box_class.def("inside", (tribool(Box::*)(const Box&)const) &Box::inside);
    box_class.def("empty", (bool(Box::*)()const) &Box::empty);
    box_class.def("widen", (Box(Box::*)()const) &Box::widen);
    box_class.def("split", (std::pair<Box,Box>(Box::*)()const) &Box::split);
    box_class.def("split", (std::pair<Box,Box>(Box::*)(uint)const) &Box::split);
    box_class.def("split", (std::pair<Box,Box>(Box::*)()const) &Box::split);
    box_class.def(self_ns::str(self));

    def("split", (std::pair<Box,Box>(*)(const Box&)) &split);
    def("disjoint", (bool(*)(const IVector&,const IVector&)) &disjoint);
    def("subset", (bool(*)(const IVector&,const IVector&)) &subset);

    def("product", (Box(*)(const Box&,const Interval&)) &product);
    def("product", (Box(*)(const Box&,const Box&)) &product);
    def("hull", (Box(*)(const Box&,const Box&)) &hull);
    def("intersection", (Box(*)(const Box&,const Box&)) &intersection);

    from_python<Box>();
    to_python< std::pair<Box,Box> >();
    implicitly_convertible<Vector<Interval>,Box>();

}

std::pair<Zonotope,Zonotope> split_pair(const Zonotope& z) {
    ListSet<Zonotope> split_list=split(z);
    ARIADNE_ASSERT(split_list.size()==2);
    return std::pair<Zonotope,Zonotope>(split_list[0],split_list[1]);
}

void export_zonotope()
{
    class_<Zonotope,bases<CompactSetInterface,OpenSetInterface,DrawableInterface> > zonotope_class("Zonotope",init<Zonotope>());
    zonotope_class.def(init< Point, Matrix<Float>, Vector<Float> >());
    zonotope_class.def(init< Point, Matrix<Float> >());
    zonotope_class.def(init< Box >());
    zonotope_class.def("centre",&Zonotope::centre,return_value_policy<copy_const_reference>());
    zonotope_class.def("generators",&Zonotope::generators,return_value_policy<copy_const_reference>());
    zonotope_class.def("error",&Zonotope::error,return_value_policy<copy_const_reference>());
    zonotope_class.def("contains",&Zonotope::contains);
    zonotope_class.def("split", (ListSet<Zonotope>(*)(const Zonotope&)) &split);
    zonotope_class.def("__str__",&__cstr__<Zonotope>);

    def("contains", (tribool(*)(const Zonotope&,const Point&)) &contains);
    def("separated", (tribool(*)(const Zonotope&,const Box&)) &separated);
    def("overlaps", (tribool(*)(const Zonotope&,const Box&)) &overlaps);
    def("separated", (tribool(*)(const Zonotope&,const Zonotope&)) &separated);

    def("polytope", (Polytope(*)(const Zonotope&)) &polytope);
    def("orthogonal_approximation", (Zonotope(*)(const Zonotope&)) &orthogonal_approximation);
    def("orthogonal_over_approximation", (Zonotope(*)(const Zonotope&)) &orthogonal_over_approximation);
    def("error_free_over_approximation", (Zonotope(*)(const Zonotope&)) &error_free_over_approximation);

    def("apply", (Zonotope(*)(const ValidatedVectorFunction&, const Zonotope&)) &apply);

    to_python< ListSet<Zonotope> >();

}

void export_polytope()
{
    class_<Polytope,bases<LocatedSetInterface,DrawableInterface> > polytope_class("Polytope",init<Polytope>());
    polytope_class.def(init<int>());
    polytope_class.def("new_vertex",&Polytope::new_vertex);
    polytope_class.def("__iter__",boost::python::range(&Polytope::vertices_begin,&Polytope::vertices_end));
    polytope_class.def(self_ns::str(self));

}

void export_curve()
{
    to_python< std::pair<const Float,Point> >();

    class_<InterpolatedCurve,bases<DrawableInterface> > interpolated_curve_class("InterpolatedCurve",init<InterpolatedCurve>());
    interpolated_curve_class.def(init<Float,Point>());
    interpolated_curve_class.def("insert", &InterpolatedCurve::insert);
    interpolated_curve_class.def("__iter__",boost::python::range(&InterpolatedCurve::begin,&InterpolatedCurve::end));
    interpolated_curve_class.def(self_ns::str(self));


}



void export_affine_set()
{

    class_<ValidatedAffineConstrainedImageSet,bases<CompactSetInterface,DrawableInterface> >
        affine_set_class("ValidatedAffineConstrainedImageSet",init<ValidatedAffineConstrainedImageSet>());
    affine_set_class.def(init<Vector<Interval>, Matrix<Float>, Vector<Float> >());
    affine_set_class.def(init<Matrix<Float>, Vector<Float> >());
    affine_set_class.def("new_parameter_constraint", (void(ValidatedAffineConstrainedImageSet::*)(const Constraint<Affine<Interval>,Float>&)) &ValidatedAffineConstrainedImageSet::new_parameter_constraint);
    affine_set_class.def("new_constraint", (void(ValidatedAffineConstrainedImageSet::*)(const Constraint<AffineModel<Interval>,Float>&)) &ValidatedAffineConstrainedImageSet::new_constraint);
    affine_set_class.def("dimension", &ValidatedAffineConstrainedImageSet::dimension);
    affine_set_class.def("bounded", &ValidatedAffineConstrainedImageSet::bounded);
    affine_set_class.def("empty", &ValidatedAffineConstrainedImageSet::empty);
    affine_set_class.def("bounding_box", &ValidatedAffineConstrainedImageSet::bounding_box);
    affine_set_class.def("separated", &ValidatedAffineConstrainedImageSet::separated);
    affine_set_class.def("adjoin_outer_approximation_to", &ValidatedAffineConstrainedImageSet::adjoin_outer_approximation_to);
    affine_set_class.def("outer_approximation", &ValidatedAffineConstrainedImageSet::outer_approximation);
    affine_set_class.def("boundary", &ValidatedAffineConstrainedImageSet::boundary);
    affine_set_class.def(self_ns::str(self));
}


void export_constraint_set()
{
    from_python< List<EffectiveConstraint> >();

    class_<ConstraintSet,bases<RegularSetInterface> >
        constraint_set_class("ConstraintSet",init<ConstraintSet>());
    constraint_set_class.def(init< List<EffectiveConstraint> >());

    class_<BoundedConstraintSet,bases<DrawableInterface> >
        bounded_constraint_set_class("BoundedConstraintSet",init<BoundedConstraintSet>());
    bounded_constraint_set_class.def(init< BoxSet, List<EffectiveConstraint> >());

    class_<BoxSet>
        box_set_class("BoxSet");
    box_set_class.def(init<IntervalVector>());
    box_set_class.def(init< List<IntervalSet> >());

    def("intersection", (BoundedConstraintSet(*)(const ConstraintSet&,const BoxSet&)) &intersection);

}


void export_constrained_image_set()
{
    from_python< List<ValidatedConstraint> >();

    class_<ValidatedConstrainedImageSet,bases<CompactSetInterface,DrawableInterface> >
        constrained_image_set_class("ValidatedConstrainedImageSet",init<ValidatedConstrainedImageSet>());
    constrained_image_set_class.def(init<Box>());
    constrained_image_set_class.def(init<Box,EffectiveVectorFunction>());
    constrained_image_set_class.def(init<Box,ValidatedVectorFunction>());
    constrained_image_set_class.def(init<Box,ValidatedVectorFunction,List<ValidatedConstraint> >());
    constrained_image_set_class.def(init<Box,ValidatedVectorFunctionModel>());
    constrained_image_set_class.def("domain", &ValidatedConstrainedImageSet::domain,return_value_policy<copy_const_reference>());
    constrained_image_set_class.def("function", &ValidatedConstrainedImageSet::function,return_value_policy<copy_const_reference>());
    constrained_image_set_class.def("constraint", &ValidatedConstrainedImageSet::constraint);
    constrained_image_set_class.def("number_of_parameters", &ValidatedConstrainedImageSet::number_of_parameters);
    constrained_image_set_class.def("number_of_constraints", &ValidatedConstrainedImageSet::number_of_constraints);
    constrained_image_set_class.def("apply", &ValidatedConstrainedImageSet::apply);
    constrained_image_set_class.def("new_space_constraint", (void(ValidatedConstrainedImageSet::*)(const EffectiveConstraint&))&ValidatedConstrainedImageSet::new_space_constraint);
    constrained_image_set_class.def("new_parameter_constraint", (void(ValidatedConstrainedImageSet::*)(const EffectiveConstraint&))&ValidatedConstrainedImageSet::new_parameter_constraint);
    //constrained_image_set_class.def("outer_approximation", &ValidatedConstrainedImageSet::outer_approximation);
    constrained_image_set_class.def("affine_approximation", &ValidatedConstrainedImageSet::affine_approximation);
    	constrained_image_set_class.def("affine_over_approximation", &ValidatedConstrainedImageSet::affine_over_approximation);
    constrained_image_set_class.def("adjoin_outer_approximation_to", &ValidatedConstrainedImageSet::adjoin_outer_approximation_to);
    constrained_image_set_class.def("bounding_box", &ValidatedConstrainedImageSet::bounding_box);
    constrained_image_set_class.def("inside", &ValidatedConstrainedImageSet::inside);
    constrained_image_set_class.def("separated", &ValidatedConstrainedImageSet::separated);
    constrained_image_set_class.def("overlaps", &ValidatedConstrainedImageSet::overlaps);
    constrained_image_set_class.def("split", (Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet>(ValidatedConstrainedImageSet::*)()const) &ValidatedConstrainedImageSet::split);
    constrained_image_set_class.def("split", (Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet>(ValidatedConstrainedImageSet::*)(uint)const) &ValidatedConstrainedImageSet::split);
    constrained_image_set_class.def(self_ns::str(self));
    constrained_image_set_class.def("__repr__", &__cstr__<ValidatedConstrainedImageSet>);

    def("product", (ValidatedConstrainedImageSet(*)(const ValidatedConstrainedImageSet&,const IntervalVector&)) &product);
}




void geometry_submodule() {
    export_set_interface();
    export_point();
    export_box();
    export_zonotope();
    export_polytope();
    export_curve();

    export_affine_set();

    export_constraint_set();
    export_constrained_image_set();

}

