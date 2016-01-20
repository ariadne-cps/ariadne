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

#include "boost_python.h"
#include "utilities.h"

#include "config.h"

#include <boost/python.hpp>

#include "geometry/geometry.h"
#include "output/geometry2d.h"
#include "geometry/point.h"
#include "geometry/curve.h"
#include "geometry/interval.h"
#include "geometry/box.h"
#include "geometry/grid_set.h"
#include "geometry/function_set.h"
#include "geometry/affine_set.h"

#include "hybrid/discrete_location.h"
#include "hybrid/hybrid_set.h"



using namespace boost::python;
using namespace Ariadne;

namespace Ariadne {

typedef RealBox BoxSet;

template<>
struct from_python_dict<ExactIntervalType> {
    from_python_dict() { converter::registry::push_back(&convertible,&construct,type_id<ExactIntervalType>()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyDict_Check(obj_ptr) || len(boost::python::extract<boost::python::dict>(obj_ptr))!=1) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::dict dct = boost::python::extract<boost::python::dict>(obj_ptr);
        boost::python::list lst=dct.items();
        assert(boost::python::len(lst)==1);
        Void* storage = ((converter::rvalue_from_python_storage<ExactIntervalType>*)data)->storage.bytes;
        double lb=boost::python::extract<double>(lst[0][0]); double ub=boost::python::extract<double>(lst[0][1]);
        new (storage) ExactIntervalType(lb,ub);
        data->convertible = storage;
    }
};


template<>
struct from_python_list<ExactIntervalType> {
    from_python_list() { converter::registry::push_back(&convertible,&construct,type_id<ExactIntervalType>()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) || len(boost::python::extract<boost::python::list>(obj_ptr))!=2) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        assert(boost::python::len(lst)==2);
        Void* storage = ((converter::rvalue_from_python_storage<ExactIntervalType>*)data)->storage.bytes;
        double lb=boost::python::extract<double>(lst[0]); double ub=boost::python::extract<double>(lst[1]);
        new (storage) ExactIntervalType(lb,ub);
        data->convertible = storage;
    }
};

template<>
struct from_python<ExactPoint> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<ExactPoint>()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::extract<boost::python::tuple> xtup(obj_ptr);
        boost::python::extract<boost::python::list> xlst(obj_ptr);
        ExactPoint pt;
        if(xtup.check()) {
            boost::python::tuple tup=xtup(); pt=ExactPoint(len(tup));
            for(Int i=0; i!=len(tup); ++i) { pt[i]=Float64Value(boost::python::extract<Float64>(tup[i])); }
        } else if(xlst.check()) {
            boost::python::list lst=xlst(); pt=ExactPoint(len(lst));
            for(Int i=0; i!=len(lst); ++i) { pt[i]=Float64Value(boost::python::extract<Float64>(lst[i])); }
        }
        Void* storage = ((converter::rvalue_from_python_storage<ExactIntervalType>*)data)->storage.bytes;
        new (storage) ExactPoint(pt);
        data->convertible = storage;
    }
};

template<>
struct from_python<ExactBoxType> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<ExactBoxType>()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        Void* storage = ((converter::rvalue_from_python_storage<ExactIntervalType>*)data)->storage.bytes;
        boost::python::list lst=boost::python::extract<boost::python::list>(obj_ptr);
        ExactBoxType* bx_ptr = new (storage) ExactBoxType(len(lst));
        for(Int i=0; i!=len(lst); ++i) { (*bx_ptr)[i]=boost::python::extract<ExactIntervalType>(lst[i]); }
        data->convertible = storage;
    }
};


template<class ES>
struct to_python< ListSet<ES> > {
    to_python() { boost::python::to_python_converter< ListSet<ES>, to_python< ListSet<ES> > >(); }

    static PyObject* convert(const ListSet<ES>& ls) {
        boost::python::list result;
        for(typename ListSet<ES>::ConstIterator iter=ls.begin(); iter!=ls.end(); ++iter) {
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
        for(typename SetType::LocationsConstIterator iter=hls.locations_begin(); iter!=hls.locations_end(); ++iter) {
            result[iter->first]=iter->second;
        }
        return boost::python::incref(boost::python::dict(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyDict_Type; }
};

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Float64Bounds>& x);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ExactIntervalType>& x) {
    ExactIntervalType const& ivl=x.reference(); return os << PythonRepresentation<Float64Bounds>(Float64Bounds(ivl.lower(),ivl.upper()));
}



class OpenSetWrapper
  : public virtual OpenSetInterface, public wrapper< OpenSetInterface >
{
  public:
    OpenSetInterface* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    Sierpinskian covers(const ExactBoxType& r) const { return this->get_override("covers")(); }
    Sierpinskian overlaps(const ExactBoxType& r) const { return this->get_override("overlaps")(); }
    OutputStream& write(OutputStream&) const { return this->get_override("write")(); }
};

class ClosedSetWrapper
  : public virtual ClosedSetInterface, public wrapper< ClosedSetInterface >
{
  public:
    ClosedSetInterface* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    Sierpinskian separated(const ExactBoxType& r) const { return this->get_override("separated")(); }
    OutputStream& write(OutputStream&) const { return this->get_override("write")(); }
};


class OvertSetWrapper
  : public virtual OvertSetInterface, public wrapper< OvertSetInterface >
{
  public:
    OvertSetInterface* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    Sierpinskian overlaps(const ExactBoxType& r) const { return this->get_override("overlaps")(); }
    OutputStream& write(OutputStream&) const { return this->get_override("write")(); }
};


class CompactSetWrapper
  : public virtual CompactSetInterface, public wrapper< CompactSetInterface >
{
  public:
    CompactSetInterface* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    Sierpinskian separated(const ExactBoxType& r) const { return this->get_override("separated")(); }
    Sierpinskian inside(const ExactBoxType& r) const { return this->get_override("inside")(); }
    Sierpinskian is_bounded() const { return this->get_override("is_bounded")(); }
    UpperBoxType bounding_box() const { return this->get_override("bounding_box")(); }
    OutputStream& write(OutputStream&) const { return this->get_override("write")(); }
};

class RegularSetWrapper
  : public virtual LocatedSetInterface, public wrapper< RegularSetWrapper >
{
  public:
    RegularSetWrapper* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    Sierpinskian overlaps(const ExactBoxType& r) const { return this->get_override("overlaps")(); }
    Sierpinskian covers(const ExactBoxType& r) const { return this->get_override("covers")(); }
    Sierpinskian separated(const ExactBoxType& r) const { return this->get_override("separated")(); }
    OutputStream& write(OutputStream&) const { return this->get_override("write")(); }
};

class LocatedSetWrapper
  : public virtual LocatedSetInterface, public wrapper< LocatedSetInterface >
{
  public:
    LocatedSetInterface* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    Sierpinskian overlaps(const ExactBoxType& r) const { return this->get_override("overlaps")(); }
    Sierpinskian separated(const ExactBoxType& r) const { return this->get_override("separated")(); }
    Sierpinskian inside(const ExactBoxType& r) const { return this->get_override("inside")(); }
    Sierpinskian is_bounded() const { return this->get_override("is_bounded")(); }
    UpperBoxType bounding_box() const { return this->get_override("bounding_box")(); }
    OutputStream& write(OutputStream&) const { return this->get_override("write")(); }
};

}


Void export_set_interface() {
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


Void export_point()
{
    class_<ExactPoint,bases<DrawableInterface>> point_class("ExactPoint",init<ExactPoint>());
    point_class.def(init<Nat>());
    point_class.def("__getitem__", &__getitem__<ExactPoint,Int,Float64Value>);
    point_class.def(self_ns::str(self));

    from_python<ExactPoint>();
    implicitly_convertible<Vector<Float64Value>,ExactPoint>();

}

typedef Logical<ExactTag> ExactLogicalType;


template<class IVL> Void export_interval(std::string name) {
    typedef IVL IntervalType;
    typedef typename IntervalType::LowerBoundType LowerBoundType;
    typedef typename IntervalType::UpperBoundType UpperBoundType;
    typedef typename IntervalType::MidpointType MidpointType;

    typedef decltype(contains(declval<IntervalType>(),declval<MidpointType>())) ContainsType;
    typedef decltype(disjoint(declval<IntervalType>(),declval<IntervalType>())) DisjointType;
    typedef decltype(subset(declval<IntervalType>(),declval<IntervalType>())) SubsetType;

    class_< IntervalType > interval_class(name.c_str(),init<IntervalType>());
    //interval_class.def(init<MidpointType>());
    //interval_class.def(init<LowerBoundType,UpperBoundType>());

    interval_class.def(self == self);
    interval_class.def(self != self);
    interval_class.def("lower", &IntervalType::lower, return_value_policy<copy_const_reference>());
    interval_class.def("upper", &IntervalType::upper, return_value_policy<copy_const_reference>());
    interval_class.def("midpoint", &IntervalType::midpoint);
    interval_class.def("radius", &IntervalType::radius);
    interval_class.def("width", &IntervalType::width);
    interval_class.def("contains", (ContainsType(*)(IntervalType const&,MidpointType const&)) &contains);
    interval_class.def("empty", &IntervalType::is_empty);
    interval_class.def(boost::python::self_ns::str(self));

    //from_python_list<IntervalType>();
    //from_python_str<ExactIntervalType>();

    def("midpoint", &IntervalType::midpoint);
    def("radius", &IntervalType::radius);
    def("width", &IntervalType::width);

    def("contains", (ContainsType(*)(IntervalType const&,MidpointType const&)) &contains);
    def("disjoint", (DisjointType(*)(IntervalType const&,IntervalType const&)) &disjoint);
    def("subset", (SubsetType(*)(IntervalType const&,IntervalType const&)) &subset);

    def("intersection", (IntervalType(*)(IntervalType const&,IntervalType const&)) &intersection);
    def("hull", (IntervalType(*)(IntervalType const&, IntervalType const&)) &hull);
}

Void export_intervals() {
    export_interval<ExactIntervalType>("ExactInterval");
    export_interval<UpperIntervalType>("UpperInterval");
    from_python_dict<ExactIntervalType>();
}

template<class BX> Void export_box()
{
    typedef Vector<ExactIntervalType> ExactIntervalVectorType;
    class_<Vector<ExactIntervalType>> interval_vector_class("ExactIntervalVectorType");

//    class_<ExactBoxType,bases<CompactSetInterface,OpenSetInterface,Vector<ExactIntervalType>,DrawableInterface > >
    class_<ExactBoxType,bases< > > box_class("ExactBox",init<ExactBoxType>());
    box_class.def(init<DimensionType>());
    box_class.def(init< Vector<ExactIntervalType> >());
    box_class.def("__eq__", (ExactLogicalType(*)(const Vector<ExactIntervalType>&,const Vector<ExactIntervalType>&)) &operator==);
    box_class.def("dimension", (DimensionType(ExactBoxType::*)()const) &ExactBoxType::dimension);
    box_class.def("centre", (ExactBoxType::CentreType(ExactBoxType::*)()const) &ExactBoxType::centre);
    box_class.def("radius", (ExactBoxType::RadiusType(ExactBoxType::*)()const) &ExactBoxType::radius);
    box_class.def("separated", (ExactLogicalType(ExactBoxType::*)(const ExactBoxType&)const) &ExactBoxType::separated);
    box_class.def("overlaps", (ExactLogicalType(ExactBoxType::*)(const ExactBoxType&)const) &ExactBoxType::overlaps);
    box_class.def("covers", (ExactLogicalType(ExactBoxType::*)(const ExactBoxType&)const) &ExactBoxType::covers);
    box_class.def("inside", (ExactLogicalType(ExactBoxType::*)(const ExactBoxType&)const) &ExactBoxType::inside);
    box_class.def("is_empty", (Bool(ExactBoxType::*)()const) &ExactBoxType::is_empty);
    box_class.def("split", (Pair<ExactBoxType,ExactBoxType>(ExactBoxType::*)()const) &ExactBoxType::split);
    box_class.def("split", (Pair<ExactBoxType,ExactBoxType>(ExactBoxType::*)(SizeType)const) &ExactBoxType::split);
    box_class.def("split", (Pair<ExactBoxType,ExactBoxType>(ExactBoxType::*)()const) &ExactBoxType::split);
    box_class.def(self_ns::str(self));

    def("disjoint", (ExactLogicalType(*)(const ExactBoxType&,const ExactBoxType&)) &disjoint);
    def("subset", (ExactLogicalType(*)(const ExactBoxType&,const ExactBoxType&)) &subset);

    def("product", (ExactBoxType(*)(const ExactBoxType&,const ExactIntervalType&)) &product);
    def("product", (ExactBoxType(*)(const ExactBoxType&,const ExactBoxType&)) &product);
    def("hull", (ExactBoxType(*)(const ExactBoxType&,const ExactBoxType&)) &hull);
    def("intersection", (ExactBoxType(*)(const ExactBoxType&,const ExactBoxType&)) &intersection);

    def("widen", (UpperBoxType(*)(ExactBoxType const&, Float64Value eps)) &widen);

    from_python<ExactBoxType>();
    to_python< Pair<ExactBoxType,ExactBoxType> >();
    implicitly_convertible<Vector<ExactIntervalType>,ExactBoxType>();
}

Void export_boxes() {
    export_box<ExactBoxType>();
}

/*
Pair<Zonotope,Zonotope> split_pair(const Zonotope& z) {
    ListSet<Zonotope> split_list=split(z);
    ARIADNE_ASSERT(split_list.size()==2);
    return Pair<Zonotope,Zonotope>(split_list[0],split_list[1]);
}

Void export_zonotope()
{
    class_<Zonotope,bases<CompactSetInterface,OpenSetInterface,DrawableInterface> > zonotope_class("Zonotope",init<Zonotope>());
    zonotope_class.def(init< Vector<Float64Value>, Matrix<Float64Value>, Vector<Float64Error> >());
    zonotope_class.def(init< Vector<Float64Value>, Matrix<Float64Value> >());
    zonotope_class.def(init< ExactBoxType >());
    zonotope_class.def("centre",&Zonotope::centre,return_value_policy<copy_const_reference>());
    zonotope_class.def("generators",&Zonotope::generators,return_value_policy<copy_const_reference>());
    zonotope_class.def("error",&Zonotope::error,return_value_policy<copy_const_reference>());
    zonotope_class.def("contains",&Zonotope::contains);
    zonotope_class.def("split", (ListSet<Zonotope>(*)(const Zonotope&)) &split);
    zonotope_class.def("__str__",&__cstr__<Zonotope>);

    def("contains", (Kleenean(*)(const Zonotope&,const ExactPoint&)) &contains);
    def("separated", (Kleenean(*)(const Zonotope&,const ExactBoxType&)) &separated);
    def("overlaps", (Kleenean(*)(const Zonotope&,const ExactBoxType&)) &overlaps);
    def("separated", (Kleenean(*)(const Zonotope&,const Zonotope&)) &separated);

    def("polytope", (Polytope(*)(const Zonotope&)) &polytope);
    def("orthogonal_approximation", (Zonotope(*)(const Zonotope&)) &orthogonal_approximation);
    def("orthogonal_over_approximation", (Zonotope(*)(const Zonotope&)) &orthogonal_over_approximation);
    def("error_free_over_approximation", (Zonotope(*)(const Zonotope&)) &error_free_over_approximation);

//    def("apply", (Zonotope(*)(const ValidatedVectorFunction&, const Zonotope&)) &apply);

    to_python< ListSet<Zonotope> >();
}


Void export_polytope()
{
    class_<Polytope,bases<LocatedSetInterface,DrawableInterface> > polytope_class("Polytope",init<Polytope>());
    polytope_class.def(init<Int>());
    polytope_class.def("new_vertex",&Polytope::new_vertex);
    polytope_class.def("__iter__",boost::python::range(&Polytope::vertices_begin,&Polytope::vertices_end));
    polytope_class.def(self_ns::str(self));

}
*/

Void export_curve()
{
    to_python< Pair<const Float64Value,ExactPoint> >();

    class_<InterpolatedCurve,bases<DrawableInterface> > interpolated_curve_class("InterpolatedCurve",init<InterpolatedCurve>());
    interpolated_curve_class.def(init<Float64Value,ExactPoint>());
    interpolated_curve_class.def("insert", (Void(InterpolatedCurve::*)(const Float64Value&, const Point<Float64Approximation>&)) &InterpolatedCurve::insert);
    interpolated_curve_class.def("__iter__",boost::python::range(&InterpolatedCurve::begin,&InterpolatedCurve::end));
    interpolated_curve_class.def(self_ns::str(self));


}



Void export_affine_set()
{

    class_<ValidatedAffineConstrainedImageSet,bases<CompactSetInterface,DrawableInterface> >
        affine_set_class("ValidatedAffineConstrainedImageSet",init<ValidatedAffineConstrainedImageSet>());
    affine_set_class.def(init<Vector<ExactIntervalType>, Matrix<Float64Value>, Vector<Float64Value> >());
    affine_set_class.def(init<Matrix<Float64Value>, Vector<Float64Value> >());
    affine_set_class.def("new_parameter_constraint", (Void(ValidatedAffineConstrainedImageSet::*)(const Constraint<Affine<Float64Bounds>,Float64Bounds>&)) &ValidatedAffineConstrainedImageSet::new_parameter_constraint);
    affine_set_class.def("new_constraint", (Void(ValidatedAffineConstrainedImageSet::*)(const Constraint<AffineModel<Float64Bounds>,Float64Bounds>&)) &ValidatedAffineConstrainedImageSet::new_constraint);
    affine_set_class.def("dimension", &ValidatedAffineConstrainedImageSet::dimension);
    affine_set_class.def("is_bounded", &ValidatedAffineConstrainedImageSet::is_bounded);
    affine_set_class.def("is_empty", &ValidatedAffineConstrainedImageSet::is_empty);
    affine_set_class.def("bounding_box", &ValidatedAffineConstrainedImageSet::bounding_box);
    affine_set_class.def("separated", &ValidatedAffineConstrainedImageSet::separated);
    affine_set_class.def("adjoin_outer_approximation_to", &ValidatedAffineConstrainedImageSet::adjoin_outer_approximation_to);
    affine_set_class.def("outer_approximation", &ValidatedAffineConstrainedImageSet::outer_approximation);
    affine_set_class.def("boundary", &ValidatedAffineConstrainedImageSet::boundary);
    affine_set_class.def(self_ns::str(self));
}


Void export_constraint_set()
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

    def("intersection", (BoundedConstraintSet(*)(const ConstraintSet&,const BoxSet&)) &intersection);

}


Void export_constrained_image_set()
{
    from_python< List<ValidatedConstraint> >();

    class_<ValidatedConstrainedImageSet,bases<CompactSetInterface,DrawableInterface> >
        constrained_image_set_class("ValidatedConstrainedImageSet",init<ValidatedConstrainedImageSet>());
    constrained_image_set_class.def(init<ExactBoxType>());
    constrained_image_set_class.def(init<ExactBoxType,EffectiveVectorFunction>());
    constrained_image_set_class.def(init<ExactBoxType,ValidatedVectorFunction>());
    constrained_image_set_class.def(init<ExactBoxType,ValidatedVectorFunction,List<ValidatedConstraint> >());
    constrained_image_set_class.def(init<ExactBoxType,ValidatedVectorFunctionModel>());
    constrained_image_set_class.def("domain", &ValidatedConstrainedImageSet::domain,return_value_policy<copy_const_reference>());
    constrained_image_set_class.def("function", &ValidatedConstrainedImageSet::function,return_value_policy<copy_const_reference>());
    constrained_image_set_class.def("constraint", &ValidatedConstrainedImageSet::constraint);
    constrained_image_set_class.def("number_of_parameters", &ValidatedConstrainedImageSet::number_of_parameters);
    constrained_image_set_class.def("number_of_constraints", &ValidatedConstrainedImageSet::number_of_constraints);
    constrained_image_set_class.def("apply", &ValidatedConstrainedImageSet::apply);
    constrained_image_set_class.def("new_space_constraint", (Void(ValidatedConstrainedImageSet::*)(const EffectiveConstraint&))&ValidatedConstrainedImageSet::new_space_constraint);
    constrained_image_set_class.def("new_parameter_constraint", (Void(ValidatedConstrainedImageSet::*)(const EffectiveConstraint&))&ValidatedConstrainedImageSet::new_parameter_constraint);
    //constrained_image_set_class.def("outer_approximation", &ValidatedConstrainedImageSet::outer_approximation);
    constrained_image_set_class.def("affine_approximation", &ValidatedConstrainedImageSet::affine_approximation);
    	constrained_image_set_class.def("affine_over_approximation", &ValidatedConstrainedImageSet::affine_over_approximation);
    constrained_image_set_class.def("adjoin_outer_approximation_to", &ValidatedConstrainedImageSet::adjoin_outer_approximation_to);
    constrained_image_set_class.def("bounding_box", &ValidatedConstrainedImageSet::bounding_box);
    constrained_image_set_class.def("inside", &ValidatedConstrainedImageSet::inside);
    constrained_image_set_class.def("separated", &ValidatedConstrainedImageSet::separated);
    constrained_image_set_class.def("overlaps", &ValidatedConstrainedImageSet::overlaps);
    constrained_image_set_class.def("split", (Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet>(ValidatedConstrainedImageSet::*)()const) &ValidatedConstrainedImageSet::split);
    constrained_image_set_class.def("split", (Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet>(ValidatedConstrainedImageSet::*)(Nat)const) &ValidatedConstrainedImageSet::split);
    constrained_image_set_class.def(self_ns::str(self));
    constrained_image_set_class.def("__repr__", &__cstr__<ValidatedConstrainedImageSet>);

//    def("product", (ValidatedConstrainedImageSet(*)(const ValidatedConstrainedImageSet&,const ExactBoxType&)) &product);
}




Void geometry_submodule() {
    export_set_interface();
    export_point();

    export_intervals();
    export_boxes();
//    export_zonotope();
//    export_polytope();
    export_curve();

    export_affine_set();

    export_constraint_set();
    export_constrained_image_set();

}

