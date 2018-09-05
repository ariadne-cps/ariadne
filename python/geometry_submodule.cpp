/***************************************************************************
 *            geometry_submodule.cpp
 *
 *  Copyright 2008--17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11.hpp"
#include "utilities.hpp"

#include "config.hpp"

#include "geometry/geometry.hpp"
#include "output/geometry2d.hpp"
#include "geometry/point.hpp"
#include "geometry/curve.hpp"
#include "geometry/interval.hpp"
#include "geometry/box.hpp"
#include "geometry/grid_paving.hpp"
#include "geometry/function_set.hpp"
#include "geometry/affine_set.hpp"

#include "hybrid/discrete_location.hpp"
#include "hybrid/hybrid_set.hpp"

namespace Ariadne {

//static constexpr auto self = boost::python::self_ns::self;
static constexpr auto self = pybind11::detail::self;

/*
template<class UB>
struct from_python_dict<Interval<UB>> {
    from_python_dict() { converter::registry::push_back(&convertible,&construct,type_id<Interval<UB>>()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyDict_Check(obj_ptr) || len(boost::python::extract<boost::python::dict>(obj_ptr))!=1) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        typedef typename Interval<UB>::LowerBoundType LB;
        boost::python::dict dct = boost::python::extract<boost::python::dict>(obj_ptr);
        boost::python::list lst=dct.items();
        assert(boost::python::len(lst)==1);
        Void* storage = ((converter::rvalue_from_python_storage<Interval<UB>>*)data)->storage.bytes;
        LB lb=boost::python::extract<LB>(lst[0][0]); UB ub=boost::python::extract<UB>(lst[0][1]);
        new (storage) Interval<UB>(lb,ub);
        data->convertible = storage;
    }
};


template<class UB>
struct from_python_list<Interval<UB>> {
    from_python_list() { converter::registry::push_back(&convertible,&construct,type_id<ExactIntervalType>()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) || len(boost::python::extract<boost::python::list>(obj_ptr))!=2) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        typedef typename Interval<UB>::LowerBoundType LB;
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        assert(boost::python::len(lst)==2);
        Void* storage = ((converter::rvalue_from_python_storage<ExactIntervalType>*)data)->storage.bytes;
        LB lb=boost::python::extract<LB>(lst[0]); UB ub=boost::python::extract<UB>(lst[1]);
        new (storage) Interval<UB>(lb,ub);
        data->convertible = storage;
    }
};

template<class X>
struct from_python<Point<X>> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<Point<X>>()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::extract<boost::python::tuple> xtup(obj_ptr);
        boost::python::extract<boost::python::list> xlst(obj_ptr);
        Point<X> pt;
        if(xtup.check()) {
            boost::python::tuple tup=xtup(); pt=Point<X>(len(tup));
            for(Nat i=0; i!=static_cast<Nat>(len(tup)); ++i) { pt[i]=boost::python::extract<X>(tup[i]); }
        } else if(xlst.check()) {
            boost::python::list lst=xlst(); pt=Point<X>(len(lst));
            for(Nat i=0; i!=static_cast<Nat>(len(lst)); ++i) { pt[i]=boost::python::extract<X>(lst[i]); }
        }
        Void* storage = ((converter::rvalue_from_python_storage<X>*)data)->storage.bytes;
        new (storage) Point<X>(pt);
        data->convertible = storage;
    }
};

template<class IVL>
struct from_python<Box<IVL>> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<Box<IVL>>()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        Void* storage = ((converter::rvalue_from_python_storage<IVL>*)data)->storage.bytes;
        boost::python::list lst=boost::python::extract<boost::python::list>(obj_ptr);
        Box<IVL>* bx_ptr = new (storage) Box<IVL>(static_cast<SizeType>(len(lst)));
        for(Int i=0; i!=len(lst); ++i) { (*bx_ptr)[static_cast<SizeType>(i)]=boost::python::extract<IVL>(lst[i]); }
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

*/


OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPBounds>& x);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ExactIntervalType>& x) {
    ExactIntervalType const& ivl=x.reference(); return os << PythonRepresentation<FloatDPBounds>(FloatDPBounds(ivl.lower(),ivl.upper()));
}


class DrawableWrapper
  : public virtual DrawableInterface, public pybind11::wrapper< DrawableInterface >
{
  public:
    virtual DrawableInterface* clone() const { return this->get_override("clone"); }
    virtual Void draw(CanvasInterface& c, const Projection2d& p) const { this->get_override("draw"); }
    virtual DimensionType dimension() const { return this->get_override("dimension"); }
    virtual OutputStream& write(OutputStream& os) const { return this->get_override("write"); }
};

class OpenSetWrapper
  : public virtual OpenSetInterface, public pybind11::wrapper< OpenSetInterface >
{
  public:
    OpenSetInterface* clone() const { return this->get_override("clone"); }
    SizeType dimension() const { return this->get_override("dimension"); }
    LowerKleenean covers(const ExactBoxType& r) const { return this->get_override("covers"); }
    LowerKleenean overlaps(const ExactBoxType& r) const { return this->get_override("overlaps"); }
    OutputStream& write(OutputStream&) const { return this->get_override("write"); }
};

class ClosedSetWrapper
  : public virtual ClosedSetInterface, public pybind11::wrapper< ClosedSetInterface >
{
  public:
    ClosedSetInterface* clone() const { return this->get_override("clone"); }
    SizeType dimension() const { return this->get_override("dimension"); }
    LowerKleenean separated(const ExactBoxType& r) const { return this->get_override("separated"); }
    OutputStream& write(OutputStream&) const { return this->get_override("write"); }
};


class OvertSetWrapper
  : public virtual OvertSetInterface, public pybind11::wrapper< OvertSetInterface >
{
  public:
    OvertSetInterface* clone() const { return this->get_override("clone"); }
    SizeType dimension() const { return this->get_override("dimension"); }
    LowerKleenean overlaps(const ExactBoxType& r) const { return this->get_override("overlaps"); }
    OutputStream& write(OutputStream&) const { return this->get_override("write"); }
};


class CompactSetWrapper
  : public virtual CompactSetInterface, public pybind11::wrapper< CompactSetInterface >
{
  public:
    CompactSetInterface* clone() const { return this->get_override("clone"); }
    SizeType dimension() const { return this->get_override("dimension"); }
    LowerKleenean separated(const ExactBoxType& r) const { return this->get_override("separated"); }
    LowerKleenean inside(const ExactBoxType& r) const { return this->get_override("inside"); }
    LowerKleenean is_bounded() const { return this->get_override("is_bounded"); }
    UpperBoxType bounding_box() const { return this->get_override("bounding_box"); }
    OutputStream& write(OutputStream&) const { return this->get_override("write"); }
};

class RegularSetWrapper
  : public virtual LocatedSetInterface, public pybind11::wrapper< RegularSetWrapper >
{
  public:
    RegularSetWrapper* clone() const { return this->get_override("clone"); }
    SizeType dimension() const { return this->get_override("dimension"); }
    LowerKleenean overlaps(const ExactBoxType& r) const { return this->get_override("overlaps"); }
    LowerKleenean covers(const ExactBoxType& r) const { return this->get_override("covers"); }
    LowerKleenean separated(const ExactBoxType& r) const { return this->get_override("separated"); }
    OutputStream& write(OutputStream&) const { return this->get_override("write"); }
};

class LocatedSetWrapper
  : public virtual LocatedSetInterface, public pybind11::wrapper< LocatedSetInterface >
{
  public:
    LocatedSetInterface* clone() const { return this->get_override("clone"); }
    SizeType dimension() const { return this->get_override("dimension"); }
    LowerKleenean overlaps(const ExactBoxType& r) const { return this->get_override("overlaps"); }
    LowerKleenean separated(const ExactBoxType& r) const { return this->get_override("separated"); }
    LowerKleenean inside(const ExactBoxType& r) const { return this->get_override("inside"); }
    LowerKleenean is_bounded() const { return this->get_override("is_bounded"); }
    UpperBoxType bounding_box() const { return this->get_override("bounding_box"); }
    OutputStream& write(OutputStream&) const { return this->get_override("write"); }
};


} // namespace Ariadne

using namespace Ariadne;



Void export_drawable_interface(pybind11::module& module) {
    pybind11::class_<DrawableInterface,DrawableWrapper> drawable_class(module, "DrawableInterface");
    drawable_class.def("clone", &DrawableInterface::clone);
    drawable_class.def("draw", &DrawableInterface::draw);
    drawable_class.def("dimension", &DrawableInterface::dimension);
}


Void export_set_interface(pybind11::module& module) {
    pybind11::class_<OpenSetInterface, OpenSetWrapper> open_set_wrapper_class(module,"OpenSetInterface");
    open_set_wrapper_class.def("covers",(LowerKleenean(OpenSetInterface::*)(const ExactBoxType& bx)const) &OpenSetInterface::covers);
#warning Export inherited members
    //    open_set_wrapper_class.def("overlaps",(LowerKleenean(OvertSetInterface::*)(const ExactBoxType& bx)const) &OvertSetInterface::overlaps);

    pybind11::class_<ClosedSetInterface, ClosedSetWrapper> closed_set_wrapper_class(module,"ClosedSetInterface");
    closed_set_wrapper_class.def("separated",(LowerKleenean(ClosedSetInterface::*)(const ExactBoxType& bx)const) &ClosedSetInterface::separated);

    pybind11::class_<OvertSetInterface, OvertSetWrapper> overt_set_wrapper_class(module,"OvertSetInterface");
    overt_set_wrapper_class.def("overlaps",(LowerKleenean(OvertSetInterface::*)(const ExactBoxType& bx)const) &OvertSetInterface::overlaps);

    pybind11::class_<CompactSetInterface, CompactSetWrapper> compact_set_wrapper_class(module,"CompactSetInterface");
//    compact_set_wrapper_class.def("separated",(LowerKleenean(ClosedSetInterface::*)(const ExactBoxType& bx)const) &CompactSetInterface::separated);
//    compact_set_wrapper_class.def("inside",(LowerKleenean(BoundedSetInterface::*)(const ExactBoxType& bx)const) &CompactSetInterface::inside);
//    compact_set_wrapper_class.def("bounding_box",&CompactSetInterface::bounding_box);

    pybind11::class_<LocatedSetInterface, pybind11::noncopyable> located_set_wrapper_class(module,"LocatedSetInterface");
    pybind11::class_<RegularSetInterface, pybind11::noncopyable> regular_set_wrapper_class(module,"RegularSetInterface");
}


Void export_point(pybind11::module& module)
{
    //static_assert(pybind11::class_<ExactPoint>::is_valid_class_option<DrawableInterface>::value);
    //static_assert(not pybind11::class_<ExactPoint>::is_valid_class_option<DrawableWrapper>::value);

//    pybind11::class_<ExactPoint,pybind11::bases<DrawableWrapper>> point_class(module,"ExactPoint");
//    pybind11::class_<ExactPoint, DrawableWrapper> point_class(module,"ExactPoint", pybind11::multiple_inheritance());
    pybind11::class_<ExactPoint, DrawableInterface> point_class(module,"ExactPoint", pybind11::multiple_inheritance());
    point_class.def(pybind11::init<ExactPoint>());
    point_class.def(pybind11::init<Nat>());
    point_class.def("__getitem__", &__getitem__<ExactPoint,Int,FloatDPValue>);
    point_class.def("__str__", &__cstr__<ExactPoint>);
    
//    from_python<ExactPoint>();
    pybind11::implicitly_convertible<Vector<FloatDPValue>,ExactPoint>();

}

typedef LogicalType<ExactTag> ExactLogicalType;


template<class IVL> Void export_interval(pybind11::module& module, std::string name) {
    typedef IVL IntervalType;
    typedef typename IntervalType::LowerBoundType LowerBoundType;
    typedef typename IntervalType::UpperBoundType UpperBoundType;
    typedef typename IntervalType::MidpointType MidpointType;

    typedef decltype(contains(declval<IntervalType>(),declval<MidpointType>())) ContainsType;
    typedef decltype(disjoint(declval<IntervalType>(),declval<IntervalType>())) DisjointType;
    typedef decltype(subset(declval<IntervalType>(),declval<IntervalType>())) SubsetType;

//    from_python_dict<IVL>();

    pybind11::class_< IntervalType > interval_class(module,name.c_str());
    interval_class.def(pybind11::init<IntervalType>());
    //interval_class.def(pybind11::init<MidpointType>());
    interval_class.def(pybind11::init<LowerBoundType,UpperBoundType>());

    // FIXME: Only export this if constructor exists
//    if constexpr (IsConstructibleGivenDefaultPrecision<UB,Dyadic>::value and not IsConstructible<UB,Dyadic>::value) {
        interval_class.def(pybind11::init<Interval<Dyadic>>());
//    }

    interval_class.def(self == self);
    interval_class.def(self != self);
    interval_class.def("lower", &IntervalType::lower);
    interval_class.def("upper", &IntervalType::upper);
    interval_class.def("midpoint", &IntervalType::midpoint);
    interval_class.def("radius", &IntervalType::radius);
    interval_class.def("width", &IntervalType::width);
    interval_class.def("contains", (ContainsType(*)(IntervalType const&,MidpointType const&)) &contains);
    interval_class.def("empty", &IntervalType::is_empty);
    interval_class.def("__str__",&__cstr__<IntervalType>);
//    interval_class.def("__repr__",&__repr__<IntervalType>);

    //from_python_list<IntervalType>();
    //from_python_str<ExactIntervalType>();

    module.def("midpoint", &IntervalType::midpoint);
    module.def("radius", &IntervalType::radius);
    module.def("width", &IntervalType::width);

    module.def("contains", (ContainsType(*)(IntervalType const&,MidpointType const&)) &contains);
    module.def("disjoint", (DisjointType(*)(IntervalType const&,IntervalType const&)) &disjoint);
    module.def("subset", (SubsetType(*)(IntervalType const&,IntervalType const&)) &subset);

    module.def("intersection", (IntervalType(*)(IntervalType const&,IntervalType const&)) &intersection);
    module.def("hull", (IntervalType(*)(IntervalType const&, IntervalType const&)) &hull);
}

Void export_intervals(pybind11::module& module) {
    export_interval<ExactIntervalType>(module,"ExactInterval");
    export_interval<UpperIntervalType>(module,"UpperInterval");
    export_interval<ApproximateIntervalType>(module,"ApproximateInterval");
    export_interval<DyadicInterval>(module,"DyadicInterval");
    export_interval<RealInterval>(module,"RealInterval");
}

template<class BX> Void export_box(pybind11::module& module, std::string name)
{
    using IVL = typename BX::IntervalType;
    using UB = typename IVL::UpperBoundType;
    //class_<Vector<ExactIntervalType>> interval_vector_class(module,"ExactIntervalVectorType");

    typedef decltype(disjoint(declval<BX>(),declval<BX>())) DisjointType;
    typedef decltype(subset(declval<BX>(),declval<BX>())) SubsetType;
    typedef decltype(separated(declval<BX>(),declval<BX>())) SeparatedType;
    typedef decltype(overlap(declval<BX>(),declval<BX>())) OverlapType;
    typedef decltype(covers(declval<BX>(),declval<BX>())) CoversType;
    typedef decltype(inside(declval<BX>(),declval<BX>())) InsideType;

//    pybind11::class_<ExactBoxType,pybind11::bases<CompactSetWrapper,OpenSetWrapper,Vector<ExactIntervalType>,DrawableWrapper > >
    pybind11::class_<BX > box_class(module,name.c_str());
    box_class.def(pybind11::init<BX>());
    if constexpr (IsConstructibleGivenDefaultPrecision<UB,Dyadic>::value and not IsConstructible<UB,Dyadic>::value) {
        box_class.def(pybind11::init<Box<Interval<Dyadic>>>());
    }

    static_assert(IsConstructibleGivenDefaultPrecision<FloatDPValue,Dyadic>::value and not IsConstructible<FloatDPValue,Dyadic>::value);

    box_class.def(pybind11::init<DimensionType>());
    box_class.def(pybind11::init< Vector<IVL> >());
    //box_class.def("__eq__", (ExactLogicalType(*)(const Vector<ExactIntervalType>&,const Vector<ExactIntervalType>&)) &operator==);
    box_class.def("dimension", (DimensionType(BX::*)()const) &BX::dimension);
    box_class.def("centre", (typename BX::CentreType(BX::*)()const) &BX::centre);
    box_class.def("radius", (typename BX::RadiusType(BX::*)()const) &BX::radius);
    box_class.def("separated", (SeparatedType(BX::*)(const BX&)const) &BX::separated);
    box_class.def("overlaps", (OverlapType(BX::*)(const BX&)const) &BX::overlaps);
    box_class.def("covers", (CoversType(BX::*)(const BX&)const) &BX::covers);
    box_class.def("inside", (InsideType(BX::*)(const BX&)const) &BX::inside);
    box_class.def("is_empty", (SeparatedType(BX::*)()const) &BX::is_empty);
    box_class.def("split", (Pair<BX,BX>(BX::*)()const) &BX::split);
    box_class.def("split", (Pair<BX,BX>(BX::*)(SizeType)const) &BX::split);
    box_class.def("split", (Pair<BX,BX>(BX::*)()const) &BX::split);
    box_class.def("__str__",&__cstr__<BX>);

    module.def("disjoint", (DisjointType(*)(const BX&,const BX&)) &disjoint);
    module.def("subset", (SubsetType(*)(const BX&,const BX&)) &subset);

    module.def("product", (BX(*)(const BX&,const IVL&)) &product);
    module.def("product", (BX(*)(const BX&,const BX&)) &product);
    module.def("hull", (BX(*)(const BX&,const BX&)) &hull);
    module.def("intersection", (BX(*)(const BX&,const BX&)) &intersection);

//    from_python<BX>();
//    to_python< Pair<BX,BX> >();
}

template<> Void export_box<DyadicBox>(pybind11::module& module, std::string name)
{
    using BX=DyadicBox;
    pybind11::class_<BX> box_class(module,name.c_str());
    box_class.def(pybind11::init<BX>());
//    from_python<BX>();
}

Void export_boxes(pybind11::module& module) {
    export_box<RealBox>(module,"RealBox");
    export_box<DyadicBox>(module,"DyadicBox");
    export_box<ExactBoxType>(module,"ExactBox");
    export_box<UpperBoxType>(module,"UpperBox");
    export_box<ApproximateBoxType>(module,"ApproximateBox");

    pybind11::implicitly_convertible<ExactBoxType,UpperBoxType>();
    pybind11::implicitly_convertible<ExactBoxType,ApproximateBoxType>();
    pybind11::implicitly_convertible<UpperBoxType,ApproximateBoxType>();

    module.def("widen", (UpperBoxType(*)(ExactBoxType const&, FloatDPValue eps)) &widen);
    module.def("image", (UpperBoxType(*)(UpperBoxType const&, ValidatedVectorFunction const&)) &_image_);
}

/*
Pair<Zonotope,Zonotope> split_pair(const Zonotope& z) {
    ListSet<Zonotope> split_list=split(z);
    ARIADNE_ASSERT(split_list.size()==2);
    return Pair<Zonotope,Zonotope>(split_list[0],split_list[1]);
}

Void export_zonotope(pybind11::module& module)
{
    pybind11::class_<Zonotope,pybind11::bases<CompactSetWrapper,OpenSetWrapper,DrawableWrapper> > zonotope_class(module,"Zonotope",pybind11::init<Zonotope>());
    zonotope_class.def(pybind11::init< Vector<FloatDPValue>, Matrix<FloatDPValue>, Vector<FloatDPError> >());
    zonotope_class.def(pybind11::init< Vector<FloatDPValue>, Matrix<FloatDPValue> >());
    zonotope_class.def(pybind11::init< ExactBoxType >());
    zonotope_class.def("centre",&Zonotope::centre,return_value_policy<copy_const_reference>());
    zonotope_class.def("generators",&Zonotope::generators,return_value_policy<copy_const_reference>());
    zonotope_class.def("error",&Zonotope::error,return_value_policy<copy_const_reference>());
    zonotope_class.def("contains",&Zonotope::contains);
    zonotope_class.def("split", (ListSet<Zonotope>(*)(const Zonotope&)) &split);
    zonotope_class.def("__str__",&__cstr__<Zonotope>);

    module.def("contains", (ValidatedKleenean(*)(const Zonotope&,const ExactPoint&)) &contains);
    module.def("separated", (ValidatedKleenean(*)(const Zonotope&,const ExactBoxType&)) &separated);
    module.def("overlaps", (ValidatedKleenean(*)(const Zonotope&,const ExactBoxType&)) &overlaps);
    module.def("separated", (ValidatedKleenean(*)(const Zonotope&,const Zonotope&)) &separated);

    module.def("polytope", (Polytope(*)(const Zonotope&)) &polytope);
    module.def("orthogonal_approximation", (Zonotope(*)(const Zonotope&)) &orthogonal_approximation);
    module.def("orthogonal_over_approximation", (Zonotope(*)(const Zonotope&)) &orthogonal_over_approximation);
    module.def("error_free_over_approximation", (Zonotope(*)(const Zonotope&)) &error_free_over_approximation);

//    module.def("apply", (Zonotope(*)(const ValidatedVectorFunction&, const Zonotope&)) &apply);

    to_python< ListSet<Zonotope> >();
}


Void export_polytope(pybind11::module& module)
{
    pybind11::class_<Polytope,pybind11::bases<LocatedSetWrapper,DrawableWrapper> > polytope_class(module,"Polytope",pybind11::init<Polytope>());
    polytope_class.def(pybind11::init<Int>());
    polytope_class.def("new_vertex",&Polytope::new_vertex);
    polytope_class.def("__iter__",boost::python::range(&Polytope::vertices_begin,&Polytope::vertices_end));
    polytope_class.def(self_ns::str(self));

}
*/

Void export_curve(pybind11::module& module)
{
//    to_python< Pair<const FloatDPValue,ExactPoint> >();

    //pybind11::class_<InterpolatedCurve> interpolated_curve_class(module,"InterpolatedCurve", drawable_class);    
    pybind11::class_<InterpolatedCurve, DrawableInterface> interpolated_curve_class(module,"InterpolatedCurve", pybind11::multiple_inheritance());
    interpolated_curve_class.def(pybind11::init<InterpolatedCurve>());
    interpolated_curve_class.def(pybind11::init<FloatDPValue,ExactPoint>());
    interpolated_curve_class.def("insert", (Void(InterpolatedCurve::*)(const FloatDPValue&, const Point<FloatDPApproximation>&)) &InterpolatedCurve::insert);
//    interpolated_curve_class.def("__iter__",boost::python::range(&InterpolatedCurve::begin,&InterpolatedCurve::end));
    interpolated_curve_class.def("__str__", &__cstr__<InterpolatedCurve>);


}



Void export_affine_set(pybind11::module& module)
{
//    pybind11::class_<ValidatedAffineConstrainedImageSet,pybind11::bases<CompactSetInterface,DrawableInterface> >
    pybind11::class_<ValidatedAffineConstrainedImageSet,pybind11::bases<DrawableInterface> >
        affine_set_class(module,"ValidatedAffineConstrainedImageSet", pybind11::multiple_inheritance());
    affine_set_class.def(pybind11::init<ValidatedAffineConstrainedImageSet>());
    affine_set_class.def(pybind11::init<RealBox>());
    affine_set_class.def(pybind11::init<ExactBoxType>());
    affine_set_class.def(pybind11::init<Vector<ExactIntervalType>, Matrix<FloatDPValue>, Vector<FloatDPValue> >());
    affine_set_class.def(pybind11::init<Matrix<FloatDPValue>, Vector<FloatDPValue> >());
    affine_set_class.def("new_parameter_constraint", (Void(ValidatedAffineConstrainedImageSet::*)(const Constraint<Affine<FloatDPBounds>,FloatDPBounds>&)) &ValidatedAffineConstrainedImageSet::new_parameter_constraint);
    affine_set_class.def("new_constraint", (Void(ValidatedAffineConstrainedImageSet::*)(const Constraint<AffineModel<ValidatedTag,FloatDP>,FloatDPBounds>&)) &ValidatedAffineConstrainedImageSet::new_constraint);
    affine_set_class.def("dimension", &ValidatedAffineConstrainedImageSet::dimension);
    affine_set_class.def("is_bounded", &ValidatedAffineConstrainedImageSet::is_bounded);
    affine_set_class.def("is_empty", &ValidatedAffineConstrainedImageSet::is_empty);
    affine_set_class.def("bounding_box", &ValidatedAffineConstrainedImageSet::bounding_box);
    affine_set_class.def("separated", &ValidatedAffineConstrainedImageSet::separated);
    affine_set_class.def("adjoin_outer_approximation_to", &ValidatedAffineConstrainedImageSet::adjoin_outer_approximation_to);
    affine_set_class.def("outer_approximation", &ValidatedAffineConstrainedImageSet::outer_approximation);
    affine_set_class.def("boundary", &ValidatedAffineConstrainedImageSet::boundary);
    affine_set_class.def("__str__",&__cstr__<ValidatedAffineConstrainedImageSet>);
    
    module.def("image", (ValidatedAffineConstrainedImageSet(*)(ValidatedAffineConstrainedImageSet const&,ValidatedVectorFunction const&)) &_image_);
}

Void export_constraint_set(pybind11::module& module)
{    
//    from_python< List<EffectiveConstraint> >();

    pybind11::class_<ConstraintSet,pybind11::bases<RegularSetInterface,OpenSetInterface> >
        constraint_set_class(module,"ConstraintSet", pybind11::multiple_inheritance());
    constraint_set_class.def(pybind11::init<ConstraintSet>());
    constraint_set_class.def(pybind11::init< List<EffectiveConstraint> >());
    constraint_set_class.def("dimension", &ConstraintSet::dimension);
    constraint_set_class.def("__str__", &__cstr__<ConstraintSet>);
    
//    pybind11::class_<BoundedConstraintSet,pybind11::bases<DrawableWrapper> >
    pybind11::class_<BoundedConstraintSet,pybind11::bases<RegularSetInterface,CompactSetInterface,DrawableInterface> >
        bounded_constraint_set_class(module,"BoundedConstraintSet", pybind11::multiple_inheritance());
    bounded_constraint_set_class.def(pybind11::init<BoundedConstraintSet>());
    bounded_constraint_set_class.def(pybind11::init< RealBox, List<EffectiveConstraint> >());
    bounded_constraint_set_class.def("dimension", &BoundedConstraintSet::dimension);
    bounded_constraint_set_class.def("__str__", &__cstr__<BoundedConstraintSet>);

    module.def("intersection", (ConstraintSet(*)(ConstraintSet const&,ConstraintSet const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(ConstraintSet const&, RealBox const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(RealBox const&, ConstraintSet const&)) &_intersection_);

    module.def("intersection", (BoundedConstraintSet(*)(BoundedConstraintSet const&, BoundedConstraintSet const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(BoundedConstraintSet const&, RealBox const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(RealBox const&, BoundedConstraintSet const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(ConstraintSet const&, BoundedConstraintSet const&)) &_intersection_);
    module.def("intersection", (BoundedConstraintSet(*)(BoundedConstraintSet const&, ConstraintSet const&)) &_intersection_);

    module.def("image", (ConstrainedImageSet(*)(BoundedConstraintSet const&, EffectiveVectorFunction const&)) &_image_);

}


Void export_constrained_image_set(pybind11::module& module)
{
//    from_python< List<ValidatedConstraint> >();

    pybind11::class_<ConstrainedImageSet,pybind11::bases<CompactSetInterface,DrawableInterface> >
        constrained_image_set_class(module,"ConstrainedImageSet");
    constrained_image_set_class.def(pybind11::init<ConstrainedImageSet>());
    constrained_image_set_class.def("dimension", &ConstrainedImageSet::dimension);
//    	constrained_image_set_class.def("affine_over_approximation", &ValidatedConstrainedImageSet::affine_over_approximation);
    constrained_image_set_class.def("__str__",&__cstr__<ConstrainedImageSet>);
    
//    pybind11::class_<ValidatedConstrainedImageSet,pybind11::bases<CompactSetInterface,DrawableInterface> >
    pybind11::class_<ValidatedConstrainedImageSet,pybind11::bases<DrawableInterface> >
        validated_constrained_image_set_class(module,"ValidatedConstrainedImageSet", pybind11::multiple_inheritance());
    validated_constrained_image_set_class.def(pybind11::init<ValidatedConstrainedImageSet>());
    validated_constrained_image_set_class.def(pybind11::init<ExactBoxType>());
    validated_constrained_image_set_class.def(pybind11::init<ExactBoxType,EffectiveVectorFunction>());
    validated_constrained_image_set_class.def(pybind11::init<ExactBoxType,ValidatedVectorFunction>());
    validated_constrained_image_set_class.def(pybind11::init<ExactBoxType,ValidatedVectorFunction,List<ValidatedConstraint> >());
    validated_constrained_image_set_class.def(pybind11::init<ExactBoxType,ValidatedVectorFunctionModelDP>());
    validated_constrained_image_set_class.def("domain", &ValidatedConstrainedImageSet::domain);
    validated_constrained_image_set_class.def("function", &ValidatedConstrainedImageSet::function);
    validated_constrained_image_set_class.def("constraint", &ValidatedConstrainedImageSet::constraint);
    validated_constrained_image_set_class.def("number_of_parameters", &ValidatedConstrainedImageSet::number_of_parameters);
    validated_constrained_image_set_class.def("number_of_constraints", &ValidatedConstrainedImageSet::number_of_constraints);
    validated_constrained_image_set_class.def("apply", &ValidatedConstrainedImageSet::apply);
    validated_constrained_image_set_class.def("new_space_constraint", (Void(ValidatedConstrainedImageSet::*)(const ValidatedConstraint&))&ValidatedConstrainedImageSet::new_space_constraint);
    validated_constrained_image_set_class.def("new_parameter_constraint", (Void(ValidatedConstrainedImageSet::*)(const ValidatedConstraint&))&ValidatedConstrainedImageSet::new_parameter_constraint);
    //constrained_image_set_class.def("outer_approximation", &ValidatedConstrainedImageSet::outer_approximation);
    validated_constrained_image_set_class.def("affine_approximation", &ValidatedConstrainedImageSet::affine_approximation);
    validated_constrained_image_set_class.def("adjoin_outer_approximation_to", &ValidatedConstrainedImageSet::adjoin_outer_approximation_to);
    validated_constrained_image_set_class.def("bounding_box", &ValidatedConstrainedImageSet::bounding_box);
    validated_constrained_image_set_class.def("inside", &ValidatedConstrainedImageSet::inside);
    validated_constrained_image_set_class.def("separated", &ValidatedConstrainedImageSet::separated);
    validated_constrained_image_set_class.def("overlaps", &ValidatedConstrainedImageSet::overlaps);
    validated_constrained_image_set_class.def("split", (Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet>(ValidatedConstrainedImageSet::*)()const) &ValidatedConstrainedImageSet::split);
    validated_constrained_image_set_class.def("split", (Pair<ValidatedConstrainedImageSet,ValidatedConstrainedImageSet>(ValidatedConstrainedImageSet::*)(Nat)const) &ValidatedConstrainedImageSet::split);
    validated_constrained_image_set_class.def("__str__", &__cstr__<ValidatedConstrainedImageSet>);
    validated_constrained_image_set_class.def("__repr__", &__cstr__<ValidatedConstrainedImageSet>);

//    module.def("product", (ValidatedConstrainedImageSet(*)(const ValidatedConstrainedImageSet&,const ExactBoxType&)) &product);
}




Void geometry_submodule(pybind11::module& module) {
    export_drawable_interface(module);
    export_set_interface(module);
    export_point(module);

    export_intervals(module);
    export_boxes(module);
//    export_zonotope(module);
//    export_polytope(module);
    export_curve(module);

    export_affine_set(module);

    export_constraint_set(module);
    export_constrained_image_set(module);

}

