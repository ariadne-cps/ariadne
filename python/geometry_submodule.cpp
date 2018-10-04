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


OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPBounds>& x);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ExactIntervalType>& x) {
    ExactIntervalType const& ivl=x.reference(); return os << PythonRepresentation<FloatDPBounds>(FloatDPBounds(ivl.lower(),ivl.upper()));
}


class DrawableWrapper
  : public pybind11::wrapper< DrawableInterface >
{
  public:
    virtual DrawableInterface* clone() const { return this->get_override("clone")(); }
    virtual Void draw(CanvasInterface& c, const Projection2d& p) const { this->get_override("draw")(c,p); }
    virtual DimensionType dimension() const { return this->get_override("dimension")(); }
    virtual OutputStream& write(OutputStream& os) const { return this->get_override("write")(os); }
};

class OpenSetWrapper
  : public pybind11::wrapper<OpenSetInterface>
{
  public:
    OpenSetInterface* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean covers(const ExactBoxType& r) const { return this->get_override("covers")(r); }
    LowerKleenean overlaps(const ExactBoxType& r) const { return this->get_override("overlaps")(r); }
    OutputStream& write(OutputStream& os) const { return this->get_override("write")(os); }
};

class ClosedSetWrapper
  : public pybind11::wrapper<ClosedSetInterface>
{
  public:
    ClosedSetInterface* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean separated(const ExactBoxType& r) const { return this->get_override("separated")(r); }
    OutputStream& write(OutputStream& os) const { return this->get_override("write")(os); }
};


class OvertSetWrapper
  : public pybind11::wrapper<OvertSetInterface>
{
  public:
    OvertSetInterface* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean overlaps(const ExactBoxType& r) const { return this->get_override("overlaps")(r); }
    OutputStream& write(OutputStream& os) const { return this->get_override("write")(os); }
};


class CompactSetWrapper
  : public pybind11::wrapper<CompactSetInterface>
{
  public:
    CompactSetInterface* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean separated(const ExactBoxType& r) const { return this->get_override("separated")(); }
    LowerKleenean inside(const ExactBoxType& r) const { return this->get_override("inside")(); }
    LowerKleenean is_bounded() const { return this->get_override("is_bounded")(); }
    UpperBoxType bounding_box() const { return this->get_override("bounding_box")(); }
    OutputStream& write(OutputStream& os) const { return this->get_override("write")(); }
};

class RegularSetWrapper
  : public pybind11::wrapper<RegularSetInterface>
{
  public:
    RegularSetWrapper* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean overlaps(const ExactBoxType& r) const { return this->get_override("overlaps")(r); }
    LowerKleenean covers(const ExactBoxType& r) const { return this->get_override("covers")(r); }
    LowerKleenean separated(const ExactBoxType& r) const { return this->get_override("separated")(r); }
    OutputStream& write(OutputStream& os) const { return this->get_override("write")(os); }
};

class LocatedSetWrapper
  : public pybind11::wrapper<LocatedSetInterface>
{
  public:
    LocatedSetInterface* clone() const { return this->get_override("clone")(); }
    SizeType dimension() const { return this->get_override("dimension")(); }
    LowerKleenean overlaps(const ExactBoxType& r) const { return this->get_override("overlaps")(r); }
    LowerKleenean separated(const ExactBoxType& r) const { return this->get_override("separated")(r); }
    LowerKleenean inside(const ExactBoxType& r) const { return this->get_override("inside")(r); }
    LowerKleenean is_bounded() const { return this->get_override("is_bounded")(); }
    UpperBoxType bounding_box() const { return this->get_override("bounding_box")(); }
    OutputStream& write(OutputStream& os) const { return this->get_override("write")(os); }
};


} // namespace Ariadne


using namespace Ariadne;


template<class IVL> IVL interval_from_dict(pybind11::dict dct) {
    typedef typename IVL::LowerBoundType LB;
    typedef typename IVL::UpperBoundType UB;
    assert(dct.size()==1);
    pybind11::detail::dict_iterator::reference item = *dct.begin();
    pybind11::handle lh = item.first;
    pybind11::handle uh = item.second;
    if constexpr (IsConstructibleGivenDefaultPrecision<UB,Dyadic>::value) {
        typedef PrecisionType<UB> PR; PR pr;
        try {
            LB lb(pybind11::cast<Dyadic>(lh),pr);
            UB ub(pybind11::cast<Dyadic>(uh),pr);
            return IVL(lb,ub);
        }
        catch(pybind11::cast_error) {
        }
    }
    LB lb = pybind11::cast<LB>(lh);
    UB ub = pybind11::cast<UB>(uh);
    return IVL(lb,ub);

}

template<class BX> BX box_from_list(pybind11::list lst) {
    typedef typename BX::IntervalType IVL;
    std::vector<IVL> vec=pybind11::cast<std::vector<IVL>>(lst);
    Array<IVL> ary(vec.begin(),vec.end());
    return BX(ary);
}


Void export_drawable_interface(pybind11::module& module) {
    pybind11::class_<DrawableInterface,DrawableWrapper> drawable_class(module, "DrawableInterface");
    drawable_class.def("clone", &DrawableInterface::clone);
    drawable_class.def("draw", &DrawableInterface::draw);
    drawable_class.def("dimension", &DrawableInterface::dimension);
}


Void export_set_interface(pybind11::module& module) {
    pybind11::class_<OvertSetInterface, OvertSetWrapper> overt_set_interface_class(module,"OvertSetInterface");
    overt_set_interface_class.def("overlaps",(LowerKleenean(OvertSetInterface::*)(const ExactBoxType& bx)const) &OvertSetInterface::overlaps);

    pybind11::class_<OpenSetInterface, OpenSetWrapper> open_set_interface_class(module,"OpenSetInterface",overt_set_interface_class);
    open_set_interface_class.def("covers",(LowerKleenean(OpenSetInterface::*)(const ExactBoxType& bx)const) &OpenSetInterface::covers);

    pybind11::class_<ClosedSetInterface, ClosedSetWrapper> closed_set_interface_class(module,"ClosedSetInterface");
    closed_set_interface_class.def("separated",(LowerKleenean(ClosedSetInterface::*)(const ExactBoxType& bx)const) &ClosedSetInterface::separated);

    pybind11::class_<BoundedSetInterface>bounded_set_interface_class(module,"BoundedSetInterface",bounded_set_interface_class);
    bounded_set_interface_class.def("inside",(LowerKleenean(BoundedSetInterface::*)(const ExactBoxType& bx)const) &BoundedSetInterface::inside);
    bounded_set_interface_class.def("bounding_box", (UpperBoxType(BoundedSetInterface::*)()const)&BoundedSetInterface::bounding_box);

    pybind11::class_<CompactSetInterface, CompactSetWrapper, ClosedSetInterface, BoundedSetInterface> compact_set_interface_class(module,"CompactSetInterface", pybind11::multiple_inheritance());

    pybind11::class_<LocatedSetInterface, OvertSetInterface,CompactSetInterface> located_set_interface_class(module,"LocatedSetInterface", pybind11::multiple_inheritance());
    pybind11::class_<RegularSetInterface, OpenSetInterface,ClosedSetInterface> regular_set_interface_class(module,"RegularSetInterface", pybind11::multiple_inheritance());
}


template<class PT> PT point_from_python(pybind11::list pylst) {
    typedef typename PT::ValueType X;
    std::vector<X> lst=pybind11::cast<std::vector<X>>(pylst);
    Array<X> ary(lst.begin(),lst.end());
    return PT(Vector<X>(ary));
}

template<class X> Void export_point(pybind11::module& module, std::string name)
{
    pybind11::class_<Point<X>, DrawableInterface> point_class(module,name.c_str());
    point_class.def(pybind11::init(&point_from_python<Point<X>>));
    point_class.def(pybind11::init<Point<X>>());
    point_class.def(pybind11::init<Nat>());
    point_class.def("__getitem__", &__getitem__<Point<X>,Int,X>);
    point_class.def("__str__", &__cstr__<Point<X>>);
}

Void export_points(pybind11::module& module) {
    export_point<FloatDPValue>(module,"ExactPoint");
}

template<class IVL> Void export_interval(pybind11::module& module, std::string name) {
    typedef IVL IntervalType;
    typedef typename IntervalType::LowerBoundType LowerBoundType;
    typedef typename IntervalType::UpperBoundType UpperBoundType;
    typedef typename IntervalType::MidpointType MidpointType;

    typedef decltype(contains(declval<IntervalType>(),declval<MidpointType>())) ContainsType;
    typedef decltype(disjoint(declval<IntervalType>(),declval<IntervalType>())) DisjointType;
    typedef decltype(subset(declval<IntervalType>(),declval<IntervalType>())) SubsetType;

    pybind11::class_< IntervalType > interval_class(module,name.c_str());
    interval_class.def(pybind11::init<IntervalType>());
    interval_class.def(pybind11::init<MidpointType>());
    interval_class.def(pybind11::init<LowerBoundType,UpperBoundType>());
    interval_class.def(pybind11::init([](pybind11::dict pydct){return interval_from_dict<IntervalType>(pydct);}));
    pybind11::implicitly_convertible<pybind11::dict, IntervalType>();

    if constexpr (IsConstructible<IntervalType,DyadicInterval>::value and not IsSame<IntervalType,DyadicInterval>::value) {
        interval_class.def(pybind11::init<DyadicInterval>());
    }
    if constexpr (IsConstructible<IntervalType,DyadicInterval>::value and not IsSame<IntervalType,DyadicInterval>::value) {
        interval_class.def(pybind11::init([](Dyadic l, Dyadic u){return IntervalType(DyadicInterval(l,u));}));
    }

    if constexpr (HasEquality<IVL,IVL>::value) {
        interval_class.def("__eq__",  &__eq__<IVL,IVL , Return<EqualityType<IVL,IVL>> >);
        interval_class.def("__ne__",  &__ne__<IVL,IVL , Return<InequalityType<IVL,IVL>> >);
    }

    interval_class.def("lower", &IntervalType::lower);
    interval_class.def("upper", &IntervalType::upper);
    interval_class.def("midpoint", &IntervalType::midpoint);
    interval_class.def("radius", &IntervalType::radius);
    interval_class.def("width", &IntervalType::width);
    interval_class.def("contains", (ContainsType(*)(IntervalType const&,MidpointType const&)) &contains);
    interval_class.def("empty", &IntervalType::is_empty);
    interval_class.def("__str__",&__cstr__<IntervalType>);
    //interval_class.def("__repr__",&__repr__<IntervalType>);

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
    using IVL=typename BX::IntervalType;

    using BoxType = BX;
    using IntervalType = IVL;

    typedef decltype(disjoint(declval<BX>(),declval<BX>())) DisjointType;
    typedef decltype(subset(declval<BX>(),declval<BX>())) SubsetType;
    typedef decltype(separated(declval<BX>(),declval<BX>())) SeparatedType;
    typedef decltype(overlap(declval<BX>(),declval<BX>())) OverlapType;
    typedef decltype(covers(declval<BX>(),declval<BX>())) CoversType;
    typedef decltype(inside(declval<BX>(),declval<BX>())) InsideType;

    //NOTE: Boxes do not inherit SetInterfaces or DrawableInterface in C++ API
    //pybind11::class_<ExactBoxType,pybind11::bases<CompactSetInterface,OpenSetInterface,DrawableInterface>>
    pybind11::class_<BoxType> box_class(module,name.c_str());
    box_class.def(pybind11::init<BoxType>());
    box_class.def(pybind11::init<DimensionType>());
    box_class.def(pybind11::init<Array<IntervalType>>());
    box_class.def(pybind11::init(&box_from_list<BoxType>));
    pybind11::implicitly_convertible<pybind11::list,BoxType>();

    if constexpr (IsConstructible<BoxType,DyadicBox>::value and not IsSame<BoxType,DyadicBox>::value) {
         box_class.def(pybind11::init<Box<Interval<Dyadic>>>());
    }

    if constexpr (HasEquality<BX,BX>::value) {
        box_class.def("__eq__",  __eq__<BX,BX , Return<EqualityType<BX,BX>> >);
        box_class.def("__ne__",  __ne__<BX,BX , Return<InequalityType<BX,BX>> >);
    }

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
}


template<> Void export_box<DyadicBox>(pybind11::module& module, std::string name)
{
    using BX=DyadicBox;
    using BoxType = BX;
    pybind11::class_<BoxType> box_class(module,name.c_str());
    box_class.def(pybind11::init(&box_from_list<BoxType>));
    box_class.def(pybind11::init<BoxType>());
    box_class.def("__str__",&__cstr__<BoxType>);

    pybind11::implicitly_convertible<pybind11::list,BoxType>();

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

Void export_zonotope(pybind11::module& module)
{
    pybind11::class_<Zonotope,pybind11::bases<CompactSetInterface,OpenSetInterface,DrawableInterface>> zonotope_class(module,"Zonotope");
    zonotope_class.def(pybind11::init<Zonotope>());
    zonotope_class.def(pybind11::init<Vector<FloatDPValue>,Matrix<FloatDPValue>,Vector<FloatDPError>>());
    zonotope_class.def(pybind11::init<Vector<FloatDPValue>,Matrix<FloatDPValue>>());
    zonotope_class.def(pybind11::init<ExactBoxType>());
    zonotope_class.def("centre",&Zonotope::centre);
    zonotope_class.def("generators",&Zonotope::generators);
    zonotope_class.def("error",&Zonotope::error);
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

//    module.def("image", (Zonotope(*)(const Zonotope&, const ValidatedVectorFunction&)) &image);
}

Void export_polytope(pybind11::module& module)
{
    pybind11::class_<Polytope,pybind11::bases<LocatedSetInterface,DrawableInterface>> polytope_class(module,"Polytope");
    polytope_class.def(pybind11::init<Polytope>());
    polytope_class.def(pybind11::init<Int>());
    polytope_class.def("new_vertex",&Polytope::new_vertex);
    polytope_class.def("__iter__",boost::python::range(&Polytope::vertices_begin,&Polytope::vertices_end));
    polytope_class.def(self_ns::str(self));
}

*/

Void export_curve(pybind11::module& module)
{
    pybind11::class_<InterpolatedCurve, DrawableInterface> interpolated_curve_class(module,"InterpolatedCurve");
    interpolated_curve_class.def(pybind11::init<InterpolatedCurve>());
    interpolated_curve_class.def(pybind11::init<FloatDPValue,ExactPoint>());
    interpolated_curve_class.def("insert", (Void(InterpolatedCurve::*)(const FloatDPValue&, const Point<FloatDPApproximation>&)) &InterpolatedCurve::insert);
    interpolated_curve_class.def("__iter__", [](InterpolatedCurve const& c){return pybind11::make_iterator(c.begin(),c.end());});
    interpolated_curve_class.def("__str__", &__cstr__<InterpolatedCurve>);


}



Void export_affine_set(pybind11::module& module)
{
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
    validated_constrained_image_set_class.def("outer_approximation", &ValidatedConstrainedImageSet::outer_approximation);
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

    //module.def("product", (ValidatedConstrainedImageSet(*)(const ValidatedConstrainedImageSet&,const ExactBoxType&)) &product);
}



Void geometry_submodule(pybind11::module& module) {
    export_drawable_interface(module);
    export_set_interface(module);

    export_points(module);
    export_intervals(module);
    export_boxes(module);
//    export_zonotope(module);
//    export_polytope(module);
    export_curve(module);

    export_affine_set(module);

    export_constraint_set(module);
    export_constrained_image_set(module);

}

