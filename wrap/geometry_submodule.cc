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
struct from_python_dict<ExactInterval> {
    from_python_dict() { converter::registry::push_back(&convertible,&construct,type_id<ExactInterval>()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyDict_Check(obj_ptr) || len(boost::python::extract<boost::python::dict>(obj_ptr))!=1) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::dict dct = boost::python::extract<boost::python::dict>(obj_ptr);
        boost::python::list lst=dct.items();
        assert(boost::python::len(lst)==1);
        void* storage = ((converter::rvalue_from_python_storage<ExactInterval>*)data)->storage.bytes;
        new (storage) ExactInterval(boost::python::extract<Float>(lst[0][0]),boost::python::extract<Float>(lst[0][1]));
        data->convertible = storage;
    }
};


template<>
struct from_python_list<ExactInterval> {
    from_python_list() { converter::registry::push_back(&convertible,&construct,type_id<ExactInterval>()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) || len(boost::python::extract<boost::python::list>(obj_ptr))!=2) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        assert(boost::python::len(lst)==2);
        void* storage = ((converter::rvalue_from_python_storage<ExactInterval>*)data)->storage.bytes;
        new (storage) ExactInterval(boost::python::extract<Float>(lst[0]),boost::python::extract<Float>(lst[1]));
        data->convertible = storage;
    }
};

template<>
struct from_python<ExactPoint> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<ExactPoint>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        extract<boost::python::tuple> xtup(obj_ptr);
        extract<boost::python::list> xlst(obj_ptr);
        ExactPoint pt;
        if(xtup.check()) {
            boost::python::tuple tup=xtup(); pt=ExactPoint(len(tup));
            for(int i=0; i!=len(tup); ++i) { pt[i]=ExactFloat(extract<double>(tup[i])); }
        } else if(xlst.check()) {
            boost::python::list lst=xlst(); pt=ExactPoint(len(lst));
            for(int i=0; i!=len(lst); ++i) { pt[i]=ExactFloat(extract<double>(lst[i])); }
        }
        void* storage = ((converter::rvalue_from_python_storage<ExactInterval>*)data)->storage.bytes;
        new (storage) ExactPoint(pt);
        data->convertible = storage;
    }
};

template<>
struct from_python<ExactBox> {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<ExactBox>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((converter::rvalue_from_python_storage<ExactInterval>*)data)->storage.bytes;
        boost::python::list lst=extract<boost::python::list>(obj_ptr);
        ExactBox* bx_ptr = new (storage) ExactBox(len(lst));
        for(int i=0; i!=len(lst); ++i) { (*bx_ptr)[i]=extract<ExactInterval>(lst[i]); }
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
std::ostream& operator<<(std::ostream& os, const PythonRepresentation<ExactInterval>& x) {
    return os << PythonRepresentation<ValidatedFloat>(ValidatedFloat(x.reference()));
}



class OpenSetWrapper
  : public virtual OpenSetInterface, public wrapper< OpenSetInterface >
{
  public:
    OpenSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool covers(const ExactBox& r) const { return this->get_override("covers")(); }
    tribool overlaps(const ExactBox& r) const { return this->get_override("overlaps")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

class ClosedSetWrapper
  : public virtual ClosedSetInterface, public wrapper< ClosedSetInterface >
{
  public:
    ClosedSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool separated(const ExactBox& r) const { return this->get_override("separated")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


class OvertSetWrapper
  : public virtual OvertSetInterface, public wrapper< OvertSetInterface >
{
  public:
    OvertSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool overlaps(const ExactBox& r) const { return this->get_override("overlaps")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


class CompactSetWrapper
  : public virtual CompactSetInterface, public wrapper< CompactSetInterface >
{
  public:
    CompactSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool separated(const ExactBox& r) const { return this->get_override("separated")(); }
    tribool inside(const ExactBox& r) const { return this->get_override("inside")(); }
    tribool bounded() const { return this->get_override("bounded")(); }
    UpperBox bounding_box() const { return this->get_override("bounding_box")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

class RegularSetWrapper
  : public virtual LocatedSetInterface, public wrapper< RegularSetWrapper >
{
  public:
    RegularSetWrapper* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool overlaps(const ExactBox& r) const { return this->get_override("overlaps")(); }
    tribool covers(const ExactBox& r) const { return this->get_override("covers")(); }
    tribool separated(const ExactBox& r) const { return this->get_override("separated")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

class LocatedSetWrapper
  : public virtual LocatedSetInterface, public wrapper< LocatedSetInterface >
{
  public:
    LocatedSetInterface* clone() const { return this->get_override("clone")(); }
    uint dimension() const { return this->get_override("dimension")(); }
    tribool overlaps(const ExactBox& r) const { return this->get_override("overlaps")(); }
    tribool separated(const ExactBox& r) const { return this->get_override("separated")(); }
    tribool inside(const ExactBox& r) const { return this->get_override("inside")(); }
    tribool bounded() const { return this->get_override("bounded")(); }
    UpperBox bounding_box() const { return this->get_override("bounding_box")(); }
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
    class_<ExactPoint,bases<DrawableInterface>> point_class("ExactPoint",init<ExactPoint>());
    point_class.def(init<uint>());
    point_class.def("__getitem__", &__getitem__<ExactPoint,int,ExactFloat>);
    point_class.def(self_ns::str(self));

    from_python<ExactPoint>();
    implicitly_convertible<Vector<Float>,ExactPoint>();

}


void export_interval()
{
    using boost::python::class_;
    using boost::python::init;
    using boost::python::self;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::def;

    class_< ExactInterval > interval_class("ExactInterval");
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
    interval_class.def("lower", &ExactInterval::lower, return_value_policy<copy_const_reference>());
    interval_class.def("upper", &ExactInterval::upper, return_value_policy<copy_const_reference>());
    interval_class.def("midpoint", &ExactInterval::midpoint);
    interval_class.def("radius", &ExactInterval::radius);
    interval_class.def("width", &ExactInterval::width);
    interval_class.def("contains", (bool(*)(ExactInterval,Float)) &contains);
    interval_class.def("empty", (bool(ExactInterval::*)()const) &ExactInterval::empty);
    interval_class.def(boost::python::self_ns::str(self));

    from_python_dict<ExactInterval>();
    from_python_list<ExactInterval>();
    //from_python_str<ExactInterval>();

    def("midpoint", &ExactInterval::midpoint);
    def("radius", &ExactInterval::radius);
    def("width", &ExactInterval::width);

    def("disjoint", (bool(*)(ExactInterval,ExactInterval)) &disjoint);
    def("subset", (bool(*)(ExactInterval,ExactInterval)) &subset);
    def("intersection", (ExactInterval(*)(ExactInterval,ExactInterval)) &intersection);

    def("hull", (ExactInterval(*)(ExactInterval,ExactInterval)) &hull);
}

void export_box()
{
    typedef Vector<ExactInterval> IVector;
    class_<Vector<ExactInterval>> interval_vector_class("ExactIntervalVector");

    class_<ExactBox,bases<CompactSetInterface,OpenSetInterface,Vector<ExactInterval>,DrawableInterface > > box_class("ExactBox",init<ExactBox>());
    box_class.def(init<uint>());
    box_class.def(init< Vector<ExactInterval> >());
    box_class.def("__eq__", (bool(*)(const Vector<ExactInterval>&,const Vector<ExactInterval>&)) &operator==);
    box_class.def("dimension", (uint(ExactBox::*)()const) &ExactBox::dimension);
    box_class.def("centre", (ExactPoint(ExactBox::*)()const) &ExactBox::centre);
    box_class.def("radius", (Float(ExactBox::*)()const) &ExactBox::radius);
    box_class.def("separated", (tribool(ExactBox::*)(const ExactBox&)const) &ExactBox::separated);
    box_class.def("overlaps", (tribool(ExactBox::*)(const ExactBox&)const) &ExactBox::overlaps);
    box_class.def("covers", (tribool(ExactBox::*)(const ExactBox&)const) &ExactBox::covers);
    box_class.def("inside", (tribool(ExactBox::*)(const ExactBox&)const) &ExactBox::inside);
    box_class.def("empty", (bool(ExactBox::*)()const) &ExactBox::empty);
    box_class.def("widen", (ExactBox(ExactBox::*)()const) &ExactBox::widen);
    box_class.def("split", (std::pair<ExactBox,ExactBox>(ExactBox::*)()const) &ExactBox::split);
    box_class.def("split", (std::pair<ExactBox,ExactBox>(ExactBox::*)(uint)const) &ExactBox::split);
    box_class.def("split", (std::pair<ExactBox,ExactBox>(ExactBox::*)()const) &ExactBox::split);
    box_class.def(self_ns::str(self));

    def("split", (std::pair<ExactBox,ExactBox>(*)(const ExactBox&)) &split);
    def("disjoint", (bool(*)(const IVector&,const IVector&)) &disjoint);
    def("subset", (bool(*)(const IVector&,const IVector&)) &subset);

    def("product", (ExactBox(*)(const ExactBox&,const ExactInterval&)) &product);
    def("product", (ExactBox(*)(const ExactBox&,const ExactBox&)) &product);
    def("hull", (ExactBox(*)(const ExactBox&,const ExactBox&)) &hull);
    def("intersection", (ExactBox(*)(const ExactBox&,const ExactBox&)) &intersection);

    from_python<ExactBox>();
    to_python< std::pair<ExactBox,ExactBox> >();
    implicitly_convertible<Vector<ExactInterval>,ExactBox>();

}

/*
std::pair<Zonotope,Zonotope> split_pair(const Zonotope& z) {
    ListSet<Zonotope> split_list=split(z);
    ARIADNE_ASSERT(split_list.size()==2);
    return std::pair<Zonotope,Zonotope>(split_list[0],split_list[1]);
}

void export_zonotope()
{
    class_<Zonotope,bases<CompactSetInterface,OpenSetInterface,DrawableInterface> > zonotope_class("Zonotope",init<Zonotope>());
    zonotope_class.def(init< Vector<ExactFloat>, Matrix<ExactFloat>, Vector<ErrorFloat> >());
    zonotope_class.def(init< Vector<ExactFloat>, Matrix<ExactFloat> >());
    zonotope_class.def(init< ExactBox >());
    zonotope_class.def("centre",&Zonotope::centre,return_value_policy<copy_const_reference>());
    zonotope_class.def("generators",&Zonotope::generators,return_value_policy<copy_const_reference>());
    zonotope_class.def("error",&Zonotope::error,return_value_policy<copy_const_reference>());
    zonotope_class.def("contains",&Zonotope::contains);
    zonotope_class.def("split", (ListSet<Zonotope>(*)(const Zonotope&)) &split);
    zonotope_class.def("__str__",&__cstr__<Zonotope>);

    def("contains", (tribool(*)(const Zonotope&,const ExactPoint&)) &contains);
    def("separated", (tribool(*)(const Zonotope&,const ExactBox&)) &separated);
    def("overlaps", (tribool(*)(const Zonotope&,const ExactBox&)) &overlaps);
    def("separated", (tribool(*)(const Zonotope&,const Zonotope&)) &separated);

    def("polytope", (Polytope(*)(const Zonotope&)) &polytope);
    def("orthogonal_approximation", (Zonotope(*)(const Zonotope&)) &orthogonal_approximation);
    def("orthogonal_over_approximation", (Zonotope(*)(const Zonotope&)) &orthogonal_over_approximation);
    def("error_free_over_approximation", (Zonotope(*)(const Zonotope&)) &error_free_over_approximation);

//    def("apply", (Zonotope(*)(const ValidatedVectorFunction&, const Zonotope&)) &apply);

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
*/

void export_curve()
{
    to_python< std::pair<const ExactFloat,ExactPoint> >();

    class_<InterpolatedCurve,bases<DrawableInterface> > interpolated_curve_class("InterpolatedCurve",init<InterpolatedCurve>());
    interpolated_curve_class.def(init<ExactFloat,ExactPoint>());
    interpolated_curve_class.def("insert", (void(InterpolatedCurve::*)(const ApproximateFloat&, const Point<ApproximateFloat>&)) &InterpolatedCurve::insert);
    interpolated_curve_class.def("__iter__",boost::python::range(&InterpolatedCurve::begin,&InterpolatedCurve::end));
    interpolated_curve_class.def(self_ns::str(self));


}



void export_affine_set()
{

    class_<ValidatedAffineConstrainedImageSet,bases<CompactSetInterface,DrawableInterface> >
        affine_set_class("ValidatedAffineConstrainedImageSet",init<ValidatedAffineConstrainedImageSet>());
    affine_set_class.def(init<Vector<ExactInterval>, Matrix<ExactFloat>, Vector<ExactFloat> >());
    affine_set_class.def(init<Matrix<ExactFloat>, Vector<ExactFloat> >());
    affine_set_class.def("new_parameter_constraint", (void(ValidatedAffineConstrainedImageSet::*)(const Constraint<Affine<ValidatedFloat>,ValidatedFloat>&)) &ValidatedAffineConstrainedImageSet::new_parameter_constraint);
    affine_set_class.def("new_constraint", (void(ValidatedAffineConstrainedImageSet::*)(const Constraint<AffineModel<ValidatedFloat>,ValidatedFloat>&)) &ValidatedAffineConstrainedImageSet::new_constraint);
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
    box_set_class.def(init<ExactIntervalVector>());
    box_set_class.def(init< List<IntervalSet> >());

    def("intersection", (BoundedConstraintSet(*)(const ConstraintSet&,const BoxSet&)) &intersection);

}


void export_constrained_image_set()
{
    from_python< List<ValidatedConstraint> >();

    class_<ValidatedConstrainedImageSet,bases<CompactSetInterface,DrawableInterface> >
        constrained_image_set_class("ValidatedConstrainedImageSet",init<ValidatedConstrainedImageSet>());
    constrained_image_set_class.def(init<ExactBox>());
    constrained_image_set_class.def(init<ExactBox,EffectiveVectorFunction>());
    constrained_image_set_class.def(init<ExactBox,ValidatedVectorFunction>());
    constrained_image_set_class.def(init<ExactBox,ValidatedVectorFunction,List<ValidatedConstraint> >());
    constrained_image_set_class.def(init<ExactBox,ValidatedVectorFunctionModel>());
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

    def("product", (ValidatedConstrainedImageSet(*)(const ValidatedConstrainedImageSet&,const ExactIntervalVector&)) &product);
}




void geometry_submodule() {
    export_set_interface();
    export_point();
    export_box();
//    export_zonotope();
//    export_polytope();
    export_curve();

    export_affine_set();

    export_constraint_set();
    export_constrained_image_set();

}

