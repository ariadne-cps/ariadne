/***************************************************************************
 *            calculus_submodule.cc
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

#include<type_traits>

#include "boost_python.h"
#include "utilities.h"

#include "algebra/expansion.tpl.h"
#include "algebra/algebra.h"
#include "function/function_interface.h"
#include "function/polynomial.h"
#include "function/function.h"
#include "function/procedure.h"

#include "function/function_model.h"
#include "function/taylor_model.h"
#include "function/taylor_function.h"

using namespace boost::python;
using namespace Ariadne;

namespace Ariadne {

ValidatedNumericType evaluate(const ValidatedScalarFunctionModel& f, const Vector<ValidatedNumericType>& x) { return f(x); }
Vector<ValidatedNumericType> evaluate(const ValidatedVectorFunctionModel& f, const Vector<ValidatedNumericType>& x) { return f(x); }

ValidatedScalarFunctionModel partial_evaluate(const ValidatedScalarFunctionModel&, SizeType, const ValidatedNumericType&);
ValidatedVectorFunctionModel partial_evaluate(const ValidatedVectorFunctionModel& f, SizeType j, const ValidatedNumericType& c);

ValidatedScalarFunctionModel compose(const ValidatedScalarFunctionModel&, const ValidatedVectorFunctionModel&);
ValidatedScalarFunctionModel compose(const ValidatedScalarFunction&, const ValidatedVectorFunctionModel&);
ValidatedVectorFunctionModel compose(const ValidatedVectorFunctionModel&, const ValidatedVectorFunctionModel&);
ValidatedVectorFunctionModel compose(const ValidatedVectorFunction&, const ValidatedVectorFunctionModel&);

ValidatedVectorFunctionModel join(const ValidatedScalarFunctionModel&, const ValidatedScalarFunctionModel&);
ValidatedVectorFunctionModel join(const ValidatedScalarFunctionModel&, const ValidatedVectorFunctionModel&);
ValidatedVectorFunctionModel join(const ValidatedVectorFunctionModel&, const ValidatedScalarFunctionModel&);
ValidatedVectorFunctionModel join(const ValidatedVectorFunctionModel&, const ValidatedVectorFunctionModel&);

template<class X> OutputStream& operator<<(OutputStream& os, const Representation< ScalarFunctionModel<X> >& frepr) {
    static_cast<const ScalarFunctionInterface<X>&>(frepr.reference()).repr(os); return os;
}

template<class X> OutputStream& operator<<(OutputStream& os, const Representation< VectorFunctionModel<X> >& frepr) {
    static_cast<const VectorFunctionInterface<X>&>(frepr.reference()).write(os); return os;
}

VectorTaylorFunction __getslice__(const VectorTaylorFunction& tf, Int start, Int stop) {
    if(start<0) { start+=tf.result_size(); }
    if(stop<0) { stop+=tf.result_size(); }
    ARIADNE_ASSERT_MSG(0<=start&&start<=stop&&SizeType(stop)<=tf.result_size(),
            "result_size="<<tf.result_size()<<", start="<<start<<", stop="<<stop);
    return VectorTaylorFunction(tf.domain(),Vector<ValidatedTaylorModel>(project(tf.models(),range(start,stop))));
}


template<>
struct from_python<MultiIndex> {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<MultiIndex>()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        Void* storage = ((boost::python::converter::rvalue_from_python_storage<MultiIndex>*)data)->storage.bytes;
        boost::python::extract<boost::python::tuple> xtup(obj_ptr);
        boost::python::extract<boost::python::list> xlst(obj_ptr);
        if(xlst.check()) {
            MultiIndex a(len(xlst)); for(SizeType i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<SizeType>(xlst()[i]); }
            new (storage) MultiIndex(a);
        } else {
            MultiIndex a(len(xtup)); for(SizeType i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<SizeType>(xtup()[i]); }
            new (storage) MultiIndex(a);
        }
        data->convertible = storage;
    }
};


template<class T>
struct from_python< Expansion<T> > {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< Expansion<T> >()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyDict_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        Void* storage = ((boost::python::converter::rvalue_from_python_storage< Expansion<T> >*)data)->storage.bytes;
        Expansion<T> r(0);
        boost::python::dict dct=boost::python::extract<boost::python::dict>(obj_ptr);
        boost::python::list lst=dct.items();
        MultiIndex a;
        if(len(lst)!=0) {
            boost::python::tuple tup=boost::python::extract<boost::python::tuple>(lst[0]);
            a=boost::python::extract<MultiIndex>(tup[0]);
            r=Expansion<T>(a.size());
            r.reserve(len(lst));
        }
        for(Int i=0; i!=len(lst); ++i) {
            boost::python::tuple tup=boost::python::extract<boost::python::tuple>(lst[i]);
            MultiIndex a=boost::python::extract<MultiIndex>(tup[0]);
            T c=boost::python::extract<T>(tup[1]);
            r.append(a,c);
        }
        new (storage) Expansion<T>(r);
        //r.unique_sort();
        data->convertible = storage;
    }
};


template<class X>
struct from_python< Vector<X> >
{
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< Vector<X> >()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        list lst=boost::python::extract<list>(obj_ptr);
        Void* storage = ((converter::rvalue_from_python_storage< Vector<X> >*) data)->storage.bytes;
        Array<X> ary(len(lst),Uninitialised());
        for(SizeType i=0; i!=ary.size(); ++i) { new (&ary[i]) X(boost::python::extract<X>(lst[i])); }
        new (storage) Vector<X>(std::move(ary));
        data->convertible = storage;
    }
};


template<>
struct from_python<VectorTaylorFunction> {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<VectorTaylorFunction>()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        Void* storage = ((boost::python::converter::rvalue_from_python_storage<VectorTaylorFunction>*)data)->storage.bytes;
        list lst=boost::python::extract<boost::python::list>(obj_ptr);
        VectorTaylorFunction* tf_ptr = new (storage) VectorTaylorFunction(len(lst));
        for(SizeType i=0; i!=tf_ptr->result_size(); ++i) { tf_ptr->set(i,boost::python::extract<ScalarTaylorFunction>(lst[i])); }
        data->convertible = storage;
    }
};


template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< X >& repr) {
    return os << repr.reference();
}

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Expansion<X> >& repr) {
    const Expansion<X>& exp=repr.reference();
    for(typename Expansion<X>::ConstIterator iter=exp.begin(); iter!=exp.end(); ++iter) {
        os << (iter==exp.begin()?'{':',') << "(";
        for(SizeType j=0; j!=iter->key().size(); ++j) {
            if(j!=0) { os << ','; } os << Int(iter->key()[j]);
        }
        os << "):" << python_representation(iter->data());
    }
    os << "}";
    return os;
}

template<class X, class CMP> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< SortedExpansion<X,CMP> >& repr) {
    return os << python_representation(static_cast<const Expansion<X>&>(repr.reference()));
}

template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<RawFloat64> >&);
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<Float64Approximation> >&);
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<Float64Bounds> >&);
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<Float64Value> >&);

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Vector<X> >& repr) {
    const Vector<X>& vec=repr.reference();
    os << "[";
    for(SizeType i=0; i!=vec.size(); ++i) {
        if(i!=0) { os << ","; }
        os << python_representation(vec[i]);
    }
    os << "]";
    return os;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation< ExactBoxType >& bx) {
    return os << PythonRepresentation< Vector<ExactIntervalType> >(bx.reference()); }

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Sweeper>& repr) {
    const Sweeper& swp=repr.reference();
    const SweeperInterface* swp_ptr = &static_cast<const SweeperInterface&>(swp);
    if(dynamic_cast<const ThresholdSweeper*>(swp_ptr)) {
        os << "ThresholdSweeper(" << dynamic_cast<const ThresholdSweeper*>(swp_ptr)->sweep_threshold() << ")";
    } else {
        os << swp;
    }
    return os;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ScalarTaylorFunction>& repr) {
    const ScalarTaylorFunction& stf=repr.reference();
    os << std::setprecision(17);
    os << "ScalarTaylorFunction"
       << "(" << python_representation(stf.domain())
       << "," << python_representation(stf.expansion())
       << "," << python_representation(stf.error())
       << "," << python_representation(stf.sweeper())
       << ")";
    return os;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<VectorTaylorFunction>& repr) {
    const VectorTaylorFunction& vtf=repr.reference();
    os << std::setprecision(17);
    os << "VectorTaylorFunction"
       << "(" << python_representation(vtf.domain())
       << "," << python_representation(vtf.expansions())
       << "," << python_representation(vtf.errors())
       << "," << python_representation(vtf.sweeper())
       << ")";
    return os;
}

List<MultiIndex> keys(const ValidatedTaylorModel& tm) {
    List<MultiIndex> r;
    for(ValidatedTaylorModel::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        r.append(iter->key());
    }
    return r;
}

ValidatedScalarFunction unrestrict(const ValidatedScalarFunctionModel& fm) {
    return ValidatedScalarFunction(fm.raw_pointer()->_clone());
}

ValidatedVectorFunction unrestrict(const ValidatedVectorFunctionModel& fm) {
    return ValidatedVectorFunction(fm.raw_pointer()->_clone());
}


ExactIntervalType _range1(const ValidatedTaylorModel&);
ExactIntervalType _range2(const ValidatedTaylorModel&);
ExactIntervalType _range3(const ValidatedTaylorModel&);

} // namespace Ariadne

Sweeper make_threshold_sweeper(double x) { return new ThresholdSweeper(x); }
Sweeper make_graded_sweeper(SizeType n) { return new GradedSweeper(n); }

Void export_expansion()
{
    from_python< Expansion<Float64Approximation> >();
    from_python< Expansion<Float64Bounds> >();
    from_python< Vector< Expansion<Float64Approximation> > >();

    class_< ExpansionValue<Float64Approximation> > expansion_value_class("ExpansionValue", init<MultiIndex,Float64Approximation>());
    // TODO: Add get/set for data
    // TODO: Use property for key
    //expansion_value_class.add_property("key", (MultiIndex const&(ExpansionValue<Float64Approximation>::*)()const)&ExpansionValue<Float64Approximation>::key);
    expansion_value_class.def("key", (const MultiIndex&(ExpansionValue<Float64Approximation>::*)()const)&ExpansionValue<Float64Approximation>::key, return_value_policy<copy_const_reference>());
    expansion_value_class.def(self_ns::str(self));

}


Void export_sweeper()
{
    class_<Sweeper> sweeper_class("Sweeper", init<Sweeper>());
    def("ThresholdSweeper", &make_threshold_sweeper );
    def("GradedSweeper", &make_graded_sweeper );
    sweeper_class.def(self_ns::str(self));

}

Expansion<Float64Value>const& get_expansion(ValidatedTaylorModel const& tm) { return tm.expansion(); }

Void export_validated_taylor_model()
{
    typedef SizeType SizeType;
    typedef ValidatedTaylorModel ValidatedTaylorModel;
    typedef VectorTaylorFunction VectorTaylorFunction;

    class_<ValidatedTaylorModel> taylor_model_class("ValidatedTaylorModel", init<ValidatedTaylorModel>());
    taylor_model_class.def( init< SizeType,Sweeper >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ValidatedTaylorModel&))&keys);
    taylor_model_class.def("value", (const Float64Value&(ValidatedTaylorModel::*)()const) &ValidatedTaylorModel::value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("gradient", (const Float64Value&(ValidatedTaylorModel::*)(SizeType)const) &ValidatedTaylorModel::gradient_value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("error", (const Float64Error&(ValidatedTaylorModel::*)()const) &ValidatedTaylorModel::error, return_value_policy<copy_const_reference>());
    taylor_model_class.def("expansion", (const Expansion<Float64Value>&(*)(ValidatedTaylorModel const&)) &get_expansion, return_value_policy<copy_const_reference>());
    taylor_model_class.def("set_error", (Void(ValidatedTaylorModel::*)(const Float64Error&)) &ValidatedTaylorModel::set_error);
    taylor_model_class.def("argument_size", &ValidatedTaylorModel::argument_size);
    taylor_model_class.def("domain", &ValidatedTaylorModel::domain);
    taylor_model_class.def("range", &ValidatedTaylorModel::range);
    taylor_model_class.def("set_sweeper", &ValidatedTaylorModel::set_sweeper);
    taylor_model_class.def("sweeper", &ValidatedTaylorModel::sweeper);
    taylor_model_class.def("sweep", (ValidatedTaylorModel&(ValidatedTaylorModel::*)()) &ValidatedTaylorModel::sweep, return_value_policy<reference_existing_object>());
    taylor_model_class.def("__getitem__", &__getitem__<ValidatedTaylorModel,MultiIndex,Float64Value>);
    taylor_model_class.def("__setitem__",&__setitem__<ValidatedTaylorModel,MultiIndex,Float64Value>);
    taylor_model_class.def(+self);
    taylor_model_class.def(-self);
    taylor_model_class.def(self+self);
    taylor_model_class.def(self-self);
    taylor_model_class.def(self*self);
    taylor_model_class.def(self/self);
    taylor_model_class.def(ValidatedNumericType()+self);
    taylor_model_class.def(ValidatedNumericType()-self);
    taylor_model_class.def(ValidatedNumericType()*self);
    taylor_model_class.def(ValidatedNumericType()/self);
    taylor_model_class.def(self+ValidatedNumericType());
    taylor_model_class.def(self-ValidatedNumericType());
    taylor_model_class.def(self*ValidatedNumericType());
    taylor_model_class.def(self/ValidatedNumericType());
    taylor_model_class.def(self+=ValidatedNumericType());
    taylor_model_class.def(self-=ValidatedNumericType());
    taylor_model_class.def(self*=ValidatedNumericType());
    taylor_model_class.def(self/=ValidatedNumericType());
    taylor_model_class.def(self+=self);
    taylor_model_class.def(self-=self);
    taylor_model_class.def(self_ns::str(self));

    taylor_model_class.def("constant",(ValidatedTaylorModel(*)(SizeType, const ValidatedNumericType&,Sweeper))&ValidatedTaylorModel::constant);
    taylor_model_class.def("coordinate",(ValidatedTaylorModel(*)(SizeType, SizeType,Sweeper))&ValidatedTaylorModel::coordinate);

    taylor_model_class.staticmethod("constant");
    taylor_model_class.staticmethod("coordinate");

    //def("max",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&,const ValidatedTaylorModel&))&max);
    //def("min",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&,const ValidatedTaylorModel&))&min);
    //def("abs",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&abs);

    typedef AlgebraOperations<ValidatedTaylorModel> Operations;
    def("pos",&Operations::_pos);
    def("neg",&Operations::_neg);
    def("rec",&Operations::_rec);
///    def("pow",&Operations::_pow);
    def("sqrt",&Operations::_sqrt);
    def("exp",&Operations::_exp);
    def("log",&Operations::_log);
    def("sin",&Operations::_sin);
    def("cos",&Operations::_cos);
    def("tan",&Operations::_tan);
    def("atan",&Operations::_atan);

    taylor_model_class.def("range", (UpperIntervalType(ValidatedTaylorModel::*)()const) &ValidatedTaylorModel::range);

    //def("evaluate", (ValidatedNumericType(*)(const ValidatedTaylorModel&, const Vector<ValidatedNumericType>&))&evaluate);
    //def("split",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&,SizeType,SplitPart)) &split);

    from_python< Vector<ValidatedTaylorModel> >();
    to_python< Vector<ValidatedTaylorModel> >();

}

Void export_approximate_taylor_model()
{
    typedef SizeType SizeType;
    typedef ApproximateTaylorModel ApproximateTaylorModel;

    class_<ApproximateTaylorModel> taylor_model_class("ApproximateTaylorModel", init<ApproximateTaylorModel>());
    taylor_model_class.def( init< SizeType,Sweeper >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ApproximateTaylorModel&))&keys);
    taylor_model_class.def("value", (const Float64Approximation&(ApproximateTaylorModel::*)()const) &ApproximateTaylorModel::value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("gradient", (const Float64Approximation&(ApproximateTaylorModel::*)(SizeType)const) &ApproximateTaylorModel::gradient_value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("expansion", (const Expansion<Float64Approximation>&(*)(ApproximateTaylorModel const&)) &get_expansion, return_value_policy<copy_const_reference>());
    taylor_model_class.def("argument_size", &ApproximateTaylorModel::argument_size);
    taylor_model_class.def("domain", &ApproximateTaylorModel::domain);
    taylor_model_class.def("range", &ApproximateTaylorModel::range);
    taylor_model_class.def("set_sweeper", &ApproximateTaylorModel::set_sweeper);
    taylor_model_class.def("sweeper", &ApproximateTaylorModel::sweeper);
    taylor_model_class.def("sweep", (ApproximateTaylorModel&(ApproximateTaylorModel::*)()) &ApproximateTaylorModel::sweep, return_value_policy<reference_existing_object>());
    taylor_model_class.def("__getitem__", &__getitem__<ApproximateTaylorModel,MultiIndex,Float64Approximation>);
    taylor_model_class.def("__setitem__",&__setitem__<ApproximateTaylorModel,MultiIndex,Float64Approximation>);
    taylor_model_class.def(+self);
    taylor_model_class.def(-self);
    taylor_model_class.def(self+self);
    taylor_model_class.def(self-self);
    taylor_model_class.def(self*self);
    taylor_model_class.def(self/self);
    taylor_model_class.def(ApproximateNumericType()+self);
    taylor_model_class.def(ApproximateNumericType()-self);
    taylor_model_class.def(ApproximateNumericType()*self);
    taylor_model_class.def(ApproximateNumericType()/self);
    taylor_model_class.def(self+ApproximateNumericType());
    taylor_model_class.def(self-ApproximateNumericType());
    taylor_model_class.def(self*ApproximateNumericType());
    taylor_model_class.def(self/ApproximateNumericType());
    taylor_model_class.def(self+=ApproximateNumericType());
    taylor_model_class.def(self-=ApproximateNumericType());
    taylor_model_class.def(self*=ApproximateNumericType());
    taylor_model_class.def(self/=ApproximateNumericType());
    taylor_model_class.def(self+=self);
    taylor_model_class.def(self-=self);
    taylor_model_class.def(self_ns::str(self));

    taylor_model_class.def("constant",(ApproximateTaylorModel(*)(SizeType, const ApproximateNumericType&,Sweeper))&ApproximateTaylorModel::constant);
    taylor_model_class.def("coordinate",(ApproximateTaylorModel(*)(SizeType, SizeType,Sweeper))&ApproximateTaylorModel::coordinate);

    taylor_model_class.staticmethod("constant");
    taylor_model_class.staticmethod("coordinate");

    //def("max",(ApproximateTaylorModel(*)(const ApproximateTaylorModel&,const ApproximateTaylorModel&))&max);
    //def("min",(ApproximateTaylorModel(*)(const ApproximateTaylorModel&,const ApproximateTaylorModel&))&min);
    //def("abs",(ApproximateTaylorModel(*)(const ApproximateTaylorModel&))&abs);

    typedef AlgebraOperations<ApproximateTaylorModel> Operations;
    def("pos",&Operations::_pos);
    def("neg",&Operations::_neg);
    def("rec",&Operations::_rec);
///    def("pow",&Operations::_pow);
    def("sqrt",&Operations::_sqrt);
    def("exp",&Operations::_exp);
    def("log",&Operations::_log);
    def("sin",&Operations::_sin);
    def("cos",&Operations::_cos);
    def("tan",&Operations::_tan);
    def("atan",&Operations::_atan);

    //def("evaluate", (ApproximateNumericType(*)(const ApproximateTaylorModel&, const Vector<ApproximateNumericType>&))&evaluate);
    //def("split",(ApproximateTaylorModel(*)(const ApproximateTaylorModel&,SizeType,SplitPart)) &split);

    from_python< Vector<ApproximateTaylorModel> >();
    to_python< Vector<ApproximateTaylorModel> >();
}


Void export_scalar_function_model()
{
     class_<ValidatedScalarFunctionModel> scalar_function_model_class("ValidatedScalarFunctionModel",init<ValidatedScalarFunctionModel>());
    scalar_function_model_class.def(init<ScalarTaylorFunction>());
    scalar_function_model_class.def("argument_size", &ValidatedScalarFunctionModel::argument_size);
    scalar_function_model_class.def("domain", &ValidatedScalarFunctionModel::domain);
    scalar_function_model_class.def("codomain", &ValidatedScalarFunctionModel::codomain);
    scalar_function_model_class.def("range", &ValidatedScalarFunctionModel::range);
    scalar_function_model_class.def("clobber", &ValidatedScalarFunctionModel::clobber);
    scalar_function_model_class.def("error", &ValidatedScalarFunctionModel::error);
    scalar_function_model_class.def("__call__", (Float64Bounds(ValidatedScalarFunctionModel::*)(const Vector<Float64Bounds>&)const) &ValidatedScalarFunctionModel::operator());
    scalar_function_model_class.def(self+self);
    scalar_function_model_class.def(self-self);
    scalar_function_model_class.def(self*self);
    scalar_function_model_class.def(self/self);
    scalar_function_model_class.def(self+ValidatedNumericType());
    scalar_function_model_class.def(self-ValidatedNumericType());
    scalar_function_model_class.def(self*ValidatedNumericType());
    scalar_function_model_class.def(self/ValidatedNumericType());
    scalar_function_model_class.def(ValidatedNumericType()+self);
    scalar_function_model_class.def(ValidatedNumericType()-self);
    scalar_function_model_class.def(ValidatedNumericType()*self);
    scalar_function_model_class.def(ValidatedNumericType()/self);
    scalar_function_model_class.def("__str__", &__cstr__<ValidatedScalarFunctionModel>);
    scalar_function_model_class.def("__repr__", &__crepr__<ValidatedScalarFunctionModel>);
    //scalar_function_model_class.def("__repr__",&__repr__<ValidatedScalarFunctionModel>);

    def("evaluate", (ValidatedNumericType(*)(const ValidatedScalarFunctionModel&,const Vector<ValidatedNumericType>&)) &evaluate);
    def("partial_evaluate", (ValidatedScalarFunctionModel(*)(const ValidatedScalarFunctionModel&,SizeType,const ValidatedNumericType&)) &partial_evaluate);

    def("compose", (ValidatedScalarFunctionModel(*)(const ValidatedScalarFunctionModel&, const ValidatedVectorFunctionModel&)) &compose);
    def("compose", (ValidatedScalarFunctionModel(*)(const ValidatedScalarFunction&, const ValidatedVectorFunctionModel&)) &compose);

    def("unrestrict", (ValidatedScalarFunction(*)(const ValidatedScalarFunctionModel&)) &unrestrict);


}

Void export_vector_function_model()
{
    class_<ValidatedVectorFunctionModel> vector_function_model_class("ValidatedVectorFunctionModel",init<ValidatedVectorFunctionModel>());
    vector_function_model_class.def(init<VectorTaylorFunction>());
    vector_function_model_class.def("result_size", &ValidatedVectorFunctionModel::result_size);
    vector_function_model_class.def("argument_size", &ValidatedVectorFunctionModel::argument_size);
    vector_function_model_class.def("domain", &ValidatedVectorFunctionModel::domain);
    vector_function_model_class.def("codomain", &ValidatedVectorFunctionModel::codomain);
    vector_function_model_class.def("range", &ValidatedVectorFunctionModel::range);
    //vector_function_model_class.def("__getslice__", (VectorTaylorFunction(*)(const VectorTaylorFunction&,Int,Int))&__getslice__);
    vector_function_model_class.def("__getitem__", &__getitem__<ValidatedVectorFunctionModel,SizeType,ValidatedScalarFunctionModel>);
    vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorFunctionModel,SizeType,ValidatedScalarFunctionModel>);
    //vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorFunctionModel,SizeType,ValidatedScalarFunction>);
    vector_function_model_class.def("__call__", (Vector<Float64Bounds>(ValidatedVectorFunctionModel::*)(const Vector<Float64Bounds>&)const) &ValidatedVectorFunctionModel::operator());
    vector_function_model_class.def(self*ValidatedNumericType());
    vector_function_model_class.def("__str__", &__cstr__<ValidatedVectorFunctionModel>);
    vector_function_model_class.def("__repr__", &__crepr__<ValidatedVectorFunctionModel>);
    //export_vector_function_model.def("__repr__",&__repr__<ValidatedVectorFunctionModel>);

    def("evaluate", (Vector<ValidatedNumericType>(*)(const ValidatedVectorFunctionModel&,const Vector<ValidatedNumericType>&)) &evaluate);

    def("compose", (ValidatedVectorFunctionModel(*)(const ValidatedVectorFunctionModel&,const ValidatedVectorFunctionModel&)) &compose);
    def("compose", (ValidatedVectorFunctionModel(*)(const ValidatedVectorFunction&,const ValidatedVectorFunctionModel&)) &compose);

    def("unrestrict", (ValidatedVectorFunction(*)(const ValidatedVectorFunctionModel&)) &unrestrict);

    def("join", (ValidatedVectorFunctionModel(*)(const ValidatedScalarFunctionModel&,const ValidatedScalarFunctionModel&)) &join);
    def("join", (ValidatedVectorFunctionModel(*)(const ValidatedScalarFunctionModel&,const ValidatedVectorFunctionModel&)) &join);
    def("join", (ValidatedVectorFunctionModel(*)(const ValidatedVectorFunctionModel&,const ValidatedScalarFunctionModel&)) &join);
    def("join", (ValidatedVectorFunctionModel(*)(const ValidatedVectorFunctionModel&,const ValidatedVectorFunctionModel&)) &join);

    to_python< List<ValidatedVectorFunctionModel> >();
}



Void export_scalar_taylor_function()
{
    class_<ScalarTaylorFunction> scalar_taylor_function_class("ScalarTaylorFunction",init<ScalarTaylorFunction>());
    scalar_taylor_function_class.def(init<ExactBoxType,ValidatedTaylorModel>());
    scalar_taylor_function_class.def(init< ExactBoxType,Sweeper >());
    scalar_taylor_function_class.def(init< ExactBoxType, const EffectiveScalarFunction&,Sweeper >());
    scalar_taylor_function_class.def(init< ExactBoxType, Expansion<Float64Value>, Float64Error, Sweeper >());
    scalar_taylor_function_class.def("error", (const Float64Error&(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::error, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("set_error", (Void(ScalarTaylorFunction::*)(const Float64Error&)) &ScalarTaylorFunction::set_error);
    scalar_taylor_function_class.def("argument_size", &ScalarTaylorFunction::argument_size);
    scalar_taylor_function_class.def("domain", &ScalarTaylorFunction::domain);
    scalar_taylor_function_class.def("codomain", &ScalarTaylorFunction::codomain);
    scalar_taylor_function_class.def("range", &ScalarTaylorFunction::range);
    scalar_taylor_function_class.def("model", (const ValidatedTaylorModel&(ScalarTaylorFunction::*)()const)&ScalarTaylorFunction::model, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("polynomial", (Polynomial<ExactIntervalType>(ScalarTaylorFunction::*)()const)&ScalarTaylorFunction::polynomial);
    scalar_taylor_function_class.def("number_of_nonzeros", (SizeType(ScalarTaylorFunction::*)()const)&ScalarTaylorFunction::number_of_nonzeros);
    scalar_taylor_function_class.def("set_sweeper", &ScalarTaylorFunction::set_sweeper);
    scalar_taylor_function_class.def("sweeper", &ScalarTaylorFunction::sweeper);
    scalar_taylor_function_class.def("sweep", (ScalarTaylorFunction&(ScalarTaylorFunction::*)()) &ScalarTaylorFunction::sweep, return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("__getitem__", &__getitem__<ScalarTaylorFunction,MultiIndex,Float64Value>);
    scalar_taylor_function_class.def("__setitem__",&__setitem__<ScalarTaylorFunction,MultiIndex,Float64Value>);
    scalar_taylor_function_class.def(+self);
    scalar_taylor_function_class.def(-self);
    scalar_taylor_function_class.def(self+self);
    scalar_taylor_function_class.def(self-self);
    scalar_taylor_function_class.def(self*self);
    scalar_taylor_function_class.def(self/self);
    scalar_taylor_function_class.def(self+ValidatedNumericType());
    scalar_taylor_function_class.def(self-ValidatedNumericType());
    scalar_taylor_function_class.def(self*ValidatedNumericType());
    scalar_taylor_function_class.def(self/ValidatedNumericType());
    scalar_taylor_function_class.def(ValidatedNumericType()+self);
    scalar_taylor_function_class.def(ValidatedNumericType()-self);
    scalar_taylor_function_class.def(ValidatedNumericType()*self);
    scalar_taylor_function_class.def(ValidatedNumericType()/self);
    scalar_taylor_function_class.def(self+=ValidatedNumericType());
    scalar_taylor_function_class.def(self-=ValidatedNumericType());
    scalar_taylor_function_class.def(self*=ValidatedNumericType());
    scalar_taylor_function_class.def(self/=ValidatedNumericType());
    scalar_taylor_function_class.def(self+ValidatedScalarFunction());
    scalar_taylor_function_class.def(self-ValidatedScalarFunction());
    scalar_taylor_function_class.def(self*ValidatedScalarFunction());
    scalar_taylor_function_class.def(self/ValidatedScalarFunction());
    scalar_taylor_function_class.def(ValidatedScalarFunction()+self);
    scalar_taylor_function_class.def(ValidatedScalarFunction()-self);
    scalar_taylor_function_class.def(ValidatedScalarFunction()*self);
    scalar_taylor_function_class.def(ValidatedScalarFunction()/self);
    scalar_taylor_function_class.def(self+=self);
    scalar_taylor_function_class.def(self-=self);
    scalar_taylor_function_class.def("__str__", &__cstr__<ScalarTaylorFunction>);
    scalar_taylor_function_class.def("__repr__", &__crepr__<ScalarTaylorFunction>);
    scalar_taylor_function_class.def("__mul__",&__mul__< VectorTaylorFunction, ScalarTaylorFunction, Vector<ValidatedNumericType> >);

    //scalar_taylor_function_class.def("__str__",(StringType(*)(const ScalarTaylorFunction&)) &__str__);
    //scalar_taylor_function_class.def("__cstr__",(StringType(*)(const ScalarTaylorFunction&)) &__cstr__);
    //scalar_taylor_function_class.def("__repr__",(StringType(*)(const ScalarTaylorFunction&)) &__repr__);
    scalar_taylor_function_class.def("value", (const Float64Value&(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::value,return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("sweep", (ScalarTaylorFunction&(ScalarTaylorFunction::*)())&ScalarTaylorFunction::sweep,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("clobber", (ScalarTaylorFunction&(ScalarTaylorFunction::*)()) &ScalarTaylorFunction::clobber,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("set_sweeper",&ScalarTaylorFunction::set_sweeper);
    scalar_taylor_function_class.def("sweeper",&ScalarTaylorFunction::sweeper);
    //scalar_taylor_function_class.def("set_default_maximum_degree",&ScalarTaylorFunction::set_default_maximum_degree);
    //scalar_taylor_function_class.def("set_default_sweep_threshold",&ScalarTaylorFunction::set_default_sweep_threshold);
    //scalar_taylor_function_class.def("default_maximum_degree",&ScalarTaylorFunction::default_maximum_degree);
    //scalar_taylor_function_class.def("default_sweep_threshold",&ScalarTaylorFunction::default_sweep_threshold);
    scalar_taylor_function_class.def("__call__", (Float64Approximation(ScalarTaylorFunction::*)(const Vector<Float64Approximation>&)const) &ScalarTaylorFunction::operator());
    scalar_taylor_function_class.def("__call__", (Float64Bounds(ScalarTaylorFunction::*)(const Vector<Float64Bounds>&)const) &ScalarTaylorFunction::operator());
    //scalar_taylor_function_class.def("gradient", (Covector<Float64Bounds>(ScalarTaylorFunction::*)(const Vector<Float64Bounds>&)const) &ScalarTaylorFunction::gradient);
    scalar_taylor_function_class.def("function", (EffectiveScalarFunction(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::function);
    scalar_taylor_function_class.def("polynomial", (Polynomial<Float64Bounds>(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::polynomial);
    scalar_taylor_function_class.def("set", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,SizeType j, const Float64Bounds&)) &partial_evaluate);
    scalar_taylor_function_class.def("restriction", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ExactBoxType&)) &restriction);
    scalar_taylor_function_class.def("extension", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ExactBoxType&)) &extension);
    //scalar_taylor_function_class.staticmethod("set_default_maximum_degree");
    //scalar_taylor_function_class.staticmethod("set_default_sweep_threshold");
    //scalar_taylor_function_class.staticmethod("default_maximum_degree");
    //scalar_taylor_function_class.staticmethod("default_sweep_threshold");



    scalar_taylor_function_class.def("zero",(ScalarTaylorFunction(*)(const ExactBoxType&,Sweeper))&ScalarTaylorFunction::zero);
    scalar_taylor_function_class.def("constant",(ScalarTaylorFunction(*)(const ExactBoxType&,const ValidatedNumericType&,Sweeper))&ScalarTaylorFunction::constant);
    scalar_taylor_function_class.def("coordinate",(ScalarTaylorFunction(*)(const ExactBoxType&,SizeType,Sweeper))&ScalarTaylorFunction::coordinate);


    scalar_taylor_function_class.staticmethod("constant");
    scalar_taylor_function_class.staticmethod("coordinate");

    def("split", (Pair<ScalarTaylorFunction,ScalarTaylorFunction>(*)(const ScalarTaylorFunction&,SizeType)) &Ariadne::split);
    def("evaluate",(Float64Bounds(*)(const ScalarTaylorFunction&,const Vector<Float64Bounds>&)) &evaluate);
    def("partial_evaluate",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,SizeType,const Float64Bounds&)) &partial_evaluate);

    def("midpoint",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&)) &midpoint);
    def("evaluate",(Float64Bounds(*)(const ScalarTaylorFunction&,const Vector<Float64Bounds>&)) &evaluate);
    def("derivative",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,SizeType)) &derivative);
    def("antiderivative",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,SizeType)) &antiderivative);
    def("antiderivative",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,SizeType,const Float64Bounds&)) &antiderivative);

    def("inconsistent",(Bool(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &inconsistent);
    def("refines",(Bool(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &refines);
    def("refinement",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &refinement);

    def("restriction",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ExactBoxType&)) &restriction);

//    def("embed",(ScalarTaylorFunction(*)(const ExactBoxType&,const ScalarTaylorFunction&,const ExactBoxType)) &embed);

    def("max",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&))&max<>);
    def("min",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&))&min<>);
    def("abs",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&abs<>);

/*
    typedef ScalarTaylorFunction SF; typedef SF const& SFcr;
    SF neg(SFcr); SF rec(SFcr); SF sqr(SFcr); SF pow(SFcr,Int);
    SF sqrt(SFcr); SF exp(SFcr); SF log(SFcr); SF atan(SFcr);
    SF sin(SFcr); SF cos(SFcr); SF tan(SFcr);


    def("neg",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&neg);
    def("rec",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&rec);
    def("sqr",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&sqr);
    def("pow",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&, Int))&pow);

    def("sqrt", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&sqrt);
    def("exp", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&exp);
    def("log", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&log);
    def("sin", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&sin);
    def("cos", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&cos);
    def("tan", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&tan);
    def("atan", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&atan);
*/

    to_python< Vector<ScalarTaylorFunction> >();
}


Void export_vector_taylor_function()
{
    typedef SizeType SizeType;

    class_<VectorTaylorFunction> vector_taylor_function_class("VectorTaylorFunction", init<VectorTaylorFunction>());
    vector_taylor_function_class.def( init< SizeType, ExactBoxType, Sweeper >());
    vector_taylor_function_class.def( init< ExactBoxType,const EffectiveVectorFunction&,Sweeper >());
    vector_taylor_function_class.def(init< ExactBoxType, Vector< Expansion<Float64Value> >, Vector<Float64Error>, Sweeper >());
    vector_taylor_function_class.def( init< Vector<ScalarTaylorFunction> >());
    vector_taylor_function_class.def("__len__", &VectorTaylorFunction::result_size);
    vector_taylor_function_class.def("result_size", &VectorTaylorFunction::result_size);
    vector_taylor_function_class.def("argument_size", &VectorTaylorFunction::argument_size);
    vector_taylor_function_class.def("domain", &VectorTaylorFunction::domain);
    vector_taylor_function_class.def("codomain", &VectorTaylorFunction::codomain);
    // FIXME: Omitted since const and non-const versions
    // vector_taylor_function_class.def("models", &VectorTaylorFunction::models, return_value_policy<copy_const_reference>());
    vector_taylor_function_class.def("centre", &VectorTaylorFunction::centre);
    vector_taylor_function_class.def("range", &VectorTaylorFunction::range);
    vector_taylor_function_class.def("errors", &VectorTaylorFunction::errors);
    vector_taylor_function_class.def("sweep", (VectorTaylorFunction&(VectorTaylorFunction::*)())&VectorTaylorFunction::sweep,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("clobber", (VectorTaylorFunction&(VectorTaylorFunction::*)()) &VectorTaylorFunction::clobber,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("set_sweeper",&VectorTaylorFunction::set_sweeper);
    vector_taylor_function_class.def("sweeper",&VectorTaylorFunction::sweeper);
    vector_taylor_function_class.def("sweep", (VectorTaylorFunction&(VectorTaylorFunction::*)()) &VectorTaylorFunction::sweep, return_value_policy<reference_existing_object>());
    //vector_taylor_function_class.def("__getslice__", &__getslice__<VectorTaylorFunction,SizeType,SizeType,ScalarTaylorFunction>);
    vector_taylor_function_class.def("__getslice__", (VectorTaylorFunction(*)(const VectorTaylorFunction&,Int,Int))&__getslice__);
    vector_taylor_function_class.def("__getitem__", &__getitem__<VectorTaylorFunction,SizeType,ScalarTaylorFunction>);
    vector_taylor_function_class.def("__setitem__",&__setitem__<VectorTaylorFunction,SizeType,ScalarTaylorFunction>);
    vector_taylor_function_class.def(-self);
    vector_taylor_function_class.def(self+self);
    vector_taylor_function_class.def(self-self);
    vector_taylor_function_class.def(self+Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self-Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self*ValidatedNumericType());
    vector_taylor_function_class.def(self/ValidatedNumericType());
    vector_taylor_function_class.def(self*ScalarTaylorFunction());
    vector_taylor_function_class.def(self+=Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self-=Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self*=ValidatedNumericType());
    vector_taylor_function_class.def(self/=ValidatedNumericType());
    vector_taylor_function_class.def(self+=self);
    vector_taylor_function_class.def(self-=self);
    vector_taylor_function_class.def("__str__", &__cstr__<VectorTaylorFunction>);
    vector_taylor_function_class.def("__repr__", &__crepr__<VectorTaylorFunction>);
    vector_taylor_function_class.def("clobber", (VectorTaylorFunction&(VectorTaylorFunction::*)()) &VectorTaylorFunction::clobber,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("__call__", (Vector<Float64Approximation>(VectorTaylorFunction::*)(const Vector<Float64Approximation>&)const) &VectorTaylorFunction::operator());
    vector_taylor_function_class.def("__call__", (Vector<Float64Bounds>(VectorTaylorFunction::*)(const Vector<Float64Bounds>&)const) &VectorTaylorFunction::operator());
     //vector_taylor_function_class.def("jacobian", (Vector<Float64Bounds>(VectorTaylorFunction::*)(const Vector<Float64Bounds>&)const) &VectorTaylorFunction::jacobian);
    vector_taylor_function_class.def("polynomials", (Vector< Polynomial<Float64Bounds> >(VectorTaylorFunction::*)()const) &VectorTaylorFunction::polynomials);
    vector_taylor_function_class.def("function", (EffectiveVectorFunction(VectorTaylorFunction::*)()const) &VectorTaylorFunction::function);


    vector_taylor_function_class.def("constant",(VectorTaylorFunction(*)(const ExactBoxType&, const Vector<ValidatedNumericType>&,Sweeper))&VectorTaylorFunction::constant);
    vector_taylor_function_class.def("identity",(VectorTaylorFunction(*)(const ExactBoxType&,Sweeper))&VectorTaylorFunction::identity);

    vector_taylor_function_class.staticmethod("constant");
    vector_taylor_function_class.staticmethod("identity");

    def("inconsistent", (Bool(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &inconsistent);
    def("refinement", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &refinement);
    def("refines", (Bool(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &refines);

    def("join", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &join);
    def("join", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const ScalarTaylorFunction&)) &join);
    def("join", (VectorTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &join);
    def("join", (VectorTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &join);

    def("combine", (VectorTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &combine);
    def("combine", (VectorTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &combine);
    def("combine", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const ScalarTaylorFunction&)) &combine);
    def("combine", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &combine);

    def("embed",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const ExactIntervalType&)) &embed);
    def("embed",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const ExactBoxType&)) &embed);
    def("embed",(VectorTaylorFunction(*)(const ExactBoxType&,const VectorTaylorFunction&)) &embed);

    def("restriction", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const ExactBoxType&)) &restriction);
    def("restriction", (VectorTaylorFunction(*)(const VectorTaylorFunction&,SizeType,const ExactIntervalType&)) &restriction);

    //def("split", (Pair<VectorTaylorFunction,VectorTaylorFunction>(*)(const VectorTaylorFunction&,SizeType)) &Ariadne::split);

//    def("evaluate",(Vector<Float64Approximation>(VectorTaylorFunction::*)(const Vector<Float64Approximation>&)const) &VectorTaylorFunction::evaluate);
    def("evaluate",(Vector<Float64Bounds>(*)(const VectorTaylorFunction&,const Vector<Float64Bounds>&)) &evaluate);
    def("partial_evaluate",(VectorTaylorFunction(*)(const VectorTaylorFunction&,SizeType,const Float64Bounds&)) &partial_evaluate);
    //def("compose",(ScalarTaylorFunction(*)(const RP&,const VectorTaylorFunction&)) &compose);
    def("compose",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(ScalarTaylorFunction(*)(const ValidatedScalarFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(VectorTaylorFunction(*)(const ValidatedVectorFunction&,const VectorTaylorFunction&)) &compose);
    def("derivative",(VectorTaylorFunction(*)(const VectorTaylorFunction&,SizeType)) &derivative);
    def("antiderivative",(VectorTaylorFunction(*)(const VectorTaylorFunction&,SizeType)) &antiderivative);
    def("antiderivative",(VectorTaylorFunction(*)(const VectorTaylorFunction&,SizeType,Float64Bounds)) &antiderivative);

    def("unchecked_compose",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &unchecked_compose);
    def("unchecked_compose",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &unchecked_compose);

    from_python<VectorTaylorFunction>();
}

Void calculus_submodule()
{
    export_expansion();
    export_sweeper();
//    export_approximate_taylor_model();
//    export_validated_taylor_model();
//    export_scalar_function_model();
//    export_vector_function_model();
//    export_scalar_taylor_function();
//    export_vector_taylor_function();
}


