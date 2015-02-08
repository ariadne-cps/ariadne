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

#include "boost_python.h"
#include "utilities.h"

#include <boost/python.hpp>

#include "function/function_interface.h"
#include "function/polynomial.h"
#include "function/function.h"
#include "function/function_model.h"
#include "function/taylor_function.h"

using namespace boost::python;
using namespace Ariadne;

namespace Ariadne {


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
        Expansion<T> r;
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
            T c=extract<T>(tup[1]);
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
        list lst=extract<list>(obj_ptr);
        Void* storage = ((converter::rvalue_from_python_storage< Vector<X> >*) data)->storage.bytes;
        Vector<X> res(len(lst));
        for(SizeType i=0; i!=res.size(); ++i) { res[i]=extract<X>(lst[i]); }
        new (storage) Vector<X>(res);
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
        list lst=extract<boost::python::list>(obj_ptr);
        VectorTaylorFunction* tf_ptr = new (storage) VectorTaylorFunction(len(lst));
        for(SizeType i=0; i!=tf_ptr->result_size(); ++i) { tf_ptr->set(i,extract<ScalarTaylorFunction>(lst[i])); }
        data->convertible = storage;
    }
};


OutputStream& operator<<(OutputStream& os, const PythonRepresentation<RawFloat64>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ApproximateFloat64>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ValidatedFloat64>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ExactFloat64>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<UpperFloat64>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ExactInterval>& repr);

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
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<ApproximateFloat64> >&);
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<ValidatedFloat64> >&);
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<ExactFloat64> >&);

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

OutputStream& operator<<(OutputStream& os, const PythonRepresentation< ExactBox >& bx) {
    return os << PythonRepresentation< Vector<ExactInterval> >(bx.reference()); }

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


ExactInterval _range1(const ValidatedTaylorModel&);
ExactInterval _range2(const ValidatedTaylorModel&);
ExactInterval _range3(const ValidatedTaylorModel&);

} // namespace Ariadne

Sweeper make_threshold_sweeper(double x) { return new ThresholdSweeper(x); }
Sweeper make_graded_sweeper(SizeType n) { return new GradedSweeper(n); }

Void export_expansion()
{
    from_python< Expansion<ApproximateFloat64> >();
    from_python< Expansion<ValidatedFloat64> >();
    from_python< Vector< Expansion<ApproximateFloat64> > >();

    class_< ExpansionValue<ApproximateFloat64> > expansion_value_class("ExpansionValue", init<MultiIndex,ApproximateFloat64>());
    // TODO: Add get/set for data
    // TODO: Use property for key
    //expansion_value_class.add_property("key", (MultiIndex const&(ExpansionValue<ApproximateFloat64>::*)()const)&ExpansionValue<ApproximateFloat64>::key);
    expansion_value_class.def("key", (const MultiIndex&(ExpansionValue<ApproximateFloat64>::*)()const)&ExpansionValue<ApproximateFloat64>::key, return_value_policy<copy_const_reference>());
    expansion_value_class.def(self_ns::str(self));

}


Void export_sweeper()
{
    class_<Sweeper> sweeper_class("Sweeper", init<Sweeper>());
    def("ThresholdSweeper", &make_threshold_sweeper );
    def("GradedSweeper", &make_graded_sweeper );
    sweeper_class.def(self_ns::str(self));

}

Expansion<ExactFloat64>const& get_expansion(ValidatedTaylorModel const& tm) { return tm.expansion(); }

Void export_taylor_model()
{
    typedef SizeType SizeType;
    typedef ValidatedTaylorModel ValidatedTaylorModel;
    typedef VectorTaylorFunction VectorTaylorFunction;


    class_<ValidatedTaylorModel> taylor_model_class("ValidatedTaylorModel", init<ValidatedTaylorModel>());
    taylor_model_class.def( init< SizeType,Sweeper >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ValidatedTaylorModel&))&keys);
    taylor_model_class.def("value", (const ExactFloat64&(ValidatedTaylorModel::*)()const) &ValidatedTaylorModel::value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("gradient", (const ExactFloat64&(ValidatedTaylorModel::*)(SizeType)const) &ValidatedTaylorModel::gradient_value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("error", (const ErrorFloat64&(ValidatedTaylorModel::*)()const) &ValidatedTaylorModel::error, return_value_policy<copy_const_reference>());
    taylor_model_class.def("expansion", (const Expansion<ExactFloat64>&(*)(ValidatedTaylorModel const&)) &get_expansion, return_value_policy<copy_const_reference>());
    taylor_model_class.def("set_error", (Void(ValidatedTaylorModel::*)(const ErrorFloat64&)) &ValidatedTaylorModel::set_error);
    taylor_model_class.def("argument_size", &ValidatedTaylorModel::argument_size);
    taylor_model_class.def("domain", &ValidatedTaylorModel::domain);
    taylor_model_class.def("range", &ValidatedTaylorModel::range);
    taylor_model_class.def("set_sweeper", &ValidatedTaylorModel::set_sweeper);
    taylor_model_class.def("sweeper", &ValidatedTaylorModel::sweeper);
    taylor_model_class.def("sweep", (ValidatedTaylorModel&(ValidatedTaylorModel::*)()) &ValidatedTaylorModel::sweep, return_value_policy<reference_existing_object>());
    taylor_model_class.def("__getitem__", &__getitem__<ValidatedTaylorModel,MultiIndex,ExactFloat64>);
    taylor_model_class.def("__setitem__",&__setitem__<ValidatedTaylorModel,MultiIndex,ExactFloat64>);
    taylor_model_class.def(+self);
    taylor_model_class.def(-self);
    taylor_model_class.def(self+self);
    taylor_model_class.def(self-self);
    taylor_model_class.def(self*self);
    taylor_model_class.def(self/self);
    taylor_model_class.def(ValidatedNumber()+self);
    taylor_model_class.def(ValidatedNumber()-self);
    taylor_model_class.def(ValidatedNumber()*self);
    taylor_model_class.def(ValidatedNumber()/self);
    taylor_model_class.def(self+ValidatedNumber());
    taylor_model_class.def(self-ValidatedNumber());
    taylor_model_class.def(self*ValidatedNumber());
    taylor_model_class.def(self/ValidatedNumber());
    taylor_model_class.def(self+=ValidatedNumber());
    taylor_model_class.def(self-=ValidatedNumber());
    taylor_model_class.def(self*=ValidatedNumber());
    taylor_model_class.def(self/=ValidatedNumber());
    taylor_model_class.def(self+=self);
    taylor_model_class.def(self-=self);
    taylor_model_class.def(self_ns::str(self));

    taylor_model_class.def("constant",(ValidatedTaylorModel(*)(SizeType, const ValidatedNumber&,Sweeper))&ValidatedTaylorModel::constant);
    taylor_model_class.def("coordinate",(ValidatedTaylorModel(*)(SizeType, SizeType,Sweeper))&ValidatedTaylorModel::coordinate);

    taylor_model_class.staticmethod("constant");
    taylor_model_class.staticmethod("coordinate");

    //def("max",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&,const ValidatedTaylorModel&))&max);
    //def("min",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&,const ValidatedTaylorModel&))&min);
    //def("abs",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&abs);

    def("neg",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&neg);
    def("rec",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&rec);
    def("sqr",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&sqr);
    def("pow",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&, Int))&pow);

    def("sqrt", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&sqrt);
    def("exp", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&exp);
    def("log", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&log);
    def("sin", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&sin);
    def("cos", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&cos);
    def("tan", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&tan);

    taylor_model_class.def("range", (UpperInterval(ValidatedTaylorModel::*)()const) &ValidatedTaylorModel::range);

    //def("evaluate", (ValidatedNumber(*)(const ValidatedTaylorModel&, const Vector<ValidatedNumber>&))&evaluate);
    //def("split",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&,SizeType,SplitPart)) &split);

    from_python< Vector<ValidatedTaylorModel> >();
    to_python< Vector<ValidatedTaylorModel> >();

/*
    class_< TMV > taylor_model_vector_class("TaylorModelVector");
    taylor_model_vector_class.def("__getitem__", &__getitem__<TMV,Int,ValidatedTaylorModel>);
    taylor_model_vector_class.def("__setitem__", &__setitem__<TMV,Int,ValidatedTaylorModel>);
    taylor_model_vector_class.def(self_ns::str(self));
*/
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
    scalar_function_model_class.def("__call__", (ValidatedFloat64(ValidatedScalarFunctionModel::*)(const Vector<ValidatedFloat64>&)const) &ValidatedScalarFunctionModel::operator());
    scalar_function_model_class.def(self+self);
    scalar_function_model_class.def(self-self);
    scalar_function_model_class.def(self*self);
    scalar_function_model_class.def(self/self);
    scalar_function_model_class.def(self+ValidatedNumber());
    scalar_function_model_class.def(self-ValidatedNumber());
    scalar_function_model_class.def(self*ValidatedNumber());
    scalar_function_model_class.def(self/ValidatedNumber());
    scalar_function_model_class.def(ValidatedNumber()+self);
    scalar_function_model_class.def(ValidatedNumber()-self);
    scalar_function_model_class.def(ValidatedNumber()*self);
    scalar_function_model_class.def(ValidatedNumber()/self);
    scalar_function_model_class.def("__str__", &__cstr__<ValidatedScalarFunctionModel>);
    scalar_function_model_class.def("__repr__", &__crepr__<ValidatedScalarFunctionModel>);
    //scalar_function_model_class.def("__repr__",&__repr__<ValidatedScalarFunctionModel>);

    //def("evaluate", (ValidatedNumber(*)(const ValidatedScalarFunctionModel&,const Vector<ValidatedNumber>&)) &evaluate);
    //def("partial_evaluate", (ValidatedScalarFunctionModel(*)(const ValidatedScalarFunctionModel&,SizeType,const ValidatedNumber&)) &partial_evaluate);
    def("compose", (ValidatedScalarFunctionModel(*)(const ValidatedScalarFunctionModel&,const ValidatedVectorFunctionModel&)) &compose);
    def("compose", (ValidatedScalarFunctionModel(*)(const ValidatedScalarFunction&,const ValidatedVectorFunctionModel&)) &compose);

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
    vector_function_model_class.def("__call__", (Vector<ValidatedFloat64>(ValidatedVectorFunctionModel::*)(const Vector<ValidatedFloat64>&)const) &ValidatedVectorFunctionModel::operator());
    vector_function_model_class.def(self*ValidatedNumber());
    vector_function_model_class.def("__str__", &__cstr__<ValidatedVectorFunctionModel>);
    vector_function_model_class.def("__repr__", &__crepr__<ValidatedVectorFunctionModel>);
    //export_vector_function_model.def("__repr__",&__repr__<ValidatedVectorFunctionModel>);

    //def("evaluate", (Vector<ValidatedNumber>(*)(const ValidatedVectorFunctionModel&,const Vector<ValidatedNumber>&)) &evaluate);
    def("compose", (ValidatedVectorFunctionModel(*)(const ValidatedVectorFunctionModel&,const ValidatedVectorFunctionModel&)) &compose);
    def("compose", (ValidatedVectorFunctionModel(*)(const ValidatedVectorFunction&,const ValidatedVectorFunctionModel&)) &compose);

    def("unrestrict", (ValidatedVectorFunction(*)(const ValidatedVectorFunctionModel&)) &unrestrict);

    def("join", (ValidatedVectorFunctionModel(*)(const ValidatedVectorFunctionModel&,const ValidatedScalarFunctionModel&)) &join);

    to_python< List<ValidatedVectorFunctionModel> >();

}



Void export_scalar_taylor_function()
{
    class_<ScalarTaylorFunction> scalar_taylor_function_class("ScalarTaylorFunction",init<ScalarTaylorFunction>());
    scalar_taylor_function_class.def(init<ExactBox,ValidatedTaylorModel>());
    scalar_taylor_function_class.def(init< ExactBox,Sweeper >());
    scalar_taylor_function_class.def(init< ExactBox, const EffectiveScalarFunction&,Sweeper >());
    scalar_taylor_function_class.def(init< ExactBox, Expansion<ExactFloat64>, ErrorFloat64, Sweeper >());
    scalar_taylor_function_class.def("error", (const ErrorFloat64&(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::error, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("set_error", (Void(ScalarTaylorFunction::*)(const ErrorFloat64&)) &ScalarTaylorFunction::set_error);
    scalar_taylor_function_class.def("argument_size", &ScalarTaylorFunction::argument_size);
    scalar_taylor_function_class.def("domain", &ScalarTaylorFunction::domain, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("codomain", &ScalarTaylorFunction::codomain);
    scalar_taylor_function_class.def("range", &ScalarTaylorFunction::range);
    scalar_taylor_function_class.def("model", (const ValidatedTaylorModel&(ScalarTaylorFunction::*)()const)&ScalarTaylorFunction::model, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("polynomial", (Polynomial<ExactInterval>(ScalarTaylorFunction::*)()const)&ScalarTaylorFunction::polynomial);
    scalar_taylor_function_class.def("number_of_nonzeros", (SizeType(ScalarTaylorFunction::*)()const)&ScalarTaylorFunction::number_of_nonzeros);
    scalar_taylor_function_class.def("set_sweeper", &ScalarTaylorFunction::set_sweeper);
    scalar_taylor_function_class.def("sweeper", &ScalarTaylorFunction::sweeper);
    scalar_taylor_function_class.def("sweep", (ScalarTaylorFunction&(ScalarTaylorFunction::*)()) &ScalarTaylorFunction::sweep, return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("__getitem__", &__getitem__<ScalarTaylorFunction,MultiIndex,ExactFloat64>);
    scalar_taylor_function_class.def("__setitem__",&__setitem__<ScalarTaylorFunction,MultiIndex,ExactFloat64>);
    scalar_taylor_function_class.def(+self);
    scalar_taylor_function_class.def(-self);
    scalar_taylor_function_class.def(self+self);
    scalar_taylor_function_class.def(self-self);
    scalar_taylor_function_class.def(self*self);
    scalar_taylor_function_class.def(self/self);
    scalar_taylor_function_class.def(self+ValidatedNumber());
    scalar_taylor_function_class.def(self-ValidatedNumber());
    scalar_taylor_function_class.def(self*ValidatedNumber());
    scalar_taylor_function_class.def(self/ValidatedNumber());
    scalar_taylor_function_class.def(ValidatedNumber()+self);
    scalar_taylor_function_class.def(ValidatedNumber()-self);
    scalar_taylor_function_class.def(ValidatedNumber()*self);
    scalar_taylor_function_class.def(ValidatedNumber()/self);
    scalar_taylor_function_class.def(self+=ValidatedNumber());
    scalar_taylor_function_class.def(self-=ValidatedNumber());
    scalar_taylor_function_class.def(self*=ValidatedNumber());
    scalar_taylor_function_class.def(self/=ValidatedNumber());
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
    scalar_taylor_function_class.def("__mul__",&__mul__< VectorTaylorFunction, ScalarTaylorFunction, Vector<ValidatedNumber> >);

    //scalar_taylor_function_class.def("__str__",(StringType(*)(const ScalarTaylorFunction&)) &__str__);
    //scalar_taylor_function_class.def("__cstr__",(StringType(*)(const ScalarTaylorFunction&)) &__cstr__);
    //scalar_taylor_function_class.def("__repr__",(StringType(*)(const ScalarTaylorFunction&)) &__repr__);
    scalar_taylor_function_class.def("value", (const ExactFloat64&(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::value,return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("sweep", (ScalarTaylorFunction&(ScalarTaylorFunction::*)())&ScalarTaylorFunction::sweep,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("clobber", (ScalarTaylorFunction&(ScalarTaylorFunction::*)()) &ScalarTaylorFunction::clobber,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("set_sweeper",&ScalarTaylorFunction::set_sweeper);
    scalar_taylor_function_class.def("sweeper",&ScalarTaylorFunction::sweeper);
    //scalar_taylor_function_class.def("set_default_maximum_degree",&ScalarTaylorFunction::set_default_maximum_degree);
    //scalar_taylor_function_class.def("set_default_sweep_threshold",&ScalarTaylorFunction::set_default_sweep_threshold);
    //scalar_taylor_function_class.def("default_maximum_degree",&ScalarTaylorFunction::default_maximum_degree);
    //scalar_taylor_function_class.def("default_sweep_threshold",&ScalarTaylorFunction::default_sweep_threshold);
    scalar_taylor_function_class.def("__call__", (ApproximateFloat64(ScalarTaylorFunction::*)(const Vector<ApproximateFloat64>&)const) &ScalarTaylorFunction::evaluate);
    scalar_taylor_function_class.def("__call__", (ValidatedFloat64(ScalarTaylorFunction::*)(const Vector<ValidatedFloat64>&)const) &ScalarTaylorFunction::evaluate);
    scalar_taylor_function_class.def("evaluate", (ApproximateFloat64(ScalarTaylorFunction::*)(const Vector<ApproximateFloat64>&)const) &ScalarTaylorFunction::evaluate);
    scalar_taylor_function_class.def("evaluate", (ValidatedFloat64(ScalarTaylorFunction::*)(const Vector<ValidatedFloat64>&)const) &ScalarTaylorFunction::evaluate);
    //scalar_taylor_function_class.def("gradient", (Covector<ValidatedFloat64>(ScalarTaylorFunction::*)(const Vector<ValidatedFloat64>&)const) &ScalarTaylorFunction::gradient);
    scalar_taylor_function_class.def("function", (EffectiveScalarFunction(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::function);
    scalar_taylor_function_class.def("polynomial", (Polynomial<ValidatedFloat64>(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::polynomial);
    scalar_taylor_function_class.def("set", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,SizeType j, const ValidatedFloat64&)) &partial_evaluate);
    scalar_taylor_function_class.def("restriction", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ExactBox&)) &restriction);
    scalar_taylor_function_class.def("extension", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ExactBox&)) &extension);
    //scalar_taylor_function_class.staticmethod("set_default_maximum_degree");
    //scalar_taylor_function_class.staticmethod("set_default_sweep_threshold");
    //scalar_taylor_function_class.staticmethod("default_maximum_degree");
    //scalar_taylor_function_class.staticmethod("default_sweep_threshold");



    scalar_taylor_function_class.def("zero",(ScalarTaylorFunction(*)(const ExactBox&,Sweeper))&ScalarTaylorFunction::zero);
    scalar_taylor_function_class.def("constant",(ScalarTaylorFunction(*)(const ExactBox&,const ValidatedNumber&,Sweeper))&ScalarTaylorFunction::constant);
    scalar_taylor_function_class.def("coordinate",(ScalarTaylorFunction(*)(const ExactBox&,SizeType,Sweeper))&ScalarTaylorFunction::coordinate);


    scalar_taylor_function_class.staticmethod("constant");
    scalar_taylor_function_class.staticmethod("coordinate");

    def("split", (Pair<ScalarTaylorFunction,ScalarTaylorFunction>(*)(const ScalarTaylorFunction&,SizeType)) &Ariadne::split);
    def("evaluate",(ValidatedFloat64(*)(const ScalarTaylorFunction&,const Vector<ValidatedFloat64>&)) &evaluate);
    def("partial_evaluate",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,SizeType,const ValidatedFloat64&)) &partial_evaluate);

    def("midpoint",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&)) &midpoint);
    def("evaluate",(ValidatedFloat64(*)(const ScalarTaylorFunction&,const Vector<ValidatedFloat64>&)) &evaluate);
    def("derivative",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,SizeType)) &derivative);
    def("antiderivative",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,SizeType)) &antiderivative);
    def("antiderivative",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,SizeType,ExactFloat64)) &antiderivative);

    def("inconsistent",(Bool(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &inconsistent);
    def("refines",(Bool(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &refines);
    def("refinement",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &refinement);

    def("restriction",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ExactBox&)) &restriction);

//    def("embed",(ScalarTaylorFunction(*)(const ExactBox&,const ScalarTaylorFunction&,const ExactBox)) &embed);

    def("max",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&))&max<>);
    def("min",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&))&min<>);
    def("abs",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&abs<>);

    def("neg",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&neg<>);
    def("rec",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&rec<>);
    def("sqr",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&sqr<>);
    def("pow",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&, Int))&pow<>);

    def("sqrt", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&sqrt<>);
    def("exp", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&exp<>);
    def("log", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&log<>);
    def("sin", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&sin<>);
    def("cos", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&cos<>);
    def("tan", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&tan<>);

    to_python< Vector<ScalarTaylorFunction> >();
}


Void export_vector_taylor_function()
{
    typedef SizeType SizeType;

    class_<VectorTaylorFunction> vector_taylor_function_class("VectorTaylorFunction", init<VectorTaylorFunction>());
    vector_taylor_function_class.def( init< SizeType, ExactBox, Sweeper >());
    vector_taylor_function_class.def( init< ExactBox,const EffectiveVectorFunction&,Sweeper >());
    vector_taylor_function_class.def(init< ExactBox, Vector< Expansion<ExactFloat64> >, Vector<ErrorFloat64>, Sweeper >());
    vector_taylor_function_class.def( init< Vector<ScalarTaylorFunction> >());
    vector_taylor_function_class.def("__len__", &VectorTaylorFunction::result_size);
    vector_taylor_function_class.def("result_size", &VectorTaylorFunction::result_size);
    vector_taylor_function_class.def("argument_size", &VectorTaylorFunction::argument_size);
    vector_taylor_function_class.def("domain", &VectorTaylorFunction::domain, return_value_policy<copy_const_reference>());
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
    vector_taylor_function_class.def(self+Vector<ValidatedNumber>());
    vector_taylor_function_class.def(self-Vector<ValidatedNumber>());
    vector_taylor_function_class.def(self*ValidatedNumber());
    vector_taylor_function_class.def(self/ValidatedNumber());
    vector_taylor_function_class.def(self*ScalarTaylorFunction());
    vector_taylor_function_class.def(self+=Vector<ValidatedNumber>());
    vector_taylor_function_class.def(self-=Vector<ValidatedNumber>());
    vector_taylor_function_class.def(self*=ValidatedNumber());
    vector_taylor_function_class.def(self/=ValidatedNumber());
    vector_taylor_function_class.def(self+=self);
    vector_taylor_function_class.def(self-=self);
    vector_taylor_function_class.def("__str__", &__cstr__<VectorTaylorFunction>);
    vector_taylor_function_class.def("__repr__", &__crepr__<VectorTaylorFunction>);
    vector_taylor_function_class.def("clobber", (VectorTaylorFunction&(VectorTaylorFunction::*)()) &VectorTaylorFunction::clobber,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("__call__", (Vector<ApproximateFloat64>(VectorTaylorFunction::*)(const Vector<ApproximateFloat64>&)const) &VectorTaylorFunction::evaluate);
    vector_taylor_function_class.def("__call__", (Vector<ValidatedFloat64>(VectorTaylorFunction::*)(const Vector<ValidatedFloat64>&)const) &VectorTaylorFunction::evaluate);
    vector_taylor_function_class.def("evaluate", (Vector<ApproximateFloat64>(VectorTaylorFunction::*)(const Vector<ApproximateFloat64>&)const) &VectorTaylorFunction::evaluate);
    vector_taylor_function_class.def("evaluate", (Vector<ValidatedFloat64>(VectorTaylorFunction::*)(const Vector<ValidatedFloat64>&)const) &VectorTaylorFunction::evaluate);
    //vector_taylor_function_class.def("jacobian", (Vector<ValidatedFloat64>(VectorTaylorFunction::*)(const Vector<ValidatedFloat64>&)const) &VectorTaylorFunction::jacobian);
    vector_taylor_function_class.def("polynomials", (Vector< Polynomial<ValidatedFloat64> >(VectorTaylorFunction::*)()const) &VectorTaylorFunction::polynomials);
    vector_taylor_function_class.def("function", (EffectiveVectorFunction(VectorTaylorFunction::*)()const) &VectorTaylorFunction::function);


    vector_taylor_function_class.def("constant",(VectorTaylorFunction(*)(const ExactBox&, const Vector<ValidatedNumber>&,Sweeper))&VectorTaylorFunction::constant);
    vector_taylor_function_class.def("identity",(VectorTaylorFunction(*)(const ExactBox&,Sweeper))&VectorTaylorFunction::identity);

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

    def("embed",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const ExactInterval&)) &embed);
    def("embed",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const ExactBox&)) &embed);
    def("embed",(VectorTaylorFunction(*)(const ExactBox&,const VectorTaylorFunction&)) &embed);

    def("restriction", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const ExactBox&)) &restriction);
    def("restriction", (VectorTaylorFunction(*)(const VectorTaylorFunction&,SizeType,const ExactInterval&)) &restriction);

    //def("split", (Pair<VectorTaylorFunction,VectorTaylorFunction>(*)(const VectorTaylorFunction&,SizeType)) &Ariadne::split);

    def("evaluate",(Vector<ApproximateFloat64>(VectorTaylorFunction::*)(const Vector<ApproximateFloat64>&)const) &VectorTaylorFunction::evaluate);
    def("evaluate",(Vector<ValidatedFloat64>(VectorTaylorFunction::*)(const Vector<ValidatedFloat64>&)const) &VectorTaylorFunction::evaluate);
    def("partial_evaluate",(VectorTaylorFunction(*)(const VectorTaylorFunction&,SizeType,const ValidatedFloat64&)) &partial_evaluate);
    //def("compose",(ScalarTaylorFunction(*)(const RP&,const VectorTaylorFunction&)) &compose);
    def("compose",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(ScalarTaylorFunction(*)(const ValidatedScalarFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(VectorTaylorFunction(*)(const ValidatedVectorFunction&,const VectorTaylorFunction&)) &compose);
    def("derivative",(VectorTaylorFunction(*)(const VectorTaylorFunction&,SizeType)) &derivative);
    def("antiderivative",(VectorTaylorFunction(*)(const VectorTaylorFunction&,SizeType)) &antiderivative);
    def("antiderivative",(VectorTaylorFunction(*)(const VectorTaylorFunction&,SizeType,ExactFloat64)) &antiderivative);

    def("unchecked_compose",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &unchecked_compose);
    def("unchecked_compose",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &unchecked_compose);

    from_python<VectorTaylorFunction>();
}

Void calculus_submodule()
{
    export_expansion();
    export_sweeper();
    export_taylor_model();
    export_scalar_function_model();
    export_vector_function_model();
    export_scalar_taylor_function();
    export_vector_taylor_function();
}


