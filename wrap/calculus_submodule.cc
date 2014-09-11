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

#include "utilities.h"

#include <boost/python.hpp>

#include "function_interface.h"
#include "polynomial.h"
#include "function.h"
#include "function_model.h"
#include "taylor_function.h"

using namespace boost::python;
using namespace Ariadne;

namespace Ariadne {

ScalarTaylorFunction max(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
ScalarTaylorFunction min(const ScalarTaylorFunction& x, const ScalarTaylorFunction& y);
ScalarTaylorFunction neg(const ScalarTaylorFunction&);
ScalarTaylorFunction rec(const ScalarTaylorFunction&);
ScalarTaylorFunction sqr(const ScalarTaylorFunction&);
ScalarTaylorFunction pow(const ScalarTaylorFunction&, int);
ScalarTaylorFunction sqrt(const ScalarTaylorFunction&);
ScalarTaylorFunction exp(const ScalarTaylorFunction&);
ScalarTaylorFunction log(const ScalarTaylorFunction&);
ScalarTaylorFunction sin(const ScalarTaylorFunction&);
ScalarTaylorFunction cos(const ScalarTaylorFunction&);
ScalarTaylorFunction tan(const ScalarTaylorFunction&);

template<class X> std::ostream& operator<<(std::ostream& os, const Representation< ScalarFunctionModel<X> >& frepr) {
    static_cast<const ScalarFunctionInterface<X>&>(frepr.reference()).repr(os); return os;
}

template<class X> std::ostream& operator<<(std::ostream& os, const Representation< VectorFunctionModel<X> >& frepr) {
    static_cast<const VectorFunctionInterface<X>&>(frepr.reference()).write(os); return os;
}


VectorTaylorFunction __getslice__(const VectorTaylorFunction& tf, int start, int stop) {
    if(start<0) { start+=tf.result_size(); }
    if(stop<0) { stop+=tf.result_size(); }
    ARIADNE_ASSERT_MSG(0<=start&&start<=stop&&uint(stop)<=tf.result_size(),
            "result_size="<<tf.result_size()<<", start="<<start<<", stop="<<stop);
    return VectorTaylorFunction(tf.domain(),Vector<ValidatedTaylorModel>(project(tf.models(),range(start,stop))));
}


template<>
struct from_python<MultiIndex> {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<MultiIndex>()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((boost::python::converter::rvalue_from_python_storage<MultiIndex>*)data)->storage.bytes;
        boost::python::extract<boost::python::tuple> xtup(obj_ptr);
        boost::python::extract<boost::python::list> xlst(obj_ptr);
        if(xlst.check()) {
            MultiIndex a(len(xlst)); for(uint i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<uint>(xlst()[i]); }
            new (storage) MultiIndex(a);
        } else {
            MultiIndex a(len(xtup)); for(uint i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<uint>(xtup()[i]); }
            new (storage) MultiIndex(a);
        }
        data->convertible = storage;
    }
};


template<class T>
struct from_python< Expansion<T> > {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< Expansion<T> >()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyDict_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((boost::python::converter::rvalue_from_python_storage< Expansion<T> >*)data)->storage.bytes;
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
        for(int i=0; i!=len(lst); ++i) {
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
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        list lst=extract<list>(obj_ptr);
        void* storage = ((converter::rvalue_from_python_storage< Vector<X> >*) data)->storage.bytes;
        Vector<X> res(len(lst));
        for(uint i=0; i!=res.size(); ++i) { res[i]=extract<X>(lst[i]); }
        new (storage) Vector<X>(res);
        data->convertible = storage;
    }
};


template<>
struct from_python<VectorTaylorFunction> {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<VectorTaylorFunction>()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((boost::python::converter::rvalue_from_python_storage<VectorTaylorFunction>*)data)->storage.bytes;
        list lst=extract<boost::python::list>(obj_ptr);
        VectorTaylorFunction* tf_ptr = new (storage) VectorTaylorFunction(len(lst));
        for(uint i=0; i!=tf_ptr->result_size(); ++i) { tf_ptr->set(i,extract<ScalarTaylorFunction>(lst[i])); }
        data->convertible = storage;
    }
};


std::ostream& operator<<(std::ostream& os, const PythonRepresentation<Float>& repr);

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<Interval>& repr);

template<class X> std::ostream& operator<<(std::ostream& os, const PythonRepresentation< Expansion<X> >& repr) {
    const Expansion<X>& exp=repr.reference();
    for(typename Expansion<X>::const_iterator iter=exp.begin(); iter!=exp.end(); ++iter) {
        os << (iter==exp.begin()?'{':',') << "(";
        for(uint j=0; j!=iter->key().size(); ++j) {
            if(j!=0) { os << ','; } os << int(iter->key()[j]);
        }
        os << "):" << python_representation(iter->data());
    }
    os << "}";
    return os;
}

template std::ostream& operator<<(std::ostream&, const PythonRepresentation< Expansion<Float> >&);
template std::ostream& operator<<(std::ostream&, const PythonRepresentation< Expansion<Interval> >&);

template<class X> std::ostream& operator<<(std::ostream& os, const PythonRepresentation< Vector<X> >& repr) {
    const Vector<X>& vec=repr.reference();
    os << "[";
    for(uint i=0; i!=vec.size(); ++i) {
        if(i!=0) { os << ","; }
        os << python_representation(vec[i]);
    }
    os << "]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const PythonRepresentation< Box >& bx) {
    return os << PythonRepresentation< Vector<Interval> >(bx.reference()); }

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<Sweeper>& repr) {
    const Sweeper& swp=repr.reference();
    const SweeperInterface* swp_ptr = &static_cast<const SweeperInterface&>(swp);
    if(dynamic_cast<const ThresholdSweeper*>(swp_ptr)) {
        os << "ThresholdSweeper(" << dynamic_cast<const ThresholdSweeper*>(swp_ptr)->sweep_threshold() << ")";
    } else {
        os << swp;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<ScalarTaylorFunction>& repr) {
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

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<VectorTaylorFunction>& repr) {
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
    for(ValidatedTaylorModel::const_iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
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


Interval _range1(const ValidatedTaylorModel&);
Interval _range2(const ValidatedTaylorModel&);
Interval _range3(const ValidatedTaylorModel&);

} // namespace Ariadne

Sweeper make_threshold_sweeper(double x) { return new ThresholdSweeper(x); }
Sweeper make_graded_sweeper(uint n) { return new GradedSweeper(n); }

void export_expansion()
{
    from_python< Expansion<Float> >();
    from_python< Expansion<Interval> >();
    from_python< Vector< Expansion<Float> > >();

    class_< ExpansionValue<Float> > expansion_value_class("ExpansionValue", init<MultiIndex,Float>());
    // TODO: Add get/set for data
    // TODO: Use property for key
    //expansion_value_class.add_property("key", (MultiIndex const&(ExpansionValue<Float>::*)()const)&ExpansionValue<Float>::key);
    expansion_value_class.def("key", (const MultiIndex&(ExpansionValue<Float>::*)()const)&ExpansionValue<Float>::key, return_value_policy<copy_const_reference>());
    expansion_value_class.def(self_ns::str(self));

}


void export_sweeper()
{
    class_<Sweeper> sweeper_class("Sweeper", init<Sweeper>());
    def("ThresholdSweeper", &make_threshold_sweeper );
    def("GradedSweeper", &make_graded_sweeper );
    sweeper_class.def(self_ns::str(self));

}

void export_taylor_model()
{
    typedef uint Nat;
    typedef ValidatedTaylorModel ValidatedTaylorModel;
    typedef VectorTaylorFunction VectorTaylorFunction;


    class_<ValidatedTaylorModel> taylor_model_class("ValidatedTaylorModel", init<ValidatedTaylorModel>());
    taylor_model_class.def( init< Nat,Sweeper >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ValidatedTaylorModel&))&keys);
    taylor_model_class.def("value", (const Float&(ValidatedTaylorModel::*)()const) &ValidatedTaylorModel::value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("gradient", (const Float&(ValidatedTaylorModel::*)(uint)const) &ValidatedTaylorModel::gradient, return_value_policy<copy_const_reference>());
    taylor_model_class.def("error", (const Float&(ValidatedTaylorModel::*)()const) &ValidatedTaylorModel::error, return_value_policy<copy_const_reference>());
    taylor_model_class.def("expansion", (const Expansion<Float>&(ValidatedTaylorModel::*)()const) &ValidatedTaylorModel::expansion, return_value_policy<copy_const_reference>());
    taylor_model_class.def("set_error", (void(ValidatedTaylorModel::*)(const Float&)) &ValidatedTaylorModel::set_error);
    taylor_model_class.def("argument_size", &ValidatedTaylorModel::argument_size);
    taylor_model_class.def("domain", &ValidatedTaylorModel::domain);
    taylor_model_class.def("range", &ValidatedTaylorModel::range);
    taylor_model_class.def("set_sweeper", &ValidatedTaylorModel::set_sweeper);
    taylor_model_class.def("sweeper", &ValidatedTaylorModel::sweeper);
    taylor_model_class.def("sweep", (ValidatedTaylorModel&(ValidatedTaylorModel::*)()) &ValidatedTaylorModel::sweep, return_value_policy<reference_existing_object>());
    taylor_model_class.def("__getitem__", &__getitem__<ValidatedTaylorModel,MultiIndex,Float>);
    taylor_model_class.def("__setitem__",&__setitem__<ValidatedTaylorModel,MultiIndex,Float>);
    taylor_model_class.def("__setitem__",&__setitem__<ValidatedTaylorModel,MultiIndex,double>);
    taylor_model_class.def("__iter__", iterator<ValidatedTaylorModel>()); // TODO: ValidatedTaylorModel iterator does not fit in general Python iterator scheme
    taylor_model_class.def(+self);
    taylor_model_class.def(-self);
    taylor_model_class.def(self+self);
    taylor_model_class.def(self-self);
    taylor_model_class.def(self*self);
    taylor_model_class.def(self/self);
    taylor_model_class.def(self+Float());
    taylor_model_class.def(self-Float());
    taylor_model_class.def(self*Float());
    taylor_model_class.def(self/Float());
    taylor_model_class.def(Float()+self);
    taylor_model_class.def(Float()-self);
    taylor_model_class.def(Float()*self);
    taylor_model_class.def(Float()/self);
    taylor_model_class.def(self+=Float());
    taylor_model_class.def(self-=Float());
    taylor_model_class.def(self*=Float());
    taylor_model_class.def(self/=Float());
    taylor_model_class.def(Interval()+self);
    taylor_model_class.def(Interval()-self);
    taylor_model_class.def(Interval()*self);
    taylor_model_class.def(Interval()/self);
    taylor_model_class.def(self+Interval());
    taylor_model_class.def(self-Interval());
    taylor_model_class.def(self*Interval());
    taylor_model_class.def(self/Interval());
    taylor_model_class.def(self+=Interval());
    taylor_model_class.def(self-=Interval());
    taylor_model_class.def(self*=Interval());
    taylor_model_class.def(self/=Interval());
    taylor_model_class.def(self+=self);
    taylor_model_class.def(self-=self);
    taylor_model_class.def(self>double());
    taylor_model_class.def(self>self);
    taylor_model_class.def(self<self);
    taylor_model_class.def(self_ns::str(self));

    taylor_model_class.def("constant",(ValidatedTaylorModel(*)(Nat, const Interval&,Sweeper))&ValidatedTaylorModel::constant);
    taylor_model_class.def("variable",(ValidatedTaylorModel(*)(Nat, Nat,Sweeper))&ValidatedTaylorModel::variable);

    taylor_model_class.staticmethod("constant");
    taylor_model_class.staticmethod("variable");
    //taylor_model_class.staticmethod("variables");

    taylor_model_class.def("restrict", (ValidatedTaylorModel&(ValidatedTaylorModel::*)(const Vector<Interval>&))&ValidatedTaylorModel::restrict, return_value_policy<reference_existing_object>());
    taylor_model_class.def("restrict", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&,uint,const Interval&))&restrict);
    taylor_model_class.def("preaffine", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&,uint,const Interval&,const Interval&))&preaffine);

    taylor_model_class.def("evaluate", (Interval(*)(const ValidatedTaylorModel&, const Vector<Interval>&))&evaluate);
    taylor_model_class.def("set",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&,uint,Interval))&partial_evaluate);

    def("max",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&,const ValidatedTaylorModel&))&max);
    def("min",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&,const ValidatedTaylorModel&))&min);
    def("abs",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&abs);

    def("neg",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&neg);
    def("rec",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&rec);
    def("sqr",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&sqr);
    def("pow",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&, int))&pow);

    def("sqrt", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&sqrt);
    def("exp", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&exp);
    def("log", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&log);
    def("sin", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&sin);
    def("cos", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&cos);
    def("tan", (ValidatedTaylorModel(*)(const ValidatedTaylorModel&))&tan);

    taylor_model_class.def("range1", (Interval(*)(const ValidatedTaylorModel&)) &_range1);
    taylor_model_class.def("range2", (Interval(*)(const ValidatedTaylorModel&)) &_range2);
    taylor_model_class.def("range3", (Interval(*)(const ValidatedTaylorModel&)) &_range3);

    def("split",(ValidatedTaylorModel(*)(const ValidatedTaylorModel&,uint,tribool)) &split);

    from_python< Vector<ValidatedTaylorModel> >();
    to_python< Vector<ValidatedTaylorModel> >();

/*
    class_< TMV > taylor_model_vector_class("TaylorModelVector");
    taylor_model_vector_class.def("__getitem__", &__getitem__<TMV,int,ValidatedTaylorModel>);
    taylor_model_vector_class.def("__setitem__", &__setitem__<TMV,int,ValidatedTaylorModel>);
    taylor_model_vector_class.def(self_ns::str(self));
*/
}

void export_scalar_function_model()
{
    class_<ValidatedScalarFunctionModel> scalar_function_model_class("ValidatedScalarFunctionModel",init<ValidatedScalarFunctionModel>());
    scalar_function_model_class.def(init<ScalarTaylorFunction>());
    scalar_function_model_class.def("argument_size", &ValidatedScalarFunctionModel::argument_size);
    scalar_function_model_class.def("domain", &ValidatedScalarFunctionModel::domain);
    scalar_function_model_class.def("codomain", &ValidatedScalarFunctionModel::codomain);
    scalar_function_model_class.def("range", &ValidatedScalarFunctionModel::range);
    scalar_function_model_class.def("clobber", &ValidatedScalarFunctionModel::clobber);
    scalar_function_model_class.def("error", &ValidatedScalarFunctionModel::error);
    scalar_function_model_class.def("__call__", (Interval(ValidatedScalarFunctionModel::*)(const IntervalVector&)const) &ValidatedScalarFunctionModel::operator());
    scalar_function_model_class.def(self+self);
    scalar_function_model_class.def(self-self);
    scalar_function_model_class.def(self*self);
    scalar_function_model_class.def(self/self);
    scalar_function_model_class.def(self+Interval());
    scalar_function_model_class.def(self-Interval());
    scalar_function_model_class.def(self*Interval());
    scalar_function_model_class.def(self/Interval());
    scalar_function_model_class.def(Interval()+self);
    scalar_function_model_class.def(Interval()-self);
    scalar_function_model_class.def(Interval()*self);
    scalar_function_model_class.def(Interval()/self);
    scalar_function_model_class.def("__str__", &__cstr__<ValidatedScalarFunctionModel>);
    scalar_function_model_class.def("__repr__", &__crepr__<ValidatedScalarFunctionModel>);
    //scalar_function_model_class.def("__repr__",&__repr__<ValidatedScalarFunctionModel>);

    def("evaluate", (Interval(*)(const ValidatedScalarFunctionModel&,const IntervalVector&)) &evaluate);
    //def("partial_evaluate", (ValidatedScalarFunctionModel(*)(const ValidatedScalarFunctionModel&,Nat,const Interval&)) &partial_evaluate);
    def("compose", (ValidatedScalarFunctionModel(*)(const ValidatedScalarFunctionModel&,const ValidatedVectorFunctionModel&)) &compose);
    def("compose", (ValidatedScalarFunctionModel(*)(const ValidatedScalarFunction&,const ValidatedVectorFunctionModel&)) &compose);

    def("unrestrict", (ValidatedScalarFunction(*)(const ValidatedScalarFunctionModel&)) &unrestrict);
}

void export_vector_function_model()
{
    class_<ValidatedVectorFunctionModel> vector_function_model_class("ValidatedVectorFunctionModel",init<ValidatedVectorFunctionModel>());
    vector_function_model_class.def(init<VectorTaylorFunction>());
    vector_function_model_class.def("result_size", &ValidatedVectorFunctionModel::result_size);
    vector_function_model_class.def("argument_size", &ValidatedVectorFunctionModel::argument_size);
    vector_function_model_class.def("domain", &ValidatedVectorFunctionModel::domain);
    vector_function_model_class.def("codomain", &ValidatedVectorFunctionModel::codomain);
    vector_function_model_class.def("range", &ValidatedVectorFunctionModel::range);
    //vector_function_model_class.def("__getslice__", (VectorTaylorFunction(*)(const VectorTaylorFunction&,int,int))&__getslice__);
    vector_function_model_class.def("__getitem__", &__getitem__<ValidatedVectorFunctionModel,uint,ValidatedScalarFunctionModel>);
    vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorFunctionModel,uint,ValidatedScalarFunctionModel>);
    //vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorFunctionModel,uint,ValidatedScalarFunction>);
    vector_function_model_class.def("__call__", (IntervalVector(ValidatedVectorFunctionModel::*)(const IntervalVector&)const) &ValidatedVectorFunctionModel::operator());
    vector_function_model_class.def(self*Interval());
    vector_function_model_class.def("__str__", &__cstr__<ValidatedVectorFunctionModel>);
    vector_function_model_class.def("__repr__", &__crepr__<ValidatedVectorFunctionModel>);
    //export_vector_function_model.def("__repr__",&__repr__<ValidatedVectorFunctionModel>);

    def("evaluate", (IntervalVector(*)(const ValidatedVectorFunctionModel&,const IntervalVector&)) &evaluate);
    def("partial_evaluate", (ValidatedVectorFunctionModel(*)(const ValidatedVectorFunctionModel&,Nat,const Interval&)) &partial_evaluate);
    def("compose", (ValidatedVectorFunctionModel(*)(const ValidatedVectorFunctionModel&,const ValidatedVectorFunctionModel&)) &compose);
    def("compose", (ValidatedVectorFunctionModel(*)(const ValidatedVectorFunction&,const ValidatedVectorFunctionModel&)) &compose);

    def("unrestrict", (ValidatedVectorFunction(*)(const ValidatedVectorFunctionModel&)) &unrestrict);

    def("join", (ValidatedVectorFunctionModel(*)(const ValidatedVectorFunctionModel&,const ValidatedScalarFunctionModel&)) &join);

    to_python< List<ValidatedVectorFunctionModel> >();

}



void export_scalar_taylor_function()
{
    typedef Vector<Float> FloatVector;
    typedef Vector<Interval> IntervalVector;

    class_<ScalarTaylorFunction> scalar_taylor_function_class("ScalarTaylorFunction",init<ScalarTaylorFunction>());
    scalar_taylor_function_class.def(init<Box,ValidatedTaylorModel>());
    scalar_taylor_function_class.def(init< Box,Sweeper >());
    scalar_taylor_function_class.def(init< Box, const EffectiveScalarFunction&,Sweeper >());
    scalar_taylor_function_class.def(init< Box, Expansion<Float>, Float, Sweeper >());
    scalar_taylor_function_class.def("error", (const Float&(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::error, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("set_error", (void(ScalarTaylorFunction::*)(const Float&)) &ScalarTaylorFunction::set_error);
    scalar_taylor_function_class.def("argument_size", &ScalarTaylorFunction::argument_size);
    scalar_taylor_function_class.def("domain", &ScalarTaylorFunction::domain, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("codomain", &ScalarTaylorFunction::codomain);
    scalar_taylor_function_class.def("centre", &ScalarTaylorFunction::centre);
    scalar_taylor_function_class.def("range", &ScalarTaylorFunction::range);
    scalar_taylor_function_class.def("model", (const ValidatedTaylorModel&(ScalarTaylorFunction::*)()const)&ScalarTaylorFunction::model, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("polynomial", (Polynomial<Interval>(ScalarTaylorFunction::*)()const)&ScalarTaylorFunction::polynomial);
    scalar_taylor_function_class.def("number_of_nonzeros", (uint(ScalarTaylorFunction::*)()const)&ScalarTaylorFunction::number_of_nonzeros);
    scalar_taylor_function_class.def("set_sweeper", &ScalarTaylorFunction::set_sweeper);
    scalar_taylor_function_class.def("sweeper", &ScalarTaylorFunction::sweeper);
    scalar_taylor_function_class.def("sweep", (ScalarTaylorFunction&(ScalarTaylorFunction::*)()) &ScalarTaylorFunction::sweep, return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("__getitem__", &__getitem__<ScalarTaylorFunction,MultiIndex,Float>);
    scalar_taylor_function_class.def("__setitem__",&__setitem__<ScalarTaylorFunction,MultiIndex,Float>);
    scalar_taylor_function_class.def(+self);
    scalar_taylor_function_class.def(-self);
    scalar_taylor_function_class.def(self+self);
    scalar_taylor_function_class.def(self-self);
    scalar_taylor_function_class.def(self*self);
    scalar_taylor_function_class.def(self/self);
    scalar_taylor_function_class.def(self+Float());
    scalar_taylor_function_class.def(self-Float());
    scalar_taylor_function_class.def(self*Float());
    scalar_taylor_function_class.def(self/Float());
    scalar_taylor_function_class.def(Float()+self);
    scalar_taylor_function_class.def(Float()-self);
    scalar_taylor_function_class.def(Float()*self);
    scalar_taylor_function_class.def(Float()/self);
    scalar_taylor_function_class.def(self+=Float());
    scalar_taylor_function_class.def(self-=Float());
    scalar_taylor_function_class.def(self*=Float());
    scalar_taylor_function_class.def(self/=Float());
    scalar_taylor_function_class.def(Interval()+self);
    scalar_taylor_function_class.def(Interval()-self);
    scalar_taylor_function_class.def(Interval()*self);
    scalar_taylor_function_class.def(Interval()/self);
    scalar_taylor_function_class.def(self+Interval());
    scalar_taylor_function_class.def(self-Interval());
    scalar_taylor_function_class.def(self*Interval());
    scalar_taylor_function_class.def(self/Interval());
    scalar_taylor_function_class.def(self+=Interval());
    scalar_taylor_function_class.def(self-=Interval());
    scalar_taylor_function_class.def(self*=Interval());
    scalar_taylor_function_class.def(self/=Interval());
    scalar_taylor_function_class.def(self+EffectiveScalarFunction());
    scalar_taylor_function_class.def(self-EffectiveScalarFunction());
    scalar_taylor_function_class.def(self*EffectiveScalarFunction());
    scalar_taylor_function_class.def(self/EffectiveScalarFunction());
    scalar_taylor_function_class.def(EffectiveScalarFunction()+self);
    scalar_taylor_function_class.def(EffectiveScalarFunction()-self);
    scalar_taylor_function_class.def(EffectiveScalarFunction()*self);
    scalar_taylor_function_class.def(EffectiveScalarFunction()/self);
    scalar_taylor_function_class.def(self+=self);
    scalar_taylor_function_class.def(self-=self);
    scalar_taylor_function_class.def(self>Float());
    scalar_taylor_function_class.def(self>self);
    scalar_taylor_function_class.def(self<self);
    scalar_taylor_function_class.def("__str__", &__cstr__<ScalarTaylorFunction>);
    scalar_taylor_function_class.def("__repr__", &__crepr__<ScalarTaylorFunction>);
    scalar_taylor_function_class.def("__mul__",&__mul__< VectorTaylorFunction, ScalarTaylorFunction, Vector<Float> >);

    //scalar_taylor_function_class.def("__str__",(std::string(*)(const ScalarTaylorFunction&)) &__str__);
    //scalar_taylor_function_class.def("__cstr__",(std::string(*)(const ScalarTaylorFunction&)) &__cstr__);
    //scalar_taylor_function_class.def("__repr__",(std::string(*)(const ScalarTaylorFunction&)) &__repr__);
    scalar_taylor_function_class.def("value", (const Float&(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::value,return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("sweep", (ScalarTaylorFunction&(ScalarTaylorFunction::*)())&ScalarTaylorFunction::sweep,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("clobber", (ScalarTaylorFunction&(ScalarTaylorFunction::*)()) &ScalarTaylorFunction::clobber,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("set_sweeper",&ScalarTaylorFunction::set_sweeper);
    scalar_taylor_function_class.def("sweeper",&ScalarTaylorFunction::sweeper);
    //scalar_taylor_function_class.def("set_default_maximum_degree",&ScalarTaylorFunction::set_default_maximum_degree);
    //scalar_taylor_function_class.def("set_default_sweep_threshold",&ScalarTaylorFunction::set_default_sweep_threshold);
    //scalar_taylor_function_class.def("default_maximum_degree",&ScalarTaylorFunction::default_maximum_degree);
    //scalar_taylor_function_class.def("default_sweep_threshold",&ScalarTaylorFunction::default_sweep_threshold);
    scalar_taylor_function_class.def("__call__", (Float(ScalarTaylorFunction::*)(const Vector<Float>&)const) &ScalarTaylorFunction::evaluate);
    scalar_taylor_function_class.def("__call__", (Interval(ScalarTaylorFunction::*)(const Vector<Interval>&)const) &ScalarTaylorFunction::evaluate);
    scalar_taylor_function_class.def("evaluate", (Float(ScalarTaylorFunction::*)(const Vector<Float>&)const) &ScalarTaylorFunction::evaluate);
    scalar_taylor_function_class.def("evaluate", (Interval(ScalarTaylorFunction::*)(const Vector<Interval>&)const) &ScalarTaylorFunction::evaluate);
    //scalar_taylor_function_class.def("gradient", (Vector<Interval>(ScalarTaylorFunction::*)(const Vector<Interval>&)const) &ScalarTaylorFunction::gradient);
    scalar_taylor_function_class.def("function", (EffectiveScalarFunction(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::function);
    scalar_taylor_function_class.def("polynomial", (Polynomial<Interval>(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::polynomial);
    scalar_taylor_function_class.def("set", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,uint j, const Interval&)) &partial_evaluate);
    scalar_taylor_function_class.def("restrict", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const Box&)) &restrict);
    scalar_taylor_function_class.def("restrict", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,uint,const Interval&)) &restrict);
    scalar_taylor_function_class.def("extend", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const Box&)) &extend);
    //scalar_taylor_function_class.staticmethod("set_default_maximum_degree");
    //scalar_taylor_function_class.staticmethod("set_default_sweep_threshold");
    //scalar_taylor_function_class.staticmethod("default_maximum_degree");
    //scalar_taylor_function_class.staticmethod("default_sweep_threshold");



    scalar_taylor_function_class.def("zero",(ScalarTaylorFunction(*)(const Box&,Sweeper))&ScalarTaylorFunction::zero);
    scalar_taylor_function_class.def("constant",(ScalarTaylorFunction(*)(const Box&,const Float&,Sweeper))&ScalarTaylorFunction::constant);
    scalar_taylor_function_class.def("constant",(ScalarTaylorFunction(*)(const Box&,const Interval&,Sweeper))&ScalarTaylorFunction::constant);
    scalar_taylor_function_class.def("coordinate",(ScalarTaylorFunction(*)(const Box&,uint,Sweeper))&ScalarTaylorFunction::coordinate);
    scalar_taylor_function_class.def("variable",(ScalarTaylorFunction(*)(const Box&,uint,Sweeper))&ScalarTaylorFunction::variable);
    scalar_taylor_function_class.def("variables",(Vector<ScalarTaylorFunction>(*)(const Box&,Sweeper)) &ScalarTaylorFunction::variables);


    scalar_taylor_function_class.staticmethod("constant");
    scalar_taylor_function_class.staticmethod("coordinate");
    scalar_taylor_function_class.staticmethod("variable");
    scalar_taylor_function_class.staticmethod("variables");

    def("split", (std::pair<ScalarTaylorFunction,ScalarTaylorFunction>(*)(const ScalarTaylorFunction&,uint)) &Ariadne::split);
    def("evaluate",(Interval(*)(const ScalarTaylorFunction&,const IntervalVector&)) &evaluate);
    def("partial_evaluate",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,uint,const Interval&)) &partial_evaluate);

    def("midpoint",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&)) &midpoint);
    def("evaluate",(Interval(*)(const ScalarTaylorFunction&,const IntervalVector&)) &evaluate);
    def("derivative",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,uint)) &derivative);
    def("antiderivative",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,uint)) &antiderivative);
    def("antiderivative",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,uint,Float)) &antiderivative);

    def("refines",(bool(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &refines);
    def("disjoint",(bool(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &disjoint);
    def("intersection",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &intersection);

    def("restrict",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const Box&)) &restrict);
    def("restrict",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,uint,const Interval&)) &restrict);

    def("embed",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const Interval&)) &embed);
    def("embed",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const Box&)) &embed);
    def("embed",(ScalarTaylorFunction(*)(const Box&,const ScalarTaylorFunction&)) &embed);

    def("max",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&))&max);
    def("min",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&))&min);
    def("abs",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&abs);

    def("neg",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&neg);
    def("rec",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&rec);
    def("sqr",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&sqr);
    def("pow",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&, int))&pow);

    def("sqrt", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&sqrt);
    def("exp", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&exp);
    def("log", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&log);
    def("sin", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&sin);
    def("cos", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&cos);
    def("tan", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&))&tan);

    to_python< Vector<ScalarTaylorFunction> >();
}


void export_vector_taylor_function()
{
    typedef uint Nat;
    typedef Vector<Float> FloatVector;
    typedef Vector<Interval> IntervalVector;
    typedef Matrix<Float> FloatMatrix;
    typedef ScalarTaylorFunction ScalarTaylorFunction;
    typedef VectorTaylorFunction VectorTaylorFunction;

    class_<VectorTaylorFunction> vector_taylor_function_class("VectorTaylorFunction", init<VectorTaylorFunction>());
    vector_taylor_function_class.def( init< Nat, Box, Sweeper >());
    vector_taylor_function_class.def( init< Box,const EffectiveVectorFunction&,Sweeper >());
    vector_taylor_function_class.def(init< Box, Vector< Expansion<Float> >, Vector<Float>, Sweeper >());
    vector_taylor_function_class.def( init< Vector<ScalarTaylorFunction> >());
    vector_taylor_function_class.def("__len__", &VectorTaylorFunction::result_size);
    vector_taylor_function_class.def("result_size", &VectorTaylorFunction::result_size);
    vector_taylor_function_class.def("argument_size", &VectorTaylorFunction::argument_size);
    vector_taylor_function_class.def("domain", &VectorTaylorFunction::domain, return_value_policy<copy_const_reference>());
    vector_taylor_function_class.def("codomain", &VectorTaylorFunction::codomain);
    vector_taylor_function_class.def("models", &VectorTaylorFunction::models, return_value_policy<copy_const_reference>());
    vector_taylor_function_class.def("centre", &VectorTaylorFunction::centre);
    vector_taylor_function_class.def("range", &VectorTaylorFunction::range);
    vector_taylor_function_class.def("errors", &VectorTaylorFunction::errors);
    vector_taylor_function_class.def("sweep", (VectorTaylorFunction&(VectorTaylorFunction::*)())&VectorTaylorFunction::sweep,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("clobber", (VectorTaylorFunction&(VectorTaylorFunction::*)()) &VectorTaylorFunction::clobber,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("set_sweeper",&VectorTaylorFunction::set_sweeper);
    vector_taylor_function_class.def("sweeper",&VectorTaylorFunction::sweeper);
    vector_taylor_function_class.def("sweep", (VectorTaylorFunction&(VectorTaylorFunction::*)()) &VectorTaylorFunction::sweep, return_value_policy<reference_existing_object>());
    //vector_taylor_function_class.def("__getslice__", &__getslice__<VectorTaylorFunction,Nat,Nat,ScalarTaylorFunction>);
    vector_taylor_function_class.def("__getslice__", (VectorTaylorFunction(*)(const VectorTaylorFunction&,int,int))&__getslice__);
    vector_taylor_function_class.def("__getitem__", &__getitem__<VectorTaylorFunction,uint,ScalarTaylorFunction>);
    vector_taylor_function_class.def("__setitem__",&__setitem__<VectorTaylorFunction,uint,ScalarTaylorFunction>);
    vector_taylor_function_class.def(-self);
    vector_taylor_function_class.def(self+self);
    vector_taylor_function_class.def(self-self);
    vector_taylor_function_class.def(self+IntervalVector());
    vector_taylor_function_class.def(self-IntervalVector());
    vector_taylor_function_class.def(self*Interval());
    vector_taylor_function_class.def(self/Interval());
    vector_taylor_function_class.def(self*ScalarTaylorFunction());
    vector_taylor_function_class.def(self+=IntervalVector());
    vector_taylor_function_class.def(self-=IntervalVector());
    vector_taylor_function_class.def(self*=Interval());
    vector_taylor_function_class.def(self/=Interval());
    vector_taylor_function_class.def(self+=self);
    vector_taylor_function_class.def(self-=self);
    vector_taylor_function_class.def("__str__", &__cstr__<VectorTaylorFunction>);
    vector_taylor_function_class.def("__repr__", &__crepr__<VectorTaylorFunction>);
    vector_taylor_function_class.def("clobber", (VectorTaylorFunction&(VectorTaylorFunction::*)()) &VectorTaylorFunction::clobber,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("__call__", (FloatVector(VectorTaylorFunction::*)(const Vector<Float>&)const) &VectorTaylorFunction::evaluate);
    vector_taylor_function_class.def("__call__", (IntervalVector(VectorTaylorFunction::*)(const IntervalVector&)const) &VectorTaylorFunction::evaluate);
    vector_taylor_function_class.def("evaluate", (FloatVector(VectorTaylorFunction::*)(const Vector<Float>&)const) &VectorTaylorFunction::evaluate);
    vector_taylor_function_class.def("evaluate", (IntervalVector(VectorTaylorFunction::*)(const IntervalVector&)const) &VectorTaylorFunction::evaluate);
    //vector_taylor_function_class.def("jacobian", (Vector<Interval>(ScalarTaylorFunction::*)(const Vector<Interval>&)const) &ScalarTaylorFunction::gradient);
    vector_taylor_function_class.def("polynomials", (Vector< Polynomial<Interval> >(VectorTaylorFunction::*)()const) &VectorTaylorFunction::polynomials);
    vector_taylor_function_class.def("function", (EffectiveVectorFunction(VectorTaylorFunction::*)()const) &VectorTaylorFunction::function);


    vector_taylor_function_class.def("constant",(VectorTaylorFunction(*)(const Box&, const FloatVector&,Sweeper))&VectorTaylorFunction::constant);
    vector_taylor_function_class.def("constant",(VectorTaylorFunction(*)(const Box&, const IntervalVector&,Sweeper))&VectorTaylorFunction::constant);
    vector_taylor_function_class.def("identity",(VectorTaylorFunction(*)(const Box&,Sweeper))&VectorTaylorFunction::identity);

    vector_taylor_function_class.staticmethod("constant");
    vector_taylor_function_class.staticmethod("identity");

    def("refines", (bool(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &refines);
    def("disjoint", (bool(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &disjoint);
    def("intersection", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &intersection);

    def("join", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &join);
    def("join", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const ScalarTaylorFunction&)) &join);
    def("join", (VectorTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &join);
    def("join", (VectorTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &join);

    def("combine", (VectorTaylorFunction(*)(const ScalarTaylorFunction&,const ScalarTaylorFunction&)) &combine);
    def("combine", (VectorTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &combine);
    def("combine", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const ScalarTaylorFunction&)) &combine);
    def("combine", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &combine);

    def("embed",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const Interval&)) &embed);
    def("embed",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const Box&)) &embed);
    def("embed",(VectorTaylorFunction(*)(const Box&,const VectorTaylorFunction&)) &embed);

    def("restrict", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const Box&)) &restrict);
    def("restrict", (VectorTaylorFunction(*)(const VectorTaylorFunction&,uint,const Interval&)) &restrict);

    def("split", (std::pair<VectorTaylorFunction,VectorTaylorFunction>(*)(const VectorTaylorFunction&,uint)) &Ariadne::split);

    def("evaluate",(FloatVector(VectorTaylorFunction::*)(const FloatVector&)const) &VectorTaylorFunction::evaluate);
    def("evaluate",(IntervalVector(VectorTaylorFunction::*)(const IntervalVector&)const) &VectorTaylorFunction::evaluate);
    def("partial_evaluate",(VectorTaylorFunction(*)(const VectorTaylorFunction&,uint,const Interval&)) &partial_evaluate);
    //def("compose",(ScalarTaylorFunction(*)(const RP&,const VectorTaylorFunction&)) &compose);
    def("compose",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(ScalarTaylorFunction(*)(const EffectiveScalarFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(VectorTaylorFunction(*)(const EffectiveVectorFunction&,const VectorTaylorFunction&)) &compose);
    def("derivative",(VectorTaylorFunction(*)(const VectorTaylorFunction&,Nat)) &derivative);
    def("antiderivative",(VectorTaylorFunction(*)(const VectorTaylorFunction&,Nat)) &antiderivative);
    def("antiderivative",(VectorTaylorFunction(*)(const VectorTaylorFunction&,Nat,Float)) &antiderivative);

    def("unchecked_compose",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &unchecked_compose);
    def("unchecked_compose",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &unchecked_compose);

    from_python<VectorTaylorFunction>();
}

void calculus_submodule()
{
    export_expansion();
    export_sweeper();
    export_taylor_model();
    export_scalar_function_model();
    export_vector_function_model();
    export_scalar_taylor_function();
    export_vector_taylor_function();
}


