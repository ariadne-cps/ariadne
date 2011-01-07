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

#include <boost/python.hpp>

#include "function_interface.h"
#include "polynomial.h"
#include "function.h"
#include "taylor_function.h"

#include "utilities.h"

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


VectorTaylorFunction __getslice__(const VectorTaylorFunction& tf, int start, int stop) {
    if(start<0) { start+=tf.result_size(); }
    if(stop<0) { stop+=tf.result_size(); }
    ARIADNE_ASSERT_MSG(0<=start&&start<=stop&&uint(stop)<=tf.result_size(),
            "result_size="<<tf.result_size()<<", start="<<start<<", stop="<<stop);
    return VectorTaylorFunction(tf.domain(),Vector<IntervalTaylorModel>(project(tf.models(),range(start,stop))));
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
        if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((boost::python::converter::rvalue_from_python_storage<MultiIndex>*)data)->storage.bytes;
        Expansion<T> r;
        boost::python::dict dct=boost::python::extract<boost::python::dict>(obj_ptr);
        boost::python::list lst=dct.values();
        for(int i=0; i!=len(lst); ++i) {
            boost::python::tuple tup=boost::python::extract<boost::python::tuple>(lst[i]);
            MultiIndex a=boost::python::extract<MultiIndex>(tup[0]);
            T c=extract<T>(tup[1]);
            r.append(a,c);
        }
        new (storage) Expansion<T>(r);
        r.unique_sort();
        data->convertible = storage;
    }
};


template<>
struct from_python< Vector<IntervalTaylorModel> > {
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< Vector<IntervalTaylorModel> >()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((converter::rvalue_from_python_storage<Interval>*)data)->storage.bytes;
        boost::python::list lst=extract<boost::python::list>(obj_ptr);
        Vector<IntervalTaylorModel>* tms_ptr = new (storage) Vector<IntervalTaylorModel>(len(lst));
        for(int i=0; i!=len(lst); ++i) { (*tms_ptr)[i]=extract<IntervalTaylorModel>(lst[i]); }
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
        void* storage = ((boost::python::converter::rvalue_from_python_storage<MultiIndex>*)data)->storage.bytes;
        list lst=extract<boost::python::list>(obj_ptr);
        VectorTaylorFunction* tf_ptr = new (storage) VectorTaylorFunction(len(lst));
        for(uint i=0; i!=tf_ptr->result_size(); ++i) { tf_ptr->set(i,extract<ScalarTaylorFunction>(lst[i])); }
        data->convertible = storage;
    }
};

List<MultiIndex> keys(const IntervalTaylorModel& tm) {
    List<MultiIndex> r;
    for(IntervalTaylorModel::const_iterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        r.append(iter->key());
    }
    return r;
}



Interval _range1(const IntervalTaylorModel&);
Interval _range2(const IntervalTaylorModel&);
Interval _range3(const IntervalTaylorModel&);


} // namespace Ariadne



void export_taylor_model()
{
    typedef int Z;
    typedef uint N;
    typedef double D;
    typedef Float Float;
    typedef Interval Interval;
    typedef MultiIndex A;
    typedef Vector<Float> RV;
    typedef Vector<Interval> IV;
    typedef Matrix<Float> RMx;
    typedef IntervalTaylorModel IntervalTaylorModel;
    typedef Vector<IntervalTaylorModel> TMV;
    typedef VectorTaylorFunction VectorTaylorFunction;
    typedef Polynomial<Float> RP;
    typedef Polynomial<Interval> IP;
    typedef RealScalarFunction FI;

    def("set_default_sweep_threshold",&TaylorModelAccuracy::set_default_sweep_threshold);

    class_< ExpansionValue<Float> > expansion_value_class("ExpansionValue", init<MultiIndex,Float>());
    // TODO: Add get/set for data
    // TODO: Use property for key
    //expansion_value_class.add_property("key", (MultiIndex const&(ExpansionValue<Float>::*)()const)&ExpansionValue<Float>::key);
    expansion_value_class.def("key", (const MultiIndex&(ExpansionValue<Float>::*)()const)&ExpansionValue<Float>::key, return_value_policy<copy_const_reference>());
    expansion_value_class.def(self_ns::str(self));

    class_<IntervalTaylorModel> taylor_model_class("IntervalTaylorModel", init<IntervalTaylorModel>());
    taylor_model_class.def( init< N >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const IntervalTaylorModel&))&keys);
    taylor_model_class.def("error", (const Float&(IntervalTaylorModel::*)()const) &IntervalTaylorModel::error, return_value_policy<copy_const_reference>());
    taylor_model_class.def("set_error", (void(IntervalTaylorModel::*)(const Float&)) &IntervalTaylorModel::set_error);
    taylor_model_class.def("argument_size", &IntervalTaylorModel::argument_size);
    taylor_model_class.def("domain", &IntervalTaylorModel::domain);
    taylor_model_class.def("range", &IntervalTaylorModel::range);
    taylor_model_class.def("set_sweep_threshold", &IntervalTaylorModel::set_sweep_threshold);
    taylor_model_class.def("get_sweep_threshold", &IntervalTaylorModel::sweep_threshold);
    taylor_model_class.def("set_maximum_degree", &IntervalTaylorModel::set_maximum_degree);
    taylor_model_class.def("get_maximum_degree", &IntervalTaylorModel::maximum_degree);
    taylor_model_class.def("__getitem__", &__getitem__<IntervalTaylorModel,MultiIndex,Float>);
    taylor_model_class.def("__setitem__",&__setitem__<IntervalTaylorModel,MultiIndex,Float>);
    taylor_model_class.def("__setitem__",&__setitem__<IntervalTaylorModel,MultiIndex,double>);
    taylor_model_class.def("__iter__", iterator<IntervalTaylorModel>()); // TODO: IntervalTaylorModel iterator does not fit in general Python iterator scheme
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

    taylor_model_class.def("constant",(IntervalTaylorModel(*)(N, const Interval&))&IntervalTaylorModel::constant);
    taylor_model_class.def("variable",(IntervalTaylorModel(*)(N, N))&IntervalTaylorModel::variable);

    taylor_model_class.staticmethod("constant");
    taylor_model_class.staticmethod("variable");
    //taylor_model_class.staticmethod("variables");

    taylor_model_class.def("restrict", (IntervalTaylorModel&(IntervalTaylorModel::*)(const Vector<Interval>&))&IntervalTaylorModel::restrict, return_value_policy<reference_existing_object>());
    taylor_model_class.def("restrict", (IntervalTaylorModel(*)(const IntervalTaylorModel&,uint,const Interval&))&restrict);
    taylor_model_class.def("preaffine", (IntervalTaylorModel(*)(const IntervalTaylorModel&,uint,const Interval&,const Interval&))&preaffine);

    taylor_model_class.def("evaluate", (Interval(IntervalTaylorModel::*)(const Vector<Interval>&)const)&IntervalTaylorModel::evaluate);
    taylor_model_class.def("set",(IntervalTaylorModel(*)(const IntervalTaylorModel&,uint,Interval))&partial_evaluate);

    def("max",(IntervalTaylorModel(*)(const IntervalTaylorModel&,const IntervalTaylorModel&))&max);
    def("min",(IntervalTaylorModel(*)(const IntervalTaylorModel&,const IntervalTaylorModel&))&min);
    def("abs",(IntervalTaylorModel(*)(const IntervalTaylorModel&))&abs);

    def("neg",(IntervalTaylorModel(*)(const IntervalTaylorModel&))&neg);
    def("rec",(IntervalTaylorModel(*)(const IntervalTaylorModel&))&rec);
    def("sqr",(IntervalTaylorModel(*)(const IntervalTaylorModel&))&sqr);
    def("pow",(IntervalTaylorModel(*)(const IntervalTaylorModel&, int))&pow);

    def("sqrt", (IntervalTaylorModel(*)(const IntervalTaylorModel&))&sqrt);
    def("exp", (IntervalTaylorModel(*)(const IntervalTaylorModel&))&exp);
    def("log", (IntervalTaylorModel(*)(const IntervalTaylorModel&))&log);
    def("sin", (IntervalTaylorModel(*)(const IntervalTaylorModel&))&sin);
    def("cos", (IntervalTaylorModel(*)(const IntervalTaylorModel&))&cos);
    def("tan", (IntervalTaylorModel(*)(const IntervalTaylorModel&))&tan);

    taylor_model_class.def("range1", (Interval(*)(const IntervalTaylorModel&)) &_range1);
    taylor_model_class.def("range2", (Interval(*)(const IntervalTaylorModel&)) &_range2);
    taylor_model_class.def("range3", (Interval(*)(const IntervalTaylorModel&)) &_range3);

    def("split",(IntervalTaylorModel(*)(const IntervalTaylorModel&,uint,tribool)) &split);

    from_python< Vector<IntervalTaylorModel> >();
    to_python< Vector<IntervalTaylorModel> >();

/*
    class_< TMV > taylor_model_vector_class("TaylorModelVector");
    taylor_model_vector_class.def("__getitem__", &__getitem__<TMV,int,IntervalTaylorModel>);
    taylor_model_vector_class.def("__setitem__", &__setitem__<TMV,int,IntervalTaylorModel>);
    taylor_model_vector_class.def(self_ns::str(self));
*/
}

void export_scalar_taylor_function()
{
    typedef Vector<Float> FloatVector;
    typedef Vector<Interval> IntervalVector;

    class_<ScalarTaylorFunction> scalar_taylor_function_class("ScalarTaylorFunction",init<ScalarTaylorFunction>());
    scalar_taylor_function_class.def(init<IntervalVector,IntervalTaylorModel>());
    scalar_taylor_function_class.def(init< IntervalVector >());
    scalar_taylor_function_class.def(init< IntervalVector, const RealScalarFunction& >());
    scalar_taylor_function_class.def("error", (const Float&(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::error, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("set_error", (void(ScalarTaylorFunction::*)(const Float&)) &ScalarTaylorFunction::set_error);
    scalar_taylor_function_class.def("argument_size", &ScalarTaylorFunction::argument_size);
    scalar_taylor_function_class.def("domain", &ScalarTaylorFunction::domain, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("codomain", &ScalarTaylorFunction::codomain);
    scalar_taylor_function_class.def("centre", &ScalarTaylorFunction::centre);
    scalar_taylor_function_class.def("range", &ScalarTaylorFunction::range);
    scalar_taylor_function_class.def("model", (const IntervalTaylorModel&(ScalarTaylorFunction::*)()const)&ScalarTaylorFunction::model, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("number_of_nonzeros", (uint(ScalarTaylorFunction::*)()const)&ScalarTaylorFunction::number_of_nonzeros);
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
    scalar_taylor_function_class.def(self+=self);
    scalar_taylor_function_class.def(self-=self);
    scalar_taylor_function_class.def(self>Float());
    scalar_taylor_function_class.def(self>self);
    scalar_taylor_function_class.def(self<self);
    scalar_taylor_function_class.def(self_ns::str(self));
    scalar_taylor_function_class.def("__mul__",&__mul__< VectorTaylorFunction, ScalarTaylorFunction, Vector<Float> >);

    //scalar_taylor_function_class.def("__str__",(std::string(*)(const ScalarTaylorFunction&)) &__str__);
    //scalar_taylor_function_class.def("__cstr__",(std::string(*)(const ScalarTaylorFunction&)) &__cstr__);
    //scalar_taylor_function_class.def("__repr__",(std::string(*)(const ScalarTaylorFunction&)) &__repr__);
    scalar_taylor_function_class.def("value", (const Float&(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::value,return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("truncate", (ScalarTaylorFunction&(ScalarTaylorFunction::*)(uint)) &ScalarTaylorFunction::truncate,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("sweep", (ScalarTaylorFunction&(ScalarTaylorFunction::*)(double))&ScalarTaylorFunction::sweep,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("truncate", (ScalarTaylorFunction&(ScalarTaylorFunction::*)()) &ScalarTaylorFunction::truncate,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("sweep", (ScalarTaylorFunction&(ScalarTaylorFunction::*)())&ScalarTaylorFunction::sweep,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("clobber", (ScalarTaylorFunction&(ScalarTaylorFunction::*)()) &ScalarTaylorFunction::clobber,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("set_maximum_degree",&ScalarTaylorFunction::set_maximum_degree);
    scalar_taylor_function_class.def("set_sweep_threshold",&ScalarTaylorFunction::set_sweep_threshold);
    scalar_taylor_function_class.def("maximum_degree",&ScalarTaylorFunction::maximum_degree);
    scalar_taylor_function_class.def("sweep_threshold",&ScalarTaylorFunction::sweep_threshold);
    //scalar_taylor_function_class.def("set_default_maximum_degree",&ScalarTaylorFunction::set_default_maximum_degree);
    //scalar_taylor_function_class.def("set_default_sweep_threshold",&ScalarTaylorFunction::set_default_sweep_threshold);
    //scalar_taylor_function_class.def("default_maximum_degree",&ScalarTaylorFunction::default_maximum_degree);
    //scalar_taylor_function_class.def("default_sweep_threshold",&ScalarTaylorFunction::default_sweep_threshold);
    scalar_taylor_function_class.def("__call__", (Float(ScalarTaylorFunction::*)(const Vector<Float>&)const) &ScalarTaylorFunction::evaluate);
    scalar_taylor_function_class.def("__call__", (Interval(ScalarTaylorFunction::*)(const Vector<Interval>&)const) &ScalarTaylorFunction::evaluate);
    scalar_taylor_function_class.def("evaluate", (Float(ScalarTaylorFunction::*)(const Vector<Float>&)const) &ScalarTaylorFunction::evaluate);
    scalar_taylor_function_class.def("evaluate", (Interval(ScalarTaylorFunction::*)(const Vector<Interval>&)const) &ScalarTaylorFunction::evaluate);
    scalar_taylor_function_class.def("function", (RealScalarFunction(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::function);
    scalar_taylor_function_class.def("polynomial", (Polynomial<Interval>(ScalarTaylorFunction::*)()const) &ScalarTaylorFunction::polynomial);
    scalar_taylor_function_class.def("set", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,uint j, const Interval&)) &partial_evaluate);
    scalar_taylor_function_class.def("restrict", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const Vector<Interval>&)) &restrict);
    scalar_taylor_function_class.def("restrict", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,uint,const Interval&)) &restrict);
    scalar_taylor_function_class.def("extend", (ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const Vector<Interval>&)) &extend);
    //scalar_taylor_function_class.staticmethod("set_default_maximum_degree");
    //scalar_taylor_function_class.staticmethod("set_default_sweep_threshold");
    //scalar_taylor_function_class.staticmethod("default_maximum_degree");
    //scalar_taylor_function_class.staticmethod("default_sweep_threshold");



    scalar_taylor_function_class.def("constant",(ScalarTaylorFunction(*)(const IntervalVector&, const Float&))&ScalarTaylorFunction::constant);
    scalar_taylor_function_class.def("constant",(ScalarTaylorFunction(*)(const IntervalVector&, const Interval&))&ScalarTaylorFunction::constant);
    scalar_taylor_function_class.def("coordinate",(ScalarTaylorFunction(*)(const IntervalVector&, uint))&ScalarTaylorFunction::coordinate);
    scalar_taylor_function_class.def("variable",(ScalarTaylorFunction(*)(const IntervalVector&, uint))&ScalarTaylorFunction::variable);
    scalar_taylor_function_class.def("variables",(Vector<ScalarTaylorFunction>(*)(const IntervalVector&)) &ScalarTaylorFunction::variables);


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

    def("restrict",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const IntervalVector&)) &restrict);
    def("restrict",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,uint,const Interval&)) &restrict);

    def("embed",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const Interval&)) &embed);
    def("embed",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const IntervalVector&)) &embed);
    def("embed",(ScalarTaylorFunction(*)(const IntervalVector&,const ScalarTaylorFunction&)) &embed);

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
    typedef Polynomial<Float> RP;
    typedef Polynomial<Interval> IP;
    typedef Vector<Polynomial<Float> > RPV;
    typedef Vector<Polynomial<Interval> > IPV;
    typedef Vector<RealScalarFunction> EV;
    typedef ScalarTaylorFunction ScalarTaylorFunction;
    typedef VectorTaylorFunction VectorTaylorFunction;

    class_<VectorTaylorFunction> vector_taylor_function_class("VectorTaylorFunction", init<VectorTaylorFunction>());
    vector_taylor_function_class.def( init< IntervalVector,const RealVectorFunction& >());
    vector_taylor_function_class.def( init< Vector<ScalarTaylorFunction> >());
    vector_taylor_function_class.def("__len__", &VectorTaylorFunction::result_size);
    vector_taylor_function_class.def("result_size", &VectorTaylorFunction::result_size);
    vector_taylor_function_class.def("argument_size", &VectorTaylorFunction::argument_size);
    vector_taylor_function_class.def("domain", &VectorTaylorFunction::domain, return_value_policy<copy_const_reference>());
    vector_taylor_function_class.def("codomain", &VectorTaylorFunction::codomain);
    vector_taylor_function_class.def("centre", &VectorTaylorFunction::centre);
    vector_taylor_function_class.def("range", &VectorTaylorFunction::range);
    vector_taylor_function_class.def("truncate", (VectorTaylorFunction&(VectorTaylorFunction::*)(uint)) &VectorTaylorFunction::truncate,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("sweep", (VectorTaylorFunction&(VectorTaylorFunction::*)(double))&VectorTaylorFunction::sweep,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("truncate", (VectorTaylorFunction&(VectorTaylorFunction::*)()) &VectorTaylorFunction::truncate,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("sweep", (VectorTaylorFunction&(VectorTaylorFunction::*)())&VectorTaylorFunction::sweep,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("clobber", (VectorTaylorFunction&(VectorTaylorFunction::*)()) &VectorTaylorFunction::clobber,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("set_maximum_degree",&VectorTaylorFunction::set_maximum_degree);
    vector_taylor_function_class.def("set_sweep_threshold",&VectorTaylorFunction::set_sweep_threshold);
//    vector_taylor_function_class.def("maximum_degree",&VectorTaylorFunction::maximum_degree);
//    vector_taylor_function_class.def("sweep_threshold",&VectorTaylorFunction::sweep_threshold);
    //vector_taylor_function_class.def("__getslice__", &__getslice__<VectorTaylorFunction,Nat,Nat,ScalarTaylorFunction>);
    vector_taylor_function_class.def("__getslice__", (VectorTaylorFunction(*)(const VectorTaylorFunction&,int,int))&__getslice__);
    vector_taylor_function_class.def("__getitem__", &__getitem__<VectorTaylorFunction,uint,ScalarTaylorFunction>);
    vector_taylor_function_class.def("__setitem__",&__setitem__<VectorTaylorFunction,uint,ScalarTaylorFunction>);
    vector_taylor_function_class.def(-self);
    vector_taylor_function_class.def(self+self);
    vector_taylor_function_class.def(self-self);
    vector_taylor_function_class.def(self+FloatVector());
    vector_taylor_function_class.def(self-FloatVector());
    vector_taylor_function_class.def(self*Float());
    vector_taylor_function_class.def(self/Float());
    vector_taylor_function_class.def(self+=FloatVector());
    vector_taylor_function_class.def(self-=FloatVector());
    vector_taylor_function_class.def(self*=Float());
    vector_taylor_function_class.def(self/=Float());
    vector_taylor_function_class.def(self+=self);
    vector_taylor_function_class.def(self-=self);
    vector_taylor_function_class.def(self_ns::str(self));
    vector_taylor_function_class.def("clobber", (VectorTaylorFunction&(VectorTaylorFunction::*)()) &VectorTaylorFunction::clobber,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("__call__", (FloatVector(VectorTaylorFunction::*)(const Vector<Float>&)const) &VectorTaylorFunction::evaluate);
    vector_taylor_function_class.def("__call__", (IntervalVector(VectorTaylorFunction::*)(const IntervalVector&)const) &VectorTaylorFunction::evaluate);
    vector_taylor_function_class.def("evaluate", (FloatVector(VectorTaylorFunction::*)(const Vector<Float>&)const) &VectorTaylorFunction::evaluate);
    vector_taylor_function_class.def("evaluate", (IntervalVector(VectorTaylorFunction::*)(const IntervalVector&)const) &VectorTaylorFunction::evaluate);
    vector_taylor_function_class.def("polynomial", (Vector< Polynomial<Interval> >(VectorTaylorFunction::*)()const) &VectorTaylorFunction::polynomial);
    vector_taylor_function_class.def("function", (RealVectorFunction(VectorTaylorFunction::*)()const) &VectorTaylorFunction::function);


    vector_taylor_function_class.def("constant",(VectorTaylorFunction(*)(const IntervalVector&, const FloatVector&))&VectorTaylorFunction::constant);
    vector_taylor_function_class.def("constant",(VectorTaylorFunction(*)(const IntervalVector&, const IntervalVector&))&VectorTaylorFunction::constant);
    vector_taylor_function_class.def("identity",(VectorTaylorFunction(*)(const IntervalVector&))&VectorTaylorFunction::identity);

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
    def("embed",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const IntervalVector&)) &embed);
    def("embed",(VectorTaylorFunction(*)(const IntervalVector&,const VectorTaylorFunction&)) &embed);

    def("restrict", (VectorTaylorFunction(*)(const VectorTaylorFunction&,const Vector<Interval>&)) &restrict);
    def("restrict", (VectorTaylorFunction(*)(const VectorTaylorFunction&,uint,const Interval&)) &restrict);

    def("split", (std::pair<VectorTaylorFunction,VectorTaylorFunction>(*)(const VectorTaylorFunction&,uint)) &Ariadne::split);

    def("evaluate",(FloatVector(VectorTaylorFunction::*)(const FloatVector&)const) &VectorTaylorFunction::evaluate);
    def("evaluate",(IntervalVector(VectorTaylorFunction::*)(const IntervalVector&)const) &VectorTaylorFunction::evaluate);
    def("partial_evaluate",(VectorTaylorFunction(*)(const VectorTaylorFunction&,uint,const Interval&)) &partial_evaluate);
    //def("compose",(ScalarTaylorFunction(*)(const RP&,const VectorTaylorFunction&)) &compose);
    def("compose",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(ScalarTaylorFunction(*)(const RealScalarFunction&,const VectorTaylorFunction&)) &compose);
    def("compose",(VectorTaylorFunction(*)(const RealVectorFunction&,const VectorTaylorFunction&)) &compose);
    def("antiderivative",(VectorTaylorFunction(*)(const VectorTaylorFunction&,Nat)) &antiderivative);
    def("antiderivative",(VectorTaylorFunction(*)(const VectorTaylorFunction&,Nat,Float)) &antiderivative);

    def("unchecked_compose",(ScalarTaylorFunction(*)(const ScalarTaylorFunction&,const VectorTaylorFunction&)) &unchecked_compose);
    def("unchecked_compose",(VectorTaylorFunction(*)(const VectorTaylorFunction&,const VectorTaylorFunction&)) &unchecked_compose);

    from_python<VectorTaylorFunction>();
}

void calculus_submodule()
{
    export_taylor_model();
    export_scalar_taylor_function();
    export_vector_taylor_function();
}


