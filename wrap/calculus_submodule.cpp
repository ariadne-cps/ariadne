/***************************************************************************
 *            calculus_submodule.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "boost_python.hpp"
#include "utilities.hpp"

#include "algebra/expansion.tpl.hpp"
#include "algebra/algebra.hpp"
#include "function/function_interface.hpp"
#include "function/polynomial.hpp"
#include "function/function.hpp"
#include "function/procedure.hpp"

#include "function/function_model.hpp"
#include "function/taylor_model.hpp"
#include "function/taylor_function.hpp"

using namespace boost::python;
using namespace Ariadne;

namespace Ariadne {

template<class... AS> inline decltype(auto) _evaluate_(AS... as) { return evaluate(as...); }
template<class... AS> inline decltype(auto) _partial_evaluate_(AS... as) { return partial_evaluate(as...); }
template<class... AS> inline decltype(auto) _unchecked_evaluate_(AS... as) { return unchecked_evaluate(as...); }
template<class... AS> inline decltype(auto) _compose_(AS... as) { return compose(as...); }
template<class... AS> inline decltype(auto) _unchecked_compose_(AS... as) { return unchecked_compose(as...); }

template<class... AS> inline decltype(auto) _join_(AS... as) { return join(as...); }
template<class... AS> inline decltype(auto) _combine_(AS... as) { return combine(as...); }

template<class... AS> inline decltype(auto) _midpoint_(AS... as) { return midpoint(as...); }
template<class... AS> inline decltype(auto) _embed_(AS... as) { return embed(as...); }
template<class... AS> inline decltype(auto) _extension_(AS... as) { return extension(as...); }
template<class... AS> inline decltype(auto) _restriction_(AS... as) { return restriction(as...); }
template<class... AS> inline decltype(auto) _split_(AS... as) { return split(as...); }
template<class... AS> inline decltype(auto) _derivative_(AS... as) { return derivative(as...); }
template<class... AS> inline decltype(auto) _antiderivative_(AS... as) { return antiderivative(as...); }

template<class... AS> inline decltype(auto) _refinement_(AS... as) { return refinement(as...); }
template<class... AS> inline decltype(auto) _refines_(AS... as) { return refines(as...); }
template<class... AS> inline decltype(auto) _inconsistent_(AS... as) { return inconsistent(as...); }



//ValidatedNumericType evaluate(const ValidatedScalarFunctionModel64& f, const Vector<ValidatedNumericType>& x) { return f(x); }
//Vector<ValidatedNumericType> evaluate(const ValidatedVectorFunctionModel64& f, const Vector<ValidatedNumericType>& x) { return f(x); }

ValidatedScalarFunctionModel64 partial_evaluate(const ValidatedScalarFunctionModel64&, SizeType, const ValidatedNumericType&);
ValidatedVectorFunctionModel64 partial_evaluate(const ValidatedVectorFunctionModel64& f, SizeType j, const ValidatedNumericType& c);

ValidatedScalarFunctionModel64 compose(const ValidatedScalarFunctionModel64&, const ValidatedVectorFunctionModel64&);
ValidatedScalarFunctionModel64 compose(const ValidatedScalarFunction&, const ValidatedVectorFunctionModel64&);
ValidatedVectorFunctionModel64 compose(const ValidatedVectorFunctionModel64&, const ValidatedVectorFunctionModel64&);
ValidatedVectorFunctionModel64 compose(const ValidatedVectorFunction&, const ValidatedVectorFunctionModel64&);

ValidatedVectorFunctionModel64 join(const ValidatedScalarFunctionModel64&, const ValidatedScalarFunctionModel64&);
ValidatedVectorFunctionModel64 join(const ValidatedScalarFunctionModel64&, const ValidatedVectorFunctionModel64&);
ValidatedVectorFunctionModel64 join(const ValidatedVectorFunctionModel64&, const ValidatedScalarFunctionModel64&);
ValidatedVectorFunctionModel64 join(const ValidatedVectorFunctionModel64&, const ValidatedVectorFunctionModel64&);

template<class P, class PR, class PRE> OutputStream& operator<<(OutputStream& os, const Representation< ScalarFunctionModel<P,PR,PRE> >& frepr) {
    static_cast<const ScalarFunctionInterface<P>&>(frepr.reference()).repr(os); return os;
}

template<class P, class PR, class PRE> OutputStream& operator<<(OutputStream& os, const Representation< VectorFunctionModel<P,PR,PRE> >& frepr) {
    static_cast<const VectorFunctionInterface<P>&>(frepr.reference()).write(os); return os;
}

ValidatedVectorTaylorFunctionModel64 __getslice__(const ValidatedVectorTaylorFunctionModel64& tf, Int start, Int stop) {
    if(start<0) { start+=tf.result_size(); }
    if(stop<0) { stop+=tf.result_size(); }
    ARIADNE_ASSERT_MSG(0<=start&&start<=stop&&SizeType(stop)<=tf.result_size(),
            "result_size="<<tf.result_size()<<", start="<<start<<", stop="<<stop);
    return ValidatedVectorTaylorFunctionModel64(tf.domain(),Vector<ValidatedTaylorModel64>(project(tf.models(),range(start,stop))));
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
struct from_python<ValidatedVectorTaylorFunctionModel64> {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<ValidatedVectorTaylorFunctionModel64>()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data) {
        Void* storage = ((boost::python::converter::rvalue_from_python_storage<ValidatedVectorTaylorFunctionModel64>*)data)->storage.bytes;
        list lst=boost::python::extract<boost::python::list>(obj_ptr);
        ValidatedVectorTaylorFunctionModel64* tf_ptr = new (storage) ValidatedVectorTaylorFunctionModel64(len(lst));
        for(SizeType i=0; i!=tf_ptr->result_size(); ++i) { tf_ptr->set(i,boost::python::extract<ValidatedScalarTaylorFunctionModel64>(lst[i])); }
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

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Sweeper<Float64>>& repr) {
    const Sweeper<Float64>& swp=repr.reference();
    auto swp_ptr = &static_cast<const SweeperInterface<Float64>&>(swp);
    auto thresh_swp_ptr = dynamic_cast<const ThresholdSweeper<Float64>*>(swp_ptr);
    if(thresh_swp_ptr) {
        os << "ThresholdSweeper64(" << thresh_swp_ptr->sweep_threshold() << ")";
    } else {
        os << swp;
    }
    return os;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ValidatedScalarTaylorFunctionModel64>& repr) {
    const ValidatedScalarTaylorFunctionModel64& stf=repr.reference();
    os << std::setprecision(17);
    os << "ValidatedScalarTaylorFunctionModel64"
       << "(" << python_representation(stf.domain())
       << "," << python_representation(stf.expansion())
       << "," << python_representation(stf.error())
       << "," << python_representation(stf.properties())
       << ")";
    return os;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ValidatedVectorTaylorFunctionModel64>& repr) {
    const ValidatedVectorTaylorFunctionModel64& vtf=repr.reference();
    os << std::setprecision(17);
    os << "ValidatedVectorTaylorFunctionModel64"
       << "(" << python_representation(vtf.domain())
       << "," << python_representation(vtf.expansions())
       << "," << python_representation(vtf.errors())
       << "," << python_representation(vtf.properties())
       << ")";
    return os;
}

List<MultiIndex> keys(const ValidatedTaylorModel64& tm) {
    List<MultiIndex> r;
    for(ValidatedTaylorModel64::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        r.append(iter->key());
    }
    return r;
}

ValidatedScalarFunction unrestrict(const ValidatedScalarFunctionModel64& fm) {
    return ValidatedScalarFunction(fm.raw_pointer()->_clone());
}

ValidatedVectorFunction unrestrict(const ValidatedVectorFunctionModel64& fm) {
    return ValidatedVectorFunction(fm.raw_pointer()->_clone());
}


ExactIntervalType _range1(const ValidatedTaylorModel64&);
ExactIntervalType _range2(const ValidatedTaylorModel64&);
ExactIntervalType _range3(const ValidatedTaylorModel64&);

} // namespace Ariadne

Sweeper<Float64> make_threshold_sweeper(Precision64 pr, double x) { return new ThresholdSweeper<Float64>(pr,x); }
Sweeper<Float64> make_graded_sweeper(Precision64 pr, SizeType n) { return new GradedSweeper<Float64>(pr,n); }

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
    class_<Sweeper<Float64>> sweeper_class("Sweeper64", init<Sweeper<Float64>>());
    def("ThresholdSweeper64", &make_threshold_sweeper );
    def("GradedSweeper64", &make_graded_sweeper );
    sweeper_class.def(self_ns::str(self));

}

/*
Expansion<Float64Value>const& get_expansion(ValidatedTaylorModel64 const& tm) { return tm.expansion(); }

template<class F> Void export_validated_taylor_model()
{
    typedef SizeType SizeType;
    typedef TaylorModel<ValidatedTag,F> ValidatedTaylorModelType;

    class_<ValidatedTaylorModel64> taylor_model_class("ValidatedTaylorModel64", init<ValidatedTaylorModel64>());
    taylor_model_class.def( init< SizeType,Sweeper64 >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ValidatedTaylorModel64&))&keys);
    taylor_model_class.def("value", (const Float64Value&(ValidatedTaylorModel64::*)()const) &ValidatedTaylorModel64::value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("gradient", (const Float64Value&(ValidatedTaylorModel64::*)(SizeType)const) &ValidatedTaylorModel64::gradient_value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("error", (const Float64Error&(ValidatedTaylorModel64::*)()const) &ValidatedTaylorModel64::error, return_value_policy<copy_const_reference>());
    taylor_model_class.def("expansion", (const Expansion<Float64Value>&(*)(ValidatedTaylorModel64 const&)) &get_expansion, return_value_policy<copy_const_reference>());
    taylor_model_class.def("set_error", (Void(ValidatedTaylorModel64::*)(const Float64Error&)) &ValidatedTaylorModel64::set_error);
    taylor_model_class.def("argument_size", &ValidatedTaylorModel64::argument_size);
    taylor_model_class.def("domain", &ValidatedTaylorModel64::domain);
    taylor_model_class.def("range", &ValidatedTaylorModel64::range);
    taylor_model_class.def("set_sweeper", &ValidatedTaylorModel64::set_sweeper);
    taylor_model_class.def("sweeper", &ValidatedTaylorModel64::sweeper);
    taylor_model_class.def("sweep", (ValidatedTaylorModel64&(ValidatedTaylorModel64::*)()) &ValidatedTaylorModel64::sweep, return_value_policy<reference_existing_object>());
    taylor_model_class.def("__getitem__", &__getitem__<ValidatedTaylorModel64,MultiIndex,Float64Value>);
    taylor_model_class.def("__setitem__",&__setitem__<ValidatedTaylorModel64,MultiIndex,Float64Value>);
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

    taylor_model_class.def("constant",(ValidatedTaylorModel64(*)(SizeType, const ValidatedNumericType&,Sweeper64))&ValidatedTaylorModel64::constant);
    taylor_model_class.def("coordinate",(ValidatedTaylorModel64(*)(SizeType, SizeType,Sweeper64))&ValidatedTaylorModel64::coordinate);

    taylor_model_class.staticmethod("constant");
    taylor_model_class.staticmethod("coordinate");

    //def("max",(ValidatedTaylorModel64(*)(const ValidatedTaylorModel64&,const ValidatedTaylorModel64&))&max);
    //def("min",(ValidatedTaylorModel64(*)(const ValidatedTaylorModel64&,const ValidatedTaylorModel64&))&min);
    //def("abs",(ValidatedTaylorModel64(*)(const ValidatedTaylorModel64&))&abs);

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

    taylor_model_class.def("range", (UpperIntervalType(ValidatedTaylorModel64::*)()const) &ValidatedTaylorModel64::range);

    //def("evaluate", (ValidatedNumericType(*)(const ValidatedTaylorModel64&, const Vector<ValidatedNumericType>&))&evaluate);
    //def("split",(ValidatedTaylorModel64(*)(const ValidatedTaylorModel64&,SizeType,SplitPart)) &split);

    from_python< Vector<ValidatedTaylorModel64> >();
    to_python< Vector<ValidatedTaylorModel64> >();

}

Void export_approximate_taylor_model()
{
    typedef SizeType SizeType;
    typedef ApproximateTaylorModel64 ApproximateTaylorModel64;

    class_<ApproximateTaylorModel64> taylor_model_class("ApproximateTaylorModel64", init<ApproximateTaylorModel64>());
    taylor_model_class.def( init< SizeType,Sweeper64 >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ApproximateTaylorModel64&))&keys);
    taylor_model_class.def("value", (const Float64Approximation&(ApproximateTaylorModel64::*)()const) &ApproximateTaylorModel64::value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("gradient", (const Float64Approximation&(ApproximateTaylorModel64::*)(SizeType)const) &ApproximateTaylorModel64::gradient_value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("expansion", (const Expansion<Float64Approximation>&(*)(ApproximateTaylorModel64 const&)) &get_expansion, return_value_policy<copy_const_reference>());
    taylor_model_class.def("argument_size", &ApproximateTaylorModel64::argument_size);
    taylor_model_class.def("domain", &ApproximateTaylorModel64::domain);
    taylor_model_class.def("range", &ApproximateTaylorModel64::range);
    taylor_model_class.def("set_sweeper", &ApproximateTaylorModel64::set_sweeper);
    taylor_model_class.def("sweeper", &ApproximateTaylorModel64::sweeper);
    taylor_model_class.def("sweep", (ApproximateTaylorModel64&(ApproximateTaylorModel64::*)()) &ApproximateTaylorModel64::sweep, return_value_policy<reference_existing_object>());
    taylor_model_class.def("__getitem__", &__getitem__<ApproximateTaylorModel64,MultiIndex,Float64Approximation>);
    taylor_model_class.def("__setitem__",&__setitem__<ApproximateTaylorModel64,MultiIndex,Float64Approximation>);
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

    taylor_model_class.def("constant",(ApproximateTaylorModel64(*)(SizeType, const ApproximateNumericType&,Sweeper64))&ApproximateTaylorModel64::constant);
    taylor_model_class.def("coordinate",(ApproximateTaylorModel64(*)(SizeType, SizeType,Sweeper64))&ApproximateTaylorModel64::coordinate);

    taylor_model_class.staticmethod("constant");
    taylor_model_class.staticmethod("coordinate");

    //def("max",(ApproximateTaylorModel64(*)(const ApproximateTaylorModel64&,const ApproximateTaylorModel64&))&max);
    //def("min",(ApproximateTaylorModel64(*)(const ApproximateTaylorModel64&,const ApproximateTaylorModel64&))&min);
    //def("abs",(ApproximateTaylorModel64(*)(const ApproximateTaylorModel64&))&abs);

    typedef AlgebraOperations<ApproximateTaylorModel64> Operations;
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

    //def("evaluate", (ApproximateNumericType(*)(const ApproximateTaylorModel64&, const Vector<ApproximateNumericType>&))&evaluate);
    //def("split",(ApproximateTaylorModel64(*)(const ApproximateTaylorModel64&,SizeType,SplitPart)) &split);

    from_python< Vector<ApproximateTaylorModel64> >();
    to_python< Vector<ApproximateTaylorModel64> >();
}

*/

/*

Void export_scalar_function_model()
{
    class_<ValidatedScalarFunctionModel64> scalar_function_model_class("ValidatedScalarFunctionModel64",init<ValidatedScalarFunctionModel64>());
    scalar_function_model_class.def(init<ValidatedScalarTaylorFunctionModel64>());
    scalar_function_model_class.def("argument_size", &ValidatedScalarFunctionModel64::argument_size);
    scalar_function_model_class.def("domain", &ValidatedScalarFunctionModel64::domain);
    scalar_function_model_class.def("codomain", &ValidatedScalarFunctionModel64::codomain);
    scalar_function_model_class.def("range", &ValidatedScalarFunctionModel64::range);
    scalar_function_model_class.def("clobber", &ValidatedScalarFunctionModel64::clobber);
    scalar_function_model_class.def("error", &ValidatedScalarFunctionModel64::error);
    scalar_function_model_class.def("__call__", (Float64Bounds(ValidatedScalarFunctionModel64::*)(const Vector<Float64Bounds>&)const) &ValidatedScalarFunctionModel64::operator());
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
    scalar_function_model_class.def("__str__", &_cstr_<ValidatedScalarFunctionModel64>);
    scalar_function_model_class.def("__repr__", &_crepr_<ValidatedScalarFunctionModel64>);
    //scalar_function_model_class.def("__repr__",&__repr__<ValidatedScalarFunctionModel64>);

    def("evaluate", (ValidatedNumericType(*)(const ValidatedScalarFunctionModel64&,const Vector<ValidatedNumericType>&)) &evaluate);
    def("partial_evaluate", (ValidatedScalarFunctionModel64(*)(const ValidatedScalarFunctionModel64&,SizeType,const ValidatedNumericType&)) &partial_evaluate);

    def("compose", (ValidatedScalarFunctionModel64(*)(const ValidatedScalarFunctionModel64&, const ValidatedVectorFunctionModel64&)) &compose);
    def("compose", (ValidatedScalarFunctionModel64(*)(const ValidatedScalarFunction&, const ValidatedVectorFunctionModel64&)) &compose);

    def("unrestrict", (ValidatedScalarFunction(*)(const ValidatedScalarFunctionModel64&)) &unrestrict);


}

Void export_vector_function_model()
{
    class_<ValidatedVectorFunctionModel64> vector_function_model_class("ValidatedVectorFunctionModel64",init<ValidatedVectorFunctionModel64>());
    vector_function_model_class.def(init<ValidatedVectorTaylorFunctionModel64>());
    vector_function_model_class.def("result_size", &ValidatedVectorFunctionModel64::result_size);
    vector_function_model_class.def("argument_size", &ValidatedVectorFunctionModel64::argument_size);
    vector_function_model_class.def("domain", &ValidatedVectorFunctionModel64::domain);
    vector_function_model_class.def("codomain", &ValidatedVectorFunctionModel64::codomain);
    vector_function_model_class.def("range", &ValidatedVectorFunctionModel64::range);
    //vector_function_model_class.def("__getslice__", (ValidatedVectorTaylorFunctionModel64(*)(const ValidatedVectorTaylorFunctionModel64&,Int,Int))&__getslice__);
    vector_function_model_class.def("__getitem__", &__getitem__<ValidatedVectorFunctionModel64,SizeType,ValidatedScalarFunctionModel64>);
    vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorFunctionModel64,SizeType,ValidatedScalarFunctionModel64>);
    //vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorFunctionModel64,SizeType,ValidatedScalarFunction>);
    vector_function_model_class.def("__call__", (Vector<Float64Bounds>(ValidatedVectorFunctionModel64::*)(const Vector<Float64Bounds>&)const) &ValidatedVectorFunctionModel64::operator());
    vector_function_model_class.def(self*ValidatedNumericType());
    vector_function_model_class.def("__str__", &_cstr_<ValidatedVectorFunctionModel64>);
    vector_function_model_class.def("__repr__", &_crepr_<ValidatedVectorFunctionModel64>);
    //export_vector_function_model.def("__repr__",&__repr__<ValidatedVectorFunctionModel64>);

    def("evaluate", (Vector<ValidatedNumericType>(*)(const ValidatedVectorFunctionModel64&,const Vector<ValidatedNumericType>&)) &evaluate);

    def("compose", (ValidatedVectorFunctionModel64(*)(const ValidatedVectorFunctionModel64&,const ValidatedVectorFunctionModel64&)) &compose);
    def("compose", (ValidatedVectorFunctionModel64(*)(const ValidatedVectorFunction&,const ValidatedVectorFunctionModel64&)) &compose);

    def("unrestrict", (ValidatedVectorFunction(*)(const ValidatedVectorFunctionModel64&)) &unrestrict);

    def("join", (ValidatedVectorFunctionModel64(*)(const ValidatedScalarFunctionModel64&,const ValidatedScalarFunctionModel64&)) &join);
    def("join", (ValidatedVectorFunctionModel64(*)(const ValidatedScalarFunctionModel64&,const ValidatedVectorFunctionModel64&)) &join);
    def("join", (ValidatedVectorFunctionModel64(*)(const ValidatedVectorFunctionModel64&,const ValidatedScalarFunctionModel64&)) &join);
    def("join", (ValidatedVectorFunctionModel64(*)(const ValidatedVectorFunctionModel64&,const ValidatedVectorFunctionModel64&)) &join);

    to_python< List<ValidatedVectorFunctionModel64> >();
}

*/

Void export_scalar_taylor_function()
{
    typedef ValidatedScalarTaylorFunctionModel64 F;
    typedef ValidatedVectorTaylorFunctionModel64 VF;
    typedef ValidatedTaylorModel64 M;
    typedef typename F::DomainType D;
    typedef typename F::NumericType X;
    typedef Vector<X> VX;
    typedef SizeType I;
    typedef ValidatedNumericType Y;
    typedef Vector<Y> VY;

    class_<ValidatedScalarTaylorFunctionModel64> scalar_taylor_function_class("ValidatedScalarTaylorFunctionModel64",init<ValidatedScalarTaylorFunctionModel64>());
    scalar_taylor_function_class.def(init<ExactBoxType,ValidatedTaylorModel64>());
    scalar_taylor_function_class.def(init< ExactBoxType,Sweeper64 >());
    scalar_taylor_function_class.def(init< ExactBoxType, const EffectiveScalarFunction&,Sweeper64 >());
    scalar_taylor_function_class.def(init< ExactBoxType, Expansion<Float64Value>, Float64Error, Sweeper64 >());
    scalar_taylor_function_class.def("error", (const Float64Error&(ValidatedScalarTaylorFunctionModel64::*)()const) &ValidatedScalarTaylorFunctionModel64::error, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("set_error", (Void(ValidatedScalarTaylorFunctionModel64::*)(const Float64Error&)) &ValidatedScalarTaylorFunctionModel64::set_error);
    scalar_taylor_function_class.def("argument_size", &F::argument_size);
    scalar_taylor_function_class.def("domain", &F::domain);
    scalar_taylor_function_class.def("codomain", &F::codomain);
    scalar_taylor_function_class.def("range", &F::range);
    scalar_taylor_function_class.def("model", (const ValidatedTaylorModel64&(ValidatedScalarTaylorFunctionModel64::*)()const)&ValidatedScalarTaylorFunctionModel64::model, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("polynomial", (Polynomial<ExactIntervalType>(ValidatedScalarTaylorFunctionModel64::*)()const)&ValidatedScalarTaylorFunctionModel64::polynomial);
    scalar_taylor_function_class.def("number_of_nonzeros", (SizeType(ValidatedScalarTaylorFunctionModel64::*)()const)&ValidatedScalarTaylorFunctionModel64::number_of_nonzeros);
//    scalar_taylor_function_class.def("set_sweeper", &ValidatedScalarTaylorFunctionModel64::set_sweeper);
//    scalar_taylor_function_class.def("sweeper", &ValidatedScalarTaylorFunctionModel64::sweeper);
//    scalar_taylor_function_class.def("sweep", (ValidatedScalarTaylorFunctionModel64&(ValidatedScalarTaylorFunctionModel64::*)()) &ValidatedScalarTaylorFunctionModel64::sweep, return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("__getitem__", &__getitem__<ValidatedScalarTaylorFunctionModel64,MultiIndex,Float64Value>);
    scalar_taylor_function_class.def("__setitem__",&__setitem__<ValidatedScalarTaylorFunctionModel64,MultiIndex,Float64Value>);
    scalar_taylor_function_class.def(+self);
    scalar_taylor_function_class.def(-self);
    scalar_taylor_function_class.def(self+self);
    scalar_taylor_function_class.def(self-self);
    scalar_taylor_function_class.def(self*self);
    scalar_taylor_function_class.def(self/self);
    scalar_taylor_function_class.def(self+Y());
    scalar_taylor_function_class.def(self-Y());
    scalar_taylor_function_class.def(self*Y());
    scalar_taylor_function_class.def(self/Y());
    scalar_taylor_function_class.def(Y()+self);
    scalar_taylor_function_class.def(Y()-self);
    scalar_taylor_function_class.def(Y()*self);
    scalar_taylor_function_class.def(Y()/self);
    scalar_taylor_function_class.def(self+=Y());
    scalar_taylor_function_class.def(self-=Y());
    scalar_taylor_function_class.def(self*=Y());
    scalar_taylor_function_class.def(self/=Y());
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
    scalar_taylor_function_class.def("__str__", &__cstr__<F>);
    scalar_taylor_function_class.def("__repr__", &__crepr__<F>);
    scalar_taylor_function_class.def("__mul__",&__mul__< VF, F, VY >);

    //scalar_taylor_function_class.def("__str__",(StringType(*)(const ValidatedScalarTaylorFunctionModel64&)) &__str__);
    //scalar_taylor_function_class.def("_cstr_",(StringType(*)(const ValidatedScalarTaylorFunctionModel64&)) &_cstr_);
    //scalar_taylor_function_class.def("__repr__",(StringType(*)(const ValidatedScalarTaylorFunctionModel64&)) &__repr__);
    scalar_taylor_function_class.def("value", (const Float64Value&(ValidatedScalarTaylorFunctionModel64::*)()const) &ValidatedScalarTaylorFunctionModel64::value,return_value_policy<copy_const_reference>());
//    scalar_taylor_function_class.def("sweep", (ValidatedScalarTaylorFunctionModel64&(ValidatedScalarTaylorFunctionModel64::*)())&ValidatedScalarTaylorFunctionModel64::sweep,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("clobber", (ValidatedScalarTaylorFunctionModel64&(ValidatedScalarTaylorFunctionModel64::*)()) &ValidatedScalarTaylorFunctionModel64::clobber,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("set_properties",&ValidatedScalarTaylorFunctionModel64::set_properties);
    scalar_taylor_function_class.def("properties",&ValidatedScalarTaylorFunctionModel64::properties);
    scalar_taylor_function_class.def("__call__", (Float64Approximation(ValidatedScalarTaylorFunctionModel64::*)(const Vector<Float64Approximation>&)const) &ValidatedScalarTaylorFunctionModel64::operator());
    scalar_taylor_function_class.def("__call__", (Float64Bounds(ValidatedScalarTaylorFunctionModel64::*)(const Vector<Float64Bounds>&)const) &ValidatedScalarTaylorFunctionModel64::operator());
    scalar_taylor_function_class.def("gradient", (Covector<Float64Bounds>(ValidatedScalarTaylorFunctionModel64::*)(const Vector<Float64Bounds>&)const) &ValidatedScalarTaylorFunctionModel64::gradient);
    scalar_taylor_function_class.def("function", (EffectiveScalarFunction(ValidatedScalarTaylorFunctionModel64::*)()const) &ValidatedScalarTaylorFunctionModel64::function);
    scalar_taylor_function_class.def("polynomial", (Polynomial<Float64Bounds>(ValidatedScalarTaylorFunctionModel64::*)()const) &ValidatedScalarTaylorFunctionModel64::polynomial);
    scalar_taylor_function_class.def("restriction",&_restriction_<F,D>);
//    scalar_taylor_function_class.def("extension",&_extension_<F,D>);

    scalar_taylor_function_class.def("zero",(ValidatedScalarTaylorFunctionModel64(*)(const ExactBoxType&,Sweeper64))&ValidatedScalarTaylorFunctionModel64::zero);
    scalar_taylor_function_class.def("constant",(ValidatedScalarTaylorFunctionModel64(*)(const ExactBoxType&,const ValidatedNumericType&,Sweeper64))&ValidatedScalarTaylorFunctionModel64::constant);
    scalar_taylor_function_class.def("coordinate",(ValidatedScalarTaylorFunctionModel64(*)(const ExactBoxType&,SizeType,Sweeper64))&ValidatedScalarTaylorFunctionModel64::coordinate);


    scalar_taylor_function_class.staticmethod("constant");
    scalar_taylor_function_class.staticmethod("coordinate");

    def("restriction",&_restriction_<F,D>);
//    def("extension",&_extension_<F,D>);
    def("embed",&_embed_<D,F,D>);
    def("split",&_split_<F,I>);
    def("evaluate",&_evaluate_<F,VX>);
    def("evaluate",&_partial_evaluate_<F,I,X>);
    def("midpoint",&_midpoint_<F>);
    def("derivative",&_derivative_<F,I>);
    def("antiderivative",&_antiderivative_<F,I>);
    def("antiderivative",&_antiderivative_<F,I,X>);

    def("inconsistent",&_inconsistent_<F,F>);
    def("refines",&_refines_<F,F>);
    def("refinement",&_refinement_<F,F>);

    def("max",&_max_<F,F>); def("min",&_min_<F,F>); def("abs",&_abs_<F>);

    def("neg",&_neg_<F>); def("rec",&_rec_<F>); def("sqr",&_sqr_<F>); def("pow",&_pow_<F,Int>);
    def("sqrt",&_sqrt_<F>); def("exp",&_exp_<F>); def("log",&_log_<F>); def("atan",&_atan_<F>);
    def("sin",&_sin_<F>); def("cos",&_cos_<F>); def("tan",&_tan_<F>);

    to_python< Vector<ValidatedScalarTaylorFunctionModel64> >();
}

Void export_vector_taylor_function()
{
    typedef SizeType I;
    typedef ValidatedScalarFunction SFN;
    typedef ValidatedVectorFunction VFN;
    typedef ValidatedScalarTaylorFunctionModel64 SF;
    typedef ValidatedVectorTaylorFunctionModel64 VF;
    typedef typename VF::DomainType D;
    typedef typename D::ScalarType Di;
    typedef typename VF::NumericType X;
    typedef Vector<X> VX;

    class_<ValidatedVectorTaylorFunctionModel64> vector_taylor_function_class("ValidatedVectorTaylorFunctionModel64", init<ValidatedVectorTaylorFunctionModel64>());
    vector_taylor_function_class.def( init< SizeType, ExactBoxType, Sweeper64 >());
    vector_taylor_function_class.def( init< ExactBoxType,const EffectiveVectorFunction&,Sweeper64 >());
    vector_taylor_function_class.def(init< ExactBoxType, Vector< Expansion<Float64Value> >, Vector<Float64Error>, Sweeper64 >());
    vector_taylor_function_class.def( init< Vector<ValidatedScalarTaylorFunctionModel64> >());
    vector_taylor_function_class.def("_len_", &ValidatedVectorTaylorFunctionModel64::result_size);
    vector_taylor_function_class.def("result_size", &ValidatedVectorTaylorFunctionModel64::result_size);
    vector_taylor_function_class.def("argument_size", &ValidatedVectorTaylorFunctionModel64::argument_size);
    vector_taylor_function_class.def("domain", &ValidatedVectorTaylorFunctionModel64::domain);
    vector_taylor_function_class.def("codomain", &ValidatedVectorTaylorFunctionModel64::codomain);
    // FIXME: Omitted since const and non-const versions
    // vector_taylor_function_class.def("models", &ValidatedVectorTaylorFunctionModel64::models, return_value_policy<copy_const_reference>());
    vector_taylor_function_class.def("centre", &ValidatedVectorTaylorFunctionModel64::centre);
    vector_taylor_function_class.def("range", &ValidatedVectorTaylorFunctionModel64::range);
    vector_taylor_function_class.def("errors", &ValidatedVectorTaylorFunctionModel64::errors);
//    vector_taylor_function_class.def("sweep", (ValidatedVectorTaylorFunctionModel64&(ValidatedVectorTaylorFunctionModel64::*)())&ValidatedVectorTaylorFunctionModel64::sweep,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("clobber", (ValidatedVectorTaylorFunctionModel64&(ValidatedVectorTaylorFunctionModel64::*)()) &ValidatedVectorTaylorFunctionModel64::clobber,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("set_properties",&ValidatedVectorTaylorFunctionModel64::set_properties);
    vector_taylor_function_class.def("properties",&ValidatedVectorTaylorFunctionModel64::properties);
//    vector_taylor_function_class.def("sweep", (ValidatedVectorTaylorFunctionModel64&(ValidatedVectorTaylorFunctionModel64::*)()) &ValidatedVectorTaylorFunctionModel64::sweep, return_value_policy<reference_existing_object>());
    //vector_taylor_function_class.def("__getslice__", &__getslice__<ValidatedVectorTaylorFunctionModel64,SizeType,SizeType,ValidatedScalarTaylorFunctionModel64>);
    vector_taylor_function_class.def("__getslice__", (ValidatedVectorTaylorFunctionModel64(*)(const ValidatedVectorTaylorFunctionModel64&,Int,Int))&__getslice__);
    vector_taylor_function_class.def("__getitem__", &__getitem__<ValidatedVectorTaylorFunctionModel64,SizeType,ValidatedScalarTaylorFunctionModel64>);
    vector_taylor_function_class.def("__setitem__",&__setitem__<ValidatedVectorTaylorFunctionModel64,SizeType,ValidatedScalarTaylorFunctionModel64>);
    vector_taylor_function_class.def(-self);
    vector_taylor_function_class.def(self+self);
    vector_taylor_function_class.def(self-self);
    vector_taylor_function_class.def(self+Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self-Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self*ValidatedNumericType());
    vector_taylor_function_class.def(self/ValidatedNumericType());
//    vector_taylor_function_class.def(self*ValidatedScalarTaylorFunctionModel64());
    vector_taylor_function_class.def(self+=Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self-=Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self*=ValidatedNumericType());
    vector_taylor_function_class.def(self/=ValidatedNumericType());
    vector_taylor_function_class.def(self+=self);
    vector_taylor_function_class.def(self-=self);
    vector_taylor_function_class.def("__str__", &__cstr__<ValidatedVectorTaylorFunctionModel64>);
    vector_taylor_function_class.def("__repr__", &__crepr__<ValidatedVectorTaylorFunctionModel64>);
    vector_taylor_function_class.def("clobber", (ValidatedVectorTaylorFunctionModel64&(ValidatedVectorTaylorFunctionModel64::*)()) &ValidatedVectorTaylorFunctionModel64::clobber,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("__call__", (Vector<Float64Approximation>(ValidatedVectorTaylorFunctionModel64::*)(const Vector<Float64Approximation>&)const) &ValidatedVectorTaylorFunctionModel64::operator());
    vector_taylor_function_class.def("__call__", (Vector<Float64Bounds>(ValidatedVectorTaylorFunctionModel64::*)(const Vector<Float64Bounds>&)const) &ValidatedVectorTaylorFunctionModel64::operator());
     //vector_taylor_function_class.def("jacobian", (Vector<Float64Bounds>(ValidatedVectorTaylorFunctionModel64::*)(const Vector<Float64Bounds>&)const) &ValidatedVectorTaylorFunctionModel64::jacobian);
    vector_taylor_function_class.def("polynomials", (Vector< Polynomial<Float64Bounds> >(ValidatedVectorTaylorFunctionModel64::*)()const) &ValidatedVectorTaylorFunctionModel64::polynomials);
    vector_taylor_function_class.def("function", (EffectiveVectorFunction(ValidatedVectorTaylorFunctionModel64::*)()const) &ValidatedVectorTaylorFunctionModel64::function);


    vector_taylor_function_class.def("constant",(ValidatedVectorTaylorFunctionModel64(*)(const ExactBoxType&, const Vector<ValidatedNumericType>&,Sweeper64))&ValidatedVectorTaylorFunctionModel64::constant);
    vector_taylor_function_class.def("identity",(ValidatedVectorTaylorFunctionModel64(*)(const ExactBoxType&,Sweeper64))&ValidatedVectorTaylorFunctionModel64::identity);

    vector_taylor_function_class.staticmethod("constant");
    vector_taylor_function_class.staticmethod("identity");

    def("inconsistent", &_inconsistent_<VF,VF>);
    def("refinement", &_refinement_<VF,VF>);
    def("refines", &_refines_<VF,VF>);

    def("join", &_join_<VF,VF>); def("join", &_join_<VF,SF>); def("join", &_join_<SF,VF>); // def("join", &_join_<SF,SF>);
    def("combine", &_combine_<VF,VF>); def("combine", &_combine_<VF,SF>); def("combine", &_combine_<SF,VF>); //def("combine", &_combine_<SF,SF>);
    def("embed", &_embed_<VF,Di>); def("embed", &_embed_<VF,D>); def("embed", &_embed_<D,VF>); def("embed", &_embed_<D,VF,D>);

    def("restriction", &_restriction_<VF,D>); def("restriction", &_restriction_<VF,I,Di>);
//    def("split", &_split_<VF,I>);

    def("evaluate", &_evaluate_<VF,VX>);
    def("partial_evaluate", &_partial_evaluate_<VF,I,X>);
    def("compose", &_compose_<VF,VF>);
    def("compose", &_compose_<SF,VF>);
    def("compose", &_compose_<SFN,VF>);
    def("compose", &_compose_<VFN,VF>);
    def("unchecked_compose", &_compose_<SF,VF>);
    def("unchecked_compose", &_compose_<VF,VF>);
    def("derivative", &_derivative_<VF,I>);
    def("antiderivative", &_antiderivative_<VF,I>);
    def("antiderivative", &_antiderivative_<VF,I,X>);

//    def("evaluate",(Vector<Float64Approximation>(ValidatedVectorTaylorFunctionModel64::*)(const Vector<Float64Approximation>&)const) &ValidatedVectorTaylorFunctionModel64::evaluate);

    from_python<ValidatedVectorTaylorFunctionModel64>();
}

Void calculus_submodule()
{
    export_expansion();
    export_sweeper();
//    export_approximate_taylor_model();
//    export_validated_taylor_model();
//    export_scalar_function_model();
//    export_vector_function_model();
    export_scalar_taylor_function();
    export_vector_taylor_function();
}


