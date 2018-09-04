/***************************************************************************
 *            calculus_submodule.cpp
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

#include <type_traits>

#include "pybind11.hpp"
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

using namespace Ariadne;

namespace Ariadne {



//ValidatedNumericType evaluate(const ValidatedScalarFunctionModelDP& f, const Vector<ValidatedNumericType>& x) { return f(x); }
//Vector<ValidatedNumericType> evaluate(const ValidatedVectorFunctionModelDP& f, const Vector<ValidatedNumericType>& x) { return f(x); }

ValidatedScalarFunctionModelDP partial_evaluate(const ValidatedScalarFunctionModelDP&, SizeType, const ValidatedNumericType&);
ValidatedVectorFunctionModelDP partial_evaluate(const ValidatedVectorFunctionModelDP&, SizeType, const ValidatedNumericType&);

//ValidatedScalarFunctionModelDP antiderivative(const ValidatedScalarFunctionModelDP&, SizeType, ValidatedNumericType);
//ValidatedVectorFunctionModelDP antiderivative(const ValidatedVectorFunctionModelDP&, SizeType, ValidatedNumericType);

ValidatedScalarFunctionModelDP compose(const ValidatedScalarFunctionModelDP&, const ValidatedVectorFunctionModelDP&);
ValidatedScalarFunctionModelDP compose(const ValidatedScalarFunction&, const ValidatedVectorFunctionModelDP&);
ValidatedVectorFunctionModelDP compose(const ValidatedVectorFunctionModelDP&, const ValidatedVectorFunctionModelDP&);
ValidatedVectorFunctionModelDP compose(const ValidatedVectorFunction&, const ValidatedVectorFunctionModelDP&);

ValidatedVectorFunctionModelDP join(const ValidatedScalarFunctionModelDP&, const ValidatedScalarFunctionModelDP&);
ValidatedVectorFunctionModelDP join(const ValidatedScalarFunctionModelDP&, const ValidatedVectorFunctionModelDP&);
ValidatedVectorFunctionModelDP join(const ValidatedVectorFunctionModelDP&, const ValidatedScalarFunctionModelDP&);
ValidatedVectorFunctionModelDP join(const ValidatedVectorFunctionModelDP&, const ValidatedVectorFunctionModelDP&);

template<class P, class PR, class PRE> OutputStream& operator<<(OutputStream& os, const Representation< ScalarFunctionModel<P,PR,PRE> >& frepr) {
    static_cast<const ScalarFunctionInterface<P>&>(frepr.reference()).repr(os); return os;
}

template<class P, class PR, class PRE> OutputStream& operator<<(OutputStream& os, const Representation< VectorFunctionModel<P,PR,PRE> >& frepr) {
    static_cast<const VectorFunctionInterface<P>&>(frepr.reference()).write(os); return os;
}

ValidatedVectorTaylorFunctionModelDP __getslice__(const ValidatedVectorTaylorFunctionModelDP& tf, Int start, Int stop) {
    Int rs = tf.result_size();
    if(start<0) { start+=rs; }
    if(stop<0) { stop+=rs; }
    ARIADNE_ASSERT_MSG(0<=start&&start<=stop&&stop<=rs,
            "result_size="<<rs<<", start="<<start<<", stop="<<stop);
    return ValidatedVectorTaylorFunctionModelDP(tf.domain(),Vector<ValidatedTaylorModelDP>(project(tf.models(),range(static_cast<SizeType>(start),static_cast<SizeType>(stop)))));
}

/*

template<>
struct from_python<MultiIndex> {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<MultiIndex>()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* coefficient) {
        Void* storage = ((boost::python::converter::rvalue_from_python_storage<MultiIndex>*)coefficient)->storage.bytes;
        boost::python::extract<boost::python::tuple> xtup(obj_ptr);
        boost::python::extract<boost::python::list> xlst(obj_ptr);
        if(xlst.check()) {
            MultiIndex a(static_cast<SizeType>(len(xlst))); for(SizeType i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<SizeType>(xlst()[i]); }
            new (storage) MultiIndex(a);
        } else {
            MultiIndex a(static_cast<SizeType>(len(xtup))); for(SizeType i=0; i!=a.size(); ++i) { (a)[i]=boost::python::extract<SizeType>(xtup()[i]); }
            new (storage) MultiIndex(a);
        }
        coefficient->convertible = storage;
    }
};


template<class T>
struct from_python< Expansion<MultiIndex,T> > {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< Expansion<MultiIndex,T> >()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyDict_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* coefficient) {
        Void* storage = ((boost::python::converter::rvalue_from_python_storage< Expansion<MultiIndex,T> >*)coefficient)->storage.bytes;
        Expansion<MultiIndex,T> r(0);
        boost::python::dict dct=boost::python::extract<boost::python::dict>(obj_ptr);
        boost::python::list lst=dct.items();
        MultiIndex a;
        if(len(lst)!=0) {
            boost::python::tuple tup=boost::python::extract<boost::python::tuple>(lst[0]);
            a=boost::python::extract<MultiIndex>(tup[0]);
            r=Expansion<MultiIndex,T>(a.size());
            r.reserve(static_cast<SizeType>(len(lst)));
        }
        for(Int i=0; i!=len(lst); ++i) {
            boost::python::tuple tup=boost::python::extract<boost::python::tuple>(lst[i]);
            MultiIndex a=boost::python::extract<MultiIndex>(tup[0]);
            T c=boost::python::extract<T>(tup[1]);
            r.append(a,c);
        }
        new (storage) Expansion<MultiIndex,T>(r);
        //r.unique_sort();
        coefficient->convertible = storage;
    }
};


template<class X>
struct from_python< Vector<X> >
{
    from_python() { boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< Vector<X> >()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* coefficient) {
        boost::python::list lst=boost::python::extract<boost::python::list>(obj_ptr);
        Void* storage = ((boost::python::converter::rvalue_from_python_storage< Vector<X> >*) coefficient)->storage.bytes;
        Array<X> ary(static_cast<SizeType>(len(lst)),Uninitialised());
        for(SizeType i=0; i!=ary.size(); ++i) { new (&ary[i]) X(boost::python::extract<X>(lst[i])); }
        new (storage) Vector<X>(std::move(ary));
        coefficient->convertible = storage;
    }
};


template<>
struct from_python<ValidatedVectorTaylorFunctionModelDP> {
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<ValidatedVectorTaylorFunctionModelDP>()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) && !PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* coefficient) {
        Void* storage = ((boost::python::converter::rvalue_from_python_storage<ValidatedVectorTaylorFunctionModelDP>*)coefficient)->storage.bytes;
        boost::python::list lst=boost::python::extract<boost::python::list>(obj_ptr);
        ValidatedVectorTaylorFunctionModelDP* tf_ptr = new (storage) ValidatedVectorTaylorFunctionModelDP(static_cast<SizeType>(len(lst)));
        for(SizeType i=0; i!=tf_ptr->result_size(); ++i) { tf_ptr->set(i,boost::python::extract<ValidatedScalarTaylorFunctionModelDP>(lst[i])); }
        coefficient->convertible = storage;
    }
};

*/


template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< X >& repr) {
    return os << repr.reference();
}

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Expansion<MultiIndex,X> >& repr) {
    const Expansion<MultiIndex,X>& exp=repr.reference();
    for(typename Expansion<MultiIndex,X>::ConstIterator iter=exp.begin(); iter!=exp.end(); ++iter) {
        os << (iter==exp.begin()?'{':',') << "(";
        for(SizeType j=0; j!=iter->index().size(); ++j) {
            if(j!=0) { os << ','; } os << Int(iter->index()[j]);
        }
        os << "):" << python_representation(iter->coefficient());
    }
    os << "}";
    return os;
}

template<class X, class CMP> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< SortedExpansion<MultiIndex,X,CMP> >& repr) {
    return os << python_representation(static_cast<const Expansion<MultiIndex,X>&>(repr.reference()));
}

template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<MultiIndex,RawFloatDP> >&);
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<MultiIndex,FloatDPApproximation> >&);
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<MultiIndex,FloatDPBounds> >&);
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<MultiIndex,FloatDPValue> >&);

template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<MultiIndex,RawFloatMP> >&);
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<MultiIndex,FloatMPApproximation> >&);
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<MultiIndex,FloatMPBounds> >&);
template OutputStream& operator<<(OutputStream&, const PythonRepresentation< Expansion<MultiIndex,FloatMPValue> >&);


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

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Sweeper<FloatDP>>& repr) {
    const Sweeper<FloatDP>& swp=repr.reference();
    auto swp_ptr = &static_cast<const SweeperInterface<FloatDP>&>(swp);
    auto thresh_swp_ptr = dynamic_cast<const ThresholdSweeper<FloatDP>*>(swp_ptr);
    if(thresh_swp_ptr) {
        os << "ThresholdSweeperDP(" << thresh_swp_ptr->sweep_threshold() << ")";
    } else {
        os << swp;
    }
    return os;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ValidatedScalarTaylorFunctionModelDP>& repr) {
    const ValidatedScalarTaylorFunctionModelDP& stf=repr.reference();
    os << std::setprecision(17);
    os << "ValidatedScalarTaylorFunctionModelDP"
       << "(" << python_representation(stf.domain())
       << "," << python_representation(stf.expansion())
       << "," << python_representation(stf.error())
       << "," << python_representation(stf.properties())
       << ")";
    return os;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ValidatedVectorTaylorFunctionModelDP>& repr) {
    const ValidatedVectorTaylorFunctionModelDP& vtf=repr.reference();
    os << std::setprecision(17);
    os << "ValidatedVectorTaylorFunctionModelDP"
       << "(" << python_representation(vtf.domain())
       << "," << python_representation(vtf.expansions())
       << "," << python_representation(vtf.errors())
       << "," << python_representation(vtf.properties())
       << ")";
    return os;
}

List<MultiIndex> keys(const ValidatedTaylorModelDP& tm) {
    List<MultiIndex> r(tm.argument_size());
    for(ValidatedTaylorModelDP::ConstIterator iter=tm.begin(); iter!=tm.end(); ++iter) {
        r.append(iter->index());
    }
    return r;
}

ValidatedScalarFunction unrestrict(const ValidatedScalarFunctionModelDP& fm) {
    return ValidatedScalarFunction(fm.raw_pointer()->_clone());
}

ValidatedVectorFunction unrestrict(const ValidatedVectorFunctionModelDP& fm) {
    return ValidatedVectorFunction(fm.raw_pointer()->_clone());
}


ExactIntervalType _range1(const ValidatedTaylorModelDP&);
ExactIntervalType _range2(const ValidatedTaylorModelDP&);
ExactIntervalType _range3(const ValidatedTaylorModelDP&);

} // namespace Ariadne

static constexpr auto self = pybind11::detail::self;

Sweeper<FloatDP> make_threshold_sweeper(DoublePrecision pr, double x) { return new ThresholdSweeper<FloatDP>(pr,x); }
Sweeper<FloatDP> make_graded_sweeper(DoublePrecision pr, SizeType n) { return new GradedSweeper<FloatDP>(pr,n); }

Void export_expansion(pybind11::module& module)
{
//    from_python< Expansion<MultiIndex,FloatDPApproximation> >();
//    from_python< Expansion<MultiIndex,FloatDPBounds> >();
//    from_python< Vector< Expansion<MultiIndex,FloatDPApproximation> > >();

    pybind11::class_< ExpansionValue<MultiIndex,FloatDPApproximation> > expansion_value_class(module,"ExpansionValue");
    expansion_value_class.def(pybind11::init<MultiIndex,FloatDPApproximation>());
    // TODO: Add get/set for coefficient
    // TODO: Use property for index
    //expansion_value_class.add_property("index", (MultiIndex const&(ExpansionValue<MultiIndex,FloatDPApproximation>::*)()const)&ExpansionValue<MultiIndex,FloatDPApproximation>::index);
    expansion_value_class.def("index", (const MultiIndex&(ExpansionValue<MultiIndex,FloatDPApproximation>::*)()const)&ExpansionValue<MultiIndex,FloatDPApproximation>::index);
    expansion_value_class.def("__str__", &__cstr__<ExpansionValue<MultiIndex,FloatDPApproximation>>);

}


Void export_sweeper(pybind11::module& module)
{
    pybind11::class_<Sweeper<FloatDP>> sweeper_class(module,"Sweeper");
    sweeper_class.def(pybind11::init<Sweeper<FloatDP>>());
    module.def("ThresholdSweeper", &make_threshold_sweeper );
    module.def("GradedSweeper", &make_graded_sweeper );
    sweeper_class.def("__str__", &__cstr__<Sweeper<FloatDP>>);
}

/*
Expansion<MultiIndex,FloatDPValue>const& get_expansion(ValidatedTaylorModelDP const& tm) { return tm.expansion(); }

template<class F> Void export_validated_taylor_model(pybind11::module& module)
{
    typedef SizeType SizeType;
    typedef TaylorModel<ValidatedTag,F> ValidatedTaylorModelType;

    pybind11::class_<ValidatedTaylorModelDP> taylor_model_class(module,"ValidatedTaylorModelDP", pybind11::init<ValidatedTaylorModelDP>());
    taylor_model_class.def( pybind11::init< SizeType,SweeperDP >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ValidatedTaylorModelDP&))&keys);
    taylor_model_class.def("value", (const FloatDPValue&(ValidatedTaylorModelDP::*)()const) &ValidatedTaylorModelDP::value);
    taylor_model_class.def("gradient", (const FloatDPValue&(ValidatedTaylorModelDP::*)(SizeType)const) &ValidatedTaylorModelDP::gradient_value);
    taylor_model_class.def("error", (const FloatDPError&(ValidatedTaylorModelDP::*)()const) &ValidatedTaylorModelDP::error);
    taylor_model_class.def("expansion", (const Expansion<MultiIndex,FloatDPValue>&(*)(ValidatedTaylorModelDP const&)) &get_expansion);
    taylor_model_class.def("set_error", (Void(ValidatedTaylorModelDP::*)(const FloatDPError&)) &ValidatedTaylorModelDP::set_error);
    taylor_model_class.def("argument_size", &ValidatedTaylorModelDP::argument_size);
    taylor_model_class.def("domain", &ValidatedTaylorModelDP::domain);
    taylor_model_class.def("range", &ValidatedTaylorModelDP::range);
    taylor_model_class.def("set_sweeper", &ValidatedTaylorModelDP::set_sweeper);
    taylor_model_class.def("sweeper", &ValidatedTaylorModelDP::sweeper);
    taylor_model_class.def("sweep", (ValidatedTaylorModelDP&(ValidatedTaylorModelDP::*)()) &ValidatedTaylorModelDP::sweep, pybind11::return_value_policy::reference_existing_object);
    taylor_model_class.def("__getitem__", &__getitem__<ValidatedTaylorModelDP,MultiIndex,FloatDPValue>);
    taylor_model_class.def("__setitem__",&__setitem__<ValidatedTaylorModelDP,MultiIndex,FloatDPValue>);
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

    taylor_model_class.def("constant",(ValidatedTaylorModelDP(*)(SizeType, const ValidatedNumericType&,SweeperDP))&ValidatedTaylorModelDP::constant);
    taylor_model_class.def("coordinate",(ValidatedTaylorModelDP(*)(SizeType, SizeType,SweeperDP))&ValidatedTaylorModelDP::coordinate);

    taylor_model_class.staticmethod("constant");
    taylor_model_class.staticmethod("coordinate");

    //def("max",(ValidatedTaylorModelDP(*)(const ValidatedTaylorModelDP&,const ValidatedTaylorModelDP&))&max);
    //def("min",(ValidatedTaylorModelDP(*)(const ValidatedTaylorModelDP&,const ValidatedTaylorModelDP&))&min);
    //def("abs",(ValidatedTaylorModelDP(*)(const ValidatedTaylorModelDP&))&abs);

    typedef AlgebraOperations<ValidatedTaylorModel> Operations;
    module.def("pos",&Operations::_pos);
    module.def("neg",&Operations::_neg);
    module.def("rec",&Operations::_rec);
///    module.def("pow",&Operations::_pow);
    module.def("sqrt",&Operations::_sqrt);
    module.def("exp",&Operations::_exp);
    module.def("log",&Operations::_log);
    module.def("sin",&Operations::_sin);
    module.def("cos",&Operations::_cos);
    module.def("tan",&Operations::_tan);
    module.def("atan",&Operations::_atan);

    taylor_model_class.def("range", (UpperIntervalType(ValidatedTaylorModelDP::*)()const) &ValidatedTaylorModelDP::range);

    //def("evaluate", (ValidatedNumericType(*)(const ValidatedTaylorModelDP&, const Vector<ValidatedNumericType>&))&evaluate);
    //def("split",(ValidatedTaylorModelDP(*)(const ValidatedTaylorModelDP&,SizeType,SplitPart)) &split);

    from_python< Vector<ValidatedTaylorModelDP> >();
    to_python< Vector<ValidatedTaylorModelDP> >();

}

Void export_approximate_taylor_model(pybind11::module& module)
{
    typedef SizeType SizeType;
    typedef ApproximateTaylorModelDP ApproximateTaylorModelDP;

    pybind11::class_<ApproximateTaylorModelDP> taylor_model_class(module,"ApproximateTaylorModelDP", pybind11::init<ApproximateTaylorModelDP>());
    taylor_model_class.def( pybind11::init< SizeType,SweeperDP >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ApproximateTaylorModelDP&))&keys);
    taylor_model_class.def("value", (const FloatDPApproximation&(ApproximateTaylorModelDP::*)()const) &ApproximateTaylorModelDP::value);
    taylor_model_class.def("gradient", (const FloatDPApproximation&(ApproximateTaylorModelDP::*)(SizeType)const) &ApproximateTaylorModelDP::gradient_value);
    taylor_model_class.def("expansion", (const Expansion<MultiIndex,FloatDPApproximation>&(*)(ApproximateTaylorModelDP const&)) &get_expansion);
    taylor_model_class.def("argument_size", &ApproximateTaylorModelDP::argument_size);
    taylor_model_class.def("domain", &ApproximateTaylorModelDP::domain);
    taylor_model_class.def("range", &ApproximateTaylorModelDP::range);
    taylor_model_class.def("set_sweeper", &ApproximateTaylorModelDP::set_sweeper);
    taylor_model_class.def("sweeper", &ApproximateTaylorModelDP::sweeper);
    taylor_model_class.def("sweep", (ApproximateTaylorModelDP&(ApproximateTaylorModelDP::*)()) &ApproximateTaylorModelDP::sweep, pybind11::return_value_policy,reference_existing_object);
    taylor_model_class.def("__getitem__", &__getitem__<ApproximateTaylorModelDP,MultiIndex,FloatDPApproximation>);
    taylor_model_class.def("__setitem__",&__setitem__<ApproximateTaylorModelDP,MultiIndex,FloatDPApproximation>);
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

    taylor_model_class.def("constant",(ApproximateTaylorModelDP(*)(SizeType, const ApproximateNumericType&,SweeperDP))&ApproximateTaylorModelDP::constant);
    taylor_model_class.def("coordinate",(ApproximateTaylorModelDP(*)(SizeType, SizeType,SweeperDP))&ApproximateTaylorModelDP::coordinate);

    taylor_model_class.staticmethod("constant");
    taylor_model_class.staticmethod("coordinate");

    //def("max",(ApproximateTaylorModelDP(*)(const ApproximateTaylorModelDP&,const ApproximateTaylorModelDP&))&max);
    //def("min",(ApproximateTaylorModelDP(*)(const ApproximateTaylorModelDP&,const ApproximateTaylorModelDP&))&min);
    //def("abs",(ApproximateTaylorModelDP(*)(const ApproximateTaylorModelDP&))&abs);

    typedef AlgebraOperations<ApproximateTaylorModelDP> Operations;
    module.def("pos",&Operations::_pos);
    module.def("neg",&Operations::_neg);
    module.def("rec",&Operations::_rec);
///    module.def("pow",&Operations::_pow);
    module.def("sqrt",&Operations::_sqrt);
    module.def("exp",&Operations::_exp);
    module.def("log",&Operations::_log);
    module.def("sin",&Operations::_sin);
    module.def("cos",&Operations::_cos);
    module.def("tan",&Operations::_tan);
    module.def("atan",&Operations::_atan);

    //def("evaluate", (ApproximateNumericType(*)(const ApproximateTaylorModelDP&, const Vector<ApproximateNumericType>&))&evaluate);
    //def("split",(ApproximateTaylorModelDP(*)(const ApproximateTaylorModelDP&,SizeType,SplitPart)) &split);

    from_python< Vector<ApproximateTaylorModelDP> >();
    to_python< Vector<ApproximateTaylorModelDP> >();
}

*/



Void export_scalar_function_model(pybind11::module& module)
{
    pybind11::class_<ValidatedScalarFunctionModelDP> scalar_function_model_class(module,"ValidatedScalarFunctionModel");
    scalar_function_model_class.def(pybind11::init<ValidatedScalarFunctionModelDP>());
    scalar_function_model_class.def(pybind11::init<ValidatedScalarTaylorFunctionModelDP>());
    scalar_function_model_class.def("argument_size", &ValidatedScalarFunctionModelDP::argument_size);
    scalar_function_model_class.def("domain", &ValidatedScalarFunctionModelDP::domain);
    scalar_function_model_class.def("codomain", &ValidatedScalarFunctionModelDP::codomain);
    scalar_function_model_class.def("range", &ValidatedScalarFunctionModelDP::range);
    scalar_function_model_class.def("clobber", &ValidatedScalarFunctionModelDP::clobber);
    scalar_function_model_class.def("error", &ValidatedScalarFunctionModelDP::error);
    scalar_function_model_class.def("__call__", (FloatDPBounds(ValidatedScalarFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedScalarFunctionModelDP::operator());
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
    scalar_function_model_class.def("__str__", &__cstr__<ValidatedScalarFunctionModelDP>);
    scalar_function_model_class.def("__repr__", &__crepr__<ValidatedScalarFunctionModelDP>);
    //scalar_function_model_class.def("__repr__",&__repr__<ValidatedScalarFunctionModelDP>);

    module.def("evaluate", (ValidatedNumericType(*)(const ValidatedScalarFunctionModelDP&,const Vector<ValidatedNumericType>&)) &evaluate);
//    module.def("partial_evaluate", (ValidatedScalarFunctionModelDP(*)(const ValidatedScalarFunctionModelDP&,SizeType,const ValidatedNumericType&)) &partial_evaluate);

    module.def("compose", (ValidatedScalarFunctionModelDP(*)(const ValidatedScalarFunctionModelDP&, const ValidatedVectorFunctionModelDP&)) &compose);
    module.def("compose", (ValidatedScalarFunctionModelDP(*)(const ValidatedScalarFunction&, const ValidatedVectorFunctionModelDP&)) &compose);

    module.def("unrestrict", (ValidatedScalarFunction(*)(const ValidatedScalarFunctionModelDP&)) &unrestrict);

    module.def("antiderivative", (ValidatedScalarFunctionModelDP(*)(ValidatedScalarFunctionModelDP,SizeType,ValidatedNumericType)) &_antiderivative_<ValidatedScalarFunctionModelDP,SizeType,ValidatedNumericType>);

}

Void export_vector_function_model(pybind11::module& module)
{
    pybind11::class_<ValidatedVectorFunctionModelDP> vector_function_model_class(module,"ValidatedVectorFunctionModel");
    vector_function_model_class.def(pybind11::init<ValidatedVectorFunctionModelDP>());
    vector_function_model_class.def(pybind11::init<ValidatedVectorTaylorFunctionModelDP>());
    vector_function_model_class.def("result_size", &ValidatedVectorFunctionModelDP::result_size);
    vector_function_model_class.def("argument_size", &ValidatedVectorFunctionModelDP::argument_size);
    vector_function_model_class.def("domain", &ValidatedVectorFunctionModelDP::domain);
    vector_function_model_class.def("codomain", &ValidatedVectorFunctionModelDP::codomain);
    vector_function_model_class.def("range", &ValidatedVectorFunctionModelDP::range);
    //vector_function_model_class.def("__getslice__", (ValidatedVectorTaylorFunctionModelDP(*)(const ValidatedVectorTaylorFunctionModelDP&,Int,Int))&__getslice__);
    vector_function_model_class.def("__getitem__", &__getitem__<ValidatedVectorFunctionModelDP,SizeType,ValidatedScalarFunctionModelDP>);
    vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorFunctionModelDP,SizeType,ValidatedScalarFunctionModelDP>);
    //vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorFunctionModelDP,SizeType,ValidatedScalarFunction>);
    vector_function_model_class.def("__call__", (Vector<FloatDPBounds>(ValidatedVectorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedVectorFunctionModelDP::operator());
    vector_function_model_class.def(-self);
    vector_function_model_class.def(self+self);
    vector_function_model_class.def(self-self);
    vector_function_model_class.def(ValidatedNumericType()*self);
    vector_function_model_class.def(self*ValidatedNumericType());
    vector_function_model_class.def("__str__", &__cstr__<ValidatedVectorFunctionModelDP>);
    vector_function_model_class.def("__repr__", &__crepr__<ValidatedVectorFunctionModelDP>);
    //export_vector_function_model.def(pybind11::module& module, "__repr__",&__repr__<ValidatedVectorFunctionModelDP>);


//    module.def("evaluate", (Vector<ValidatedNumericType>(*)(const ValidatedVectorFunctionModelDP&,const Vector<ValidatedNumericType>&)) &evaluate);

    module.def("compose", (ValidatedVectorFunctionModelDP(*)(const ValidatedVectorFunctionModelDP&,const ValidatedVectorFunctionModelDP&)) &compose);
    module.def("compose", (ValidatedVectorFunctionModelDP(*)(const ValidatedVectorFunction&,const ValidatedVectorFunctionModelDP&)) &compose);

    module.def("unrestrict", (ValidatedVectorFunction(*)(const ValidatedVectorFunctionModelDP&)) &unrestrict);

    module.def("join", (ValidatedVectorFunctionModelDP(*)(const ValidatedScalarFunctionModelDP&,const ValidatedScalarFunctionModelDP&)) &join);
    module.def("join", (ValidatedVectorFunctionModelDP(*)(const ValidatedScalarFunctionModelDP&,const ValidatedVectorFunctionModelDP&)) &join);
    module.def("join", (ValidatedVectorFunctionModelDP(*)(const ValidatedVectorFunctionModelDP&,const ValidatedScalarFunctionModelDP&)) &join);
    module.def("join", (ValidatedVectorFunctionModelDP(*)(const ValidatedVectorFunctionModelDP&,const ValidatedVectorFunctionModelDP&)) &join);

    module.def("antiderivative", (ValidatedVectorFunctionModelDP(*)(ValidatedVectorFunctionModelDP,SizeType,ValidatedNumericType)) &_antiderivative_<ValidatedVectorFunctionModelDP,SizeType,ValidatedNumericType>);
    module.def("antiderivative", (ValidatedVectorFunctionModelDP(*)(ValidatedVectorFunctionModelDP,SizeType,ValidatedNumber)) &_antiderivative_<ValidatedVectorFunctionModelDP,SizeType,ValidatedNumber>);

//    to_python< List<ValidatedVectorFunctionModelDP> >();
}



Void export_scalar_taylor_function(pybind11::module& module)
{
    typedef ValidatedScalarTaylorFunctionModelDP F;
    typedef ValidatedVectorTaylorFunctionModelDP VF;
    typedef typename F::DomainType D;
    typedef typename F::NumericType X;
    typedef Vector<X> VX;
    typedef SizeType I;
    typedef typename X::GenericType Y;
    typedef Vector<Y> VY;

    typedef F const& Fcr;
    
    pybind11::class_<ValidatedScalarTaylorFunctionModelDP> scalar_taylor_function_class(module,"ValidatedScalarTaylorFunctionModel");
    scalar_taylor_function_class.def(pybind11::init<ValidatedScalarTaylorFunctionModelDP>());
    scalar_taylor_function_class.def(pybind11::init<ExactBoxType,ValidatedTaylorModelDP>());
    scalar_taylor_function_class.def(pybind11::init< ExactBoxType,SweeperDP >());
    scalar_taylor_function_class.def(pybind11::init< ExactBoxType, const EffectiveScalarFunction&,SweeperDP >());
    scalar_taylor_function_class.def(pybind11::init< ExactBoxType, Expansion<MultiIndex,FloatDPValue>, FloatDPError, SweeperDP >());
    scalar_taylor_function_class.def("error", (const FloatDPError&(ValidatedScalarTaylorFunctionModelDP::*)()const) &ValidatedScalarTaylorFunctionModelDP::error);
    scalar_taylor_function_class.def("set_error", (Void(ValidatedScalarTaylorFunctionModelDP::*)(const FloatDPError&)) &ValidatedScalarTaylorFunctionModelDP::set_error);
    scalar_taylor_function_class.def("argument_size", &F::argument_size);
    scalar_taylor_function_class.def("domain", &F::domain);
    scalar_taylor_function_class.def("codomain", &F::codomain);
    scalar_taylor_function_class.def("range", &F::range);
    scalar_taylor_function_class.def("model", (const ValidatedTaylorModelDP&(ValidatedScalarTaylorFunctionModelDP::*)()const)&ValidatedScalarTaylorFunctionModelDP::model);
    scalar_taylor_function_class.def("polynomial", (Polynomial<ValidatedNumericType>(ValidatedScalarTaylorFunctionModelDP::*)()const)&ValidatedScalarTaylorFunctionModelDP::polynomial);
    scalar_taylor_function_class.def("number_of_nonzeros", (SizeType(ValidatedScalarTaylorFunctionModelDP::*)()const)&ValidatedScalarTaylorFunctionModelDP::number_of_nonzeros);
//    scalar_taylor_function_class.def("set_sweeper", &ValidatedScalarTaylorFunctionModelDP::set_sweeper);
//    scalar_taylor_function_class.def("sweeper", &ValidatedScalarTaylorFunctionModelDP::sweeper);
//    scalar_taylor_function_class.def("sweep", (ValidatedScalarTaylorFunctionModelDP&(ValidatedScalarTaylorFunctionModelDP::*)()) &ValidatedScalarTaylorFunctionModelDP::sweep, pybind11::return_value_policy,reference_existing_object);
    scalar_taylor_function_class.def("__getitem__", &__getitem__<ValidatedScalarTaylorFunctionModelDP,MultiIndex,FloatDPValue>);
    scalar_taylor_function_class.def("__setitem__",&__setitem__<ValidatedScalarTaylorFunctionModelDP,MultiIndex,FloatDPValue>);
    scalar_taylor_function_class.def(+self);
    scalar_taylor_function_class.def(-self);
    scalar_taylor_function_class.def(self+self);
    scalar_taylor_function_class.def(self-self);
    scalar_taylor_function_class.def(self*self);
    scalar_taylor_function_class.def(self/self);
    scalar_taylor_function_class.def(self+X());
    scalar_taylor_function_class.def(self-X());
    scalar_taylor_function_class.def(self*X());
    scalar_taylor_function_class.def(self/X());
    scalar_taylor_function_class.def(X()+self);
    scalar_taylor_function_class.def(X()-self);
    scalar_taylor_function_class.def(X()*self);
    scalar_taylor_function_class.def(X()/self);
    scalar_taylor_function_class.def(self+=X());
    scalar_taylor_function_class.def(self-=X());
    scalar_taylor_function_class.def(self*=X());
    scalar_taylor_function_class.def(self/=X());
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

    //scalar_taylor_function_class.def("__str__",(StringType(*)(const ValidatedScalarTaylorFunctionModelDP&)) &__str__);
    //scalar_taylor_function_class.def("_cstr_",(StringType(*)(const ValidatedScalarTaylorFunctionModelDP&)) &_cstr_);
    //scalar_taylor_function_class.def("__repr__",(StringType(*)(const ValidatedScalarTaylorFunctionModelDP&)) &__repr__);
    scalar_taylor_function_class.def("value", (const FloatDPValue&(ValidatedScalarTaylorFunctionModelDP::*)()const) &ValidatedScalarTaylorFunctionModelDP::value);
//    scalar_taylor_function_class.def("sweep", (ValidatedScalarTaylorFunctionModelDP&(ValidatedScalarTaylorFunctionModelDP::*)())&ValidatedScalarTaylorFunctionModelDP::sweep,pybind11::return_value_policy::reference_internal);
    scalar_taylor_function_class.def("clobber", (Void(ValidatedScalarTaylorFunctionModelDP::*)()) &ValidatedScalarTaylorFunctionModelDP::clobber);
    scalar_taylor_function_class.def("set_properties",&ValidatedScalarTaylorFunctionModelDP::set_properties);
    scalar_taylor_function_class.def("properties",&ValidatedScalarTaylorFunctionModelDP::properties);
    scalar_taylor_function_class.def("__call__", (FloatDPApproximation(ValidatedScalarTaylorFunctionModelDP::*)(const Vector<FloatDPApproximation>&)const) &ValidatedScalarTaylorFunctionModelDP::operator());
    scalar_taylor_function_class.def("__call__", (FloatDPBounds(ValidatedScalarTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedScalarTaylorFunctionModelDP::operator());
    scalar_taylor_function_class.def("gradient", (Covector<FloatDPBounds>(ValidatedScalarTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedScalarTaylorFunctionModelDP::gradient);
    scalar_taylor_function_class.def("function", (ValidatedScalarFunction(ValidatedScalarTaylorFunctionModelDP::*)()const) &ValidatedScalarTaylorFunctionModelDP::function);
    scalar_taylor_function_class.def("polynomial", (Polynomial<FloatDPBounds>(ValidatedScalarTaylorFunctionModelDP::*)()const) &ValidatedScalarTaylorFunctionModelDP::polynomial);
    scalar_taylor_function_class.def("restriction", (F(*)(F,D)) &_restriction_<F,D>);
//    scalar_taylor_function_class.def("extension",&_extension_<F,D>);

    scalar_taylor_function_class.def_static("zero",(ValidatedScalarTaylorFunctionModelDP(*)(const ExactBoxType&,SweeperDP))&ValidatedScalarTaylorFunctionModelDP::zero);
    scalar_taylor_function_class.def_static("constant",(ValidatedScalarTaylorFunctionModelDP(*)(const ExactBoxType&,const ValidatedNumericType&,SweeperDP))&ValidatedScalarTaylorFunctionModelDP::constant);
    scalar_taylor_function_class.def_static("coordinate",(ValidatedScalarTaylorFunctionModelDP(*)(const ExactBoxType&,SizeType,SweeperDP))&ValidatedScalarTaylorFunctionModelDP::coordinate);


    module.def("restriction", (F(*)(F,D)) &_restriction_<F,D>);
    module.def("join",(VF(*)(F,F)) &_join_<F,F>);
//    module.def("extension",&_extension_<F,D>);
    module.def("embed",(F(*)(D,F,D)) &_embed_<D,F,D>);
    module.def("split",(Pair<F,F>(*)(F,I)) &_split_<F,I>);
    module.def("evaluate",(X(*)(F,VX)) &_evaluate_<F,VX>);
    module.def("evaluate",(F(*)(F,I,X)) &_partial_evaluate_<F,I,X>);
    module.def("midpoint",(F(*)(F)) &_midpoint_<F>);
    module.def("derivative", (F(*)(F,I)) &_derivative_<F,I>);
    module.def("antiderivative", (F(*)(F,I)) &_antiderivative_<F,I>);
    module.def("antiderivative", (F(*)(F,I,X)) &_antiderivative_<F,I,X>);

    module.def("inconsistent", (Bool(*)(F,F)) &_inconsistent_<F,F>);
    module.def("refines", (Bool(*)(F,F)) &_refines_<F,F>);
    module.def("refinement", (F(*)(F,F)) &_refinement_<F,F>);

    module.def("max", (F(*)(Fcr,Fcr)) &_max_<F,F>); module.def("min", (F(*)(Fcr,Fcr)) &_min_<F,F>); module.def("abs", (F(*)(Fcr)) &_abs_<F>);

    module.def("neg",(F(*)(Fcr))&_neg_<F>); module.def("rec",(F(*)(Fcr))&_rec_<F>); module.def("sqr",(F(*)(Fcr))&_sqr_<F>); 
    module.def("pow",(F(*)(Fcr,Int const&))&_pow_<F,Int>);
    module.def("sqrt",(F(*)(Fcr))&_sqrt_<F>); module.def("exp",(F(*)(Fcr))&_exp_<F>); module.def("log",(F(*)(Fcr))&_log_<F>); module.def("atan",(F(*)(Fcr))&_atan_<F>);
    module.def("sin",(F(*)(Fcr))&_sin_<F>); module.def("cos",(F(*)(Fcr))&_cos_<F>); module.def("tan",(F(*)(Fcr))&_tan_<F>);

//    to_python< Vector<ValidatedScalarTaylorFunctionModelDP> >();
}

Void export_vector_taylor_function(pybind11::module& module)
{

    typedef SizeType I;
    typedef ValidatedScalarFunction SFN;
    typedef ValidatedVectorFunction VFN;
    typedef ValidatedScalarTaylorFunctionModelDP SF;
    typedef ValidatedVectorTaylorFunctionModelDP VF;
    typedef typename VF::DomainType D;
    typedef typename D::ScalarType Di;
    typedef typename VF::NumericType X;
    typedef Vector<X> VX;

    pybind11::class_<ValidatedVectorTaylorFunctionModelDP> vector_taylor_function_class(module,"ValidatedVectorTaylorFunctionModel");
    vector_taylor_function_class.def( pybind11::init<ValidatedVectorTaylorFunctionModelDP>());
    vector_taylor_function_class.def( pybind11::init< SizeType, ExactBoxType, SweeperDP >());
    vector_taylor_function_class.def( pybind11::init< ExactBoxType,const EffectiveVectorFunction&,SweeperDP >());
    vector_taylor_function_class.def(pybind11::init< ExactBoxType, Vector< Expansion<MultiIndex,FloatDPValue> >, Vector<FloatDPError>, SweeperDP >());
    vector_taylor_function_class.def( pybind11::init< Vector<ValidatedScalarTaylorFunctionModelDP> >());
    vector_taylor_function_class.def("_len_", &ValidatedVectorTaylorFunctionModelDP::result_size);
    vector_taylor_function_class.def("result_size", &ValidatedVectorTaylorFunctionModelDP::result_size);
    vector_taylor_function_class.def("argument_size", &ValidatedVectorTaylorFunctionModelDP::argument_size);
    vector_taylor_function_class.def("domain", &ValidatedVectorTaylorFunctionModelDP::domain);
    vector_taylor_function_class.def("codomain", &ValidatedVectorTaylorFunctionModelDP::codomain);
    // FIXME: Omitted since const and non-const versions
    // vector_taylor_function_class.def("models", &ValidatedVectorTaylorFunctionModelDP::models);
    vector_taylor_function_class.def("centre", &ValidatedVectorTaylorFunctionModelDP::centre);
    vector_taylor_function_class.def("range", &ValidatedVectorTaylorFunctionModelDP::range);
    vector_taylor_function_class.def("errors", &ValidatedVectorTaylorFunctionModelDP::errors);
//    vector_taylor_function_class.def("sweep", (ValidatedVectorTaylorFunctionModelDP&(ValidatedVectorTaylorFunctionModelDP::*)())&ValidatedVectorTaylorFunctionModelDP::sweep,pybind11::return_value_policy::reference_existing_object);
    vector_taylor_function_class.def("clobber", (Void(ValidatedVectorTaylorFunctionModelDP::*)()) &ValidatedVectorTaylorFunctionModelDP::clobber);
    vector_taylor_function_class.def("set_properties",&ValidatedVectorTaylorFunctionModelDP::set_properties);
    vector_taylor_function_class.def("properties",&ValidatedVectorTaylorFunctionModelDP::properties);
//    vector_taylor_function_class.def("sweep", (ValidatedVectorTaylorFunctionModelDP&(ValidatedVectorTaylorFunctionModelDP::*)()) &ValidatedVectorTaylorFunctionModelDP::sweep, pybind11::return_value_policy::reference_existing_object);
    //vector_taylor_function_class.def("__getslice__", &__getslice__<ValidatedVectorTaylorFunctionModelDP,SizeType,SizeType,ValidatedScalarTaylorFunctionModelDP>);
    vector_taylor_function_class.def("__getslice__", (ValidatedVectorTaylorFunctionModelDP(*)(const ValidatedVectorTaylorFunctionModelDP&,Int,Int))&__getslice__);
    vector_taylor_function_class.def("__getitem__", &__getitem__<ValidatedVectorTaylorFunctionModelDP,SizeType,ValidatedScalarTaylorFunctionModelDP>);
    vector_taylor_function_class.def("__setitem__",&__setitem__<ValidatedVectorTaylorFunctionModelDP,SizeType,ValidatedScalarTaylorFunctionModelDP>);
    vector_taylor_function_class.def(-self);
    vector_taylor_function_class.def(self+self);
    vector_taylor_function_class.def(self-self);
    vector_taylor_function_class.def(self+Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self-Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self*ValidatedNumericType());
    vector_taylor_function_class.def(self/ValidatedNumericType());
//    vector_taylor_function_class.def(self*ValidatedScalarTaylorFunctionModelDP());
    vector_taylor_function_class.def(self+=Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self-=Vector<ValidatedNumericType>());
    vector_taylor_function_class.def(self*=ValidatedNumericType());
    vector_taylor_function_class.def(self/=ValidatedNumericType());
    vector_taylor_function_class.def(self+=self);
    vector_taylor_function_class.def(self-=self);
    vector_taylor_function_class.def("__str__", &__cstr__<ValidatedVectorTaylorFunctionModelDP>);
    vector_taylor_function_class.def("__repr__", &__crepr__<ValidatedVectorTaylorFunctionModelDP>);
    vector_taylor_function_class.def("clobber", (Void(ValidatedVectorTaylorFunctionModelDP::*)()) &ValidatedVectorTaylorFunctionModelDP::clobber);
    vector_taylor_function_class.def("__call__", (Vector<FloatDPApproximation>(ValidatedVectorTaylorFunctionModelDP::*)(const Vector<FloatDPApproximation>&)const) &ValidatedVectorTaylorFunctionModelDP::operator());
    vector_taylor_function_class.def("__call__", (Vector<FloatDPBounds>(ValidatedVectorTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedVectorTaylorFunctionModelDP::operator());
     //vector_taylor_function_class.def("jacobian", (Vector<FloatDPBounds>(ValidatedVectorTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedVectorTaylorFunctionModelDP::jacobian);
    vector_taylor_function_class.def("polynomials", (Vector< Polynomial<FloatDPBounds> >(ValidatedVectorTaylorFunctionModelDP::*)()const) &ValidatedVectorTaylorFunctionModelDP::polynomials);
    vector_taylor_function_class.def("function", (ValidatedVectorFunction(ValidatedVectorTaylorFunctionModelDP::*)()const) &ValidatedVectorTaylorFunctionModelDP::function);


    vector_taylor_function_class.def_static("constant",(ValidatedVectorTaylorFunctionModelDP(*)(const ExactBoxType&, const Vector<ValidatedNumericType>&,SweeperDP))&ValidatedVectorTaylorFunctionModelDP::constant);
    vector_taylor_function_class.def_static("identity",(ValidatedVectorTaylorFunctionModelDP(*)(const ExactBoxType&,SweeperDP))&ValidatedVectorTaylorFunctionModelDP::identity);

    module.def("inconsistent", (Bool(*)(VF,VF)) &_inconsistent_<VF,VF>);
    module.def("refinement", (VF(*)(VF,VF)) &_refinement_<VF,VF>);
    module.def("refines", (Bool(*)(VF,VF)) &_refines_<VF,VF>);

    module.def("join", (VF(*)(VF,VF)) &_join_<VF,VF>); module.def("join", (VF(*)(VF,SF)) &_join_<VF,SF>); 
        module.def("join", (VF(*)(SF,VF)) &_join_<SF,VF>); // module.def("join", (VF(*)(SF,SF)) &_join_<SF,SF>);
    module.def("combine", (VF(*)(VF,VF)) &_combine_<VF,VF>); module.def("combine", (VF(*)(VF,SF)) &_combine_<VF,SF>); 
        module.def("combine", (VF(*)(SF,VF)) &_combine_<SF,VF>); module.def("combine", (VF(*)(SF,SF)) &_combine_<SF,SF>);
    module.def("embed", (VF(*)(VF,Di)) &_embed_<VF,Di>); module.def("embed", (VF(*)(VF,D)) &_embed_<VF,D>); 
        module.def("embed", (VF(*)(D,VF)) &_embed_<D,VF>); module.def("embed", (VF(*)(D,VF,D)) &_embed_<D,VF,D>);

    module.def("restriction", (VF(*)(VF,D)) &_restriction_<VF,D>); module.def("restriction", (VF(*)(VF,I,Di)) &_restriction_<VF,I,Di>);
//    module.def("split", &_split_<VF,I>);

    module.def("evaluate", (VX(*)(VF,VX)) &_evaluate_<VF,VX>);
    module.def("partial_evaluate", (VF(*)(VF,I,X)) &_partial_evaluate_<VF,I,X>);
    module.def("compose", (VF(*)(VF,VF)) &_compose_<VF,VF>);
    module.def("compose", (SF(*)(SF,VF)) &_compose_<SF,VF>);
    module.def("compose", (SF(*)(SFN,VF)) &_compose_<SFN,VF>);
    module.def("compose", (VF(*)(VFN,VF)) &_compose_<VFN,VF>);
    module.def("unchecked_compose", (SF(*)(SF,VF)) &_compose_<SF,VF>);
    module.def("unchecked_compose", (VF(*)(VF,VF)) &_compose_<VF,VF>);
    module.def("derivative", (VF(*)(VF,I)) &_derivative_<VF,I>);
    module.def("antiderivative", (VF(*)(VF,I)) &_antiderivative_<VF,I>);
    module.def("antiderivative", (VF(*)(VF,I,X)) &_antiderivative_<VF,I,X>);

//    module.def("evaluate",(Vector<FloatDPApproximation>(ValidatedVectorTaylorFunctionModelDP::*)(const Vector<FloatDPApproximation>&)const) &ValidatedVectorTaylorFunctionModelDP::evaluate);

//    from_python<ValidatedVectorTaylorFunctionModelDP>();


}

Void calculus_submodule(pybind11::module& module)
{
    export_expansion(module);
    export_sweeper(module);
//    export_approximate_taylor_model(module);
//    export_validated_taylor_model(module);
    export_scalar_function_model(module);
    export_vector_function_model(module);
    export_scalar_taylor_function(module);
    export_vector_taylor_function(module);
}


