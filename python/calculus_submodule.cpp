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
    if(start<0) { start+=tf.result_size(); }
    if(stop<0) { stop+=tf.result_size(); }
    ARIADNE_ASSERT_MSG(0<=start&&start<=stop&&SizeType(stop)<=tf.result_size(),
            "result_size="<<tf.result_size()<<", start="<<start<<", stop="<<stop);
    return ValidatedVectorTaylorFunctionModelDP(tf.domain(),Vector<ValidatedTaylorModelDP>(project(tf.models(),range(static_cast<SizeType>(start),static_cast<SizeType>(stop)))));
}


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
    from_python() { converter::registry::push_back(&convertible,&construct,type_id< Vector<X> >()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* coefficient) {
        list lst=boost::python::extract<list>(obj_ptr);
        Void* storage = ((converter::rvalue_from_python_storage< Vector<X> >*) coefficient)->storage.bytes;
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
        list lst=boost::python::extract<boost::python::list>(obj_ptr);
        ValidatedVectorTaylorFunctionModelDP* tf_ptr = new (storage) ValidatedVectorTaylorFunctionModelDP(static_cast<SizeType>(len(lst)));
        for(SizeType i=0; i!=tf_ptr->result_size(); ++i) { tf_ptr->set(i,boost::python::extract<ValidatedScalarTaylorFunctionModelDP>(lst[i])); }
        coefficient->convertible = storage;
    }
};


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

Sweeper<FloatDP> make_threshold_sweeper(DoublePrecision pr, double x) { return new ThresholdSweeper<FloatDP>(pr,x); }
Sweeper<FloatDP> make_graded_sweeper(DoublePrecision pr, SizeType n) { return new GradedSweeper<FloatDP>(pr,n); }

Void export_expansion()
{
    from_python< Expansion<MultiIndex,FloatDPApproximation> >();
    from_python< Expansion<MultiIndex,FloatDPBounds> >();
    from_python< Vector< Expansion<MultiIndex,FloatDPApproximation> > >();

    class_< ExpansionValue<MultiIndex,FloatDPApproximation> > expansion_value_class("ExpansionValue", init<MultiIndex,FloatDPApproximation>());
    // TODO: Add get/set for coefficient
    // TODO: Use property for index
    //expansion_value_class.add_property("index", (MultiIndex const&(ExpansionValue<MultiIndex,FloatDPApproximation>::*)()const)&ExpansionValue<MultiIndex,FloatDPApproximation>::index);
    expansion_value_class.def("index", (const MultiIndex&(ExpansionValue<MultiIndex,FloatDPApproximation>::*)()const)&ExpansionValue<MultiIndex,FloatDPApproximation>::index, return_value_policy<copy_const_reference>());
    expansion_value_class.def(self_ns::str(self));

}


Void export_sweeper()
{
    class_<Sweeper<FloatDP>> sweeper_class("Sweeper", init<Sweeper<FloatDP>>());
    def("ThresholdSweeper", &make_threshold_sweeper );
    def("GradedSweeper", &make_graded_sweeper );
    sweeper_class.def(self_ns::str(self));
}

/*
Expansion<MultiIndex,FloatDPValue>const& get_expansion(ValidatedTaylorModelDP const& tm) { return tm.expansion(); }

template<class F> Void export_validated_taylor_model()
{
    typedef SizeType SizeType;
    typedef TaylorModel<ValidatedTag,F> ValidatedTaylorModelType;

    class_<ValidatedTaylorModelDP> taylor_model_class("ValidatedTaylorModelDP", init<ValidatedTaylorModelDP>());
    taylor_model_class.def( init< SizeType,SweeperDP >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ValidatedTaylorModelDP&))&keys);
    taylor_model_class.def("value", (const FloatDPValue&(ValidatedTaylorModelDP::*)()const) &ValidatedTaylorModelDP::value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("gradient", (const FloatDPValue&(ValidatedTaylorModelDP::*)(SizeType)const) &ValidatedTaylorModelDP::gradient_value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("error", (const FloatDPError&(ValidatedTaylorModelDP::*)()const) &ValidatedTaylorModelDP::error, return_value_policy<copy_const_reference>());
    taylor_model_class.def("expansion", (const Expansion<MultiIndex,FloatDPValue>&(*)(ValidatedTaylorModelDP const&)) &get_expansion, return_value_policy<copy_const_reference>());
    taylor_model_class.def("set_error", (Void(ValidatedTaylorModelDP::*)(const FloatDPError&)) &ValidatedTaylorModelDP::set_error);
    taylor_model_class.def("argument_size", &ValidatedTaylorModelDP::argument_size);
    taylor_model_class.def("domain", &ValidatedTaylorModelDP::domain);
    taylor_model_class.def("range", &ValidatedTaylorModelDP::range);
    taylor_model_class.def("set_sweeper", &ValidatedTaylorModelDP::set_sweeper);
    taylor_model_class.def("sweeper", &ValidatedTaylorModelDP::sweeper);
    taylor_model_class.def("sweep", (ValidatedTaylorModelDP&(ValidatedTaylorModelDP::*)()) &ValidatedTaylorModelDP::sweep, return_value_policy<reference_existing_object>());
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

    taylor_model_class.def("range", (UpperIntervalType(ValidatedTaylorModelDP::*)()const) &ValidatedTaylorModelDP::range);

    //def("evaluate", (ValidatedNumericType(*)(const ValidatedTaylorModelDP&, const Vector<ValidatedNumericType>&))&evaluate);
    //def("split",(ValidatedTaylorModelDP(*)(const ValidatedTaylorModelDP&,SizeType,SplitPart)) &split);

    from_python< Vector<ValidatedTaylorModelDP> >();
    to_python< Vector<ValidatedTaylorModelDP> >();

}

Void export_approximate_taylor_model()
{
    typedef SizeType SizeType;
    typedef ApproximateTaylorModelDP ApproximateTaylorModelDP;

    class_<ApproximateTaylorModelDP> taylor_model_class("ApproximateTaylorModelDP", init<ApproximateTaylorModelDP>());
    taylor_model_class.def( init< SizeType,SweeperDP >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ApproximateTaylorModelDP&))&keys);
    taylor_model_class.def("value", (const FloatDPApproximation&(ApproximateTaylorModelDP::*)()const) &ApproximateTaylorModelDP::value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("gradient", (const FloatDPApproximation&(ApproximateTaylorModelDP::*)(SizeType)const) &ApproximateTaylorModelDP::gradient_value, return_value_policy<copy_const_reference>());
    taylor_model_class.def("expansion", (const Expansion<MultiIndex,FloatDPApproximation>&(*)(ApproximateTaylorModelDP const&)) &get_expansion, return_value_policy<copy_const_reference>());
    taylor_model_class.def("argument_size", &ApproximateTaylorModelDP::argument_size);
    taylor_model_class.def("domain", &ApproximateTaylorModelDP::domain);
    taylor_model_class.def("range", &ApproximateTaylorModelDP::range);
    taylor_model_class.def("set_sweeper", &ApproximateTaylorModelDP::set_sweeper);
    taylor_model_class.def("sweeper", &ApproximateTaylorModelDP::sweeper);
    taylor_model_class.def("sweep", (ApproximateTaylorModelDP&(ApproximateTaylorModelDP::*)()) &ApproximateTaylorModelDP::sweep, return_value_policy<reference_existing_object>());
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

    //def("evaluate", (ApproximateNumericType(*)(const ApproximateTaylorModelDP&, const Vector<ApproximateNumericType>&))&evaluate);
    //def("split",(ApproximateTaylorModelDP(*)(const ApproximateTaylorModelDP&,SizeType,SplitPart)) &split);

    from_python< Vector<ApproximateTaylorModelDP> >();
    to_python< Vector<ApproximateTaylorModelDP> >();
}

*/



Void export_scalar_function_model()
{
    class_<ValidatedScalarFunctionModelDP> scalar_function_model_class("ValidatedScalarFunctionModel",init<ValidatedScalarFunctionModelDP>());
    scalar_function_model_class.def(init<ValidatedScalarTaylorFunctionModelDP>());
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

    def("evaluate", (ValidatedNumericType(*)(const ValidatedScalarFunctionModelDP&,const Vector<ValidatedNumericType>&)) &evaluate);
//    def("partial_evaluate", (ValidatedScalarFunctionModelDP(*)(const ValidatedScalarFunctionModelDP&,SizeType,const ValidatedNumericType&)) &partial_evaluate);

    def("compose", (ValidatedScalarFunctionModelDP(*)(const ValidatedScalarFunctionModelDP&, const ValidatedVectorFunctionModelDP&)) &compose);
    def("compose", (ValidatedScalarFunctionModelDP(*)(const ValidatedScalarFunction&, const ValidatedVectorFunctionModelDP&)) &compose);

    def("unrestrict", (ValidatedScalarFunction(*)(const ValidatedScalarFunctionModelDP&)) &unrestrict);

    def("antiderivative", &_antiderivative_<ValidatedScalarFunctionModelDP,SizeType,ValidatedNumericType>);

}

Void export_vector_function_model()
{
    class_<ValidatedVectorFunctionModelDP> vector_function_model_class("ValidatedVectorFunctionModel",init<ValidatedVectorFunctionModelDP>());
    vector_function_model_class.def(init<ValidatedVectorTaylorFunctionModelDP>());
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
    //export_vector_function_model.def("__repr__",&__repr__<ValidatedVectorFunctionModelDP>);


//    def("evaluate", (Vector<ValidatedNumericType>(*)(const ValidatedVectorFunctionModelDP&,const Vector<ValidatedNumericType>&)) &evaluate);

    def("compose", (ValidatedVectorFunctionModelDP(*)(const ValidatedVectorFunctionModelDP&,const ValidatedVectorFunctionModelDP&)) &compose);
    def("compose", (ValidatedVectorFunctionModelDP(*)(const ValidatedVectorFunction&,const ValidatedVectorFunctionModelDP&)) &compose);

    def("unrestrict", (ValidatedVectorFunction(*)(const ValidatedVectorFunctionModelDP&)) &unrestrict);

    def("join", (ValidatedVectorFunctionModelDP(*)(const ValidatedScalarFunctionModelDP&,const ValidatedScalarFunctionModelDP&)) &join);
    def("join", (ValidatedVectorFunctionModelDP(*)(const ValidatedScalarFunctionModelDP&,const ValidatedVectorFunctionModelDP&)) &join);
    def("join", (ValidatedVectorFunctionModelDP(*)(const ValidatedVectorFunctionModelDP&,const ValidatedScalarFunctionModelDP&)) &join);
    def("join", (ValidatedVectorFunctionModelDP(*)(const ValidatedVectorFunctionModelDP&,const ValidatedVectorFunctionModelDP&)) &join);

    def("antiderivative", &_antiderivative_<ValidatedVectorFunctionModelDP,SizeType,ValidatedNumericType>);
    def("antiderivative", &_antiderivative_<ValidatedVectorFunctionModelDP,SizeType,ValidatedNumber>);

    to_python< List<ValidatedVectorFunctionModelDP> >();
}



Void export_scalar_taylor_function()
{
    typedef ValidatedScalarTaylorFunctionModelDP F;
    typedef ValidatedVectorTaylorFunctionModelDP VF;
    typedef typename F::DomainType D;
    typedef typename F::NumericType X;
    typedef Vector<X> VX;
    typedef SizeType I;
    typedef typename X::GenericType Y;
    typedef Vector<Y> VY;

    class_<ValidatedScalarTaylorFunctionModelDP> scalar_taylor_function_class("ValidatedScalarTaylorFunctionModel",init<ValidatedScalarTaylorFunctionModelDP>());
    scalar_taylor_function_class.def(init<ExactBoxType,ValidatedTaylorModelDP>());
    scalar_taylor_function_class.def(init< ExactBoxType,SweeperDP >());
    scalar_taylor_function_class.def(init< ExactBoxType, const EffectiveScalarFunction&,SweeperDP >());
    scalar_taylor_function_class.def(init< ExactBoxType, Expansion<MultiIndex,FloatDPValue>, FloatDPError, SweeperDP >());
    scalar_taylor_function_class.def("error", (const FloatDPError&(ValidatedScalarTaylorFunctionModelDP::*)()const) &ValidatedScalarTaylorFunctionModelDP::error, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("set_error", (Void(ValidatedScalarTaylorFunctionModelDP::*)(const FloatDPError&)) &ValidatedScalarTaylorFunctionModelDP::set_error);
    scalar_taylor_function_class.def("argument_size", &F::argument_size);
    scalar_taylor_function_class.def("domain", &F::domain);
    scalar_taylor_function_class.def("codomain", &F::codomain);
    scalar_taylor_function_class.def("range", &F::range);
    scalar_taylor_function_class.def("model", (const ValidatedTaylorModelDP&(ValidatedScalarTaylorFunctionModelDP::*)()const)&ValidatedScalarTaylorFunctionModelDP::model, return_value_policy<copy_const_reference>());
    scalar_taylor_function_class.def("polynomial", (Polynomial<ExactIntervalType>(ValidatedScalarTaylorFunctionModelDP::*)()const)&ValidatedScalarTaylorFunctionModelDP::polynomial);
    scalar_taylor_function_class.def("number_of_nonzeros", (SizeType(ValidatedScalarTaylorFunctionModelDP::*)()const)&ValidatedScalarTaylorFunctionModelDP::number_of_nonzeros);
//    scalar_taylor_function_class.def("set_sweeper", &ValidatedScalarTaylorFunctionModelDP::set_sweeper);
//    scalar_taylor_function_class.def("sweeper", &ValidatedScalarTaylorFunctionModelDP::sweeper);
//    scalar_taylor_function_class.def("sweep", (ValidatedScalarTaylorFunctionModelDP&(ValidatedScalarTaylorFunctionModelDP::*)()) &ValidatedScalarTaylorFunctionModelDP::sweep, return_value_policy<reference_existing_object>());
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
    scalar_taylor_function_class.def("value", (const FloatDPValue&(ValidatedScalarTaylorFunctionModelDP::*)()const) &ValidatedScalarTaylorFunctionModelDP::value,return_value_policy<copy_const_reference>());
//    scalar_taylor_function_class.def("sweep", (ValidatedScalarTaylorFunctionModelDP&(ValidatedScalarTaylorFunctionModelDP::*)())&ValidatedScalarTaylorFunctionModelDP::sweep,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("clobber", (ValidatedScalarTaylorFunctionModelDP&(ValidatedScalarTaylorFunctionModelDP::*)()) &ValidatedScalarTaylorFunctionModelDP::clobber,return_value_policy<reference_existing_object>());
    scalar_taylor_function_class.def("set_properties",&ValidatedScalarTaylorFunctionModelDP::set_properties);
    scalar_taylor_function_class.def("properties",&ValidatedScalarTaylorFunctionModelDP::properties);
    scalar_taylor_function_class.def("__call__", (FloatDPApproximation(ValidatedScalarTaylorFunctionModelDP::*)(const Vector<FloatDPApproximation>&)const) &ValidatedScalarTaylorFunctionModelDP::operator());
    scalar_taylor_function_class.def("__call__", (FloatDPBounds(ValidatedScalarTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedScalarTaylorFunctionModelDP::operator());
    scalar_taylor_function_class.def("gradient", (Covector<FloatDPBounds>(ValidatedScalarTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedScalarTaylorFunctionModelDP::gradient);
    scalar_taylor_function_class.def("function", (EffectiveScalarFunction(ValidatedScalarTaylorFunctionModelDP::*)()const) &ValidatedScalarTaylorFunctionModelDP::function);
    scalar_taylor_function_class.def("polynomial", (Polynomial<FloatDPBounds>(ValidatedScalarTaylorFunctionModelDP::*)()const) &ValidatedScalarTaylorFunctionModelDP::polynomial);
    scalar_taylor_function_class.def("restriction",&_restriction_<F,D>);
//    scalar_taylor_function_class.def("extension",&_extension_<F,D>);

    scalar_taylor_function_class.def("zero",(ValidatedScalarTaylorFunctionModelDP(*)(const ExactBoxType&,SweeperDP))&ValidatedScalarTaylorFunctionModelDP::zero);
    scalar_taylor_function_class.def("constant",(ValidatedScalarTaylorFunctionModelDP(*)(const ExactBoxType&,const ValidatedNumericType&,SweeperDP))&ValidatedScalarTaylorFunctionModelDP::constant);
    scalar_taylor_function_class.def("coordinate",(ValidatedScalarTaylorFunctionModelDP(*)(const ExactBoxType&,SizeType,SweeperDP))&ValidatedScalarTaylorFunctionModelDP::coordinate);


    scalar_taylor_function_class.staticmethod("constant");
    scalar_taylor_function_class.staticmethod("coordinate");

    def("restriction",&_restriction_<F,D>);
    def("join",&_join_<F,F>);
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

    to_python< Vector<ValidatedScalarTaylorFunctionModelDP> >();
}

Void export_vector_taylor_function()
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

    class_<ValidatedVectorTaylorFunctionModelDP> vector_taylor_function_class("ValidatedVectorTaylorFunctionModel", init<ValidatedVectorTaylorFunctionModelDP>());
    vector_taylor_function_class.def( init< SizeType, ExactBoxType, SweeperDP >());
    vector_taylor_function_class.def( init< ExactBoxType,const EffectiveVectorFunction&,SweeperDP >());
    vector_taylor_function_class.def(init< ExactBoxType, Vector< Expansion<MultiIndex,FloatDPValue> >, Vector<FloatDPError>, SweeperDP >());
    vector_taylor_function_class.def( init< Vector<ValidatedScalarTaylorFunctionModelDP> >());
    vector_taylor_function_class.def("_len_", &ValidatedVectorTaylorFunctionModelDP::result_size);
    vector_taylor_function_class.def("result_size", &ValidatedVectorTaylorFunctionModelDP::result_size);
    vector_taylor_function_class.def("argument_size", &ValidatedVectorTaylorFunctionModelDP::argument_size);
    vector_taylor_function_class.def("domain", &ValidatedVectorTaylorFunctionModelDP::domain);
    vector_taylor_function_class.def("codomain", &ValidatedVectorTaylorFunctionModelDP::codomain);
    // FIXME: Omitted since const and non-const versions
    // vector_taylor_function_class.def("models", &ValidatedVectorTaylorFunctionModelDP::models, return_value_policy<copy_const_reference>());
    vector_taylor_function_class.def("centre", &ValidatedVectorTaylorFunctionModelDP::centre);
    vector_taylor_function_class.def("range", &ValidatedVectorTaylorFunctionModelDP::range);
    vector_taylor_function_class.def("errors", &ValidatedVectorTaylorFunctionModelDP::errors);
//    vector_taylor_function_class.def("sweep", (ValidatedVectorTaylorFunctionModelDP&(ValidatedVectorTaylorFunctionModelDP::*)())&ValidatedVectorTaylorFunctionModelDP::sweep,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("clobber", (ValidatedVectorTaylorFunctionModelDP&(ValidatedVectorTaylorFunctionModelDP::*)()) &ValidatedVectorTaylorFunctionModelDP::clobber,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("set_properties",&ValidatedVectorTaylorFunctionModelDP::set_properties);
    vector_taylor_function_class.def("properties",&ValidatedVectorTaylorFunctionModelDP::properties);
//    vector_taylor_function_class.def("sweep", (ValidatedVectorTaylorFunctionModelDP&(ValidatedVectorTaylorFunctionModelDP::*)()) &ValidatedVectorTaylorFunctionModelDP::sweep, return_value_policy<reference_existing_object>());
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
    vector_taylor_function_class.def("clobber", (ValidatedVectorTaylorFunctionModelDP&(ValidatedVectorTaylorFunctionModelDP::*)()) &ValidatedVectorTaylorFunctionModelDP::clobber,return_value_policy<reference_existing_object>());
    vector_taylor_function_class.def("__call__", (Vector<FloatDPApproximation>(ValidatedVectorTaylorFunctionModelDP::*)(const Vector<FloatDPApproximation>&)const) &ValidatedVectorTaylorFunctionModelDP::operator());
    vector_taylor_function_class.def("__call__", (Vector<FloatDPBounds>(ValidatedVectorTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedVectorTaylorFunctionModelDP::operator());
     //vector_taylor_function_class.def("jacobian", (Vector<FloatDPBounds>(ValidatedVectorTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedVectorTaylorFunctionModelDP::jacobian);
    vector_taylor_function_class.def("polynomials", (Vector< Polynomial<FloatDPBounds> >(ValidatedVectorTaylorFunctionModelDP::*)()const) &ValidatedVectorTaylorFunctionModelDP::polynomials);
    vector_taylor_function_class.def("function", (EffectiveVectorFunction(ValidatedVectorTaylorFunctionModelDP::*)()const) &ValidatedVectorTaylorFunctionModelDP::function);


    vector_taylor_function_class.def("constant",(ValidatedVectorTaylorFunctionModelDP(*)(const ExactBoxType&, const Vector<ValidatedNumericType>&,SweeperDP))&ValidatedVectorTaylorFunctionModelDP::constant);
    vector_taylor_function_class.def("identity",(ValidatedVectorTaylorFunctionModelDP(*)(const ExactBoxType&,SweeperDP))&ValidatedVectorTaylorFunctionModelDP::identity);

    vector_taylor_function_class.staticmethod("constant");
    vector_taylor_function_class.staticmethod("identity");

    def("inconsistent", &_inconsistent_<VF,VF>);
    def("refinement", &_refinement_<VF,VF>);
    def("refines", &_refines_<VF,VF>);

    def("join", &_join_<VF,VF>); def("join", &_join_<VF,SF>); def("join", &_join_<SF,VF>); // def("join", &_join_<SF,SF>);
    def("combine", &_combine_<VF,VF>); def("combine", &_combine_<VF,SF>); def("combine", &_combine_<SF,VF>); def("combine", &_combine_<SF,SF>);
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

//    def("evaluate",(Vector<FloatDPApproximation>(ValidatedVectorTaylorFunctionModelDP::*)(const Vector<FloatDPApproximation>&)const) &ValidatedVectorTaylorFunctionModelDP::evaluate);

    from_python<ValidatedVectorTaylorFunctionModelDP>();


}

Void calculus_submodule()
{
    export_expansion();
    export_sweeper();
//    export_approximate_taylor_model();
//    export_validated_taylor_model();
    export_scalar_function_model();
    export_vector_function_model();
    export_scalar_taylor_function();
    export_vector_taylor_function();
}


