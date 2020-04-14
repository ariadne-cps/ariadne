/***************************************************************************
 *            calculus_submodule.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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



//ValidatedNumericType evaluate(const ValidatedScalarMultivariateFunctionModelDP& f, const Vector<ValidatedNumericType>& x) { return f(x); }
//Vector<ValidatedNumericType> evaluate(const ValidatedVectorMultivariateFunctionModelDP& f, const Vector<ValidatedNumericType>& x) { return f(x); }

//ValidatedScalarMultivariateFunctionModelDP partial_evaluate(const ValidatedScalarMultivariateFunctionModelDP&, SizeType, const ValidatedNumericType&);
//ValidatedVectorMultivariateFunctionModelDP partial_evaluate(const ValidatedVectorMultivariateFunctionModelDP&, SizeType, const ValidatedNumericType&);

//ValidatedScalarMultivariateFunctionModelDP antiderivative(const ValidatedScalarMultivariateFunctionModelDP&, SizeType, ValidatedNumericType);
//ValidatedVectorMultivariateFunctionModelDP antiderivative(const ValidatedVectorMultivariateFunctionModelDP&, SizeType, ValidatedNumericType);

//ValidatedScalarMultivariateFunctionModelDP compose(const ValidatedScalarMultivariateFunctionModelDP&, const ValidatedVectorMultivariateFunctionModelDP&);
//ValidatedScalarMultivariateFunctionModelDP compose(const ValidatedScalarMultivariateFunction&, const ValidatedVectorMultivariateFunctionModelDP&);
//ValidatedVectorMultivariateFunctionModelDP compose(const ValidatedVectorMultivariateFunctionModelDP&, const ValidatedVectorMultivariateFunctionModelDP&);
//ValidatedVectorMultivariateFunctionModelDP compose(const ValidatedVectorMultivariateFunction&, const ValidatedVectorMultivariateFunctionModelDP&);

//ValidatedVectorMultivariateFunctionModelDP join(const ValidatedScalarMultivariateFunctionModelDP&, const ValidatedScalarMultivariateFunctionModelDP&);
//ValidatedVectorMultivariateFunctionModelDP join(const ValidatedScalarMultivariateFunctionModelDP&, const ValidatedVectorMultivariateFunctionModelDP&);
//ValidatedVectorMultivariateFunctionModelDP join(const ValidatedVectorMultivariateFunctionModelDP&, const ValidatedScalarMultivariateFunctionModelDP&);
//ValidatedVectorMultivariateFunctionModelDP join(const ValidatedVectorMultivariateFunctionModelDP&, const ValidatedVectorMultivariateFunctionModelDP&);

template<class P, class D, class C, class PR, class PRE> OutputStream& operator<<(OutputStream& os, const Representation< FunctionModel<P,D,C,PR,PRE> >& frepr) {
    static_cast<const FunctionInterface<P,D,C>&>(frepr.reference()).repr(os); return os;
}

ValidatedVectorMultivariateTaylorFunctionModelDP __getslice__(const ValidatedVectorMultivariateTaylorFunctionModelDP& tf, Int start, Int stop) {
    Int rs = tf.result_size();
    if(start<0) { start+=rs; }
    if(stop<0) { stop+=rs; }
    ARIADNE_ASSERT_MSG(0<=start&&start<=stop&&stop<=rs,
            "result_size="<<rs<<", start="<<start<<", stop="<<stop);
    return ValidatedVectorMultivariateTaylorFunctionModelDP(tf.domain(),Vector<ValidatedTaylorModelDP>(project(tf.models(),range(static_cast<SizeType>(start),static_cast<SizeType>(stop)))));
}



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
    return os << PythonRepresentation< Vector<ExactIntervalType> >(cast_vector(bx.reference())); }

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

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ValidatedScalarMultivariateTaylorFunctionModelDP>& repr) {
    const ValidatedScalarMultivariateTaylorFunctionModelDP& stf=repr.reference();
    os << std::setprecision(17);
    os << "ValidatedScalarMultivariateTaylorFunctionModelDP"
       << "(" << python_representation(stf.domain())
       << "," << python_representation(stf.expansion())
       << "," << python_representation(stf.error())
       << "," << python_representation(stf.properties())
       << ")";
    return os;
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ValidatedVectorMultivariateTaylorFunctionModelDP>& repr) {
    const ValidatedVectorMultivariateTaylorFunctionModelDP& vtf=repr.reference();
    os << std::setprecision(17);
    os << "ValidatedVectorMultivariateTaylorFunctionModelDP"
       << "(" << python_representation(vtf.domain())
       << "," << python_representation(vtf.expansions())
       << "," << python_representation(vtf.errors())
       << "," << python_representation(vtf.properties())
       << ")";
    return os;
}

template<class P, class F> List<MultiIndex> keys(const TaylorModel<P,F>& tm) {
    List<MultiIndex> r(tm.argument_size());
    for(auto iter=tm.expansion().begin(); iter!=tm.expansion().end(); ++iter) {
        r.append(iter->index());
    }
    return r;
}

ValidatedScalarMultivariateFunction unrestrict(const ValidatedScalarMultivariateFunctionModelDP& fm) {
    return ValidatedScalarMultivariateFunction(fm.raw_pointer()->_clone());
}

ValidatedVectorMultivariateFunction unrestrict(const ValidatedVectorMultivariateFunctionModelDP& fm) {
    return ValidatedVectorMultivariateFunction(fm.raw_pointer()->_clone());
}

} // namespace Ariadne

static constexpr auto self = pybind11::detail::self;

Sweeper<FloatDP> make_threshold_sweeper(DoublePrecision pr, double x) { return new ThresholdSweeper<FloatDP>(pr,x); }
Sweeper<FloatDP> make_graded_sweeper(DoublePrecision pr, SizeType n) { return new GradedSweeper<FloatDP>(pr,n); }



Void export_expansion(pybind11::module& module)
{
//    from_python< Expansion<MultiIndex,FloatDPApproximation> >();
//    from_python< Expansion<MultiIndex,FloatDPBounds> >();
//    from_python< Vector< Expansion<MultiIndex,FloatDPApproximation> > >();
//    from_python< Vector< Expansion<MultiIndex,FloatDPBounds> > >();
}


Void export_sweeper(pybind11::module& module)
{
    pybind11::class_<Sweeper<FloatDP>> sweeper_class(module,"Sweeper");
    sweeper_class.def(pybind11::init<Sweeper<FloatDP>>());
    sweeper_class.def("__str__", &__cstr__<Sweeper<FloatDP>>);

    pybind11::class_<ThresholdSweeper<FloatDP>> threshold_sweeper_class(module,"ThresholdSweeper");
    threshold_sweeper_class.def(pybind11::init<DP,double>());
    threshold_sweeper_class.def("__str__", &__cstr__<ThresholdSweeper<FloatDP>>);
    sweeper_class.def(pybind11::init<ThresholdSweeper<FloatDP>>());
    pybind11::implicitly_convertible<ThresholdSweeper<FloatDP>,Sweeper<FloatDP>>();

    pybind11::class_<GradedSweeper<FloatDP>> graded_sweeper_class(module,"GradedSweeper");
    graded_sweeper_class.def(pybind11::init<DP,int>());
    graded_sweeper_class.def("__str__", &__cstr__<GradedSweeper<FloatDP>>);
    sweeper_class.def(pybind11::init<GradedSweeper<FloatDP>>());
    pybind11::implicitly_convertible<GradedSweeper<FloatDP>,Sweeper<FloatDP>>();
}


Expansion<MultiIndex,FloatDPValue>const& get_expansion(ValidatedTaylorModelDP const& tm) { return tm.expansion(); }
Expansion<MultiIndex,FloatDPApproximation>const& get_expansion(ApproximateTaylorModelDP const& tm) { return tm.expansion(); }


Void export_validated_taylor_model(pybind11::module& module)
{
    typedef ValidatedTaylorModelDP ModelType;
    typedef NumericType<ModelType> NumericType;
    typedef GenericType<NumericType> GenericNumericType;

    Tag<GenericNumericType> generic_number_tag;

    pybind11::class_<ValidatedTaylorModelDP> taylor_model_class(module,"ValidatedTaylorModelDP");
    taylor_model_class.def(pybind11::init<ValidatedTaylorModelDP>());
    taylor_model_class.def(pybind11::init< SizeType,SweeperDP >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ValidatedTaylorModelDP&))&keys);
    taylor_model_class.def("value", (const FloatDPValue(ValidatedTaylorModelDP::*)()const) &ValidatedTaylorModelDP::value);
    taylor_model_class.def("error", (const FloatDPError&(ValidatedTaylorModelDP::*)()const) &ValidatedTaylorModelDP::error);
    taylor_model_class.def("expansion", (const Expansion<MultiIndex,FloatDPValue>&(*)(ValidatedTaylorModelDP const&)) &get_expansion);
    taylor_model_class.def("set_error", (Void(ValidatedTaylorModelDP::*)(const FloatDPError&)) &ValidatedTaylorModelDP::set_error);
    taylor_model_class.def("argument_size", &ValidatedTaylorModelDP::argument_size);
    taylor_model_class.def("domain", &ValidatedTaylorModelDP::domain);
    taylor_model_class.def("range", &ValidatedTaylorModelDP::range);
    taylor_model_class.def("set_sweeper", &ValidatedTaylorModelDP::set_sweeper);
    taylor_model_class.def("sweeper", &ValidatedTaylorModelDP::sweeper);
    taylor_model_class.def("sweep", (ValidatedTaylorModelDP&(ValidatedTaylorModelDP::*)()) &ValidatedTaylorModelDP::sweep, pybind11::return_value_policy::reference);
    taylor_model_class.def("__getitem__", &__getitem__<ValidatedTaylorModelDP,MultiIndex,FloatDPValue>);
    taylor_model_class.def("__setitem__",&__setitem__<ValidatedTaylorModelDP,MultiIndex,FloatDPValue>);

    define_elementary_algebra(module,taylor_model_class);
    define_inplace_algebra(module,taylor_model_class);
    define_mixed_arithmetic(module, taylor_model_class, generic_number_tag);

    define_lattice(module,taylor_model_class);


    taylor_model_class.def("__str__", &__cstr__<ValidatedTaylorModelDP>);

    taylor_model_class.def_static("constant",(ValidatedTaylorModelDP(*)(SizeType, const FloatDPBounds&,SweeperDP))&ValidatedTaylorModelDP::constant);
    taylor_model_class.def_static("coordinate",(ValidatedTaylorModelDP(*)(SizeType, SizeType,SweeperDP))&ValidatedTaylorModelDP::coordinate);

    taylor_model_class.def("range", (UpperIntervalType(ValidatedTaylorModelDP::*)()const) &ValidatedTaylorModelDP::range);

    module.def("evaluate",&_evaluate_<ValidatedTaylorModelDP,Vector<FloatDPBounds>>);
    //module.def("split",&_split_<ValidatedTaylorModelDP,SizeType,SplitPart>);

    //to_python< Vector<ValidatedTaylorModelDP> >();
    //from_python< Vector<ValidatedTaylorModelDP> >();

}

Void export_approximate_taylor_model(pybind11::module& module)
{
    typedef ApproximateTaylorModelDP ModelType;

    pybind11::class_<ModelType> taylor_model_class(module,"ApproximateTaylorModelDP");
    taylor_model_class.def(pybind11::init<ModelType>());
    taylor_model_class.def(pybind11::init< SizeType,SweeperDP >());
    taylor_model_class.def("keys", (List<MultiIndex>(*)(const ModelType&))&keys);
    taylor_model_class.def("value", (const FloatDPApproximation(ModelType::*)()const) &ModelType::value);
    taylor_model_class.def("expansion", (const Expansion<MultiIndex,FloatDPApproximation>&(*)(ModelType const&)) &get_expansion);
    taylor_model_class.def("argument_size", &ModelType::argument_size);
    taylor_model_class.def("domain", &ModelType::domain);
    taylor_model_class.def("range", &ModelType::range);
    taylor_model_class.def("set_sweeper", &ModelType::set_sweeper);
    taylor_model_class.def("sweeper", &ModelType::sweeper);
    taylor_model_class.def("sweep", (ModelType&(ModelType::*)()) &ModelType::sweep, pybind11::return_value_policy::reference);
    taylor_model_class.def("__getitem__", &__getitem__<ModelType,MultiIndex,FloatDPApproximation>);
    taylor_model_class.def("__setitem__",&__setitem__<ModelType,MultiIndex,FloatDPApproximation>);

    define_elementary_algebra(module,taylor_model_class);
    define_inplace_algebra(module,taylor_model_class);
    define_lattice(module,taylor_model_class);

    taylor_model_class.def("__str__",&__cstr__<ModelType>);

    taylor_model_class.def_static("constant",(ModelType(*)(SizeType, const FloatDPApproximation&,SweeperDP))&ModelType::constant);
    taylor_model_class.def_static("coordinate",(ModelType(*)(SizeType, SizeType,SweeperDP))&ModelType::coordinate);

    //FIXME: Not in C++ API
    //module.def("evaluate",&_evaluate_<ModelType,Vector<ApproximateNumericType>>);
    //module.def("split",&_split_<ModelType,SizeType,SplitPart>);

    //from_python< Vector<ModelType> >();
    //to_python< Vector<ModelType> >();
}



Void export_scalar_function_model(pybind11::module& module)
{
    typedef ValidatedScalarMultivariateFunctionModelDP::NumericType NumericType;

    pybind11::class_<ValidatedScalarMultivariateFunctionModelDP> scalar_function_model_class(module,"ValidatedScalarMultivariateFunctionModelDP");
    scalar_function_model_class.def(pybind11::init<ValidatedScalarMultivariateFunctionModelDP>());
    scalar_function_model_class.def(pybind11::init<ValidatedScalarMultivariateTaylorFunctionModelDP>());
    scalar_function_model_class.def("argument_size", &ValidatedScalarMultivariateFunctionModelDP::argument_size);
    scalar_function_model_class.def("domain", &ValidatedScalarMultivariateFunctionModelDP::domain);
    scalar_function_model_class.def("codomain", &ValidatedScalarMultivariateFunctionModelDP::codomain);
    scalar_function_model_class.def("range", &ValidatedScalarMultivariateFunctionModelDP::range);
    scalar_function_model_class.def("clobber", &ValidatedScalarMultivariateFunctionModelDP::clobber);
    scalar_function_model_class.def("error", &ValidatedScalarMultivariateFunctionModelDP::error);
    scalar_function_model_class.def("__call__", (FloatDPBounds(ValidatedScalarMultivariateFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedScalarMultivariateFunctionModelDP::operator());
    scalar_function_model_class.def(self+self);
    scalar_function_model_class.def(self-self);
    scalar_function_model_class.def(self*self);
    scalar_function_model_class.def(self/self);
    scalar_function_model_class.def(self+NumericType());
    scalar_function_model_class.def(self-NumericType());
    scalar_function_model_class.def(self*NumericType());
    scalar_function_model_class.def(self/NumericType());
    scalar_function_model_class.def(NumericType()+self);
    scalar_function_model_class.def(NumericType()-self);
    scalar_function_model_class.def(NumericType()*self);
    scalar_function_model_class.def(NumericType()/self);
    scalar_function_model_class.def("__str__", &__cstr__<ValidatedScalarMultivariateFunctionModelDP>);
    scalar_function_model_class.def("__repr__", &__crepr__<ValidatedScalarMultivariateFunctionModelDP>);

    module.def("evaluate", &_evaluate_<ValidatedScalarMultivariateFunctionModelDP,Vector<NumericType>>);
    module.def("partial_evaluate", &_partial_evaluate_<ValidatedScalarMultivariateFunctionModelDP,SizeType,NumericType>);

    module.def("compose", _compose_<ValidatedScalarMultivariateFunctionModelDP,ValidatedVectorMultivariateFunctionModelDP>);
    module.def("compose", _compose_<ValidatedScalarMultivariateFunction,ValidatedVectorMultivariateFunctionModelDP>);

    module.def("unrestrict", (ValidatedScalarMultivariateFunction(*)(const ValidatedScalarMultivariateFunctionModelDP&)) &unrestrict);

    module.def("antiderivative",  &_antiderivative_<ValidatedScalarMultivariateFunctionModelDP,SizeType,NumericType>);

}

Void export_vector_function_model(pybind11::module& module)
{
    typedef ValidatedScalarMultivariateFunctionModelDP::NumericType NumericType;
    //using VectorMultivariateFunctionModelType = ValidatedVectorMultivariateFunctionModelDP;
    //using ScalarMultivariateFunctionModelType = ValidatedScalarMultivariateFunctionModelDP;

    pybind11::class_<ValidatedVectorMultivariateFunctionModelDP> vector_function_model_class(module,"ValidatedVectorMultivariateFunctionModelDP");
    vector_function_model_class.def(pybind11::init<ValidatedVectorMultivariateFunctionModelDP>());
//    vector_function_model_class.def(pybind11::init([](Array<ValidatedScalarMultivariateFunctionModelDP> ary){return ValidatedVectorMultivariateFunctionModelDP(Vector<ValidatedScalarMultivariateFunctionModelDP>(ary));}));
    vector_function_model_class.def(pybind11::init<ValidatedVectorMultivariateTaylorFunctionModelDP>());
    vector_function_model_class.def("result_size", &ValidatedVectorMultivariateFunctionModelDP::result_size);
    vector_function_model_class.def("argument_size", &ValidatedVectorMultivariateFunctionModelDP::argument_size);
    vector_function_model_class.def("domain", &ValidatedVectorMultivariateFunctionModelDP::domain);
    vector_function_model_class.def("codomain", &ValidatedVectorMultivariateFunctionModelDP::codomain);
    vector_function_model_class.def("range", &ValidatedVectorMultivariateFunctionModelDP::range);
    //vector_function_model_class.def("__getslice__", (ValidatedVectorMultivariateTaylorFunctionModelDP(*)(const ValidatedVectorMultivariateTaylorFunctionModelDP&,Int,Int))&__getslice__);
    vector_function_model_class.def("__getitem__", &__getitem__<ValidatedVectorMultivariateFunctionModelDP,SizeType,ValidatedScalarMultivariateFunctionModelDP>);
    vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorMultivariateFunctionModelDP,SizeType,ValidatedScalarMultivariateFunctionModelDP>);
    //vector_function_model_class.def("__setitem__",&__setitem__<ValidatedVectorMultivariateFunctionModelDP,SizeType,ValidatedScalarMultivariateFunction>);
    vector_function_model_class.def("__call__", (Vector<FloatDPBounds>(ValidatedVectorMultivariateFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedVectorMultivariateFunctionModelDP::operator());

    // NOTE: Not all operations are exported in C++ API.
    //define_vector_algebra_arithmetic<VectorMultivariateFunctionModelType,ScalarMultivariateFunctionModelType,NumericType>(module,vector_function_model_class);
    //define_vector_arithmetic<VectorMultivariateFunctionModelType,ScalarMultivariateFunctionModelType>(module,vector_function_model_class);

    vector_function_model_class.def("__str__", &__cstr__<ValidatedVectorMultivariateFunctionModelDP>);
    vector_function_model_class.def("__repr__", &__crepr__<ValidatedVectorMultivariateFunctionModelDP>);
    //export_vector_function_model.def(pybind11::module& module, "__repr__",&__repr__<ValidatedVectorMultivariateFunctionModelDP>);


//    module.def("evaluate", (Vector<ValidatedNumericType>(*)(const ValidatedVectorMultivariateFunctionModelDP&,const Vector<ValidatedNumericType>&)) &evaluate);
    module.def("partial_evaluate", &_partial_evaluate_<ValidatedVectorMultivariateFunctionModelDP,SizeType,NumericType>);

    module.def("compose", &_compose_<ValidatedVectorMultivariateFunctionModelDP,ValidatedVectorMultivariateFunctionModelDP>);
    module.def("compose", &_compose_<ValidatedVectorMultivariateFunction,ValidatedVectorMultivariateFunctionModelDP>);

    module.def("unrestrict", (ValidatedVectorMultivariateFunction(*)(const ValidatedVectorMultivariateFunctionModelDP&)) &unrestrict);

    module.def("join", &_join_<ValidatedScalarMultivariateFunctionModelDP,ValidatedScalarMultivariateFunctionModelDP>);
    module.def("join", &_join_<ValidatedScalarMultivariateFunctionModelDP,ValidatedVectorMultivariateFunctionModelDP>);
    module.def("join", &_join_<ValidatedVectorMultivariateFunctionModelDP,ValidatedScalarMultivariateFunctionModelDP>);
    module.def("join", &_join_<ValidatedVectorMultivariateFunctionModelDP,ValidatedVectorMultivariateFunctionModelDP>);

    module.def("combine", &_combine_<ValidatedVectorMultivariateFunctionModelDP,ValidatedVectorMultivariateFunctionModelDP>);
        module.def("combine", &_combine_<ValidatedVectorMultivariateFunctionModelDP,ValidatedScalarMultivariateFunctionModelDP>);

    module.def("antiderivative", &_antiderivative_<ValidatedVectorMultivariateFunctionModelDP,SizeType,NumericType>);
    module.def("antiderivative", &_antiderivative_<ValidatedVectorMultivariateFunctionModelDP,SizeType,ValidatedNumber>);

//    to_python< List<ValidatedVectorMultivariateFunctionModelDP> >();
}


Void export_scalar_taylor_function(pybind11::module& module)
{
    using FunctionModelType = ValidatedScalarMultivariateTaylorFunctionModelDP;
//    using FunctionType = ScalarMultivariateFunction<Paradigm<FunctionModelType>>;
    using GenericFunctionType = ScalarMultivariateFunction<ValidatedTag>;
    using NumericType = NumericType<FunctionModelType>;
    using GenericNumericType = GenericType<NumericType>;

    typedef ValidatedScalarMultivariateTaylorFunctionModelDP F;
    typedef ValidatedVectorMultivariateTaylorFunctionModelDP VF;
    typedef typename F::DomainType D;
    typedef typename F::NumericType X;
    typedef Vector<X> VX;
    typedef SizeType I;

    Tag<GenericNumericType> generic_number_tag;
    Tag<GenericFunctionType> generic_function_tag;

    pybind11::class_<ValidatedScalarMultivariateTaylorFunctionModelDP> scalar_taylor_function_class(module,"ValidatedScalarMultivariateTaylorFunctionModelDP");
    scalar_taylor_function_class.def(pybind11::init<ValidatedScalarMultivariateTaylorFunctionModelDP>());
    scalar_taylor_function_class.def(pybind11::init<ExactBoxType,ValidatedTaylorModelDP>());
    scalar_taylor_function_class.def(pybind11::init< ExactBoxType,SweeperDP >());
    scalar_taylor_function_class.def(pybind11::init< ExactBoxType, const EffectiveScalarMultivariateFunction&,SweeperDP >());
    scalar_taylor_function_class.def(pybind11::init< ExactBoxType, Expansion<MultiIndex,FloatDPValue>, FloatDPError, SweeperDP >());
    scalar_taylor_function_class.def("error", (const FloatDPError(ValidatedScalarMultivariateTaylorFunctionModelDP::*)()const) &ValidatedScalarMultivariateTaylorFunctionModelDP::error);
    scalar_taylor_function_class.def("set_error", (Void(ValidatedScalarMultivariateTaylorFunctionModelDP::*)(const FloatDPError&)) &ValidatedScalarMultivariateTaylorFunctionModelDP::set_error);
    scalar_taylor_function_class.def("argument_size", &F::argument_size);
    scalar_taylor_function_class.def("domain", &F::domain);
    scalar_taylor_function_class.def("codomain", &F::codomain);
    scalar_taylor_function_class.def("range", &F::range);
    scalar_taylor_function_class.def("model", (const ValidatedTaylorModelDP&(ValidatedScalarMultivariateTaylorFunctionModelDP::*)()const)&ValidatedScalarMultivariateTaylorFunctionModelDP::model);
    scalar_taylor_function_class.def("polynomial", (MultivariatePolynomial<NumericType>(ValidatedScalarMultivariateTaylorFunctionModelDP::*)()const)&ValidatedScalarMultivariateTaylorFunctionModelDP::polynomial);
    scalar_taylor_function_class.def("number_of_nonzeros", (SizeType(ValidatedScalarMultivariateTaylorFunctionModelDP::*)()const)&ValidatedScalarMultivariateTaylorFunctionModelDP::number_of_nonzeros);
    scalar_taylor_function_class.def("__getitem__", &__getitem__<ValidatedScalarMultivariateTaylorFunctionModelDP,MultiIndex,FloatDPValue>);
    scalar_taylor_function_class.def("__setitem__",&__setitem__<ValidatedScalarMultivariateTaylorFunctionModelDP,MultiIndex,FloatDPValue>);

    define_elementary_algebra(module,scalar_taylor_function_class);
    define_inplace_algebra(module,scalar_taylor_function_class);
    define_mixed_arithmetic(module,scalar_taylor_function_class,generic_number_tag);
    define_mixed_arithmetic(module,scalar_taylor_function_class,generic_function_tag);
    define_transcendental(module,scalar_taylor_function_class);
    define_lattice(module,scalar_taylor_function_class);

    scalar_taylor_function_class.def("__str__", &__cstr__<F>);
    scalar_taylor_function_class.def("__repr__", &__crepr__<F>);

    scalar_taylor_function_class.def("value", (const FloatDPValue(ValidatedScalarMultivariateTaylorFunctionModelDP::*)()const) &ValidatedScalarMultivariateTaylorFunctionModelDP::value);
    scalar_taylor_function_class.def("clobber", (Void(ValidatedScalarMultivariateTaylorFunctionModelDP::*)()) &ValidatedScalarMultivariateTaylorFunctionModelDP::clobber);
    scalar_taylor_function_class.def("set_properties",&ValidatedScalarMultivariateTaylorFunctionModelDP::set_properties);
    scalar_taylor_function_class.def("properties",&ValidatedScalarMultivariateTaylorFunctionModelDP::properties);
    scalar_taylor_function_class.def("__call__", (FloatDPApproximation(ValidatedScalarMultivariateTaylorFunctionModelDP::*)(const Vector<FloatDPApproximation>&)const) &ValidatedScalarMultivariateTaylorFunctionModelDP::operator());
    scalar_taylor_function_class.def("__call__", (FloatDPBounds(ValidatedScalarMultivariateTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedScalarMultivariateTaylorFunctionModelDP::operator());
//      scalar_taylor_function_class.def("gradient", (Covector<FloatDPBounds>(ValidatedScalarMultivariateTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedScalarMultivariateTaylorFunctionModelDP::gradient);
    scalar_taylor_function_class.def("function", (ValidatedScalarMultivariateFunction(ValidatedScalarMultivariateTaylorFunctionModelDP::*)()const) &ValidatedScalarMultivariateTaylorFunctionModelDP::function);
    scalar_taylor_function_class.def("polynomial", (MultivariatePolynomial<FloatDPBounds>(ValidatedScalarMultivariateTaylorFunctionModelDP::*)()const) &ValidatedScalarMultivariateTaylorFunctionModelDP::polynomial);
    scalar_taylor_function_class.def("restriction", &_restriction_<F,D>);
    module.def("restrict", &_restriction_<F,D>);
//    scalar_taylor_function_class.def("extension",&_extension_<F,D>);

    scalar_taylor_function_class.def_static("zero",(ValidatedScalarMultivariateTaylorFunctionModelDP(*)(const ExactBoxType&,SweeperDP))&ValidatedScalarMultivariateTaylorFunctionModelDP::zero);
    scalar_taylor_function_class.def_static("constant",(ValidatedScalarMultivariateTaylorFunctionModelDP(*)(const ExactBoxType&,const NumericType&,SweeperDP))&ValidatedScalarMultivariateTaylorFunctionModelDP::constant);
    scalar_taylor_function_class.def_static("constant",[](const ExactBoxType& bx, const ValidatedNumber& c,SweeperDP swp){return ValidatedScalarMultivariateTaylorFunctionModelDP::constant(bx,NumericType(c,swp.precision()),swp);});
    scalar_taylor_function_class.def_static("coordinate",(ValidatedScalarMultivariateTaylorFunctionModelDP(*)(const ExactBoxType&,SizeType,SweeperDP))&ValidatedScalarMultivariateTaylorFunctionModelDP::coordinate);


    module.def("restriction", &_restriction_<F,D>);
    module.def("join",(VF(*)(F,F)) &_join_<F,F>);
//    module.def("extension",&_extension_<F,D>);
    module.def("embed", &_embed_<D,F,D>);
//    module.def("split", &_split_<F,I>);
    module.def("evaluate", &_evaluate_<F,VX>);
    module.def("partial_evaluate", &_partial_evaluate_<F,I,X>);
    module.def("midpoint", &_midpoint_<F>);
    module.def("derivative", &_derivative_<F,I>);
    module.def("antiderivative", &_antiderivative_<F,I>);
    module.def("antiderivative", &_antiderivative_<F,I,X>);

    module.def("inconsistent", &_inconsistent_<F,F>);
    module.def("refines", &_refines_<F,F>);
    module.def("refinement", &_refinement_<F,F>);


//    to_python< Vector<ValidatedScalarMultivariateTaylorFunctionModelDP> >();
}


Void export_vector_taylor_function(pybind11::module& module)
{

    //using VectorMultivariateFunctionModelType = ValidatedVectorMultivariateFunctionModelDP;
    using ScalarMultivariateFunctionModelType = ValidatedScalarMultivariateFunctionModelDP;
    using NumericType = NumericType<ScalarMultivariateFunctionModelType>;

    typedef SizeType I;
    typedef ValidatedScalarMultivariateFunction SFN;
    typedef ValidatedVectorMultivariateFunction VFN;
    typedef ValidatedScalarMultivariateTaylorFunctionModelDP SF;
    typedef ValidatedVectorMultivariateTaylorFunctionModelDP VF;
    typedef typename VF::DomainType D;
    typedef typename D::IntervalType Di;
    typedef typename VF::NumericType X;
    typedef Vector<X> VX;

    Tag<ValidatedScalarMultivariateFunctionModelDP> scalar_taylor_function_tag;
    Tag<Vector<NumericType>> number_vector_tag;

    pybind11::class_<ValidatedVectorMultivariateTaylorFunctionModelDP> vector_taylor_function_class(module,"ValidatedVectorMultivariateTaylorFunctionModelDP");
    vector_taylor_function_class.def( pybind11::init<ValidatedVectorMultivariateTaylorFunctionModelDP>());
//    vector_taylor_function_class.def( pybind11::init<Vector<ValidatedScalarMultivariateTaylorFunctionModelDP>>());
    vector_taylor_function_class.def( pybind11::init([](Array<ValidatedScalarMultivariateTaylorFunctionModelDP> ary){return ValidatedVectorMultivariateTaylorFunctionModelDP(Vector<ValidatedScalarMultivariateTaylorFunctionModelDP>(ary));}));
    vector_taylor_function_class.def( pybind11::init< SizeType, ExactBoxType, SweeperDP >());
    vector_taylor_function_class.def( pybind11::init< ExactBoxType,const EffectiveVectorMultivariateFunction&,SweeperDP >());
    vector_taylor_function_class.def(pybind11::init< ExactBoxType, Vector< Expansion<MultiIndex,FloatDPValue> >, Vector<FloatDPError>, SweeperDP >());
    vector_taylor_function_class.def( pybind11::init< Vector<ValidatedScalarMultivariateTaylorFunctionModelDP> >());
    vector_taylor_function_class.def("_len_", &ValidatedVectorMultivariateTaylorFunctionModelDP::result_size);
    vector_taylor_function_class.def("result_size", &ValidatedVectorMultivariateTaylorFunctionModelDP::result_size);
    vector_taylor_function_class.def("argument_size", &ValidatedVectorMultivariateTaylorFunctionModelDP::argument_size);
    vector_taylor_function_class.def("domain", &ValidatedVectorMultivariateTaylorFunctionModelDP::domain);
    vector_taylor_function_class.def("codomain", &ValidatedVectorMultivariateTaylorFunctionModelDP::codomain);
    // FIXME: Omitted since const and non-const versions
    // vector_taylor_function_class.def("models", &ValidatedVectorMultivariateTaylorFunctionModelDP::models);
    vector_taylor_function_class.def("centre", &ValidatedVectorMultivariateTaylorFunctionModelDP::centre);
    vector_taylor_function_class.def("range", &ValidatedVectorMultivariateTaylorFunctionModelDP::range);
    vector_taylor_function_class.def("errors", &ValidatedVectorMultivariateTaylorFunctionModelDP::errors);
    vector_taylor_function_class.def("clobber", (Void(ValidatedVectorMultivariateTaylorFunctionModelDP::*)()) &ValidatedVectorMultivariateTaylorFunctionModelDP::clobber);
    vector_taylor_function_class.def("set_properties",&ValidatedVectorMultivariateTaylorFunctionModelDP::set_properties);
    vector_taylor_function_class.def("properties",&ValidatedVectorMultivariateTaylorFunctionModelDP::properties);
    //FIXME: Omitted since project(...) fails
    //vector_taylor_function_class.def("__getslice__", &__getslice__<ValidatedVectorMultivariateTaylorFunctionModelDP,SizeType,SizeType,ValidatedVectorMultivariateTaylorFunctionModelDP>);
    vector_taylor_function_class.def("__getitem__", &__getitem__<ValidatedVectorMultivariateTaylorFunctionModelDP,SizeType,ValidatedScalarMultivariateTaylorFunctionModelDP>);
    vector_taylor_function_class.def("__setitem__",&__setitem__<ValidatedVectorMultivariateTaylorFunctionModelDP,SizeType,ValidatedScalarMultivariateTaylorFunctionModelDP>);

    define_vector_arithmetic(module,vector_taylor_function_class, scalar_taylor_function_tag);
    //FIXME: Not completely defined in C++ API
    //define_inplace_vector_arithmetic(module,vector_taylor_function_class, scalar_taylor_function_tag);
    define_inplace_mixed_vector_arithmetic(module,vector_taylor_function_class, number_vector_tag);

    vector_taylor_function_class.def("__str__", &__cstr__<ValidatedVectorMultivariateTaylorFunctionModelDP>);
    vector_taylor_function_class.def("__repr__", &__crepr__<ValidatedVectorMultivariateTaylorFunctionModelDP>);
    vector_taylor_function_class.def("clobber", (Void(ValidatedVectorMultivariateTaylorFunctionModelDP::*)()) &ValidatedVectorMultivariateTaylorFunctionModelDP::clobber);
    vector_taylor_function_class.def("__call__", (Vector<FloatDPApproximation>(ValidatedVectorMultivariateTaylorFunctionModelDP::*)(const Vector<FloatDPApproximation>&)const) &ValidatedVectorMultivariateTaylorFunctionModelDP::operator());
    vector_taylor_function_class.def("__call__", (Vector<FloatDPBounds>(ValidatedVectorMultivariateTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedVectorMultivariateTaylorFunctionModelDP::operator());
     vector_taylor_function_class.def("jacobian", (Matrix<FloatDPBounds>(ValidatedVectorMultivariateTaylorFunctionModelDP::*)(const Vector<FloatDPBounds>&)const) &ValidatedVectorMultivariateTaylorFunctionModelDP::jacobian);
    vector_taylor_function_class.def("polynomials", (Vector< MultivariatePolynomial<FloatDPBounds> >(ValidatedVectorMultivariateTaylorFunctionModelDP::*)()const) &ValidatedVectorMultivariateTaylorFunctionModelDP::polynomials);
    vector_taylor_function_class.def("function", (ValidatedVectorMultivariateFunction(ValidatedVectorMultivariateTaylorFunctionModelDP::*)()const) &ValidatedVectorMultivariateTaylorFunctionModelDP::function);

    vector_taylor_function_class.def_static("constant",(ValidatedVectorMultivariateTaylorFunctionModelDP(*)(const ExactBoxType&, const Vector<NumericType>&,SweeperDP))&ValidatedVectorMultivariateTaylorFunctionModelDP::constant);
    vector_taylor_function_class.def_static("constant",[](const ExactBoxType& bx, const Vector<ValidatedNumber>& c,SweeperDP swp){return ValidatedVectorMultivariateTaylorFunctionModelDP::constant(bx,Vector<NumericType>(c,swp.precision()),swp);});
    vector_taylor_function_class.def_static("identity",(ValidatedVectorMultivariateTaylorFunctionModelDP(*)(const ExactBoxType&,SweeperDP))&ValidatedVectorMultivariateTaylorFunctionModelDP::identity);

    module.def("inconsistent", &_inconsistent_<VF,VF>);
    module.def("refinement", &_refinement_<VF,VF>);
    module.def("refines", &_refines_<VF,VF>);

    module.def("join", &_join_<VF,SF>);
        module.def("join", &_join_<SF,SF>);
    module.def("combine", &_combine_<VF,SF>);
        module.def("combine", &_combine_<SF,SF>);
    module.def("embed", &_embed_<VF,D>);
        module.def("embed", &_embed_<D,VF,D>);

    module.def("restriction", &_restriction_<VF,I,Di>);
    module.def("restrict", &_restriction_<VF,I,Di>);
//    module.def("split", &_split_<VF,I>);

    module.def("evaluate", &_evaluate_<VF,VX>);
    module.def("partial_evaluate", &_partial_evaluate_<VF,I,X>);
    module.def("compose", &_compose_<VF,VF>);
    module.def("compose", &_compose_<SF,VF>);
    module.def("compose", &_compose_<SFN,VF>);
    module.def("compose", &_compose_<VFN,VF>);
    module.def("unchecked_compose", &_compose_<SF,VF>);
    module.def("unchecked_compose", &_compose_<VF,VF>);
    module.def("derivative", &_derivative_<VF,I>);
    module.def("antiderivative", &_antiderivative_<VF,I>);
    module.def("antiderivative", &_antiderivative_<VF,I,X>);

    module.def("compose", &_compose_<EffectiveScalarMultivariateFunction,VF>);
    module.def("compose", &_compose_<EffectiveVectorMultivariateFunction,VF>);

}


Void calculus_submodule(pybind11::module& module)
{
    export_expansion(module);
    export_sweeper(module);
/*
    export_approximate_taylor_model(module);
    export_validated_taylor_model(module);
*/
    export_scalar_function_model(module);
    export_vector_function_model(module);
    export_scalar_taylor_function(module);
    export_vector_taylor_function(module);
}


