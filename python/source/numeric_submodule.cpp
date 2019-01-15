/***************************************************************************
 *            numeric_submodule.cpp
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
#include <pybind11/operators.h>

#include "utilities.hpp"

#if defined(__GNUG__) && !defined(__clang__)
#  pragma GCC diagnostic ignored "-Wattributes"
#endif

#include "numeric/logical.hpp"
#include "numeric/integer.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "numeric/real.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/float-user.hpp"

namespace Ariadne {

extern template Nat Error<FloatDP>::output_places;
extern template Nat Error<FloatMP>::output_places;

template<> String class_name<DoublePrecision>() { return "DoublePrecision"; }
template<> String class_name<MultiplePrecision>() { return "MultiplePrecision"; }

template<class T> String numeric_class_tag();
template<> String numeric_class_tag<DoublePrecision>() { return "DP"; }
template<> String numeric_class_tag<MultiplePrecision>() { return "MP"; }

template<class T> struct PythonName { const char* get() const { return class_name<T>().c_str(); } };
template<class T> inline const char* python_name() { return PythonName<T>().get(); }

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Integer>& repr) {
    return os << "Integer("<<repr.reference()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Dyadic>& repr) {
    return os << "Dyadic("<<repr.reference().mantissa()<<","<<repr.reference().exponent()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Rational>& repr) {
    return os << "Rational("<<repr.reference().numerator()<<","<<repr.reference().denominator()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Real>& repr) {
    return os << "Real("<<repr.reference()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<RawFloatDP>& repr) {
    return os << "FloatDP("<<repr.reference()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPApproximation>& repr) {
    return os << "FloatDPApproximation("<<repr.reference().raw()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPBounds>& repr) {
    return os << "FloatDPBounds("<<repr.reference().lower().raw()<<","<<repr.reference().upper().raw()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPValue>& repr) {
    return os << "FloatDPValue("<<repr.reference().raw()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPUpperBound>& repr) {
    return os << "FloatDPUpperBound("<<repr.reference().raw()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<PositiveFloatDPUpperBound>& repr) {
    return os << "FloatDPError("<<repr.reference().raw()<<")"; }

ExactDouble exact(double d) { return ExactDouble(d); }
Dyadic cast_exact(double d) { return Dyadic(ExactDouble(d)); }

inline Dyadic operator/(Dyadic x, Two w) { return Dyadic(x/(w^1)); }


template<class L> Bool _decide_(L l) { return decide(l); }
template<class L> Bool _definitely_(L l) { return definitely(l); }
template<class L> Bool _possibly_(L l) { return possibly(l); }
template<class L> Bool _probably_(L l) { return probably(l); }

template<class L, class E=Effort> auto _check_(L const& l, E const& e) -> decltype(l.check(e)) { return l.check(e); }

template<class OP, class... TS> auto py_apply(TS const& ... ts) -> decltype(OP()(ts...)){ OP op; return op(ts...); }

template<class A> auto _is_nan_(A const& a) -> Bool { return is_nan(a); }
template<class A> auto _is_inf_(A const& a) -> Bool { return is_inf(a); }
template<class A> auto _is_finite_(A const& a) -> Bool { return is_finite(a); }
template<class A> auto _is_zero_(A const& a) -> Bool { return is_zero(a); }

template<class X> Void define_infinitary_checks(pybind11::module& module, pybind11::class_<X>& pyclass) {
    module.def("is_nan", &_is_nan_<X>);
    module.def("is_inf", &_is_inf_<X>);
    module.def("is_finite", &_is_finite_<X>);
    module.def("is_zero", &_is_zero_<X>);
}

template<class X,EnableIf<Not<Or<IsSame<X,FloatDP>,IsSame<X,FloatMP>>>> =dummy> Void define_infinitary(pybind11::module& module, pybind11::class_<X>& pyclass) {
    pyclass.def_static("nan", (X(*)()) &X::nan);
    pyclass.def_static("inf", (X(*)()) &X::inf);
    pyclass.def_static("inf", (X(*)(Sign)) &X::inf);

    define_infinitary_checks(module,pyclass);
}
template<class X,EnableIf<Or<IsSame<X,FloatDP>,IsSame<X,FloatMP>>> =dummy> Void define_infinitary(pybind11::module& module, pybind11::class_<X>& pyclass) {
    pyclass.def_static("nan", (X(*)(typename X::PrecisionType)) &X::nan);
    pyclass.def_static("inf", (X(*)(typename X::PrecisionType)) &X::inf);
    pyclass.def_static("inf", (X(*)(Sign,typename X::PrecisionType)) &X::inf);

    define_infinitary_checks(module,pyclass);
}

template<class X> X bounds_from_dict(pybind11::dict dct) {
    typedef decltype(declval<X>().lower_raw()) L;
    typedef decltype(declval<X>().upper_raw()) U;
    assert(dct.size()==1);
    pybind11::detail::dict_iterator::reference item = *dct.begin();
    pybind11::handle lh = item.first;
    pybind11::handle uh = item.second;
    L l = pybind11::cast<L>(lh);
    U u = pybind11::cast<U>(uh);
    return X(l,u);
}


} // namespace Ariadne


using namespace Ariadne;
using pymodule = pybind11::module;
using pybind11::init;
using pybind11::detail::self;
using pybind11::implicitly_convertible;

template<class P> void export_effective_logical(pymodule& module, std::string name)
{
    OutputStream& operator<<(OutputStream& os, LogicalType<P> l);

    pybind11::class_<LogicalType<P>> logical_class(module,name.c_str());
    logical_class.def(init<bool>());
    logical_class.def(init<LogicalType<P>>());
    logical_class.def("__str__", &__cstr__<LogicalType<P>>);
    logical_class.def("__repr__", &__cstr__<LogicalType<P>>);
    logical_class.def("check", &LogicalType<P>::check);
    define_logical(module,logical_class);
    module.def("check", &_check_<LogicalType<P>,Effort>);
}

template<class P> void export_logical(pymodule& module, std::string name)
{
    OutputStream& operator<<(OutputStream& os, LogicalType<P> l);
    pybind11::class_<LogicalType<P>> logical_class(module,name.c_str());
    logical_class.def(init<bool>());
    logical_class.def(init<LogicalType<P>>());
    logical_class.def("__str__", &__cstr__<LogicalType<P>>);
    logical_class.def("__repr__", &__cstr__<LogicalType<P>>);
    define_logical(module,logical_class);

    module.def("decide", (bool(*)(LogicalType<P>)) &_decide_);
    module.def("probably", (bool(*)(LogicalType<P>)) &_probably_);
    module.def("possibly", (bool(*)(LogicalType<P>)) &_possibly_);
    module.def("definitely", (bool(*)(LogicalType<P>)) &_definitely_);
}

template<> void export_logical<ExactTag>(pymodule& module, std::string name)
{
    typedef ExactTag P;
    OutputStream& operator<<(OutputStream& os, LogicalType<P> l);
    pybind11::class_<LogicalType<P>> logical_class(module,name.c_str());
    logical_class.def(init<bool>());
    logical_class.def(init<LogicalType<P>>());
    logical_class.def("__str__", &__cstr__<LogicalType<P>>);
    logical_class.def("__repr__", &__cstr__<LogicalType<P>>);
    define_logical(module,logical_class);

//    implicitly_convertible<LogicalType<ExactTag>,bool>();
}



Void export_logicals(pymodule& module) {
    export_logical<ExactTag>(module,"Boolean");
    export_effective_logical<EffectiveTag>(module,"Kleenean");
    export_effective_logical<EffectiveUpperTag>(module,"Sierpinskian");
    export_effective_logical<EffectiveLowerTag>(module,"NegatedSierpinskian");
    export_logical<ValidatedTag>(module,"ValidatedKleenean");
    export_logical<UpperTag>(module,"ValidatedUpperKleenean");
    export_logical<LowerTag>(module,"ValidatedLowerKleenean");
    export_logical<ApproximateTag>(module,"ApproximateKleenean");
}


void export_builtins(pymodule& module)
{
    pybind11::class_<ExactDouble> exact_double_class(module,"ExactDouble");
    exact_double_class.def(init<double>());
    exact_double_class.def(init<ExactDouble>());
    exact_double_class.def("__str__", &__cstr__<ExactDouble>);

    module.def("exact", (ExactDouble(*)(double)) &exact);

}


void export_integer(pymodule& module)
{
    pybind11::class_<Integer> integer_class(module,"Integer");
    integer_class.def(init<int>());
    integer_class.def(init<Integer>());
    integer_class.def("__str__", &__cstr__<Integer>);
    integer_class.def("__repr__", &__repr__<Integer>);

    define_arithmetic(module,integer_class);
    define_lattice(module,integer_class);
    define_comparisons(module,integer_class);
    module.def("pow", &_pow_<Integer,Nat>);

    implicitly_convertible<int,Integer>();
}

void export_dyadic(pymodule& module)
{
    pybind11::class_<Dyadic> dyadic_class(module,"Dyadic");
    dyadic_class.def(init<Int,Nat>());
    dyadic_class.def(init<Integer,Natural>());
    dyadic_class.def(init<ExactDouble>());
    dyadic_class.def(init<Integer>());
    dyadic_class.def(init<Dyadic>());
    dyadic_class.def(init<FloatDPValue>());
    dyadic_class.def(init<FloatMPValue>());
    dyadic_class.def("__str__", &__cstr__<Dyadic>);
    dyadic_class.def("__repr__", &__repr__<Dyadic>);
    dyadic_class.def("__rmul__", &__rmul__<Dyadic,Dyadic>);

    define_infinitary(module,dyadic_class);
    define_arithmetic(module,dyadic_class);
    define_lattice(module,dyadic_class);
    define_comparisons(module,dyadic_class);
    module.def("hlf", &_hlf_<Dyadic>);

    implicitly_convertible<int,Dyadic>();
    implicitly_convertible<Integer,Dyadic>();

    pybind11::class_<Two> two_class(module,"Two");
    two_class.def("__pow__", &__pow__<Two,Int>);
    two_class.def(__py_rdiv__, &__rdiv__<Two,Dyadic, Return<Dyadic>>);
    dyadic_class.def(__py_div__, &__div__<Dyadic,Two, Return<Dyadic>>);
    module.attr("two") = Two();

    pybind11::class_<TwoExp> two_exp_class(module,"TwoExp");
    two_exp_class.def(__py_rdiv__, &__rdiv__<TwoExp,Dyadic, Return<Dyadic>>);
    dyadic_class.def(__py_div__, &__div__<Dyadic,TwoExp, Return<Dyadic>>);

    module.def("cast_exact",(Dyadic(*)(double)) &cast_exact);
}

void export_decimal(pymodule& module)
{
    pybind11::class_<Decimal> decimal_class(module,"Decimal");
    decimal_class.def(init<Decimal>());
    decimal_class.def(init<Dyadic>());
    decimal_class.def(init<std::string>());
    decimal_class.def(init<double>());
    decimal_class.def("__str__", &__cstr__<Decimal>);
    decimal_class.def("__repr__", &__cstr__<Decimal>);
    implicitly_convertible<std::string,Decimal>();
}

void export_rational(pymodule& module)
{
    pybind11::class_<Rational> rational_class(module,"Rational");
    rational_class.def(init<int,int>());
    rational_class.def(init<Integer,Integer>());
    rational_class.def(init<int>());
    rational_class.def(init<Integer>());
    rational_class.def(init<Dyadic>());
    rational_class.def(init<Decimal>());
    rational_class.def(init<Rational>());
    rational_class.def("__str__", &__cstr__<Rational>);
    rational_class.def("__repr__", &__repr__<Rational>);

    define_infinitary(module,rational_class);

    define_arithmetic(module,rational_class);
    define_lattice(module,rational_class);
    define_comparisons(module,rational_class);
    module.def("rec", &_rec_<Rational>);

    implicitly_convertible<int,Rational>();
    implicitly_convertible<Integer,Rational>();
    implicitly_convertible<Dyadic,Rational>();
    implicitly_convertible<Decimal,Rational>();
}

void export_real(pymodule& module)
{
    pybind11::class_<Real> real_class(module,"Real");
    real_class.def(init<int>());
    real_class.def(init<Integer>());
    real_class.def(init<Dyadic>());
    real_class.def(init<Rational>());
    real_class.def(init<Real>());
    real_class.def("__str__", &__cstr__<Real>);
    real_class.def("__repr__", &__cstr__<Real>);

    real_class.def("get", (FloatDPBounds(Real::*)(DoublePrecision)const) &Real::get);
    real_class.def("get", (FloatMPBounds(Real::*)(MultiplePrecision)const) &Real::get);
    real_class.def("compute", (ValidatedReal(Real::*)(Effort)const) &Real::compute);
    real_class.def("compute", (ValidatedReal(Real::*)(Accuracy)const) &Real::compute);

    define_elementary(module,real_class);
    define_lattice(module,real_class);
    define_comparisons(module,real_class);

    implicitly_convertible<Int,Real>();
    implicitly_convertible<Integer,Real>();
    implicitly_convertible<Dyadic,Real>();
    implicitly_convertible<Decimal,Real>();
    implicitly_convertible<Rational,Real>();


    pybind11::class_<ValidatedReal> validated_real_class(module,"ValidatedReal");
    validated_real_class.def(init<DyadicBounds>());
    validated_real_class.def("__str__", &__cstr__<ValidatedReal>);
    validated_real_class.def("__repr__", &__cstr__<ValidatedReal>);

    validated_real_class.def("get", (DyadicBounds(ValidatedReal::*)()const) &ValidatedReal::get);
    validated_real_class.def("get", (FloatDPBounds(ValidatedReal::*)(DoublePrecision)const) &ValidatedReal::get);
    validated_real_class.def("get", (FloatMPBounds(ValidatedReal::*)(MultiplePrecision)const) &ValidatedReal::get);
}


template<class P> void export_number(pybind11::module& module)
{
    pybind11::class_<Number<P>> number_class(module,class_name<P>()+"Number");
    number_class.def(init<Rational>());

    define_self_arithmetic(module,number_class);
//    number_class.def("get", (FloatDPBounds(Number<P>::*)(ValidatedTag,DoublePrecision)const) &Number<P>::get);
//    number_class.def("get", (FloatDPApproximation(Number<P>::*)(ApproximateTag,DoublePrecision)const) &Number<P>::get);
//    number_class.def("get", (FloatMPBounds(Number<P>::*)(ValidatedTag,MultiplePrecision)const) &Number<P>::get);
//    number_class.def("get", (FloatMPApproximation(Number<P>::*)(ApproximateTag,MultiplePrecision)const) &Number<P>::get);
}

FloatDPBounds get(ExactNumber const& n, DoublePrecision const& pr) { return n.get(BoundedTag(),pr); }
FloatMPBounds get(ExactNumber const& n, MultiplePrecision const& pr) { return n.get(BoundedTag(),pr); }
FloatDPBounds get(EffectiveNumber const& n, DoublePrecision const& pr) { return n.get(BoundedTag(),pr); }
FloatMPBounds get(EffectiveNumber const& n, MultiplePrecision const& pr) { return n.get(BoundedTag(),pr); }
FloatDPBounds get(ValidatedNumber const& n, DoublePrecision const& pr) { return n.get(BoundedTag(),pr); }
FloatMPBounds get(ValidatedNumber const& n, MultiplePrecision const& pr) { return n.get(BoundedTag(),pr); }
FloatDPUpperBound get(ValidatedUpperNumber const& n, DoublePrecision const& pr) { return n.get(UpperTag(),pr); }
FloatMPUpperBound get(ValidatedUpperNumber const& n, MultiplePrecision const& pr) { return n.get(UpperTag(),pr); }
FloatDPLowerBound get(ValidatedLowerNumber const& n, DoublePrecision const& pr) { return n.get(LowerTag(),pr); }
FloatMPLowerBound get(ValidatedLowerNumber const& n, MultiplePrecision const& pr) { return n.get(LowerTag(),pr); }
FloatDPApproximation get(ApproximateNumber const& n, DoublePrecision const& pr) { return n.get(ApproximateTag(),pr); }
FloatMPApproximation get(ApproximateNumber const& n, MultiplePrecision const& pr) { std::cerr<<"get(AN,MP)\n";return n.get(ApproximateTag(),pr); }

void export_numbers(pymodule& module)
{
    pybind11::class_<ApproximateNumber> approximate_number_class(module,(class_name<ApproximateTag>()+"Number").c_str());
    approximate_number_class.def(init<double>());
    approximate_number_class.def(init<Rational>());
    approximate_number_class.def(init<Real>());
    approximate_number_class.def(init<ExactNumber>());
    approximate_number_class.def(init<EffectiveNumber>());
    approximate_number_class.def(init<ValidatedNumber>());
    approximate_number_class.def(init<ApproximateNumber>());
    approximate_number_class.def("get", (FloatDPApproximation(*)(ApproximateNumber const&, DoublePrecision const&)) &get);
    approximate_number_class.def("get", (FloatMPApproximation(*)(ApproximateNumber const&, MultiplePrecision const&)) &get);
    approximate_number_class.def("__str__", &__cstr__<ApproximateNumber>);
    approximate_number_class.def("__repr__", &__cstr__<ApproximateNumber>);

    // TODO: These exports should be with FloatApproximation
    approximate_number_class.def(init<FloatDPApproximation>());
    approximate_number_class.def(init<FloatMPApproximation>());

    define_elementary(module,approximate_number_class);
    define_lattice(module,approximate_number_class);


    pybind11::class_<ValidatedLowerNumber> lower_number_class(module,(class_name<LowerTag>()+"Number").c_str());
    lower_number_class.def(init<ValidatedNumber>());
    lower_number_class.def("get", (FloatDPLowerBound(*)(ValidatedLowerNumber const&, DoublePrecision const&)) &get);
    lower_number_class.def("get", (FloatMPLowerBound(*)(ValidatedLowerNumber const&, MultiplePrecision const&)) &get);
    lower_number_class.def("__str__", &__cstr__<ValidatedLowerNumber>);
    lower_number_class.def("__repr__", &__cstr__<ValidatedLowerNumber>);

    define_monotonic(module,lower_number_class);


    pybind11::class_<ValidatedUpperNumber> upper_number_class(module,(class_name<UpperTag>()+"Number").c_str());
    upper_number_class.def(init<ValidatedNumber>());
    upper_number_class.def("get", (FloatDPUpperBound(*)(ValidatedUpperNumber const&, DoublePrecision const&)) &get);
    upper_number_class.def("get", (FloatMPUpperBound(*)(ValidatedUpperNumber const&, MultiplePrecision const&)) &get);
    upper_number_class.def("__str__", &__cstr__<ValidatedUpperNumber>);
    upper_number_class.def("__repr__", &__cstr__<ValidatedUpperNumber>);

    define_monotonic(module,upper_number_class);


    pybind11::class_<ValidatedNumber> validated_number_class(module,(class_name<ValidatedTag>()+"Number").c_str());
    validated_number_class.def(init<Rational>());
    validated_number_class.def(init<Real>());
    validated_number_class.def(init<DyadicBounds>());
    validated_number_class.def(init<ExactNumber>());
    validated_number_class.def(init<EffectiveNumber>());
    validated_number_class.def(init<ValidatedNumber>());
    validated_number_class.def("get", (FloatDPBounds(*)(ValidatedNumber const&, DoublePrecision const&)) &get);
    validated_number_class.def("get", (FloatMPBounds(*)(ValidatedNumber const&, MultiplePrecision const&)) &get);
    validated_number_class.def("__str__", &__cstr__<ValidatedNumber>);
    validated_number_class.def("__repr__", &__cstr__<ValidatedNumber>);

    define_elementary(module,validated_number_class);
    validated_number_class.def(pybind11::init([](pybind11::dict pydct){return ValidatedNumber(bounds_from_dict<DyadicBounds>(pydct));}));
    implicitly_convertible<pybind11::dict,ValidatedNumber>();

// TODO: These exports should be with FloatBounds
    validated_number_class.def(init<FloatDPBounds>());
    validated_number_class.def(init<FloatMPBounds>());

    pybind11::class_<EffectiveNumber> effective_number_class(module,(class_name<EffectiveTag>()+"Number").c_str());
    effective_number_class.def(init<Rational>());
    effective_number_class.def(init<Real>());
    effective_number_class.def(init<ExactNumber>());
    effective_number_class.def(init<EffectiveNumber>());
    effective_number_class.def("get", (FloatDPBounds(*)(EffectiveNumber const&, DoublePrecision const&)) &get);
    effective_number_class.def("get", (FloatMPBounds(*)(EffectiveNumber const&, MultiplePrecision const&)) &get);
    effective_number_class.def("__str__", &__cstr__<EffectiveNumber>);
    effective_number_class.def("__repr__", &__cstr__<EffectiveNumber>);

    define_elementary(module,effective_number_class);

    pybind11::class_<ExactNumber> exact_number_class(module,(class_name<ExactTag>()+"Number").c_str());
    exact_number_class.def(init<Rational>());
    exact_number_class.def(init<ExactNumber>());
    exact_number_class.def("get", (FloatDPBounds(*)(ExactNumber const&, DoublePrecision const&)) &get);
    exact_number_class.def("get", (FloatMPBounds(*)(ExactNumber const&, MultiplePrecision const&)) &get);
    exact_number_class.def("__str__", &__cstr__<ExactNumber>);
    exact_number_class.def("__repr__", &__cstr__<ExactNumber>);


    implicitly_convertible<ValidatedNumber,ApproximateNumber>();
    implicitly_convertible<ValidatedNumber,ValidatedLowerNumber>();
    implicitly_convertible<ValidatedNumber,ValidatedUpperNumber>();
    implicitly_convertible<EffectiveNumber,ValidatedNumber>();
    implicitly_convertible<ExactNumber,ValidatedNumber>();
    implicitly_convertible<ExactNumber,EffectiveNumber>();

    implicitly_convertible<DyadicBounds,ValidatedNumber>();

    implicitly_convertible<Real,EffectiveNumber>();
    implicitly_convertible<Rational,ExactNumber>();
    implicitly_convertible<Rational,ValidatedNumber>();
    implicitly_convertible<Decimal,ExactNumber>();
    implicitly_convertible<Decimal,ValidatedNumber>();
    implicitly_convertible<Dyadic,ExactNumber>();
    implicitly_convertible<Dyadic,ValidatedNumber>();
    implicitly_convertible<Integer,ExactNumber>();
    implicitly_convertible<Integer,ValidatedNumber>();

    implicitly_convertible<Int,ApproximateNumber>();
    implicitly_convertible<Int,ValidatedLowerNumber>();
    implicitly_convertible<Int,ValidatedUpperNumber>();
    implicitly_convertible<Int,ValidatedNumber>();
    implicitly_convertible<Int,EffectiveNumber>();
    implicitly_convertible<Int,ExactNumber>();

    implicitly_convertible<double,ApproximateNumber>();
}




void export_dyadic_bounds(pymodule& module)
{
    pybind11::class_<DyadicBounds> dyadic_bounds_class(module,"DyadicBounds");
    dyadic_bounds_class.def(init<DyadicBounds>());
    dyadic_bounds_class.def(init<Dyadic,Dyadic>());
    dyadic_bounds_class.def(pybind11::init([](pybind11::dict pydct){return bounds_from_dict<DyadicBounds>(pydct);}));
    implicitly_convertible<pybind11::dict,DyadicBounds>();

    dyadic_bounds_class.def("__str__", &__cstr__<DyadicBounds>);
    dyadic_bounds_class.def("__repr__", &__cstr__<DyadicBounds>);

}

namespace Ariadne {
const Rounding round_upward = Rounding(upward);
const Rounding round_downward = Rounding(downward);
const Rounding round_to_nearest = Rounding(to_nearest);
}

void export_rounding_mode(pymodule& module) {
    pybind11::class_<Rounding> rounding_mode_class(module,"Rounding");
    rounding_mode_class.def(init<Rounding>());
    rounding_mode_class.def("__str__", &__cstr__<Rounding>);
    module.attr("upward") = round_upward;
    module.attr("downward") = round_downward;
    module.attr("to_near") = round_to_nearest;
    module.attr("up") = round_upward;
    module.attr("down") = round_downward;
    module.attr("near") = round_to_nearest;
}

template<class PR> void export_raw_float(pymodule& module)
{
    typedef RawFloat<PR> F;

//    typedef typename F::RoundingModeType RND;
//    implicitly_convertible<Rounding,RND>();
    typedef Rounding RND;

    pybind11::class_<F> raw_float_class(module,("Float"+numeric_class_tag<PR>()).c_str());
    raw_float_class.def(init<double,PR>());
    raw_float_class.def(init<Dyadic,PR>());
    raw_float_class.def(init<Rational,RND,PR>());

    define_infinitary(module,raw_float_class);

    raw_float_class.def_static("eps", (F(*)(PR)) &F::eps);
    raw_float_class.def_static("max", (F(*)(PR)) &F::max);
    raw_float_class.def_static("min", (F(*)(PR)) &F::min);

    raw_float_class.def("__str__", &__cstr__<RawFloat<PR>>);
    raw_float_class.def("__repr__", &__cstr__<RawFloat<PR>>);

    module.def("nul", &_nul_<F>);
    module.def("pos", &_pos_<F>);
    module.def("neg", &_neg_<F>);
    module.def("hlf", &_hlf_<F>);

    module.def("nul", &_nul_<RND,F>);
    module.def("pos", &_pos_<RND,F>);
    module.def("neg", &_neg_<RND,F>);
    module.def("hlf", &_hlf_<RND,F>);
    module.def("sqr", &_sqr_<RND,F>);
    module.def("rec", &_rec_<RND,F>);
    module.def("add", &_add_<RND,F,F>);
    module.def("sub", &_sub_<RND,F,F>);
    module.def("mul", &_mul_<RND,F,F>);
    module.def("div", &_div_<RND,F,F>);
    module.def("fma", &_fma_<RND,F,F,F>);
    module.def("pow", &_pow_<RND,F,Int>);
    module.def("sqrt", &_sqrt_<RND,F>);
    module.def("exp", &_exp_<RND,F>);
    module.def("log", &_log_<RND,F>);
    module.def("sin", &_sin_<RND,F>);
    module.def("cos", &_cos_<RND,F>);
    module.def("tan", &_tan_<RND,F>);
    module.def("atan", &_atan_<RND,F>);

    module.def("abs", &_abs_<F>);
    module.def("max", &_min_<F,F>);
    module.def("min", &_max_<F,F>);

    define_comparisons<F>(module,raw_float_class);

}

template<class PR> void export_float_value(pymodule& module)
{
    pybind11::class_<FloatValue<PR>> float_value_class(module,("Float"+numeric_class_tag<PR>()+"Value").c_str());
    float_value_class.def(init<PR>());
    float_value_class.def(init<RawFloat<PR>>());
    float_value_class.def(init<int>());
    float_value_class.def(init<double>());
//    float_value_class.def(init<ExactDouble>());
    float_value_class.def(init<Integer,PR>());
    float_value_class.def(init<Dyadic,PR>());
    float_value_class.def(init<FloatValue<PR>>());

    float_value_class.def("__pos__", &__pos__<FloatValue<PR>>);
    float_value_class.def("__neg__", &__neg__<FloatValue<PR>>);

    float_value_class.def("precision", &FloatValue<PR>::precision);
    float_value_class.def("get_d",&FloatValue<PR>::get_d);

    define_mixed_arithmetic<FloatValue<PR>,Int>(module,float_value_class);
    define_arithmetic<FloatValue<PR>>(module,float_value_class);
    //float_value_class.define_mixed_arithmetic<ApproximateNumericType>();
    //float_value_class.define_mixed_arithmetic<LowerNumericType>();
    //float_value_class.define_mixed_arithmetic<UpperNumericType>();
    //float_value_class.define_mixed_arithmetic<ValidatedNumericType>();
    define_lattice(module,float_value_class);
    define_comparisons(module,float_value_class);

    float_value_class.def("precision", &FloatValue<PR>::precision);
    float_value_class.def("raw", (RawFloat<PR>const&(FloatValue<PR>::*)()const)&FloatValue<PR>::raw);
    float_value_class.def("get_d",&FloatValue<PR>::get_d);

    float_value_class.def("__str__", &__cstr__<FloatValue<PR>>);
    float_value_class.def("__repr__", &__cstr__<FloatValue<PR>>);

    float_value_class.def_static("set_output_places",&FloatValue<PR>::set_output_places);

    implicitly_convertible<int,FloatValue<PR>>();
    implicitly_convertible<double,FloatValue<PR>>();
    implicitly_convertible<FloatValue<PR>,ExactNumber>();
    implicitly_convertible<FloatValue<PR>,ValidatedNumber>();
    implicitly_convertible<FloatValue<PR>,ApproximateNumber>();
}

template<class PRE> void export_float_error(pymodule& module)
{
    pybind11::class_<FloatError<PRE>> float_error_class(module,("Float"+numeric_class_tag<PRE>()+"Error").c_str());
    float_error_class.def(init<PRE>());
    float_error_class.def(init<RawFloat<PRE>>());
    float_error_class.def(init<uint>());
    float_error_class.def(init<double>());
    float_error_class.def(init<FloatError<PRE>>());

    float_error_class.def("__pos__", &__pos__<FloatError<PRE>>, pybind11::is_operator());
    float_error_class.def("__add__", &__add__<FloatError<PRE>,FloatError<PRE>>, pybind11::is_operator());
    float_error_class.def("__mul__", &__mul__<FloatError<PRE>,FloatError<PRE>>, pybind11::is_operator());

    module.def("max", &_max_<FloatError<PRE>,FloatError<PRE>>);
    module.def("min", &_min_<FloatError<PRE>,FloatError<PRE>>);

    float_error_class.def("raw",(RawFloat<PRE>const&(FloatError<PRE>::*)()const)&FloatError<PRE>::raw);
    float_error_class.def("__str__", &__cstr__<FloatError<PRE>>);
    float_error_class.def("__repr__", &__cstr__<FloatError<PRE>>);

    float_error_class.def("precision", &FloatError<PRE>::precision);

    module.def("log2", (FloatUpperBound<PRE>(*)(FloatError<PRE>const&)) &_log2_);
    float_error_class.def_static("set_output_places",&FloatError<PRE>::set_output_places);
}

template<class PR, class PRE=PR> void export_float_ball(pymodule& module)
{
    String tag = IsSame<PR,PRE>::value ? numeric_class_tag<PR>() : numeric_class_tag<PR>()+numeric_class_tag<PRE>();
    pybind11::class_<FloatBall<PR,PRE>> float_ball_class(module,("Float"+tag+"Ball").c_str());
    float_ball_class.def(init<PR,PRE>());
    float_ball_class.def(init<RawFloat<PR>,RawFloat<PRE>>());
    float_ball_class.def(init<FloatValue<PR>,FloatError<PRE>>());
    float_ball_class.def(init<Real,PR>());

    float_ball_class.def(init<ExactDouble,PR>());
    float_ball_class.def(init<ValidatedNumber,PR>());
    float_ball_class.def(init<FloatValue<PR>>());
    float_ball_class.def(init<FloatBall<PR,PRE>>());
    float_ball_class.def(init<FloatBounds<PR>>());

    float_ball_class.def("value", &FloatBall<PR,PRE>::value);
    float_ball_class.def("error", &FloatBall<PR,PRE>::error);
    float_ball_class.def("lower", &FloatBall<PR,PRE>::lower);
    float_ball_class.def("upper", &FloatBall<PR,PRE>::upper);

    float_ball_class.def("precision", &FloatBall<PR,PRE>::precision);
    float_ball_class.def("error_precision", &FloatBall<PR,PRE>::error_precision);

    define_mixed_arithmetic<FloatBall<PR,PRE>,ValidatedNumber>(module,float_ball_class);
    define_mixed_arithmetic<FloatBall<PR,PRE>,Int>(module,float_ball_class);
    define_elementary<FloatBall<PR,PRE>>(module,float_ball_class);
    define_lattice<FloatBall<PR,PRE>>(module,float_ball_class);
    define_mixed_lattice<FloatBall<PR,PRE>,ValidatedNumber>(module,float_ball_class);
    define_comparisons<FloatBall<PR,PRE>>(module,float_ball_class);
    define_mixed_comparisons<FloatBall<PR,PRE>,ValidatedNumber>(module,float_ball_class);
//    float_ball_class.define_mixed_arithmetic<FloatBall<PR,PRE>>();
//    float_ball_class.define_mixed_arithmetic<ApproximateNumericType>();
//    float_ball_class.define_mixed_arithmetic<LowerNumericType>();
//    float_ball_class.define_mixed_arithmetic<UpperNumericType>();
//    float_ball_class.define_mixed_arithmetic<ValidatedNumericType>();

    float_ball_class.def("__str__", &__cstr__<FloatBall<PR,PRE>>);
    float_ball_class.def("__repr__", &__cstr__<FloatBall<PR,PRE>>);
}



template<class PR> void export_float_bounds(pymodule& module)
{
    pybind11::class_<FloatBounds<PR>> float_bounds_class(module,("Float"+numeric_class_tag<PR>()+"Bounds").c_str());
    float_bounds_class.def(init<PR>());
    float_bounds_class.def(init<RawFloat<PR>>());
    float_bounds_class.def(init<RawFloat<PR>,RawFloat<PR>>());
    float_bounds_class.def(init<double,double>());
    float_bounds_class.def(init<FloatLowerBound<PR>,FloatUpperBound<PR>>());
    float_bounds_class.def(init<Real,Real,PR>());
    float_bounds_class.def(init<Real,PR>());

    float_bounds_class.def(init<ExactDouble,PR>());
    float_bounds_class.def(init<ValidatedNumber,PR>());
    float_bounds_class.def(init<FloatValue<PR>>());
    float_bounds_class.def(init<FloatBall<PR>>());
    float_bounds_class.def(init<FloatBounds<PR>>());

    float_bounds_class.def("__str__", &__cstr__<FloatBounds<PR>>);
    float_bounds_class.def("__repr__", &__cstr__<FloatBounds<PR>>);

    float_bounds_class.def("precision", &FloatBounds<PR>::precision);
    float_bounds_class.def("lower", &FloatBounds<PR>::lower);
    float_bounds_class.def("upper", &FloatBounds<PR>::upper);
    float_bounds_class.def("value", &FloatBounds<PR>::value);
    float_bounds_class.def("error", &FloatBounds<PR>::error);


    define_mixed_arithmetic<FloatBounds<PR>,ValidatedNumber>(module,float_bounds_class);
    define_mixed_arithmetic<FloatBounds<PR>,Int>(module,float_bounds_class);
    define_mixed_arithmetic<FloatBounds<PR>,FloatBall<PR>>(module,float_bounds_class);
    define_elementary<FloatBounds<PR>>(module,float_bounds_class);
    define_lattice<FloatBounds<PR>>(module,float_bounds_class);
    define_mixed_lattice<FloatBounds<PR>,ValidatedNumber>(module,float_bounds_class);
    define_comparisons<FloatBounds<PR>>(module,float_bounds_class);
    define_mixed_comparisons<FloatBounds<PR>,ValidatedNumber>(module,float_bounds_class);
//    float_bounds_class.define_mixed_arithmetic<ApproximateNumericType>();
//    float_bounds_class.define_mixed_arithmetic<LowerNumericType>();
//    float_bounds_class.define_mixed_arithmetic<UpperNumericType>();
//    float_bounds_class.define_mixed_arithmetic<ValidatedNumericType>();

    module.def("refinement", &_refinement_<FloatBounds<PR>,FloatBounds<PR>>);
    module.def("refines", &_refines_<FloatBounds<PR>,FloatBounds<PR>>);
    module.def("inconsistent", &_inconsistent_<FloatBounds<PR>,FloatBounds<PR>>);

    implicitly_convertible<FloatValue<PR>,FloatBounds<PR>>();
    implicitly_convertible<FloatBall<PR>,FloatBounds<PR>>();
    implicitly_convertible<FloatBounds<PR>,ValidatedNumber>();
    implicitly_convertible<FloatBounds<PR>,ApproximateNumber>();

    float_bounds_class.def_static("set_output_places",&FloatBounds<PR>::set_output_places);
}

template<class PR> void export_float_upper_bound(pymodule& module)
{
    pybind11::class_<FloatUpperBound<PR>> float_upper_bound_class(module,("Float"+numeric_class_tag<PR>()+"UpperBound").c_str());
    float_upper_bound_class.def(init<PR>());
    float_upper_bound_class.def(init<RawFloat<PR>>());
    float_upper_bound_class.def(init<int>());
    float_upper_bound_class.def(init<double>());
    float_upper_bound_class.def(init<Real,PR>());

    float_upper_bound_class.def(init<FloatValue<PR>>());
    float_upper_bound_class.def(init<FloatBall<PR>>());
    float_upper_bound_class.def(init<FloatBounds<PR>>());
    float_upper_bound_class.def(init<FloatUpperBound<PR>>());
    float_upper_bound_class.def(init<ValidatedUpperNumber,PR>());
    float_upper_bound_class.def("__str__", &__cstr__<FloatUpperBound<PR>>);
    float_upper_bound_class.def("__repr__", &__cstr__<FloatUpperBound<PR>>);

    float_upper_bound_class.def("precision", &FloatUpperBound<PR>::precision);
    float_upper_bound_class.def("raw", (RawFloat<PR>const&(FloatUpperBound<PR>::*)()const)&FloatUpperBound<PR>::raw);

    define_monotonic(module,float_upper_bound_class);
    define_mixed_monotonic<FloatUpperBound<PR>,Int>(module,float_upper_bound_class);
    define_lattice<FloatUpperBound<PR>>(module,float_upper_bound_class);
    define_mixed_lattice<FloatUpperBound<PR>,ValidatedUpperNumber>(module,float_upper_bound_class);

    float_upper_bound_class.def("__lt__", &__lt__<FloatUpperBound<PR>,FloatLowerBound<PR>>);
    float_upper_bound_class.def("__le__", &__le__<FloatUpperBound<PR>,FloatLowerBound<PR>>);

//    float_upper_bound_class.def(self > FloatBounds<PR>());
//    float_upper_bound_class.def(self > FloatApproximation<PR>());
//    float_upper_bound_class.def(self < FloatLowerBound<PR>());
//    float_upper_bound_class.def(self >= FloatBounds<PR>());
//    float_upper_bound_class.def(self >= FloatApproximation<PR>());
//    float_upper_bound_class.def(self <= FloatLowerBound<PR>());

//    float_upper_bound_class.define_mixed_arithmetic<ApproximateNumericType>();
//    float_upper_bound_class.def(UpperNumericType() + self);
//    float_upper_bound_class.def(LowerNumericType() - self);
//    float_upper_bound_class.def(self + UpperNumericType());
//    float_upper_bound_class.def(self - LowerNumericType());

    implicitly_convertible<FloatBounds<PR>,FloatUpperBound<PR>>();
}

template<class PR> void export_float_lower_bound(pymodule& module)
{
    pybind11::class_<FloatLowerBound<PR>> float_lower_bound_class(module,("Float"+numeric_class_tag<PR>()+"LowerBound").c_str());
    float_lower_bound_class.def(init<PR>());
    float_lower_bound_class.def(init<RawFloat<PR>>());
    float_lower_bound_class.def(init<int>());
    float_lower_bound_class.def(init<double>());
    float_lower_bound_class.def(init<Real,PR>());

    float_lower_bound_class.def(init<FloatValue<PR>>());
    float_lower_bound_class.def(init<FloatBall<PR>>());
    float_lower_bound_class.def(init<FloatBounds<PR>>());
    float_lower_bound_class.def(init<FloatLowerBound<PR>>());
    float_lower_bound_class.def(init<ValidatedLowerNumber,PR>());
    float_lower_bound_class.def("__str__", &__cstr__<FloatLowerBound<PR>>);
    float_lower_bound_class.def("__repr__", &__cstr__<FloatLowerBound<PR>>);

    float_lower_bound_class.def("precision", &FloatLowerBound<PR>::precision);

    define_monotonic(module,float_lower_bound_class);
    define_mixed_monotonic<FloatLowerBound<PR>,Int>(module,float_lower_bound_class);
    define_lattice<FloatLowerBound<PR>>(module,float_lower_bound_class);
    define_mixed_lattice<FloatLowerBound<PR>,ValidatedLowerNumber>(module,float_lower_bound_class);

    float_lower_bound_class.def("__lt__", &__lt__<FloatLowerBound<PR>,FloatUpperBound<PR>>);
    float_lower_bound_class.def("__le__", &__le__<FloatLowerBound<PR>,FloatUpperBound<PR>>);

//    float_lower_bound_class.def(self < FloatBounds<PR>());
//    float_lower_bound_class.def(self < FloatApproximation<PR>());
//    float_lower_bound_class.def(self > FloatUpperBound<PR>());
//    float_lower_bound_class.def(self <= FloatBounds<PR>());
//    float_lower_bound_class.def(self <= FloatApproximation<PR>());
//    float_lower_bound_class.def(self >= FloatUpperBound<PR>());
//    float_lower_bound_class.def("raw", (RawFloat<PR>const&(FloatLowerBound<PR>::*)()const)&FloatLowerBound<PR>::raw);

//    float_lower_bound_class.define_mixed_arithmetic<ApproximateNumericType>();

//    float_lower_bound_class.def(LowerNumericType() + self);
//    float_lower_bound_class.def(UpperNumericType() - self);
//    float_lower_bound_class.def(self + LowerNumericType());
//    float_lower_bound_class.def(self - UpperNumericType());

    implicitly_convertible<FloatBounds<PR>,FloatLowerBound<PR>>();
}



template<class PR> void export_float_approximation(pymodule& module)
{
    pybind11::class_<FloatApproximation<PR>> float_approximation_class(module,("Float"+numeric_class_tag<PR>()+"Approximation").c_str());

    if(IsSame<PR,DoublePrecision>::value) {
        float_approximation_class.def(init<double>());
    }
    float_approximation_class.def(init<PR>());
    float_approximation_class.def(init<double,PR>());
    float_approximation_class.def(init<RawFloat<PR>>());
    float_approximation_class.def(init<Real,PR>());

    float_approximation_class.def(init<FloatValue<PR>>());
    float_approximation_class.def(init<FloatBall<PR>>());
    float_approximation_class.def(init<FloatBounds<PR>>());
    float_approximation_class.def(init<FloatLowerBound<PR>>());
    float_approximation_class.def(init<FloatUpperBound<PR>>());
    float_approximation_class.def(init<FloatApproximation<PR>>());
    float_approximation_class.def(init<ApproximateNumber,PR>());

    float_approximation_class.def("precision", &FloatApproximation<PR>::precision);
    float_approximation_class.def("raw", (RawFloat<PR>const&(FloatApproximation<PR>::*)()const)&FloatApproximation<PR>::raw);
    float_approximation_class.def("get_d", &FloatApproximation<PR>::get_d);

    define_mixed_arithmetic<FloatApproximation<PR>,double>(module,float_approximation_class);
    define_mixed_arithmetic<FloatApproximation<PR>,ApproximateNumber>(module,float_approximation_class);
    define_elementary<FloatApproximation<PR>>(module,float_approximation_class);
    define_lattice<FloatApproximation<PR>>(module,float_approximation_class);
    define_mixed_lattice<FloatApproximation<PR>,ApproximateNumber>(module,float_approximation_class);
    define_comparisons<FloatApproximation<PR>>(module,float_approximation_class);
    define_mixed_comparisons<FloatApproximation<PR>,ApproximateNumber>(module,float_approximation_class);

    float_approximation_class.def("__str__", &__cstr__<FloatApproximation<PR>>);
    float_approximation_class.def("__repr__", &__cstr__<FloatApproximation<PR>>);

    //NOTE: Conversion to FloatMP needs precision, so disallow conversion to FloatDP as well
    //implicitly_convertible<double,FloatApproximation<PR>>();
    implicitly_convertible<FloatValue<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatBall<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatBounds<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatUpperBound<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatLowerBound<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatApproximation<PR>,ApproximateNumber>();

    module.def("cast_exact",(FloatValue<PR>const&(*)(FloatApproximation<PR> const&)) &cast_exact);
    module.def("cast_exact",(FloatValue<PR>const&(*)(RawFloat<PR> const&)) &cast_exact);
}

Void export_effort(pymodule& module) {
    pybind11::class_<Effort> effort_class(module,"Effort");
    effort_class.def(init<Nat>());
    effort_class.def("work",&Effort::work);
    effort_class.def("__str__", &__cstr__<Effort>);
}

Void export_accuracy(pymodule& module) {
    pybind11::class_<Accuracy> accuracy_class(module,"Accuracy");
    accuracy_class.def(init<Nat>());
    accuracy_class.def("bits",&Accuracy::bits);
    accuracy_class.def("__str__", &__cstr__<Accuracy>);
    accuracy_class.def("__repr__", &__cstr__<Accuracy>);
}

template<class PR> Void export_precision(pymodule& module);

template<> Void export_precision<DoublePrecision>(pymodule& module) {
    pybind11::class_<DoublePrecision> precision_class(module,"DoublePrecision");
    precision_class.def(init<>());
    precision_class.def("__str__", &__cstr__<DoublePrecision>);
    module.attr("double_precision") = double_precision;
}

template<> Void export_precision<MultiplePrecision>(pymodule& module) {
    pybind11::class_<MultiplePrecision> precision_class(module,"MultiplePrecision");
    precision_class.def(init<Nat>());
    precision_class.def("bits",&MultiplePrecision::bits);
    precision_class.def("__str__", &__cstr__<MultiplePrecision>);
    module.def("multiple_precision", &multiple_precision);
    module.def("precision", &precision);

}

template<class PR> Void export_user_floats(pymodule& module) {
    export_float_approximation<PR>(module);
    export_float_upper_bound<PR>(module);
    export_float_lower_bound<PR>(module);
    export_float_bounds<PR>(module);
    export_float_ball<PR>(module);
    export_float_value<PR>(module);
    export_float_error<PR>(module);

    module.def("Value", [](Dyadic const& y, PR pr){return Value(y,pr);});
    module.def("Value", [](RawFloatType<PR> const& v){return Value(v);});

    if constexpr (IsSame<PR,MultiplePrecision>::value) {
        module.def("Ball", [](ValidatedNumber const& y, MultiplePrecision pr, DoublePrecision pre){return Ball(y,pr,pre);});
        module.def("Ball", [](RawFloatType<MultiplePrecision> const& v, RawFloatType<DoublePrecision> const& e){return Ball(v,e);});
    }

    module.def("Ball", [](ValidatedNumber const& y, PR pr){return Ball(y,pr);});
    module.def("Ball", [](ValidatedNumber const& y, PR pr, PR pre){return Ball(y,pr,pre);});
    module.def("Ball", [](RawFloatType<PR> const& v, RawFloatType<PR> const& e){return Ball(v,e);});
    module.def("Bounds", [](ValidatedNumber const& y, PR pr){return Bounds(y,pr);});
    module.def("Bounds", [](ValidatedLowerNumber const& yl, ValidatedUpperNumber const& yu, PR pr){return Bounds(yl,yu,pr);});
    module.def("Bounds", [](RawFloatType<PR> const& l, RawFloatType<PR> const& u){return Bounds(l,u);});
    module.def("UpperBound", [](ValidatedUpperNumber const& y, PR pr){return UpperBound(y,pr);});
    module.def("UpperBound", [](RawFloatType<PR> const& u){return UpperBound(u);});
    module.def("LowerBound", [](ValidatedLowerNumber const& y, PR pr){return LowerBound(y,pr);});
    module.def("LowerBound", [](RawFloatType<PR> const& l){return LowerBound(l);});
    module.def("Approximation", [](ApproximateNumber const& y, PR pr){return Approximation(y,pr);});
    module.def("Approximation", [](RawFloatType<PR> const& a){return Approximation(a);});
}


Void numeric_submodule(pymodule& module) {
    export_effort(module);
    export_accuracy(module);

    export_logicals(module);

    export_builtins(module);

    export_integer(module);
    export_dyadic(module);
    export_decimal(module);
    export_rational(module);
    export_real(module);

    export_dyadic_bounds(module);

    export_numbers(module);

    export_rounding_mode(module);
    export_precision<DoublePrecision>(module);
    export_precision<MultiplePrecision>(module);
    export_raw_float<DoublePrecision>(module);
    export_raw_float<MultiplePrecision>(module);

    export_user_floats<DoublePrecision>(module);
    export_user_floats<MultiplePrecision>(module);
    export_float_ball<MultiplePrecision,DoublePrecision>(module);
}

