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

template<class T> std::string to_cppstring(T const& t) { std::stringstream ss; ss<<t; return ss.str(); }


template<class T> struct PythonRepresentation {
    const T* pointer;
    PythonRepresentation(const T& t) : pointer(&t) { }
    const T& reference() const { return *pointer; }
};

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Rational>& repr) {
    return os << "Rational("<<repr.reference().numerator()<<","<<repr.reference().denominator()<<")"; }
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

Dyadic cast_exact(double d) { return Dyadic(ExactDouble(d)); }


template<class L> Bool decide(L l) { return Ariadne::decide(l); }
template<class L> Bool definitely(L l) { return Ariadne::definitely(l); }
template<class L> Bool possibly(L l) { return Ariadne::possibly(l); }

template<class L> decltype(auto) check(L const& l, Effort const& e) { return l.check(e); }
template<class L> decltype(auto) operator&(L const& l1, L const& l2) { return l1 and l2; }
template<class L> decltype(auto) operator|(L const& l1, L const& l2) { return l1 or l2; }
template<class L> decltype(auto) operator~(L const& l) { return not l; }


template<class OP, class... TS> auto py_apply(TS const& ... ts) -> decltype(OP()(ts...)){ OP op; return op(ts...); }

template<class... T> struct Tag { };

template<class OP> struct PythonOperator { };
PythonOperator<Pos> pos_op(pybind11::detail::self_t) { return PythonOperator<Pos>(); }
PythonOperator<Neg> neg_op(pybind11::detail::self_t) { return PythonOperator<Neg>(); }
PythonOperator<Sqr> sqr_op(pybind11::detail::self_t) { return PythonOperator<Sqr>(); }
PythonOperator<Hlf> hlf_op(pybind11::detail::self_t) { return PythonOperator<Hlf>(); }
PythonOperator<Rec> rec_op(pybind11::detail::self_t) { return PythonOperator<Rec>(); }
PythonOperator<Sqrt> sqrt_op(pybind11::detail::self_t) { return PythonOperator<Sqrt>(); }
PythonOperator<Exp> exp_op(pybind11::detail::self_t) { return PythonOperator<Exp>(); }
PythonOperator<Log> log_op(pybind11::detail::self_t) { return PythonOperator<Log>(); }
PythonOperator<Sin> sin_op(pybind11::detail::self_t) { return PythonOperator<Sin>(); }
PythonOperator<Cos> cos_op(pybind11::detail::self_t) { return PythonOperator<Cos>(); }
PythonOperator<Tan> tan_op(pybind11::detail::self_t) { return PythonOperator<Tan>(); }
PythonOperator<Atan> atan_op(pybind11::detail::self_t) { return PythonOperator<Atan>(); }

PythonOperator<Abs> abs_op(pybind11::detail::self_t) { return PythonOperator<Abs>(); }



template<class A> decltype(auto) _nul_(A const& a) { return nul(a); }
template<class A> decltype(auto) _pos_(A const& a) { return pos(a); }
template<class A> decltype(auto) _neg_(A const& a) { return neg(a); }
template<class A> decltype(auto) _hlf_(A const& a) { return hlf(a); }
template<class A> decltype(auto) _sqr_(A const& a) { return sqr(a); }
template<class A> decltype(auto) _rec_(A const& a) { return rec(a); }
template<class A1, class A2> decltype(auto) _add_(A1 const& a1, A2 const& a2) { return add(a1,a2); }
template<class A1, class A2> decltype(auto) _sub_(A1 const& a1, A2 const& a2) { return sub(a1,a2); }
template<class A1, class A2> decltype(auto) _mul_(A1 const& a1, A2 const& a2) { return mul(a1,a2); }
template<class A1, class A2> decltype(auto) _div_(A1 const& a1, A2 const& a2) { return div(a1,a2); }
template<class A1, class A2, class A3> decltype(auto) _fma_(A1 const& a1, A2 const& a2, A3 const& a3) { return fma(a1,a2,a3); }
template<class A, class N> decltype(auto) _pow_(A const& a, N n) { return pow(a,n); }
template<class A> decltype(auto) _sqrt_(A const& a) { return sqrt(a); }
template<class A> decltype(auto) _exp_(A const& a) { return exp(a); }
template<class A> decltype(auto) _log_(A const& a) { return log(a); }
template<class A> decltype(auto) _sin_(A const& a) { return sin(a); }
template<class A> decltype(auto) _cos_(A const& a) { return cos(a); }
template<class A> decltype(auto) _tan_(A const& a) { return tan(a); }
template<class A> decltype(auto) _atan_(A const& a) { return atan(a); }

template<class A> decltype(auto) _log2_(A const& a) { return log2(a); }

template<class A> decltype(auto) _abs_(A const& a) { return abs(a); }
template<class A1, class A2> decltype(auto) _max_(A1 const& a1, A2 const& a2) { return max(a1,a2); }
template<class A1, class A2> decltype(auto) _min_(A1 const& a1, A2 const& a2) { return min(a1,a2); }


template<class A> decltype(auto) _not_(A const& a) { return !a; }
template<class A1, class A2> decltype(auto) _or_(A1 const& a1, A2 const& a2) { return a1 || a2; }
template<class A1, class A2> decltype(auto) _and_(A1 const& a1, A2 const& a2) { return a1 && a2; }

template<class R, class A> decltype(auto) _nul_rnd_(R r, A const& a) { return nul(r,a); }
template<class R, class A> decltype(auto) _pos_rnd_(R r, A const& a) { return pos(r,a); }
template<class R, class A> decltype(auto) _neg_rnd_(R r, A const& a) { return neg(r,a); }
template<class R, class A> decltype(auto) _hlf_rnd_(R r, A const& a) { return hlf(r,a); }
template<class R, class A> decltype(auto) _sqr_rnd_(R r, A const& a) { return sqr(r,a); }
template<class R, class A> decltype(auto) _rec_rnd_(R r, A const& a) { return rec(r,a); }
template<class R, class A1, class A2> decltype(auto) _add_rnd_(R r, A1 const& a1, A2 const& a2) { return add(r,a1,a2); }
template<class R, class A1, class A2> decltype(auto) _sub_rnd_(R r, A1 const& a1, A2 const& a2) { return sub(r,a1,a2); }
template<class R, class A1, class A2> decltype(auto) _mul_rnd_(R r, A1 const& a1, A2 const& a2) { return mul(r,a1,a2); }
template<class R, class A1, class A2> decltype(auto) _div_rnd_(R r, A1 const& a1, A2 const& a2) { return div(r,a1,a2); }
template<class R, class A1, class A2, class A3> decltype(auto) _fma_rnd_(R r, A1 const& a1, A2 const& a2, A3 const& a3) { return fma(r,a1,a2,a3); }
template<class R, class A, class N> decltype(auto) _pow_rnd_(R r, A const& a, N n) { return pow(r,a,n); }
template<class R, class A> decltype(auto) _sqrt_rnd_(R r, A const& a) { return sqrt(r,a); }
template<class R, class A> decltype(auto) _exp_rnd_(R r, A const& a) { return exp(r,a); }
template<class R, class A> decltype(auto) _log_rnd_(R r, A const& a) { return log(r,a); }
template<class R, class A> decltype(auto) _sin_rnd_(R r, A const& a) { return sin(r,a); }
template<class R, class A> decltype(auto) _cos_rnd_(R r, A const& a) { return cos(r,a); }
template<class R, class A> decltype(auto) _tan_rnd_(R r, A const& a) { return tan(r,a); }
template<class R, class A> decltype(auto) _atan_rnd_(R r, A const& a) { return atan(r,a); }


template<class T> class numeric_class_ : public pybind11::class_<T> {
    typedef T const& Tcr;
  public:
    numeric_class_(pybind11::module const& module, String const& name)
        : pybind11::class_<T>(module,name.c_str()) { }
    using pybind11::class_<T>::def;
    static constexpr auto self = pybind11::detail::self;
    template<class OP> void def(PythonOperator<OP>) { this->pybind11::class_<T>::def(to_str(OP()).c_str(),&py_apply<OP,T>); }
    void define_self_arithmetic() {
        using pybind11::self;
        T const* other_ptr=nullptr; T const& other=*other_ptr;
        this->def(+self); this->def(-self);
        this->def(self+other); this->def(self-other); this->def(self*other); this->def(self/other);
        this->def(other+self); this->def(other-self); this->def(other*self); this->def(other/self);
    }
    template<class X> void define_mixed_arithmetic() {
        X const* other_ptr=nullptr; X const& other=*other_ptr;
        this->def(self+other); this->def(self-other); this->def(self*other); this->def(self/other);
        this->def(other+self); this->def(other-self); this->def(other*self); this->def(other/self);
    }
    void define_unary_arithmetic() {
        this->def(+self); this->def(-self);
        this->def(pos_op(self)); this->def(neg_op(self));
        this->def(sqr_op(self)); this->def(rec_op(self));
    }
    void define_transcendental_functions() {
        this->define_unary_arithmetic();
        this->def(sqrt_op(self)); this->def(exp_op(self)); this->def(log_op(self));
        this->def(sin_op(self)); this->def(cos_op(self)); this->def(tan_op(self)); this->def(atan_op(self));
    }
    void define_monotonic_functions() {
        this->def(pos_op(self)); this->def(neg_op(self));
        this->def(sqrt_op(self)); this->def(exp_op(self)); this->def(log_op(self)); this->def(atan_op(self));
    }
    void define_self_comparisons() {
        this->def(self==self); this->def(self!=self);
        this->def(self<=self); this->def(self>=self); this->def(self<self); this->def(self>self);
    }
    void define_self_logical() {
        //T const* self_ptr=nullptr; T const& self=*self_ptr;
//        //using pybind11::self_ns::self;
        this->def("__and__", (T(*)(T const&, T const&)) &operator&);
        //self & self); this->def(self | self); this->def(~self);
    }
    void convert_from() { }
    template<class A, class... AS> void convert_from() {
        this->def(pybind11::init<A>()); pybind11::implicitly_convertible<A,T>; this->convert_from<AS...>();
    }
    void inits() { }
    template<class A, class... AS> void inits() {
        this->def(pybind11::init<A>()); this->inits<AS...>();
    }
};

template<class T> class logical_class_ : public pybind11::class_<T> {
    typedef T const& Tcr;
  public:
    logical_class_(pybind11::module const& module, String const& name)
        : pybind11::class_<T>(module,name.c_str()) { }
    using pybind11::class_<T>::def;
    void define_self_logical() {
        typedef decltype(! declval<T>()) NT;
        this->def("__and__", (T(*)(T const&, T const&)) _and_);
        this->def("__or__", (T(*)(T const&, T const&)) _or_);
        this->def("__not__", (NT(*)(T const&)) _not_);
    }
};


} // namespace Ariadne 


using namespace Ariadne;
using pymodule = pybind11::module;
using pybind11::init;
using pybind11::detail::self;
using pybind11::implicitly_convertible;

template<class P> void export_effective_logical(pymodule& module, std::string name)
{
    typedef decltype(declval<LogicalType<P>>().check(declval<Effort>())) CheckType;
    OutputStream& operator<<(OutputStream& os, LogicalType<P> l);

    logical_class_<LogicalType<P>> logical_class(module,name);
    logical_class.def(init<bool>());
    logical_class.def(init<LogicalType<P>>());
    logical_class.def("check", (CheckType(LogicalType<P>::*)(Effort)) &LogicalType<P>::check);
    logical_class.define_self_logical();
//    logical_class.def(pybind11::detail::str_op(self));
//    logical_class.def(pybind11::detail::repr_op(self));
    module.def("check", (CheckType(*)(LogicalType<P> const&,Effort const&)) &Ariadne::check<LogicalType<P>>);
}

template<class P> void export_logical(pymodule& module, std::string name)
{
    OutputStream& operator<<(OutputStream& os, LogicalType<P> l);
    logical_class_<LogicalType<P>> logical_class(module,name);
    logical_class.def(init<bool>());
    logical_class.def(init<LogicalType<P>>());
    logical_class.define_self_logical();
    logical_class.def("__str__", &to_cppstring<LogicalType<P>>);
    logical_class.def("__repr__", &to_cppstring<LogicalType<P>>);

    module.def("decide", (bool(*)(LogicalType<P>)) &decide);
    module.def("possibly", (bool(*)(LogicalType<P>)) &possibly);
    module.def("definitely", (bool(*)(LogicalType<P>)) &definitely);

}

template<> void export_logical<ExactTag>(pymodule& module, std::string name)
{
    typedef ExactTag P;
    OutputStream& operator<<(OutputStream& os, LogicalType<P> l);
    logical_class_<LogicalType<P>> logical_class(module,name);
    logical_class.def(init<bool>());
    logical_class.def(init<LogicalType<P>>());
    logical_class.define_self_logical();
    logical_class.def("__str__", &to_cppstring<LogicalType<P>>);
    logical_class.def("__repr__", &to_cppstring<LogicalType<P>>);

//    implicitly_convertible<LogicalType<ExactTag>,bool>();
}



Void export_logicals(pymodule& module) {
    export_logical<ExactTag>(module,"Boolean");
    export_effective_logical<EffectiveTag>(module,"Kleenean");
    export_effective_logical<EffectiveUpperTag>(module,"Sierpinskian");
    export_effective_logical<EffectiveLowerTag>(module,"NegatedSierpinskian");
    export_logical<ValidatedTag>(module,"Tribool");
    export_logical<UpperTag>(module,"Verified");
    export_logical<LowerTag>(module,"Falsified");
    export_logical<ApproximateTag>(module,"Fuzzy");
}


void export_integer(pymodule& module)
{
    numeric_class_<Integer> integer_class(module,"Integer");
    integer_class.def(init<int>());
    integer_class.def(init<Integer>());

    integer_class.def("__str__", &to_cppstring<Integer>);
    integer_class.def("__repr__", &to_cppstring<Integer>);

    // Allows mixed operations with integral types
    integer_class.define_self_arithmetic();
    integer_class.define_self_comparisons();

    module.def("pow", (Integer(*)(const Integer&,Nat)) &_pow_) ;
    integer_class.def(abs_op(self));


    implicitly_convertible<int,Integer>();

}

void export_dyadic(pymodule& module)
{

    numeric_class_<Dyadic> dyadic_class(module,"Dyadic");
    dyadic_class.def(init<Int,Nat>());
//    dyadic_class.def(init<Integer,Natural>());
//    dyadic_class.def(init<Int>());
    dyadic_class.def(init<Integer>());
    dyadic_class.def(init<Dyadic>());
    dyadic_class.def(init<FloatDPValue>());
    dyadic_class.def(init<FloatMPValue>());

    dyadic_class.define_self_arithmetic();
    dyadic_class.define_self_comparisons();
    module.def("hlf", (Dyadic(*)(Dyadic const&)) &_hlf_);

    dyadic_class.def("__str__", &to_cppstring<Dyadic>);
    dyadic_class.def("__repr__", &to_cppstring<Dyadic>);

    dyadic_class.def("get_d", &Dyadic::get_d);

    implicitly_convertible<int,Dyadic>();
    implicitly_convertible<Integer,Dyadic>();

    module.def("cast_exact",(Dyadic(*)(double)) &cast_exact);
}

void export_rational(pymodule& module)
{
    numeric_class_<Rational> rational_class(module,"Rational");
    rational_class.def(init<int,int>());
    rational_class.def(init<Integer,Integer>());
    rational_class.def(init<int>());
    //rational_class.def(init<double>());
    rational_class.def(init<Integer>());
    rational_class.def(init<Dyadic>());
    rational_class.def(init<Decimal>());
    rational_class.def(init<Rational>());

    rational_class.define_self_arithmetic();
    rational_class.define_self_comparisons();

    rational_class.def("__str__", &to_cppstring<Rational>);
    rational_class.def("__repr__", &to_cppstring<Rational>);

    rational_class.def("get_d", &Rational::get_d);

    implicitly_convertible<int,Rational>();
    implicitly_convertible<Integer,Rational>();
    implicitly_convertible<Dyadic,Rational>();
}

void export_decimal(pymodule& module)
{
    numeric_class_<Decimal> decimal_class(module,"Decimal");
    decimal_class.def(init<Decimal>());
    decimal_class.def(init<StringType>());
    decimal_class.def(init<double>());
    decimal_class.def("__str__", &to_cppstring<Decimal>);
    decimal_class.def("__repr__", &to_cppstring<Decimal>);
    implicitly_convertible<Decimal,Rational>();
}

void export_real(pymodule& module)
{
    numeric_class_<Real> real_class(module,"Real");
    real_class.def(init<int>());
    real_class.def(init<Integer>());
    real_class.def(init<Dyadic>());
    real_class.def(init<Rational>());
    real_class.def(init<Real>());

    real_class.define_self_arithmetic();
    real_class.define_transcendental_functions();
    real_class.define_self_comparisons();

    module.def("sqrt", (Real(*)(Real const&)) &_sqrt_);
    module.def("exp", (Real(*)(Real const&)) &_exp_);
    module.def("log", (Real(*)(Real const&)) &_log_);
    module.def("sin", (Real(*)(Real const&)) &_sin_);
    module.def("cos", (Real(*)(Real const&)) &_cos_);
    module.def("tan", (Real(*)(Real const&)) &_tan_);

    real_class.def("__str__", &to_cppstring<Real>);
    real_class.def("__repr__", &to_cppstring<Real>);

    real_class.def("get", (FloatDPBounds(Real::*)(DoublePrecision)const) &Real::get);
    real_class.def("get", (FloatMPBounds(Real::*)(MultiplePrecision)const) &Real::get);
    real_class.def("compute", (ValidatedReal(Real::*)(Effort)const) &Real::compute);
    real_class.def("compute", (ValidatedReal(Real::*)(Accuracy)const) &Real::compute);
    real_class.def("get_d", &Real::get_d);

    implicitly_convertible<Rational,Real>();


    numeric_class_<ValidatedReal> validated_real_class(module,"ValidatedReal");
    validated_real_class.def(init<DyadicBounds>());
    validated_real_class.def("get", (DyadicBounds(ValidatedReal::*)()const) &ValidatedReal::get);
    validated_real_class.def("get", (FloatDPBounds(ValidatedReal::*)(DoublePrecision)const) &ValidatedReal::get);
    validated_real_class.def("get", (FloatMPBounds(ValidatedReal::*)(MultiplePrecision)const) &ValidatedReal::get);

    validated_real_class.def("__str__", &to_cppstring<ValidatedReal>);
    validated_real_class.def("__repr__", &to_cppstring<ValidatedReal>);
}


template<class P> void export_number()
{
    numeric_class_<Number<P>> number_class(class_name<P>()+"Number");
    number_class.def(init<Rational>());
    number_class.define_self_arithmetic();
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
    numeric_class_<ApproximateNumber> approximate_number_class(module,class_name<ApproximateTag>()+"Number");
    approximate_number_class.def(init<Rational>());
    approximate_number_class.def(init<Real>());
    approximate_number_class.def(init<ExactNumber>());
    approximate_number_class.def(init<EffectiveNumber>());
    approximate_number_class.def(init<ValidatedNumber>());
    approximate_number_class.def(init<ApproximateNumber>());
    approximate_number_class.def("get", (FloatDPApproximation(*)(ApproximateNumber const&, DoublePrecision const&)) &get);
    approximate_number_class.def("get", (FloatMPApproximation(*)(ApproximateNumber const&, MultiplePrecision const&)) &get);
    approximate_number_class.define_self_arithmetic();
    approximate_number_class.define_mixed_arithmetic<ApproximateNumber>();
    approximate_number_class.define_transcendental_functions();
    approximate_number_class.def("__str__", &to_cppstring<ApproximateNumber>);
    approximate_number_class.def("__repr__", &to_cppstring<ApproximateNumber>);

    numeric_class_<ValidatedLowerNumber> lower_number_class(module,class_name<LowerTag>()+"Number");
    lower_number_class.def(init<ValidatedNumber>());
    lower_number_class.def("get", (FloatDPLowerBound(*)(ValidatedLowerNumber const&, DoublePrecision const&)) &get);
    lower_number_class.def("get", (FloatMPLowerBound(*)(ValidatedLowerNumber const&, MultiplePrecision const&)) &get);
    lower_number_class.define_monotonic_functions();
    lower_number_class.def("__str__", &to_cppstring<ValidatedLowerNumber>);
    lower_number_class.def("__repr__", &to_cppstring<ValidatedLowerNumber>);

    numeric_class_<ValidatedUpperNumber> upper_number_class(module,class_name<UpperTag>()+"Number");
    upper_number_class.def(init<ValidatedNumber>());
    upper_number_class.def("get", (FloatDPUpperBound(*)(ValidatedUpperNumber const&, DoublePrecision const&)) &get);
    upper_number_class.def("get", (FloatMPUpperBound(*)(ValidatedUpperNumber const&, MultiplePrecision const&)) &get);
    upper_number_class.define_monotonic_functions();
    upper_number_class.def("__str__", &to_cppstring<ValidatedUpperNumber>);
    upper_number_class.def("__repr__", &to_cppstring<ValidatedUpperNumber>);

    numeric_class_<ValidatedNumber> validated_number_class(module,class_name<ValidatedTag>()+"Number");
    validated_number_class.def(init<Rational>());
    validated_number_class.def(init<Real>());
    validated_number_class.def(init<ExactNumber>());
    validated_number_class.def(init<EffectiveNumber>());
    validated_number_class.def(init<ValidatedNumber>());
    validated_number_class.def("get", (FloatDPBounds(*)(ValidatedNumber const&, DoublePrecision const&)) &get);
    validated_number_class.def("get", (FloatMPBounds(*)(ValidatedNumber const&, MultiplePrecision const&)) &get);
    validated_number_class.define_self_arithmetic();
    validated_number_class.define_mixed_arithmetic<ValidatedNumber>();
    validated_number_class.define_transcendental_functions();
    validated_number_class.def("__str__", &to_cppstring<ValidatedNumber>);
    validated_number_class.def("__repr__", &to_cppstring<ValidatedNumber>);

    numeric_class_<EffectiveNumber> effective_number_class(module,class_name<EffectiveTag>()+"Number");
    effective_number_class.def(init<Rational>());
    effective_number_class.def(init<Real>());
    effective_number_class.def(init<ExactNumber>());
    effective_number_class.def(init<EffectiveNumber>());
    effective_number_class.def("get", (FloatDPBounds(*)(EffectiveNumber const&, DoublePrecision const&)) &get);
    effective_number_class.def("get", (FloatMPBounds(*)(EffectiveNumber const&, MultiplePrecision const&)) &get);
    effective_number_class.define_self_arithmetic();
    effective_number_class.define_mixed_arithmetic<EffectiveNumber>();
    effective_number_class.define_transcendental_functions();
    effective_number_class.def("__str__", &to_cppstring<EffectiveNumber>);
    effective_number_class.def("__repr__", &to_cppstring<EffectiveNumber>);

    numeric_class_<ExactNumber> exact_number_class(module,class_name<ExactTag>()+"Number");
    exact_number_class.def(init<Rational>());
    exact_number_class.def(init<ExactNumber>());
    exact_number_class.def("get", (FloatDPBounds(*)(ExactNumber const&, DoublePrecision const&)) &get);
    exact_number_class.def("get", (FloatMPBounds(*)(ExactNumber const&, MultiplePrecision const&)) &get);
    exact_number_class.def("__str__", &to_cppstring<ExactNumber>);
    exact_number_class.def("__repr__", &to_cppstring<ExactNumber>);


    implicitly_convertible<ValidatedNumber,ApproximateNumber>();
    implicitly_convertible<ValidatedNumber,ValidatedLowerNumber>();
    implicitly_convertible<ValidatedNumber,ValidatedUpperNumber>();
    implicitly_convertible<EffectiveNumber,ValidatedNumber>();
    implicitly_convertible<ExactNumber,EffectiveNumber>();

    implicitly_convertible<Rational,ExactNumber>();
    implicitly_convertible<Real,EffectiveNumber>();
}




void export_dyadic_bounds(pymodule& module)
{
    numeric_class_<DyadicBounds> dyadic_bounds_class(module,"DyadicBounds");
    dyadic_bounds_class.def(init<DyadicBounds>());
    dyadic_bounds_class.def(init<Dyadic,Dyadic>());

    dyadic_bounds_class.def("__str__", &to_cppstring<DyadicBounds>);
    dyadic_bounds_class.def("__repr__", &to_cppstring<DyadicBounds>);

}

void export_rounding_mode(pymodule& module) {
    numeric_class_<Rounding> rounding_mode_class(module,"Rounding");
    rounding_mode_class.def(init<Rounding>());
//    boost::python::scope().attr("up") = Rounding(up);
//    boost::python::scope().attr("down") = Rounding(down);
//    boost::python::scope().attr("near") = Rounding(near);
}

template<class PR> void export_raw_float(pymodule& module)
{
    typedef RawFloat<PR> F;
    typedef F const& Fcr;
    typedef typename F::RoundingModeType RND;
//    implicitly_convertible<Rounding,RND>();

    FloatMP const& arg_type(FloatMP);
    FloatDP arg_type(FloatDP);

    numeric_class_<F> raw_float_class(module,"Float"+numeric_class_tag<PR>());
    raw_float_class.def(init<double,PR>());
    raw_float_class.def(init<Dyadic,PR>());
    raw_float_class.def(init<Rational,RND,PR>());
    raw_float_class.def("__str__", &to_cppstring<RawFloat<PR>>);
    raw_float_class.def("__repr__", &to_cppstring<RawFloat<PR>>);

    module.def("nul", (F(*)(Fcr)) &_nul_);
    module.def("pos", (F(*)(Fcr)) &_pos_);
    module.def("neg", (F(*)(Fcr)) &_neg_);
    module.def("hlf", (F(*)(Fcr)) &_hlf_);

    module.def("nul", (F(*)(RND,Fcr)) &_nul_rnd_);
    module.def("pos", (F(*)(RND,Fcr)) &_pos_rnd_);
    module.def("neg", (F(*)(RND,Fcr)) _neg_rnd_);
    module.def("hlf", (F(*)(RND,Fcr)) _hlf_rnd_);
    module.def("sqr", (F(*)(RND,Fcr)) _sqr_rnd_);
    module.def("rec", (F(*)(RND,Fcr)) _rec_rnd_);
    module.def("add", (F(*)(RND,Fcr,Fcr)) _add_rnd_);
    module.def("sub", (F(*)(RND,Fcr,Fcr)) _sub_rnd_);
    module.def("mul", (F(*)(RND,Fcr,Fcr)) _mul_rnd_);
    module.def("div", (F(*)(RND,Fcr,Fcr)) _div_rnd_);
    module.def("fma", (F(*)(RND,Fcr,Fcr,Fcr)) _fma_rnd_);
    module.def("pow", (F(*)(RND,Fcr,Int)) _pow_rnd_);
    module.def("sqrt", (F(*)(RND,Fcr)) _sqrt_rnd_);
    module.def("exp", (F(*)(RND,Fcr)) _exp_rnd_);
    module.def("log", (F(*)(RND,Fcr)) _log_rnd_);
    module.def("sin", (F(*)(RND,Fcr)) _sin_rnd_);
    module.def("cos", (F(*)(RND,Fcr)) _cos_rnd_);
    module.def("tan", (F(*)(RND,Fcr)) _tan_rnd_);
    module.def("atan", (F(*)(RND,Fcr)) _atan_rnd_);

    module.def("abs", (F(*)(Fcr)) _abs_);
    module.def("max", (F(*)(Fcr,Fcr)) _min_);
    module.def("min", (F(*)(Fcr,Fcr)) _max_);

}

template<class PR> void export_float_value(pymodule& module)
{
    numeric_class_<FloatValue<PR>> float_value_class(module,"Float"+numeric_class_tag<PR>()+"Value");
    float_value_class.def(init<RawFloat<PR>>());
    float_value_class.def(init<int>());
    float_value_class.def(init<double>());
    float_value_class.def(init<Integer,PR>());
    float_value_class.def(init<Dyadic,PR>());
    float_value_class.def(init<FloatValue<PR>>());

    float_value_class.def(pos_op(self));
    float_value_class.def(neg_op(self));

    float_value_class.def("precision", &FloatValue<PR>::precision);
    float_value_class.def("get_d",&FloatValue<PR>::get_d);

    float_value_class.template define_mixed_arithmetic<Int>();
    float_value_class.define_self_arithmetic();
    //float_value_class.define_mixed_arithmetic<ApproximateNumericType>();
    //float_value_class.define_mixed_arithmetic<LowerNumericType>();
    //float_value_class.define_mixed_arithmetic<UpperNumericType>();
    //float_value_class.define_mixed_arithmetic<ValidatedNumericType>();
    float_value_class.define_self_comparisons();

//    float_value_class.def(pos_op(self));
//    float_value_class.def(neg_op(self));

    float_value_class.def("precision", &FloatValue<PR>::precision);
    float_value_class.def("raw", (RawFloat<PR>const&(FloatValue<PR>::*)()const)&FloatValue<PR>::raw);
    float_value_class.def("get_d",&FloatValue<PR>::get_d);

    float_value_class.def("__str__", &to_cppstring<FloatValue<PR>>);
    float_value_class.def("__repr__", &to_cppstring<FloatValue<PR>>);

    float_value_class.def_static("set_output_places",&FloatValue<PR>::set_output_places);
}

template<class PRE> void export_float_error(pymodule& module)
{
    numeric_class_<FloatError<PRE>> float_error_class(module,"Float"+numeric_class_tag<PRE>()+"Error");
    float_error_class.def(init<RawFloat<PRE>>());
    float_error_class.def(init<uint>());
    float_error_class.def(init<double>());
    float_error_class.def(init<FloatError<PRE>>());

    float_error_class.def(+self);
    float_error_class.def(self+self);
    float_error_class.def(self*self);

    //float_error_class.def("raw",(RawFloat<PRE>&(FloatError<PRE>::*)())&FloatError<PRE>::raw,return_value_policy<Reference_existing_object>());
    float_error_class.def("raw",(RawFloat<PRE>const&(FloatError<PRE>::*)()const)&FloatError<PRE>::raw);
    float_error_class.def("__str__", &to_cppstring<FloatError<PRE>>);
    float_error_class.def("__repr__", &to_cppstring<FloatError<PRE>>);

    float_error_class.def("precision", &FloatError<PRE>::precision);

    module.def("log2", (FloatUpperBound<PRE>(*)(FloatError<PRE>const&)) &_log2_);
    float_error_class.def_static("set_output_places",&FloatError<PRE>::set_output_places);
}

template<class PR, class PRE=PR> void export_float_ball(pymodule& module)
{
    String tag = IsSame<PR,PRE>::value ? numeric_class_tag<PR>() : numeric_class_tag<PR>()+numeric_class_tag<PRE>();
    numeric_class_<FloatBall<PR,PRE>> float_ball_class(module,"Float"+tag+"Ball");
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

    float_ball_class.template define_mixed_arithmetic<ValidatedNumber>();
    float_ball_class.template define_mixed_arithmetic<Int>();
    float_ball_class.define_self_arithmetic();
//    float_ball_class.define_mixed_arithmetic<FloatBall<PR>>();
//    float_ball_class.define_mixed_arithmetic<ApproximateNumericType>();
//    float_ball_class.define_mixed_arithmetic<LowerNumericType>();
//    float_ball_class.define_mixed_arithmetic<UpperNumericType>();
//    float_ball_class.define_mixed_arithmetic<ValidatedNumericType>();
    float_ball_class.define_self_comparisons();

    float_ball_class.define_unary_arithmetic();
    float_ball_class.define_transcendental_functions();

    float_ball_class.def("__str__", &to_cppstring<FloatBall<PR>>);
    float_ball_class.def("__repr__", &to_cppstring<FloatBall<PR>>);
}



template<class PR> void export_float_bounds(pymodule& module)
{
    numeric_class_<FloatBounds<PR>> float_bounds_class(module,"Float"+numeric_class_tag<PR>()+"Bounds");
    float_bounds_class.def(init<RawFloat<PR>,RawFloat<PR>>());
    float_bounds_class.def(init<double,double>());
    float_bounds_class.def(init<FloatLowerBound<PR>,FloatUpperBound<PR>>());
    float_bounds_class.def(init<Real,PR>());

    float_bounds_class.def(init<ExactDouble,PR>());
    float_bounds_class.def(init<ValidatedNumber,PR>());
    float_bounds_class.def(init<FloatValue<PR>>());
    float_bounds_class.def(init<FloatBall<PR>>());
    float_bounds_class.def(init<FloatBounds<PR>>());

    float_bounds_class.def("lower", &FloatBounds<PR>::lower);
    float_bounds_class.def("upper", &FloatBounds<PR>::upper);
    float_bounds_class.def("value", &FloatBounds<PR>::value);
    float_bounds_class.def("error", &FloatBounds<PR>::error);

    float_bounds_class.def("precision", &FloatBounds<PR>::precision);

    float_bounds_class.template define_mixed_arithmetic<ValidatedNumber>();
    float_bounds_class.template define_mixed_arithmetic<Int>();
    float_bounds_class.template define_mixed_arithmetic<FloatBall<PR>>();
    float_bounds_class.define_self_arithmetic();
//    float_bounds_class.define_mixed_arithmetic<ApproximateNumericType>();
//    float_bounds_class.define_mixed_arithmetic<LowerNumericType>();
//    float_bounds_class.define_mixed_arithmetic<UpperNumericType>();
//    float_bounds_class.define_mixed_arithmetic<ValidatedNumericType>();
    float_bounds_class.define_self_comparisons();

    float_bounds_class.define_transcendental_functions();

    float_bounds_class.def("__str__", &to_cppstring<FloatBounds<PR>>);
    float_bounds_class.def("__repr__", &to_cppstring<FloatBounds<PR>>);

    implicitly_convertible<FloatValue<PR>,FloatBounds<PR>>();
    implicitly_convertible<FloatBall<PR>,FloatBounds<PR>>();

    float_bounds_class.def_static("set_output_places",&FloatBounds<PR>::set_output_places);
}

template<class PR> void export_float_upper_bound(pymodule& module)
{
    numeric_class_<FloatUpperBound<PR>> float_upper_bound_class(module,"Float"+numeric_class_tag<PR>()+"UpperBound");
    float_upper_bound_class.def(init<RawFloat<PR>>());
    float_upper_bound_class.def(init<int>());
    float_upper_bound_class.def(init<double>());
    float_upper_bound_class.def(init<Real,PR>());

    float_upper_bound_class.def(init<FloatValue<PR>>());
    float_upper_bound_class.def(init<FloatBall<PR>>());
    float_upper_bound_class.def(init<FloatBounds<PR>>());
    float_upper_bound_class.def(init<FloatUpperBound<PR>>());
    float_upper_bound_class.def(init<ValidatedUpperNumber,PR>());
    float_upper_bound_class.def("__str__", &to_cppstring<FloatUpperBound<PR>>);
    float_upper_bound_class.def("__repr__", &to_cppstring<FloatUpperBound<PR>>);

    float_upper_bound_class.def("precision", &FloatUpperBound<PR>::precision);
    float_upper_bound_class.def("raw", (RawFloat<PR>const&(FloatUpperBound<PR>::*)()const)&FloatUpperBound<PR>::raw);

    float_upper_bound_class.def(+self);
    float_upper_bound_class.def(-self);
    float_upper_bound_class.def(self + self);
    float_upper_bound_class.def(self + FloatUpperBound<PR>());
    float_upper_bound_class.def(FloatUpperBound<PR>() + self);
    float_upper_bound_class.def(self - self);
    float_upper_bound_class.def(self - FloatLowerBound<PR>());
    float_upper_bound_class.def(FloatLowerBound<PR>() - self);

    float_upper_bound_class.def(self + Int());
    float_upper_bound_class.def(Int() + self);
    float_upper_bound_class.def(self - Int());
    float_upper_bound_class.def(Int() - self);

    float_upper_bound_class.def(self > FloatBounds<PR>());
    float_upper_bound_class.def(self > FloatApproximation<PR>());
    float_upper_bound_class.def(self < FloatLowerBound<PR>());
    float_upper_bound_class.def(self >= FloatBounds<PR>());
    float_upper_bound_class.def(self >= FloatApproximation<PR>());
    float_upper_bound_class.def(self <= FloatLowerBound<PR>());

//    float_upper_bound_class.define_mixed_arithmetic<ApproximateNumericType>();
//    float_upper_bound_class.def(UpperNumericType() + self);
//    float_upper_bound_class.def(LowerNumericType() - self);
//    float_upper_bound_class.def(self + UpperNumericType());
//    float_upper_bound_class.def(self - LowerNumericType());

    float_upper_bound_class.define_monotonic_functions();

    implicitly_convertible<FloatBounds<PR>,FloatUpperBound<PR>>();
}

template<class PR> void export_float_lower_bound(pymodule& module)
{
    numeric_class_<FloatLowerBound<PR>> float_lower_bound_class(module,"Float"+numeric_class_tag<PR>()+"LowerBound");
    float_lower_bound_class.def(init<RawFloat<PR>>());
    float_lower_bound_class.def(init<int>());
    float_lower_bound_class.def(init<double>());
    float_lower_bound_class.def(init<Real,PR>());

    float_lower_bound_class.def(init<FloatValue<PR>>());
    float_lower_bound_class.def(init<FloatBall<PR>>());
    float_lower_bound_class.def(init<FloatBounds<PR>>());
    float_lower_bound_class.def(init<FloatLowerBound<PR>>());
    float_lower_bound_class.def(init<ValidatedLowerNumber,PR>());
    float_lower_bound_class.def("__str__", &to_cppstring<FloatLowerBound<PR>>);
    float_lower_bound_class.def("__repr__", &to_cppstring<FloatLowerBound<PR>>);

    float_lower_bound_class.def("precision", &FloatLowerBound<PR>::precision);

    float_lower_bound_class.def(+self);
    float_lower_bound_class.def(-self);
    float_lower_bound_class.def(self + self);
    float_lower_bound_class.def(self + FloatLowerBound<PR>());
    float_lower_bound_class.def(FloatLowerBound<PR>() + self);
    float_lower_bound_class.def(self - self);
    float_lower_bound_class.def(self - FloatUpperBound<PR>());
    float_lower_bound_class.def(FloatUpperBound<PR>() - self);

    float_lower_bound_class.def(self + Int());
    float_lower_bound_class.def(Int() + self);
    float_lower_bound_class.def(self - Int());
    float_lower_bound_class.def(Int() - self);

    float_lower_bound_class.def(self < FloatBounds<PR>());
    float_lower_bound_class.def(self < FloatApproximation<PR>());
    float_lower_bound_class.def(self > FloatUpperBound<PR>());
    float_lower_bound_class.def(self <= FloatBounds<PR>());
    float_lower_bound_class.def(self <= FloatApproximation<PR>());
    float_lower_bound_class.def(self >= FloatUpperBound<PR>());
    float_lower_bound_class.def("raw", (RawFloat<PR>const&(FloatLowerBound<PR>::*)()const)&FloatLowerBound<PR>::raw);

    //    float_lower_bound_class.define_mixed_arithmetic<ApproximateNumericType>();

//    float_lower_bound_class.def(LowerNumericType() + self);
//    float_lower_bound_class.def(UpperNumericType() - self);
//    float_lower_bound_class.def(self + LowerNumericType());
//    float_lower_bound_class.def(self - UpperNumericType());

    float_lower_bound_class.define_monotonic_functions();

    implicitly_convertible<FloatBounds<PR>,FloatLowerBound<PR>>();
}



template<class PR> void export_float_approximation(pymodule& module)
{
    numeric_class_<FloatApproximation<PR>> float_approximation_class(module,"Float"+numeric_class_tag<PR>()+"Approximation");

    if(IsSame<PR,DoublePrecision>::value) {
        float_approximation_class.def(init<double>());
    }
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

    float_approximation_class.template define_mixed_arithmetic<double>();
    float_approximation_class.template define_mixed_arithmetic<ApproximateNumber>();
    float_approximation_class.define_self_arithmetic();
    float_approximation_class.define_self_comparisons();

    float_approximation_class.def("precision", &FloatApproximation<PR>::precision);
    float_approximation_class.def("raw", (RawFloat<PR>const&(FloatApproximation<PR>::*)()const)&FloatApproximation<PR>::raw);
    float_approximation_class.def("get_d", &FloatApproximation<PR>::get_d);

    float_approximation_class.def("__str__", &to_cppstring<FloatApproximation<PR>>);
    float_approximation_class.def("__repr__", &to_cppstring<FloatApproximation<PR>>);

    implicitly_convertible<double,FloatApproximation<PR>>();
    //implicitly_convertible<Integer,FloatBounds<PR>>();
    implicitly_convertible<FloatValue<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatBall<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatBounds<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatUpperBound<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatLowerBound<PR>,FloatApproximation<PR>>();

    float_approximation_class.def(abs_op(self));

    float_approximation_class.define_unary_arithmetic();
    float_approximation_class.define_transcendental_functions();
}

Void export_effort(pymodule& module) {
    numeric_class_<Effort> effort_class(module,"Effort");
    effort_class.def(init<Nat>());
    effort_class.def("work",&Effort::work);
    effort_class.def("__str__", &to_cppstring<Effort>);
}

Void export_accuracy(pymodule& module) {
    numeric_class_<Accuracy> accuracy_class(module,"Accuracy");
    accuracy_class.def(init<Nat>());
    accuracy_class.def("bits",&Accuracy::bits);
    accuracy_class.def("__str__", &to_cppstring<Accuracy>);
    accuracy_class.def("__repr__", &to_cppstring<Accuracy>);
}

template<class PR> Void export_precision(pymodule& module);

template<> Void export_precision<DoublePrecision>(pymodule& module) {
    numeric_class_<DoublePrecision> precision_class(module,"DoublePrecision");
    precision_class.def(init<>());
    precision_class.def("__str__", &to_cppstring<DoublePrecision>);
//    boost::python::scope().attr("double_precision") = double_precision;
}

template<> Void export_precision<MultiplePrecision>(pymodule& module) {
    numeric_class_<MultiplePrecision> precision_class(module,"MultiplePrecision");
    precision_class.def(init<Nat>());
    precision_class.def("bits",&MultiplePrecision::bits);
    precision_class.def("__str__", &to_cppstring<MultiplePrecision>);
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
}


Void numeric_submodule(pymodule& module) {
    export_effort(module);
    export_accuracy(module);

    export_logicals(module);

    export_rounding_mode(module);
    export_precision<DoublePrecision>(module);
    export_precision<MultiplePrecision>(module);
    export_raw_float<DoublePrecision>(module);
    export_raw_float<MultiplePrecision>(module);

    export_user_floats<DoublePrecision>(module);
    export_user_floats<MultiplePrecision>(module);
    export_float_ball<MultiplePrecision,DoublePrecision>(module);

    export_dyadic_bounds(module);

    export_numbers(module);
    export_real(module);
    export_rational(module);
    export_decimal(module);
    export_dyadic(module);
    export_integer(module);
}

