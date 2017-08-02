/***************************************************************************
 *            numeric_submodule.cpp
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


#include "boost_python.hpp"
#include "utilities.hpp"

#include "numeric/logical.hpp"
#include "numeric/integer.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "numeric/real.hpp"
#include "numeric/float64.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/float-user.hpp"

#define DECLARE_NUMERIC_OPERATIONS(Xcr,X,PX) \
    X add(Xcr,Xcr); X sub(Xcr,Xcr); X mul(Xcr,Xcr); X div(Xcr,Xcr); \
    X pos(Xcr); X neg(Xcr); X sqr(Xcr); X rec(Xcr); X pow(Xcr,Int); \
    X max(Xcr,Xcr); X min(Xcr,Xcr); PX abs(Xcr); \
    X sqrt(Xcr); X exp(Xcr); X log(Xcr); X sin(Xcr); X cos(Xcr); X tan(Xcr); X atan(Xcr); \


namespace Ariadne {

// Declare friend functions in namespace
Rational operator/(Integer const&, Integer const&);
Integer pow(Integer const& z, Nat m);
//Integer abs(Integer const& z);

Rational operator/(Dyadic const& w1, Dyadic const& w2) { return Rational(w1) / Rational(w2); }
Dyadic hlf(Dyadic const& w);

DECLARE_NUMERIC_OPERATIONS(Real const&,Real,PositiveReal);
//DECLARE_NUMERIC_OPERATIONS(FloatDPApproximation const&,FloatDPApproximation,PositiveFloatDPApproximation);
//DECLARE_NUMERIC_OPERATIONS(FloatBall<PR>);
//DECLARE_NUMERIC_OPERATIONS(FloatBounds<PR>);
//DECLARE_NUMERIC_OPERATIONS(FloatApproximation<PR>);

template<> String class_name<DoublePrecision>() { return "DoublePrecision"; }
template<> String class_name<MultiplePrecision>() { return "MultiplePrecision"; }

template<class T> String class_tag();
template<> String class_tag<DoublePrecision>() { return "DP"; }
template<> String class_tag<MultiplePrecision>() { return "MP"; }

template<class T> struct PythonName { const char* get() const { return class_name<T>().c_str(); } };
template<class T> inline const char* python_name() { return PythonName<T>().get(); }

} // namespace Ariadne


using namespace boost::python;


namespace Ariadne {

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

template<class OP> struct PythonOperator { };
PythonOperator<Pos> pos(boost::python::self_ns::self_t) { return PythonOperator<Pos>(); }
PythonOperator<Neg> neg(boost::python::self_ns::self_t) { return PythonOperator<Neg>(); }
PythonOperator<Sqr> sqr(boost::python::self_ns::self_t) { return PythonOperator<Sqr>(); }
PythonOperator<Hlf> hlf(boost::python::self_ns::self_t) { return PythonOperator<Hlf>(); }
PythonOperator<Rec> rec(boost::python::self_ns::self_t) { return PythonOperator<Rec>(); }
PythonOperator<Sqrt> sqrt(boost::python::self_ns::self_t) { return PythonOperator<Sqrt>(); }
PythonOperator<Exp> exp(boost::python::self_ns::self_t) { return PythonOperator<Exp>(); }
PythonOperator<Log> log(boost::python::self_ns::self_t) { return PythonOperator<Log>(); }
PythonOperator<Sin> sin(boost::python::self_ns::self_t) { return PythonOperator<Sin>(); }
PythonOperator<Cos> cos(boost::python::self_ns::self_t) { return PythonOperator<Cos>(); }
PythonOperator<Tan> tan(boost::python::self_ns::self_t) { return PythonOperator<Tan>(); }
PythonOperator<Atan> atan(boost::python::self_ns::self_t) { return PythonOperator<Atan>(); }


template<class OP, class T> auto py_apply(T const& t) -> decltype(OP()(t)){ OP op; return op(t); }

template<class... T> struct Tag { };

template<class T, class B = boost::python::bases<> > class class_ : public boost::python::class_<T,B> {
    typedef T const& Tcr;
  public:
    class_(String const& name)
        : boost::python::class_<T,B>(name.c_str()) { }
    template<class... U> class_(String const& name, boost::python::init<U...> const& initialiser)
        : boost::python::class_<T,B>(name.c_str(),initialiser) { }
    using boost::python::class_<T,B>::def;
    template<class OP> void def(PythonOperator<OP>) { boost::python::def(to_str(OP()).c_str(),&py_apply<OP,T>); }
    void define_self_arithmetic() {
        using boost::python::self;
        this->def(+self); this->def(-self);
        this->def(self+self); this->def(self-self); this->def(self*self); this->def(self/self);
    }
    template<class X> void define_mixed_arithmetic() {
        using boost::python::self_ns::self;
        X const* other_ptr=nullptr; X const& other=*other_ptr;
        this->def(self+other); this->def(self-other); this->def(self*other); this->def(self/other);
        this->def(other+self); this->def(other-self); this->def(other*self); this->def(other/self);
    }
    void define_unary_arithmetic() {
        using boost::python::self_ns::self;
        this->def(+self); this->def(-self);
        using boost::python::def;
        this->def(pos(self)); this->def(neg(self));
        this->def(sqr(self)); this->def(rec(self));
    }
    void define_transcendental_functions() {
        this->define_unary_arithmetic();
        this->def(sqrt(self)); this->def(exp(self)); this->def(log(self));
        this->def(sin(self)); this->def(cos(self)); this->def(tan(self)); this->def(atan(self));
    }
    void define_monotonic_functions() {
        typedef decltype(-declval<T>()) NT;
        using boost::python::def;
        this->def(pos(self)); this->def(neg(self));
        this->def(sqrt(self)); this->def(exp(self)); this->def(log(self)); this->def(atan(self));
    }
    void define_self_comparisons() {
        using boost::python::self_ns::self;
        this->def(self==self); this->def(self!=self);
        this->def(self<=self); this->def(self>=self); this->def(self<self); this->def(self>self);
    }
    void define_self_logical() {
        using boost::python::self_ns::self;
        this->def(self & self); this->def(self | self); this->def(~self);
    }
    void convert_from() { }
    template<class A, class... AS> void convert_from() {
        this->def(boost::python::init<A>()); boost::python::implicitly_convertible<A,T>; this->convert_from<AS...>();
    }
    void inits() { }
    template<class A, class... AS> void inits() {
        this->def(boost::python::init<A>()); this->inits<AS...>();
    }
};


template<class P> Bool decide(Logical<P> l) { return Ariadne::decide(l); }
template<class P> Bool definitely(Logical<P> l) { return Ariadne::definitely(l); }
template<class P> Bool possibly(Logical<P> l) { return Ariadne::possibly(l); }

template<class P> auto check(Logical<P> l, Effort e) -> decltype(l.check(e)) { return l.check(e); }
template<class P1, class P2> Logical<Weaker<P1,P2>> operator&(Logical<P1> l1, Logical<P2> l2) { return l1 and l2; }
template<class P1, class P2> Logical<Weaker<P1,P2>> operator|(Logical<P1> l1, Logical<P2> l2) { return l1 or l2; }
template<class P> Logical<Opposite<P>> operator~(Logical<P> l) { return not l; }

template<class P> void export_effective_logical(std::string name)
{
    typedef decltype(declval<Logical<P>>().check(declval<Effort>())) CheckType;
    OutputStream& operator<<(OutputStream& os, Logical<P> l);

    class_<Logical<P>> logical_class(name,init<bool>());
    logical_class.def(init<Logical<P>>());
    logical_class.def("check", (CheckType(Logical<P>::*)(Effort)) &Logical<P>::check);
    logical_class.define_self_logical();
    logical_class.def(self_ns::str(self));
    logical_class.def(self_ns::repr(self));
    def("check", (CheckType(*)(Logical<P> const&,Effort)) &Ariadne::check<P>);
};

template<class P> void export_logical(std::string name)
{
    typedef decltype(~declval<Logical<P>>()) NotType;
    OutputStream& operator<<(OutputStream& os, Logical<P> l);
    class_<Logical<P>> logical_class(name,init<bool>());
    logical_class.def(init<Logical<P>>());
    logical_class.define_self_logical();
    logical_class.def(self_ns::str(self));
    logical_class.def(self_ns::repr(self));

    def("decide", (bool(*)(Logical<P>)) &decide);
    def("possibly", (bool(*)(Logical<P>)) &possibly);
    def("definitely", (bool(*)(Logical<P>)) &definitely);

};

template<> void export_logical<ExactTag>(std::string name)
{
    typedef ExactTag P;
    OutputStream& operator<<(OutputStream& os, Logical<P> l);
    class_<Logical<P>> logical_class(name,init<bool>());
    logical_class.def(init<Logical<P>>());
    logical_class.define_self_logical();
    logical_class.def(self_ns::str(self));
    logical_class.def(self_ns::repr(self));

    implicitly_convertible<Logical<ExactTag>,bool>();

}

void export_integer()
{
    class_<Integer> integer_class("Integer");
    integer_class.def(init<int>());
    integer_class.def(init<Integer>());

    integer_class.def(self_ns::str(self));
    integer_class.def(self_ns::repr(self));

//    integer_class.define_self_arithmetic();

    // Required for mixed operations with integral types
    integer_class.define_mixed_arithmetic<Integer>();
    integer_class.define_self_comparisons();

    def("pow", (Integer(*)(const Integer&,Nat)) &pow) ;
    integer_class.def(abs(self));


    implicitly_convertible<int,Integer>();

}

void export_dyadic()
{
    class_<Dyadic> dyadic_class("Dyadic");
//    dyadic_class.def(init<Int,Nat>());
//    dyadic_class.def(init<Integer,Natural>());
//    dyadic_class.def(init<Int>());
    dyadic_class.def(init<Integer>());
    dyadic_class.def(init<Dyadic>());
    dyadic_class.def(init<FloatDPValue>());
    dyadic_class.def(init<FloatMPValue>());

    dyadic_class.define_self_arithmetic();
    dyadic_class.define_self_comparisons();
    def("hlf", (Dyadic(*)(Dyadic const&)) &hlf);

    dyadic_class.def(self_ns::str(self));
    dyadic_class.def(self_ns::repr(self));

    dyadic_class.def("get_d", &Dyadic::get_d);

    implicitly_convertible<Integer,Dyadic>();
}

void export_rational()
{
    class_<Rational> rational_class("Rational");
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

    rational_class.def(self_ns::str(self));
    rational_class.def(self_ns::repr(self));

    rational_class.def("get_d", &Rational::get_d);

    implicitly_convertible<Integer,Rational>();
    implicitly_convertible<Dyadic,Rational>();
}

void export_decimal()
{
    class_<Decimal> decimal_class("Decimal", init<Decimal>());
    decimal_class.def(self_ns::str(self));
    decimal_class.def(self_ns::repr(self));
    implicitly_convertible<Decimal,Rational>();
}

void export_real()
{
    class_<Real> real_class("Real");
    real_class.def(init<int>());
    real_class.def(init<Integer>());
    real_class.def(init<Dyadic>());
    real_class.def(init<Rational>());
    real_class.def(init<Real>());

    real_class.define_self_arithmetic();
    real_class.template define_mixed_arithmetic<Int>();
    real_class.template define_mixed_arithmetic<Rational>();
    real_class.define_transcendental_functions();
    real_class.define_self_comparisons();

    def("sqrt", (Real(*)(Real const&)) &sqrt);
    def("exp", (Real(*)(Real const&)) &exp);
    def("log", (Real(*)(Real const&)) &log);
    def("sin", (Real(*)(Real const&)) &sin);
    def("cos", (Real(*)(Real const&)) &cos);
    def("tan", (Real(*)(Real const&)) &tan);

    real_class.def(self_ns::str(self));
    real_class.def(self_ns::repr(self));

    real_class.def("get", (FloatDPBounds(Real::*)(DoublePrecision)const) &Real::get);
    real_class.def("get", (FloatMPBounds(Real::*)(MultiplePrecision)const) &Real::get);
    real_class.def("get", (FloatMPBall(Real::*)(Accuracy)const) &Real::get);
    real_class.def("evaluate", (FloatMPBall(Real::*)(Accuracy)const) &Real::evaluate);
    real_class.def("get_d", &Real::get_d);

    implicitly_convertible<Rational,Real>();
}


template<class P> void export_number()
{
    class_<Number<P>> number_class(class_name<P>()+"Number");
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

void export_numbers()
{
    class_<ApproximateNumber> approximate_number_class(class_name<ApproximateTag>()+"Number");
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
    approximate_number_class.def(self_ns::str(self));
    approximate_number_class.def(self_ns::repr(self));

    class_<ValidatedLowerNumber> lower_number_class(class_name<LowerTag>()+"Number");
    lower_number_class.def(init<ValidatedNumber>());
    lower_number_class.def("get", (FloatDPLowerBound(*)(ValidatedLowerNumber const&, DoublePrecision const&)) &get);
    lower_number_class.def("get", (FloatMPLowerBound(*)(ValidatedLowerNumber const&, MultiplePrecision const&)) &get);
    lower_number_class.define_monotonic_functions();
    lower_number_class.def(self_ns::str(self));
    lower_number_class.def(self_ns::repr(self));

    class_<ValidatedUpperNumber> upper_number_class(class_name<UpperTag>()+"Number");
    upper_number_class.def(init<ValidatedNumber>());
    upper_number_class.def("get", (FloatDPUpperBound(*)(ValidatedUpperNumber const&, DoublePrecision const&)) &get);
    upper_number_class.def("get", (FloatMPUpperBound(*)(ValidatedUpperNumber const&, MultiplePrecision const&)) &get);
    upper_number_class.define_monotonic_functions();
    upper_number_class.def(self_ns::str(self));
    upper_number_class.def(self_ns::repr(self));

    class_<ValidatedNumber> validated_number_class(class_name<ValidatedTag>()+"Number");
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
    validated_number_class.def(self_ns::str(self));
    validated_number_class.def(self_ns::repr(self));

    class_<EffectiveNumber> effective_number_class(class_name<EffectiveTag>()+"Number");
    effective_number_class.def(init<Rational>());
    effective_number_class.def(init<Real>());
    effective_number_class.def(init<ExactNumber>());
    effective_number_class.def(init<EffectiveNumber>());
    effective_number_class.def("get", (FloatDPBounds(*)(EffectiveNumber const&, DoublePrecision const&)) &get);
    effective_number_class.def("get", (FloatMPBounds(*)(EffectiveNumber const&, MultiplePrecision const&)) &get);
    effective_number_class.define_self_arithmetic();
    effective_number_class.define_mixed_arithmetic<EffectiveNumber>();
    effective_number_class.define_transcendental_functions();
    effective_number_class.def(self_ns::str(self));
    effective_number_class.def(self_ns::repr(self));

    class_<ExactNumber> exact_number_class(class_name<ExactTag>()+"Number");
    exact_number_class.def(init<Rational>());
    exact_number_class.def(init<ExactNumber>());
    exact_number_class.def("get", (FloatDPBounds(*)(ExactNumber const&, DoublePrecision const&)) &get);
    exact_number_class.def("get", (FloatMPBounds(*)(ExactNumber const&, MultiplePrecision const&)) &get);
    exact_number_class.def(self_ns::str(self));
    exact_number_class.def(self_ns::repr(self));


    implicitly_convertible<ValidatedNumber,ApproximateNumber>();
    implicitly_convertible<ValidatedNumber,ValidatedLowerNumber>();
    implicitly_convertible<ValidatedNumber,ValidatedUpperNumber>();
    implicitly_convertible<EffectiveNumber,ValidatedNumber>();
    implicitly_convertible<ExactNumber,EffectiveNumber>();

    implicitly_convertible<Rational,ExactNumber>();
    implicitly_convertible<Real,EffectiveNumber>();
}

template<class PR> void export_raw_float()
{
    class_<RawFloat<PR>> raw_float_class("Float"+class_tag<PR>());
    raw_float_class.def(init<double,PR>());
    raw_float_class.def(self_ns::str(self));
    raw_float_class.def(self_ns::repr(self));
}

template<class PR> void export_float_value()
{
    class_<FloatValue<PR>> float_value_class("Float"+class_tag<PR>()+"Value");
    float_value_class.def(init<RawFloat<PR>>());
    float_value_class.def(init<int>());
    float_value_class.def(init<double>());
    float_value_class.def(init<Integer,PR>());
    float_value_class.def(init<FloatValue<PR>>());

    float_value_class.def(pos(self));
    float_value_class.def(neg(self));

    float_value_class.def("precision", &FloatValue<PR>::precision);
    float_value_class.def("get_d",&FloatValue<PR>::get_d);

    float_value_class.template define_mixed_arithmetic<Int>();
    float_value_class.define_self_arithmetic();
    //float_value_class.define_mixed_arithmetic<ApproximateNumericType>();
    //float_value_class.define_mixed_arithmetic<LowerNumericType>();
    //float_value_class.define_mixed_arithmetic<UpperNumericType>();
    //float_value_class.define_mixed_arithmetic<ValidatedNumericType>();
    float_value_class.define_self_comparisons();

    float_value_class.def(pos(self));
    float_value_class.def(neg(self));

    float_value_class.def("precision", &FloatValue<PR>::precision);
    float_value_class.def("raw", (RawFloat<PR>const&(FloatValue<PR>::*)()const)&FloatValue<PR>::raw, return_value_policy<copy_const_reference>());
    float_value_class.def("get_d",&FloatValue<PR>::get_d);

    float_value_class.def(self_ns::str(self));
    float_value_class.def(self_ns::repr(self));

    float_value_class.def("set_output_places",&FloatValue<PR>::set_output_places).staticmethod("set_output_places");
}

FloatUpperBound<DoublePrecision> log2(FloatError<DoublePrecision> const& x);
FloatUpperBound<MultiplePrecision> log2(FloatError<MultiplePrecision> const& x);

template<class PRE> void export_float_error()
{
    class_<FloatError<PRE>> float_error_class("Float"+class_tag<PRE>()+"Error");
    float_error_class.def(init<RawFloat<PRE>>());
    float_error_class.def(init<uint>());
    float_error_class.def(init<double>());
    float_error_class.def(init<FloatError<PRE>>());

    float_error_class.def(+self);
    float_error_class.def(self+self);
    float_error_class.def(self*self);

    //float_error_class.def("raw",(RawFloat<PRE>&(FloatError<PRE>::*)())&FloatError<PRE>::raw,return_value_policy<reference_existing_object>());
    float_error_class.def("raw",(RawFloat<PRE>const&(FloatError<PRE>::*)()const)&FloatError<PRE>::raw,return_value_policy<copy_const_reference>());
    float_error_class.def(self_ns::str(self));
    float_error_class.def(self_ns::repr(self));

    float_error_class.def("precision", &FloatError<PRE>::precision);

    def("log2", (FloatUpperBound<PRE>(*)(FloatError<PRE>const&)) &log2);
    //    float_error_class.def("set_output_places",&FloatError<PRE>::set_output_places).staticmethod("set_output_places");
}

template<class PR, class PRE=PR> void export_float_ball()
{
    String tag = IsSame<PR,PRE>::value ? class_tag<PR>() : class_tag<PR>()+class_tag<PRE>();
    class_<FloatBall<PR,PRE>> float_ball_class("Float"+tag+"Ball");
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

    float_ball_class.def(self_ns::str(self));
    float_ball_class.def(self_ns::repr(self));
};



template<class PR> void export_float_bounds()
{
    class_<FloatBounds<PR>> float_bounds_class("Float"+class_tag<PR>()+"Bounds");
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

    float_bounds_class.def(self_ns::str(self));
    float_bounds_class.def(self_ns::repr(self));

    implicitly_convertible<FloatBall<PR>,FloatBounds<PR>>();

    float_bounds_class.def("set_output_places",&FloatBounds<PR>::set_output_places).staticmethod("set_output_places");

}

template<class PR> void export_float_upper_bound()
{
    class_<FloatUpperBound<PR>> float_upper_bound_class("Float"+class_tag<PR>()+"UpperBound");
    float_upper_bound_class.def(init<RawFloat<PR>>());
    float_upper_bound_class.def(init<int>());
    float_upper_bound_class.def(init<double>());
    float_upper_bound_class.def(init<Real,PR>());

    float_upper_bound_class.def(init<FloatValue<PR>>());
    float_upper_bound_class.def(init<FloatBall<PR>>());
    float_upper_bound_class.def(init<FloatBounds<PR>>());
    float_upper_bound_class.def(init<FloatUpperBound<PR>>());
    float_upper_bound_class.def(init<ValidatedUpperNumber,PR>());
    float_upper_bound_class.def(self_ns::str(self));
    float_upper_bound_class.def(self_ns::repr(self));

    float_upper_bound_class.def("precision", &FloatUpperBound<PR>::precision);
    float_upper_bound_class.def("raw", (RawFloat<PR>const&(FloatUpperBound<PR>::*)()const)&FloatUpperBound<PR>::raw, return_value_policy<copy_const_reference>());

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

template<class PR> void export_float_lower_bound()
{
    class_<FloatLowerBound<PR>> float_lower_bound_class("Float"+class_tag<PR>()+"LowerBound");
    float_lower_bound_class.def(init<RawFloat<PR>>());
    float_lower_bound_class.def(init<int>());
    float_lower_bound_class.def(init<double>());
    float_lower_bound_class.def(init<Real,PR>());

    float_lower_bound_class.def(init<FloatValue<PR>>());
    float_lower_bound_class.def(init<FloatBall<PR>>());
    float_lower_bound_class.def(init<FloatBounds<PR>>());
    float_lower_bound_class.def(init<FloatLowerBound<PR>>());
    float_lower_bound_class.def(init<ValidatedLowerNumber,PR>());
    float_lower_bound_class.def(self_ns::str(self));
    float_lower_bound_class.def(self_ns::repr(self));

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
    float_lower_bound_class.def("raw", (RawFloat<PR>const&(FloatLowerBound<PR>::*)()const)&FloatLowerBound<PR>::raw, return_value_policy<copy_const_reference>());

    //    float_lower_bound_class.define_mixed_arithmetic<ApproximateNumericType>();

//    float_lower_bound_class.def(LowerNumericType() + self);
//    float_lower_bound_class.def(UpperNumericType() - self);
//    float_lower_bound_class.def(self + LowerNumericType());
//    float_lower_bound_class.def(self - UpperNumericType());

    float_lower_bound_class.define_monotonic_functions();

    implicitly_convertible<FloatBounds<PR>,FloatLowerBound<PR>>();
}



template<class PR> void export_float_approximation()
{
    class_<FloatApproximation<PR>> float_approximation_class("Float"+class_tag<PR>()+"Approximation");

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
    float_approximation_class.def("raw", (RawFloat<PR>const&(FloatApproximation<PR>::*)()const)&FloatApproximation<PR>::raw, return_value_policy<copy_const_reference>());
    float_approximation_class.def("get_d", &FloatApproximation<PR>::get_d);

    float_approximation_class.def(self_ns::str(self));
    float_approximation_class.def(self_ns::repr(self));

    implicitly_convertible<double,FloatApproximation<PR>>();
    //implicitly_convertible<Integer,FloatBounds<PR>>();
    implicitly_convertible<FloatValue<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatBall<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatBounds<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatUpperBound<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatLowerBound<PR>,FloatApproximation<PR>>();

    float_approximation_class.def(abs(self));

    float_approximation_class.define_unary_arithmetic();
    float_approximation_class.define_transcendental_functions();
}

Void export_effort() {
    class_<Effort> effort_class("Effort",init<Nat>());
    effort_class.def(init<>());
    effort_class.def(self_ns::str(self));
}

Void export_accuracy() {
    class_<Accuracy> accuracy_class("Accuracy",init<Nat>());
    accuracy_class.def("bits",&Accuracy::bits);
    accuracy_class.def(self_ns::str(self));
    accuracy_class.def(self_ns::repr(self));
}

template<class PR> Void export_precision();

template<> Void export_precision<DoublePrecision>() {
    class_<DoublePrecision> precision_class("DoublePrecision",init<>());
    precision_class.def(self_ns::str(self));
}

template<> Void export_precision<MultiplePrecision>() {
    class_<MultiplePrecision> precision_class("MultiplePrecision",init<Nat>());
    precision_class.def("bits",&MultiplePrecision::bits);
    precision_class.def(self_ns::str(self));
}

template<class PR> Void export_user_floats() {
    export_float_approximation<PR>();
    export_float_upper_bound<PR>();
    export_float_lower_bound<PR>();
    export_float_bounds<PR>();
    export_float_ball<PR>();
    export_float_value<PR>();
    export_float_error<PR>();
}


} // namespace Ariadne;

void
numeric_submodule()
{
    using namespace Ariadne;
    export_effort();
    export_accuracy();

    export_logical<ExactTag>("Boolean");
    export_effective_logical<EffectiveTag>("Kleenean");
    export_effective_logical<EffectiveUpperTag>("Sierpinskian");
    export_effective_logical<EffectiveLowerTag>("NegatedSierpinskian");
    export_logical<ValidatedTag>("Tribool");
    export_logical<UpperTag>("Verified");
    export_logical<LowerTag>("Falsified");
    export_logical<ApproximateTag>("Fuzzy");

    export_precision<DoublePrecision>();
    export_precision<MultiplePrecision>();
    export_raw_float<DoublePrecision>();
    export_raw_float<MultiplePrecision>();
    export_user_floats<DoublePrecision>();
    export_user_floats<MultiplePrecision>();

    export_float_ball<MultiplePrecision,DoublePrecision>();

    export_numbers();
    export_real();
    export_rational();
    export_decimal();
    export_dyadic();
    export_integer();
}
