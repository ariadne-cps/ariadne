/***************************************************************************
 *            numeric_submodule.cc
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

#include "numeric/logical.h"
#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/real.h"
#include "numeric/float64.h"
#include "numeric/floatmp.h"
#include "numeric/float-user.h"

#define DECLARE_NUMERIC_OPERATIONS(X,PX) \
    X add(X,X); X sub(X,X); X mul(X,X); X div(X,X); \
    X pos(X); X neg(X); X sqr(X); X rec(X); X pow(X,Int); \
    X max(X,X); X min(X,X); PX abs(X); \
    X sqrt(X); X exp(X); X log(X); X sin(X); X cos(X); X tan(X); X atan(X); \


namespace Ariadne {

// Declare friend functions in namespace
Rational operator/(Integer const&, Integer const&);
Integer pow(Integer const& z, Nat m);
//Integer abs(Integer const& z);

DECLARE_NUMERIC_OPERATIONS(Real,PositiveReal);
//DECLARE_NUMERIC_OPERATIONS(FloatBall<PR>);
//DECLARE_NUMERIC_OPERATIONS(FloatBounds<PR>);
//DECLARE_NUMERIC_OPERATIONS(FloatApproximation<PR>64);

template<> String class_name<Precision64>() { return "Precision64"; }
template<> String class_name<PrecisionMP>() { return "PrecisionMP"; }

template<class T> String class_tag();
template<> String class_tag<Precision64>() { return "64"; }
template<> String class_tag<PrecisionMP>() { return "MP"; }

template<class T> struct PythonName { const char* get() const { return class_name<T>().c_str(); } };
template<class T> inline const char* python_name() { return PythonName<T>().get(); }

} // namespace Ariadne


using namespace boost::python;


namespace Ariadne {

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Rational>& repr) {
    return os << "Rational("<<repr.reference().numerator()<<","<<repr.reference().denominator()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<RawFloat64>& repr) {
    return os << "Float64("<<repr.reference()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Float64Approximation>& repr) {
    return os << "Float64Approximation("<<repr.reference().raw()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Float64Bounds>& repr) {
    return os << "Float64Bounds("<<repr.reference().lower().raw()<<","<<repr.reference().upper().raw()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Float64Value>& repr) {
    return os << "Float64Value("<<repr.reference().raw()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Float64UpperBound>& repr) {
    return os << "Float64UpperBound("<<repr.reference().raw()<<")"; }
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<PositiveFloat64UpperBound>& repr) {
    return os << "Float64Error("<<repr.reference().raw()<<")"; }

template<class OP> struct PythonOperator { };
PythonOperator<Pos> pos(boost::python::self_ns::self_t) { return PythonOperator<Pos>(); }
PythonOperator<Neg> neg(boost::python::self_ns::self_t) { return PythonOperator<Neg>(); }
PythonOperator<Sqr> sqr(boost::python::self_ns::self_t) { return PythonOperator<Sqr>(); }
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

void export_rational()
{
    class_<Rational> rational_class("Rational");
    rational_class.def(init<int,int>());
    rational_class.def(init<Integer,Integer>());
    rational_class.def(init<int>());
    //rational_class.def(init<double>());
    rational_class.def(init<Integer>());
    rational_class.def(init<Rational>());

    rational_class.define_self_arithmetic();
    rational_class.define_self_comparisons();
    rational_class.define_mixed_arithmetic<double>();
    rational_class.define_mixed_arithmetic<Rational>();

    rational_class.def(self_ns::str(self));
    rational_class.def(self_ns::repr(self));

    rational_class.def("get_d", &Rational::get_d);

    implicitly_convertible<Integer,Rational>();

}

void export_real()
{
    class_<Real> real_class("Real");
    real_class.def(init<int>());
    real_class.def(init<Integer>());
    real_class.def(init<Rational>());
    real_class.def(init<Real>());

    real_class.define_self_arithmetic();
    real_class.template define_mixed_arithmetic<Int>();
    real_class.template define_mixed_arithmetic<Rational>();
    real_class.define_transcendental_functions();
    real_class.define_self_comparisons();

    def("exp", (Real(*)(Real)) &exp);

    real_class.def(self_ns::str(self));
    real_class.def(self_ns::repr(self));

    real_class.def("get", (Float64Bounds(Real::*)(Precision64)const) &Real::get);
    real_class.def("get", (FloatMPBounds(Real::*)(PrecisionMP)const) &Real::get);
    real_class.def("get", (FloatMPBounds(Real::*)(Accuracy)const) &Real::evaluate);
    real_class.def("evaluate", (FloatMPBounds(Real::*)(Accuracy)const) &Real::evaluate);
    real_class.def("get_d", &Real::get_d);

    implicitly_convertible<Rational,Real>();
}


template<class P> void export_number()
{
    class_<Number<P>> number_class(class_name<P>()+"Number");
    number_class.def(init<Rational>());
    number_class.define_self_arithmetic();
//    number_class.def("get", (Float64Bounds(Number<P>::*)(ValidatedTag,Precision64)const) &Number<P>::get);
//    number_class.def("get", (Float64Approximation(Number<P>::*)(ApproximateTag,Precision64)const) &Number<P>::get);
//    number_class.def("get", (FloatMPBounds(Number<P>::*)(ValidatedTag,PrecisionMP)const) &Number<P>::get);
//    number_class.def("get", (FloatMPApproximation(Number<P>::*)(ApproximateTag,PrecisionMP)const) &Number<P>::get);
}

Float64Bounds get(ExactNumber const& n, Precision64 const& pr) { return n.get(BoundedTag(),pr); }
FloatMPBounds get(ExactNumber const& n, PrecisionMP const& pr) { return n.get(BoundedTag(),pr); }
Float64Bounds get(EffectiveNumber const& n, Precision64 const& pr) { return n.get(BoundedTag(),pr); }
FloatMPBounds get(EffectiveNumber const& n, PrecisionMP const& pr) { return n.get(BoundedTag(),pr); }
Float64Bounds get(ValidatedNumber const& n, Precision64 const& pr) { return n.get(BoundedTag(),pr); }
FloatMPBounds get(ValidatedNumber const& n, PrecisionMP const& pr) { return n.get(BoundedTag(),pr); }
Float64UpperBound get(ValidatedUpperNumber const& n, Precision64 const& pr) { return n.get(UpperTag(),pr); }
FloatMPUpperBound get(ValidatedUpperNumber const& n, PrecisionMP const& pr) { return n.get(UpperTag(),pr); }
Float64LowerBound get(ValidatedLowerNumber const& n, Precision64 const& pr) { return n.get(LowerTag(),pr); }
FloatMPLowerBound get(ValidatedLowerNumber const& n, PrecisionMP const& pr) { return n.get(LowerTag(),pr); }
Float64Approximation get(ApproximateNumber const& n, Precision64 const& pr) { return n.get(ApproximateTag(),pr); }
FloatMPApproximation get(ApproximateNumber const& n, PrecisionMP const& pr) { std::cerr<<"get(AN,MP)\n";return n.get(ApproximateTag(),pr); }

void export_numbers()
{
    class_<ApproximateNumber> approximate_number_class(class_name<ApproximateTag>()+"Number");
    approximate_number_class.def(init<Rational>());
    approximate_number_class.def(init<Real>());
    approximate_number_class.def(init<ExactNumber>());
    approximate_number_class.def(init<EffectiveNumber>());
    approximate_number_class.def(init<ValidatedNumber>());
    approximate_number_class.def(init<ApproximateNumber>());
    approximate_number_class.def("get", (Float64Approximation(*)(ApproximateNumber const&, Precision64 const&)) &get);
    approximate_number_class.def("get", (FloatMPApproximation(*)(ApproximateNumber const&, PrecisionMP const&)) &get);
    approximate_number_class.define_self_arithmetic();
    approximate_number_class.define_mixed_arithmetic<ApproximateNumber>();
    approximate_number_class.define_transcendental_functions();
    approximate_number_class.def(self_ns::str(self));
    approximate_number_class.def(self_ns::repr(self));

    class_<ValidatedLowerNumber> lower_number_class(class_name<LowerTag>()+"Number");
    lower_number_class.def(init<ValidatedNumber>());
    lower_number_class.def("get", (Float64LowerBound(*)(ValidatedLowerNumber const&, Precision64 const&)) &get);
    lower_number_class.def("get", (FloatMPLowerBound(*)(ValidatedLowerNumber const&, PrecisionMP const&)) &get);
    lower_number_class.define_monotonic_functions();
    lower_number_class.def(self_ns::str(self));
    lower_number_class.def(self_ns::repr(self));

    class_<ValidatedUpperNumber> upper_number_class(class_name<UpperTag>()+"Number");
    upper_number_class.def(init<ValidatedNumber>());
    upper_number_class.def("get", (Float64UpperBound(*)(ValidatedUpperNumber const&, Precision64 const&)) &get);
    upper_number_class.def("get", (FloatMPUpperBound(*)(ValidatedUpperNumber const&, PrecisionMP const&)) &get);
    upper_number_class.define_monotonic_functions();
    upper_number_class.def(self_ns::str(self));
    upper_number_class.def(self_ns::repr(self));

    class_<ValidatedNumber> validated_number_class(class_name<ValidatedTag>()+"Number");
    validated_number_class.def(init<Rational>());
    validated_number_class.def(init<Real>());
    validated_number_class.def(init<ExactNumber>());
    validated_number_class.def(init<EffectiveNumber>());
    validated_number_class.def(init<ValidatedNumber>());
    validated_number_class.def("get", (Float64Bounds(*)(ValidatedNumber const&, Precision64 const&)) &get);
    validated_number_class.def("get", (FloatMPBounds(*)(ValidatedNumber const&, PrecisionMP const&)) &get);
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
    effective_number_class.def("get", (Float64Bounds(*)(EffectiveNumber const&, Precision64 const&)) &get);
    effective_number_class.def("get", (FloatMPBounds(*)(EffectiveNumber const&, PrecisionMP const&)) &get);
    effective_number_class.define_self_arithmetic();
    effective_number_class.define_mixed_arithmetic<EffectiveNumber>();
    effective_number_class.define_transcendental_functions();
    effective_number_class.def(self_ns::str(self));
    effective_number_class.def(self_ns::repr(self));

    class_<ExactNumber> exact_number_class(class_name<ExactTag>()+"Number");
    exact_number_class.def(init<Rational>());
    exact_number_class.def(init<ExactNumber>());
    exact_number_class.def("get", (Float64Bounds(*)(ExactNumber const&, Precision64 const&)) &get);
    exact_number_class.def("get", (FloatMPBounds(*)(ExactNumber const&, PrecisionMP const&)) &get);
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

template<class PR> void export_exact_float()
{
    class_<FloatValue<PR>> exact_float_class("Float"+class_tag<PR>()+"Value");
    exact_float_class.def(init<int>());
    exact_float_class.def(init<double>());
    exact_float_class.def(init<Integer,PR>());
    exact_float_class.def(init<FloatValue<PR>>());

    exact_float_class.template define_mixed_arithmetic<Int>();
    exact_float_class.template define_mixed_arithmetic<FloatValue<PR>>();
    exact_float_class.define_self_arithmetic();
    //exact_float_class.define_mixed_arithmetic<ApproximateNumericType>();
    //exact_float_class.define_mixed_arithmetic<LowerNumericType>();
    //exact_float_class.define_mixed_arithmetic<UpperNumericType>();
    //exact_float_class.define_mixed_arithmetic<ValidatedNumericType>();
    exact_float_class.define_self_comparisons();

    def("pos", (FloatValue<PR>(*)(FloatValue<PR> const&)) &pos);
    def("neg", (FloatValue<PR>(*)(FloatValue<PR> const&)) &neg);

    exact_float_class.def("precision", &FloatValue<PR>::precision);
    exact_float_class.def("get_d",&FloatValue<PR>::get_d);

    exact_float_class.def(self_ns::str(self));
    exact_float_class.def(self_ns::repr(self));

    implicitly_convertible<int,FloatValue<PR>>();

    exact_float_class.def("set_output_precision",&FloatValue<PR>::set_output_precision).staticmethod("set_output_precision");
}

template<class PR> void export_error_float()
{
    class_<FloatError<PR>> error_float_class("Float"+class_tag<PR>()+"Error");
    error_float_class.def(init<uint>());
    error_float_class.def(init<double>());
    error_float_class.def(init<FloatError<PR>>());

    error_float_class.def(+self);
    error_float_class.def(self+self);
    error_float_class.def(self*self);
    error_float_class.def("get_d",&FloatError<PR>::get_d);
    error_float_class.def(self_ns::str(self));
    error_float_class.def(self_ns::repr(self));

    error_float_class.def("precision", &FloatError<PR>::precision);

    //    error_float_class.def("set_output_precision",&FloatError<PR>::set_output_precision).staticmethod("set_output_precision");
}

template<class PR> void export_metric_float()
{
    class_<FloatBall<PR>> metric_float_class("Float"+class_tag<PR>()+"Ball");
    metric_float_class.def(init<double,double>());
    metric_float_class.def(init<FloatValue<PR>,FloatError<PR>>());
    metric_float_class.def(init<Real,PR>());

    metric_float_class.def(init<double>());
    metric_float_class.def(init<ValidatedNumericType>());
    metric_float_class.def(init<FloatValue<PR>>());
    metric_float_class.def(init<FloatBall<PR>>());
    metric_float_class.def(init<FloatBounds<PR>>());

    metric_float_class.def("value", &FloatBall<PR>::value);
    metric_float_class.def("error", &FloatBall<PR>::error);
    metric_float_class.def("lower", &FloatBall<PR>::lower);
    metric_float_class.def("upper", &FloatBall<PR>::upper);

    metric_float_class.def("precision", &FloatBall<PR>::precision);

    metric_float_class.template define_mixed_arithmetic<Int>();
    metric_float_class.template define_mixed_arithmetic<FloatBall<PR>>();
    metric_float_class.template define_mixed_arithmetic<ValidatedNumber>();
    metric_float_class.define_self_arithmetic();
//    metric_float_class.define_mixed_arithmetic<FloatBall<PR>>();
//    metric_float_class.define_mixed_arithmetic<ApproximateNumericType>();
//    metric_float_class.define_mixed_arithmetic<LowerNumericType>();
//    metric_float_class.define_mixed_arithmetic<UpperNumericType>();
//    metric_float_class.define_mixed_arithmetic<ValidatedNumericType>();
    metric_float_class.define_self_comparisons();

    metric_float_class.define_unary_arithmetic();
    metric_float_class.define_transcendental_functions();

    metric_float_class.def(self_ns::str(self));
    metric_float_class.def(self_ns::repr(self));

    implicitly_convertible<FloatValue<PR>,FloatBall<PR>>();

}


template<class PR> void export_bounded_float()
{
    class_<FloatBounds<PR>> bounded_float_class("Float"+class_tag<PR>()+"Bounds");
    bounded_float_class.def(init<double,double>());
    bounded_float_class.def(init<FloatLowerBound<PR>,FloatUpperBound<PR>>());
    bounded_float_class.def(init<Real,PR>());

    bounded_float_class.def(init<double>());
    bounded_float_class.def(init<ValidatedNumericType>());
    bounded_float_class.def(init<FloatValue<PR>>());
    bounded_float_class.def(init<FloatBall<PR>>());
    bounded_float_class.def(init<FloatBounds<PR>>());

    bounded_float_class.def("lower", &FloatBounds<PR>::lower);
    bounded_float_class.def("upper", &FloatBounds<PR>::upper);
    bounded_float_class.def("value", &FloatBounds<PR>::value);
    bounded_float_class.def("error", &FloatBounds<PR>::error);

    bounded_float_class.def("precision", &FloatBounds<PR>::precision);

    bounded_float_class.template define_mixed_arithmetic<Int>();
    bounded_float_class.template define_mixed_arithmetic<FloatBall<PR>>();
    bounded_float_class.template define_mixed_arithmetic<FloatBounds<PR>>();
    bounded_float_class.template define_mixed_arithmetic<ValidatedNumber>();
    bounded_float_class.define_self_arithmetic();
//    bounded_float_class.define_mixed_arithmetic<ApproximateNumericType>();
//    bounded_float_class.define_mixed_arithmetic<LowerNumericType>();
//    bounded_float_class.define_mixed_arithmetic<UpperNumericType>();
//    bounded_float_class.define_mixed_arithmetic<ValidatedNumericType>();
    bounded_float_class.define_self_comparisons();

    bounded_float_class.define_unary_arithmetic();
    bounded_float_class.define_transcendental_functions();

    bounded_float_class.def(self_ns::str(self));
    bounded_float_class.def(self_ns::repr(self));

    implicitly_convertible<FloatBall<PR>,FloatBounds<PR>>();

    bounded_float_class.def("set_output_precision",&FloatBounds<PR>::set_output_precision).staticmethod("set_output_precision");

}

template<class PR> void export_upper_float()
{
    class_<FloatUpperBound<PR>> upper_float_class("Float"+class_tag<PR>()+"UpperBound");
    upper_float_class.def(init<int>());
    upper_float_class.def(init<double>());
    upper_float_class.def(init<Real,PR>());

    upper_float_class.def(init<FloatValue<PR>>());
    upper_float_class.def(init<FloatBall<PR>>());
    upper_float_class.def(init<FloatBounds<PR>>());
    upper_float_class.def(init<FloatUpperBound<PR>>());
    upper_float_class.def(init<UpperNumericType>());
    upper_float_class.def(self_ns::str(self));
    upper_float_class.def(self_ns::repr(self));

    upper_float_class.def("precision", &FloatUpperBound<PR>::precision);

    upper_float_class.def(+self);
    upper_float_class.def(-self);
    upper_float_class.def(self + self);
    upper_float_class.def(self + FloatUpperBound<PR>());
    upper_float_class.def(FloatUpperBound<PR>() + self);
    upper_float_class.def(self - self);
    upper_float_class.def(self - FloatLowerBound<PR>());
    upper_float_class.def(FloatLowerBound<PR>() - self);

    upper_float_class.def(self + Int());
    upper_float_class.def(Int() + self);
    upper_float_class.def(self - Int());
    upper_float_class.def(Int() - self);

    upper_float_class.def(self > FloatBounds<PR>());
    upper_float_class.def(self > FloatApproximation<PR>());
    upper_float_class.def(self < FloatLowerBound<PR>());
    upper_float_class.def(self >= FloatBounds<PR>());
    upper_float_class.def(self >= FloatApproximation<PR>());
    upper_float_class.def(self <= FloatLowerBound<PR>());

//    upper_float_class.define_mixed_arithmetic<ApproximateNumericType>();
//    upper_float_class.def(UpperNumericType() + self);
//    upper_float_class.def(LowerNumericType() - self);
//    upper_float_class.def(self + UpperNumericType());
//    upper_float_class.def(self - LowerNumericType());

    upper_float_class.define_monotonic_functions();

    implicitly_convertible<FloatBounds<PR>,FloatUpperBound<PR>>();
}

template<class PR> void export_lower_float()
{
    class_<FloatLowerBound<PR>> lower_float_class("Float"+class_tag<PR>()+"LowerBound");
    lower_float_class.def(init<int>());
    lower_float_class.def(init<double>());
    lower_float_class.def(init<Real,PR>());

    lower_float_class.def(init<FloatValue<PR>>());
    lower_float_class.def(init<FloatBall<PR>>());
    lower_float_class.def(init<FloatBounds<PR>>());
    lower_float_class.def(init<FloatLowerBound<PR>>());
    lower_float_class.def(init<LowerNumericType>());
    lower_float_class.def(self_ns::str(self));
    lower_float_class.def(self_ns::repr(self));

    lower_float_class.def("precision", &FloatLowerBound<PR>::precision);

    lower_float_class.def(+self);
    lower_float_class.def(-self);
    lower_float_class.def(self + self);
    lower_float_class.def(self + FloatLowerBound<PR>());
    lower_float_class.def(FloatLowerBound<PR>() + self);
    lower_float_class.def(self - self);
    lower_float_class.def(self - FloatUpperBound<PR>());
    lower_float_class.def(FloatUpperBound<PR>() - self);

    lower_float_class.def(self + Int());
    lower_float_class.def(Int() + self);
    lower_float_class.def(self - Int());
    lower_float_class.def(Int() - self);

    lower_float_class.def(self < FloatBounds<PR>());
    lower_float_class.def(self < FloatApproximation<PR>());
    lower_float_class.def(self > FloatUpperBound<PR>());
    lower_float_class.def(self <= FloatBounds<PR>());
    lower_float_class.def(self <= FloatApproximation<PR>());
    lower_float_class.def(self >= FloatUpperBound<PR>());

    //    lower_float_class.define_mixed_arithmetic<ApproximateNumericType>();

//    lower_float_class.def(LowerNumericType() + self);
//    lower_float_class.def(UpperNumericType() - self);
//    lower_float_class.def(self + LowerNumericType());
//    lower_float_class.def(self - UpperNumericType());

    lower_float_class.define_monotonic_functions();

    implicitly_convertible<FloatBounds<PR>,FloatLowerBound<PR>>();
}



template<class PR> void export_approximate_float()
{
    class_<FloatApproximation<PR>> approximate_float_class("Float"+class_tag<PR>()+"Approximation");
    approximate_float_class.def(init<double>());
    approximate_float_class.def(init<Real,PR>());

    approximate_float_class.def(init<FloatValue<PR>>());
    approximate_float_class.def(init<FloatBall<PR>>());
    approximate_float_class.def(init<FloatBounds<PR>>());
    approximate_float_class.def(init<FloatLowerBound<PR>>());
    approximate_float_class.def(init<FloatUpperBound<PR>>());
    approximate_float_class.def(init<FloatApproximation<PR>>());
    approximate_float_class.def(init<ApproximateNumericType>());

    approximate_float_class.define_self_arithmetic();
    approximate_float_class.template define_mixed_arithmetic<FloatApproximation<PR>>();
    approximate_float_class.template define_mixed_arithmetic<ApproximateNumber>();
    approximate_float_class.template define_mixed_arithmetic<double>();
    approximate_float_class.define_self_comparisons();

    approximate_float_class.def("precision", &FloatApproximation<PR>::precision);
    approximate_float_class.def("get_d", &FloatApproximation<PR>::get_d);


    approximate_float_class.def(self_ns::str(self));
    approximate_float_class.def(self_ns::repr(self));

    implicitly_convertible<double,FloatApproximation<PR>>();
    //implicitly_convertible<Integer,FloatBounds<PR>>();
    implicitly_convertible<FloatValue<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatBall<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatBounds<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatUpperBound<PR>,FloatApproximation<PR>>();
    implicitly_convertible<FloatLowerBound<PR>,FloatApproximation<PR>>();

    approximate_float_class.def(abs(self));

    approximate_float_class.define_unary_arithmetic();
    approximate_float_class.define_transcendental_functions();
}

Void export_effort() {
    class_<Effort> effort_class("Effort",init<Nat>());
    effort_class.def(init<>());
    effort_class.def(self_ns::str(self));
}

Void export_precision() {
    class_<Precision64> precision64_class("Precision64",init<>());
    precision64_class.def(self_ns::str(self));
    class_<PrecisionMP> precisionmp_class("PrecisionMP",init<Nat>());
    precisionmp_class.def("bits",&PrecisionMP::bits);
    precisionmp_class.def(self_ns::str(self));
    class_<Accuracy> accuracy_class("Accuracy",init<Nat>());
    accuracy_class.def("bits",&Accuracy::bits);
    accuracy_class.def(self_ns::str(self));
}

template<class PR> Void export_user_floats() {
    class_<RawFloat<PR>> raw_float_class("RawFloat"+class_tag<PR>(),init<double>());
    raw_float_class.def(self_ns::str(self));

    export_approximate_float<PR>();
    export_upper_float<PR>();
    export_lower_float<PR>();
    export_bounded_float<PR>();
    export_metric_float<PR>();
    export_exact_float<PR>();
    export_error_float<PR>();
}


} // namespace Ariadne;

void
numeric_submodule()
{
    using namespace Ariadne;
    export_effort();

    export_logical<ExactTag>("Boolean");
    export_effective_logical<EffectiveTag>("Kleenean");
    export_effective_logical<EffectiveUpperTag>("Sierpinskian");
    export_effective_logical<EffectiveLowerTag>("NegatedSierpinskian");
    export_logical<ValidatedTag>("Tribool");
    export_logical<UpperTag>("Verified");
    export_logical<LowerTag>("Falsified");
    export_logical<ApproximateTag>("Fuzzy");


    export_precision();
    export_user_floats<Precision64>();
    export_user_floats<PrecisionMP>();

    export_numbers();
    export_real();
    export_rational();
    export_integer();
}
