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
template<class PR> using ApproximateFloat=Float<Approximate,PR>;
template<class PR> using LowerFloat=Float<Lower,PR>;
template<class PR> using UpperFloat=Float<Upper,PR>;
template<class PR> using BoundedFloat=Float<Bounded,PR>;
template<class PR> using MetricFloat=Float<Metric,PR>;
template<class PR> using ExactFloat=Float<Exact,PR>;
template<class PR> using ErrorFloat=Float<Error,PR>;

// Declare friend functions in namespace
Rational operator/(Integer const&, Integer const&);
Integer pow(Integer const& z, Nat m);
//Integer abs(Integer const& z);

DECLARE_NUMERIC_OPERATIONS(Real,PositiveReal);
//DECLARE_NUMERIC_OPERATIONS(MetricFloat<PR>);
//DECLARE_NUMERIC_OPERATIONS(BoundedFloat<PR>);
//DECLARE_NUMERIC_OPERATIONS(ApproximateFloat<PR>64);

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

template<class T> struct Tag { };

template<class T> class class_ : public boost::python::class_<T> {
    typedef T const& Tcr;
  public:
    using boost::python::class_<T>::class_;
    class_(String const& name) : boost::python::class_<T>(name.c_str()) { }
    void define_self_arithmetic() {
        this->def(+self); this->def(-self);
        this->def(self+self); this->def(self-self); this->def(self*self); this->def(self/self);
    }
    template<class X> void define_mixed_arithmetic(Tag<X> tag = Tag<X>()) {
        X const* other_ptr=nullptr; X const& other=*other_ptr;
        this->def(self+other); this->def(self-other); this->def(self*other); this->def(self/other);
        this->def(other+self); this->def(other-self); this->def(other*self); this->def(other/self);
    }
    void define_unary_arithmetic() {
        this->def(+self); this->def(-self);
        using boost::python::def;
        def("pos",(T(*)(Tcr))&pos); def("neg",(T(*)(Tcr))&neg);
        def("sqr",(T(*)(Tcr))&sqr); def("rec",(T(*)(Tcr))&rec);
    }
    void define_transcendental_functions();
    void define_monotonic_functions() {
        typedef decltype(-declval<T>()) NT;
        using boost::python::def;
        def("pos",(T(*)(Tcr))&pos); def("neg",(NT(*)(Tcr))&neg);
        def("sqrt",(T(*)(Tcr))&sqrt); def("exp",(T(*)(Tcr))&exp); def("log",(T(*)(Tcr))&log); def("atan",(T(*)(Tcr))&atan);
    }
};

template<> void class_<Real>::define_transcendental_functions() {
    typedef Real T;
    using boost::python::def;
    T sqrt(T); T exp(T); T log(T); T sin(T); T cos(T); T tan(T); T atan(T);
    def("sqrt",(T(*)(T))&sqrt); def("exp",(T(*)(T))&exp); def("log",(T(*)(T))&log);
    def("sin",(T(*)(T))&sin); def("cos",(T(*)(T))&cos); def("tan",(T(*)(T))&tan); def("atan",(T(*)(T))&atan);
}

template<class T> void class_<T>::define_transcendental_functions() {
    using boost::python::def;
    def("sqrt",(T(*)(Tcr))&sqrt); def("exp",(T(*)(Tcr))&exp); def("log",(T(*)(Tcr))&log);
    def("sin",(T(*)(Tcr))&sin); def("cos",(T(*)(Tcr))&cos); def("tan",(T(*)(Tcr))&tan); def("atan",(T(*)(Tcr))&atan);
}

void export_integer()
{
    class_<Integer> integer_class("Integer");
    integer_class.def(init<int>());
    integer_class.def(init<Integer>());

    integer_class.def(boost::python::self_ns::str(self));
    integer_class.def(boost::python::self_ns::repr(self));

//    integer_class.define_self_arithmetic();

    // Required for mixed operations with integral types
    integer_class.define_mixed_arithmetic<Integer>();

    def("pow", (Integer(*)(const Integer&,Nat)) &pow) ;
    integer_class.def(abs(self));

    integer_class.def(self == self);
    integer_class.def(self < self);

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
    rational_class.define_mixed_arithmetic<double>();
    rational_class.define_mixed_arithmetic<Rational>();

    rational_class.def(boost::python::self_ns::str(self));
    rational_class.def(boost::python::self_ns::repr(self));
    rational_class.def(self < self);

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
//    real_class.define_mixed_arithmetic<double>();
    // Needed to dispatch Real # ExactNumericType operations
    //real_class.define_mixed_arithmetic<EffectiveNumericType>();
//    real_class.define_mixed_arithmetic<Real>();
    real_class.define_transcendental_functions();

    real_class.def(boost::python::self_ns::str(self));
    real_class.def(boost::python::self_ns::repr(self));
    real_class.def(self < self);

    real_class.def("get", (BoundedFloat64(Real::*)(Precision64)const) &Real::get);
    real_class.def("get", (BoundedFloatMP(Real::*)(PrecisionMP)const) &Real::get);
    real_class.def("get", (BoundedFloatMP(Real::*)(Accuracy)const) &Real::evaluate);
    real_class.def("evaluate", (BoundedFloatMP(Real::*)(Accuracy)const) &Real::evaluate);
    real_class.def("get_d", &Real::get_d);

    implicitly_convertible<Rational,Real>();
}




template<class PR> void export_exact_float()
{
    class_<ExactFloat<PR>> exact_float_class("ExactFloat"+class_tag<PR>());
    exact_float_class.def(init<int>());
    exact_float_class.def(init<double>());
    exact_float_class.def(init<Integer,PR>());
    exact_float_class.def(init<ExactFloat<PR>>());

    exact_float_class.define_mixed_arithmetic(Tag<Int>());
    exact_float_class.define_mixed_arithmetic(Tag<ExactFloat<PR>>());
    exact_float_class.define_self_arithmetic();
    //exact_float_class.define_mixed_arithmetic<ApproximateNumericType>();
    //exact_float_class.define_mixed_arithmetic<LowerNumericType>();
    //exact_float_class.define_mixed_arithmetic<UpperNumericType>();
    //exact_float_class.define_mixed_arithmetic<ValidatedNumericType>();

    def("pos", (ExactFloat<PR>(*)(ExactFloat<PR> const&)) &pos);
    def("neg", (ExactFloat<PR>(*)(ExactFloat<PR> const&)) &neg);

    exact_float_class.def("precision", &ExactFloat<PR>::precision);
    exact_float_class.def("get_d",&ExactFloat<PR>::get_d);

    exact_float_class.def(boost::python::self_ns::str(self));
    exact_float_class.def(boost::python::self_ns::repr(self));

    implicitly_convertible<int,ExactFloat<PR>>();

    exact_float_class.def("set_output_precision",&ExactFloat<PR>::set_output_precision).staticmethod("set_output_precision");
}

template<class PR> void export_error_float()
{
    class_<ErrorFloat<PR>> error_float_class("ErrorFloat"+class_tag<PR>());
    error_float_class.def(init<uint>());
    error_float_class.def(init<double>());
    error_float_class.def(init<ErrorFloat<PR>>());
    error_float_class.def(+self);
    error_float_class.def(self+self);
    error_float_class.def(self*self);
    error_float_class.def("get_d",&ErrorFloat<PR>::get_d);
    error_float_class.def(boost::python::self_ns::str(self));
    error_float_class.def(boost::python::self_ns::repr(self));

    error_float_class.def("precision", &ErrorFloat<PR>::precision);

    //    error_float_class.def("set_output_precision",&ErrorFloat<PR>::set_output_precision).staticmethod("set_output_precision");
}

template<class PR> void export_metric_float()
{
    class_<MetricFloat<PR>> metric_float_class("MetricFloat"+class_tag<PR>());
    metric_float_class.def(init<double>());
    metric_float_class.def(init<Real,PR>());
    metric_float_class.def(init<double,double>());
    metric_float_class.def(init<ExactFloat<PR>,ErrorFloat<PR>>());
    metric_float_class.def(init<ValidatedNumericType>());
    metric_float_class.def(init<ExactFloat<PR>>());
    metric_float_class.def(init<MetricFloat<PR>>());
    metric_float_class.def(init<BoundedFloat<PR>>());

    metric_float_class.def("value", &MetricFloat<PR>::value);
    metric_float_class.def("error", &MetricFloat<PR>::error);
    metric_float_class.def("lower", &MetricFloat<PR>::lower);
    metric_float_class.def("upper", &MetricFloat<PR>::upper);

    metric_float_class.def("precision", &MetricFloat<PR>::precision);

    metric_float_class.define_mixed_arithmetic( Tag<Int>() );
    metric_float_class.define_mixed_arithmetic( Tag<MetricFloat<PR>>() );
    metric_float_class.define_self_arithmetic();
//    metric_float_class.define_mixed_arithmetic<MetricFloat<PR>>();
//    metric_float_class.define_mixed_arithmetic<ApproximateNumericType>();
//    metric_float_class.define_mixed_arithmetic<LowerNumericType>();
//    metric_float_class.define_mixed_arithmetic<UpperNumericType>();
//    metric_float_class.define_mixed_arithmetic<ValidatedNumericType>();

    metric_float_class.define_unary_arithmetic();
    metric_float_class.define_transcendental_functions();

    metric_float_class.def(boost::python::self_ns::str(self));
    metric_float_class.def(boost::python::self_ns::repr(self));

    implicitly_convertible<ExactFloat<PR>,MetricFloat<PR>>();
}


template<class PR> void export_bounded_float()
{
    class_<BoundedFloat<PR>> bounded_float_class("BoundedFloat"+class_tag<PR>());
    bounded_float_class.def(init<double>());
    bounded_float_class.def(init<Real,PR>());
    bounded_float_class.def(init<double,double>());
    bounded_float_class.def(init<LowerFloat<PR>,UpperFloat<PR>>());
    bounded_float_class.def(init<ValidatedNumericType>());
    bounded_float_class.def(init<ExactFloat<PR>>());
    bounded_float_class.def(init<MetricFloat<PR>>());
    bounded_float_class.def(init<BoundedFloat<PR>>());

    bounded_float_class.def("lower", &BoundedFloat<PR>::lower);
    bounded_float_class.def("upper", &BoundedFloat<PR>::upper);
    bounded_float_class.def("value", &BoundedFloat<PR>::value);
    bounded_float_class.def("error", &BoundedFloat<PR>::error);

    bounded_float_class.def("precision", &BoundedFloat<PR>::precision);

    bounded_float_class.define_mixed_arithmetic( Tag<Int>() );
    bounded_float_class.define_mixed_arithmetic( Tag<MetricFloat<PR>>() );
    bounded_float_class.define_mixed_arithmetic( Tag<BoundedFloat<PR>>() );
    bounded_float_class.define_self_arithmetic();
//    bounded_float_class.define_mixed_arithmetic<ApproximateNumericType>();
//    bounded_float_class.define_mixed_arithmetic<LowerNumericType>();
//    bounded_float_class.define_mixed_arithmetic<UpperNumericType>();
//    bounded_float_class.define_mixed_arithmetic<ValidatedNumericType>();

    bounded_float_class.define_unary_arithmetic();
    bounded_float_class.define_transcendental_functions();

    bounded_float_class.def(boost::python::self_ns::str(self));
    bounded_float_class.def(boost::python::self_ns::repr(self));

    implicitly_convertible<MetricFloat<PR>,BoundedFloat<PR>>();

    bounded_float_class.def("set_output_precision",&BoundedFloat<PR>::set_output_precision).staticmethod("set_output_precision");

}

template<class PR> void export_upper_float()
{
    class_<UpperFloat<PR>> upper_float_class("UpperFloat"+class_tag<PR>());
    upper_float_class.def(init<int>());
    upper_float_class.def(init<double>());
    upper_float_class.def(init<Real,PR>());
    upper_float_class.def(init<ExactFloat<PR>>());
    upper_float_class.def(init<MetricFloat<PR>>());
    upper_float_class.def(init<BoundedFloat<PR>>());
    upper_float_class.def(init<UpperFloat<PR>>());
    upper_float_class.def(init<UpperNumericType>());
    upper_float_class.def(boost::python::self_ns::str(self));
    upper_float_class.def(boost::python::self_ns::repr(self));

    upper_float_class.def("precision", &UpperFloat<PR>::precision);

    upper_float_class.def(+self);
    upper_float_class.def(-self);
    upper_float_class.def(self + self);
    upper_float_class.def(self + UpperFloat<PR>());
    upper_float_class.def(UpperFloat<PR>() + self);
    upper_float_class.def(self - self);
    upper_float_class.def(self - LowerFloat<PR>());
    upper_float_class.def(LowerFloat<PR>() - self);

    upper_float_class.def(self + Int());
    upper_float_class.def(Int() + self);
    upper_float_class.def(self - Int());
    upper_float_class.def(Int() - self);

//    upper_float_class.define_mixed_arithmetic<ApproximateNumericType>();
//    upper_float_class.def(UpperNumericType() + self);
//    upper_float_class.def(LowerNumericType() - self);
//    upper_float_class.def(self + UpperNumericType());
//    upper_float_class.def(self - LowerNumericType());

    upper_float_class.define_monotonic_functions();

    implicitly_convertible<BoundedFloat<PR>,UpperFloat<PR>>();
}

template<class PR> void export_lower_float()
{
    class_<LowerFloat<PR>> lower_float_class("LowerFloat"+class_tag<PR>());
    lower_float_class.def(init<int>());
    lower_float_class.def(init<double>());
    lower_float_class.def(init<Real,PR>());
    lower_float_class.def(init<ExactFloat<PR>>());
    lower_float_class.def(init<MetricFloat<PR>>());
    lower_float_class.def(init<BoundedFloat<PR>>());
    lower_float_class.def(init<LowerFloat<PR>>());
    lower_float_class.def(init<LowerNumericType>());
    lower_float_class.def(boost::python::self_ns::str(self));
    lower_float_class.def(boost::python::self_ns::repr(self));

    lower_float_class.def("precision", &LowerFloat<PR>::precision);

    lower_float_class.def(+self);
    lower_float_class.def(-self);
    lower_float_class.def(self + self);
    lower_float_class.def(self + LowerFloat<PR>());
    lower_float_class.def(LowerFloat<PR>() + self);
    lower_float_class.def(self - self);
    lower_float_class.def(self - UpperFloat<PR>());
    lower_float_class.def(UpperFloat<PR>() - self);

    lower_float_class.def(self + Int());
    lower_float_class.def(Int() + self);
    lower_float_class.def(self - Int());
    lower_float_class.def(Int() - self);

//    lower_float_class.define_mixed_arithmetic<ApproximateNumericType>();

//    lower_float_class.def(LowerNumericType() + self);
//    lower_float_class.def(UpperNumericType() - self);
//    lower_float_class.def(self + LowerNumericType());
//    lower_float_class.def(self - UpperNumericType());

    lower_float_class.define_monotonic_functions();

    implicitly_convertible<BoundedFloat<PR>,LowerFloat<PR>>();
}



template<class PR> void export_approximate_float()
{
    class_<ApproximateFloat<PR>> approximate_float_class("ApproximateFloat"+class_tag<PR>());
    approximate_float_class.def(init<double>());
    approximate_float_class.def(init<Real,PR>());
    approximate_float_class.def(init<ExactFloat<PR>>());
    approximate_float_class.def(init<MetricFloat<PR>>());
    approximate_float_class.def(init<BoundedFloat<PR>>());
    approximate_float_class.def(init<LowerFloat<PR>>());
    approximate_float_class.def(init<UpperFloat<PR>>());
    approximate_float_class.def(init<ApproximateFloat<PR>>());
    approximate_float_class.def(init<ApproximateNumericType>());

    approximate_float_class.define_self_arithmetic();
    approximate_float_class.define_mixed_arithmetic(Tag<ApproximateFloat<PR>>());
    approximate_float_class.define_mixed_arithmetic(Tag<double>());
//
    approximate_float_class.def("precision", &ApproximateFloat<PR>::precision);
    approximate_float_class.def("get_d", &ApproximateFloat<PR>::get_d);

    approximate_float_class.def(boost::python::self_ns::str(self));
    approximate_float_class.def(boost::python::self_ns::repr(self));

    implicitly_convertible<double,ApproximateFloat<PR>>();
    //implicitly_convertible<Integer,BoundedFloat<PR>>();
    implicitly_convertible<ExactFloat<PR>,ApproximateFloat<PR>>();
    implicitly_convertible<MetricFloat<PR>,ApproximateFloat<PR>>();
    implicitly_convertible<BoundedFloat<PR>,ApproximateFloat<PR>>();
    implicitly_convertible<UpperFloat<PR>,ApproximateFloat<PR>>();
    implicitly_convertible<LowerFloat<PR>,ApproximateFloat<PR>>();

    approximate_float_class.def(abs(self));

    approximate_float_class.define_unary_arithmetic();
    approximate_float_class.define_transcendental_functions();
}

Void export_precision() {
    class_<Precision64> precision64_class("Precision64",init<>());
    precision64_class.def(boost::python::self_ns::str(self));
    class_<PrecisionMP> precisionmp_class("PrecisionMP",init<Nat>());
    precisionmp_class.def("bits",&PrecisionMP::bits);
    precisionmp_class.def(boost::python::self_ns::str(self));
    class_<Accuracy> accuracy_class("Accuracy",init<Nat>());
    accuracy_class.def("bits",&Accuracy::bits);
    accuracy_class.def(boost::python::self_ns::str(self));
}

template<class PR> Void export_user_floats() {
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
    export_integer();
    export_rational();
    export_real();

    export_precision();
    export_user_floats<Precision64>();
    export_user_floats<PrecisionMP>();
}
