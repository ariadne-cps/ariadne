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

#include <iostream>
#include <iomanip>

#include "config.h"

#include <boost/python.hpp>

#include "utility/tribool.h"
#include "numeric/numeric.h"
#include "utility/exceptions.h"

using namespace boost::python;
using namespace Ariadne;


namespace Ariadne {

Void set_output_precision(Nat p) { std::cout << std::setprecision(p); }

Dyadic operator+(Dyadic const& x1, Dyadic const& x2) {
    Dyadic r(x1.value()+x2.value()); ARIADNE_ASSERT(Rational(r)==Rational(x1)+Rational(x2)); return r; }
Dyadic operator-(Dyadic const& x1, Dyadic const& x2) {
    Dyadic r(x1.value()-x2.value()); ARIADNE_ASSERT(Rational(r)==Rational(x1)-Rational(x2)); return r; }
Dyadic operator*(Dyadic const& x1, Dyadic const& x2) {
    Dyadic r(x1.value()*x2.value()); ARIADNE_ASSERT(Rational(r)==Rational(x1)*Rational(x2)); return r; }
Rational operator/(Dyadic const& x1, Dyadic const& x2) {
    return Rational(x1)/Rational(x2); }

class PythonRational : public Rational {
  public:
    using Rational::Rational;
    PythonRational() = default;
    PythonRational(double x) { throw std::runtime_error("Cannot construct Rational from builtin floating-point number"); }
};

template<>
struct from_python_dict<ValidatedFloat64> {
    from_python_dict() { converter::registry::push_back(&convertible,&construct,type_id<ValidatedFloat64>()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyDict_Check(obj_ptr) || len(extract<dict>(obj_ptr))!=1) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::dict dct = boost::python::extract<boost::python::dict>(obj_ptr);
        boost::python::list lst=dct.items();
        assert(boost::python::len(lst)==1);
        Void* storage = ((converter::rvalue_from_python_storage<ValidatedFloat64>*)data)->storage.bytes;
        new (storage) ValidatedFloat64(boost::python::extract<Float64>(lst[0][0]),boost::python::extract<Float64>(lst[0][1]));
        data->convertible = storage;
    }
};


template<>
struct from_python_list<ValidatedFloat64> {
    from_python_list() { converter::registry::push_back(&convertible,&construct,type_id<ValidatedFloat64>()); }
    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) || len(extract<list>(obj_ptr))!=2) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        assert(boost::python::len(lst)==2);
        Void* storage = ((converter::rvalue_from_python_storage<ValidatedFloat64>*)data)->storage.bytes;
        new (storage) ValidatedFloat64(boost::python::extract<Float64>(lst[0]),boost::python::extract<Float64>(lst[1]));
        data->convertible = storage;
    }
};

/*
struct interval_from_python_str {
    interval_from_python_str() { converter::registry::push_back(&convertible,&construct,type_id<ValidatedFloat64>()); }
    static Void* convertible(PyObject* obj_ptr) { if (!PyString_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        StringType str = boost::python::extract<StringType>(obj_ptr);
        Void* storage = ((converter::rvalue_from_python_storage<ValidatedFloat64>*)data)->storage.bytes;
        storage = new ValidatedFloat64(str);
        data->convertible = storage;
    }
};
*/


OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Float64>& x) {
    return os << std::showpoint << std::setprecision(18) << x.reference().get_d();
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ApproximateFloat64>& x) {
    return os << std::showpoint << std::setprecision(18) << x.reference().get_d();
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<UpperFloat64>& x) {
    return os << std::showpoint << std::setprecision(18) << x.reference().get_d();
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ValidatedFloat64>& x) {
    rounding_mode_t rnd=Float64::get_rounding_mode();
    os << '{';
    Float64::set_rounding_downward();
    os << std::showpoint << std::setprecision(18) << x.reference().lower().get_d();
    os << ':';
    Float64::set_rounding_upward();
    os << std::showpoint << std::setprecision(18) << x.reference().upper().get_d();
    Float64::set_rounding_mode(rnd);
    os << '}';
    return os;

}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<ExactFloat64>& x) {
    return os << std::showpoint << std::setprecision(18) << x.reference().get_d();
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Integer>& x) {
    return os << x.reference();
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Rational>& x) {
    return os << "Rational(" << x.reference().get_num() << "," << x.reference().get_den() << ")";
}

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Real>& x) {
    return os << "Real(" << x.reference() << ")";
}

} // namespace Ariadne


StringType __str__(Kleenean tb) {
  if(is_indeterminate(tb)) { return "Indeterminate"; }
  if(definitely(tb)) { return "True"; }
  else if(not possibly(tb)) { return "False"; }
  else { return "Indeterminate"; }
}

StringType __repr__(Kleenean tb) {
  if(definitely(tb)) { return "Kleenean(True)"; }
  else if(not possibly(tb)) { return "Kleenean(False)"; }
  else { return "Kleenean(Indeterminate)"; }
}

Bool __nonzero__(Kleenean tb) { return definitely(tb); }
Kleenean indeterminate_const() { return indeterminate; }

namespace Ariadne {
Bool possibly(Logical<Validated>);
Bool definitely(Logical<Validated>);
Bool is_determinate(Logical<Validated>);
}

Void export_tribool() {

    class_<Kleenean> tribool_class("Kleenean",init<Bool>());
    tribool_class.def(init<Int>());
    tribool_class.def(init<Kleenean>());
    tribool_class.def("__eq__", &__eq__<Kleenean,Kleenean,Kleenean>);
    tribool_class.def("__ne__", &__ne__<Kleenean,Kleenean,Kleenean>);
    tribool_class.def("__and__", &__and__<Kleenean,Kleenean,Kleenean>);
    tribool_class.def("__or__", &__or__<Kleenean,Kleenean,Kleenean>);
    // WARNING: __not__ is not a special method!
    tribool_class.def("__not__", &__not__<Kleenean,Kleenean>);
    tribool_class.def("__nonzero__", (Bool(*)(Kleenean))&__nonzero__);

    //tribool_class.def("__eq__", (Logical<Validated>(*)(Logical<Validated>,Bool))(&operator==));
    //tribool_class.def("__neq__", (Logical<Validated>(*)(Logical<Validated>,Bool))(&operator!=));
    //tribool_class.def("__and__", (Logical<Validated>(*)(Logical<Validated>,Bool))(&operator!=));
    //tribool_class.def("__or__", (Logical<Validated>(*)(Logical<Validated>,Bool))(&operator!=));

    tribool_class.def("__str__", (StringType(*)(Kleenean))&__str__);
    tribool_class.def("__repr__", (StringType(*)(Kleenean))&__repr__);

    implicitly_convertible<Bool,Kleenean>();

    def("indeterminate",(Kleenean(*)())&indeterminate_const);
    def("possibly",(Bool(*)(Logical<Validated>))&possibly);
    def("definitely",(Bool(*)(Logical<Validated>))&definitely);
    def("is_determinate",(Bool(*)(Logical<Validated>))&is_determinate);
    // no facility for wrapping C++ constants
    // def("Indeterminate",tribool_indeterminate_constant);

}


Void export_integer()
{
    class_<Integer> integer_class("Integer");
    integer_class.def(init<Int>());
    integer_class.def(init<Integer>());
    integer_class.def(boost::python::self_ns::str(self));
    integer_class.def("__repr__", &__repr__<Integer>);
    integer_class.def("__lt__",&__lt__<Bool,Integer,Integer>);

    integer_class.def("__pos__", &__pos__<Integer,Integer>);
    integer_class.def("__neg__", &__neg__<Integer,Integer>);
    integer_class.def("__add__", &__add__<Integer,Integer,Integer>);
    integer_class.def("__sub__", &__sub__<Integer,Integer,Integer>);
    integer_class.def("__mul__", &__mul__<Integer,Integer,Integer>);

    implicitly_convertible<Int,Integer>();
}


namespace Ariadne {
Rational sqr(Rational const&);
}

Void export_rational()
{
    class_<Rational> rational_class("Rational");
    rational_class.def(init<Integer,Integer>());
    rational_class.def(init<Int>());
    rational_class.def(init<double>());
    rational_class.def(init<StringType>());
    rational_class.def(init<Rational>());
    rational_class.def(boost::python::self_ns::str(self));
    rational_class.def("__repr__", &__repr__<Rational>);
    rational_class.def("__lt__",&__lt__<Bool,Rational,Rational>);

    rational_class.def("__pos__", &__pos__<Rational,Rational>);
    rational_class.def("__neg__", &__neg__<Rational,Rational>);
    rational_class.def("__add__", &__add__<Rational,Rational,Rational>);
    rational_class.def("__sub__", &__sub__<Rational,Rational,Rational>);
    rational_class.def("__mul__", &__mul__<Rational,Rational,Rational>);
    rational_class.def("__div__", &__div__<Rational,Rational,Rational>);

    def("sqr",(Rational(*)(Rational const&)) &sqr);

    implicitly_convertible<Integer,Rational>();

}

Void export_dyadic()
{
    class_< Dyadic > dyadic_class("Dyadic",init<Dyadic>());
    dyadic_class.def(init<>());
    dyadic_class.def(init<Int>());
    dyadic_class.def(init<double>());
    dyadic_class.def(+self);
    dyadic_class.def(-self);
    dyadic_class.def(self+self);
    dyadic_class.def(self-self);
    dyadic_class.def(self*self);
    dyadic_class.def(self/self);

    dyadic_class.def("__str__", &__cstr__<Dyadic>);
    dyadic_class.def("__repr__", &__cstr__<Dyadic>);
}

Void export_decimal()
{
    class_< Decimal > decimal_class("Decimal");
    decimal_class.def(init<double>());
    decimal_class.def(init<StringType>());
    decimal_class.def(boost::python::self_ns::str(self));
}

Real pi_function() { return pi; }

namespace Ariadne {
Real pow(Real, Int);
Real sqr(Real);
Real rec(Real);
Real sqrt(Real x);
Real exp(Real);
Real log(Real);
Real sin(Real);
Real cos(Real);
Real tan(Real);
Real atan(Real);
}

Void export_real()
{
    class_<Real> real_class("Real",init<Real>());
    real_class.def(init<Int>());
    real_class.def(init<Integer>());
    real_class.def(init<Rational>());
    real_class.def(init<Dyadic>());
    real_class.def(init<Decimal>());

    real_class.def("radius", &ValidatedFloat64::radius);
    real_class.def(boost::python::self_ns::str(self));
    real_class.def("__repr__", &__repr__<Real>);

    real_class.def(+self);
    real_class.def(-self);
    real_class.def(self + self);
    real_class.def(self - self);
    real_class.def(self * self);
    real_class.def(self / self);

    real_class.def(Int() + self);
    real_class.def(Int() - self);
    real_class.def(Int() * self);
    real_class.def(Int() / self);

    real_class.def(self == self);
    real_class.def(self != self);
    real_class.def(self >= self);
    real_class.def(self <= self);
    real_class.def(self > self);
    real_class.def(self < self);

    def("pi", (Real(*)()) &pi_function);

    def("pow",  (Real(*)(Real, Int)) &pow);
    def("sqr", (Real(*)(Real)) &sqr);
    def("rec", (Real(*)(Real)) &rec);
    def("sqrt", (Real(*)(Real)) &sqrt);
    def("exp", (Real(*)(Real)) &exp);
    def("log", (Real(*)(Real)) &log);

    def("sin", (Real(*)(Real)) &sin);
    def("cos", (Real(*)(Real)) &cos);
    def("tan", (Real(*)(Real)) &tan);
    def("atan", (Real(*)(Real)) &atan);

    implicitly_convertible<Int,Real>();
    implicitly_convertible<Integer,Real>();
    implicitly_convertible<Rational,Real>();
    implicitly_convertible<Decimal,Real>();
    implicitly_convertible<Dyadic,Real>();
}




Void export_exact_float()
{
    class_< ExactFloat64 > exact_float_class("ExactFloat64",init<ExactFloat64>());
    exact_float_class.def(init<>());
    exact_float_class.def(init<Float64>());
    exact_float_class.def(+self);
    exact_float_class.def(-self);
    exact_float_class.def(self+self);
    exact_float_class.def(self-self);
    exact_float_class.def(self*self);
    exact_float_class.def(self/self);

    exact_float_class.def("set_output_precision", &ExactFloat64::set_output_precision);
    exact_float_class.staticmethod("set_output_precision");
    exact_float_class.def("__str__", &__cstr__<ExactFloat64>);
    exact_float_class.def("__repr__", &__cstr__<ExactFloat64>);
}

Void export_validated_float()
{
    using boost::python::class_;
    using boost::python::init;
    using boost::python::self;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::def;

    def("down",&down);
    def("up",&up);

    class_< ValidatedFloat64 > validated_float_class("ValidatedFloat64");
    validated_float_class.def(init<double>());
    validated_float_class.def(init<double,double>());
    validated_float_class.def(init<Float64,Float64>());
    validated_float_class.def(init<Real>());
    validated_float_class.def(init<Decimal>());
    validated_float_class.def(init<Dyadic>());
    validated_float_class.def(init<ValidatedFloat64>());
    validated_float_class.def(init<Float64>());
    validated_float_class.def(init<Rational>());

    validated_float_class.def(+self);
    validated_float_class.def(-self);
    validated_float_class.def(self + self);
    validated_float_class.def(self - self);
    validated_float_class.def(self * self);
    validated_float_class.def(self / self);

    validated_float_class.def(Real() + self);
    validated_float_class.def(Real() - self);
    validated_float_class.def(Real() * self);
    validated_float_class.def(Real() / self);

    validated_float_class.def(self == self);
    validated_float_class.def(self != self);
    validated_float_class.def(self >= self);
    validated_float_class.def(self <= self);
    validated_float_class.def(self > self);
    validated_float_class.def(self < self);
    validated_float_class.def("lower", &ValidatedFloat64::lower);
    validated_float_class.def("upper", &ValidatedFloat64::upper);
    validated_float_class.def("lower_value", &ValidatedFloat64::lower_value, return_value_policy<copy_const_reference>());
    validated_float_class.def("upper_value", &ValidatedFloat64::upper_value, return_value_policy<copy_const_reference>());
    validated_float_class.def("midpoint", &ValidatedFloat64::midpoint);
    validated_float_class.def("radius", &ValidatedFloat64::radius);
    validated_float_class.def("width", &ValidatedFloat64::width);
    validated_float_class.def(boost::python::self_ns::str(self));
    validated_float_class.def("__repr__", &__repr__<ValidatedFloat64>);

    validated_float_class.def("set_output_precision", &ValidatedFloat64::set_output_precision);
    validated_float_class.staticmethod("set_output_precision");

    implicitly_convertible<Real,ValidatedFloat64>();
    implicitly_convertible<ExactFloat64,ValidatedFloat64>();

    from_python_dict<ValidatedFloat64>();
    from_python_list<ValidatedFloat64>();
    //from_python_str<ValidatedFloat64>();

    def("mag", (UpperFloat64(*)(ValidatedFloat64)) &mag);

    def("max", (ValidatedFloat64(*)(ValidatedFloat64,ValidatedFloat64)) &max);
    def("min", (ValidatedFloat64(*)(ValidatedFloat64,ValidatedFloat64)) &min);

    def("trunc", (ValidatedFloat64(*)(ValidatedFloat64,Nat)) &trunc, "truncate to n binary digits");
    def("abs", (ValidatedFloat64(*)(ValidatedFloat64)) &abs, "validated_float absolute value function");
    def("pow",  (ValidatedFloat64(*)(ValidatedFloat64,Int)) &pow, "validated_float power function");
    def("sqr", (ValidatedFloat64(*)(ValidatedFloat64)) &sqr, "validated_float square function");
    def("rec", (ValidatedFloat64(*)(ValidatedFloat64)) &rec);
    def("sqrt", (ValidatedFloat64(*)(ValidatedFloat64)) &sqrt);
    def("exp", (ValidatedFloat64(*)(ValidatedFloat64)) &exp);
    def("log", (ValidatedFloat64(*)(ValidatedFloat64)) &log);

    def("sin", (ValidatedFloat64(*)(ValidatedFloat64)) &sin);
    def("cos", (ValidatedFloat64(*)(ValidatedFloat64)) &cos);
    def("tan", (ValidatedFloat64(*)(ValidatedFloat64)) &tan);
    def("asin", (ValidatedFloat64(*)(ValidatedFloat64)) &asin);
    def("acos", (ValidatedFloat64(*)(ValidatedFloat64)) &acos);
    def("atan", (ValidatedFloat64(*)(ValidatedFloat64)) &atan);
}

Void export_upper_float()
{
    class_< UpperFloat64 > upper_float_class("UpperFloat64");
    upper_float_class.def(init<double>());
    upper_float_class.def(init<ValidatedFloat64>());
    upper_float_class.def(init<ExactFloat64>());

    upper_float_class.def(-self);
    upper_float_class.def(self + self);
    upper_float_class.def(self - LowerFloat64());

    upper_float_class.def(boost::python::self_ns::str(self));
    implicitly_convertible<ValidatedFloat64,UpperFloat64>();
}

Void export_lower_float()
{
    class_< LowerFloat64 > lower_float_class("LowerFloat64");
    lower_float_class.def(init<double>());
    lower_float_class.def(init<Float64>());
    lower_float_class.def(init<ValidatedFloat64>());
    lower_float_class.def(init<ExactFloat64>());

    lower_float_class.def(-self);
    lower_float_class.def(self + self);
    lower_float_class.def(self - UpperFloat64());

    lower_float_class.def(boost::python::self_ns::str(self));
    implicitly_convertible<ValidatedFloat64,LowerFloat64>();
}

Void export_approximate_float()
{
    using boost::python::class_;
    using boost::python::init;
    using boost::python::self;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::def;

    def("down",&down);
    def("up",&up);

    class_< ApproximateFloat64 > approximate_float_class("ApproximateFloat64");
    approximate_float_class.def(init<double>());
    approximate_float_class.def(init<Float64>());
    approximate_float_class.def(init<Decimal>());
    approximate_float_class.def(init<Dyadic>());
    approximate_float_class.def(init<ApproximateFloat64>());
    approximate_float_class.def(init<LowerFloat64>());
    approximate_float_class.def(init<UpperFloat64>());
    approximate_float_class.def(init<ValidatedFloat64>());
    approximate_float_class.def(init<Rational>());

    approximate_float_class.def(+self);
    approximate_float_class.def(-self);
    approximate_float_class.def(self + self);
    approximate_float_class.def(self - self);
    approximate_float_class.def(self * self);
    approximate_float_class.def(self / self);

    approximate_float_class.def(Real() + self);
    approximate_float_class.def(Real() - self);
    approximate_float_class.def(Real() * self);
    approximate_float_class.def(Real() / self);

    approximate_float_class.def(self == self);
    approximate_float_class.def(self != self);
    approximate_float_class.def(self >= self);
    approximate_float_class.def(self <= self);
    approximate_float_class.def(self > self);
    approximate_float_class.def(self < self);
    approximate_float_class.def(boost::python::self_ns::str(self));
    approximate_float_class.def("__repr__", &__repr__<ApproximateFloat64>);

    approximate_float_class.def("set_output_precision", &ApproximateFloat64::set_output_precision);
    approximate_float_class.staticmethod("set_output_precision");

    implicitly_convertible<LowerFloat64,ApproximateFloat64>();
    implicitly_convertible<UpperFloat64,ApproximateFloat64>();
    implicitly_convertible<ValidatedFloat64,ApproximateFloat64>();

    def("max", (ApproximateFloat64(*)(ApproximateFloat64,ApproximateFloat64)) &max);
    def("min", (ApproximateFloat64(*)(ApproximateFloat64,ApproximateFloat64)) &min);

    def("abs", (ApproximateFloat64(*)(ApproximateFloat64)) &abs);
    def("pow",  (ApproximateFloat64(*)(ApproximateFloat64,Int)) &pow);
    def("sqr", (ApproximateFloat64(*)(ApproximateFloat64)) &sqr);
    def("rec", (ApproximateFloat64(*)(ApproximateFloat64)) &rec);
    def("sqrt", (ApproximateFloat64(*)(ApproximateFloat64)) &sqrt);
    def("exp", (ApproximateFloat64(*)(ApproximateFloat64)) &exp);
    def("log", (ApproximateFloat64(*)(ApproximateFloat64)) &log);

    def("sin", (ApproximateFloat64(*)(ApproximateFloat64)) &sin);
    def("cos", (ApproximateFloat64(*)(ApproximateFloat64)) &cos);
    def("tan", (ApproximateFloat64(*)(ApproximateFloat64)) &tan);
    def("asin", (ApproximateFloat64(*)(ApproximateFloat64)) &asin);
    def("acos", (ApproximateFloat64(*)(ApproximateFloat64)) &acos);
    def("atan", (ApproximateFloat64(*)(ApproximateFloat64)) &atan);
}

Void export_float()
{
    class_< Float64 > float_class("Float64");
    float_class.def(init<double>());
    float_class.def(init<Float64>());
    float_class.def(boost::python::self_ns::str(self));
    float_class.def("__repr__", &__repr__<Float64>);

    float_class.def(+self);
    float_class.def(-self);
    float_class.def(self + self);
    float_class.def(self - self);
    float_class.def(self * self);
    float_class.def(self / self);

    float_class.def(double() + self);
    float_class.def(double() - self);
    float_class.def(double() * self);
    float_class.def(double() / self);

    float_class.def(self == self);
    float_class.def(self != self);
    float_class.def(self >= self);
    float_class.def(self <= self);
    float_class.def(self > self);
    float_class.def(self < self);

    implicitly_convertible<double,Float64>();

    def("min",(Float64(*)(Float64,Float64)) &Ariadne::min);
    def("max",(Float64(*)(Float64,Float64)) &Ariadne::max);
    def("abs", (Float64(*)(Float64)) &Ariadne::abs);

    def("pow",(Float64(*)(Float64,Int)) &Ariadne::pow);

    def("rec",(Float64(*)(Float64)) &Ariadne::rec);
    def("sqr",(Float64(*)(Float64)) &Ariadne::sqr);
    def("sqrt", (Float64(*)(Float64)) &Ariadne::sqrt);

    def("exp", (Float64(*)(Float64)) &Ariadne::exp);
    def("log", (Float64(*)(Float64)) &Ariadne::log);

    def("sin", (Float64(*)(Float64)) &Ariadne::sin);
    def("cos", (Float64(*)(Float64)) &Ariadne::cos);
    def("tan", (Float64(*)(Float64)) &Ariadne::tan);
    def("asin", (Float64(*)(Float64)) &Ariadne::asin);
    def("acos", (Float64(*)(Float64)) &Ariadne::acos);
    def("atan", (Float64(*)(Float64)) &Ariadne::atan);

}


Void
numeric_submodule()
{
    export_tribool();
    export_float();
    export_approximate_float();
    export_lower_float();
    export_upper_float();
    export_validated_float();
    export_exact_float();
    export_real();
    export_rational();
    export_decimal();
    export_dyadic();
    export_integer();
}
