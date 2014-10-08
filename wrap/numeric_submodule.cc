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

bool definitely(bool b) { return b; }
bool possibly(bool b) { return b; }

void set_output_precision(uint p) { std::cout << std::setprecision(p); }

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
struct from_python_dict<ValidatedFloat> {
    from_python_dict() { converter::registry::push_back(&convertible,&construct,type_id<ValidatedFloat>()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyDict_Check(obj_ptr) || len(extract<dict>(obj_ptr))!=1) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::dict dct = boost::python::extract<boost::python::dict>(obj_ptr);
        boost::python::list lst=dct.items();
        assert(boost::python::len(lst)==1);
        void* storage = ((converter::rvalue_from_python_storage<ValidatedFloat>*)data)->storage.bytes;
        new (storage) ValidatedFloat(boost::python::extract<Float>(lst[0][0]),boost::python::extract<Float>(lst[0][1]));
        data->convertible = storage;
    }
};


template<>
struct from_python_list<ValidatedFloat> {
    from_python_list() { converter::registry::push_back(&convertible,&construct,type_id<ValidatedFloat>()); }
    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr) || len(extract<list>(obj_ptr))!=2) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::list lst = boost::python::extract<boost::python::list>(obj_ptr);
        assert(boost::python::len(lst)==2);
        void* storage = ((converter::rvalue_from_python_storage<ValidatedFloat>*)data)->storage.bytes;
        new (storage) ValidatedFloat(boost::python::extract<Float>(lst[0]),boost::python::extract<Float>(lst[1]));
        data->convertible = storage;
    }
};

/*
struct interval_from_python_str {
    interval_from_python_str() { converter::registry::push_back(&convertible,&construct,type_id<ValidatedFloat>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyString_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        std::string str = boost::python::extract<std::string>(obj_ptr);
        void* storage = ((converter::rvalue_from_python_storage<ValidatedFloat>*)data)->storage.bytes;
        storage = new ValidatedFloat(str);
        data->convertible = storage;
    }
};
*/


std::ostream& operator<<(std::ostream& os, const PythonRepresentation<Float>& x) {
    return os << std::showpoint << std::setprecision(18) << x.reference().get_d();
}

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<ApproximateFloat>& x) {
    return os << std::showpoint << std::setprecision(18) << x.reference().get_d();
}

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<UpperFloat>& x) {
    return os << std::showpoint << std::setprecision(18) << x.reference().get_d();
}

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<ValidatedFloat>& x) {
    rounding_mode_t rnd=get_rounding_mode();
    os << '{';
    set_rounding_downward();
    os << std::showpoint << std::setprecision(18) << x.reference().lower().get_d();
    os << ':';
    set_rounding_upward();
    os << std::showpoint << std::setprecision(18) << x.reference().upper().get_d();
    set_rounding_mode(rnd);
    os << '}';
    return os;

}

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<ExactFloat>& x) {
    return os << std::showpoint << std::setprecision(18) << x.reference().get_d();
}

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<Integer>& x) {
    return os << x.reference();
}

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<Rational>& x) {
    return os << "Rational(" << x.reference().get_num() << "," << x.reference().get_den() << ")";
}

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<Real>& x) {
    return os << "Real(" << x.reference() << ")";
}

} // namespace Ariadne


std::string __str__(tribool tb) {
  if(indeterminate(tb)) { return "Indeterminate"; }
  if(tb==true) { return "True"; }
  else if(tb==false) { return "False"; }
  else { return "Indeterminate"; }
}

std::string __repr__(tribool tb) {
  if(tb==true) { return "tribool(True)"; }
  else if(tb==false) { return "tribool(False)"; }
  else { return "tribool(Indeterminate)"; }
}

bool __nonzero__(tribool tb) { return tb==true; }
tribool indeterminate_const() { return indeterminate; }

void export_tribool() {

    class_<tribool> tribool_class("tribool",init<bool>());
    tribool_class.def(init<int>());
    tribool_class.def(init<tribool>());
    tribool_class.def("__eq__", (tribool(*)(tribool,tribool))(&boost::logic::operator==));
    tribool_class.def("__eq__", (tribool(*)(tribool,bool))(&boost::logic::operator==));
    tribool_class.def("__neq__", (tribool(*)(tribool,tribool))(&boost::logic::operator!=));
    tribool_class.def("__neq__", (tribool(*)(tribool,bool))(&boost::logic::operator!=));
    tribool_class.def("__and__", (tribool(*)(tribool,tribool))(&boost::logic::operator!=));
    tribool_class.def("__and__", (tribool(*)(tribool,bool))(&boost::logic::operator!=));
    tribool_class.def("__or__", (tribool(*)(tribool,tribool))(&boost::logic::operator!=));
    tribool_class.def("__or__", (tribool(*)(tribool,bool))(&boost::logic::operator!=));
    // WARNING: __not__ is not a special method!
    tribool_class.def("__not__", (tribool(*)(tribool))(&boost::logic::operator!));
    tribool_class.def("__nonzero__", (bool(*)(tribool))&__nonzero__);
    tribool_class.def("__str__", (std::string(*)(tribool))&__str__);
    tribool_class.def("__repr__", (std::string(*)(tribool))&__repr__);

    implicitly_convertible<bool,tribool>();

    def("indeterminate",(tribool(*)(void))&indeterminate_const);
    def("possibly",(bool(*)(bool))&possibly);
    def("possibly",(bool(*)(tribool))&possibly);
    def("definitely",(bool(*)(bool))&definitely);
    def("definitely",(bool(*)(tribool))&definitely);
    // no facility for wrapping C++ constants
    // def("Indeterminate",tribool_indeterminate_constant);

}


#ifdef HAVE_GMPXX_H
void export_integer()
{
    class_<Integer> integer_class("Integer");
    integer_class.def(init<int>());
    integer_class.def(init<Integer>());
    integer_class.def(boost::python::self_ns::str(self));
    integer_class.def("__repr__", &__repr__<Integer>);
    integer_class.def("__less__",(bool(*)(const mpz_class&, const mpz_class&)) &operator<);

    integer_class.def("__pos__", &__pos__<Integer,Integer>);
    integer_class.def("__neg__", &__neg__<Integer,Integer>);
    integer_class.def("__add__", &__add__<Integer,Integer,Integer>);
    integer_class.def("__sub__", &__sub__<Integer,Integer,Integer>);
    integer_class.def("__mul__", &__mul__<Integer,Integer,Integer>);

    implicitly_convertible<int,Integer>();
}
#endif

#ifdef HAVE_GMPXX_H
void export_rational()
{
    class_<Rational> rational_class("Rational");
    rational_class.def(init<Integer,Integer>());
    rational_class.def(init<int>());
    rational_class.def(init<double>());
    rational_class.def(init<std::string>());
    rational_class.def(init<Rational>());
    rational_class.def(boost::python::self_ns::str(self));
    rational_class.def("__repr__", &__repr__<Rational>);
    rational_class.def("__less__",(bool(*)(const mpq_class&, const mpq_class&)) &operator<);

    rational_class.def("__pos__", &__pos__<Rational,Rational>);
    rational_class.def("__neg__", &__neg__<Rational,Rational>);
    rational_class.def("__add__", &__add__<Rational,Rational,Rational>);
    rational_class.def("__sub__", &__sub__<Rational,Rational,Rational>);
    rational_class.def("__mul__", &__mul__<Rational,Rational,Rational>);
    rational_class.def("__div__", &__div__<Rational,Rational,Rational>);

    def("sqr",(Rational(*)(Rational const&)) &sqr);

    implicitly_convertible<Integer,Rational>();

}
#endif

void export_dyadic()
{
    class_< Dyadic > dyadic_class("Dyadic",init<Dyadic>());
    dyadic_class.def(init<>());
    dyadic_class.def(init<int>());
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

void export_decimal()
{
    class_< Decimal > decimal_class("Decimal");
    decimal_class.def(init<double>());
    decimal_class.def(init<std::string>());
    decimal_class.def(boost::python::self_ns::str(self));
}

Real pi_function() { return pi; }

void export_real()
{
    class_<Real> real_class("Real",init<Real>());
    real_class.def(init<int>());
#ifdef HAVE_GMPXX_H
    real_class.def(init<Integer>());
    real_class.def(init<Rational>());
#endif
    real_class.def(init<Dyadic>());
    real_class.def(init<Decimal>());

    real_class.def("radius", &ValidatedFloat::radius);
    real_class.def(boost::python::self_ns::str(self));
    real_class.def("__repr__", &__repr__<Real>);

    real_class.def(+self);
    real_class.def(-self);
    real_class.def(self + self);
    real_class.def(self - self);
    real_class.def(self * self);
    real_class.def(self / self);

    real_class.def(int() + self);
    real_class.def(int() - self);
    real_class.def(int() * self);
    real_class.def(int() / self);

    real_class.def(self == self);
    real_class.def(self != self);
    real_class.def(self >= self);
    real_class.def(self <= self);
    real_class.def(self > self);
    real_class.def(self < self);

    def("pi", (Real(*)()) &pi_function);

    def("pow",  (Real(*)(Real const&, int)) &pow);
    def("sqr", (Real(*)(Real const&)) &sqr);
    def("rec", (Real(*)(Real const&)) &rec);
    def("sqrt", (Real(*)(Real const&)) &sqrt);
    def("exp", (Real(*)(Real const&)) &exp);
    def("log", (Real(*)(Real const&)) &log);

    def("sin", (Real(*)(Real const&)) &sin);
    def("cos", (Real(*)(Real const&)) &cos);
    def("tan", (Real(*)(Real const&)) &tan);
    def("atan", (Real(*)(Real const&)) &atan);

    implicitly_convertible<int,Real>();
#ifdef HAVE_GMPXX_H
    implicitly_convertible<Integer,Real>();
    implicitly_convertible<Rational,Real>();
#endif
    implicitly_convertible<Decimal,Real>();
    implicitly_convertible<Dyadic,Real>();
}




void export_exact_float()
{
    class_< ExactFloat > exact_float_class("ExactFloat",init<ExactFloat>());
    exact_float_class.def(init<>());
    exact_float_class.def(init<Float>());
    exact_float_class.def(+self);
    exact_float_class.def(-self);
    exact_float_class.def(self+self);
    exact_float_class.def(self-self);
    exact_float_class.def(self*self);
    exact_float_class.def(self/self);

    exact_float_class.def("set_output_precision", &ExactFloat::set_output_precision);
    exact_float_class.staticmethod("set_output_precision");
    exact_float_class.def("__str__", &__cstr__<ExactFloat>);
    exact_float_class.def("__repr__", &__cstr__<ExactFloat>);
}

void export_validated_float()
{
    using boost::python::class_;
    using boost::python::init;
    using boost::python::self;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::def;

    def("down",&down);
    def("up",&up);

    class_< ValidatedFloat > validated_float_class("ValidatedFloat");
    validated_float_class.def(init<double>());
    validated_float_class.def(init<double,double>());
    validated_float_class.def(init<Float,Float>());
    validated_float_class.def(init<Real>());
    validated_float_class.def(init<Decimal>());
    validated_float_class.def(init<Dyadic>());
    validated_float_class.def(init<ValidatedFloat>());
    validated_float_class.def(init<Float>());
#ifdef HAVE_GMPXX_H
    validated_float_class.def(init<Rational>());
#endif
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
    validated_float_class.def("lower", &ValidatedFloat::lower);
    validated_float_class.def("upper", &ValidatedFloat::upper);
    validated_float_class.def("lower_value", &ValidatedFloat::lower_value, return_value_policy<copy_const_reference>());
    validated_float_class.def("upper_value", &ValidatedFloat::upper_value, return_value_policy<copy_const_reference>());
    validated_float_class.def("midpoint", &ValidatedFloat::midpoint);
    validated_float_class.def("radius", &ValidatedFloat::radius);
    validated_float_class.def("width", &ValidatedFloat::width);
    validated_float_class.def(boost::python::self_ns::str(self));
    validated_float_class.def("__repr__", &__repr__<ValidatedFloat>);

    validated_float_class.def("set_output_precision", &ValidatedFloat::set_output_precision);
    validated_float_class.staticmethod("set_output_precision");

    implicitly_convertible<Real,ValidatedFloat>();
    implicitly_convertible<ExactFloat,ValidatedFloat>();

    from_python_dict<ValidatedFloat>();
    from_python_list<ValidatedFloat>();
    //from_python_str<ValidatedFloat>();

    def("mag", (Float(*)(ValidatedFloat)) &mag);

    def("max", (ValidatedFloat(*)(ValidatedFloat,ValidatedFloat)) &max);
    def("min", (ValidatedFloat(*)(ValidatedFloat,ValidatedFloat)) &min);

    def("trunc", (ValidatedFloat(*)(ValidatedFloat,uint)) &trunc, "truncate to n binary digits");
    def("abs", (ValidatedFloat(*)(ValidatedFloat)) &abs, "validated_float absolute value function");
    def("pow",  (ValidatedFloat(*)(ValidatedFloat,int)) &pow, "validated_float power function");
    def("sqr", (ValidatedFloat(*)(ValidatedFloat)) &sqr, "validated_float square function");
    def("rec", (ValidatedFloat(*)(ValidatedFloat)) &rec);
    def("sqrt", (ValidatedFloat(*)(ValidatedFloat)) &sqrt);
    def("exp", (ValidatedFloat(*)(ValidatedFloat)) &exp);
    def("log", (ValidatedFloat(*)(ValidatedFloat)) &log);

    def("sin", (ValidatedFloat(*)(ValidatedFloat)) &sin);
    def("cos", (ValidatedFloat(*)(ValidatedFloat)) &cos);
    def("tan", (ValidatedFloat(*)(ValidatedFloat)) &tan);
    def("asin", (ValidatedFloat(*)(ValidatedFloat)) &asin);
    def("acos", (ValidatedFloat(*)(ValidatedFloat)) &acos);
    def("atan", (ValidatedFloat(*)(ValidatedFloat)) &atan);
}

void export_upper_float()
{
    class_< UpperFloat > upper_float_class("UpperFloat");
    upper_float_class.def(init<double>());
    upper_float_class.def(init<ValidatedFloat>());
    upper_float_class.def(init<ExactFloat>());

    upper_float_class.def(-self);
    upper_float_class.def(self + self);
    upper_float_class.def(self - LowerFloat());

    upper_float_class.def(boost::python::self_ns::str(self));
    implicitly_convertible<ValidatedFloat,UpperFloat>();
}

void export_lower_float()
{
    class_< LowerFloat > lower_float_class("LowerFloat");
    lower_float_class.def(init<double>());
    lower_float_class.def(init<Float>());
    lower_float_class.def(init<ValidatedFloat>());
    lower_float_class.def(init<ExactFloat>());

    lower_float_class.def(-self);
    lower_float_class.def(self + self);
    lower_float_class.def(self - UpperFloat());

    lower_float_class.def(boost::python::self_ns::str(self));
    implicitly_convertible<ValidatedFloat,LowerFloat>();
}

void export_approximate_float()
{
    using boost::python::class_;
    using boost::python::init;
    using boost::python::self;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::def;

    def("down",&down);
    def("up",&up);

    class_< ApproximateFloat > approximate_float_class("ApproximateFloat");
    approximate_float_class.def(init<double>());
    approximate_float_class.def(init<Float>());
    approximate_float_class.def(init<Decimal>());
    approximate_float_class.def(init<Dyadic>());
    approximate_float_class.def(init<ApproximateFloat>());
    approximate_float_class.def(init<LowerFloat>());
    approximate_float_class.def(init<UpperFloat>());
    approximate_float_class.def(init<ValidatedFloat>());
#ifdef HAVE_GMPXX_H
    approximate_float_class.def(init<Rational>());
#endif
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
    approximate_float_class.def("__repr__", &__repr__<ApproximateFloat>);

    approximate_float_class.def("set_output_precision", &ApproximateFloat::set_output_precision);
    approximate_float_class.staticmethod("set_output_precision");

    implicitly_convertible<LowerFloat,ApproximateFloat>();
    implicitly_convertible<UpperFloat,ApproximateFloat>();
    implicitly_convertible<ValidatedFloat,ApproximateFloat>();

    def("max", (ApproximateFloat(*)(ApproximateFloat,ApproximateFloat)) &max);
    def("min", (ApproximateFloat(*)(ApproximateFloat,ApproximateFloat)) &min);

    def("abs", (ApproximateFloat(*)(ApproximateFloat)) &abs);
    def("pow",  (ApproximateFloat(*)(ApproximateFloat,int)) &pow);
    def("sqr", (ApproximateFloat(*)(ApproximateFloat)) &sqr);
    def("rec", (ApproximateFloat(*)(ApproximateFloat)) &rec);
    def("sqrt", (ApproximateFloat(*)(ApproximateFloat)) &sqrt);
    def("exp", (ApproximateFloat(*)(ApproximateFloat)) &exp);
    def("log", (ApproximateFloat(*)(ApproximateFloat)) &log);

    def("sin", (ApproximateFloat(*)(ApproximateFloat)) &sin);
    def("cos", (ApproximateFloat(*)(ApproximateFloat)) &cos);
    def("tan", (ApproximateFloat(*)(ApproximateFloat)) &tan);
    def("asin", (ApproximateFloat(*)(ApproximateFloat)) &asin);
    def("acos", (ApproximateFloat(*)(ApproximateFloat)) &acos);
    def("atan", (ApproximateFloat(*)(ApproximateFloat)) &atan);
}

void export_float()
{
    class_< Float > float_class("Float");
    float_class.def(init<double>());
    float_class.def(init<Float>());
    float_class.def(boost::python::self_ns::str(self));
    float_class.def("__repr__", &__repr__<Float>);

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

    implicitly_convertible<double,Float>();

    def("min",(Float(*)(Float,Float)) &Ariadne::min);
    def("max",(Float(*)(Float,Float)) &Ariadne::max);
    def("abs", (Float(*)(Float)) &Ariadne::abs);

    def("pow",(Float(*)(Float,int)) &Ariadne::pow);

    def("rec",(Float(*)(Float)) &Ariadne::rec);
    def("sqr",(Float(*)(Float)) &Ariadne::sqr);
    def("sqrt", (Float(*)(Float)) &Ariadne::sqrt);

    def("exp", (Float(*)(Float)) &Ariadne::exp);
    def("log", (Float(*)(Float)) &Ariadne::log);

    def("sin", (Float(*)(Float)) &Ariadne::sin);
    def("cos", (Float(*)(Float)) &Ariadne::cos);
    def("tan", (Float(*)(Float)) &Ariadne::tan);
    def("asin", (Float(*)(Float)) &Ariadne::asin);
    def("acos", (Float(*)(Float)) &Ariadne::acos);
    def("atan", (Float(*)(Float)) &Ariadne::atan);

}


void
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
    export_dyadic();
    export_decimal();
#ifdef HAVE_GMPXX_H
    export_integer();
    export_rational();
#endif
}
