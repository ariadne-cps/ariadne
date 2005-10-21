// g++ -fpic -I../include/ -I/usr/include/python2.4 -c interval.cc
// g++ -shared interval.o -o interval.so -lboost_python -lmpfr -lgmp  -lsuperlu -lblas

#include <iostream>

#include <boost/python.hpp>

#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/io.hpp>

#include <gmp.h>
#include <mpfr.h>


struct full_rounding :
  boost::numeric::interval_lib::rounded_arith_opp<double>
{
 private:
  typedef int mpfr_func(mpfr_t, const __mpfr_struct*, mp_rnd_t);
  
  double invoke_mpfr(double x, mpfr_func f, mp_rnd_t r) {
    mpfr_t xx;
    mpfr_init_set_d(xx, x, GMP_RNDN);
    f(xx, xx, r);
    double res = mpfr_get_d(xx, r);
    mpfr_clear(xx);
    return res;
  }
 public:
  double exp_down(double x) { return invoke_mpfr(x, mpfr_exp, GMP_RNDD); } 
  double log_down(double x) { return invoke_mpfr(x, mpfr_log, GMP_RNDD); } 
  double sin_down(double x) { return invoke_mpfr(x, mpfr_sin, GMP_RNDD); } 
  double cos_down(double x) { return invoke_mpfr(x, mpfr_cos, GMP_RNDD); } 
  double tan_down(double x) { return invoke_mpfr(x, mpfr_tan, GMP_RNDD); } 
  double asin_down(double x) { return invoke_mpfr(x, mpfr_asin, GMP_RNDD); } 
  double acos_down(double x) { return invoke_mpfr(x, mpfr_acos, GMP_RNDD); } 
  double atan_down(double x) { return invoke_mpfr(x, mpfr_atan, GMP_RNDD); } 
  double sinh_down(double x) { return invoke_mpfr(x, mpfr_sinh, GMP_RNDD); } 
  double cosh_down(double x) { return invoke_mpfr(x, mpfr_cosh, GMP_RNDD); } 
  double tanh_down(double x) { return invoke_mpfr(x, mpfr_tanh, GMP_RNDD); }
  double asinh_down(double x) { return invoke_mpfr(x, mpfr_asinh, GMP_RNDD); } 
  double acosh_down(double x) { return invoke_mpfr(x, mpfr_acosh, GMP_RNDD); } 
  double atanh_down(double x) { return invoke_mpfr(x, mpfr_atanh, GMP_RNDD); }

  double exp_up (double x) { return invoke_mpfr(x, mpfr_exp, GMP_RNDU); }
  double log_up (double x) { return invoke_mpfr(x, mpfr_log, GMP_RNDU); }
  double sin_up (double x) { return invoke_mpfr(x, mpfr_sin, GMP_RNDU); }
  double cos_up (double x) { return invoke_mpfr(x, mpfr_cos, GMP_RNDU); }
  double tan_up (double x) { return invoke_mpfr(x, mpfr_tan, GMP_RNDU); }
  double asin_up (double x) { return invoke_mpfr(x, mpfr_asin, GMP_RNDU); }
  double acos_up (double x) { return invoke_mpfr(x, mpfr_acos, GMP_RNDU); }
  double atan_up (double x) { return invoke_mpfr(x, mpfr_atan, GMP_RNDU); }
  double sinh_up (double x) { return invoke_mpfr(x, mpfr_sinh, GMP_RNDU); }
  double cosh_up (double x) { return invoke_mpfr(x, mpfr_cosh, GMP_RNDU); }
  double tanh_up (double x) { return invoke_mpfr(x, mpfr_tanh, GMP_RNDU); }
  double asinh_up (double x) { return invoke_mpfr(x, mpfr_asinh, GMP_RNDU); }
  double acosh_up (double x) { return invoke_mpfr(x, mpfr_acosh, GMP_RNDU); }
  double atanh_up (double x) { return invoke_mpfr(x, mpfr_atanh, GMP_RNDU); }

};

using boost::numeric::interval_lib::policies;
using boost::numeric::interval_lib::save_state;
using boost::numeric::interval_lib::checking_strict;

typedef boost::numeric::interval<double, policies<save_state<full_rounding>,checking_strict<double> > > Interval;



#include <sla/sla.h>
typedef SLA::Vector<Interval> Vector;
typedef SLA::Matrix<Interval> Matrix;

typedef Vector vector_unary_func(const Vector& v);
typedef Vector vector_binary_func(const Vector& v, const Vector& v);

Vector operator-(const Vector& v) {
  return SLA::operator-(v); }
Vector operator+(const Vector& v1, const Vector& v2) {
  return SLA::operator+(v1,v2); }
Vector operator-(const Vector& v1, const Vector& v2) {
  return SLA::operator-(v1,v2); }
Vector operator*(const Interval& i, const Vector& v) {
  return SLA::operator*(v,i); }
Vector operator*(const Vector& v, const Interval& i) {
  return SLA::operator*(v,i); }
Vector operator/(const Vector& v, const Interval& i) {
  return SLA::operator/(v,i); }

Matrix operator-(const Matrix& a) {
  return SLA::operator-(a); }
Matrix operator+(const Matrix& a1, const Matrix& a2) {
  return SLA::operator+(a1,a2); }
Matrix operator-(const Matrix& a1, const Matrix& a2) {
  return SLA::operator-(a1,a2); }
Matrix operator*(const Interval& i, const Matrix& a) {
  return SLA::operator*(a,i); }
Matrix operator*(const Matrix& a, const Interval& i) {
  return SLA::operator*(a,i); }
Vector operator*(const Matrix& a, const Vector& v) {
  return SLA::operator*(a,v); }
Matrix operator*(const Matrix& a, const Matrix& v) {
  return SLA::operator*(a,v); }
Matrix operator/(const Matrix& a, const Interval& i) {
  return SLA::operator/(a,i); }


BOOST_PYTHON_MODULE(interval)
{
    using boost::python::class_;
    using boost::python::init;
    using boost::python::self;
    using boost::python::return_value_policy;
    using boost::python::copy_const_reference;
    using boost::python::def;

    typedef Interval (*IFUN)(const Interval&);

    class_< Interval >("Interval")
      .def(init<double,double>())
      .def(init<double>())
      .def(-self)        // __neg__
      .def(self + self)  // __add__
      .def(self - self)  // __sub__
      .def(self * self)  // __mul__
      .def(self / self)  // __div__
      .def(self + double())  // __add__
      .def(self - double())  // __sub__
      .def(self * double())  // __mul__
      .def(self / double())  // __div__
      .def(double() + self)  // __add__
      .def(double() - self)  // __sub__
      .def(double() * self)  // __mul__
      .def(double() / self)  // __div__
      .def("lower", &Interval::lower, return_value_policy<copy_const_reference>())
      .def("upper", &Interval::upper, return_value_policy<copy_const_reference>())
      .def(boost::python::self_ns::str(self))    // __str__
    ;

    IFUN iabs(&boost::numeric::abs);
    IFUN iexp(&boost::numeric::exp);
    IFUN ilog(&boost::numeric::log);
    IFUN isin(&boost::numeric::sin);
    IFUN icos(&boost::numeric::cos);
    IFUN itan(&boost::numeric::tan);
    IFUN iasin(&boost::numeric::asin);
    IFUN iacos(&boost::numeric::acos);
    IFUN iatan(&boost::numeric::atan);
    IFUN isinh(&boost::numeric::sinh);
    IFUN icosh(&boost::numeric::cosh);
    IFUN itanh(&boost::numeric::tanh);
    IFUN iasinh(&boost::numeric::asinh);
    IFUN iacosh(&boost::numeric::acosh);
    IFUN iatanh(&boost::numeric::atanh);


    def("abs", iabs, "interval absolute value function");
    def("exp", iexp);
    def("log", ilog);
    def("sin", isin);
    def("cos", icos);
    def("tan", itan);
    def("asin", iasin);
    def("acos", iacos);
    def("atan", iatan);
    def("sinh", isinh);
    def("cosh", icosh);
    def("tanh", itanh);
    def("asinh", iasinh);
    def("acosh", iacosh);
    def("atanh", iatanh);

    class_< Vector >("Vector")
      .def(init<int>())
      .def(-self)        // __neg__
      .def(self + self)  // __add__
      .def(self - self)  // __sub__
      .def(Interval() * self)  // __mul__
      .def(self * Interval())  // __mul__
      .def(self / Interval())  // __div__
      .def(boost::python::self_ns::str(self))    // __str__
      .def("__setitem__", &Vector::set)
      .def("__getitem__", &Vector::get, return_value_policy<copy_const_reference>())
      ;

    class_< Matrix >("Matrix")
      .def(init<int,int>())
      .def(-self)        // __neg__
      .def(self + self)  // __add__
      .def(self - self)  // __sub__
      .def(Interval() * self)  // __mul__
      .def(self * Interval())  // __mul__
      .def(self * Vector())    // __mul__
      .def(self * Matrix())    // __mul__
      .def(self / Interval())  // __div__
      .def(boost::python::self_ns::str(self))    // __str__
      .def("set", &Matrix::set)
      .def("get", &Matrix::get, return_value_policy<copy_const_reference>())
      .def("inverse", &Matrix::inverse)
      .def("determinant", &Matrix::determinant)
      .def("transpose", &Matrix::transpose)
      .def("solve", &Matrix::solve)
      .def("__setitem__", &Matrix::set)
      .def("__getitem__", &Matrix::get, return_value_policy<copy_const_reference>())
      ;


}
