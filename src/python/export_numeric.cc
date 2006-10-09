/***************************************************************************
 *            python/export_numeric.cc
 *
 *  22 June 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include "numeric/integer.h"
#include "numeric/rational.h"
#include "numeric/dyadic.h"
#include "numeric/float64.h"
#include "numeric/mpfloat.h"

#include "numeric/interval.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;

#include <boost/python.hpp>
using namespace boost::python;

// Give some mixed arithmetic operators for MPFloat
template<typename R> inline Interval<MPFloat> operator+(const MPFloat& x1, const R& x2) {
  return x1+MPFloat(x2); }
template<typename R> inline Interval<MPFloat> operator-(const MPFloat& x1, const R& x2) {
  return x1-MPFloat(x2); }
template<typename R> inline Interval<MPFloat> operator*(const MPFloat& x1, const R& x2) {
  return x1*MPFloat(x2); }
template<typename R> inline Interval<MPFloat> operator/(const MPFloat& x1, const R& x2) {
  return x1/MPFloat(x2); }
  
template<typename R> inline Interval<MPFloat> operator-(const R& x1, const MPFloat& x2) {
  return MPFloat(x1)-x2; }
template<typename R> inline Interval<MPFloat> operator/(const R& x1, const MPFloat& x2) {
  return MPFloat(x1)/x2; }

template<typename R> inline bool operator==(const MPFloat& x1, const R& x2) {
  return x1==MPFloat(x2); }
template<typename R> inline bool operator!=(const MPFloat& x1, const R& x2) {
  return x1!=MPFloat(x2); }
template<typename R> inline bool operator<=(const MPFloat& x1, const R& x2) {
  return x1<=MPFloat(x2); }
template<typename R> inline bool operator>=(const MPFloat& x1, const R& x2) {
  return x1>=MPFloat(x2); }
template<typename R> inline bool operator< (const MPFloat& x1, const R& x2) {
  return x1< MPFloat(x2); }
template<typename R> inline bool operator> (const MPFloat& x1, const R& x2) {
  return x1> MPFloat(x2); }

  
// Since gmpxx operators return intermediate types which need conversion, need to explicitly provide arithmetic functions
inline Integer neg_z(const Integer& z1) { return -z1; }
inline Integer add_zz(const Integer& z1, const Integer& z2) { return z1+z2; }
inline Integer add_zi(const Integer& z1, const int& z2) { return z1+z2; }
inline Integer sub_zz(const Integer& z1, const Integer& z2) { return z1-z2; }
inline Integer sub_zi(const Integer& z1, const int& z2) { return z1-z2; }
inline Integer rsub_zi(const Integer& z1, const int& z2) { return z2-z1; }
inline Integer mul_zz(const Integer& z1, const Integer& z2) { return z1*z2; }
inline Integer mul_zi(const Integer& z1, const int& z2) { return z1*z2; }
inline Rational div_zz(const Integer& z1, const Integer& z2) { return Rational(z1)/z2; }
inline bool eq_zz(const Integer& z1, const Integer& z2) { return z1==z2; }
inline bool eq_zi(const Integer& z1, const int& z2) { return z1==z2; }
inline bool eq_zx(const Integer& z1, const double& z2) { return z1==z2; }
inline bool ne_zz(const Integer& z1, const Integer& z2) { return z1!=z2; }
inline bool ne_zi(const Integer& z1, const int& z2) { return z1!=z2; }
inline bool ne_zx(const Integer& z1, const double& z2) { return z1!=z2; }
inline bool lt_zz(const Integer& z1, const Integer& z2) { return z1<z2; }
inline bool lt_zi(const Integer& z1, const int& z2) { return z1<z2; }
inline bool lt_zx(const Integer& z1, const double& z2) { return z1<z2; }
inline bool gt_zi(const Integer& z1, const int& z2) { return z1>z2; }
inline bool gt_zx(const Integer& z1, const double& z2) { return z1>z2; }
inline bool le_zz(const Integer& z1, const Integer& z2) { return z1<=z2; }
inline bool le_zi(const Integer& z1, const int& z2) { return z1<=z2; }
inline bool le_zx(const Integer& z1, const double& z2) { return z1<=z2; }
inline bool ge_zi(const Integer& z1, const int& z2) { return z1>=z2; }
inline bool ge_zx(const Integer& z1, const double& z2) { return z1>=z2; }

inline MPFloat neg_f(const MPFloat& f1) { return -f1; }
inline Interval<MPFloat> add_ff(const MPFloat& f1, const MPFloat& f2) { return f1+f2; }
inline Interval<MPFloat> add_fx(const MPFloat& f1, const double& f2) { return f1+f2; }
inline Interval<MPFloat> sub_ff(const MPFloat& f1, const MPFloat& f2) { return f1-f2; }
inline Interval<MPFloat> sub_fx(const MPFloat& f1, const double& f2) { return f1-f2; }
inline Interval<MPFloat> rsub_fx(const MPFloat& f1, const double& f2) { return f2-f1; }
inline Interval<MPFloat> mul_ff(const MPFloat& f1, const MPFloat& f2) { return f1*f2; }
inline Interval<MPFloat> mul_fx(const MPFloat& f1, const double& f2) { return f1*f2; }
inline Interval<MPFloat> div_ff(const MPFloat& f1, const MPFloat& f2) { return f1/f2; }
inline Interval<MPFloat> div_fn(const MPFloat& f1, const uint& f2) { return f1/f2; }
inline bool eq_ff(const MPFloat& f1, const MPFloat& f2) { return f1==f2; }
inline bool eq_fx(const MPFloat& f1, const double& f2) { return f1==f2; }
inline bool eq_fz(const MPFloat& f1, const Integer& f2) { return f1==f2; }
inline bool ne_ff(const MPFloat& f1, const MPFloat& f2) { return f1!=f2; }
inline bool ne_fx(const MPFloat& f1, const double& f2) { return f1!=f2; }
inline bool ne_fz(const MPFloat& f1, const Integer& f2) { return f1!=f2; }
inline bool lt_ff(const MPFloat& f1, const MPFloat& f2) { return f1<f2; }
inline bool lt_fx(const MPFloat& f1, const double& f2) { return f1<f2; }
inline bool lt_fz(const MPFloat& f1, const Integer& f2) { return f1<f2; }
inline bool gt_ff(const MPFloat& f1, const MPFloat& f2) { return f1>f2; }
inline bool gt_fx(const MPFloat& f1, const double& f2) { return f1>f2; }
inline bool gt_fz(const MPFloat& f1, const Integer& f2) { return f1>f2; }
inline bool le_ff(const MPFloat& f1, const MPFloat& f2) { return f1<=f2; }
inline bool le_fx(const MPFloat& f1, const double& f2) { return f1<=f2; }
inline bool le_fz(const MPFloat& f1, const Integer& f2) { return f1<=f2; }
inline bool ge_ff(const MPFloat& f1, const MPFloat& f2) { return f1>=f2; }
inline bool ge_fx(const MPFloat& f1, const double& f2) { return f1>=f2; }
inline bool ge_fz(const MPFloat& f1, const Integer& f2) { return f1>=f2; }
//inline MPFloat mantissa_f(const MPFloat f) { return mantissa(f); }
//inline Integer exponent_f(const MPFloat f) { return exponent(f); }
inline Integer precision_f(const MPFloat f) { return precision(f); }

inline Dyadic neg_d(const Dyadic& d1) { return -d1; }
inline Dyadic add_dd(const Dyadic& d1, const Dyadic& d2) { return d1+d2; }
inline Dyadic add_dx(const Dyadic& d1, const double& d2) { return d1+d2; }
inline Dyadic sub_dd(const Dyadic& d1, const Dyadic& d2) { return d1-d2; }
inline Dyadic sub_dx(const Dyadic& d1, const double& d2) { return d1-d2; }
inline Dyadic rsub_dx(const Dyadic& d1, const double& d2) { return d2-d1; }
inline Dyadic mul_dd(const Dyadic& d1, const Dyadic& d2) { return d1*d2; }
inline Dyadic mul_dx(const Dyadic& d1, const double& d2) { return d1*d2; }
inline Dyadic div_dn(const Dyadic& d1, const uint& d2) { return d1/d2; }
inline bool eq_dd(const Dyadic& d1, const Dyadic& d2) { return d1==d2; }
inline bool eq_dx(const Dyadic& d1, const double& d2) { return d1==d2; }
inline bool eq_dz(const Dyadic& n1, const Integer& n2) { return n1==n2; }
inline bool ne_dd(const Dyadic& d1, const Dyadic& d2) { return d1!=d2; }
inline bool ne_dx(const Dyadic& d1, const double& d2) { return d1!=d2; }
inline bool ne_dz(const Dyadic& n1, const Integer& n2) { return n1!=n2; }
inline bool lt_dd(const Dyadic& d1, const Dyadic& d2) { return d1<d2; }
inline bool lt_dx(const Dyadic& d1, const double& d2) { return d1<d2; }
inline bool lt_dz(const Dyadic& n1, const Integer& n2) { return n1<n2; }
inline bool gt_dd(const Dyadic& d1, const Dyadic& d2) { return d1>d2; }
inline bool gt_dx(const Dyadic& d1, const double& d2) { return d1>d2; }
inline bool gt_dz(const Dyadic& n1, const Integer& n2) { return n1>n2; }
inline bool le_dd(const Dyadic& d1, const Dyadic& d2) { return d1<=d2; }
inline bool le_dx(const Dyadic& d1, const double& d2) { return d1<=d2; }
inline bool le_dz(const Dyadic& n1, const Integer& n2) { return n1<=n2; }
inline bool ge_dd(const Dyadic& d1, const Dyadic& d2) { return d1>=d2; }
inline bool ge_dx(const Dyadic& d1, const double& d2) { return d1>=d2; }
inline bool ge_dz(const Dyadic& n1, const Integer& n2) { return n1>=n2; }
inline Dyadic mantissa_d(const Dyadic d) { return mantissa(d); }
inline int exponent_d(const Dyadic d) { return exponent(d); }
inline int precision_d(const Dyadic d) { return precision(d); }

inline Rational neg_q(const Rational& q1) { return -q1; }
inline Rational add_qq(const Rational& q1, const Rational& q2) { return q1+q2; }
inline Rational add_qx(const Rational& q1, const double& q2) { return q1+q2; }
inline Rational sub_qq(const Rational& q1, const Rational& q2) { return q1-q2; }
inline Rational sub_qx(const Rational& q1, const double& q2) { return q1-q2; }
inline Rational rsub_qx(const Rational& q1, const double& q2) { return q2-q1; }
inline Rational mul_qq(const Rational& q1, const Rational& q2) { return q1*q2; }
inline Rational mul_qx(const Rational& q1, const double& q2) { return q1*q2; }
inline Rational div_qq(const Rational& q1, const Rational& q2) { return q1/q2; }
inline Rational div_qx(const Rational& q1, const double& q2) { return q1/q2; }
inline Rational rdiv_qx(const Rational& q1, const double& q2) { return q2/q1; }
inline bool eq_qq(const Rational& q1, const Rational& q2) { return q1==q2; }
inline bool eq_qx(const Rational& q1, const double& q2) { return q1==q2; }
inline bool ne_qq(const Rational& q1, const Rational& q2) { return q1!=q2; }
inline bool ne_qx(const Rational& q1, const double& q2) { return q1!=q2; }
inline bool lt_qq(const Rational& q1, const Rational& q2) { return q1<q2; }
inline bool lt_qx(const Rational& q1, const double& q2) { return q1<q2; }
inline bool gt_qx(const Rational& q1, const double& q2) { return q1>q2; }
inline bool le_qq(const Rational& q1, const Rational& q2) { return q1<=q2; }
inline bool le_qx(const Rational& q1, const double& q2) { return q1<=q2; }
inline bool ge_qx(const Rational& q1, const double& q2) { return q1>=q2; }

inline bool eq_qd(const Rational& n1, const Dyadic& n2) { return n1==n2; }
inline bool eq_qz(const Rational& n1, const Integer& n2) { return n1==n2; }
inline bool ne_qd(const Rational& n1, const Dyadic& n2) { return n1!=n2; }
inline bool ne_qz(const Rational& n1, const Integer& n2) { return n1!=n2; }
inline bool lt_qd(const Rational& n1, const Dyadic& n2) { return n1<n2; }
inline bool lt_qz(const Rational& n1, const Integer& n2) { return n1<n2; }
inline bool gt_qd(const Rational& n1, const Dyadic& n2) { return n1>n2; }
inline bool gt_qz(const Rational& n1, const Integer& n2) { return n1>n2; }
inline bool le_qd(const Rational& n1, const Dyadic& n2) { return n1<=n2; }
inline bool le_qz(const Rational& n1, const Integer& n2) { return n1<=n2; }
inline bool ge_qd(const Rational& n1, const Dyadic& n2) { return n1>=n2; }
inline bool ge_qz(const Rational& n1, const Integer& n2) { return n1>=n2; }
inline Integer numerator_q(const Rational& q) { return numerator(q); }
inline Integer denominator_q(const Rational& q) { return denominator(q); }

void export_numeric() {
  class_<Integer>("Integer")
    .def(init<int>())
    .def(init<Integer>())
    .def("__neg__", &neg_z)
    .def("__add__", &add_zz)
    .def("__add__", &add_zi)
    .def("__radd__", &add_zi)
    .def("__sub__", &sub_zz)
    .def("__sub__", &sub_zi)
    .def("__rsub__", &rsub_zi)
    .def("__mul__", &mul_zz)
    .def("__mul__", &mul_zi)
    .def("__rmul__", &mul_zi)
    .def("__eq__", &eq_zz)
    .def("__eq__", &eq_zx)
    .def("__ne__", &ne_zz)
    .def("__ne__", &ne_zx)
    .def("__lt__", &lt_zz)
    .def("__lt__", &lt_zx)
    .def("__gt__", &gt_zx)
    .def("__le__", &le_zz)
    .def("__le__", &le_zx)
    .def("__ge__", &le_zx)
    .def("__div__", &div_zz)
    .def(self_ns::str(self))
  ;

  class_<MPFloat>("MPFloat")
    .def(init<int>())
    .def(init<Integer>())
    .def(init<double>())
    .def(init<MPFloat>())
    .def("__neg__", &neg_f)
    .def("__add__", &add_ff)
    .def("__add__", &add_fx)
    .def("__radd__", &add_fx)
    .def("__sub__", &sub_ff)
    .def("__sub__", &sub_fx)
    .def("__rsub__", &rsub_fx)
    .def("__mul__", &mul_ff)
    .def("__mul__", &mul_fx)
    .def("__rmul__", &mul_fx)
    .def("__div__", &div_fn)
    .def("__eq__", &eq_ff)
    .def("__eq__", &eq_fx)
    .def("__eq__", &eq_fz)
    .def("__ne__", &ne_ff)
    .def("__ne__", &ne_fx)
    .def("__ne__", &ne_fz)
    .def("__lt__", &lt_ff)
    .def("__lt__", &lt_fx)
    .def("__lt__", &lt_fz)
    .def("__gt__", &gt_ff)
    .def("__gt__", &gt_fx)
    .def("__gt__", &gt_fz)
    .def("__le__", &le_ff)
    .def("__le__", &le_fx)
    .def("__le__", &le_fz)
    .def("__ge__", &ge_ff)
    .def("__ge__", &ge_fx)
    .def("__ge__", &ge_fz)
    .def("precision", &precision_f)
//    .def("exponent", &exponent_f)
//    .def("mantissa", &mantissa_f)
    .def(self_ns::str(self))
  ;

  class_<Dyadic>("Dyadic")
    .def(init<int,int>())
    .def(init<int>())
    .def(init<Integer>())
    .def(init<double>())
//    .def(init<MPFloat>())
    .def(init<Dyadic>())
    .def("__neg__", &neg_d)
    .def("__add__", &add_dd)
    .def("__add__", &add_dx)
    .def("__radd__", &add_dx)
    .def("__sub__", &sub_dd)
    .def("__sub__", &sub_dx)
    .def("__rsub__", &rsub_dx)
    .def("__mul__", &mul_dd)
    .def("__mul__", &mul_dx)
    .def("__rmul__", &mul_dx)
    .def("__div__", &div_dn)
    .def("__eq__", &eq_dd)
    .def("__eq__", &eq_dx)
    .def("__eq__", &eq_dz)
    .def("__ne__", &ne_dd)
    .def("__ne__", &ne_dx)
    .def("__ne__", &ne_dz)
    .def("__lt__", &lt_dd)
    .def("__lt__", &lt_dx)
    .def("__lt__", &lt_dz)
    .def("__gt__", &gt_dx)
    .def("__gt__", &gt_dz)
    .def("__le__", &le_dd)
    .def("__le__", &le_dx)
    .def("__le__", &le_dz)
    .def("__ge__", &ge_dd)
    .def("__ge__", &ge_dx)
    .def("__ge__", &ge_dz)
    .def("precision", &precision_d)
    .def("exponent", &exponent_d)
    .def("mantissa", &mantissa_d)
    .def(self_ns::str(self))
  ;

  class_<Rational>("Rational")
    .def(init<int,int>())
    .def(init<int>())
    .def(init<Integer>())
    .def(init<double>())
    .def(init<MPFloat>())
    .def(init<Dyadic>())
    .def(init<Rational>())
    .def("__neg__", &neg_q)
    .def("__add__", &add_qq)
    .def("__add__", &add_qx)
    .def("__radd__", &add_qx)
    .def("__sub__", &sub_qq)
    .def("__sub__", &sub_qx)
    .def("__rsub__", &rsub_qx)
    .def("__mul__", &mul_qq)
    .def("__mul__", &mul_qx)
    .def("__rmul__", &mul_qx)
    .def("__div__", &div_qq)
    .def("__div__", &div_qx)
    .def("__rdiv__", &rdiv_qx)
    .def("__eq__", &eq_qq)
    .def("__eq__", &eq_qx)
    .def("__eq__", &eq_qz)
    .def("__eq__", &eq_qd)
    .def("__ne__", &ne_qq)
    .def("__ne__", &ne_qx)
    .def("__ne__", &ne_qz)
    .def("__ne__", &ne_qd)
    .def("__lt__", &lt_qq)
    .def("__lt__", &lt_qx)
    .def("__lt__", &lt_qz)
    .def("__lt__", &lt_qd)
    .def("__gt__", &gt_qx)
    .def("__gt__", &gt_qz)
    .def("__gt__", &gt_qd)
    .def("__le__", &le_qq)
    .def("__le__", &le_qx)
    .def("__le__", &le_qz)
    .def("__le__", &le_qd)
    .def("__ge__", &ge_qx)
    .def("__ge__", &ge_qz)
    .def("__ge__", &ge_qd)
    .def("numerator", &numerator_q)
    .def("denominator", &denominator_q)
    .def(self_ns::str(self))
  ;
  
/*
inline Real neg_r(const Real& r) { return -r; }
inline Real add_rr(const Real& r1, const Real& r2) { return r1+r2; }
inline Real add_rx(const Real& r1, const double& r2) { return r1+r2; }
inline Real sub_rr(const Real& r1, const Real& r2) { return r1-r2; }
inline Real sub_rx(const Real& r1, const double& r2) { return r1-r2; }
inline Real rsub_rx(const Real& r1, const double& r2) { return r2-r1; }
inline Real mul_rr(const Real& r1, const Real& r2) { return r1*r2; }
inline Real mul_rx(const Real& r1, const double& r2) { return r1*r2; }
inline Real div_rr(const Real& r1, const Real& r2) { return r1/r2; }
inline Real div_rn(const Real& r1, const uint& r2) { return r1/r2; }
inline bool eq_rq(const Real& r1, const Rational& r2) { return r1==r2; }
inline bool eq_rr(const Real& r1, const Real& r2) { return r1==r2; }
inline bool eq_rx(const Real& r1, const double& r2) { return r1==r2; }
inline bool eq_rz(const Real& r1, const Integer& r2) { return r1==r2; }
inline bool ne_rq(const Real& r1, const Rational& r2) { return r1!=r2; }
inline bool ne_rr(const Real& r1, const Real& r2) { return r1!=r2; }
inline bool ne_rx(const Real& r1, const double& r2) { return r1!=r2; }
inline bool ne_rz(const Real& r1, const Integer& r2) { return r1!=r2; }
inline bool lt_rq(const Real& r1, const Rational& r2) { return r1<r2; }
inline bool lt_rr(const Real& r1, const Real& r2) { return r1<r2; }
inline bool lt_rx(const Real& r1, const double& r2) { return r1<r2; }
inline bool lt_rz(const Real& r1, const Integer& r2) { return r1<r2; }
inline bool gt_rq(const Real& r1, const Rational& r2) { return r1>r2; }
inline bool gt_rr(const Real& r1, const Real& r2) { return r1>r2; }
inline bool gt_rx(const Real& r1, const double& r2) { return r1>r2; }
inline bool gt_rz(const Real& r1, const Integer& r2) { return r1>r2; }
inline bool le_rq(const Real& r1, const Rational& r2) { return r1<=r2; }
inline bool le_rr(const Real& r1, const Real& r2) { return r1<=r2; }
inline bool le_rx(const Real& r1, const double& r2) { return r1<=r2; }
inline bool le_rz(const Real& r1, const Integer& r2) { return r1<=r2; }
inline bool ge_rq(const Real& r1, const Rational& r2) { return r1>=r2; }
inline bool ge_rr(const Real& r1, const Real& r2) { return r1>=r2; }
inline bool ge_rx(const Real& r1, const double& r2) { return r1>=r2; }
inline bool ge_rz(const Real& r1, const Integer& r2) { return r1>=r2; }
inline Real div_approx_rr(const Real r1, const Real& r2, const Real& e) { return div_approx(r1,r2,e);  }
inline Real div_approx_rx(const Real r1, const double& r2, const Real& e) { return div_approx(r1,Real(r2),e);  }
inline Real div_approx_xr(const double r1, const Real& r2, const Real& e) { return div_approx(Real(r1),r2,e);  }
inline Real div_prec_rr(const Real r1, const Real& r2, const int& n) { return div_approx(r1,r2,n);  }
inline Real mantissa_r(const Real r) { return mantissa(r); }
inline Integer exponent_r(const Real r) { return exponent(r); }
inline Integer precision_r(const Real r) { return precision(r); }  class_<Real>("Real")
    .def(init<int>())
    .def(init<double>())
    .def(init<int,int>())
    .def(init<Real>())
    .def("__neg__", &neg_r)
    .def("__add__", &add_rr)
    .def("__add__", &add_rx)
    .def("__radd__", &add_rx)
    .def("__sub__", &sub_rr)
    .def("__sub__", &sub_rx)
    .def("__rsub__", &rsub_rx)
    .def("__mul__", &mul_rr)
    .def("__mul__", &mul_rx)
    .def("__rmul__", &mul_rx)
    .def("__div__", &div_rn)
    .def("__eq__", &eq_rq)
    .def("__eq__", &eq_rr)
    .def("__eq__", &eq_rx)
    .def("__eq__", &eq_rz)
    .def("__ne__", &ne_rq)
    .def("__ne__", &ne_rr)
    .def("__ne__", &ne_rx)
    .def("__ne__", &ne_rz)
    .def("__lt__", &lt_rq)
    .def("__lt__", &lt_rr)
    .def("__lt__", &lt_rx)
    .def("__lt__", &lt_rz)
    .def("__gt__", &gt_rq)
    .def("__gt__", &gt_rr)
    .def("__gt__", &gt_rx)
    .def("__gt__", &gt_rz)
    .def("__le__", &le_rq)
    .def("__le__", &le_rr)
    .def("__le__", &le_rx)
    .def("__le__", &le_rz)
    .def("__ge__", &ge_rq)
    .def("__ge__", &ge_rr)
    .def("__ge__", &ge_rx)
    .def("__ge__", &ge_rz)
    .def("precision", &precision_r)
    .def("exponent", &exponent_r)
    .def("mantissa", &mantissa_r)
    .def(self_ns::str(self))
  ;
*/

  
}
