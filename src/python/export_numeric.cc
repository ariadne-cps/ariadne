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
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "base/numerical_type.h"
#include "base/arithmetic.h"

#include <boost/python.hpp>

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::self_ns::str;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

using Ariadne::Integer;
using Ariadne::Dyadic;
using Ariadne::Rational;

using Ariadne::div_approx;
using Ariadne::precision;
using Ariadne::mantissa;
using Ariadne::exponent;
using Ariadne::numerator;
using Ariadne::denominator;

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

inline Dyadic neg_d(const Dyadic d1) { return -d1; }
inline Dyadic add_dd(const Dyadic d1, const Dyadic& d2) { return d1+d2; }
inline Dyadic add_dx(const Dyadic d1, const double& d2) { return d1+d2; }
inline Dyadic sub_dd(const Dyadic d1, const Dyadic& d2) { return d1-d2; }
inline Dyadic sub_dx(const Dyadic d1, const double& d2) { return d1-d2; }
inline Dyadic rsub_dx(const Dyadic d1, const double& d2) { return d2-d1; }
inline Dyadic mul_dd(const Dyadic d1, const Dyadic& d2) { return d1*d2; }
inline Dyadic mul_dx(const Dyadic d1, const double& d2) { return d1*d2; }
inline Dyadic div_dn(const Dyadic d1, const uint& d2) { return d1/d2; }
inline bool eq_dd(const Dyadic& d1, const Dyadic& d2) { return d1==d2; }
inline bool eq_dx(const Dyadic& d1, const double& d2) { return d1==d2; }
inline bool ne_dd(const Dyadic& d1, const Dyadic& d2) { return d1!=d2; }
inline bool ne_dx(const Dyadic& d1, const double& d2) { return d1!=d2; }
inline bool lt_dd(const Dyadic& d1, const Dyadic& d2) { return d1<d2; }
inline bool lt_dx(const Dyadic& d1, const double& d2) { return d1<d2; }
inline bool gt_dx(const Dyadic& d1, const double& d2) { return d1>d2; }
inline bool le_dd(const Dyadic& d1, const Dyadic& d2) { return d1<=d2; }
inline bool le_dx(const Dyadic& d1, const double& d2) { return d1<=d2; }
inline bool ge_dx(const Dyadic& d1, const double& d2) { return d1>=d2; }
inline Dyadic div_approx_dd(const Dyadic d1, const Dyadic& d2, const Dyadic& e) { return div_approx(d1,d2,e);  }
inline Dyadic div_approx_dx(const Dyadic d1, const double& d2, const Dyadic& e) { return div_approx(d1,Dyadic(d2),e);  }
inline Dyadic div_approx_xd(const double d1, const Dyadic& d2, const Dyadic& e) { return div_approx(Dyadic(d1),d2,e);  }
inline Dyadic div_prec_dd(const Dyadic d1, const Dyadic& d2, const int& n) { return div_approx(d1,d2,n);  }
inline Dyadic mantissa_d(const Dyadic d) { return mantissa(d); }
inline int exponent_d(const Dyadic d) { return exponent(d); }
inline int precision_d(const Dyadic d) { return precision(d); }

inline Rational neg_q(const Rational q1) { return -q1; }
inline Rational add_qq(const Rational q1, const Rational& q2) { return q1+q2; }
inline Rational add_qx(const Rational q1, const double& q2) { return q1+q2; }
inline Rational sub_qq(const Rational q1, const Rational& q2) { return q1-q2; }
inline Rational sub_qx(const Rational q1, const double& q2) { return q1-q2; }
inline Rational rsub_qx(const Rational q1, const double& q2) { return q2-q1; }
inline Rational mul_qq(const Rational q1, const Rational& q2) { return q1*q2; }
inline Rational mul_qx(const Rational q1, const double& q2) { return q1*q2; }
inline Rational div_qq(const Rational q1, const Rational& q2) { return q1/q2; }
inline Rational div_qx(const Rational q1, const double& q2) { return q1/q2; }
inline Rational rdiv_qx(const Rational q1, const double& q2) { return q2/q1; }
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
inline bool eq_dz(const Dyadic& n1, const Integer& n2) { return n1==n2; }
inline bool ne_qd(const Rational& n1, const Dyadic& n2) { return n1!=n2; }
inline bool ne_qz(const Rational& n1, const Integer& n2) { return n1!=n2; }
inline bool ne_dz(const Dyadic& n1, const Integer& n2) { return n1!=n2; }
inline bool lt_qd(const Rational& n1, const Dyadic& n2) { return n1<n2; }
inline bool lt_qz(const Rational& n1, const Integer& n2) { return n1<n2; }
inline bool lt_dz(const Dyadic& n1, const Integer& n2) { return n1<n2; }
inline bool gt_qd(const Rational& n1, const Dyadic& n2) { return n1>n2; }
inline bool gt_qz(const Rational& n1, const Integer& n2) { return n1>n2; }
inline bool gt_dz(const Dyadic& n1, const Integer& n2) { return n1>n2; }
inline bool le_qd(const Rational& n1, const Dyadic& n2) { return n1<=n2; }
inline bool le_qz(const Rational& n1, const Integer& n2) { return n1<=n2; }
inline bool le_dz(const Dyadic& n1, const Integer& n2) { return n1<=n2; }
inline bool ge_qd(const Rational& n1, const Dyadic& n2) { return n1>=n2; }
inline bool ge_qz(const Rational& n1, const Integer& n2) { return n1>=n2; }
inline bool ge_dz(const Dyadic& n1, const Integer& n2) { return n1>=n2; }
inline Integer numerator_q(const Rational& q) { return numerator(q); }
inline Integer denominator_q(const Rational& q) { return denominator(q); }

void export_numeric() {
  class_<Integer>("Integer")
    .def(init<int>())
    .def(init<double>())
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
    .def(str(self))    // __str__
  ;

  class_<Dyadic>("Dyadic")
    .def(init<int>())
    .def(init<double>())
    .def(init<int,int>())
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
    .def("__ge__", &ge_dx)
    .def("__ge__", &ge_dz)
    .def("precision", &precision_d)
    .def("exponent", &exponent_d)
    .def("mantissa", &mantissa_d)
    .def(str(self))    // __str__
  ;

  def("div_approx", &div_approx_dd);
  def("div_approx", &div_approx_dx);
  def("div_approx", &div_approx_xd);
  def("div_prec", &div_prec_dd);

  class_<Rational>("Rational")
    .def(init<int,int>())
    .def(init<int>())
    .def(init<double>())
    .def(init<Integer>())
    .def(init<Dyadic>())
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
    .def(str(self))    // __str__
  ;
}
