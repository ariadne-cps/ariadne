/***************************************************************************
 *            python/ariadne.cc
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

#include "numerical_type.h"
#include "interval.h"
#include "state.h"
#include "rectangle.h"

#include <boost/python.hpp>

using namespace Ariadne;
using namespace Ariadne::Geometry;

typedef Interval<Rational> QInterval;
typedef State<Rational> QState;
typedef Rectangle<Rational> QRectangle;

// Since gmpxx operators return intermediate types which need conversion, need to explicitly provide arithmetic functions
inline Integer neg_z(const Integer& z1) { return -z1; }
inline Integer add_z(const Integer& z1, const Integer& z2) { return z1+z2; }
inline Integer sub_z(const Integer& z1, const Integer& z2) { return z1-z2; }
inline Integer mul_z(const Integer& z1, const Integer& z2) { return z1*z2; }
inline Rational div_z(const Integer& z1, const Integer& z2) { return Rational(z1)/z2; }
inline bool eq_z(const Integer& z1, const Integer& z2) { return z1==z2; }
inline bool ne_z(const Integer& z1, const Integer& z2) { return z1!=z2; }
inline bool lt_z(const Integer& z1, const Integer& z2) { return z1<z2; }
inline bool le_z(const Integer& z1, const Integer& z2) { return z1<=z2; }

inline Dyadic neg_d(const Dyadic d1) { return -d1; }
inline Dyadic add_d(const Dyadic d1, const Dyadic& d2) { return d1+d2; }
inline Dyadic sub_d(const Dyadic d1, const Dyadic& d2) { return d1-d2; }
inline Dyadic mul_d(const Dyadic d1, const Dyadic& d2) { return d1*d2; }
inline Dyadic div_d(const Dyadic d1, const Dyadic& d2) { return d1/d2; }
inline bool eq_d(const Dyadic& d1, const Dyadic& d2) { return d1==d2; }
inline bool ne_d(const Dyadic& d1, const Dyadic& d2) { return d1!=d2; }
inline bool lt_d(const Dyadic& d1, const Dyadic& d2) { return d1<d2; }
inline bool le_d(const Dyadic& d1, const Dyadic& d2) { return d1<=d2; }

inline Rational neg_q(const Rational q1) { return -q1; }
inline Rational add_q(const Rational q1, const Rational& q2) { return q1+q2; }
inline Rational sub_q(const Rational q1, const Rational& q2) { return q1-q2; }
inline Rational mul_q(const Rational q1, const Rational& q2) { return q1*q2; }
inline Rational div_q(const Rational q1, const Rational& q2) { return q1/q2; }
inline bool eq_q(const Rational& q1, const Rational& q2) { return q1==q2; }
inline bool ne_q(const Rational& q1, const Rational& q2) { return q1!=q2; }
inline bool lt_q(const Rational& q1, const Rational& q2) { return q1<q2; }
inline bool le_q(const Rational& q1, const Rational& q2) { return q1<=q2; }



BOOST_PYTHON_MODULE(ariadne)
{
  using boost::python::class_;
  using boost::python::init;
  using boost::python::self;
  using boost::python::return_value_policy;
  using boost::python::copy_const_reference;
  using boost::python::def;
  
  typedef bool (*QRectFun) (const QRectangle&, const QRectangle&);
  QRectFun qinner_subset=&inner_subset<Rational>;

  def("regular_intersection", &regular_intersection<Rational>);
  def("interiors_intersect", &interiors_intersect<Rational>);
  def("disjoint", &disjoint<Rational>);
  def("inner_subset", qinner_subset);
  def("subset", &subset<Rational>);

  class_<Integer>("Integer")
    .def(init<int>())
    .def(init<double>())
    .def("__neg__", &neg_z)
    .def("__add__", &add_z)
    .def("__sub__", &sub_z)
    .def("__mul__", &mul_z)
    .def("__div__", &div_z)
    .def("__eq__", &eq_z)
    .def("__ne__", &ne_z)
    .def("__lt__", &lt_z)
    .def("__le__", &le_z)
    .def(boost::python::self_ns::str(self))    // __str__
    ;

  class_<Dyadic>("Dyadic")
    .def(init<int>())
    .def(init<double>())
    .def(init<Integer>())
    .def("__neg__", &neg_d)
    .def("__add__", &add_d)
    .def("__sub__", &sub_d)
    .def("__mul__", &mul_d)
    .def("__div__", &div_d)
    .def("__eq__", &eq_d)
    .def("__ne__", &ne_d)
    .def("__lt__", &lt_d)
    .def("__le__", &le_d)
    .def(boost::python::self_ns::str(self))    // __str__
    ;

  class_<Rational>("Rational")
    .def(init<int,int>())
    .def(init<int>())
    .def(init<double>())
    .def(init<Integer>())
    .def(init<Dyadic>())
    .def("__neg__", &neg_q)
    .def("__add__", &add_q)
    .def("__sub__", &sub_q)
    .def("__mul__", &mul_q)
    .def("__div__", &div_q)
    .def("__eq__", &eq_q)
    .def("__ne__", &ne_q)
    .def("__lt__", &lt_q)
    .def("__le__", &le_q)
    .def(boost::python::self_ns::str(self))    // __str__
    ;

  class_< Interval<Rational> >("QInterval")
    .def(init<Rational,Rational>())
    .def("lower", &Interval<Rational>::lower, return_value_policy<copy_const_reference>())
    .def("upper", &Interval<Rational>::upper, return_value_policy<copy_const_reference>())
    .def(boost::python::self_ns::str(self))    // __str__
    ;

  class_< State<Rational> >("QState")
    .def(init<int>())
    .def(init<int,Rational>())
    .def(init<QState>())
    .def("dimension", &State<Rational>::dimension)
    .def("__len__", &State<Rational>::dimension)
    .def("__getitem__", &State<Rational>::get)
    .def("__setitem__", &State <Rational>::set)
    .def("__eq__", &State<Rational>::operator==)
    .def("__ne__", &State<Rational>::operator!=)
    .def(boost::python::self_ns::str(self))    // __str__
    ;

  class_< QRectangle >("QRectangle")
    .def(init<int>())
    .def(init<QState,QState>())
    .def(init<QRectangle>())
    .def("dimension", &QRectangle::dimension)
    .def("lower_corner", &QRectangle::lower_corner)
    .def("upper_corner", &QRectangle::upper_corner)
    .def("set_lower", &QRectangle::set_lower)
    .def("set_upper", &QRectangle::set_upper)
    .def("__len__", &QRectangle::dimension)
    .def("__getitem__", &QRectangle::get)
    .def("__setitem__", &QRectangle::set)
    .def("__eq__", &QRectangle::operator==)
    .def("__ne__", &QRectangle::operator!=)
    .def(boost::python::self_ns::str(self))    // __str__
    ;


}

