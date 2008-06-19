/***************************************************************************
 *            python/export_constraint.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

#include "python/float.h"

#include "function/function_interface.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"

#include "geometry/constraint_interface.h"
#include "geometry/constraint.h"
#include "geometry/linear_constraint.h"


using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
class ConstraintWrapper
  : public ConstraintInterface<R>, 
    public wrapper< ConstraintInterface<R> >
{
  typedef typename traits<R>::arithmetic_type A;
 public:
  ConstraintWrapper() { }
  ConstraintWrapper<R>* clone() const { return this->get_override("clone")(); }
  dimension_type dimension() const { return this->get_override("dimension")(); }
  Comparison comparison() const { return this->get_override("comparison")(); }
  smoothness_type smoothness() const { return this->get_override("smoothness")(); }
  std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
  
  A value(const Point<A>& pt) const { return this->get_override("value")(); }
  Vector<A> gradient(const Point<A>& pt) const { return this->get_override("gradient")(); }
};
  



template<class R>
void export_constraint() 
{
  typedef Interval<R> I;

  enum_<Comparison>("Comparison")
    .value("greater", greater)
    .value("less", less)
  ;

  class_< ConstraintWrapper<R>, boost::noncopyable >
   ("ConstraintInterface",init<>())
    .def("dimension",&ConstraintWrapper<R>::dimension)
    .def("comparison",&ConstraintWrapper<R>::comparison)
    .def("smoothness",&ConstraintWrapper<R>::smoothness)
    .def("value",&ConstraintWrapper<R>::value)
    .def("gradient",&ConstraintWrapper<R>::gradient)
  ;

  class_< LinearConstraint<R>, 
    bases< ConstraintInterface<R> > >
      ("LinearConstraint",init<Vector<R>,R>())
    .def(init< Vector<R>,Comparison,R >())
    .def("dimension",&LinearConstraint<R>::dimension)
    .def("dimension",&LinearConstraint<R>::dimension)
    .def("comparison",&LinearConstraint<R>::comparison)
    .def("smoothness",&LinearConstraint<R>::smoothness)
    .def("value",&LinearConstraint<R>::value)
    .def("gradient",&LinearConstraint<R>::gradient)
    .def(self_ns::str(self))
  ;

  class_< Constraint<R>, 
    bases< ConstraintInterface<R> > >
      ("Constraint",init<const FunctionInterface<R>&>())
    .def("dimension",&Constraint<R>::dimension)
    .def("dimension",&Constraint<R>::dimension)
    .def("comparison",&Constraint<R>::comparison)
    .def("smoothness",&Constraint<R>::smoothness)
    .def("value",&Constraint<R>::value)
    .def("gradient",&Constraint<R>::gradient)
    .def(self_ns::str(self))
  ;


}

template void export_constraint<FloatPy>();
