/***************************************************************************
 *            python/export_detector.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"

#include "system/vector_field_interface.h"

#include "evaluation/integrator_interface.h"
#include "evaluation/detector_interface.h"
#include "evaluation/detector.h"


using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
class DetectorWrapper : public DetectorInterface<R>, public wrapper< DetectorInterface<R> >
{
  typedef Interval<R> I;
 public:
  DetectorWrapper() : DetectorInterface<R>() { }
  DetectorWrapper<R>* clone() const { 
    return this->get_override("clone")(); }
  Interval<R> value(const ConstraintInterface<R>& c, const Rectangle<R>& r) const {
    return this->get_override("value")(); }
  tribool forces(const ConstraintInterface<R>& c1, const ConstraintInterface<R>& c2, const Rectangle<R>& r) const {
    return this->get_override("forces")(); }
  Interval<R> normal_derivative(const VectorFieldInterface<R>& vf, const DifferentiableConstraintInterface<R>& c, const Point<I>& r) const {
    return this->get_override("normal_derivative")(); }
  Interval<R> crossing_time(const VectorFieldInterface<R>& vf, const ConstraintInterface<R>& c, const Point<I>& pt, const Rectangle<R>& bb) const {
    return this->get_override("crossing_time")(); }
  TimeModel<R> crossing_time(const VectorFieldInterface<R>& vf, const ConstraintInterface<R>& c, const Rectangle<R>& pt, const Rectangle<R>& bb) const {
    return this->get_override("crossing_time")(); }
};
  
 
template<class R>
void export_detector() 
{
  typedef Interval<R> I;

  class_< DetectorWrapper<R>, boost::noncopyable >("DetectorInterface",init<>());

  class_< Detector<R>, bases<DetectorInterface<R> > > detector_class("Detector",init<>());
  detector_class.def("satisfies", (tribool(Detector<R>::*)(const Rectangle<R>&,const ConstraintInterface<R>&)const) &Detector<R>::satisfies);
  detector_class.def("satisfies", (tribool(Detector<R>::*)(const Zonotope<R,UniformErrorTag>&,const ConstraintInterface<R>&)const) &Detector<R>::satisfies);
  detector_class.def("satisfies", (tribool(Detector<R>::*)(const Zonotope<R,IntervalTag>&,const ConstraintInterface<R>&)const) &Detector<R>::satisfies);
  detector_class.def("value", (Interval<R>(Detector<R>::*)(const ConstraintInterface<R>&,const Rectangle<R>&)const) &Detector<R>::value);
  detector_class.def("value", (Interval<R>(Detector<R>::*)(const ConstraintInterface<R>&,const Zonotope<R,UniformErrorTag>&)const) &Detector<R>::value);
  detector_class.def("value", (Interval<R>(Detector<R>::*)(const ConstraintInterface<R>&,const Zonotope<R,IntervalTag>&)const) &Detector<R>::value);
  detector_class.def("forces", (tribool(Detector<R>::*)(const ConstraintInterface<R>&,const ConstraintInterface<R>&,const Rectangle<R>&)const) &Detector<R>::forces);
  detector_class.def("normal_derivative", (Interval<R>(Detector<R>::*)(const VectorFieldInterface<R>&,const DifferentiableConstraintInterface<R>&,const Point<I>&)const) &Detector<R>::normal_derivative);
  detector_class.def("crossing_time", (Interval<R>(Detector<R>::*)(const VectorFieldInterface<R>&,const ConstraintInterface<R>&,const Point<I>&,const Rectangle<R>&)const) &Detector<R>::crossing_time);
  detector_class.def("crossing_time", (TimeModel<R>(Detector<R>::*)(const VectorFieldInterface<R>&,const ConstraintInterface<R>&,const Rectangle<R>&,const Rectangle<R>&)const) &Detector<R>::crossing_time);

}

template void export_detector<FloatPy>();
