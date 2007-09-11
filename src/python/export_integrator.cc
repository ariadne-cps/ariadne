/***************************************************************************
 *            python/export_integrator.cc
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

#include "python/python_float.h"

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/grid_set.h"
#include "geometry/list_set.h"

#include "system/vector_field.h"
#include "system/affine_vector_field.h"

#include "evaluation/integrator_interface.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/euler_integrator.h"


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
class IntegratorWrapper
  : public IntegratorInterface<R>,
    public wrapper< IntegratorInterface<R> >
{
  typedef Interval<R> I;
 public:
  IntegratorWrapper() { }
  IntegratorWrapper<R>* clone() const { return this->get_override("clone")(); }
  Point<I> flow_step(const VectorFieldInterface<R>&, const Point<I>&, const I&, const Rectangle<R>&) const {
    return this->get_override("flow_step")(); }
  Rectangle<R> integration_step(const VectorFieldInterface<R>&, const Rectangle<R>&, const I&, const Rectangle<R>&) const {
    return this->get_override("integration_step")(); }
  Rectangle<R> reachability_step(const VectorFieldInterface<R>&, const Rectangle<R>&, const I&, const Rectangle<R>&) const {
    return this->get_override("reachability_step")(); }
  Zonotope<I,R> integration_step(const VectorFieldInterface<R>&, const Zonotope<I,R>&, const I&, const Rectangle<R>&) const {
    return this->get_override("integration_step")(); }
  Zonotope<I,R> reachability_step(const VectorFieldInterface<R>&, const Zonotope<I,R>&, const I&, const Rectangle<R>&) const {
    return this->get_override("reachability_step")(); }
  Zonotope<I,I> integration_step(const VectorFieldInterface<R>&, const Zonotope<I,I>&, const I&, const Rectangle<R>&) const {
    return this->get_override("integration_step")(); }
  Zonotope<I,I> reachability_step(const VectorFieldInterface<R>&, const Zonotope<I,I>&, const I&, const Rectangle<R>&) const {
    return this->get_override("reachability_step")(); }
};

template<class R>
void export_integrator() 
{
  typedef Interval<R> I;

  class_< IntegratorWrapper<R>, boost::noncopyable >("IntegratorInterface",init<>());

  class_< C1LohnerIntegrator<R>, bases<IntegratorInterface<R> > >("C1LohnerIntegrator",init<>());

  class_< LohnerIntegrator<R>, bases<IntegratorInterface<R> > >("LohnerIntegrator",init<>());

  class_< AffineIntegrator<R>, bases<IntegratorInterface<R> > >("AffineIntegrator",init<>());

  class_< EulerIntegrator<R>, bases<IntegratorInterface<R> > >("EulerIntegrator",init<>());
}

template void export_integrator<Float>();
