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

template<class BS>
class IntegratorWrapper
  : public IntegratorInterface<BS>,
    public wrapper< IntegratorInterface<BS> >
{
  typedef typename BS::real_type R;
  typedef Interval<R> I;
 public:
  IntegratorWrapper() { }
  IntegratorWrapper<BS>* clone() const { return this->get_override("clone")(); }
  BS integration_step(const VectorFieldInterface<R>&, const BS&, const I&, const Rectangle<R>&) const {
    return this->get_override("integration_step")(); }
  BS reachability_step(const VectorFieldInterface<R>&, const BS&, const I&, const Rectangle<R>&) const {
    return this->get_override("reachability_step")(); }
  std::ostream& write(std::ostream&) const {
    return this->get_override("write")(); }
};

template<class R>
void export_integrator() 
{
  typedef Interval<R> I;

  class_< IntegratorWrapper< Rectangle<R> >, boost::noncopyable >("RectangleIntegratorInterface",init<>());
  class_< IntegratorWrapper< Zonotope<I,R> >, boost::noncopyable >("C0ZonotopeIntegratorInterface",init<>());
  class_< IntegratorWrapper< Zonotope<I,I> >, boost::noncopyable >("C1ZonotopeIntegratorInterface",init<>());

  class_< LohnerIntegrator<R>, bases<IntegratorInterface< Zonotope<I,R> > > >("C0LohnerIntegrator",init<>());
  class_< C1LohnerIntegrator<R>, bases<IntegratorInterface< Zonotope<I,I> > > >("C1LohnerIntegrator",init<>());


  //  class_< AffineIntegrator<R>, bases<IntegratorInterface< Zonotope<I,R> > > >("C0AffineIntegrator",init<>());
  //  class_< AffineIntegrator<R>, bases<IntegratorInterface< Zonotope<I,I> > > >("C1AffineIntegrator",init<>());
  class_< AffineIntegrator<R>, bases< IntegratorInterface< Zonotope<I,R> >, IntegratorInterface< Zonotope<I,I> > > >("AffineIntegrator",init<>());

  class_< EulerIntegrator<R>, bases<IntegratorInterface< Rectangle<R> > > >("EulerIntegrator",init<>());
}

template void export_integrator<Float>();
