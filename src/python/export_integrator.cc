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

#include "python/float.h"

#include "geometry/box.h"
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"

#include "system/vector_field.h"
#include "system/affine_vector_field.h"

#include "evaluation/integrator_interface.h"
#include "evaluation/kuhn_integrator.h"
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

template<class Integrator, class R>
tuple
flow_bounds(const Integrator& i, const VectorField<R>& vf, const Box<R>& bb, const Rational& t)
{
  std::pair< Rational, Box<R> > result=i.flow_bounds(vf,bb,t);
  return boost::python::make_tuple(result.first,result.second);
}

template<class BS>
class IntegratorWrapper
  : public IntegratorInterface<BS>,
    public wrapper< IntegratorInterface<BS> >
{
  typedef Numeric::Rational T;
  typedef typename BS::real_type R;
  typedef Interval<R> I;
 public:
  IntegratorWrapper() { }
  IntegratorWrapper<BS>* clone() const { return this->get_override("clone")(); }
  std::pair< T,Box<R> > flow_bounds(const VectorField<R>&, const Box<R>&, const T&) const {
    return this->get_override("flow_bounds")(); }
  BS integration_step(const VectorField<R>&, const BS&, const T&, const Box<R>&) const {
    return this->get_override("integration_step")(); }
  BS reachability_step(const VectorField<R>&, const BS&, const T&, const Box<R>&) const {
    return this->get_override("reachability_step")(); }
  std::ostream& write(std::ostream&) const {
    return this->get_override("write")(); }
};

template<class R>
void export_integrator() 
{
  typedef Interval<R> I;

  class_< IntegratorWrapper< Box<R> >, boost::noncopyable >("BoxIntegratorInterface",init<>());
  class_< IntegratorWrapper< Zonotope<R> >, boost::noncopyable >("ZonotopeIntegratorInterface",init<>());

  class_< KuhnIntegrator<R>, bases<IntegratorInterface< Zonotope<R> > > > kuhn_integrator_class("KuhnIntegrator",init<uint>());
  kuhn_integrator_class.def("flow_bounds",&flow_bounds<KuhnIntegrator<R>,R>);
  kuhn_integrator_class.def("integration_step",&KuhnIntegrator<R>::integration_step);
  kuhn_integrator_class.def("reachability_step",&KuhnIntegrator<R>::reachability_step);

  class_< LohnerIntegrator<R>, bases<IntegratorInterface< Zonotope<R> > > >
    lohner_integrator_class("LohnerIntegrator",init<>());
  lohner_integrator_class.def("flow_bounds",&flow_bounds<LohnerIntegrator<R>,R>);
  lohner_integrator_class.def("integration_step",&LohnerIntegrator<R>::integration_step);
  lohner_integrator_class.def("reachability_step",&LohnerIntegrator<R>::reachability_step);

  class_< AffineIntegrator<R>, bases< IntegratorInterface< Zonotope<R> > > >
    affine_integrator_class("AffineIntegrator",init<>());
  affine_integrator_class.def("flow_bounds",&flow_bounds<AffineIntegrator<R>,R>);
  affine_integrator_class.def("integration_step",(Zonotope<R>(AffineIntegrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,const Rational&,const Box<R>&)const) &AffineIntegrator<R>::integration_step);
  affine_integrator_class.def("reachability_step",(Zonotope<R>(AffineIntegrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,const Rational&,const Box<R>&)const) &AffineIntegrator<R>::reachability_step);

  class_< EulerIntegrator<R>, bases<IntegratorInterface< Box<R> > > >("EulerIntegrator",init<>());
}

template void export_integrator<FloatPy>();
