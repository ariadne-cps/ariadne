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

#include "evaluation/evolution_parameters.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/integrator_plugin_interface.h"
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
class IntegratorPluginWrapper
  : public IntegratorPluginInterface<BS>,
    public wrapper< IntegratorPluginInterface<BS> >
{
  typedef typename BS::real_type R;
  typedef Interval<R> I;
 public:
  IntegratorPluginWrapper() { }
  IntegratorPluginWrapper<BS>* clone() const { return this->get_override("clone")(); }
  Point<I> flow_step(const VectorFieldInterface<R>&, const Point<I>&, const I&, const Rectangle<R>&) const {
    return this->get_override("flow_step")(); }
  BS integration_step(const VectorFieldInterface<R>&, const BS&, const I&, const Rectangle<R>&) const {
    return this->get_override("integration_step")(); }
  BS reachability_step(const VectorFieldInterface<R>&, const BS&, const I&, const Rectangle<R>&) const {
    return this->get_override("reachability_step")(); }
};

template<class R>
void export_integrate() 
{
  typedef time_type T;
  typedef Interval<R> I;
  typedef Zonotope<I,R> BS;

  /*
   class_< IntegratorWrapper<R>, boost::noncopyable >("VectorFieldEvolver",init<T,T,R>())
    .def("integrate",(ListSet< Rectangle<R> >(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const ListSet< Rectangle<R> >&,const time_type&)const)
                                    (&VectorFieldEvolver<R>::integrate))
    .def("reach",(ListSet< Rectangle<R> >(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const ListSet< Rectangle<R> >&,const time_type&)const)
                                    (&VectorFieldEvolver<R>::reach))
    .def("integrate",(GridMaskSet<R>(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                                    (&VectorFieldEvolver<R>::integrate))
    .def("reach",(GridMaskSet<R>(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                              (&VectorFieldEvolver<R>::reach))
    .def("chainreach",(GridMaskSet<R>(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
         (&VectorFieldEvolver<R>::chainreach))
    .def("viable",(GridMaskSet<R>(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&)const)
         (&VectorFieldEvolver<R>::viable))
    .def("verify",(tribool(VectorFieldEvolver<R>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
         (&VectorFieldEvolver<R>::verify))
  ;
  */




  class_< VectorFieldEvolver<BS> >("VectorFieldEvolver",init<const EvolutionParameters<R>&,const IntegratorPluginInterface<BS>&>())
    .def("integrate",(ListSet< Rectangle<R> >(VectorFieldEvolver<BS>::*)(const VectorFieldInterface<R>&,const ListSet< Rectangle<R> >&,const time_type&)const)
                                    (&VectorFieldEvolver<BS>::integrate))
    .def("reach",(ListSet< Rectangle<R> >(VectorFieldEvolver<BS>::*)(const VectorFieldInterface<R>&,const ListSet< Rectangle<R> >&,const time_type&)const)
                                    (&VectorFieldEvolver<BS>::reach))
    .def("integrate",(GridMaskSet<R>(VectorFieldEvolver<BS>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                                    (&VectorFieldEvolver<BS>::integrate))
    .def("reach",(GridMaskSet<R>(VectorFieldEvolver<BS>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                              (&VectorFieldEvolver<BS>::reach))
    .def("chainreach",(GridMaskSet<R>(VectorFieldEvolver<BS>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
         (&VectorFieldEvolver<BS>::chainreach))
    .def("viable",(GridMaskSet<R>(VectorFieldEvolver<BS>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&)const)
         (&VectorFieldEvolver<BS>::viable))
    .def("verify",(tribool(VectorFieldEvolver<BS>::*)(const VectorFieldInterface<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
         (&VectorFieldEvolver<BS>::verify))
  ;


  class_< IntegratorPluginWrapper< Rectangle<R> >, boost::noncopyable >("IntegratorPluginInterface",init<>());
  class_< IntegratorPluginWrapper< Zonotope<I,R> >, boost::noncopyable >("IntegratorPluginInterface",init<>());
  class_< IntegratorPluginWrapper< Zonotope<I,I> >, boost::noncopyable >("IntegratorPluginInterface",init<>());

  class_< C1LohnerIntegrator<R>, bases<IntegratorPluginInterface< Zonotope<I,I> > > >("C1LohnerIntegrator",init<>());

  class_< LohnerIntegrator<R>, bases<IntegratorPluginInterface< Zonotope<I,R> > > >("LohnerIntegrator",init<>());

  class_< AffineIntegrator<R>, bases<IntegratorPluginInterface< Zonotope<I,I> > > >("AffineIntegrator",init<>());

  class_< EulerIntegrator<R>, bases<IntegratorPluginInterface< Rectangle<R> > > >("EulerIntegrator",init<>());
}

template void export_integrate<Float>();
