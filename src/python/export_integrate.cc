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

#include "real_typedef.h"

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/grid_set.h"
#include "geometry/list_set.h"

#include "system/vector_field.h"

#include "evaluation/integrator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"


using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
class IntegratorWrapper : public Integrator<R>, public wrapper< Integrator<R> >
{
 public:
  IntegratorWrapper() : Integrator<R>(0.125,0.125,0.125) { }
  Zonotope<R> integrate(const VectorField<R>& vf,const Zonotope<R>& z,const time_type& t) {
    if (override integrate = this->get_override("integrate")) { return this->integrate(vf,z,t); } 
    else { return this->Integrator<R>::integrate(vf,z,t); } }
  Zonotope<R> reach(const VectorField<R>& vf,const Zonotope<R>& z,const time_type& t) {
    if (override reach = this->get_override("reach")) { return this->reach(vf,z,t); } 
    else { return this->Integrator<R>::reach(vf,z,t); } }
  Zonotope<R> default_integrate(const VectorField<R>& vf,const Zonotope<R>& z,const time_type& t) {
    return this->Integrator<R>::integrate(vf,z,t); }
  
};
  
 
template<class R>
void export_integrate() 
{

 class_< IntegratorWrapper<R>, boost::noncopyable >("Integrator")
    .def("integrate",(Zonotope<R>(Integrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,const time_type&)const)
                              (&Integrator<R>::integrate))
    .def("integrate",(ListSet< Zonotope<R> >(Integrator<R>::*)(const VectorField<R>&,const ListSet< Zonotope<R> >&,const time_type&)const)
                              (&Integrator<R>::integrate))
    .def("reach",(ListSet< Zonotope<R> >(Integrator<R>::*)(const VectorField<R>&,const ListSet< Zonotope<R> >&,const time_type&)const)
                              (&Integrator<R>::reach))
    .def("integrate",(GridMaskSet<R>(Integrator<R>::*)(const VectorField<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                              (&Integrator<R>::integrate))
    .def("reach",(GridMaskSet<R>(Integrator<R>::*)(const VectorField<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                              (&Integrator<R>::reach))
    .def("chainreach",(GridMaskSet<R>(Integrator<R>::*)(const VectorField<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
         (&Integrator<R>::chainreach))
  ;

 class_< LohnerIntegrator<R>, bases<Integrator<R> > >("LohnerIntegrator",init<R,R,R>())
    .def(init<double,double,double>()) 
    .def("integrate",(Zonotope<R>(LohnerIntegrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,const time_type&)const)
                              (&LohnerIntegrator<R>::integrate))
    .def("integrate",(ListSet< Zonotope<R> >(LohnerIntegrator<R>::*)(const VectorField<R>&,const ListSet< Zonotope<R> >&,const time_type&)const)
                              (&LohnerIntegrator<R>::integrate))
    .def("reach",(ListSet< Zonotope<R> >(LohnerIntegrator<R>::*)(const VectorField<R>&,const ListSet< Zonotope<R> >&,const time_type&)const)
                              (&LohnerIntegrator<R>::reach))
  ;

 class_< AffineIntegrator<R>, bases<Integrator<R> > >("AffineIntegrator",init<R,R,R>())
    .def(init<double,double,double>()) 
    .def("integrate",(Zonotope<R>(AffineIntegrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,const time_type&)const)
                              (&AffineIntegrator<R>::integrate))
    .def("integrate",(ListSet< Zonotope<R> >(AffineIntegrator<R>::*)(const VectorField<R>&,const ListSet< Zonotope<R> >&,const time_type&)const)
                              (&AffineIntegrator<R>::integrate))
    .def("reach",(ListSet< Zonotope<R> >(AffineIntegrator<R>::*)(const VectorField<R>&,const ListSet< Zonotope<R> >&,const time_type&)const)
                              (&AffineIntegrator<R>::reach))
  ;

}

template void export_integrate<Real>();
