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

/*
template<class R>
class LohnerIntegratorWrapper : public LohnerIntegrator<R>
{
 public:
  LohnerIntegratorWrapper(const R& mst,const R& lgt, const R& mss) : LohnerIntegrator<R>(mst,lgt,mss) { }
  Zonotope<R> integrate(const VectorField<R>& vf,const Zonotope<R>& z,const time_type& t) const {
    return this->LohnerIntegrator<R>::integrate(vf,z,t);
  }
};
*/

template<class R>
class LohnerIntegratorWrapper
{
 public:
  LohnerIntegratorWrapper(const R& mst,const R& lgt, const R& mss) : _integrator(mst,lgt,mss) { }
  Zonotope<R> integrate(const VectorField<R>& vf,const Zonotope<R>& z,const time_type& t) const {
    Zonotope<R> fz=this->_integrator.integrate(vf,z,t);
    std::cerr << "LohnerIntegratorWrapper::integrate(...): result=" << fz << std::endl;
    return fz;
  }
  ListSet<R,Zonotope> reach(const VectorField<R>& vf,const Zonotope<R>& z,const time_type& t) const {
    ListSet<R,Zonotope> fz=this->_integrator.reach(vf,z,t);
    std::cerr << "LohnerIntegratorWrapper::reach(...): result=" << fz << std::endl;
    return fz;
  }
  Zonotope<R> simple_identity(const Zonotope<R>& z) const { return z; }
  Zonotope<R> identity(const Zonotope<R>& z) const { return _integrator.identity(z); }
  Vector< Interval<R> > field(const VectorField<R>& vf, const Zonotope<R>& z) const { return _integrator.field(vf,z); }

 private:
  LohnerIntegrator<R> _integrator;
};
 
template<class R>
class Identity 
{
 public:
  Identity() { }
  Rectangle<R> image(const Rectangle<R>& r) const { return r; }
  Zonotope<R> image(const Zonotope<R>& z) const { return z; }
};

  
 
template<class R>
void export_integrate() 
{

/*  class_< IntegratorWrapper<R>, boost::noncopyable>("Integrator")
    .def_readwrite("maximum_step_size", &Integrator<R>::maximum_step_size)
    .def_readwrite("maximum_basic_set_radius", &Integrator<R>::maximum_basic_set_radius)
    .def("integrate",(Rectangle<R>(Integrator<R>::*)(const VectorField<R>&,const Rectangle<R>&,const time_type&)const)
                              (&Integrator<R>::integrate))
    .def("integrate",(Parallelotope<R>(Integrator<R>::*)(const VectorField<R>&,const Parallelotope<R>&,const time_type&)const)
                              (&Integrator<R>::integrate))

    .def("integrate",
         (Zonotope<R>(Integrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,const time_type&)const)(&Integrator<R>::integrate),
         (Zonotope<R>(IntegratorWrapper<R>::*)(const VectorField<R>&,const Zonotope<R>&,const time_type&)const)(&IntegratorWrapper<R>::default_integrate))
  
    .def("integrate",(ListSet<R,Parallelotope>(Integrator<R>::*)(const VectorField<R>&,const ListSet<R,Parallelotope>&,const time_type&)const)
                              (&Integrator<R>::integrate))
    .def("integrate",(ListSet<R,Zonotope>(Integrator<R>::*)(const VectorField<R>&,const ListSet<R,Zonotope>&,const time_type&)const)
                              (&Integrator<R>::integrate))
    .def("integrate",(GridMaskSet<R>(Integrator<R>::*)(const VectorField<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                              (&Integrator<R>::integrate))
    .def("reach",(ListSet<R,Zonotope>(Integrator<R>::*)(const VectorField<R>&,const ListSet<R,Parallelotope>&,const time_type&)const)
                              (&Integrator<R>::reach))
    .def("reach",(ListSet<R,Zonotope>(Integrator<R>::*)(const VectorField<R>&,const ListSet<R,Zonotope>&,const time_type&)const)
                              (&Integrator<R>::reach))
    .def("reach",(GridMaskSet<R>(Integrator<R>::*)(const VectorField<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&,const time_type&)const)
                              (&Integrator<R>::reach))

    .def("chainreach",(GridMaskSet<R>(Integrator<R>::*)(const VectorField<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
                              (&Integrator<R>::chainreach),"chain reach of a set")
    .def("verify",(bool(Integrator<R>::*)(const VectorField<R>&,const GridMaskSet<R>&,const GridMaskSet<R>&)const)
                              (&Integrator<R>::verify),"verify a safety property")
  ;
*/
/*
  class_< LohnerIntegrator<R>, bases< Integrator<R> > >("LohnerIntegrator",init<R,R,R>())
    .def(init<double,double,double>()) 
    .def("integration_step",(Parallelotope<R>(LohnerIntegrator<R>::*)(const VectorField<R>&,const Parallelotope<R>&,time_type&)const)
                              (&LohnerIntegrator<R>::integration_step))
    .def("integration_step",(Zonotope<R>(LohnerIntegrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,time_type&)const)
                              (&LohnerIntegrator<R>::integration_step))
    .def("reach_step",(Zonotope<R>(LohnerIntegrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,time_type&)const)
                              (&LohnerIntegrator<R>::reachability_step))
  ;
*/

 class_< Identity<R> >("Identity",init<>())
    .def("image",(Rectangle<R>(Identity<R>::*)(const Rectangle<R>&)const)(&Identity<R>::image))
    .def("image",(Zonotope<R>(Identity<R>::*)(const Zonotope<R>&)const)(&Identity<R>::image))
  ;

  class_< LohnerIntegrator<R> >("LohnerIntegrator",init<R,R,R>())
    .def(init<double,double,double>()) 
    .def("integrate",(Parallelotope<R>(Integrator<R>::*)(const VectorField<R>&,const Parallelotope<R>&,const time_type&)const)
                              (&Integrator<R>::integrate))
    .def("integrate",(Zonotope<R>(Integrator<R>::*)(const VectorField<R>&,const Zonotope<R>&,const time_type&)const)
                              (&Integrator<R>::integrate))
    .def("integrate",(ListSet<R,Parallelotope>(Integrator<R>::*)(const VectorField<R>&,const ListSet<R,Parallelotope>&,const time_type&)const)
                              (&Integrator<R>::integrate))
    .def("integrate",(ListSet<R,Zonotope>(Integrator<R>::*)(const VectorField<R>&,const ListSet<R,Zonotope>&,const time_type&)const)
                              (&Integrator<R>::integrate))
    .def("reach",(ListSet<R,Zonotope>(Integrator<R>::*)(const VectorField<R>&,const ListSet<R,Zonotope>&,const time_type&)const)
                              (&Integrator<R>::reach))
  ;



/*
class_< LohnerIntegratorWrapper<R> >("LohnerIntegrator",init<R,R,R>())
    .def(init<double,double,double>()) 
    .def("integrate",&LohnerIntegratorWrapper<R>::integrate)
    .def("reach",&LohnerIntegratorWrapper<R>::reach)
  ;
*/
  
}

template void export_integrate<Real>();
