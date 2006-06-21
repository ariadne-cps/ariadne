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

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/zonotope.h"
#include "geometry/grid_set.h"
#include "geometry/list_set.h"

#include "system/vector_field.h"

#include "evaluation/integrator.h"
#include "evaluation/lohner_integrator.h"


#include "python/typedefs.h"
using namespace Ariadne;
using namespace Ariadne::Evaluation;

#include <boost/python.hpp>
using namespace boost::python;

typedef Integrator<Real> RIntegrator;
typedef C0Integrator<Real> RC0Integrator;
typedef C1Integrator<Real> RC1Integrator;
typedef C1LohnerIntegrator<Real> RC1LohnerIntegrator;

void export_integrate() {
  typedef RRectangle (RC1LohnerIntegrator::*IntStepRectFunc) (const RVectorFieldBase&, const RRectangle&, Real&) const;
  typedef RParallelotope (RC1LohnerIntegrator::*IntStepPltpFunc) (const RVectorFieldBase&, const RParallelotope&, Real&) const;
  typedef RRectangle (RC1Integrator::*RchStepRectFunc) (const RVectorFieldBase&, const RRectangle&, Real&) const;
  typedef RZonotope (RC1Integrator::*RchStepPltpFunc) (const RVectorFieldBase&, const RParallelotope&, Real&) const;
  typedef RZonotope (RC1Integrator::*RchStepZntpFunc) (const RVectorFieldBase&, const RZonotope&, Real&) const;
  typedef RRectangle (RC1LohnerIntegrator::*IntRectFunc) (const RVectorFieldBase&, const RRectangle&, const Real&) const;
  typedef RParallelotope (RC1LohnerIntegrator::*IntPltpFunc) (const RVectorFieldBase&, const RParallelotope&, const Real&) const;
  typedef RZonotope (RC1LohnerIntegrator::*IntZltzFunc) (const RVectorFieldBase&, const RZonotope&, const Real&) const;
  typedef RRectangleListSet (RC1LohnerIntegrator::*IntLSRectFunc) (const RVectorFieldBase&, const RRectangleListSet&, const Real&) const;
  typedef RZonotopeListSet (RC1LohnerIntegrator::*IntLSZNtpFunc) (const RVectorFieldBase&, const RZonotopeListSet&, const Real&) const;
  typedef RParallelotopeListSet (RC1LohnerIntegrator::*ReachLSPltpFunc) (const RVectorFieldBase&, const RParallelotopeListSet&, const Real&) const;
  typedef RZonotopeListSet (RC1LohnerIntegrator::*ReachLSPltzFunc) (const RVectorFieldBase&, const RZonotopeListSet&, const Real&) const;
  typedef RParallelotopeListSet (RC1LohnerIntegrator::*ReachLSPpFunc) (const RVectorFieldBase&, const RParallelotope&, const Real&) const;
  typedef RZonotopeListSet (RC1LohnerIntegrator::*ReachLSZzFunc) (const RVectorFieldBase&, const RZonotope&, const Real&) const;
  typedef RParallelotopeListSet (RC1LohnerIntegrator::*IntLSPltpFunc) (const RVectorFieldBase&, const RParallelotopeListSet&, const Real&) const;
  typedef RZonotopeListSet (RC1LohnerIntegrator::*IntLSZltzFunc) (const RVectorFieldBase&, const RZonotopeListSet&, const Real&) const;
  typedef RGridMaskSet (RC1Integrator::*IntGMSFunc) (const RVectorFieldBase&, const RGridMaskSet&, const RGridMaskSet&, const Real&) const;
  typedef RGridMaskSet (RC1LohnerIntegrator::*CRGMSFunc) (const RVectorFieldBase&, const RGridMaskSet&, const RGridMaskSet&) const;
 
  class_<RC1LohnerIntegrator>("C1LohnerIntegrator",init<Real,Real,Real>())
    .def(init<double,double,double>()) 
    .def_readwrite("maximum_step_size", &RC1LohnerIntegrator::maximum_step_size)
    .def_readwrite("maximum_basic_set_radius", &RC1LohnerIntegrator::maximum_basic_set_radius)
    .def("integration_step", IntStepRectFunc(&RC1Integrator::integration_step))
    .def("integration_step", IntStepPltpFunc(&RC1LohnerIntegrator::integration_step))
    .def("reach_step", RchStepPltpFunc(&RC1LohnerIntegrator::reachability_step))
    .def("reach_step", RchStepZntpFunc(&RC1LohnerIntegrator::reachability_step))
    .def("integrate", IntRectFunc(&RC1LohnerIntegrator::integrate))
    .def("integrate", IntPltpFunc(&RC1LohnerIntegrator::integrate))
    .def("integrate", IntZltzFunc(&RC1LohnerIntegrator::integrate))
    .def("integrate", IntLSPltpFunc(&RC1LohnerIntegrator::integrate))
    .def("integrate", IntLSZltzFunc(&RC1LohnerIntegrator::integrate))
    .def("integrate", IntGMSFunc(&RC1Integrator::integrate))
    .def("reach", ReachLSPltpFunc(&RC1LohnerIntegrator::reach))
    .def("reach", ReachLSPltzFunc(&RC1LohnerIntegrator::reach))
    .def("reach", ReachLSPpFunc(&RC1LohnerIntegrator::reach))
    .def("reach", ReachLSZzFunc(&RC1LohnerIntegrator::reach))
    .def("reach", IntGMSFunc(&RC1Integrator::reach))
    .def("chainreach", CRGMSFunc(&RC1LohnerIntegrator::chainreach), "chain reach of a set" )
    ;

}
