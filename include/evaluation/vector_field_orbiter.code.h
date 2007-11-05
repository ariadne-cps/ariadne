/***************************************************************************
 *            vector_field_orbiter.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include "vector_field_orbiter.h"

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../combinatoric/lattice_set.h"

#include "../geometry/rectangle.h"
#include "../geometry/zonotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid.h"
#include "../geometry/grid_set.h"
#include "../geometry/partition_tree_set.h"
#include "../geometry/grid_approximation.h"
#include "../geometry/rectangular_set.h"
#include "../geometry/orbit.h"

#include "../system/grid_multimap.h"


#include "../system/vector_field.h"

#include "../evaluation/bounder_interface.h"
#include "../evaluation/integrator_interface.h"
#include "../evaluation/lohner_integrator.h"
#include "../evaluation/bounder.h"

#include "../output/logging.h"

namespace Ariadne {
  
namespace Evaluation { static int& verbosity = applicator_verbosity; }



template<class BS>
Evaluation::VectorFieldOrbiter<BS>::VectorFieldOrbiter(const EvolutionParameters<R>& parameters, const IntegratorInterface<BS>& integrator)
  : _bounder(new Bounder<R>()),
    _integrator(integrator.clone()),
    _parameters(parameters.clone())
{
}


template<class BS>
Evaluation::VectorFieldOrbiter<BS>::VectorFieldOrbiter(const VectorFieldOrbiter<BS>& orbiter)
  : _bounder(orbiter._bounder->clone()),
    _integrator(orbiter._integrator->clone()),
    _parameters(orbiter._parameters->clone())
{
}


template<class BS>
Evaluation::VectorFieldOrbiter<BS>*
Evaluation::VectorFieldOrbiter<BS>::clone() const 
{
  return new VectorFieldOrbiter<BS>(*this);
}





template<class BS>
Geometry::ContinuousTimeOrbit< Numeric::Rational, Geometry::Rectangle<typename BS::real_type>, Geometry::Rectangle<typename BS::real_type> >
Evaluation::VectorFieldOrbiter<BS>::orbit(const System::VectorFieldInterface<R>& vf, const Geometry::Rectangle<R>& r, const Numeric::Rational& t) const
{
  ARIADNE_LOG(4,"ContinuousTimeOrbit<Rational,Rectangle> VectorFieldOrbiter::orbit(VectorFieldInterface,Rectangle,Rational)\n");
  assert(t>=0);
  Geometry::ContinuousTimeOrbit<Numeric::Rational, Geometry::Rectangle<R>, Geometry::Rectangle<R> > orbit(r);
  Geometry::Rectangle<R> bb;
  BS bs(r);
  BS rs(bs);
  const Numeric::Rational& maxh=this->_parameters->maximum_step_size();
  Numeric::Rational h=maxh;
  Numeric::Rational s=0;
  while(s<t) {
    bb=this->_bounder->flow_bounds(vf,bs.bounding_box(),h);
    if(s+h>t) { h=t-s; }
    rs=this->_integrator->reachability_step(vf,bs,h,bb);
    bs=this->_integrator->integration_step(vf,bs,h,bb);
    s=s+h;
    orbit.push_back(s,rs.bounding_box(),bs.bounding_box());
  }
  return orbit;
}



template<class BS>
Geometry::ContinuousTimeOrbit< Numeric::Rational, Geometry::GridCellListSet<typename BS::real_type>, Geometry::GridCellListSet<typename BS::real_type> >
Evaluation::VectorFieldOrbiter<BS>::orbit(const System::VectorFieldInterface<R>& vf, const Geometry::GridCell<R>& gc, const Numeric::Rational& t) const
{
  ARIADNE_LOG(4,"ContinuousTimeOrbit<Rational,Rectangle> VectorFieldOrbiter::orbit(VectorFieldInterface,Rectangle,Rational)\n");
  assert(t>=0);
  const Geometry::Grid<R>& g(gc.grid());
  Geometry::GridCellListSet<R> rgcls(g);
  Geometry::GridCellListSet<R> egcls(g);
  egcls.adjoin(gc);
  Geometry::ContinuousTimeOrbit<Numeric::Rational, Geometry::GridCellListSet<R>, Geometry::GridCellListSet<R> > orbit(egcls);
  Geometry::Rectangle<R> bb;
  BS es(gc);
  BS rs(es);
  const Numeric::Rational& maxh=this->_parameters->maximum_step_size();
  Numeric::Rational h=maxh;
  Numeric::Rational s=0;
  while(s<t) {
    h=maxh;
    bb=this->_bounder->flow_bounds(vf,es.bounding_box(),h);
    if(s+h>t) { h=t-s; }
    rs=this->_integrator->reachability_step(vf,es,h,bb);
    es=this->_integrator->integration_step(vf,es,h,bb);
    rgcls=this->_approximator->outer_approximation(rs,g);
    egcls=this->_approximator->outer_approximation(es,g);
    s=s+h;
    orbit.push_back(s,rgcls,egcls);
  }
  return orbit;
}





}
