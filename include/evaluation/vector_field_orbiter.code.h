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

#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "combinatoric/lattice_set.h"

#include "geometry/box.h"
#include "geometry/box_list_set.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/grid_approximation.h"
#include "geometry/orbit.h"

#include "system/grid_multimap.h"


#include "system/vector_field.h"

#include "evaluation/bounder_interface.h"
#include "evaluation/integrator_interface.h"

#include "evaluation/standard_bounder.h"
#include "evaluation/lohner_integrator.h"

#include "output/logging.h"

namespace Ariadne {
  
namespace Evaluation { static int& verbosity = applicator_verbosity; }



template<class BS>
Evaluation::VectorFieldOrbiter<BS>::VectorFieldOrbiter(const EvolutionParameters<R>& parameters, const IntegratorInterface<BS>& integrator, const ApproximatorInterface<BS>& approximator)
  : _bounder(new StandardBounder<R>()),
    _integrator(integrator.clone()),
    _approximator(approximator.clone()),
    _parameters(parameters.clone())
{
}


template<class BS>
Evaluation::VectorFieldOrbiter<BS>*
Evaluation::VectorFieldOrbiter<BS>::clone() const 
{
  return new VectorFieldOrbiter<BS>(*this);
}





template<class BS>
Geometry::GridCellListSet<typename BS::real_type> 
Evaluation::VectorFieldOrbiter<BS>::upper_evolve(const System::VectorField<R>& f, 
                                                 const Geometry::GridCell<R>& gc, 
                                                 const Numeric::Rational& t) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  /*
  const Geometry::Grid<R>& grid=gc.grid();
  BS bs=this->over_approximation(Geometry::Box<R>(gc));
  std::vector<BSL> orbit=this->_orbit(f,bs,n);
  GCLS result=this->outer_approximation(orbit[n],grid);
  return result;
  */
}

template<class BS>
Geometry::GridCellListSet<typename BS::real_type> 
Evaluation::VectorFieldOrbiter<BS>::upper_reach(const System::VectorField<R>& f, 
                                                const Geometry::GridCell<R>& gc, 
                                                const Numeric::Rational& t) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  /*
  const Geometry::Grid<R>& grid=gc.grid();
  BS bs=this->over_approximation(Geometry::Box<R>(gc));
  std::vector<BSL> orbit=this->_orbit(f,bs,n);
  GCLS result(grid);
  for(size_type i=0; i!=orbit.size(); ++i) {
    result.adjoin(this->outer_approximation(orbit[n],grid));
  }
  result.unique_sort();
  return result;
  */
}


template<class BS>
Geometry::Box<typename BS::real_type>
Evaluation::VectorFieldOrbiter<BS>::lower_evolve(const System::VectorField<R>& f, 
                                         const Geometry::Box<R>& bx, 
                                         const Numeric::Rational& t) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  /*
  BS bs=this->over_approximation(bx);
  std::vector<BSL> orbit=this->_orbit(f,bs,n);
  return this->bounding_box(orbit[n]);
  */
}


template<class BS>
Geometry::BoxListSet<typename BS::real_type>
Evaluation::VectorFieldOrbiter<BS>::lower_reach(const System::VectorField<R>& f, 
                                        const Geometry::Box<R>& bx, 
                                        const Numeric::Rational& t) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
  /*
  BS bs=this->over_approximation(bx);
  std::vector<BSL> orbit=this->_orbit(f,bs,n);
  Geometry::BoxListSet<R> result;
  for(size_type i=0; i!=orbit.size(); ++i) {
    result.adjoin(this->bounding_box(orbit[i]));
  }
  return result;
  */
} 



/*

template<class BS>
Geometry::Orbit<Numeric::Rational,BS,BS>*
Evaluation::VectorFieldOrbiter<BS>::orbit(const System::VectorField<R>& vf, 
                                          const Geometry::Box<R>& is, 
                                          const Numeric::Rational& t) const
{
  ARIADNE_LOG(4,"Orbit<Rational,Box> VectorFieldOrbiter::orbit(VectorField,Box,Rational)\n");
  assert(t>=0);
  BS es(is);
  BS rs(es);
  Geometry::Orbit<Numeric::Rational,BS,BS>* orbit=new Geometry::Orbit<Numeric::Rational,BS,BS>(es);
  Geometry::Box<R> bb;
  const Numeric::Rational& maxh=this->_parameters->maximum_step_size();
  Numeric::Rational h=maxh;
  Numeric::Rational s=0;
  while(s<t) {
    bb=this->_bounder->flow_bounds(vf,es.bounding_box(),h);
    if(s+h>t) { h=t-s; }
    rs=this->_integrator->reachability_step(vf,es,h,bb);
    es=this->_integrator->integration_step(vf,es,h,bb);
    s=s+h;
    orbit->push_back(s,rs,es);
  }
  return orbit;
}

*/

}
