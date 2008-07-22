/***************************************************************************
 *            model_checker.code.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *
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
 
#include "model_checker.h"

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
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/rectangular_set.h"
#include "geometry/orbit.h"

#include "system/grid_multimap.h"


#include "system/map.h"
#include "system/discrete_time_system.h"

#include "output/logging.h"

namespace { 

using namespace Ariadne;

template<class R> void ltgt(Integer& n, const EvolutionParameters<R>& params) {
  n=params.lock_to_grid_steps(); }
template<class R> void ltgt(Rational& t, const EvolutionParameters<R>& params) {
  t=params.lock_to_grid_time(); }

template<class T, class R> T lock_to_grid_time(const EvolutionParameters<R>& params) {
  T t; ltgt(t,params); return t; }

}

namespace Ariadne {


template<class T, class Aprx>
tribool
ModelChecker<T,Aprx>::verify(const TSI& f, 
                                         const PTS& initial_set, 
                                         const PTS& safe_set) const
{
  int verbosity=this->verbosity();
  ARIADNE_LOG(2,"triboolModelChecker::verify(TransitionSystemInterface map, GridMaskSet initial_set,GridMaskSet  safe_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n"<<"safe_set="<<safe_set);
  
  ARIADNE_CHECK_BOUNDED(initial_set,"ModelChecker<R>::verify(...)");
  ARIADNE_CHECK_BOUNDED(safe_set,"ModelChecker<R>::verify(...)");
  
  T time = lock_to_grid_time<T>(this->parameters());
  const Pv& g=initial_set.grid();
  LatticeBlock bd=safe_set.block();
  Bx bb=safe_set.bounding_box();
  PTS reach(g,bd);
  PLS found(g);
  PLS cell_image(g);
  PLS image(g);
  
  Box<R> r(g.dimension());
  
  found=initial_set;
  while(!subset(found,reach)) {
    found=difference(found,reach);
    reach.adjoin(found);
    image.clear();
    for(typename PLS::const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
      cell_image=f.upper_reach(*iter,time);
      if(!subset(cell_image,safe_set)) {
        return false;
      }
      image.adjoin(cell_image);
    }
    found=image;
  }
  return true;
}



template<class T, class Aprx>
tribool
ModelChecker<T,Aprx>::verify(const TransitionSystem& system, 
                                         const PartitionTreeSet& initial_set, 
                                         const TimedLogicFormula& formula) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}









}
