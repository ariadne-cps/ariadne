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
#include "geometry/grid_approximation.h"
#include "geometry/rectangular_set.h"
#include "geometry/orbit.h"

#include "system/grid_multimap.h"


#include "system/map.h"
#include "system/discrete_time_system.h"

#include "output/logging.h"

namespace Ariadne {
  
namespace Evaluation { 
static int& verbosity = applicator_verbosity; 
}





template<class R> inline
const Evaluation::EvolutionParameters<R>&
Evaluation::ModelChecker<R>::parameters() const
{
  return *this->_parameters;
}








template<class R>
Evaluation::ModelChecker<R>::ModelChecker() 
  : _parameters(new EvolutionParameters<R>())
{
}


template<class R>
Evaluation::ModelChecker<R>::ModelChecker(const EvolutionParameters<R>& parameters)
  : _parameters(new EvolutionParameters<R>(parameters))
{
}



template<class R>
Evaluation::ModelChecker<R>::ModelChecker(const ModelChecker<R>& other) 
  : _parameters(new EvolutionParameters<R>(*other._parameters))
{
}



template<class R>
Evaluation::ModelChecker<R>::~ModelChecker() 
{
  delete this->_parameters;
}


template<class R>
Evaluation::ModelChecker<R>*
Evaluation::ModelChecker<R>::clone() const 
{
  return new ModelChecker<R>(*this);
}

















template<class R>
Geometry:: BoxListSet<R>
Evaluation::ModelChecker<R>::evolve(const System::TransitionSystemInterface<R>& f, 
                                    const Geometry::BoxListSet<R>& initial_set,
                                    const Numeric::Integer& steps) const 
{
  using namespace Numeric;
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  ARIADNE_LOG(2,"ListSet<Box> ModelChecker::reach(TransitionSystemInterface,ListSet<Box>\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  BoxListSet<R> result;
  R mbsr=inf<R>();
  Box<R> r;
  for(typename BoxListSet<R>::const_iterator r_iter=initial_set.begin();
      r_iter!=initial_set.end(); ++r_iter)
  {
    ARIADNE_LOG(6,"  computing iterate for r="<<*r_iter);
    r=*r_iter;
    result.adjoin(f.lower_reach(r));
  }
  return result;
}

template<class R>
Geometry::BoxListSet<R>
Evaluation::ModelChecker<R>::reach(const System::TransitionSystemInterface<R>& f, 
                                   const Geometry::BoxListSet<R>& initial_set,
                                   const Numeric::Integer& steps) const 
{
  using namespace Numeric;
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  ARIADNE_LOG(2,"ListSet<Box> ModelChecker::reach(TransitionSystemInterface,ListSet<Box>\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  BoxListSet<R> result; 
  Box<R> r;
  for(typename BoxListSet<R>::const_iterator r_iter=initial_set.begin();
      r_iter!=initial_set.end(); ++r_iter)
  {
    r=*r_iter;
    result.adjoin(f.lower_reach(r));
  }
  return result;
}


template<class R>
Geometry::BoxListSet<R>
Evaluation::ModelChecker<R>::lower_reach(const System::TransitionSystemInterface<R>& f, 
                                       const Geometry::BoxListSet<R>& initial_set) const 
{
  using namespace Numeric;
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  ARIADNE_LOG(2,"ListSet<Box> ModelChecker::lower_reach(TransitionSystemInterface,ListSet<Box>\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  size_type n=this->_parameters->maximum_number_of_steps();
  R mbsr=this->parameters().maximum_basic_set_radius();
  BoxListSet<R> result;
  for(typename BoxListSet<R>::const_iterator r_iter=initial_set.begin();
      r_iter!=initial_set.end(); ++r_iter)
  {
    ARIADNE_LOG(6,"  computing reach for r="<<*r_iter);
    Box<R> r=*r_iter;
    result.adjoin(f.lower_reach(r));
  }
  return result;
}




template<class R>
Geometry::GridMaskSet<R> 
Evaluation::ModelChecker<R>::chainreach(const System::TransitionSystemInterface<R>& f, 
                                      const Geometry::GridMaskSet<R>& initial_set, 
                                      const Geometry::Box<R>& bounding_set) const
{
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
  ARIADNE_LOG(2,"GridMaskSet ModelChecker::chainreach(TransitionSystemInterface map, GridMaskSet initial_set, GridMaskSet bounding_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\nbounding_set="<<bounding_set<<"\n");
  ARIADNE_CHECK_BOUNDED(initial_set,"GridMaskSet ModelChecker<R>::chainreach(TransitionSystemInterface,GridMaskSet,GridMaskSet)");
  ARIADNE_CHECK_BOUNDED(bounding_set,"GridMaskSet ModelChecker<R>::chainreach(TransitionSystemInterface,GridMaskSet,GridMaskSet)");
  
  const Grid<R>& g=initial_set.grid();
  GridBlock<R> bb=over_approximation(bounding_set,g);
  GridMaskSet<R> result(g,bb.lattice_set());
  GridCellListSet<R> found(g);
  GridCellListSet<R> image(g);
  
  Box<R> r(g.dimension());
  
  uint step=0;
  found=initial_set;
  while(!subset(found,result)) {
    ARIADNE_LOG(3,"Chainreach step "<<step<<": found "<<found.size()<<" cells, ");
    found=difference(found,result);
    ARIADNE_LOG(3,found.size()<<" of which are new.\n");
    ARIADNE_LOG(3,"reached "<<result.size()<<" cells in total.\n");
    result.adjoin(found);
    image.clear(); 
    for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
      image.adjoin(f.upper_evolve(*iter));
    }
    image.unique_sort();
    found.restrict(bb);
    ++step;
  }
  return result;
}



template<class R>
Geometry::GridMaskSet<R>
Evaluation::ModelChecker<R>::viable(const System::TransitionSystemInterface<R>& f, 
                                  const Geometry::GridMaskSet<R>& bounding_set) const
{
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
  ARIADNE_LOG(2,"GridMaskSet ModelChecker::viable(TransitionSystemInterface map, GridMaskSet bounding_set)\n");
  ARIADNE_LOG(3,"bounding_set="<<bounding_set<<"\n");
  ARIADNE_CHECK_BOUNDED(bounding_set,"ModelChecker<R>::viable(TransitionSystemInterface,GridMaskSet)");
  
  const Grid<R>& g=bounding_set.grid();
  Combinatoric::LatticeBlock bd=bounding_set.block();
  GridBlock<R> bb(g,bd);
  GridMaskSet<R> result(g,bd);
  GridCellListSet<R> unsafe(g);
  
  Box<R> r(g.dimension());
  Box<R> fr(g.dimension());
  GridBlock<R> fgb(g);
  Zonotope<R> z(g.dimension());
  Zonotope<R> fz(g.dimension());
  GridCellListSet<R> fgcls(g);
  
  result=bounding_set;
  size_type step=0;

  ARIADNE_LOG(3,"Computing discretization...\n");
  System::GridMultiMap<R> discretization=this->discretize(f,bounding_set,bounding_set.grid());
  ARIADNE_LOG(3,"   Done computing discretization.\n");
  do {
    ARIADNE_LOG(3,"Viability step "<<step+1<< ": testing "<<result.size()<<" cells.\n");
    ++step;
    unsafe.clear();
    for(gms_const_iterator iter=result.begin(); iter!=result.end(); ++iter) {
      fgcls=discretization.image(*iter);
      ARIADNE_LOG(7,"cell="<<*iter<<", image.size()="<<fgcls.size()<<"\n");
      if(!overlap(result,fgcls)) {
        unsafe.adjoin(*iter);
      }
    }
    result.remove(unsafe);
  } while(!unsafe.empty());
  return result;
}


template<class R>
tribool
Evaluation::ModelChecker<R>::verify(const System::TransitionSystemInterface<R>& f, 
                                  const Geometry::GridMaskSet<R>& initial_set, 
                                  const Geometry::GridMaskSet<R>& safe_set) const
{
  ARIADNE_LOG(2,"triboolEvaluation::ModelChecker::verify(TransitionSystemInterface map, GridMaskSet initial_set,GridMaskSet  safe_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n"<<"safe_set="<<safe_set);
  typedef Numeric::Interval<R> I;
  using namespace Geometry;
  typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
  ARIADNE_CHECK_BOUNDED(initial_set,"ModelChecker<R>::verify(...)");
  ARIADNE_CHECK_BOUNDED(safe_set,"ModelChecker<R>::verify(...)");
  
  const Grid<R>& g=initial_set.grid();
  Combinatoric::LatticeBlock bd=safe_set.block();
  Box<R> bb=safe_set.bounding_box();
  GridMaskSet<R> reach(g,bd);
  GridCellListSet<R> found(g);
  GridCellListSet<R> cell_image(g);
  GridCellListSet<R> image(g);
  
  Box<R> r(g.dimension());
  
  found=initial_set;
  while(!subset(found,reach)) {
    found=difference(found,reach);
    reach.adjoin(found);
    image.clear();
    for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
      cell_image=f.upper_reach(*iter);
      if(!subset(cell_image,safe_set)) {
        return false;
      }
      image.adjoin(cell_image);
    }
    found=image;
  }
  return true;
}






template<class R>
System::GridMultiMap<R> 
Evaluation::ModelChecker<R>::discretize(const System::TransitionSystemInterface<R>& f, 
                                      const Geometry::GridMaskSet<R>& domain,
                                      const Geometry::Grid<R>& range_grid) const
{
  ARIADNE_LOG(2,"GridMultiMap*Evaluation::ModelChecker::discretize(TransitionSystemInterface map, GridMaskSet domain, Grid range_grid)\n");
  ARIADNE_LOG(3,"domain="<<domain<<"\n"<<"range_grid="<<range_grid);
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  System::GridMultiMap<R> result(domain.grid(),range_grid);
  for(typename GridMaskSet<R>::const_iterator dom_iter=domain.begin();
      dom_iter!=domain.end(); ++dom_iter)
    {
      const GridCell<R>& gc=*dom_iter;
      GridCellListSet<R> gcls=f.upper_reach(gc);
      result.adjoin_to_image(gc,gcls);
    }
  return result;
}


/*

template<class R>
System::GridMultiMap<R> 
Evaluation::ModelChecker<R>::control_synthesis(const System::DiscreteTimeSystem<R>& f, 
                                             const Geometry::SetInterface<R>& initial_set,
                                             const Geometry::SetInterface<R>& target_set,
                                             const Geometry::GridMaskSet<R>& state_bounding_set,
                                             const Geometry::GridMaskSet<R>& input_bounding_set,
                                             const Geometry::GridMaskSet<R>& noise_bounding_set) const
{
  // TODO: Use on-the-fly discretization
  // TODO: Use adaptive grid refinement
  
  ARIADNE_LOG(2,"GridMultiMap* ModelChecker::control_synthesis(...)\n");
  typedef typename Numeric::traits<R>::interval_type I;
  
  using namespace Combinatoric;
  using namespace Geometry;
  using namespace System;
  
  const Grid<R>& state_grid = state_bounding_set.grid();
  const Grid<R>& input_grid = input_bounding_set.grid();
  dimension_type state_space_dimension = f.state_space_dimension();
  dimension_type input_space_dimension = f.control_space_dimension();
  
  // Discretize function
  std::map< std::pair<LatticeCell,LatticeCell>, LatticeCellListSet > discretization;
  for(typename GridMaskSet<R>::const_iterator state_iter=state_bounding_set.begin();
      state_iter!=state_bounding_set.end(); ++state_iter)
    {
      Point<I> state=static_cast< Box<R> >(*state_iter);
      
      for(typename GridMaskSet<R>::const_iterator input_iter=input_bounding_set.begin();
          input_iter!=input_bounding_set.end(); ++input_iter)
        {
          std::pair<LatticeCell,LatticeCell> control = std::make_pair(state_iter->lattice_set(),input_iter->lattice_set());
          discretization.insert(std::make_pair(control,LatticeCellListSet(state_space_dimension)));
          
          Point<I> input = static_cast< Box<R> >(*input_iter);
          
          for(typename GridMaskSet<R>::const_iterator noise_iter=noise_bounding_set.begin();
              noise_iter!=noise_bounding_set.end(); ++noise_iter)
            {
              Point<I> noise = static_cast< Box<R> >(*noise_iter);
              
              Point<I> image = f.image(state,input,noise);
              
              GridBlock<R> image_set = outer_approximation(image,state_grid);
              discretization.find(control)->second.adjoin(image_set.lattice_set());
            }
        }
    }
  
  // Discretize target set
  GridMaskSet<R> target_approximation(state_bounding_set.grid(),state_bounding_set.block());
  target_approximation.adjoin_inner_approximation(target_set);
  Combinatoric::LatticeMaskSet target_lattice_set = target_approximation.lattice_set();
  
  // Discretize initial set
  GridMaskSet<R> initial_approximation(state_bounding_set.grid(),state_bounding_set.block());
  initial_approximation.adjoin_inner_approximation(initial_set);
  Combinatoric::LatticeMaskSet initial_lattice_set = initial_approximation.lattice_set();
  
  Combinatoric::LatticeMaskSet bounding_lattice_set = state_bounding_set.lattice_set();
  Combinatoric::LatticeMaskSet input_lattice_set = input_bounding_set.lattice_set();
  
  // Solve control problem
  Combinatoric::LatticeMultiMap lattice_control(state_space_dimension,input_space_dimension);
  Combinatoric::LatticeMaskSet controllable_lattice_set(state_bounding_set.block());
  Combinatoric::LatticeMaskSet new_controllable_lattice_set(state_bounding_set.block());;
  
  new_controllable_lattice_set.adjoin(target_lattice_set);
  controllable_lattice_set.adjoin(new_controllable_lattice_set);
  while(!new_controllable_lattice_set.empty() && !subset(initial_lattice_set,controllable_lattice_set)) {
    new_controllable_lattice_set.clear();
    for(LatticeMaskSet::const_iterator state_iter = bounding_lattice_set.begin();
        state_iter!=bounding_lattice_set.end(); ++state_iter)
      {
        if(!subset(*state_iter,controllable_lattice_set)) {
          for(LatticeMaskSet::const_iterator input_iter = input_lattice_set.begin();
              input_iter!=input_lattice_set.end(); ++input_iter)
            
            {
              if(Combinatoric::subset(discretization.find(std::make_pair(*state_iter,*input_iter))->second,controllable_lattice_set)) {
                new_controllable_lattice_set.adjoin(*state_iter);
                lattice_control.adjoin_to_image(*state_iter,*input_iter);
              }
            }
        }
      }
    controllable_lattice_set.adjoin(new_controllable_lattice_set);
  }
  
  System::GridMultiMap<R> result(state_grid,input_grid,lattice_control);
  
  return result;      
}

*/






}
