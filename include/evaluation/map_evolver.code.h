/***************************************************************************
 *            map_evolver.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
#include "map_evolver.h"
#include "map_orbiter.h"

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

#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "geometry/grid.h"
#include "geometry/grid_set.h"
#include "geometry/partition_tree_set.h"
#include "geometry/grid_approximation.h"
#include "geometry/rectangular_set.h"
#include "geometry/orbit.h"

#include "system/grid_multimap.h"


#include "system/map_interface.h"
#include "system/discrete_map.h"
#include "system/discrete_time_system.h"

#include "evaluation/applicator.h"
#include "evaluation/model_checker.h"

#include "output/logging.h"

namespace Ariadne {
  
namespace Evaluation { 
static int& verbosity = applicator_verbosity; 
static const double DEFAULT_MAXIMUM_BASIC_SET_RADIUS=0.25;
static const double DEFAULT_GRID_LENGTH=0.125;
}



template<class R> inline
Evaluation::EvolutionParameters<R>*
Evaluation::MapEvolver<R>::default_parameters() 
{
  return new EvolutionParameters<R>();
}

template<class R> inline
Evaluation::MapOrbiterInterface<R>*
Evaluation::MapEvolver<R>::default_orbiter() 
{
  typedef Geometry::Zonotope<I,R> BS;
  const EvolutionParameters<R>& parameters=*this->_parameters;
  Applicator<R> applicator;
  return new MapOrbiter<BS>(parameters,applicator);
}

template<class R> inline
Evaluation::ModelChecker<R>
Evaluation::MapEvolver<R>::model_checker() const
{
  return ModelChecker<R>(*this->_parameters);
}

template<class R> inline
System::DiscreteMap<R>
Evaluation::MapEvolver<R>::discrete_map(const System::MapInterface<R>& f) const
{
  return System::DiscreteMap<R>(f,*this->_orbiter);
}



template<class R>
Evaluation::MapEvolver<R>::MapEvolver() 
  : _parameters(default_parameters()),
    _orbiter(default_orbiter())
{
  _parameters->set_maximum_basic_set_radius(DEFAULT_MAXIMUM_BASIC_SET_RADIUS);
  _parameters->set_grid_length(DEFAULT_GRID_LENGTH);
}


template<class R>
Evaluation::MapEvolver<R>::MapEvolver(const EvolutionParameters<R>& parameters)
  : _parameters(new EvolutionParameters<R>(parameters)),
    _orbiter(default_orbiter())
{
}



template<class R>
Evaluation::MapEvolver<R>::MapEvolver(const MapEvolver<R>& other) 
  : _parameters(new EvolutionParameters<R>(other.parameters())),
    _orbiter(other._orbiter->clone())
{
}



template<class R>
Evaluation::MapEvolver<R>::~MapEvolver() 
{
  delete this->_parameters;
  delete this->_orbiter;
}


template<class R>
Evaluation::MapEvolver<R>*
Evaluation::MapEvolver<R>::clone() const 
{
  return new MapEvolver<R>(*this);
}




template<class R>
const Evaluation::EvolutionParameters<R>&
Evaluation::MapEvolver<R>::parameters() const
{
  return *this->_parameters;
}


template<class R>
Evaluation::EvolutionParameters<R>&
Evaluation::MapEvolver<R>::parameters() 
{
  return *this->_parameters;
}









template<class R>
Geometry::Rectangle<R>
Evaluation::MapEvolver<R>::apply(const System::MapInterface<R>& f, const Geometry::Rectangle<R>& r) const
{
  ARIADNE_LOG(4,"BS MapEvolver::apply(MapInterface,Rectangle)\n");
  return this->_orbiter->apply(f,r);
}

template<class R>
Geometry::GridCellListSet<R>
Evaluation::MapEvolver<R>::apply(const System::MapInterface<R>& f, const Geometry::GridCell<R>& gc) const
{
  ARIADNE_LOG(4,"GridCellListSet MapEvolver::apply(MapInterface,GridCell)\n");
  return this->_orbiter->apply(f,gc,gc.grid());
}

template<class R>
Geometry::GridCellListSet<R>
Evaluation::MapEvolver<R>::apply(const System::MapInterface<R>& f, const Geometry::GridCell<R>& gc, const Geometry::Grid<R>& g) const
{
  ARIADNE_LOG(4,"GridCellListSet MapEvolver::apply(MapInterface,GridCell)\n");
  return this->_orbiter->apply(f,gc,g);
}




template<class R>
Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::Rectangle<R> >
Evaluation::MapEvolver<R>::orbit(const System::MapInterface<R>& f, const Geometry::Rectangle<R>& r, const Numeric::Integer& n) const
{
  ARIADNE_LOG(4,"DiscreteTimeOrbit<Integer,Rectangle> MapEvolver::orbit(MapInterface,RectangleRectangle)\n");
  return this->_orbiter->orbit(f,r,n,Numeric::inf<R>());
}

template<class R>
Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::Rectangle<R> >
Evaluation::MapEvolver<R>::orbit(const System::MapInterface<R>& f, const Geometry::Rectangle<R>& r, const Numeric::Integer& n, const R& mbsr) const
{
  ARIADNE_LOG(4,"DiscreteTimeOrbit<Integer,Rectangle> MapEvolver::orbit(MapInterface,Rectangle,Integer,Float)\n");
  return this->_orbiter->orbit(f,r,n,mbsr);
}

template<class R>
Geometry::DiscreteTimeOrbit< Numeric::Integer, Geometry::GridCellListSet<R> >
Evaluation::MapEvolver<R>::orbit(const System::MapInterface<R>& f, const Geometry::GridCell<R>& gc, const Numeric::Integer& n) const
{
  ARIADNE_LOG(4,"DiscreteTimeOrbit<Integer,GridCellListSet> MapEvolver::orbit(MapInterface,GridCell,Integer)\n");
  return this->_orbiter->orbit(f,gc,n);
}








template<class R>
Geometry::ListSet< Geometry::Rectangle<R> >
Evaluation::MapEvolver<R>::image(const System::MapInterface<R>& f, const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set) const 
{
  ARIADNE_LOG(2,"GridMaskSet MapEvolver::image(MapInterface map, ListSet< Rectangle<Float> > initial_set)\n");
  return this->model_checker().image(this->discrete_map(f),initial_set);
}


template<class R>
Geometry::GridCellListSet<R> 
Evaluation::MapEvolver<R>::image(const System::MapInterface<R>& f, 
                                 const Geometry::GridCellListSet<R>& initial_set, 
                                 const Geometry::Grid<R>& image_grid) const 
{
  ARIADNE_LOG(2,"GridMaskSet MapEvolver::image(MapInterface map, GridCellListSet initial_set, Grid image_grid)\n");
  return this->model_checker().image(this->discrete_map(f),initial_set,image_grid);
}




template<class R>
Geometry::GridMaskSet<R> 
Evaluation::MapEvolver<R>::image(const System::MapInterface<R>& f, 
                                 const Geometry::GridMaskSet<R>& initial_set, 
                                 const Geometry::GridMaskSet<R>& bounding_set) const 
{
  ARIADNE_LOG(2,"GridMaskSet MapEvolver::preimage(MapInterface,PartitionTreeSet,Rectangle)\n");
  return this->model_checker().image(this->discrete_map(f),initial_set,bounding_set);
}





template<class R>
Geometry::GridMaskSet<R> 
Evaluation::MapEvolver<R>::preimage(const System::MapInterface<R>& f, 
                                    const Geometry::GridMaskSet<R>& set, 
                                    const Geometry::GridMaskSet<R>& bounding_set) const 
{
  ARIADNE_LOG(2,"GridMaskSet MapEvolver::preimage(MapInterface,GridMaskSet,GridMaskSet)\n");
  return this->model_checker().preimage(this->discrete_map(f),set,bounding_set);
}

template<class R>
Geometry::PartitionTreeSet<R> 
Evaluation::MapEvolver<R>::preimage(const System::MapInterface<R>& f, 
                                    const Geometry::PartitionTreeSet<R>& set, 
                                    const Geometry::Rectangle<R>& bound) const 
{
  ARIADNE_LOG(2,"GridMaskSet MapEvolver::preimage(MapInterface,PartitionTreeSet,Rectangle)\n");
  return this->model_checker().preimage(this->discrete_map(f),set,bound);
}



template<class R>
Geometry::ListSet< Geometry::Rectangle<R> >
Evaluation::MapEvolver<R>::iterate(const System::MapInterface<R>& f, 
                                   const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set,
                                   const Numeric::Integer& steps) const 
{
  ARIADNE_LOG(2,"ListSet<Rectangle> MapEvolver::initial(MapInterface,ListSet<Rectangle>,Integer)\n");
  return this->model_checker().iterate(this->discrete_map(f),initial_set,steps);
}

template<class R>
Geometry::ListSet< Geometry::Rectangle<R> >
Evaluation::MapEvolver<R>::reach(const System::MapInterface<R>& f, 
                                 const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set,
                                 const Numeric::Integer& steps) const 
{
  ARIADNE_LOG(2,"ListSet<Rectangle> MapEvolver::reach(MapInterface,ListSet<Rectangle>,Integer)\n");
  return this->model_checker().reach(this->discrete_map(f),initial_set,steps);
}


template<class R>
Geometry::ListSet< Geometry::Rectangle<R> >
Evaluation::MapEvolver<R>::lower_reach(const System::MapInterface<R>& f, 
                                       const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set) const 
{
  ARIADNE_LOG(2,"ListSet<Rectangle> MapEvolver::(MapInterface map, ListSet<Rectangle> initial_set\n");
  return this->model_checker().lower_reach(this->discrete_map(f),initial_set);
}




template<class R>
Geometry::GridMaskSet<R> 
Evaluation::MapEvolver<R>::chainreach(const System::MapInterface<R>& f, 
                                      const Geometry::GridMaskSet<R>& initial_set, 
                                      const Geometry::GridMaskSet<R>& bounding_set) const
{
  ARIADNE_LOG(2,"GridMaskSet MapEvolver::chainreach(MapInterface map, GridMaskSet initial_set, GridMaskSet bounding_set)\n");
  return this->model_checker().chainreach(this->discrete_map(f),initial_set,bounding_set);
}



template<class R>
Geometry::GridMaskSet<R>
Evaluation::MapEvolver<R>::viable(const System::MapInterface<R>& f, 
                                  const Geometry::GridMaskSet<R>& bounding_set) const
{
  ARIADNE_LOG(2,"GridMaskSet MapEvolver::viable(MapInterface map, GridMaskSet bounding_set)\n");
  return this->model_checker().viable(this->discrete_map(f),bounding_set);
}


template<class R>
tribool
Evaluation::MapEvolver<R>::verify(const System::MapInterface<R>& f, 
                                  const Geometry::GridMaskSet<R>& initial_set, 
                                  const Geometry::GridMaskSet<R>& safe_set) const
{
  ARIADNE_LOG(2,"triboolEvaluation::MapEvolver::verify(MapInterface map, GridMaskSet initial_set, GridMaskSet safe_set)\n");
  return this->model_checker().verify(this->discrete_map(f),initial_set,safe_set);
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n"<<"safe_set="<<safe_set);
  typedef Numeric::Interval<R> I;
  using namespace Geometry;
  typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
  ARIADNE_CHECK_BOUNDED(initial_set,"MapEvolver<R>::verify(...)");
  ARIADNE_CHECK_BOUNDED(safe_set,"MapEvolver<R>::verify(...)");
  
  const Grid<R>& g=initial_set.grid();
  Combinatoric::LatticeBlock bd=safe_set.block();
  Rectangle<R> bb=safe_set.bounding_box();
  GridMaskSet<R> reach(g,bd);
  GridCellListSet<R> found(g);
  GridCellListSet<R> cell_image(g);
  GridCellListSet<R> image(g);
  
  Rectangle<R> r(g.dimension());
  
  found=initial_set;
  while(!subset(found,reach)) {
    found=difference(found,reach);
    reach.adjoin(found);
    image.clear();
    for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
      cell_image=this->apply(f,*iter);
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
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::image(const System::MapInterface<R>& f, 
                                 const Geometry::SetInterface<R>& set) const
{
  // FIXME: Only computes an over-approximation to the image.
  using namespace Geometry;
  using namespace System;
  typedef Numeric::Interval<R> I;
  ARIADNE_LOG(2,"SetInterface* MapEvolver::image(MapInterface,SetInterface)\n");
  ARIADNE_LOG(3,"set="<<set<<"\n");
  Grid<R> argument_grid(f.argument_dimension(),this->parameters().grid_length());
  Grid<R> result_grid(f.result_dimension(),this->parameters().grid_length());
  GridCellListSet<R> grid_initial_set=Geometry::outer_approximation(set,argument_grid);
  DiscreteMap<R> map(f,*this->_orbiter);
  ModelChecker<R> model_checker(this->parameters());
  return new GridCellListSet<R>(model_checker.image(map,grid_initial_set,result_grid));
}



template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::preimage(const System::MapInterface<R>& map, 
                                    const Geometry::SetInterface<R>& set,
                                    const Geometry::SetInterface<R>& bound) const
{
  typedef Numeric::Interval<R> I;
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* MapEvolver::preimage(MapInterface map, SetInterface set)\n");
  ARIADNE_LOG(3,"set="<<set<<"\n");
  Rectangle<R> preimage_bounding_box=bound.bounding_box();
  Grid<R> preimage_grid(map.argument_dimension(),this->parameters().grid_length());
  Grid<R> image_grid(map.result_dimension(),this->parameters().grid_length());
  GridMaskSet<R> grid_image_set(image_grid,set.bounding_box());
  grid_image_set.adjoin_inner_approximation(set);
  GridMaskSet<R> grid_preimage_bounding_set(preimage_grid,preimage_bounding_box);
  grid_preimage_bounding_set.adjoin_outer_approximation(bound);
  return new GridMaskSet<R>(this->preimage(map,grid_image_set,grid_preimage_bounding_set));
}



template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::iterate(const System::MapInterface<R>& f, 
                                   const Geometry::SetInterface<R>& initial_set,
                                   const Numeric::Integer& steps) const
{
  typedef Numeric::Interval<R> I;
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* MapEvolver::iterate(MapInterface,SetInterface,Integer)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  ARIADNE_LOG(3,"steps="<<steps<<"\n");
  const ListSet< Rectangle<R> >* list_initial_set_ptr=dynamic_cast<const ListSet< Rectangle<R> >*>(&initial_set);
  ARIADNE_LOG(4,"list_initial_set_ptr="<<list_initial_set_ptr<<"\n");
  if(list_initial_set_ptr) {
    return new ListSet< Rectangle<R> >(this->iterate(f,*list_initial_set_ptr,steps));
  }
  
  ARIADNE_CHECK_BOUNDED(initial_set,"SetInterface* MapEvolver::iterate(MapInterface,SetInterface,Integer)");
  
  ListSet< Rectangle<R> > list_initial_set;
  if(dynamic_cast<const RectangularSet<R>*>(&initial_set)) {
    ARIADNE_LOG(4,"Cast to RectangularSet<R>\n");
  }
  
  Grid<R> grid(initial_set.dimension(),this->parameters().grid_length());
  Rectangle<R> bounding_box=initial_set.bounding_box();
  GridMaskSet<R> gms(grid,bounding_box);
  gms.adjoin_inner_approximation(initial_set);
  ARIADNE_LOG(4,"gms="<<gms);
  list_initial_set=ListSet< Rectangle<R> >(gms);
  
  return new ListSet< Rectangle<R> >(this->iterate(f,list_initial_set,steps));      
}

template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::reach(const System::MapInterface<R>& f, 
                                       const Geometry::SetInterface<R>& initial_set,
                                       const Numeric::Integer& steps) const
{
  typedef Numeric::Interval<R> I;
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* MapEvolver::reach(MapInterface,SetInterface,Integer)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  ARIADNE_LOG(3,"steps="<<steps<<"\n");
  const ListSet< Rectangle<R> >* list_initial_set_ptr=dynamic_cast<const ListSet< Rectangle<R> >*>(&initial_set);
  ARIADNE_LOG(4,"list_initial_set_ptr="<<list_initial_set_ptr<<"\n");
  if(list_initial_set_ptr) {
    return new ListSet< Rectangle<R> >(this->reach(f,*list_initial_set_ptr,steps));
  }
  
  ARIADNE_CHECK_BOUNDED(initial_set,"SetInterface* MapEvolver::reach(MapInterface,SetInterface,Integer)");
  
  ListSet< Rectangle<R> > list_initial_set;
  if(dynamic_cast<const RectangularSet<R>*>(&initial_set)) {
    ARIADNE_LOG(4,"Cast to RectangularSet<R>\n");
  }
  
  Grid<R> grid(initial_set.dimension(),this->parameters().grid_length());
  Rectangle<R> bounding_box=initial_set.bounding_box();
  GridMaskSet<R> gms(grid,bounding_box);
  gms.adjoin_inner_approximation(initial_set);
  ARIADNE_LOG(4,"gms="<<gms);
  list_initial_set=ListSet< Rectangle<R> >(gms);
  
  return new ListSet< Rectangle<R> >(this->reach(f,list_initial_set,steps));      
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::lower_reach(const System::MapInterface<R>& f, 
                                       const Geometry::SetInterface<R>& initial_set) const
{
  typedef Numeric::Interval<R> I;
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* MapEvolver::lower_reach(MapInterface,SetInterface)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  
  const ListSet< Rectangle<R> >* rectangle_list_initial_set_ptr=dynamic_cast<const ListSet< Rectangle<R> >*>(&initial_set);
  ARIADNE_LOG(4,"rectangle_list_initial_set_ptr="<<rectangle_list_initial_set_ptr<<"\n");
  if(rectangle_list_initial_set_ptr) {
    return new ListSet< Rectangle<R> >(this->lower_reach(f,*rectangle_list_initial_set_ptr));
  }
  
  ARIADNE_CHECK_BOUNDED(initial_set,"SetInterface* MapEvolver::lower_reach(MapInterface,SetInterface)");
  
  ListSet< Rectangle<R> > list_initial_set;
  if(dynamic_cast<const RectangularSet<R>*>(&initial_set)) {
    ARIADNE_LOG(4,"Cast to RectangularSet<R>\n");
  }
  
  Grid<R> grid(initial_set.dimension(),this->parameters().grid_length());
  ListSet< Rectangle<R> > rls=Geometry::point_approximation(initial_set,grid);
  ARIADNE_LOG(4,"rls="<<rls);
  
  return new ListSet< Rectangle<R> >(this->lower_reach(f,rls));      
}



template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::chainreach(const System::MapInterface<R>& map, 
                                      const Geometry::SetInterface<R>& initial_set) const
{
  Geometry::RectangularSet<R> bounding_set(this->_parameters->bounding_box(map.argument_dimension()));
  return this->chainreach(map,initial_set,bounding_set);
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::chainreach(const System::MapInterface<R>& map, 
                                      const Geometry::SetInterface<R>& initial_set, 
                                      const Geometry::SetInterface<R>& bounding_set) const
{
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* MapEvolver::chainreach(MapInterface,SetInterface,SetInterface)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n"<<"bounding_set="<<bounding_set<<"\n");
  
  ARIADNE_CHECK_BOUNDED(bounding_set,"SetInterface* MapEvolver::chainreach(MapInterface map, SetInterface initial_set, SetInterface bounding_set)");
  Rectangle<R> bounding_box=bounding_set.bounding_box();
  Grid<R> grid(bounding_set.dimension(),this->parameters().grid_length());
  FiniteGrid<R> finite_grid(grid,bounding_box);
  GridMaskSet<R> grid_bounding_set(finite_grid);
  grid_bounding_set.adjoin_outer_approximation(bounding_set);
  GridMaskSet<R> grid_initial_set(finite_grid);
  grid_initial_set.adjoin_outer_approximation(initial_set);
  
  return new GridMaskSet<R>(this->chainreach(map,grid_initial_set,grid_bounding_set));
}




template<class R>
Geometry::SetInterface<R>*
Evaluation::MapEvolver<R>::viable(const System::MapInterface<R>& map, 
                                  const Geometry::SetInterface<R>& bounding_set) const
{
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* MapEvolver::viable(MapInterface,SetInterface)\n");
  ARIADNE_LOG(3,"bounding_set="<<bounding_set<<"\n");
  ARIADNE_CHECK_BOUNDED(bounding_set,"SetInterface* MapEvolver::viable(MapInterface map, SetInterface bounding_set)");
  Rectangle<R> bounding_box=bounding_set.bounding_box();
  Grid<R> grid(bounding_set.dimension(),this->parameters().grid_length());
  GridMaskSet<R> grid_bounding_set(grid,bounding_box);
  grid_bounding_set.adjoin_outer_approximation(bounding_set);
  return new GridMaskSet<R>(this->viable(map,grid_bounding_set));
}



template<class R>
tribool
Evaluation::MapEvolver<R>::verify(const System::MapInterface<R>& f, 
                                  const Geometry::SetInterface<R>& initial_set, 
                                  const Geometry::SetInterface<R>& safe_set) const
{
  using namespace Geometry;
  ARIADNE_LOG(2,"triboolEvaluation::MapEvolver::verify(MapInterface,SetInterface,SetInterface)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n"<<"safe_set="<<safe_set<<"\n");
  ARIADNE_CHECK_BOUNDED(safe_set,"SetInterface* MapEvolver::verify(MapInterface map, SetInterface initial_set, SetInterface safe_set)");
  Rectangle<R> bounding_box=safe_set.bounding_box();
  Grid<R> grid(safe_set.dimension(),this->parameters().grid_length());
  FiniteGrid<R> finite_grid(grid,bounding_box);
  GridMaskSet<R> grid_inner_safe_set(finite_grid);
  grid_inner_safe_set.adjoin_inner_approximation(safe_set);
  GridMaskSet<R> grid_initial_set(finite_grid);
  grid_initial_set.adjoin_outer_approximation(initial_set);
  
  return this->verify(f,grid_initial_set,grid_inner_safe_set);
}




template<class R>
System::GridMultiMap<R> 
Evaluation::MapEvolver<R>::discretize(const System::MapInterface<R>& f, 
                                      const Geometry::GridMaskSet<R>& domain,
                                      const Geometry::Grid<R>& range_grid) const
{
  ARIADNE_LOG(2,"GridMultiMap*Evaluation::MapEvolver::discretize(MapInterface map, GridMaskSet domain, Grid range_grid)\n");
  ARIADNE_LOG(3,"domain="<<domain<<"\n"<<"range_grid="<<range_grid);
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  System::GridMultiMap<R> result(domain.grid(),range_grid);
  for(typename GridMaskSet<R>::const_iterator dom_iter=domain.begin();
      dom_iter!=domain.end(); ++dom_iter)
    {
      const GridCell<R>& gc=*dom_iter;
      GridCellListSet<R> gcls=this->apply(f,gc,range_grid);
      result.adjoin_to_image(gc,gcls);
    }
  return result;
}



template<class R>
System::GridMultiMap<R> 
Evaluation::MapEvolver<R>::control_synthesis(const System::DiscreteTimeSystem<R>& f, 
                                             const Geometry::SetInterface<R>& initial_set,
                                             const Geometry::SetInterface<R>& target_set,
                                             const Geometry::GridMaskSet<R>& state_bounding_set,
                                             const Geometry::GridMaskSet<R>& input_bounding_set,
                                             const Geometry::GridMaskSet<R>& noise_bounding_set) const
{
  // TODO: Use on-the-fly discretization
  // TODO: Use adaptive grid refinement
  
  ARIADNE_LOG(2,"GridMultiMap* MapEvolver::control_synthesis(...)\n");
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
      Point<I> state=static_cast< Rectangle<R> >(*state_iter);
      
      for(typename GridMaskSet<R>::const_iterator input_iter=input_bounding_set.begin();
          input_iter!=input_bounding_set.end(); ++input_iter)
        {
          std::pair<LatticeCell,LatticeCell> control = std::make_pair(state_iter->lattice_set(),input_iter->lattice_set());
          discretization.insert(std::make_pair(control,LatticeCellListSet(state_space_dimension)));
          
          Point<I> input = static_cast< Rectangle<R> >(*input_iter);
          
          for(typename GridMaskSet<R>::const_iterator noise_iter=noise_bounding_set.begin();
              noise_iter!=noise_bounding_set.end(); ++noise_iter)
            {
              Point<I> noise = static_cast< Rectangle<R> >(*noise_iter);
              
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

}
