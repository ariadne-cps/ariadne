/***************************************************************************
 *            vector_field_evolver.code.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
#include <cassert>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>
#include <cstring>
#include <cassert>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../base/array.h"

#include "../numeric/rational.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid.h"
#include "../geometry/grid_set.h"

#include "../system/vector_field.h"
#include "../system/affine_vector_field.h"

#include "../evaluation/evolution_parameters.h"
#include "../evaluation/bounder.h"
#include "../evaluation/integrator_interface.h"

#include "../output/logging.h"

#include "vector_field_evolver.h"

namespace {
    
using namespace Ariadne;


template<class BS, class BST>
struct pair_first_less {
  bool operator()(const std::pair<BS,BST>& ts1, const std::pair<BS,BST>& ts2) {
    return ts1.first < ts2.first; 
  }
};

/*!\ brief A class representing pre-computed bounds for an integration step. */
template<class BS>
class IntegrationStepBound {
  typedef typename BS::real_type R;
 public:
  /*!\ brief Constructor. */
  IntegrationStepBound(const Geometry::Rectangle<R>& bound, const BS& integration_time) 
    : _bound(bound), _integration_time(integration_time) { }
  /*!\ brief Constructor. */
  IntegrationStepBound(const Geometry::Rectangle<R>& bound, const Numeric::Interval<BS>& integration_time) 
    : _bound(bound), _integration_time(integration_time.upper()) { }
  /*! The spacial bound for the integrations step. */
  const Geometry::Rectangle<R>& bound() const { return _bound; }
  /*! The step size in time of the integrations step. */
  const BS& integration_time() const { return _integration_time; }
 private:
  Geometry::Rectangle<R> _bound;
  BS _integration_time;
};

}


namespace Ariadne {
    
namespace Evaluation { static int& verbosity = integrator_verbosity; }
    

using namespace Geometry;
using namespace System;

template<class R>
Evaluation::VectorFieldEvolver<R>::~VectorFieldEvolver()
{
  delete _parameters;
  delete _bounder;
  delete _integrator;
}

template<class R>
Evaluation::VectorFieldEvolver<R>::VectorFieldEvolver(const VectorFieldEvolver<R>& i)
  : _parameters(i._parameters->clone()),
    _bounder(i._bounder->clone()),
    _integrator(i._integrator->clone())
{
}


template<class R>
Evaluation::VectorFieldEvolver<R>*
  Evaluation::VectorFieldEvolver<R>::clone() const
{
  return new VectorFieldEvolver<R>(*this);
}

template<class R>
Evaluation::VectorFieldEvolver<R>::VectorFieldEvolver(const EvolutionParameters<R>& parameters, const IntegratorInterface<R>& plugin)
  : _parameters(new EvolutionParameters<R>(parameters)),
    _bounder(new Bounder<R>()),
    _integrator(plugin.clone())
{
}



template<class R>
Evaluation::EvolutionParameters<R>&
Evaluation::VectorFieldEvolver<R>::parameters() 
{
  return *this->_parameters;
}

template<class R>
const Evaluation::EvolutionParameters<R>&
Evaluation::VectorFieldEvolver<R>::parameters() const
{
  return *this->_parameters;
}

template<class R>
time_type 
Evaluation::VectorFieldEvolver<R>::minimum_step_size() const
{
  return this->_parameters->minimum_step_size();
}

template<class R>
time_type 
Evaluation::VectorFieldEvolver<R>::maximum_step_size() const
{
  return this->_parameters->maximum_step_size();
}

template<class R>
R
Evaluation::VectorFieldEvolver<R>::minimum_basic_set_radius() const
{
  return this->_parameters->minimum_basic_set_radius();
}

template<class R>
R
Evaluation::VectorFieldEvolver<R>::maximum_basic_set_radius() const
{
  return this->_parameters->maximum_basic_set_radius();
}

template<class R>
R
Evaluation::VectorFieldEvolver<R>::grid_length() const
{
  return this->_parameters->grid_length();
}

template<class R>
time_type 
Evaluation::VectorFieldEvolver<R>::lock_to_grid_time() const
{
  return this->_parameters->lock_to_grid_time();
}




template<class R>
typename Evaluation::VectorFieldEvolver<R>::list_set_type
Evaluation::VectorFieldEvolver<R>::subdivide(const basic_set_type& set) const
{
  ARIADNE_LOG(5,"VectorFieldEvolver::subdivide(BasicSet)\n");
  return set.subdivide();
}


template<class R>
typename Evaluation::VectorFieldEvolver<R>::basic_set_type
Evaluation::VectorFieldEvolver<R>::integration_step(const vector_field_type& vector_field, 
                                             const basic_set_type& initial_set, 
                                             time_type& step_size) const
{
  ARIADNE_LOG(2,"BasicSet VectorFieldEvolver::integration_step(VectorFieldInterface,BasicSet,Time)\n");
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"VectorFieldEvolver::integration_step(VectorFieldInterface,BasicSet,Time)");

  bounding_set_type bounding_set=this->_bounder->estimate_flow_bounds(vector_field,initial_set.bounding_box(),step_size);
  bounding_set=this->_bounder->refine_flow_bounds(vector_field,initial_set.bounding_box(),bounding_set,step_size);
  bounding_set=this->_bounder->refine_flow_bounds(vector_field,initial_set.bounding_box(),bounding_set,step_size);
  return this->_integrator->integration_step(vector_field,initial_set,step_size,bounding_set);
}


template<class R>
typename Evaluation::VectorFieldEvolver<R>::basic_set_type
Evaluation::VectorFieldEvolver<R>::reachability_step(const vector_field_type& vector_field, 
                                              const basic_set_type& initial_set, 
                                              time_type& step_size) const
{
  ARIADNE_LOG(2,"BasicSet VectorFieldEvolver::reachability_step(VectorFieldInterface,BasicSet,Time)\n");
  ARIADNE_CHECK_EQUAL_DIMENSIONS(vector_field,initial_set,"VectorFieldEvolver::reachability_step(VectorFieldInterface,BasicSet,Time)");
  bounding_set_type bounding_set=this->_bounder->estimate_flow_bounds(vector_field,initial_set.bounding_box(),step_size);
  bounding_set=this->_bounder->refine_flow_bounds(vector_field,initial_set.bounding_box(),bounding_set,step_size);
  bounding_set=this->_bounder->refine_flow_bounds(vector_field,initial_set.bounding_box(),bounding_set,step_size);
  return this->_integrator->reachability_step(vector_field,initial_set,step_size,bounding_set);
}




template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::integrate(const System::VectorFieldInterface<R>& vector_field, 
                                      const Geometry::SetInterface<R>& initial_set, 
                                      const Numeric::Rational& time) const
{
  ARIADNE_LOG(2,"SetInterface* VectorFieldEvolver::integrate(VectorFieldInterface,SetInterface>,Time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  const VectorFieldInterface<R>& cast_vector_field=vector_field;
  Rectangle<R> bb=initial_set.bounding_box();
  Grid<R> grid(vector_field.dimension(),this->grid_length());
  ListSet< Rectangle<R> > rectangle_list_initial_set=lower_approximation(initial_set,grid);
  ARIADNE_LOG(3,"rectangle_list_initial_set="<<rectangle_list_initial_set);
  ListSet<BS> list_initial_set(rectangle_list_initial_set);
  ARIADNE_LOG(3,"list_initial_set="<<list_initial_set);
  ListSet<BS> list_final_set=this->lower_integrate(cast_vector_field,list_initial_set,time);
  ARIADNE_LOG(3,"list_final_set="<<list_final_set);
  return new ListSet<BS>(list_final_set);
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::integrate(const System::VectorFieldInterface<R>& vector_field,
                                     const Geometry::SetInterface<R>& initial_set,
                                     const Geometry::SetInterface<R>& bounding_set,
                                     const time_type& time) const
{
  // FIXME: Only computes over-approximation;
  using namespace Geometry;
  const System::VectorFieldInterface<R>& vf=vector_field;
  
  Rectangle<R> bb;
  try {
    bb=bounding_set.bounding_box();
  }
  catch(UnboundedSet&) {
    throw UnboundedSet("integrate(Map,Set,Set,Time): bounding_set unbounded");
  }
  Grid<R> g(bounding_set.dimension(),this->grid_length());
  FiniteGrid<R> fg(g,bb);
  GridMaskSet<R> gbs(fg);
  gbs.adjoin_outer_approximation(bounding_set);
  GridMaskSet<R> gis(fg);
  gis.adjoin_outer_approximation(initial_set);
  
  GridMaskSet<R> gs=this->reach(vf,gis,gbs,time);
  return new GridMaskSet<R>(gs);
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::reach(const System::VectorFieldInterface<R>& vector_field,
                                  const Geometry::SetInterface<R>& initial_set,
                                  const time_type& time) const
{
  using namespace Geometry;
  typedef Numeric::Interval<BS> I;
  ARIADNE_LOG(2,"SetInterface* VectorFieldEvolver::reach(VectorFile,SetInterface,Time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  const VectorFieldInterface<R>& cast_vector_field=vector_field;
  Rectangle<R> bb=initial_set.bounding_box();
  Grid<R> grid(vector_field.dimension(),this->grid_length());
  ListSet< Rectangle<R> > rectangle_list_initial_set=lower_approximation(initial_set,grid);
  ListSet<BS> list_initial_set(rectangle_list_initial_set);
  ListSet<BS> list_reach_set=this->lower_reach(cast_vector_field,list_initial_set,time);
  return new ListSet<BS>(list_reach_set);
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::reach(const System::VectorFieldInterface<R>& vector_field,
                                 const Geometry::SetInterface<R>& initial_set,
                                 const Geometry::SetInterface<R>& bounding_set,
                                 const time_type& time) const
{
  using namespace Geometry;
  const System::VectorFieldInterface<R>& vf=vector_field;
  
  Rectangle<R> bb;
  try {
    bb=bounding_set.bounding_box();
  }
  catch(UnboundedSet&) {
    throw UnboundedSet("reach(Map,Set,Set,Time): bounding_set unbounded");
  }
  Grid<R> g(bounding_set.dimension(),this->grid_length());
  FiniteGrid<R> fg(g,bb);
  GridMaskSet<R> gbs(fg);
  gbs.adjoin_outer_approximation(bounding_set);
  GridMaskSet<R> gis(fg);
  gis.adjoin_outer_approximation(initial_set);
  
  GridMaskSet<R> grs=this->reach(vf,gis,gbs,time);
  return new GridMaskSet<R>(grs);
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::reach(const System::VectorFieldInterface<R>& vector_field,
                                           const Geometry::SetInterface<R>& initial_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::chainreach(const System::VectorFieldInterface<R>& vector_field,
                                      const Geometry::SetInterface<R>& initial_set,
                                      const Geometry::SetInterface<R>& bounding_set) const
{
  using namespace Geometry;
  const System::VectorFieldInterface<R>& vf=vector_field;
  
  Rectangle<R> bb;
  try {
    bb=bounding_set.bounding_box();
  }
  catch(UnboundedSet&) {
    throw UnboundedSet("chainreach(Map,Set,Set): bounding_set unbounded");
  }
  Grid<R> g(bounding_set.dimension(),this->grid_length());
  FiniteGrid<R> fg(g,bb);
  GridMaskSet<R> gbs(fg);
  gbs.adjoin_outer_approximation(bounding_set);
  GridMaskSet<R> gis(fg);
  gis.adjoin_outer_approximation(initial_set);
  
  GridMaskSet<R> gcrs=this->chainreach(vf,gis,gbs);
  return new GridMaskSet<R>(gcrs);
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::VectorFieldEvolver<R>::viable(const System::VectorFieldInterface<R>& vector_field,
                                  const Geometry::SetInterface<R>& bounding_set) const
{
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* MapEvolver::viable(VectorFieldInterface,SetInterface)\n");
  ARIADNE_LOG(3,"bounding_set="<<bounding_set<<"\n");
  ARIADNE_CHECK_BOUNDED(bounding_set,"SetInterface* VectorFieldEvolver::viable(VectorFieldInterface vector_field, SetInterface bounding_set)");
  Rectangle<R> bounding_box=bounding_set.bounding_box();
  Grid<R> grid(bounding_set.dimension(),this->grid_length());
  GridMaskSet<R> grid_bounding_set(grid,bounding_box);
  grid_bounding_set.adjoin_outer_approximation(bounding_set);
  return new GridMaskSet<R>(this->viable(vector_field,grid_bounding_set));
}


template<class R>
tribool
Evaluation::VectorFieldEvolver<R>::verify(const System::VectorFieldInterface<R>& vector_field,
                                  const Geometry::SetInterface<R>& initial_set,
                                  const Geometry::SetInterface<R>& safe_set) const
{
  using namespace Geometry;
  const System::VectorFieldInterface<R>& vf=vector_field;
  
  Rectangle<R> bb;
  try {
    bb=safe_set.bounding_box();
  }
  catch(UnboundedSet&) {
    throw UnboundedSet("verify(Map,Set,Set): safe_set unbounded");
  }
  if(!initial_set.subset(bb)) {
    return false;
  }
  Grid<R> g(safe_set.dimension(),this->grid_length());
  FiniteGrid<R> fg(g,bb);
  GridMaskSet<R> gss(fg);
  gss.adjoin_outer_approximation(safe_set);
  GridMaskSet<R> gis(fg);
  gis.adjoin_outer_approximation(initial_set);
  
  return this->verify(vf,gis,gss);
}







template<class R>
typename Evaluation::VectorFieldEvolver<R>::BS
Evaluation::VectorFieldEvolver<R>::integrate(const VectorFieldInterface<R>& vector_field, 
                                               const BS& initial_set, 
                                               const time_type& time) const
{
  ARIADNE_LOG(4,"BasicSet VectorFieldEvolver::reach(VectorFieldInterface vector_field, BasicSet initial_set, Time time)\n");
  ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");
  
  if(time==0) { 
    return initial_set;
  }
  
  const VectorFieldInterface<R>& vf=vector_field;
  BS bs=initial_set;
  time_type t=0;
  time_type h=this->maximum_step_size();
  while(t<time) {
    h=min(time_type(time-t),h);
    bs=this->integration_step(vf,bs,h);
    t=t+h;
    h=min(time_type(2*h),this->maximum_step_size());
    h=max(h,this->minimum_step_size());
  }
  if(verbosity>4) { std::clog << "  t=" << t << "  final_set=" << bs << std::endl; }
  return bs;
}




template<class R>
Geometry::ListSet<typename Evaluation::VectorFieldEvolver<R>::BS> 
Evaluation::VectorFieldEvolver<R>::reach(const VectorFieldInterface<R>& vector_field, 
                                           const BS& initial_set, 
                                           const time_type& time) const
{
  ARIADNE_LOG(4,"ListSet<BS> VectorFieldEvolver::reach(VectorFieldInterface vector_field, BasicSet initial_set, Time time)\n");
  ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");
  
  if(time==0) { 
    return initial_set;
  }
  
  ListSet<BS> reach_set(initial_set.dimension());

  const VectorFieldInterface<R>& vf=vector_field;
  BS bs=initial_set;
  BS rbs=bs;
  time_type t=0;
  time_type h=this->maximum_step_size();
  while(t<time) {
    h=min(time_type(time-t),h);
    rbs=this->reachability_step(vf,bs,h);
    reach_set.adjoin(rbs);
    bs=this->integration_step(vf,bs,h);
    t=t+h;
    h=min(time_type(2*h),this->maximum_step_size());
    h=max(h,this->minimum_step_size());
  }
  if(verbosity>4) { std::clog << "  t=" << t << "  reach_set=" << bs << std::endl; }
  return bs;
}





// Template pattern for integrating a list set
template<class R>
Geometry::ListSet<typename Evaluation::VectorFieldEvolver<R>::BS>
Evaluation::VectorFieldEvolver<R>::lower_integrate(const VectorFieldInterface<R>& vector_field, 
                                                     const ListSet<BS>& initial_set, 
                                                     const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet< Rectangle<R> > VectorFieldEvolver::lower_integrate(VectorFieldInterface vector_field, ListSet< Rectangle<R> > initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  using namespace Numeric;
  
  if(time==0) { 
    return initial_set;
  }
  
  const VectorFieldInterface<R>& vf=vector_field;
  time_type step_size=this->maximum_step_size();
  R maximum_set_radius=this->maximum_basic_set_radius();
  
  if(verbosity>7) { std::clog << "step_size=" << step_size << "  maximum_set_radius=" << maximum_set_radius << std::endl; }
  
  time_type t=0; // t is the time elapsed!
  time_type h=step_size;
  BS bs(initial_set.dimension());
  
  typedef std::pair< time_type, BS> timed_set_type;
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets;
  
  ListSet<BS> final_set(initial_set.dimension());
  
  typedef typename ListSet<BS>::const_iterator list_set_const_iterator;
  for(list_set_const_iterator bs_iter=initial_set.begin(); bs_iter!=initial_set.end(); ++bs_iter) {
    t=0;
    bs=*bs_iter;
    h=step_size;
    
    while(t!=time && bs.radius()<=maximum_set_radius) {
      if(verbosity>5) { 
        std::clog << "      t=" << conv_approx<double>(t) << "  h=" << conv_approx<double>(h) << "  c=" << bs.centre() 
                  << "  r=" << bs.radius() << std::endl;
      }
      h=min(time_type(time-t),h);
      bs=this->integration_step(vf,bs,h);
      t=t+h;  // t is the time elapsed!
      h=min(time_type(2*h),step_size);
    } 
    
    if(verbosity>5) { 
      std::clog << "      t=" << conv_approx<double>(t) << "  c=" << bs.centre() 
                << "  r=" << bs.radius() << std::endl;
    }
    
    if(t==time) {
      final_set.adjoin(bs);
    }
  }
  if(verbosity>6) { std::clog << "  final_set=" << final_set << std::endl; }
  return final_set;
}



template<class R>
Geometry::ListSet<typename Evaluation::VectorFieldEvolver<R>::BS>
Evaluation::VectorFieldEvolver<R>::lower_reach(const VectorFieldInterface<R>& vector_field, 
                                                 const ListSet<BS>& initial_set, 
                                                 const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet< Rectangle<R> > VectorFieldEvolver::lower_reach(VectorFieldInterface vector_field, ListSet< Rectangle<R> > initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  using namespace Numeric;
  
  
  if(time==0) { 
    return initial_set;
  }
  
  const VectorFieldInterface<R>& vf=vector_field;
  time_type step_size=this->maximum_step_size();
  R maximum_set_radius=this->maximum_basic_set_radius();
  
  if(verbosity>7) { std::clog << "step_size=" << step_size << "  maximum_set_radius=" << maximum_set_radius << std::endl; }
  
  time_type t=0; // t is the time elapsed!
  time_type h=step_size;
  BS bs(initial_set.dimension());
  BS rbs(initial_set.dimension());
  
  typedef std::pair< time_type, BS> timed_set_type;
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets;
  
  ListSet<BS> final_set(initial_set.dimension());
  
  typedef typename ListSet<BS>::const_iterator list_set_const_iterator;
  for(list_set_const_iterator bs_iter=initial_set.begin(); bs_iter!=initial_set.end(); ++bs_iter) {
    t=0;
    bs=*bs_iter;
    h=step_size;
    
    while(t!=time && bs.radius()<=maximum_set_radius) {
      if(verbosity>5) { 
        std::clog << "      t=" << conv_approx<double>(t) << "  h=" << conv_approx<double>(h) << "  c=" << bs.centre() 
                  << "  r=" << bs.radius() << std::endl;
      }
      h=min(time_type(time-t),h);
      rbs=this->reachability_step(vf,bs,h);
      bs=this->integration_step(vf,bs,h);
      final_set.adjoin(rbs);
      t=t+h;  // t is the time elapsed!
      h=min(time_type(2*h),step_size);
    } 
    if(verbosity>5) { 
      std::clog << "      t=" << conv_approx<double>(t) << "  c=" << bs.centre() 
                << "  r=" << bs.radius() << std::endl;
    }
  }
  if(verbosity>6) { std::clog << "  final_set=" << final_set << std::endl; }
  return final_set;
}



// Template pattern for integrating a list set
template<class R>
Geometry::ListSet<typename Evaluation::VectorFieldEvolver<R>::BS>
Evaluation::VectorFieldEvolver<R>::upper_integrate(const VectorFieldInterface<R>& vector_field, 
                                                     const ListSet<BS>& initial_set, 
                                                     const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet< Rectangle<R> > VectorFieldEvolver::upper_integrate(VectorFieldInterface vector_field, ListSet< Rectangle<R> > initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  using namespace Numeric;
  
  if(verbosity>4) { std::clog << "ListSet< Rectangle<R> > VectorFieldEvolver::integrate(VectorFieldInterface,ListSet< Rectangle<R> >,Time)" << std::endl; } 
  
  if(time==0) { 
    return initial_set;
  }
  
  const VectorFieldInterface<R>& vf=vector_field;
  time_type step_size=this->maximum_step_size();
  R maximum_set_radius=this->maximum_basic_set_radius();
  
  if(verbosity>7) { std::clog << "step_size=" << step_size << "  maximum_set_radius=" << maximum_set_radius << std::endl; }
  
  time_type t=0; // t is the time elapsed!
  time_type h=step_size;
  BS bs(initial_set.dimension());
  
  typedef std::pair< time_type, BS> timed_set_type;
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets;
  
  ListSet<BS> final_set(initial_set.dimension());
  
  typedef typename ListSet<BS>::const_iterator list_set_const_iterator;
  for(list_set_const_iterator bs_iter=initial_set.begin(); bs_iter!=initial_set.end(); ++bs_iter) {
    working_sets.push_back(timed_set_type(0,*bs_iter));
  }
  
  while(!working_sets.empty()) {
    if(verbosity>7) { std::clog << "working_sets.size()=" << working_sets.size() << "\n"; }
    
    const timed_set_type& ts=working_sets.back();
    t=ts.first;
    bs=ts.second;
    h=step_size;
    working_sets.pop_back();
    
    if(verbosity>5) { std::clog << "  t=" << t << "  bs=" << bs << std::endl; }
    if(bs.radius()>maximum_set_radius) {
      if(verbosity>5) { std::clog << "    subdividing..." << std::flush; }
      ListSet<BS> subdivisions=this->subdivide(bs);
      if(verbosity>5) { std::clog << "    subdivisions.size() =" << subdivisions.size() << std::endl; }
      for(list_set_const_iterator subdiv_iter=subdivisions.begin(); 
          subdiv_iter!=subdivisions.end(); ++subdiv_iter)
        {
          working_sets.push_back(timed_set_type(t,*subdiv_iter));
        }
      if(verbosity>5) { std::clog << " done" << std::endl; }
    }
    else {
      if(verbosity>5) { std::clog << "    integrating..." << std::endl; }
      do {
        if(verbosity>5) { 
          std::clog << "      t=" << conv_approx<double>(t) << "  h=" << conv_approx<double>(h) << "  c=" << bs.centre() 
                    << "  r=" << bs.radius() << std::endl;
        }
        h=min(time_type(time-t),h);
        bs=this->integration_step(vf,bs,h);
        t=t+h;  // t is the time remaining!
        h=min(time_type(2*h),step_size);
      } while(t!=time && bs.radius()<=maximum_set_radius);
      
      if(verbosity>5) { 
        std::clog << "      t=" << conv_approx<double>(t) << "  c=" << bs.centre() 
                  << "  r=" << bs.radius() << std::endl;
      }
      
      if(t==time) {
        final_set.adjoin(bs);
      } else {
        working_sets.push_back(timed_set_type(t,bs));
      }
    }
  }
  if(verbosity>6) { std::clog << "  final_set=" << final_set << std::endl; }
  return final_set;
}



template<class R>
Geometry::ListSet<typename Evaluation::VectorFieldEvolver<R>::BS>
Evaluation::VectorFieldEvolver<R>::upper_reach(const VectorFieldInterface<R>& vector_field, 
                                                 const ListSet<BS>& initial_set, 
                                                 const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet<BasicSet> VectorFieldEvolver::integrate(VectorFieldInterface vector_field, ListSet<BasicSet> initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  using namespace Numeric;
  
  const VectorFieldInterface<R>& vf=vector_field;
  time_type step_size=this->maximum_step_size();
  R maximum_set_radius=this->maximum_basic_set_radius();
  
  if(verbosity>4) {
    std::clog << "step_size=" << conv_approx<double>(step_size) << "  maximum_set_radius()=" << maximum_set_radius << std::endl<<std::flush;
  }
  
  time_type t=0;
  time_type h=step_size;
  BS bs(initial_set.dimension());
  BS rs(initial_set.dimension());
  
  typedef typename ListSet<BS>::const_iterator list_set_const_iterator;
  typedef typename ListSet<BS>::iterator list_set_iterator;
  typedef std::pair< time_type, BS > timed_set_type;
  
  std::vector< timed_set_type > working_sets;
  ListSet<BS> reach_set(initial_set.dimension());
  
  for(list_set_const_iterator bs_iter=initial_set.begin(); bs_iter!=initial_set.end(); ++bs_iter) {
    working_sets.push_back(timed_set_type(0,*bs_iter));
  }
  
  if(verbosity>6) { std::clog << "initial_set.size()=" << initial_set.size() << std::endl; }
  assert(working_sets.size()==initial_set.size());
  assert(reach_set.size()==0);
  
  while(!working_sets.empty()) {
    if(verbosity>6) { 
      std::clog << "  working_sets.size()=" << working_sets.size() << "  reach_set.size()=" << reach_set.size() << std::endl; 
    }
    
    const timed_set_type& ts=working_sets.back();
    t=ts.first;
    bs=ts.second;
    h=step_size;
    working_sets.pop_back();
    
    if(verbosity>6) { std::clog << "  t=" << conv_approx<double>(t) << "  bs.centre()=" << bs.centre() << "  bs.radius() = " << bs.radius() << std::endl; }
    
    if(bs.radius()>maximum_set_radius) {
      if(verbosity>5) { std::clog << "  subdividing..." << std::endl; }
      ListSet<BS> subdivisions=this->subdivide(bs);
      if(verbosity>5) { std::clog << "    subdivisions.size() =" << subdivisions.size() << std::endl; }
      for(list_set_const_iterator subdiv_iter=subdivisions.begin(); 
          subdiv_iter!=subdivisions.end(); ++subdiv_iter)
        {
          working_sets.push_back(timed_set_type(t,*subdiv_iter));
        }
    }
    else {
      if(verbosity>6) { std::clog << "  integrating..." << std::endl; }
      
      do {
        if(verbosity>6) {
          std::clog << "    t=" << conv_approx<double>(t) << "  h=" << conv_approx<double>(h)
                    << "  bs=" << bs << std::endl;
        }
        
        h=min(time_type(time-t),h);
        rs=this->reachability_step(vf,bs,h);
        reach_set.adjoin(rs);
        
        if(t<time) {
          bs=integration_step(vf,bs,h);
          t=t+h;
        }
      } while(t!=time && bs.radius()<=maximum_set_radius);
      
      if(t<time) {
        working_sets.push_back(timed_set_type(t,bs));
      }
    }
  }
  
  if(verbosity>4) {
    std::clog << "  reach_set.size()=" << reach_set.size() <<  std::endl;
  }
  return reach_set;
} 







template<class R>
Geometry::ListSet< Rectangle<R> >
Evaluation::VectorFieldEvolver<R>::integrate(const System::VectorFieldInterface<R>& vector_field, 
                                               const Geometry::ListSet< Rectangle<R> >& initial_set,
                                               const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet<Rectangle> VectorFieldEvolver::integrate(VectorFieldInterface vector_field, ListSet<Rectangle> initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  ListSet< Rectangle<R> > result;
  const VectorFieldInterface<R>& vf=vector_field;
  ListSet<BS> ils=initial_set;
  ListSet<BS> fls=this->lower_integrate(vf,ils,time);
  for(typename ListSet<BS>::const_iterator bs_iter=fls.begin();
      bs_iter!=fls.end(); ++bs_iter)
  {
    result.adjoin(bs_iter->bounding_box());
  }
  return result;
}



template<class R>
Geometry::ListSet< Rectangle<R> >
Evaluation::VectorFieldEvolver<R>::reach(const System::VectorFieldInterface<R>& vector_field, 
                                           const Geometry::ListSet< Rectangle<R> >& initial_set,
                                           const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet<Rectangle> VectorFieldEvolver::integrate(VectorFieldInterface vector_field, ListSet<Rectangle> initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");

  ListSet< Rectangle<R> > result;
  const VectorFieldInterface<R>& vf=vector_field;
  ListSet<BS> ils=initial_set;
  ListSet<BS> rls=this->lower_reach(vf,ils,time);
  for(typename ListSet<BS>::const_iterator bs_iter=rls.begin();
      bs_iter!=rls.end(); ++bs_iter)
  {
    result.adjoin(bs_iter->bounding_box());
  }
  return result;
}



template<class R>
Geometry::GridMaskSet<R>
Evaluation::VectorFieldEvolver<R>::integrate(const System::VectorFieldInterface<R>& vector_field, 
                                               const Geometry::GridMaskSet<R>& initial_set,
                                               const Geometry::GridMaskSet<R>& bounding_set,
                                               const time_type& time) const
{
  ARIADNE_LOG(2,"GridMaskSet VectorFieldEvolver::integrate(VectorFieldInterface vector_field, GridMaskSet initial_set, GridMaskSet bounding_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  const VectorFieldInterface<R>& vf=vector_field;
  
  ARIADNE_CHECK_SAME_GRID(initial_set,bounding_set,"GridMaskSet VectorFieldEvolver::integrate(VectorFieldInterface,GridMaskSet,GridMaskSet,Time)");
  using namespace System;
  using namespace Geometry;
  using namespace LinearAlgebra;
  
  if(time==0) { 
    return initial_set;
  }
  
  time_type step_size=this->maximum_step_size();
  
  Rectangle<R> bb=bounding_set.bounding_box();
  
  GridMaskSet<R> result(bounding_set);
  result.clear();
  
  Rectangle<R> tmp_rectangle;
  BS tmp_basic_set;
  
  time_type t=time;
  time_type h=step_size;
  
  //bb=regular_intersection(bb,bounding_box());
  
  R spacial_tolerance=2;
  
  ListSet<BS> start_set;
  ListSet<BS> finish_set;
  for(typename GridMaskSet<R>::const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
    tmp_rectangle=*iter;
    tmp_basic_set=tmp_rectangle;
    start_set.adjoin(tmp_basic_set);
  }
  
  while(t!=0) {
    if(verbosity > 6) { std::clog << "time left=" << t << "  stepsize=" << h << "  sets in list=" << start_set.size() << "\n"; }
    h=min(t,h);
    for(typename ListSet<BS>::const_iterator iter=start_set.begin(); iter!=start_set.end(); ++iter) {
      BS p(*iter);
      p=this->integrate(vf,p,h);
      finish_set.adjoin(p);
    }
    start_set.clear();
    GridMaskSet<R> mask_set(bounding_set);
    mask_set.clear();
    for(typename ListSet<BS>::const_iterator iter=finish_set.begin(); iter!=finish_set.end(); ++iter) {
      const BS& p=*iter;
      if(p.radius()>spacial_tolerance) {
        if(verbosity > 6) { std::clog << "Splitting, radius=" << p.radius() << "\n" << p << "\n"; }
        mask_set.adjoin_outer_approximation(p);
      }
      else {
        start_set.adjoin(p);
      }
    }
    for(typename GridMaskSet<R>::const_iterator iter=mask_set.begin(); iter!=mask_set.end(); ++iter) {
      tmp_rectangle=*iter;
      tmp_basic_set=tmp_rectangle;
      start_set.adjoin(tmp_basic_set);
    }
    finish_set.clear();
    t-=h;
  }
  
  for(typename ListSet<BS>::const_iterator iter=start_set.begin(); iter!=start_set.end(); ++iter) {
    const BS& bs=*iter;
    result.adjoin_outer_approximation(bs);
  }
  return result;
}





template<class R>
Geometry::GridMaskSet<R>
Evaluation::VectorFieldEvolver<R>::reach(const System::VectorFieldInterface<R>& vector_field, 
                                           const Geometry::GridMaskSet<R>& initial_set,
                                           const Geometry::GridMaskSet<R>& bounding_set,
                                           const time_type& time) const
{
  using namespace Numeric;
  
  ARIADNE_LOG(2,"GridMaskSet VectorFieldEvolver::reach(VectorFieldInterface vector_field, GridMaskSet initial_set, GridMaskSet bounding_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  typedef typename GridMaskSet<R>::const_iterator gms_const_iterator;
  typedef typename ListSet<BS>::const_iterator ls_const_iterator;
  ARIADNE_CHECK_BOUNDED(initial_set,"GridMaskSet VectorFieldEvolver::reach(VectorFieldInterface,GridMaskSet,GridMaskSet,Time)");
  ARIADNE_CHECK_BOUNDED(bounding_set,"GridMaskSet VectorFieldEvolver::reach(VectorFieldInterface,GridMaskSet,GridMaskSet,Time)");
  const VectorFieldInterface<R>& vf=vector_field;
  
  if(!subset(initial_set,bounding_set)) {
    throw std::runtime_error("GridMaskSet VectorFieldEvolver::reach(VectorFieldInterface,GridMaskSet,GridMaskSet,Time): Initial set must be subset of bounding set");
  }
  
  const GridMaskSet<R>& is=initial_set;
  const Rectangle<R> bb=bounding_set.bounding_box();
  
  GridMaskSet<R> result(bounding_set);
  result.clear();
  
  GridMaskSet<R> stored(result);
  GridMaskSet<R> found(result);
  GridMaskSet<R> image(result);
  found.adjoin(is);
  
  int steps=int_up<int>(time_type(time/this->lock_to_grid_time()));
  if (steps==0) { steps=1; }
  
  time_type time_step=time/steps;
  
  for(int step=0; step!=steps; ++step) {
    found=difference(found,stored);
    stored.adjoin(found);
    image.clear();
    GridMaskSet<R> image=this->integrate(vf,found,bounding_set,time_step);
    found=image;
  }
  
  ListSet<BS> input_list;
  ListSet<BS> output_list;
  for(gms_const_iterator iter=stored.begin(); iter!=stored.end(); ++iter) {
    input_list.adjoin(BS(Rectangle<R>(*iter)));
  }
  output_list=this->upper_reach(vf,input_list,time_step);
  for(ls_const_iterator iter=output_list.begin(); iter!=output_list.end(); ++iter) {
    const BS& fz=*iter;
    result.adjoin_outer_approximation(fz);
  }
  return result;
}



template<class R>
Geometry::GridMaskSet<R>
Evaluation::VectorFieldEvolver<R>::chainreach(const System::VectorFieldInterface<R>& vector_field, 
                                                const Geometry::GridMaskSet<R>& initial_set, 
                                                const Geometry::GridMaskSet<R>& bounding_set) const
{
  ARIADNE_LOG(2,"GridMaskSet VectorFieldEvolver::chainreach(VectorFieldInterface vector_field, GridMaskSet initial_set, GridMaskSet bounding_set)\n");
  ARIADNE_LOG(5,"vector_field="<<vector_field<<"\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  typedef typename GridCellListSet<R>::const_iterator gcls_const_iterator;
  typedef typename GridMaskSet<R>::const_iterator gms_const_iterator;
  typedef typename ListSet<BS>::const_iterator ls_const_iterator;
  ARIADNE_CHECK_BOUNDED(initial_set,"GridMaskSet VectorFieldEvolver::chainreach(VectorFieldInterface,GridMaskSet,GridMaskSet)");
  ARIADNE_CHECK_BOUNDED(bounding_set,"GridMaskSet VectorFieldEvolver::chainreach(VectorFieldInterface,GridMaskSet,GridMaskSet)");
  const VectorFieldInterface<R>& vf=vector_field;
  
  if(!subset(initial_set,bounding_set)) {
    throw std::runtime_error("GridMaskSet VectorFieldEvolver::chainreach(VectorFieldInterface,GridMaskSet,GridMaskSet): Initial set must be subset of bounding set");
  }
  
  const GridMaskSet<R>& is=initial_set;
  const Rectangle<R> bb=bounding_set.bounding_box();
  
  GridMaskSet<R> result(initial_set);
  result.clear();
  GridCellListSet<R> image(initial_set.grid());
  GridCellListSet<R> found(initial_set.grid());
  found.adjoin(is);
  
  time_type step_size=this->maximum_step_size();
  time_type time_step=this->lock_to_grid_time();
  
  if(verbosity>4) { std::clog << "Beginning integration phase" << std::endl; }
  while(!subset(found,result)) {
    if(verbosity>5) { std::clog << "Found " << found.size() << " cells, " << std::flush; }
    found.unique_sort();
    if(verbosity>5) { std::clog << "of which " << found.size() << " are not duplicates," << std::flush; }
    found=difference(found,result);
    if(verbosity>5) { std::clog << " and " << found.size() << " are new" << std::endl; }
    result.adjoin(found);
    image.clear();
    uint size=0;
    ListSet<BS> basic_set_list;
    for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
      ++size;
      Rectangle<R> r=*iter;
      BS z(r);
      basic_set_list.adjoin(z);
    }
    basic_set_list=this->upper_integrate(vf,basic_set_list,time_step);
    if(verbosity>5) { std::clog << "new basic sets " << basic_set_list << std::endl; }
    for(ls_const_iterator iter=basic_set_list.begin(); iter!=basic_set_list.end(); ++iter) {
      const BS& fp=*iter;
      if(!disjoint(fp.bounding_box(),bounding_set)) {
        image.adjoin_outer_approximation(fp);
      } 
    }
    if(verbosity>5) { std::clog << "image set " << image << std::endl; }
    found=regular_intersection(image,bounding_set);
  }
  if(verbosity>5) { std::clog << "Found " << result.size() << " cells, " << std::endl; }
  
  if(verbosity>4) { std::clog << "Beginning reachability phase" << std::endl; }
  ListSet<BS> reach_basic_set_list;
  for(gms_const_iterator iter=result.begin(); iter!=result.end(); ++iter) {
    Rectangle<R> r=*iter;
    BS z(r);
    reach_basic_set_list.adjoin(z);
  }
  reach_basic_set_list=this->upper_reach(vf,reach_basic_set_list,time_step);
  for(ls_const_iterator iter=reach_basic_set_list.begin(); iter!=reach_basic_set_list.end(); ++iter) {
    BS fz=*iter;
    if(!disjoint(fz.bounding_box(),bounding_set)) {
      result.adjoin_outer_approximation(fz);
    }
  }
  if(verbosity>4) { std::clog << "Reached " << result.size() << " cells, " << std::endl; }
  
  result.restrict(bounding_set);
  
  return result;
}


template<class R>
Geometry::GridMaskSet<R>
Evaluation::VectorFieldEvolver<R>::viable(const System::VectorFieldInterface<R>& vector_field, 
                                            const Geometry::GridMaskSet<R>& bounding_set) const
{
  ARIADNE_LOG(2,"GridMaskSet VectorFieldEvolver::viable(VectorFieldInterface vector_field, GridMaskSet bounding_set)\n");
  ARIADNE_LOG(3,"bounding_set="<<bounding_set<<"\n");
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
  ARIADNE_CHECK_BOUNDED(bounding_set,"VectorFieldEvolver<R>::viable(VectorFieldInterface,GridMaskSet)");
  
  const VectorFieldInterface<R>& vf=vector_field;
  const Grid<R>& g=bounding_set.grid();
  Combinatoric::LatticeBlock bd=bounding_set.block();
  GridBlock<R> bb(g,bd);
  GridMaskSet<R> result(g,bd);
  GridCellListSet<R> unsafe=bounding_set;
  
  Rectangle<R> r(g.dimension());
  Rectangle<R> fr(g.dimension());
  GridBlock<R> fgb(g);
  BS bs(g.dimension());
  BS fbs(g.dimension());
  GridCellListSet<R> fgcls(g);
  
  // FIXME: Reduct integration time properly
  time_type h=1;
  result=bounding_set;
  size_type step=0;
  size_type moved=1;
  while(moved>0) {
    moved=0;
    while(!unsafe.empty()) {
      ARIADNE_LOG(3,"Viability step "<<step+1<< ": testing "<<result.size()<<" cells with an integration time of "<<h<<"\n");
      ++step;
      unsafe.clear();
      for(gms_const_iterator iter=result.begin(); iter!=result.end(); ++iter) {
        fgcls.clear();
        r=*iter;
        bs=r;
        fbs=this->integration_step(vf,bs,h);
        if(disjoint(fbs.bounding_box(),r)) {
          ++moved;
        }
        fgcls.adjoin_outer_approximation(fbs);
        if(!overlap(result,outer_approximation(fbs,g))) {
          unsafe.adjoin(*iter);
        }
      }
      result.remove(unsafe);
    }
    h/=2;
  }
  return result;
}


template<class R>
tribool
Evaluation::VectorFieldEvolver<R>::verify(const System::VectorFieldInterface<R>& vector_field, 
                                            const Geometry::GridMaskSet<R>& initial_set, 
                                            const Geometry::GridMaskSet<R>& safe_set) const
{
  ARIADNE_LOG(2,"GridMaskSet VectorFieldEvolver::verify(VectorFieldInterface vector_field, GridMaskSet initial_set, GridMaskSet safe_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  ARIADNE_LOG(3,"safe_set="<<safe_set<<"\n");

  typedef typename GridCellListSet<R>::const_iterator gcls_const_iterator;
  typedef typename ListSet<BS>::const_iterator ls_const_iterator;
  ARIADNE_CHECK_BOUNDED(initial_set,"tribool VectorFieldEvolver::verify(VectorFieldInterface,GridMaskSet,GridMaskSet)");
  ARIADNE_CHECK_BOUNDED(safe_set,"tribool VectorFieldEvolver::verify(VectorFieldInterface,GridMaskSet,GridMaskSet)");
  const VectorFieldInterface<R>& vf=vector_field;
  
  if(!subset(initial_set,safe_set)) {
    return false;
  }
  
  const GridMaskSet<R>& is=initial_set;
  const Rectangle<R> bb=safe_set.bounding_box();
  
  GridMaskSet<R> chainreach(is);
  GridCellListSet<R> found(is.grid());
  GridCellListSet<R> image(is.grid());
  GridCellListSet<R> cellimage(is.grid());
  found.adjoin(is);
  
  time_type step_size=this->maximum_step_size();
  time_type time_step=this->lock_to_grid_time();
  
  while(!subset(found,chainreach)) {
    found=difference(found,chainreach);
    chainreach.adjoin(found);
    image.clear();
    uint size=0;
    ListSet<BS> basic_set_list;
    for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
      ++size;
      Rectangle<R> r=*iter;
      BS pp(r);
      basic_set_list.adjoin(pp);
    }
    basic_set_list=this->upper_integrate(vf,basic_set_list,time_step);
    for(ls_const_iterator iter=basic_set_list.begin(); iter!=basic_set_list.end(); ++iter) {
      const BS& fp=*iter;
      cellimage.adjoin_outer_approximation(fp);
      if(!subset(cellimage,safe_set)) {
        return false;
      }
      image.adjoin(cellimage);
      cellimage.clear();
    }
    found=image;
  }
  return indeterminate;
}



}
