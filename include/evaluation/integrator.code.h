/***************************************************************************
 *            integrator.code.h
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
 
#include <iosfwd>
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

#include "../output/logging.h"

#include "integrator.h"

namespace {
    
using namespace Ariadne;

template<class R, class BST>
struct pair_first_less {
  bool operator()(const std::pair<R,BST>& ts1, const std::pair<R,BST>& ts2) {
    return ts1.first < ts2.first; 
  }
};

/*!\ brief A class representing pre-computed bounds for an integration step. */
template<class R>
class IntegrationStepBound {
 public:
  /*!\ brief Constructor. */
  IntegrationStepBound(const Geometry::Rectangle<R>& bound, const time_type& integration_time) 
    : _bound(bound), _integration_time(integration_time) { }
  /*! The spacial bound for the integrations step. */
  const Geometry::Rectangle<R>& bound() const { return _bound; }
  /*! The step size in time of the integrations step. */
  const R& integration_time() const { return _integration_time; }
 private:
  Geometry::Rectangle<R> _bound;
  R _integration_time;
};

}


namespace Ariadne {
    
namespace Evaluation { static int& verbosity = integrator_verbosity; }
    


template<class R>
Evaluation::Integrator<R>::~Integrator()
{
}



template<class R>
Evaluation::Integrator<R>::Integrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
  : _maximum_step_size(maximum_step_size),
    _lock_to_grid_time(lock_to_grid_time),
    _grid_size(0.125),
    _maximum_basic_set_radius(maximum_basic_set_radius)
{
}



template<class R>
time_type Evaluation::Integrator<R>::minimum_step_size() const
{
  return this->_minimum_step_size;
}

template<class R>
time_type Evaluation::Integrator<R>::maximum_step_size() const
{
  return this->_maximum_step_size;
}

template<class R>
R Evaluation::Integrator<R>::minimum_basic_set_radius() const
{
  return this->_minimum_basic_set_radius;
}

template<class R>
R Evaluation::Integrator<R>::maximum_basic_set_radius() const
{
  return this->_maximum_basic_set_radius;
}

template<class R>
R Evaluation::Integrator<R>::grid_size() const
{
  return this->_grid_size;
}

template<class R>
void Evaluation::Integrator<R>::set_grid_size(const R& grid_size)
{
  this->_grid_size=grid_size;
}

template<class R>
time_type Evaluation::Integrator<R>::lock_to_grid_time() const
{
  return this->_lock_to_grid_time;
}



template<class R, class VF, class BS>
Geometry::SetInterface<R>*
Evaluation::IntegratorBase<R,VF,BS>::integrate(const System::VectorField<R>& vector_field, 
                                               const Geometry::SetInterface<R>& initial_set, 
                                               const time_type& time) const
{
  typedef Numeric::Interval<R> I;
  ARIADNE_LOG(2,"SetInterface* Integrator::integrate(VectorField,SetInterface,Time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  const VectorField& cast_vector_field=dynamic_cast<const VectorField&>(vector_field);
  Rectangle bb=initial_set.bounding_box();
  Grid grid(vector_field.dimension(),this->grid_size());
  RectangleListSet rectangle_list_initial_set=lower_approximation(initial_set,grid);
  ARIADNE_LOG(3,"rectangle_list_initial_set="<<rectangle_list_initial_set);
  ListSet list_initial_set(rectangle_list_initial_set);
  ARIADNE_LOG(3,"list_initial_set="<<list_initial_set);
  ListSet list_final_set=this->lower_integrate(cast_vector_field,list_initial_set,time);
  ARIADNE_LOG(3,"list_final_set="<<list_final_set);
  return new ListSet(list_final_set);
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::Integrator<R>::integrate(const System::VectorField<R>& vector_field,
                                     const Geometry::SetInterface<R>& initial_set,
                                     const Geometry::SetInterface<R>& bounding_set,
                                     const time_type& time) const
{
  // FIXME: Only computes over-approximation;
  using namespace Geometry;
  const System::VectorField<R>& vf=vector_field;
  
  Rectangle<R> bb;
  try {
    bb=bounding_set.bounding_box();
  }
  catch(UnboundedSet&) {
    throw UnboundedSet("integrate(Map,Set,Set,Time): bounding_set unbounded");
  }
  Grid<R> g(bounding_set.dimension(),this->grid_size());
  FiniteGrid<R> fg(g,bb);
  GridMaskSet<R> gbs(fg);
  gbs.adjoin_outer_approximation(bounding_set);
  GridMaskSet<R> gis(fg);
  gis.adjoin_outer_approximation(initial_set);
  
  GridMaskSet<R> gs=this->reach(vf,gis,gbs,time);
  return new GridMaskSet<R>(gs);
}


template<class R, class VF, class BS>
Geometry::SetInterface<R>*
Evaluation::IntegratorBase<R,VF,BS>::reach(const System::VectorField<R>& vector_field,
                                           const Geometry::SetInterface<R>& initial_set,
                                           const time_type& time) const
{
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  ARIADNE_LOG(2,"SetInterface* Integrator::reach(VectorFile,SetInterface,Time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  const VectorField& cast_vector_field=dynamic_cast<const VectorField&>(vector_field);
  Rectangle bb=initial_set.bounding_box();
  Grid grid(vector_field.dimension(),this->grid_size());
  RectangleListSet rectangle_list_initial_set=lower_approximation(initial_set,grid);
  ListSet list_initial_set(rectangle_list_initial_set);
  ListSet list_reach_set=this->lower_reach(cast_vector_field,list_initial_set,time);
  return new ListSet(list_reach_set);
}


template<class R>
Geometry::SetInterface<R>*
Evaluation::Integrator<R>::reach(const System::VectorField<R>& vector_field,
                                 const Geometry::SetInterface<R>& initial_set,
                                 const Geometry::SetInterface<R>& bounding_set,
                                 const time_type& time) const
{
  using namespace Geometry;
  const System::VectorField<R>& vf=vector_field;
  
  Rectangle<R> bb;
  try {
    bb=bounding_set.bounding_box();
  }
  catch(UnboundedSet&) {
    throw UnboundedSet("reach(Map,Set,Set,Time): bounding_set unbounded");
  }
  Grid<R> g(bounding_set.dimension(),this->grid_size());
  FiniteGrid<R> fg(g,bb);
  GridMaskSet<R> gbs(fg);
  gbs.adjoin_outer_approximation(bounding_set);
  GridMaskSet<R> gis(fg);
  gis.adjoin_outer_approximation(initial_set);
  
  GridMaskSet<R> grs=this->reach(vf,gis,gbs,time);
  return new GridMaskSet<R>(grs);
}


template<class R, class VF, class BS>
Geometry::SetInterface<R>*
Evaluation::IntegratorBase<R,VF,BS>::reach(const System::VectorField<R>& vector_field,
                                           const Geometry::SetInterface<R>& initial_set) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
Geometry::SetInterface<R>*
Evaluation::Integrator<R>::chainreach(const System::VectorField<R>& vector_field,
                                      const Geometry::SetInterface<R>& initial_set,
                                      const Geometry::SetInterface<R>& bounding_set) const
{
  using namespace Geometry;
  const System::VectorField<R>& vf=vector_field;
  
  Rectangle<R> bb;
  try {
    bb=bounding_set.bounding_box();
  }
  catch(UnboundedSet&) {
    throw UnboundedSet("chainreach(Map,Set,Set): bounding_set unbounded");
  }
  Grid<R> g(bounding_set.dimension(),this->grid_size());
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
Evaluation::Integrator<R>::viable(const System::VectorField<R>& vector_field,
                                  const Geometry::SetInterface<R>& bounding_set) const
{
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* Applicator::viable(VectorField,SetInterface)\n");
  ARIADNE_LOG(3,"bounding_set="<<bounding_set<<"\n");
  ARIADNE_CHECK_BOUNDED(bounding_set,"SetInterface* Integrator::viable(VectorField vector_field, SetInterface bounding_set)");
  Rectangle<R> bounding_box=bounding_set.bounding_box();
  Grid<R> grid(bounding_set.dimension(),this->grid_size());
  GridMaskSet<R> grid_bounding_set(grid,bounding_box);
  grid_bounding_set.adjoin_outer_approximation(bounding_set);
  return new GridMaskSet<R>(this->viable(vector_field,grid_bounding_set));
}


template<class R>
tribool
Evaluation::Integrator<R>::verify(const System::VectorField<R>& vector_field,
                                  const Geometry::SetInterface<R>& initial_set,
                                  const Geometry::SetInterface<R>& safe_set) const
{
  using namespace Geometry;
  const System::VectorField<R>& vf=vector_field;
  
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
  Grid<R> g(safe_set.dimension(),this->grid_size());
  FiniteGrid<R> fg(g,bb);
  GridMaskSet<R> gss(fg);
  gss.adjoin_outer_approximation(safe_set);
  GridMaskSet<R> gis(fg);
  gis.adjoin_outer_approximation(initial_set);
  
  return this->verify(vf,gis,gss);
}





template<class R>
bool
Evaluation::Integrator<R>::check_flow_bounds(const System::VectorField<R>& vf,
                                             const Geometry::Rectangle<R>& r,
                                             const Geometry::Rectangle<R>& b,
                                             const time_type& h) const
{
  if(verbosity>6) { std::clog << "Integrator::check_flow_bounds" << std::endl; }
  using namespace Geometry;
  using namespace Numeric;
  return subset(r+Interval<R>(0,h)*vf(b),b);
}



template<class R>
Geometry::Rectangle<R>
Evaluation::Integrator<R>::estimate_flow_bounds(const System::VectorField<R>& vf,
                                                const Geometry::Rectangle<R>& r,
                                                const time_type& h,
                                                const unsigned int& maximum_iterations) const
{
  using namespace Geometry;
  using namespace Numeric;
  
  if(verbosity>6) { std::clog << "Integrator::estimate_flow_bounds" << " (maximum_iterations=" << maximum_iterations << ")" << std::endl; }
  if(verbosity>7) { std::clog << "  h=" << conv_approx<double>(h) << "  r=" << r << "  vf(r)=" << vf(r) << std::endl; }
  
  typedef typename Numeric::traits<R>::arithmetic_type F;
  uint iteration=0;
  R multiplier=1.125;
  time_type t=h;
  Rectangle<R> reach(vf.dimension());
  Rectangle<R> bounds(vf.dimension());
  reach=r;
  
  while(t>0) {
    bounds=reach+multiplier*Interval<R>(R(0),h)*vf(reach);
    LinearAlgebra::Vector< Interval<R> > df=vf(bounds);
    
    time_type dt=t;
    for(dimension_type i=0; i!=vf.dimension(); ++i) {
      if(df(i).upper()>0) {
        dt=min(dt,time_type(div_up(sub_up(bounds[i].upper(),reach[i].upper()),df(i).upper())));
      }
      if(df(i).lower()<0) {
        dt=min(dt,time_type(div_up(sub_up(bounds[i].lower(),reach[i].lower()),df(i).lower())));
      }
    }
    reach=bounds;
    t-=dt;
    
    if(verbosity>8) { std::clog << "t=" << conv_approx<double>(t) << "  reach="  << reach << std::endl; }
    
    ++iteration;
    if(iteration==maximum_iterations) {
      throw std::runtime_error(std::string(__FUNCTION__)+": Cannot find bounding box for flow");
    }
  }
  return reach;
}



template<class R>
Geometry::Rectangle<R>
Evaluation::Integrator<R>::estimate_flow_bounds(const System::VectorField<R>& vf,
                                                const Geometry::Rectangle<R>& r,
                                                time_type& h) const
{
  using namespace Geometry;
  using namespace Numeric;
  
  if(verbosity>6) { std::clog << "Integrator::estimate_flow_bounds" << std::endl; }
  
  static const unsigned int max_tries=8;
  
  unsigned int max_iterations=4;
  unsigned int remaining_tries=max_tries;
  
  Rectangle<R> bounds(vf.dimension());
  while(bounds.empty()) {
    try {
      bounds=estimate_flow_bounds(vf,r,h,max_iterations);
    }
    catch(std::runtime_error) { 
      h/=2;
      max_iterations+=1;
      --remaining_tries;
      if(remaining_tries==0) {
        throw std::runtime_error(std::string(__FUNCTION__)+": cannnot find bounding box for flow");
      }
    }
  }
  
  if(verbosity>7) { std::clog << "  h=" << conv_approx<double>(h) << "  b=" << bounds << std::endl; }
  
  return bounds;
}



template<class R>
Geometry::Rectangle<R>
Evaluation::Integrator<R>::refine_flow_bounds(const System::VectorField<R>& vector_field,
                                              const Geometry::Rectangle<R>& initial_set,
                                              const Geometry::Rectangle<R>& estimated_bounds,
                                              const time_type& step_size) const
{
  if(verbosity>6) { std::clog << "Integrator::refine_flow_bounds" << std::endl; }
  
  using namespace System;
  using namespace Geometry;
  using namespace LinearAlgebra;
  using namespace Numeric;
  const VectorField<R>& vf=vector_field;
  Rectangle<R> rx=initial_set;
  Rectangle<R> b=estimated_bounds;
  Numeric::Interval<R> h=step_size;
  
  Rectangle<R> xb=rx+Numeric::Interval<R>(0,step_size)*vf(b);
  Rectangle<R> xxb=rx+Numeric::Interval<R>(0,step_size)*vf(xb);
  
  if(verbosity>7) { std::clog << "new_bounds " << xxb << "," << xb << " vs old_bounds " << b << "  " << subset(xb,b) << std::endl; }
  
  Vector< Interval<R> > ddphi=vf.jacobian(xb)*vf(xb);
  Vector< Interval<R> > dfx=vf(rx);
  Vector< Interval<R> > hdfx=h*dfx;
  Vector< Interval<R> > hhddphi=(h*h/R(2))*ddphi;
  Vector< Interval<R> > dx=hdfx+hhddphi;
  return rx+dx;
}






template<class R, class VF, class BS>
typename Evaluation::IntegratorBase<R,VF,BS>::BasicSet 
Evaluation::IntegratorBase<R,VF,BS>::integrate(const VectorField& vector_field, 
                                               const BasicSet& initial_set, 
                                               const time_type& time) const
{
  ARIADNE_LOG(4,"BasicSet Integrator::reach(VectorField vector_field, BasicSet initial_set, Time time)\n");
  ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");
  
  if(time==0) { 
    return initial_set;
  }
  
  const VF& vf=vector_field;
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




template<class R, class VF, class BS>
typename Evaluation::IntegratorBase<R,VF,BS>::ListSet 
Evaluation::IntegratorBase<R,VF,BS>::reach(const VectorField& vector_field, 
                                           const BasicSet& initial_set, 
                                           const time_type& time) const
{
  ARIADNE_LOG(4,"ListSet Integrator::reach(VectorField vector_field, BasicSet initial_set, Time time)\n");
  ARIADNE_LOG(5,"initial_set="<<initial_set<<"\n");
  
  if(time==0) { 
    return initial_set;
  }
  
  ListSet reach_set(initial_set.dimension());

  const VF& vf=vector_field;
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
template<class R, class VF, class BS>
typename Evaluation::IntegratorBase<R,VF,BS>::ListSet
Evaluation::IntegratorBase<R,VF,BS>::lower_integrate(const VectorField& vector_field, 
                                                     const ListSet& initial_set, 
                                                     const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet Integrator::lower_integrate(VectorField vector_field, ListSet initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  using namespace Numeric;
  
  if(time==0) { 
    return initial_set;
  }
  
  const VectorField& vf=vector_field;
  time_type step_size=this->maximum_step_size();
  R maximum_set_radius=this->maximum_basic_set_radius();
  
  if(verbosity>7) { std::clog << "step_size=" << step_size << "  maximum_set_radius=" << maximum_set_radius << std::endl; }
  
  time_type t=0; // t is the time elapsed!
  time_type h=step_size;
  BasicSet bs(initial_set.dimension());
  
  typedef std::pair< time_type, BasicSet> timed_set_type;
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets;
  
  ListSet final_set(initial_set.dimension());
  
  typedef typename ListSet::const_iterator list_set_const_iterator;
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



template<class R, class VF, class BS>
typename Evaluation::IntegratorBase<R,VF,BS>::ListSet
Evaluation::IntegratorBase<R,VF,BS>::lower_reach(const VectorField& vector_field, 
                                                 const ListSet& initial_set, 
                                                 const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet Integrator::lower_reach(VectorField vector_field, ListSet initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  using namespace Numeric;
  
  
  if(time==0) { 
    return initial_set;
  }
  
  const VectorField& vf=vector_field;
  time_type step_size=this->maximum_step_size();
  R maximum_set_radius=this->maximum_basic_set_radius();
  
  if(verbosity>7) { std::clog << "step_size=" << step_size << "  maximum_set_radius=" << maximum_set_radius << std::endl; }
  
  time_type t=0; // t is the time elapsed!
  time_type h=step_size;
  BasicSet bs(initial_set.dimension());
  BasicSet rbs(initial_set.dimension());
  
  typedef std::pair< time_type, BasicSet> timed_set_type;
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets;
  
  ListSet final_set(initial_set.dimension());
  
  typedef typename ListSet::const_iterator list_set_const_iterator;
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
template<class R, class VF, class BS>
typename Evaluation::IntegratorBase<R,VF,BS>::ListSet
Evaluation::IntegratorBase<R,VF,BS>::upper_integrate(const VectorField& vector_field, 
                                                     const ListSet& initial_set, 
                                                     const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet Integrator::upper_integrate(VectorField vector_field, ListSet initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  using namespace Numeric;
  
  if(verbosity>4) { std::clog << "ListSet IntegratorBase::integrate(VectorField,ListSet,Time)" << std::endl; } 
  
  if(time==0) { 
    return initial_set;
  }
  
  const VectorField& vf=vector_field;
  time_type step_size=this->maximum_step_size();
  R maximum_set_radius=this->maximum_basic_set_radius();
  
  if(verbosity>7) { std::clog << "step_size=" << step_size << "  maximum_set_radius=" << maximum_set_radius << std::endl; }
  
  time_type t=0; // t is the time elapsed!
  time_type h=step_size;
  BasicSet bs(initial_set.dimension());
  
  typedef std::pair< time_type, BasicSet> timed_set_type;
  
  // Working sets contains (time,set) pairs, storing the sets reached with different remaining
  std::vector< timed_set_type > working_sets;
  
  ListSet final_set(initial_set.dimension());
  
  typedef typename ListSet::const_iterator list_set_const_iterator;
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
      ListSet subdivisions=bs.subdivide();
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



template<class R, class VF, class BS>
typename Evaluation::IntegratorBase<R,VF,BS>::ListSet
Evaluation::IntegratorBase<R,VF,BS>::upper_reach(const VectorField& vector_field, 
                                                 const ListSet& initial_set, 
                                                 const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet Integrator::integrate(VectorField vector_field, ListSet initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  using namespace Numeric;
  
  const VectorField& vf=vector_field;
  time_type step_size=this->maximum_step_size();
  R maximum_set_radius=this->maximum_basic_set_radius();
  
  if(verbosity>4) {
    std::clog << "step_size=" << conv_approx<double>(step_size) << "  maximum_set_radius()=" << maximum_set_radius << std::endl<<std::flush;
  }
  
  time_type t=0;
  time_type h=step_size;
  BasicSet bs(initial_set.dimension());
  BasicSet rs(initial_set.dimension());
  
  typedef typename ListSet::const_iterator list_set_const_iterator;
  typedef typename ListSet::iterator list_set_iterator;
  typedef std::pair< time_type, BasicSet > timed_set_type;
  
  std::vector< timed_set_type > working_sets;
  ListSet reach_set(initial_set.dimension());
  
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
      if(verbosity>6) { std::clog << "  subdividing..." << std::endl; }
      ListSet subdivisions=bs.subdivide();
      if(verbosity>6) { std::clog << "    subdivisions.size() =" << subdivisions.size() << std::endl; }
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







template<class R, class VF, class BS>
Geometry::ListSet< Geometry::Rectangle<R> >
Evaluation::IntegratorBase<R,VF,BS>::integrate(const System::VectorField<R>& vector_field, 
                                               const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set,
                                               const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet<Rectangle> Integrator::integrate(VectorField vector_field, ListSet<Rectangle> initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  Geometry::ListSet< Geometry::Rectangle<R> > result;
  const VectorField& vf=dynamic_cast<const VectorField&>(vector_field);
  ListSet ils=initial_set;
  ListSet fls=this->lower_integrate(vf,ils,time);
  for(typename ListSet::const_iterator bs_iter=fls.begin();
      bs_iter!=fls.end(); ++bs_iter)
  {
    result.adjoin(bs_iter->bounding_box());
  }
  return result;
}



template<class R, class VF, class BS>
Geometry::ListSet< Geometry::Rectangle<R> >
Evaluation::IntegratorBase<R,VF,BS>::reach(const System::VectorField<R>& vector_field, 
                                           const Geometry::ListSet< Geometry::Rectangle<R> >& initial_set,
                                           const time_type& time) const
{
  ARIADNE_LOG(2,"ListSet<Rectangle> Integrator::integrate(VectorField vector_field, ListSet<Rectangle> initial_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");

  Geometry::ListSet< Geometry::Rectangle<R> > result;
  const VectorField& vf=dynamic_cast<const VectorField&>(vector_field);
  ListSet ils=initial_set;
  ListSet rls=this->lower_reach(vf,ils,time);
  for(typename ListSet::const_iterator bs_iter=rls.begin();
      bs_iter!=rls.end(); ++bs_iter)
  {
    result.adjoin(bs_iter->bounding_box());
  }
  return result;
}



template<class R, class VF, class BS>
Geometry::GridMaskSet<R>
Evaluation::IntegratorBase<R,VF,BS>::integrate(const System::VectorField<R>& vector_field, 
                                               const Geometry::GridMaskSet<R>& initial_set,
                                               const Geometry::GridMaskSet<R>& bounding_set,
                                               const time_type& time) const
{
  ARIADNE_LOG(2,"GridMaskSet Integrator::integrate(VectorField vector_field, GridMaskSet initial_set, GridMaskSet bounding_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  const VectorField& vf=dynamic_cast<const VectorField&>(vector_field);
  
  ARIADNE_CHECK_SAME_GRID(initial_set,bounding_set,"GridMaskSet IntegratorBase::integrate(VectorField,GridMaskSet,GridMaskSet,Time)");
  using namespace System;
  using namespace Geometry;
  using namespace LinearAlgebra;
  
  if(time==0) { 
    return initial_set;
  }
  
  time_type step_size=this->maximum_step_size();
  
  Rectangle bb=bounding_set.bounding_box();
  
  GridMaskSet result(bounding_set);
  result.clear();
  
  Rectangle tmp_rectangle;
  BasicSet tmp_basic_set;
  
  time_type t=time;
  time_type h=step_size;
  
  //bb=regular_intersection(bb,bounding_box());
  
  R spacial_tolerance=2;
  
  ListSet start_set;
  ListSet finish_set;
  for(typename GridMaskSet::const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
    tmp_rectangle=*iter;
    tmp_basic_set=tmp_rectangle;
    start_set.adjoin(tmp_basic_set);
  }
  
  while(t!=0) {
    if(verbosity > 6) { std::clog << "time left=" << t << "  stepsize=" << h << "  sets in list=" << start_set.size() << "\n"; }
    h=min(t,h);
    for(typename ListSet::const_iterator iter=start_set.begin(); iter!=start_set.end(); ++iter) {
      BasicSet p(*iter);
      p=this->integrate(vf,p,h);
      finish_set.adjoin(p);
    }
    start_set.clear();
    GridMaskSet mask_set(bounding_set);
    mask_set.clear();
    for(typename ListSet::const_iterator iter=finish_set.begin(); iter!=finish_set.end(); ++iter) {
      const BasicSet& p=*iter;
      if(p.radius()>spacial_tolerance) {
        if(verbosity > 6) { std::clog << "Splitting, radius=" << p.radius() << "\n" << p << "\n"; }
        mask_set.adjoin_outer_approximation(p);
      }
      else {
        start_set.adjoin(p);
      }
    }
    for(typename GridMaskSet::const_iterator iter=mask_set.begin(); iter!=mask_set.end(); ++iter) {
      tmp_rectangle=*iter;
      tmp_basic_set=tmp_rectangle;
      start_set.adjoin(tmp_basic_set);
    }
    finish_set.clear();
    t-=h;
  }
  
  for(typename ListSet::const_iterator iter=start_set.begin(); iter!=start_set.end(); ++iter) {
    const BasicSet& bs=*iter;
    result.adjoin_outer_approximation(bs);
  }
  return result;
}





template<class R, class VF, class BS>
Geometry::GridMaskSet<R>
Evaluation::IntegratorBase<R,VF,BS>::reach(const System::VectorField<R>& vector_field, 
                                           const Geometry::GridMaskSet<R>& initial_set,
                                           const Geometry::GridMaskSet<R>& bounding_set,
                                           const time_type& time) const
{
  using namespace Numeric;
  
  ARIADNE_LOG(2,"GridMaskSet Integrator::reach(VectorField vector_field, GridMaskSet initial_set, GridMaskSet bounding_set, Time time)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  typedef typename GridMaskSet::const_iterator gms_const_iterator;
  typedef typename ListSet::const_iterator ls_const_iterator;
  ARIADNE_CHECK_BOUNDED(initial_set,"GridMaskSet IntegratorBase::reach(VectorField,GridMaskSet,GridMaskSet,Time)");
  ARIADNE_CHECK_BOUNDED(bounding_set,"GridMaskSet IntegratorBase::reach(VectorField,GridMaskSet,GridMaskSet,Time)");
  const VectorField& vf=dynamic_cast<const VectorField&>(vector_field);
  
  if(!subset(initial_set,bounding_set)) {
    throw std::runtime_error("GridMaskSet IntegratorBase::reach(VectorField,GridMaskSet,GridMaskSet,Time): Initial set must be subset of bounding set");
  }
  
  const GridMaskSet& is=initial_set;
  const Rectangle bb=bounding_set.bounding_box();
  
  GridMaskSet result(bounding_set);
  result.clear();
  
  GridMaskSet stored(result);
  GridMaskSet found(result);
  GridMaskSet image(result);
  found.adjoin(is);
  
  int steps=int_up<int>(time_type(time/this->lock_to_grid_time()));
  if (steps==0) { steps=1; }
  
  time_type time_step=time/steps;
  
  for(int step=0; step!=steps; ++step) {
    found=difference(found,stored);
    stored.adjoin(found);
    image.clear();
    GridMaskSet image=this->integrate(vf,found,bounding_set,time_step);
    found=image;
  }
  
  ListSet input_list;
  ListSet output_list;
  for(gms_const_iterator iter=stored.begin(); iter!=stored.end(); ++iter) {
    input_list.adjoin(BasicSet(Rectangle(*iter)));
  }
  output_list=this->upper_reach(vf,input_list,time_step);
  for(ls_const_iterator iter=output_list.begin(); iter!=output_list.end(); ++iter) {
    const BasicSet& fz=*iter;
    result.adjoin_outer_approximation(fz);
  }
  return result;
}



template<class R, class VF, class BS>
Geometry::GridMaskSet<R>
Evaluation::IntegratorBase<R,VF,BS>::chainreach(const System::VectorField<R>& vector_field, 
                                                const Geometry::GridMaskSet<R>& initial_set, 
                                                const Geometry::GridMaskSet<R>& bounding_set) const
{
  ARIADNE_LOG(2,"GridMaskSet Integrator::chainreach(VectorField vector_field, GridMaskSet initial_set, GridMaskSet bounding_set)\n");
  ARIADNE_LOG(5,"vector_field="<<vector_field<<"\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  typedef typename GridCellListSet::const_iterator gcls_const_iterator;
  typedef typename GridMaskSet::const_iterator gms_const_iterator;
  typedef typename ListSet::const_iterator ls_const_iterator;
  ARIADNE_CHECK_BOUNDED(initial_set,"GridMaskSet IntegratorBase::chainreach(VectorField,GridMaskSet,GridMaskSet)");
  ARIADNE_CHECK_BOUNDED(bounding_set,"GridMaskSet IntegratorBase::chainreach(VectorField,GridMaskSet,GridMaskSet)");
  const VectorField& vf=dynamic_cast<const VectorField&>(vector_field);
  
  if(!subset(initial_set,bounding_set)) {
    throw std::runtime_error("GridMaskSet IntegratorBase::chainreach(VectorField,GridMaskSet,GridMaskSet): Initial set must be subset of bounding set");
  }
  
  const GridMaskSet& is=initial_set;
  const Rectangle bb=bounding_set.bounding_box();
  
  GridMaskSet result(initial_set);
  result.clear();
  GridCellListSet image(initial_set.grid());
  GridCellListSet found(initial_set.grid());
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
    ListSet basic_set_list;
    for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
      ++size;
      Rectangle r=*iter;
      BasicSet z(r);
      basic_set_list.adjoin(z);
    }
    basic_set_list=this->upper_integrate(vf,basic_set_list,time_step);
    if(verbosity>5) { std::clog << "new basic sets " << basic_set_list << std::endl; }
    for(ls_const_iterator iter=basic_set_list.begin(); iter!=basic_set_list.end(); ++iter) {
      const BasicSet& fp=*iter;
      if(!disjoint(fp.bounding_box(),bounding_set)) {
        image.adjoin_outer_approximation(fp);
      } 
    }
    if(verbosity>5) { std::clog << "image set " << image << std::endl; }
    found=regular_intersection(image,bounding_set);
  }
  if(verbosity>5) { std::clog << "Found " << result.size() << " cells, " << std::endl; }
  
  if(verbosity>4) { std::clog << "Beginning reachability phase" << std::endl; }
  ListSet reach_basic_set_list;
  for(gms_const_iterator iter=result.begin(); iter!=result.end(); ++iter) {
    Rectangle r=*iter;
    BasicSet z(r);
    reach_basic_set_list.adjoin(z);
  }
  reach_basic_set_list=this->upper_reach(vf,reach_basic_set_list,time_step);
  for(ls_const_iterator iter=reach_basic_set_list.begin(); iter!=reach_basic_set_list.end(); ++iter) {
    BasicSet fz=*iter;
    if(!disjoint(fz.bounding_box(),bounding_set)) {
      result.adjoin_outer_approximation(fz);
    }
  }
  if(verbosity>4) { std::clog << "Reached " << result.size() << " cells, " << std::endl; }
  
  result.restrict(bounding_set);
  
  return result;
}


template<class R, class VF, class BS>
Geometry::GridMaskSet<R>
Evaluation::IntegratorBase<R,VF,BS>::viable(const System::VectorField<R>& vector_field, 
                                            const Geometry::GridMaskSet<R>& bounding_set) const
{
  ARIADNE_LOG(2,"GridMaskSet Integrator::viable(VectorField vector_field, GridMaskSet bounding_set)\n");
  ARIADNE_LOG(3,"bounding_set="<<bounding_set<<"\n");
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
  ARIADNE_CHECK_BOUNDED(bounding_set,"Integrator<R>::viable(VectorField,GridMaskSet)");
  
  const VectorField& vf=dynamic_cast<const VectorField&>(vector_field);
  const Grid& g=bounding_set.grid();
  Combinatoric::LatticeBlock bd=bounding_set.block();
  GridBlock<R> bb(g,bd);
  GridMaskSet result(g,bd);
  GridCellListSet unsafe=bounding_set;
  
  Rectangle r(g.dimension());
  Rectangle fr(g.dimension());
  GridBlock<R> fgb(g);
  BS bs(g.dimension());
  BS fbs(g.dimension());
  GridCellListSet fgcls(g);
  
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


template<class R, class VF, class BS>
tribool
Evaluation::IntegratorBase<R,VF,BS>::verify(const System::VectorField<R>& vector_field, 
                                            const Geometry::GridMaskSet<R>& initial_set, 
                                            const Geometry::GridMaskSet<R>& safe_set) const
{
  ARIADNE_LOG(2,"GridMaskSet Integrator::verify(VectorField vector_field, GridMaskSet initial_set, GridMaskSet safe_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  ARIADNE_LOG(3,"safe_set="<<safe_set<<"\n");

  typedef typename GridCellListSet::const_iterator gcls_const_iterator;
  typedef typename ListSet::const_iterator ls_const_iterator;
  ARIADNE_CHECK_BOUNDED(initial_set,"tribool IntegratorBase::verify(VectorField,GridMaskSet,GridMaskSet)");
  ARIADNE_CHECK_BOUNDED(safe_set,"tribool IntegratorBase::verify(VectorField,GridMaskSet,GridMaskSet)");
  const VectorField& vf=dynamic_cast<const VectorField&>(vector_field);
  
  if(!subset(initial_set,safe_set)) {
    return false;
  }
  
  const GridMaskSet& is=initial_set;
  const Rectangle bb=safe_set.bounding_box();
  
  GridMaskSet chainreach(is);
  GridCellListSet found(is.grid());
  GridCellListSet image(is.grid());
  GridCellListSet cellimage(is.grid());
  found.adjoin(is);
  
  time_type step_size=this->maximum_step_size();
  time_type time_step=this->lock_to_grid_time();
  
  while(!subset(found,chainreach)) {
    found=difference(found,chainreach);
    chainreach.adjoin(found);
    image.clear();
    uint size=0;
    ListSet basic_set_list;
    for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
      ++size;
      Rectangle r=*iter;
      BasicSet pp(r);
      basic_set_list.adjoin(pp);
    }
    basic_set_list=this->upper_integrate(vf,basic_set_list,time_step);
    for(ls_const_iterator iter=basic_set_list.begin(); iter!=basic_set_list.end(); ++iter) {
      const BasicSet& fp=*iter;
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
