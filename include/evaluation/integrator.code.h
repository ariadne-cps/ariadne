/***************************************************************************
 *            integrator.code.h
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
 
//#define DEBUG

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

#include "../logging.h"

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

#include "integrator.h"

namespace Ariadne {
  namespace Evaluation {
    
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
    
    
    
    
    
    template<class R>
    Integrator<R>::~Integrator()
    {
    }
    
    
    
    template<class R>
    Integrator<R>::Integrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
      : _maximum_step_size(maximum_step_size),
        _lock_to_grid_time(lock_to_grid_time),
        _maximum_basic_set_radius(maximum_basic_set_radius)
    {
    }
    
    
    
    template<class R>
    time_type Integrator<R>::minimum_step_size() const
    {
      return this->_minimum_step_size;
    }

    template<class R>
    time_type Integrator<R>::maximum_step_size() const
    {
      return this->_maximum_step_size;
    }

    template<class R>
    R Integrator<R>::minimum_basic_set_radius() const
    {
      return this->_minimum_basic_set_radius;
    }

    template<class R>
    R Integrator<R>::maximum_basic_set_radius() const
    {
      return this->_maximum_basic_set_radius;
    }

    template<class R>
    R Integrator<R>::grid_size() const
    {
      return this->_grid_size;
    }

    template<class R>
    time_type Integrator<R>::lock_to_grid_time() const
    {
      return this->_lock_to_grid_time;
    }
    

    
    template<class R>
    Geometry::SetInterface<R>*
    Integrator<R>::integrate(const System::VectorField<R>& vector_field,
                             const Geometry::SetInterface<R>& initial_set,
                             const time_type& time) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    Geometry::SetInterface<R>*
    Integrator<R>::integrate(const System::VectorField<R>& vector_field,
                             const Geometry::SetInterface<R>& initial_set,
                             const Geometry::SetInterface<R>& bounding_set,
                             const time_type& time) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    Geometry::SetInterface<R>*
    Integrator<R>::reach(const System::VectorField<R>& vector_field,
                         const Geometry::SetInterface<R>& initial_set,
                         const time_type& time) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    Geometry::SetInterface<R>*
    Integrator<R>::reach(const System::VectorField<R>& vector_field,
                         const Geometry::SetInterface<R>& initial_set,
                         const Geometry::SetInterface<R>& bounding_set,
                         const time_type& time) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    Geometry::SetInterface<R>*
    Integrator<R>::reach(const System::VectorField<R>& vector_field,
                         const Geometry::SetInterface<R>& initial_set) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    Geometry::SetInterface<R>*
    Integrator<R>::reach(const System::VectorField<R>& vector_field,
                         const Geometry::SetInterface<R>& initial_set,
                         const Geometry::SetInterface<R>& bounding_set) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    Geometry::SetInterface<R>*
    Integrator<R>::chainreach(const System::VectorField<R>& vector_field,
                              const Geometry::SetInterface<R>& initial_set,
                              const Geometry::SetInterface<R>& bounding_set) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }

    template<class R>
    tribool
    Integrator<R>::verify(const System::VectorField<R>& vector_field,
                          const Geometry::SetInterface<R>& initial_set,
                          const Geometry::SetInterface<R>& bounding_set) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }



    
    
    template<class R>
    bool
    Integrator<R>::check_flow_bounds(const System::VectorField<R>& vf,
                                     const Geometry::Rectangle<R>& r,
                                     const Geometry::Rectangle<R>& b,
                                     const time_type& h) const
    {
      if(verbosity>6) { std::clog << __FUNCTION__ << std::endl; }
      using namespace Geometry;
      using namespace Numeric;
      return subset(r+Interval<R>(0,h)*vf(b),b);
    }
    
    
    
    template<class R>
    Geometry::Rectangle<R>
    Integrator<R>::estimate_flow_bounds(const System::VectorField<R>& vf,
                                        const Geometry::Rectangle<R>& r,
                                        const time_type& h,
                                        const unsigned int& maximum_iterations) const
    {
      using namespace Geometry;
      using namespace Numeric;

      if(verbosity>6) { std::clog << __FUNCTION__ << " (maximum_iterations=" << maximum_iterations << ")" << std::endl; }
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
    Integrator<R>::estimate_flow_bounds(const System::VectorField<R>& vf,
                                        const Geometry::Rectangle<R>& r,
                                        time_type& h) const
    {
      using namespace Geometry;
      using namespace Numeric;
      
      if(verbosity>6) { std::clog << __FUNCTION__ << std::endl; }

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
    Integrator<R>::refine_flow_bounds(const System::VectorField<R>& vector_field,
                                      const Geometry::Rectangle<R>& initial_set,
                                      const Geometry::Rectangle<R>& estimated_bounds,
                                      const time_type& step_size) const
    {
      if(verbosity>6) { std::clog << __FUNCTION__ << std::endl; }

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
    typename IntegratorBase<R,VF,BS>::BasicSet 
    IntegratorBase<R,VF,BS>::integrate(const VectorField& vector_field, 
                                                const BasicSet& initial_set, 
                                                const time_type& time) const
    {
      if(verbosity>4) { std::clog << "IntegratorBase::integrate(VectorField,BasicSet,time_type)" << std::endl; }
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
    
    
    
    // Template pattern for integrating a list set
    template<class R, class VF, class BS>
    typename IntegratorBase<R,VF,BS>::ListSet
    IntegratorBase<R,VF,BS>::integrate(const VectorField& vector_field, 
                                      const ListSet& initial_set, 
                                      const time_type& time) const
    {
      using namespace Numeric;
      
      if(verbosity>4) { std::clog << "IntegratorBase::integrate(VectorField,ListSet,time_type)" << std::endl; }

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
    typename IntegratorBase<R,VF,BS>::ListSet
    IntegratorBase<R,VF,BS>::reach(const VectorField& vector_field, 
                                   const ListSet& initial_set, 
                                   const time_type& time) const
    {
      using namespace Numeric;
      
      if(verbosity>4) { std::clog << "IntegratorBase::reach(VectorField,ListSet,time_type)" << std::endl; }
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
    typename IntegratorBase<R,VF,BS>::ListSet
    IntegratorBase<R,VF,BS>::reach(const VectorField& vector_field, 
                                   const ListSet& initial_set) const
    {
      using namespace Numeric;
      
      if(verbosity>4) { std::clog << "IntegratorBase::reach(VectorField,ListSet)" << std::endl; }
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

      ListSet reach_set(initial_set.dimension());
      
      if(verbosity>6) { std::clog << "initial_set.size()=" << initial_set.size() << std::endl; }
      assert(reach_set.size()==0);

      for(list_set_const_iterator bs_iter=initial_set.begin(); bs_iter!=initial_set.end(); ++bs_iter) {
        bs=*bs_iter;
        h=step_size;
        
        if(verbosity>6) { std::clog << "  t=" << conv_approx<double>(t) << "  bs.centre()=" << bs.centre() << "  bs.radius() = " << bs.radius() << std::endl; }

        while(bs.radius()<maximum_set_radius) {
          if(verbosity>6) { std::clog << "  integrating..." << std::endl; }
         
          if(verbosity>6) {
            std::clog << "  h=" << conv_approx<double>(h)
                      << "  bs=" << bs << std::endl;
          }

          rs=this->reachability_step(vf,bs,h);
          reach_set.adjoin(rs);
          
          bs=integration_step(vf,bs,h);
        } 
          
      }
      
      if(verbosity>4) {
        std::clog << "  reach_set.size()=" << reach_set.size() <<  std::endl;
      }
      return reach_set;
    } 

    
    
    
    
    
    


    template<class R, class VF, class BS>
    Geometry::GridMaskSet<R>
    IntegratorBase<R,VF,BS>::integrate(const System::VectorField<R>& vector_field, 
                                      const Geometry::GridMaskSet<R>& initial_set,
                                      const Geometry::GridMaskSet<R>& bounding_set,
                                      const time_type& time) const
    {
      if(verbosity>4) { std::clog << "IntegratorBase::integrate(VectorField,GridMaskSet,GridMaskSet,time_type)" << std::endl; }

      const VectorField& vf=dynamic_cast<const VectorField&>(vector_field);

      check_same_grid(initial_set,bounding_set,"IntegratorBase::integrate(VectorField,GridMaskSet,GridMaskSet,time_type)");
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

        if(verbosity>4) { std::clog << "time left=" << t << "  stepsize=" << h << "  sets in list=" << start_set.size() << "\n"; }

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

            if(verbosity>4) { std::clog << "Splitting, radius=" << p.radius() << "\n" << p << "\n"; }

            mask_set.adjoin_over_approximation(p);
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
        result.adjoin_over_approximation(bs);
      }
      return result;
    }
   
    
    
  
    
    template<class R, class VF, class BS>
    Geometry::GridMaskSet<R>
    IntegratorBase<R,VF,BS>::reach(const System::VectorField<R>& vector_field, 
                                  const Geometry::GridMaskSet<R>& initial_set,
                                  const Geometry::GridMaskSet<R>& bounding_set,
                                  const time_type& time) const
    {
      using namespace Numeric;
      
      if(verbosity>4) { std::clog << "IntegratorBase::reach(VectorField,GridMaskSet,GridMaskSet,time_type)" << std::endl; }

      typedef typename GridMaskSet::const_iterator gms_const_iterator;
      typedef typename ListSet::const_iterator ls_const_iterator;
      check_bounded(initial_set,"IntegratorBase::reach(VectorField,GridMaskSet,GridMaskSet,time_type)");
      check_bounded(bounding_set,"IntegratorBase::reach(VectorField,GridMaskSet,GridMaskSet,time_type)");
      const VectorField& vf=dynamic_cast<const VectorField&>(vector_field);
      
      if(!subset(initial_set,bounding_set)) {
        throw std::runtime_error("IntegratorBase::reach(VectorField,GridMaskSet,GridMaskSet,time_type): Initial set must be subset of bounding set");
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
      output_list=this->reach(vf,input_list,time_step);
      for(ls_const_iterator iter=output_list.begin(); iter!=output_list.end(); ++iter) {
        const BasicSet& fz=*iter;
        result.adjoin_over_approximation(fz);
      }
      return result;
    }
    
    
    
    template<class R, class VF, class BS>
    Geometry::GridMaskSet<R>
    IntegratorBase<R,VF,BS>::chainreach(const System::VectorField<R>& vector_field, 
                                        const Geometry::GridMaskSet<R>& initial_set, 
                                        const Geometry::GridMaskSet<R>& bounding_set) const
    {
      if(verbosity>4) { std::clog << "IntegratorBase::chainreach(VectorField,GridMaskSet,GridMaskSet)"; }
      typedef typename GridCellListSet::const_iterator gcls_const_iterator;
      typedef typename GridMaskSet::const_iterator gms_const_iterator;
      typedef typename ListSet::const_iterator ls_const_iterator;
      check_bounded(initial_set,"IntegratorBase::chainreach(VectorField,GridMaskSet,GridMaskSet)");
      check_bounded(bounding_set,"IntegratorBase::chainreach(VectorField,GridMaskSet,GridMaskSet)");
      const VectorField& vf=dynamic_cast<const VectorField&>(vector_field);
     
      if(!subset(initial_set,bounding_set)) {
        throw std::runtime_error("IntegratorBase::chainreach(VectorField,GridMaskSet,GridMaskSet): Initial set must be subset of bounding set");
      }

      if(verbosity>3) { std::clog << "initial_set=" << initial_set << std::endl; }
      if(verbosity>3) { std::clog << "bounding_set=" << bounding_set << std::endl; }

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
        basic_set_list=this->integrate(vf,basic_set_list,time_step);
        for(ls_const_iterator iter=basic_set_list.begin(); iter!=basic_set_list.end(); ++iter) {
          const BasicSet& fp=*iter;
          if(!disjoint(fp.bounding_box(),bounding_set)) {
            image.adjoin_over_approximation(fp);
          }
        }
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
      reach_basic_set_list=this->reach(vf,reach_basic_set_list,time_step);
      for(ls_const_iterator iter=reach_basic_set_list.begin(); iter!=reach_basic_set_list.end(); ++iter) {
        BasicSet fz=*iter;
        if(!disjoint(fz.bounding_box(),bounding_set)) {
          result.adjoin_over_approximation(fz);
        }
      }
      result=regular_intersection(result,bounding_set);
      if(verbosity>4) { std::clog << "Reached " << result.size() << " cells, " << std::endl; }

      if(!subset(result,bounding_set)) {
        std::cerr << "WARNING: result is not a subset of bounding_set\n"; 
      }
      return result;
    }

    template<class R, class VF, class BS>
    tribool
    IntegratorBase<R,VF,BS>::verify(const System::VectorField<R>& vector_field, 
                                   const Geometry::GridMaskSet<R>& initial_set, 
                                   const Geometry::GridMaskSet<R>& safe_set) const
    {
      if(verbosity>4) { std::clog << "IntegratorBase::verify(VectorField,GridMaskSet,GridMaskSet)"; }
      typedef typename GridCellListSet::const_iterator gcls_const_iterator;
      typedef typename ListSet::const_iterator ls_const_iterator;
      check_bounded(initial_set,"IntegratorBase::verify(VectorField,GridMaskSet,GridMaskSet)");
      check_bounded(safe_set,"IntegratorBase::verify(VectorField,GridMaskSet,GridMaskSet)");
      const VectorField& vf=dynamic_cast<const VectorField&>(vector_field);
     
      if(!subset(initial_set,safe_set)) {
        throw std::runtime_error("IntegratorBase::verify(VectorField,GridMaskSet,GridMaskSet): Initial set must be subset of bounding set");
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
        basic_set_list=this->integrate(vf,basic_set_list,time_step);
        for(ls_const_iterator iter=basic_set_list.begin(); iter!=basic_set_list.end(); ++iter) {
          const BasicSet& fp=*iter;
          cellimage.adjoin_over_approximation(fp);
          if(!subset(cellimage,safe_set)) {
            return false;
          }
          image.adjoin(cellimage);
          cellimage.clear();
        }
        found=image;
      }
      return true;
    }

  }
}
