/***************************************************************************
 *            integrator.tpl
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

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../declarations.h"

#include "../base/array.h"
#include "../base/exceptions.h"

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
    time_type Integrator<R>::lock_to_grid_time() const
    {
      return this->_lock_to_grid_time;
    }
    
    
    
    
    template<class R>
    LinearAlgebra::Matrix<R>
    Integrator<R>::symmetrize(const LinearAlgebra::Vector< Interval<R> >& iv)
    {
      LinearAlgebra::Matrix<R> A(iv.size(),iv.size()+1);
      for(size_type i=0; i!=A.number_of_rows(); ++i) {
        A(i,i)=iv(i).radius();
        A(i,iv.size())=iv(i).centre();
      }
      return A;
    }
    
    
    
    template<class R>
    bool
    Integrator<R>::check_flow_bounds(const System::VectorField<R>& vf,
                                     const Geometry::Rectangle<R>& r,
                                     const Geometry::Rectangle<R>& b,
                                     const time_type& h) const
    {
      using namespace Geometry;
      return subset(r+Interval<R>(0,h)*vf(b),b);
    }
    
    
    
    template<class R>
    Geometry::Rectangle<R>
    Integrator<R>::estimate_flow_bounds(const System::VectorField<R>& vf,
                                        const Geometry::Rectangle<R>& r,
                                        const time_type& h,
                                        const unsigned int& maximum_iterations) const
    {
#ifdef DEBUG
      std::cerr << "estimate_flow_bounds(VectorField<R>, Rectangle<R>, R, uint)\n";
      std::cerr << "r=" << r << "  h=" << h << "  max_iter=" << maximum_iterations << std::endl;
#endif

      using namespace Geometry;
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
        
#ifdef DEBUG
        std::cerr << "t=" << convert_to<double>(t) << "  reach="  << reach << std::endl;
#endif
        ++iteration;
        if(iteration==maximum_iterations) {
          throw std::runtime_error("Cannot find bounding box for flow");
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
#ifdef DEBUG
      std::cerr << "estimate_flow_bounds" << std::endl;
#endif

      unsigned int max_iterations=16;
      
      Geometry::Rectangle<R> bounds(vf.dimension());
      while(bounds.empty()) {
        try {
          bounds=estimate_flow_bounds(vf,r,h,max_iterations);
        }
        catch(std::runtime_error) { 
          h/=2;
          max_iterations*=2;
        }
      }
#ifdef DEBUG
      std::cerr << "bounds=" << bounds << std::endl;
#endif
      return bounds;
    }
    
    
    
    template<class R>
    Geometry::Rectangle<R>
    Integrator<R>::refine_flow_bounds(const System::VectorField<R>& vector_field,
                                      const Geometry::Rectangle<R>& initial_set,
                                      const Geometry::Rectangle<R>& estimated_bounds,
                                      const time_type& step_size) const
    {
      using namespace System;
      using namespace Geometry;
      using namespace LinearAlgebra;
      const VectorField<R>& vf=vector_field;
      Rectangle<R> rx=initial_set;
      Rectangle<R> b=estimated_bounds;
      Interval<R> h=step_size;
      
      Rectangle<R> xb=rx+Interval<R>(0,step_size)*vf(b);
      Rectangle<R> xxb=rx+Interval<R>(0,step_size)*vf(xb);
#ifdef DEBUG
      std::cerr << "new bounds " << xxb << "," << xb << " vs old bounds " << b << "  " << subset(xb,b) << std::endl;
#endif
      Vector< Interval<R> > ddphi=vf.jacobian(xb)*vf(xb);
      Vector< Interval<R> > dfx=vf(rx);
      Vector< Interval<R> > hdfx=h*dfx;
      Vector< Interval<R> > hhddphi=(h*h/R(2))*ddphi;
      Vector< Interval<R> > dx=hdfx+hhddphi;
      return rx+dx;
    }
    
    
    
    template<class R>
    template<template<class> class BS>
    BS<R> 
    Integrator<R>::integrate_basic_set(const System::VectorField<R>& vector_field, 
                                       const BS<R>& initial_set, 
                                       const time_type& time) const
    {
      if(time==0) { 
        return initial_set;
      }
      
      const System::VectorField<R>& vf=vector_field;
      BS<R> r(initial_set);
      time_type t=time;
      time_type h=this->maximum_step_size();
      while(t>0) {
        h=min(t,h);
        r=this->integration_step(vf,r,h);
        t=t-h;
        h=min(time_type(2*h),this->maximum_step_size());
      }
      return r;
    }
    
    
    
    template<class R>
    template<template<class> class BS>
    Geometry::ListSet<R,BS> 
    Integrator<R>::integrate_list_set(const System::VectorField<R>& vector_field, 
                                      const Geometry::ListSet<R,BS>& initial_set, 
                                      const time_type& time) const
    {
      if(time==0) { 
        return initial_set;
      }

      const System::VectorField<R>& vf=vector_field;
      time_type step_size=this->maximum_step_size();
      R maximum_set_radius=this->maximum_basic_set_radius();
#ifdef DEBUG
        std::cerr << "step_size=" << step_size << "  maximum_set_radius=" << maximum_set_radius << std::endl;
#endif
      
      time_type t=time;
      time_type h=step_size;
      BS<R> bs(initial_set.dimension());
      
      std::multiset< std::pair<R,BS<R> >, pair_first_less<R,BS<R> > > working_sets;
      Geometry::ListSet<R,BS> final_set(initial_set.dimension());
      
      typedef typename Geometry::ListSet<R,BS>::const_iterator list_set_const_iterator;
      typedef std::pair< time_type,BS<R> > timed_set;
      for(list_set_const_iterator bs_iter=initial_set.begin(); bs_iter!=initial_set.end(); ++bs_iter) {
        working_sets.insert(timed_set(time,*bs_iter));
      }
      
      while(!working_sets.empty()) {
#ifdef DEBUG
        //std::cerr << working_sets << "\n\n\n";
#endif
        timed_set ts=*working_sets.begin();
        working_sets.erase(working_sets.begin());
        t=ts.first;
        bs=ts.second;
        h=step_size;
        

        if(bs.radius()>maximum_set_radius) {
          Geometry::ListSet<R,BS> subdivisions=bs.subdivide();
          for(list_set_const_iterator subdiv_iter=subdivisions.begin(); 
              subdiv_iter!=subdivisions.end(); ++subdiv_iter)
          {
            working_sets.insert(timed_set(t,*subdiv_iter));
          }
        }
        else {
          do {
#ifdef DEBUG
        std::cerr << "time left=" << t << "  stepsize=" << h << "  centre=" << bs.centre() 
                  << "  radius=" << bs.radius() << std::endl;
#endif
            h=min(t,h);
            bs=this->integration_step(vf,bs,h);
            t=t-h;
            h=min(time_type(2*h),step_size);
          } while(t!=0 && bs.radius()<=maximum_set_radius);
          if(t==0) {
            final_set.adjoin(bs);
          } else {
            working_sets.insert(timed_set(t,bs));
          }
        }
      }
      return final_set;
    }
    
    
    
    
    
    
    
    template<class R>
    Geometry::ListSet<R,Geometry::Rectangle>
    Integrator<R>::reach(const System::VectorField<R>& vector_field,
                         const Geometry::ListSet<R,Geometry::Rectangle>& initial_set,
                         const time_type& time) const
    {
      throw std::domain_error("Integrator<R>::reach(const System::VectorField<R>& vector_field,"
                              "const Geometry::ListSet<R,Geometry::Rectangle>& initial_set,"
                              "const time_type& time) const not implemented");
    }
    
    
    
    template<class R>
    Geometry::Rectangle<R>
    Integrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Rectangle<R>& initial_set, 
                                      time_type& time) const
    {
      throw NotImplemented("Integrator<R>::integration_step(VectorField<R>,Rectangle<R>,T)");
    }
    
    template<class R>
    Geometry::Parallelotope<R>
    Integrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Parallelotope<R>& initial_set, 
                                      time_type& time) const
    {
      return Geometry::Parallelotope<R>(this->integration_step(vector_field,initial_set.bounding_box(),time));
    }
    
    template<class R>
    Geometry::Zonotope<R>
    Integrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Zonotope<R>& initial_set, 
                                      time_type& time) const
    {
      return Geometry::Zonotope<R>(this->integration_step(vector_field,initial_set.bounding_box(),time));
    }
    
    template<class R>
    Geometry::Parallelotope< Interval<R> >
    Integrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Parallelotope< Interval<R> >& initial_set, 
                                      time_type& time) const
    {
      return Geometry::Parallelotope< Interval<R> >(this->integration_step(vector_field,initial_set.bounding_box(),time));
    }
    
    template<class R>
    Geometry::Zonotope< Interval<R> >
    Integrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                      const Geometry::Zonotope< Interval<R> >& initial_set, 
                                      time_type& time) const
    {
      return Geometry::Zonotope< Interval<R> >(this->integration_step(vector_field,initial_set.bounding_box(),time));
    }
    
    
    
    template<class R>
    Geometry::Rectangle<R> 
    Integrator<R>::reachability_step(const System::VectorField<R>& vf,
                                       const Geometry::Rectangle<R>& is,
                                       time_type& h) const 
    {
      throw NotImplemented("Integrator<R>::reachability_step(VectorField<R>,Rectangle<R>,time_type)");
    }
    
    
    template<class R>
    Geometry::Zonotope<R> 
    Integrator<R>::reachability_step(const System::VectorField<R>& vf,
                                       const Geometry::Zonotope<R>& is,
                                       time_type& h) const 
    {
      return this->reachability_step(vf,is.bounding_box(),h);
    }

    template<class R>
    Geometry::Zonotope< Interval<R> > 
    Integrator<R>::reachability_step(const System::VectorField<R>& vf,
                                       const Geometry::Zonotope< Interval<R> >& is,
                                       time_type& h) const 
    {
      return Geometry::Zonotope< Interval<R> >(this->reachability_step(vf,is.bounding_box(),h));
    }




    template<class R>
    Geometry::Rectangle<R>
    Integrator<R>::integrate(const System::VectorField<R>& vector_field, 
                               const Geometry::Rectangle<R>& initial_set, 
                               const time_type& time) const
    {
      return this->integrate_basic_set(vector_field,initial_set,time);
    }
    
    
    template<class R>
    Geometry::Parallelotope<R>
    Integrator<R>::integrate(const System::VectorField<R>& vector_field, 
                               const Geometry::Parallelotope<R>& initial_set, 
                               const time_type& time) const
    {
      return this->integrate_basic_set(vector_field,initial_set,time);
    }
    
    
    template<class R>
    Geometry::Zonotope<R>
    Integrator<R>::integrate(const System::VectorField<R>& vector_field, 
                               const Geometry::Zonotope<R>& initial_set, 
                               const time_type& time) const
    {
      return this->integrate_basic_set(vector_field,initial_set,time);
    }
    
    
    
    template<class R>
    Geometry::ListSet<R,Geometry::Rectangle>
    Integrator<R>::integrate(const System::VectorField<R>& vector_field,
                               const Geometry::ListSet<R,Geometry::Rectangle>& initial_set,
                               const time_type& time) const
    {
     return this->integrate_list_set(vector_field,initial_set,time);
    }
    
    
    template<class R>
    Geometry::ListSet<R,Geometry::Parallelotope> 
    Integrator<R>::integrate(const System::VectorField<R>& vector_field, 
                               const Geometry::ListSet<R,Geometry::Parallelotope>& initial_set, 
                               const time_type& time) const
    {
      return this->integrate_list_set(vector_field,initial_set,time);
    }
    
    
    template<class R>
    Geometry::ListSet<R,Geometry::Zonotope> 
    Integrator<R>::integrate(const System::VectorField<R>& vector_field, 
                               const Geometry::ListSet<R,Geometry::Zonotope>& initial_set, 
                               const time_type& time) const
    {
      return this->integrate_list_set(vector_field,initial_set,time);
    }
   
    
    template<class R>
    Geometry::GridMaskSet<R> 
    Integrator<R>::integrate(const System::VectorField<R>& vector_field, 
                             const Geometry::GridMaskSet<R>& initial_set,
                             const Geometry::GridMaskSet<R>& bounding_set,
                             const time_type& time) const
    {
      assert(initial_set.grid()==bounding_set.grid());
      using namespace System;
      using namespace Geometry;
      using namespace LinearAlgebra;
      
      if(time==0) { 
        return initial_set;
      }
      
      time_type step_size=this->maximum_step_size();
      
      const Geometry::Grid<R>& g(initial_set.grid());
      Combinatoric::LatticeBlock lb=bounding_set.block();
      Geometry::Rectangle<R> bb=bounding_set.bounding_box();
      
      Geometry::GridMaskSet<R> result(bounding_set.grid(),bounding_set.block());
      
      time_type t=time;
      time_type h=step_size;
      
      //bb=regular_intersection(bb,bounding_box());
      
      R spacial_tolerance=2;
      
      ListSet<R,Parallelotope> start_set;
      ListSet<R,Parallelotope> finish_set;
      for(typename GridMaskSet<R>::const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
        start_set.adjoin(Parallelotope<R>(*iter));
      }
      
      while(t!=0) {
#ifdef DEBUG
        std::cerr << "time left=" << t << "  stepsize=" << h << "  sets in list=" << start_set.size() << "\n";
#endif
        h=min(t,h);
        for(typename ListSet<R,Parallelotope>::const_iterator iter=start_set.begin(); iter!=start_set.end(); ++iter) {
          const VectorField<R>& vf=vector_field;
          Geometry::Parallelotope<R> p(*iter);
          p=this->integrate(vf,p,h);
          finish_set.adjoin(p);
        }
        start_set.clear();
        GridMaskSet<R> mask_set(g,lb);
        for(typename ListSet<R,Parallelotope>::const_iterator iter=finish_set.begin(); iter!=finish_set.end(); ++iter) {
          const Parallelotope<R>& p=*iter;
          if(p.radius()>spacial_tolerance) {
#ifdef DEBUG
            std::cerr << "Splitting, radius=" << p.radius() << "\n" << p << "\n";
#endif
            GridCellListSet<R> p_approx=over_approximation(p,g);
            mask_set.adjoin(p_approx);
          }
          else {
            start_set.adjoin(p);
          }
        }
        for(typename GridMaskSet<R>::const_iterator iter=mask_set.begin(); iter!=mask_set.end(); ++iter) {
          start_set.adjoin(Geometry::Parallelotope<R>(*iter));
        }
        finish_set.clear();
        t-=h;
      }
    
      for(typename ListSet<R,Parallelotope>::const_iterator iter=start_set.begin(); iter!=start_set.end(); ++iter) {
        GridCellListSet<R> oai=over_approximation(*iter,g);
        result.adjoin(oai);
      }
      return result;
    }
   
        
    template<class R>
    Geometry::ListSet<R,Geometry::Zonotope> 
    Integrator<R>::reach(const System::VectorField<R>& vector_field, 
                         const Geometry::ListSet<R,Geometry::Parallelotope>& initial_set, 
                         const time_type& time) const
    {
      const System::VectorField<R>& vf=vector_field;
      time_type step_size=this->maximum_step_size();
      R maximum_set_radius=this->maximum_basic_set_radius();
#ifdef DEBUG
        std::cerr << "step_size=" << step_size << "  maximum_set_radius()=" << maximum_set_radius << std::endl;
#endif
      
      time_type t=time;
      time_type h=step_size;
      Geometry::Parallelotope<R> bs(initial_set.dimension());
      Geometry::Zonotope<R> rs(initial_set.dimension());
      
      std::multiset< std::pair< R,Geometry::Parallelotope<R> >, pair_first_less< R,Geometry::Parallelotope<R> > > working_sets;
      Geometry::ListSet< R,Geometry::Zonotope > reach_set(initial_set.dimension());
      
      typedef typename Geometry::ListSet<R,Geometry::Parallelotope>::const_iterator list_set_const_iterator;
      typedef std::pair< time_type,Geometry::Parallelotope<R> > timed_set;
      for(list_set_const_iterator bs_iter=initial_set.begin(); bs_iter!=initial_set.end(); ++bs_iter) {
        working_sets.insert(timed_set(time,*bs_iter));
      }
      
      while(!working_sets.empty()) {
        
        timed_set ts=*working_sets.begin();
        working_sets.erase(working_sets.begin());
        t=ts.first;
        bs=ts.second;
        h=step_size;

        if(bs.radius()>maximum_set_radius) {
          Geometry::ListSet<R,Geometry::Parallelotope> subdivisions=bs.subdivide();
          for(list_set_const_iterator subdiv_iter=subdivisions.begin(); 
              subdiv_iter!=subdivisions.end(); ++subdiv_iter)
          {
            working_sets.insert(timed_set(t,*subdiv_iter));
          }
        }
        else {
          do {
#ifdef DEBUG
        std::cerr << "time left=" << t << "  stepsize=" << h << "  centre=" 
                  << bs.centre()  << std::endl;
#endif
            h=min(t,h);
            rs=Geometry::over_approximation(Geometry::Zonotope< Interval<R> >(this->reachability_step(vf,bs,h)));
            reach_set.adjoin(rs);
#ifdef DEBUG
        std::cerr << "rs.centre=" << rs.centre() << std::endl;
#endif
            
            bs=integration_step(vf,bs,h);
            t=t-h;
            h=min(time_type(2*h),step_size);
          } while(t!=0 && bs.radius()<=maximum_set_radius);
          if(t!=0) {
            working_sets.insert(timed_set(t,bs));
          }
        }
      }
      
      return reach_set;
    }


    template<class R>
    Geometry::ListSet<R,Geometry::Zonotope> 
    Integrator<R>::reach(const System::VectorField<R>& vector_field, 
                           const Geometry::ListSet<R,Geometry::Zonotope>& initial_set, 
                           const time_type& time) const
    {
      const System::VectorField<R>& vf=vector_field;
      time_type step_size=this->maximum_step_size();
      R maximum_set_radius=this->maximum_basic_set_radius();
#ifdef DEBUG
        std::cerr << "Integrator<R>::reach(const VectorField<R>& vector_field,const Geometry::ListSet<R,Geometry::Zonotope>& initial_set, const time_type& time)" << std::endl;
        std::cerr << "step_size=" << step_size << "  maximum_set_radius()=" << maximum_set_radius << std::endl<<std::flush;
#endif
      
      time_type t=time;
      time_type h=step_size;
      Geometry::Zonotope<R> bs(initial_set.dimension(),initial_set.dimension());
      Geometry::Zonotope<R> rs(initial_set.dimension(),initial_set.dimension()+1);
      
      std::multiset< std::pair< R,Geometry::Zonotope<R> >, pair_first_less< R,Geometry::Zonotope<R> > > working_sets;
      Geometry::ListSet<R,Geometry::Zonotope> reach_set(initial_set.dimension());
      
      typedef typename Geometry::ListSet<R,Geometry::Zonotope>::const_iterator list_set_const_iterator;
      typedef std::pair< time_type,Geometry::Zonotope<R> > timed_set;
      for(list_set_const_iterator bs_iter=initial_set.begin(); bs_iter!=initial_set.end(); ++bs_iter) {
        working_sets.insert(timed_set(time,*bs_iter));
      }
      
      while(!working_sets.empty()) {
#ifdef DEBUG
        //std::cerr << working_sets << "\n\n\n";
#endif
        timed_set ts=*working_sets.begin();
        working_sets.erase(working_sets.begin());
        t=ts.first;
        bs=ts.second;
        h=step_size;
        

        if(bs.radius()>maximum_set_radius) {
          Geometry::ListSet<R,Geometry::Zonotope> subdivisions=bs.subdivide();
          for(list_set_const_iterator subdiv_iter=subdivisions.begin(); 
              subdiv_iter!=subdivisions.end(); ++subdiv_iter)
          {
            working_sets.insert(timed_set(t,*subdiv_iter));
          }
        }
        else {
          do {
#ifdef DEBUG
        std::cerr << "time left=" << t << "  stepsize=" << h << "  centre=" 
	          << bs.centre() << std::endl;
#endif
            h=min(t,h);
            rs=this->reachability_step(vf,bs,h);
            reach_set.adjoin(rs);
#ifdef DEBUG
        std::cerr << "rs.centre=" << rs.centre() << std::endl;
#endif
            
            bs=integration_step(vf,bs,h);
            t=t-h;
            h=min(time_type(2*h),step_size);
          } while(t!=0 && bs.radius()<=maximum_set_radius);
          if(t!=0) {
            working_sets.insert(timed_set(t,bs));
          }
        }
      }
      
#ifdef DEBUG
        std::cerr << "Exiting reach"<< std::endl;
#endif
      return reach_set;
    } 

    template<class R>
    Geometry::GridMaskSet<R> 
    Integrator<R>::reach(const System::VectorField<R>& vector_field, 
                           const Geometry::GridMaskSet<R>& initial_set,
                           const Geometry::GridMaskSet<R>& bounding_set,
                           const time_type& time) const
    {
      typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
      typedef typename Geometry::ListSet<R,Geometry::Zonotope>::const_iterator zls_const_iterator;
      assert(initial_set.bounded() && bounding_set.bounded());
      
      if(!subset(initial_set,bounding_set)) {
        throw std::runtime_error("chainreach: Initial set must be subset of bounding set");
      }
        
      const Geometry::Grid<R>& g=initial_set.grid();
      const Geometry::GridMaskSet<R>& is=initial_set;
      const Geometry::Rectangle<R>& bb=bounding_set.bounding_box();
      const Combinatoric::LatticeBlock lb=over_approximation(bb,g).lattice_set();
      
      Geometry::GridMaskSet<R> stored(g,lb);
      Geometry::GridMaskSet<R> found(g,lb);
      Geometry::GridMaskSet<R> image(g,lb);
      Geometry::GridMaskSet<R> result(g,lb);
      found.adjoin(is);
      
      int steps=ceil<int>(time_type(time/this->lock_to_grid_time()));
      if (steps==0) { steps=1; }

      time_type time_step=time/steps;
     
      for(int step=0; step!=steps; ++step) {
        found=difference(found,stored);
        stored.adjoin(found);
        image.clear();
        Geometry::GridMaskSet<R> image=this->integrate(vector_field,found,bounding_set,time_step);
        found=image;
      }
      
      Geometry::ListSet<R,Geometry::Zonotope> input_list;
      Geometry::ListSet<R,Geometry::Zonotope> output_list;
      for(gms_const_iterator iter=stored.begin(); iter!=stored.end(); ++iter) {
        //Geometry::Rectangle<R> r=*iter;
        //Geometry::Zonotope<R> z(r);
        input_list.adjoin(Geometry::Zonotope<R>(*iter));
      }
      output_list=this->reach(vector_field,input_list,time_step);
      for(zls_const_iterator iter=output_list.begin(); iter!=output_list.end(); ++iter) {
        Geometry::Zonotope<R> fz=*iter;
        Geometry::GridCellListSet<R> oai=over_approximation(fz,g);
        result.adjoin(oai);
      }
      return result;
    }
    
    
    
    template<class R>
    Geometry::GridMaskSet<R> 
    Integrator<R>::chainreach(const System::VectorField<R>& vf, 
                              const Geometry::GridMaskSet<R>& initial_set, 
                              const Geometry::GridMaskSet<R>& bounding_set) const
    {
      typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
      typedef typename Geometry::ListSet<R,Geometry::Parallelotope>::const_iterator pls_const_iterator;
      assert(initial_set.bounded() && bounding_set.bounded());
     
      if(!subset(initial_set,bounding_set)) {
        throw std::runtime_error("chainreach: Initial set must be subset of bounding set");
      }
        
      const Geometry::Grid<R>& g=initial_set.grid();
      const Geometry::GridMaskSet<R>& is=initial_set;
      const Geometry::Rectangle<R>& bb=bounding_set.bounding_box();
      const Combinatoric::LatticeBlock lb=over_approximation(bb,g).lattice_set();
      
      Geometry::GridMaskSet<R> result(g,lb);
      Geometry::GridCellListSet<R> found(g);
      Geometry::GridCellListSet<R> image(g);
      found.adjoin(is);
      
      time_type step_size=this->maximum_step_size();
      time_type time_step=this->lock_to_grid_time();
      
      while(!subset(found,result)) {
        found=difference(found,result);
        result.adjoin(found);
        image.clear();
        uint size=0;
        Geometry::ListSet<R,Geometry::Parallelotope> parallotope_list;
        for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
          ++size;
          Geometry::Rectangle<R> r=*iter;
          Geometry::Parallelotope<R> pp(r);
          parallotope_list.adjoin(pp);
        }
        parallotope_list=this->integrate_list_set(vf,parallotope_list,time_step);
        for(pls_const_iterator iter=parallotope_list.begin(); iter!=parallotope_list.end(); ++iter) {
          Geometry::Parallelotope<R> fp=*iter;
          if(!disjoint(fp.bounding_box(),bounding_set)) {
            image.adjoin(over_approximation(fp,g));
          }
        }
        found=regular_intersection(image,bounding_set);
      }
      return result;
    }

    template<class R>
    bool
    Integrator<R>::verify(const System::VectorField<R>& vf, 
                          const Geometry::GridMaskSet<R>& initial_set, 
                          const Geometry::GridMaskSet<R>& safe_set) const
    {
      typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
      typedef typename Geometry::ListSet<R,Geometry::Parallelotope>::const_iterator pls_const_iterator;
      assert(initial_set.bounded() && safe_set.bounded());
     
      if(!subset(initial_set,safe_set)) {
        throw std::runtime_error("chainreach: Initial set must be subset of bounding set");
      }
        
      const Geometry::Grid<R>& g=initial_set.grid();
      const Geometry::GridMaskSet<R>& is=initial_set;
      const Geometry::Rectangle<R>& bb=safe_set.bounding_box();
      const Combinatoric::LatticeBlock lb=over_approximation(bb,g).lattice_set();
      
      Geometry::GridMaskSet<R> chainreach(g,lb);
      Geometry::GridCellListSet<R> found(g);
      Geometry::GridCellListSet<R> image(g);
      Geometry::GridCellListSet<R> cellimage(g);
      found.adjoin(is);
      
      time_type step_size=this->maximum_step_size();
      time_type time_step=this->lock_to_grid_time();
      
      while(!subset(found,chainreach)) {
        found=difference(found,chainreach);
        chainreach.adjoin(found);
        image.clear();
        uint size=0;
        Geometry::ListSet<R,Geometry::Parallelotope> parallotope_list;
        for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
          ++size;
          Geometry::Rectangle<R> r=*iter;
          Geometry::Parallelotope<R> pp(r);
          parallotope_list.adjoin(pp);
        }
        parallotope_list=this->integrate_list_set(vf,parallotope_list,time_step);
        for(pls_const_iterator iter=parallotope_list.begin(); iter!=parallotope_list.end(); ++iter) {
          Geometry::Parallelotope<R> fp=*iter;
          cellimage=over_approximation(fp,g);
          if(!subset(cellimage,safe_set)) {
            return false;
          }
          image.adjoin(cellimage);
        }
        found=image;
      }
      return true;
    }

  }
}
