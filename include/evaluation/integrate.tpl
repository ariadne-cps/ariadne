/***************************************************************************
 *            integrate.tpl
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

#include "../utility/stlio.h"
#include "../base/array.h"
#include "../numeric/arithmetic.h"
#include "../numeric/function.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/interval_vector.h"
#include "../linear_algebra/zonotopic_vector.h"

#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_matrix.h"

#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid.h"
#include "../geometry/grid_set.h"

#include "../evaluation/vector_field.h"
#include "../evaluation/affine_vector_field.h"

#include "../evaluation/integration_step.h"
#include "../evaluation/integrate.h"

namespace Ariadne {
  namespace Evaluation {

    template<typename R>
    Geometry::Rectangle<R> 
    integrate(const VectorField<R>& vector_field, 
              const Geometry::Rectangle<R>& initial_set, 
              const R& time,
              const R& step_size) 
    {
      if(time==0) { 
        return initial_set;
      }
      
      const VectorField<R>& vf=vector_field;
      Geometry::Rectangle<R> r(initial_set);
      R t=time;
      R h=step_size;
      while(t>0) {
        h=min(t,h);
        r=integration_step(vf,r,h);
        t=t-h;
        //h=max(R(2*h),step_size);
      }
      return r;
    }
    

    template<typename R>
    Geometry::Parallelotope<R> 
    integrate(const VectorField<R>& vector_field, 
              const Geometry::Parallelotope<R>& initial_set, 
              const R& time,
              const R& step_size) 
    {
      if(time==0) { 
        return initial_set;
      }
      
      const VectorField<R>& vf=vector_field;
      Geometry::Parallelotope<R> p(initial_set);
      R t=time;
      R h=step_size;
      while(t>0) {
        h=min(t,h);
        p=integration_step(vf,p,h);
        t=t-h;
        //h=max(R(2*h),step_size);
      }
      return p;
    }
    
    
    
/*
    template<typename R, template<typename> class BS>
      bool operator<(const std::pair< R,BS<R> >& ts1, const std::pair< R,BS<R> >& ts2) {
        return ts1.first < ts2.first; 
      }
*/
      
    template<typename R, typename BST>
    struct pair_first_less {
      bool operator()(const std::pair<R,BST>& ts1, const std::pair<R,BST>& ts2) {
        return ts1.first < ts2.first; 
      }
    };
    
    template<typename R, template<typename> class BS>
    Geometry::ListSet<R,BS> 
    integrate(const VectorField<R>& vector_field, 
              const Geometry::ListSet<R,BS>& initial_set, 
              const R& time, 
              const IntegrationParameters<R>& parameters)
    {
      if(time==0) { 
        return initial_set;
      }

      const VectorField<R>& vf=vector_field;
      R step_size=parameters.step_size;
      R maximum_set_radius=parameters.maximum_set_radius;
#ifdef DEBUG
        std::cerr << "step_size=" << step_size << "  maximum_set_radius=" << maximum_set_radius << std::endl;
#endif
      
      R t=time;
      R h=step_size;
      BS<R> bs(initial_set.dimension());
      
      std::multiset< std::pair<R,BS<R> >, pair_first_less<R,BS<R> > > working_sets;
      Geometry::ListSet<R,BS> final_set(initial_set.dimension());
      
      typedef typename Geometry::ListSet<R,BS>::const_iterator list_set_const_iterator;
      typedef std::pair< R,BS<R> > timed_set;
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
            bs=integration_step(vf,bs,h);
            t=t-h;
            h=min(R(2*h),step_size);
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
    

    template<typename R>
    Geometry::GridMaskSet<R> 
    integrate(const VectorField<R>& vector_field, 
              const Geometry::GridMaskSet<R>& initial_set,
              const Geometry::GridMaskSet<R>& bounding_set,
              const R& time,
              const IntegrationParameters<R>& parameters) 
    {
      assert(initial_set.grid()==bounding_set.grid());
      using namespace Geometry;
      using namespace LinearAlgebra;
      
      if(time==0) { 
        return initial_set;
      }
      
      R step_size=parameters.step_size;
      
      const Geometry::Grid<R>& g(initial_set.grid());
      Geometry::LatticeRectangle lb=bounding_set.bounds();
      Geometry::Rectangle<R> bb=bounding_set.bounding_box();
      
      Geometry::GridMaskSet<R> result(bounding_set.grid(),bounding_set.bounds());
      
      R t=time;
      R h=step_size;
      
      //bb=regular_intersection(bb,bounding_box());
      
      R spacial_tolerance=2.0;
      
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
          p=integrate(vf,p,h,h);
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
        GridCellListSet<R> oai=over_approximation_of_intersection(*iter,bb,g);
        result.adjoin(oai);
      }
      return result;
    }
    
    
    template<typename R, template<typename> class BSR, template<typename> class BSA>
    Geometry::ListSet<R,BSR> 
    reach(const VectorField<R>& vector_field, 
          const Geometry::ListSet<R,BSA>& initial_set, 
          const R& time,
          const IntegrationParameters<R>& parameters) 
    {
      const VectorField<R>& vf=vector_field;
      R step_size=parameters.step_size;
      R maximum_set_radius=parameters.maximum_set_radius;
#ifdef DEBUG
        std::cerr << "step_size=" << step_size << "  maximum_set_radius=" << maximum_set_radius << std::endl;
#endif
      
      R t=time;
      R h=step_size;
      BSA<R> bs(initial_set.dimension());
      BSR<R> rs(initial_set.dimension());
      
      std::multiset< std::pair<R,BSA<R> >, pair_first_less<R,BSA<R> > > working_sets;
      Geometry::ListSet<R,BSR> reach_set(initial_set.dimension());
      
      typedef typename Geometry::ListSet<R,BSA>::const_iterator list_set_const_iterator;
      typedef std::pair< R,BSA<R> > timed_set;
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
          Geometry::ListSet<R,BSA> subdivisions=bs.subdivide();
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
            rs=reach_step(vf,BSR<R>(bs),h);
            reach_set.adjoin(rs);
#ifdef DEBUG
        std::cerr << "rs.centre=" << rs.centre() << "  radius=" << rs.radius() 
                  << "  new size=" << reach_set.size() << std::endl;
#endif
            
            bs=integration_step(vf,bs,h);
            t=t-h;
            h=min(R(2*h),step_size);
          } while(t!=0 && bs.radius()<=maximum_set_radius);
          if(t==0) {
          } else {
            working_sets.insert(timed_set(t,bs));
          }
        }
      }
      
      return reach_set;
    }

    template<typename R>
    Geometry::GridMaskSet<R> 
    reach(const VectorField<R>& vector_field, 
          const Geometry::GridMaskSet<R>& initial_set,
          const Geometry::GridMaskSet<R>& bounding_set,
          const R& time,
          const IntegrationParameters<R>& parameters) 
    {
      typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
      typedef typename Geometry::ListSet<R,Geometry::Zonotope>::const_iterator zls_const_iterator;
      
      if(!subset(initial_set,bounding_set)) {
        throw std::runtime_error("chainreach: Initial set must be subset of bounding set");
      }
        
      const Geometry::Grid<R>& g=initial_set.grid();
      const Geometry::GridMaskSet<R>& is=initial_set;
      const Geometry::Rectangle<R>& bb=bounding_set.bounding_box();
      const Geometry::LatticeRectangle lb=over_approximation(bb,g).position();
      
      Geometry::GridMaskSet<R> stored(g,lb);
      Geometry::GridMaskSet<R> found(g,lb);
      Geometry::GridMaskSet<R> image(g,lb);
      Geometry::GridMaskSet<R> result(g,lb);
      found.adjoin(is);
      
      int steps=quotient(time,parameters.lock_to_grid_time);
      R time_step=time/steps;
      
      for(int step=0; step!=steps; ++step) {
        found=difference(found,stored);
        stored.adjoin(found);
        image.clear();
        Geometry::GridMaskSet<R> image=integrate(vector_field,found,bounding_set,time_step,parameters);
        found=image;
      }
      
      Geometry::ListSet<R,Geometry::Parallelotope> parallelotope_list;
      Geometry::ListSet<R,Geometry::Zonotope> zonotope_list;
      for(gms_const_iterator iter=stored.begin(); iter!=stored.end(); ++iter) {
        Geometry::Rectangle<R> r=*iter;
        Geometry::Parallelotope<R> pp(r);
        parallelotope_list.adjoin(pp);
      }
      zonotope_list=reach<R,Geometry::Zonotope,Geometry::Parallelotope>(vector_field,parallelotope_list,time_step,parameters);
      for(zls_const_iterator iter=zonotope_list.begin(); iter!=zonotope_list.end(); ++iter) {
        Geometry::Zonotope<R> fz=*iter;
        Geometry::GridCellListSet<R> oai=over_approximation_of_intersection(fz,bb,g);
        result.adjoin(oai);
      }
      return result;
    }
    
    
    
    template<typename R>
    Ariadne::Geometry::GridMaskSet<R> 
    chainreach(const VectorField<R>& vf, 
               const Geometry::GridMaskSet<R>& initial_set, 
               const Geometry::GridMaskSet<R>& bounding_set, 
               const IntegrationParameters<R>& parameters)
    {
      typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
      typedef typename Geometry::ListSet<R,Geometry::Parallelotope>::const_iterator pls_const_iterator;
      
      if(!subset(initial_set,bounding_set)) {
        throw std::runtime_error("chainreach: Initial set must be subset of bounding set");
      }
        
      const Geometry::Grid<R>& g=initial_set.grid();
      const Geometry::GridMaskSet<R>& is=initial_set;
      const Geometry::Rectangle<R>& bb=bounding_set.bounding_box();
      const Geometry::LatticeRectangle lb=over_approximation(bb,g).position();
      
      Geometry::GridMaskSet<R> result(g,lb);
      Geometry::GridMaskSet<R> found(g,lb);
      Geometry::GridMaskSet<R> image(g,lb);
      found.adjoin(is);
      
      R step_size=parameters.step_size;
      R time_step=parameters.lock_to_grid_time;
      
      while(!subset(found,result)) {
        found=difference(found,result);
        result.adjoin(found);
        image.clear();
        uint size=0;
        Geometry::ListSet<R,Geometry::Parallelotope> parallotope_list;
        for(gms_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
          ++size;
          Geometry::Rectangle<R> r=*iter;
          Geometry::Parallelotope<R> pp(r);
          parallotope_list.adjoin(pp);
        }
        parallotope_list=integrate(vf,parallotope_list,time_step,parameters);
        for(pls_const_iterator iter=parallotope_list.begin(); iter!=parallotope_list.end(); ++iter) {
          Geometry::Parallelotope<R> fp=*iter;
          Geometry::GridCellListSet<R> oai=over_approximation_of_intersection(fp,bb,g);
          image.adjoin(oai);
        }
        found=image;
      }
      return result;
    }

  }
}
