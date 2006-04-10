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
    

    template<typename R, template<typename> class BS>
    Geometry::ListSet<R,BS> 
    integrate(const VectorField<R>& vector_field, const Geometry::ListSet<R,BS>& initial_set, const R& time, const R& step_size)
    {
      if(time==0) { 
        return initial_set;
      }

      R t=time;
      R h=step_size;
      R spacial_tolerance=0.1;
      
      Geometry::ListSet<R,BS> start=initial_set;
      Geometry::ListSet<R,BS> finish(initial_set.dimension());
      
      while(t!=0) {
        h=min(t,h);
        for(typename Geometry::ListSet<R,BS>::const_iterator iter=start.begin(); iter!=start.end(); ++iter) {
          const VectorField<R>& vf=vector_field;
          BS<R> bs(*iter);
          bs=integrate(vf,bs,h,h);
          finish.adjoin(bs);
        }
        std::cerr << "finish.size()=" << finish.size() << "  ";
        start.clear();
        for(typename Geometry::ListSet<R,BS>::const_iterator iter=finish.begin(); iter!=finish.end(); ++iter) {
          BS<R> bs=*iter;
          if(bs.radius()>spacial_tolerance) {
            start.adjoin(bs.subdivide());
          }
          else {
            start.adjoin(bs);
          }
        }
        std::cerr << "start.size()=" << start.size() << "\n";
        finish.clear();
        t-=h;
      }
      return start;
    }
    

    template<typename R>
    Geometry::GridMaskSet<R> 
    integrate(const VectorField<R>& vector_field, 
              const Geometry::GridMaskSet<R>& initial_set,
              const Geometry::GridMaskSet<R>& bounding_set,
              const R& time,
              const R& step_size) 
    {
      assert(initial_set.grid()==bounding_set.grid());
      using namespace Geometry;
      using namespace LinearAlgebra;
      
      if(time==0) { 
        return initial_set;
      }
      
      const Geometry::Grid<R>& g(initial_set.grid());
      Geometry::LatticeRectangle lb=bounding_set.bounds();
      Geometry::GridRectangle<R> bb=bounding_set.bounding_box();
      
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
        std::cerr << "time left=" << t << "  stepsize=" << h << "  sets in list=" << start_set.size() << "\n";
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
            std::cerr << "Splitting, radius=" << p.radius() << "\n" << p << "\n";
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
    
    
    template<typename R, template<typename> class BS>
    Geometry::ListSet<R,BS> 
    reach(const VectorField<R>& vector_field, 
          const Geometry::ListSet<R,BS>& initial_set, 
          const R& time,
          const R& step_size) 
    {
      throw std::runtime_error("reach(...): Not implemented.");
      
      return initial_set;
    }


    template<typename R>
    Geometry::GridMaskSet<R> 
    reach(const VectorField<R>& vector_field, 
              const Geometry::GridMaskSet<R>& initial_set, 
              const Geometry::GridMaskSet<R>& bounding_set, 
              const R& time,
              const R& step_size) 
    {
      const Geometry::Grid<R>& g(initial_set.grid());
      const Geometry::GridRectangle<R> gbb=Geometry::GridRectangle<R>(initial_set.grid(),initial_set.bounds());
      
      Geometry::GridMaskSet<R> result(vector_field.dimension());
      for(typename Geometry::GridMaskSet<R>::const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
        const VectorField<R>& vf=vector_field;
        Geometry::Parallelotope<R> p(*iter);
        R t=time;
        R h=step_size;
        while(t>0) {
          h=min(t,h);
          result.adjoin(over_approximation_of_intersection(reach_step(vf,p,h),gbb,g));
          p=integration_step(vf,p,h);
          t=t-h;
          //h=max(R(2*h),step_size);
        }
      }
      return result;
    }

/*
    template<typename R>
    Ariadne::Geometry::GridMaskSet<R> 
    chainreach(const VectorField<R>& f, 
               const Ariadne::Geometry::ListSet<R,Ariadne::Geometry::Rectangle>& is, 
               const Ariadne::Geometry::FiniteGrid<R>& g, 
               const Ariadne::Geometry::Rectangle<R>& bb) 
    {
      typedef typename Ariadne::Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
      
      Ariadne::Geometry::GridMaskSet<R> result(g);
      Ariadne::Geometry::GridMaskSet<R> found=over_approximation_of_intersection(is,bb,g);
      Ariadne::Geometry::GridMaskSet<R> image(g);
      
      R h(1); // FIXME: change stepsize
      
      while(!subset(found,result)) {
        found=difference(found,result);
        result.adjoin(found);
        image=Ariadne::Geometry::GridMaskSet<R>(g);
        uint size=0;
        for(gms_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
          ++size;
          Ariadne::Geometry::Rectangle<R> r=*iter;
          Ariadne::Geometry::Parallelotope<R> pp(r);
          Ariadne::Geometry::Parallelotope<R> fp(Ariadne::Evaluation::integrate(f,pp,h));
          Geometry::GridCellListSet<R> oai=over_approximation_of_intersection(fp,bb,g);
          image.adjoin(oai);
        }
        std::cerr << size << std::endl;
        found=image;
      }
      return result;
    }
*/
  }
}
