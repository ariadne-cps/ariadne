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

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../base/utility.h"
#include "../base/array.h"
#include "../base/arithmetic.h"
#include "../base/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../linear_algebra/interval_vector.h"
#include "../linear_algebra/interval_matrix.h"

#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_set.h"

#include "../evaluation/evaluation_declarations.h"
#include "../evaluation/vector_field.h"

namespace Ariadne {
  namespace Evaluation {
   
    template<typename R>
    Geometry::Rectangle<R> 
    move(const Geometry::Rectangle<R>& r, 
         const R& h,
         const LinearAlgebra::vector< Interval<R> >& iv)
    {
      Geometry::Rectangle<R> result(r.dimension());
      assert(r.dimension()==iv.size());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+h*iv[i]);
      }
      return result;
    }
    
    template<typename R>
    Geometry::Rectangle<R> 
    operator+(const Geometry::Rectangle<R>& r, 
              const LinearAlgebra::vector< Interval<R> >& iv)
    {
      Geometry::Rectangle<R> result(r.dimension());
      assert(r.dimension()==iv.size());
      
      for(size_type i=0; i!=result.dimension(); ++i) {
        result.set_interval(i,r[i]+iv[i]);
      }
      return result;
    }
    

    /*! An inefficient C0 algorithm for integrating forward a rectangle. */
    template<typename R>
    Geometry::Rectangle<R> 
    integrate(const VectorField<R>& vf, const Geometry::Rectangle<R>& r, const R& t) 
    {
      std::cerr << "integrate(const VectorField<R>& vf, const Geometry::Rectangle<R>& r, const R& t)" << std::endl;
      
      assert(vf.dimension()==r.dimension());
      const size_type n=r.dimension();
      Geometry::Rectangle<R> result=r;
      R time=t;
      
      while(time>0) {
        std::cerr << "time left=" << time << std::endl;
        R expansion=result.upper_bound(0)-result.lower_bound(0);
        for(dimension_type i=1; i!=n; ++i) {
          R radius=result.upper_bound(i)-result.lower_bound(i);
          if(radius<expansion) {
            expansion=radius;
          }
        }
        
        std::cerr << "expansion=" << expansion << std::endl;
        
        Geometry::Rectangle<R> q=result;
        q.expand_by(expansion);
        
        std::cerr << "bound=" << q << std::endl;

        LinearAlgebra::vector< Interval<R> > fq=vf.apply(q);

        std::cerr << "derivative=" << fq << std::endl;

        R vmax=0;
        for(dimension_type i=0; i!=n; ++i) {
          R vabs=Ariadne::abs(fq[i].lower());
          if(vabs>vmax) {
            vmax=vabs;
          }
          vabs=Ariadne::abs(fq[i].upper());
          if(vabs>vmax) {
            vmax=vabs;
          }
        }
        
        R h=std::min(R(expansion/vmax),time);
        std::cerr << "stepsize=" << h << std::endl;
        
/*
        for(dimension_type i=0; i!=n; ++i) {
          //Interval<R> ri=result[i];
          Interval<R> ri(result.lower_bound(i),result.upper_bound(i));
          Interval<R> fqi=fq(i);
          Interval<R> nri=ri+h*fqi;
          result.set_interval(i,nri);
        }
*/
        result=move(result,h,fq);
        std::cerr << "position=" << result << std::endl;

        time=time-h;
      }
      return result;
    }

    template<typename R>
    Geometry::Parallelotope<R> 
    integrate(const VectorField<R>& vf, const Geometry::Parallelotope<R>& p, const R& t) 
    {
      assert(vf.dimension()==p.dimension());
      const size_type n=p.dimension();
  
      R h=t;
      
      LinearAlgebra::vector< Interval<R> > f;
      Geometry::Rectangle<R> ibb=p.bounding_box();
      Geometry::Rectangle<R> obb=move(ibb,R(2*h),vf.apply(ibb));
      Geometry::Rectangle<R> bb=move(ibb,h,vf.apply(obb));
      while(!subset(bb,obb)) {
        std::cerr << "stepsize=" << h << std::endl;
        Geometry::Rectangle<R> obb=move(ibb,h,vf.apply(ibb));
        h=h/2;
        Geometry::Rectangle<R> bb=move(ibb,h,vf.apply(obb));
      }
      std::cerr << "stepsize=" << h << std::endl;
      
      f=vf.apply(obb);
      LinearAlgebra::matrix< Interval<R> > df=vf.derivative(obb);

      
/*
      LinearAlgebra::vector< Interval<R> > cuboid_vector(m);
      const Interval<R> unit_interval(-1,1);
      for(size_type i=0; i!=cuboid_vector.size(); ++i) {
        cuboid_vector[i]=Interval<R>(-1,1);
      }
            
      const Geometry::Point<R>& c=p.centre();
      const LinearAlgebra::matrix<R>& g=p.generators();
      
      Geometry::Point<R> img_centre=f.apply(c);
      LinearAlgebra::matrix< Interval<R> > df_on_set = f.derivative(p.bounding_box());
      LinearAlgebra::matrix<R> df_at_centre = f.derivative(c);
      
      LinearAlgebra::matrix<R> img_generators = df_at_centre*g;
      
      LinearAlgebra::matrix<R> img_generators_inverse = LinearAlgebra::inverse(img_generators);
      
      LinearAlgebra::matrix< Interval<R> > cuboid_transform = img_generators_inverse * df_on_set * g;
      
      LinearAlgebra::vector< Interval<R> > new_cuboid = cuboid_transform * cuboid_vector;
      
      R new_cuboid_sup(0);
      for(size_type j=0; j!=n; ++j) {
        new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid[j].lower())) );
        new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid[j].upper())) );
      }
*/
      
//      Geometry::Parallelotope<R> result(img_centre,img_generators);
      Geometry::Parallelotope<R> result(n);
      return result;
    }

    template<typename R, template<typename> class BS>
    Ariadne::Geometry::ListSet<R,BS> 
    integrate(const VectorField<R>& vf, const Ariadne::Geometry::ListSet<R,BS>& ds, const R& t) 
    {
      Ariadne::Geometry::ListSet<R,BS> result(vf.dimension());
      for(typename Geometry::ListSet<R,BS>::const_iterator iter=ds.begin(); iter!=ds.end(); ++iter) {
        result.push_back(integrate(vf,*iter,t));
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
