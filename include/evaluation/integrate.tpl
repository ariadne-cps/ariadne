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
#include "../linear_algebra/interval_vector.h"
#include "../linear_algebra/zonotopic_vector.h"

#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_matrix.h"

#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_set.h"

#include "../evaluation/evaluation_declarations.h"
#include "../evaluation/vector_field.h"

namespace Ariadne {
  namespace Evaluation {
    
    LinearAlgebra::matrix< Interval<double> >
    approximate_matrix(const LinearAlgebra::matrix< Interval<Rational> >& A) {
      LinearAlgebra::matrix< Interval<double> > result(A.size1(),A.size2());
      for(size_type i=0; i!=result.size1(); ++i) {
        for(size_type j=0; j!=result.size2(); ++j) {
          double l=convert_to<double>(A(i,j).lower());
          double u=convert_to<double>(A(i,j).upper());
          result(i,j)=Interval<double>(l,u);
        }
      }
      return result;
    }

    LinearAlgebra::matrix<double>
    approximate_matrix(const LinearAlgebra::matrix<Rational>& A) {
      LinearAlgebra::matrix<double> result(A.size1(),A.size2());
      for(size_type i=0; i!=result.size1(); ++i) {
        for(size_type j=0; j!=result.size2(); ++j) {
          result(i,j)=convert_to<double>(A(i,j));
        }
      }
      return result;
    }

/*
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
*/
    
    template<typename R>
    Geometry::Rectangle<R>
    estimate_flow_bounds(const Evaluation::VectorField<R>& vf,
                         const Geometry::Rectangle<R>& r,
                         R& h)
    {
      using namespace Geometry;
      Rectangle<R> estimate(vf.dimension());
      Rectangle<R> bounds(vf.dimension());
      std::cerr << h/2 << " " << r << std::endl;
      estimate=r+Interval<R>(0,2*h)*vf.apply(r);
      bounds=r+Interval<R>(0,h)*vf.apply(estimate);
      while(!subset(bounds,estimate)) {
        h=h/2;
        estimate=bounds;
        bounds=r+Interval<R>(0,h)*vf.apply(estimate);
        std::cerr << h << " " << bounds << " " << estimate << std::endl;
      };
      assert(subset(r+Interval<R>(0,h)*vf.apply(bounds),bounds));
      return bounds;
    }
    
    template<typename R>
    Geometry::Rectangle<R>
    estimate_flow_derivative_bounds(const Evaluation::VectorField<R>& vf,
                                    const Geometry::Rectangle<R>& r,
                                    R& h)
    {
      Geometry::Rectangle<R> estimate(vf.dimension());
      Geometry::Rectangle<R> bounds(vf.dimension());
      h=2*h;
      std::cerr << h/2 << " " << r << " " << vf.apply(r) << std::endl;
      do {
        estimate=r+Interval<R>(0,h)*vf.apply(r);
        h=h/2;
        bounds=r+Interval<R>(0,h)*vf.apply(estimate);
        std::cerr << h << " " << bounds << " " << estimate << std::endl;
      } while(!subset(bounds,estimate));
      return bounds;
    }
    
    template<typename R>
    Geometry::Rectangle<R>
    refine_flow_bounds(const Evaluation::VectorField<R>& vf,
                       const Geometry::Point<R>& x,
                       const Geometry::Rectangle<R>& b,
                       const R& h)
    {
      /* TODO: Use higher order method. */
      using namespace Geometry;
      using namespace LinearAlgebra;
      Rectangle<R> rx(x,x);
      Rectangle<R> xb=rx+Interval<R>(0,h)*vf.apply(b);
      Rectangle<R> xxb=rx+Interval<R>(0,h)*vf.apply(xb);
      std::cerr << "new bounds " << xxb << "," << xb << " vs old bounds " << b << "  " << subset(xb,b) << std::endl;
      vector< Interval<R> > ddphi=vf.derivative(xb)*vf.apply(xb);
      vector< Interval<R> > dfx=vf.apply(x);
      vector< Interval<R> > hdfx=(h*dfx);
      vector< Interval<R> > hhddphi=(R(h*h/2)*ddphi);
      vector< Interval<R> > dx=hdfx+hhddphi;
      return rx+dx;
    }
    
    
    /*! An inefficient C0 algorithm for integrating forward a rectangle. */
    template<typename R>
    Geometry::Rectangle<R> 
    integrate(const VectorField<R>& vector_field, 
              const Geometry::Rectangle<R>& initial_set, 
              const R& time) 
    {
      using namespace Geometry;
      using namespace LinearAlgebra;
      
      std::cerr << "integrate(const VectorField<R>& vf, const Geometry::Rectangle<R>& r, const R& t)" << std::endl;
      assert(vector_field.dimension()==initial_set.dimension());
      
      const VectorField<R>& vf(vector_field);
      Rectangle<R> r=initial_set;
      R t=time;
      R h=time;

      while(t>0) {
        std::cerr << "time left=" << t << std::endl;
        Geometry::Rectangle<R> q=estimate_flow_bounds(vf,r,h);
                
        std::cerr << "stepsize=" << h << std::endl;
        std::cerr << "bound=" << q << std::endl;

        LinearAlgebra::vector< Interval<R> > fq=vf.apply(q);

        std::cerr << "derivative=" << fq << std::endl;

        r=r+h*fq;
        std::cerr << "position=" << r << std::endl;

        t=t-h;
        h=std::min(R(2*h),t);
      }
      return r;
    }

    template<typename R>
    Geometry::Parallelotope<R> 
    integrate(const VectorField<R>& vector_field, 
              const Geometry::Parallelotope<R>& initial_set, 
              const R& time) 
    {
      typedef typename numerical_traits<R>::field_extension_type F;
      std::cerr << "integrate(VectorField<R>, Parallelotope<R>, R)\n";

      using namespace LinearAlgebra;
      using namespace Geometry;

      assert(vector_field.dimension()==initial_set.dimension());

      const VectorField<R>& vf(vector_field);
      Parallelotope<R> p=initial_set;
      const size_type n=p.dimension();
      R t=time;
      R h=time;
      const matrix<R> id=identity_matrix<R>(n);
      
      while(t>0) {
        std::cerr << "time left=" << t << std::endl;
        Rectangle<R> b=estimate_flow_bounds(vf,p.bounding_box(),h);
                
        std::cerr << "stepsize=" << h << std::endl;
        std::cerr << "bound=" << b << std::endl;
      
        vector< Interval<R> > f=vf.apply(b);
        std::cerr << "flow=" << f << std::endl;
        matrix< Interval<R> > df=vf.derivative(b);
        std::cerr << "jacobian=" << df << std::endl;
        matrix< Interval<R> > dphi=id+h*df;
        std::cerr << "flow derivative=" << dphi << std::endl;
        
        Point<R> c=p.centre();
        std::cerr << "centre=" << c << std::endl;
        Rectangle<R> phic=refine_flow_bounds(vf,c,b,h);
        std::cerr << "bounds on centre=" << phic << std::endl;
        
        /* TODO: Make this code more accurate by computing Minkowski sum zonotope. */
        matrix<R> B(n,n);
        matrix<F> BF(n,n);
        for(size_type i=0; i!=n; ++i) {
          for(size_type j=0; j!=n; ++j) {
            B(i,j)=(dphi(i,j).upper()+dphi(i,j).lower())/2;
            BF(i,j)=B(i,j);
          }
        }
        matrix<F> Binv=inverse(BF);
        matrix< Interval<F> > E(n,n);
        for(size_type i=0; i!=n; ++i) {
          for(size_type j=0; j!=n; ++j) {
            for(size_type k=0; k!=n; ++k) {
              Interval<F> invl(dphi(i,k).lower(),dphi(i,k).upper());
              E(i,j)+=invl*Binv(k,j);
            }
          }
        }
        
        std::cerr << "flow derivative centre=" << B << std::endl;
        std::cerr << "flow derivative centre inverse=" << approximate_matrix(Binv) << std::endl;
        std::cerr << "error matrix=" << approximate_matrix(E) << std::endl;
        
        matrix<F> Ainv=inverse(p.generators());
        matrix<F> Rerr(n,n);
        for(size_type i=0; i!=n; ++i) {
          const Rectangle<R>& cnstphic(phic);
          Interval<R> intvl=cnstphic[i];
          Rerr(i,i)=(intvl.upper()-intvl.lower())/2;
        }
        matrix<F> RerrAinvBinv=Rerr*Ainv*Binv;
        std::cerr << "relative rectangle error" << approximate_matrix(RerrAinvBinv) << std::endl;
        E+=RerrAinvBinv;
        std::cerr << "total error matrix=" << approximate_matrix(E) << std::endl;
        
        F excess=norm(E);
        uint precision=log_floor(2,F(1/(excess-1)));
        R approx_excess=1;
        
        p=Parallelotope<R>(phic.centre(),approx_excess*(B*p.generators()));
        
        t=t-h;
        h*=2;
        h=std::min(h,t);
      }
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
      return p;
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
