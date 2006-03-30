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
    
    template<typename R>
    Geometry::Zonotope<R>
    operator+(const Geometry::Zonotope<R>& z, const LinearAlgebra::zonotopic_vector<R>& v)
    {
      return Geometry::Zonotope<R>(z.centre()+v.centre(),LinearAlgebra::concatenate_columns(z.generators(),v.generators()));
    }
    
    template<typename R>
    Geometry::Zonotope<R>
    operator+(const Geometry::Rectangle<R>& r, const LinearAlgebra::zonotopic_vector<R>& v)
    {
      return Geometry::Zonotope<R>(r)+v;
    }
    
    template<typename R>
    LinearAlgebra::zonotopic_vector<R>
    symmetrise(const LinearAlgebra::vector< Interval<R> >& iv)
    {
      using namespace Geometry;
      using namespace LinearAlgebra;

      matrix<R> A(iv.size(),iv.size()+1);
      for(size_type i=0; i!=A.size1(); ++i) {
        A(i,i)=(iv[i].upper()-iv[i].lower())/2;
        A(i,iv.size())=(iv[i].upper()+iv[i].lower())/2;
      }
      return zonotopic_vector<R>(vector<R>(iv.size()),A);
    }
      
    template<typename R>
    LinearAlgebra::matrix<R>
    over_approximating_matrix(const LinearAlgebra::matrix< Interval<R> >& D, const LinearAlgebra::matrix<R>& A)
    {
      using namespace Geometry;
      using namespace LinearAlgebra;
      
      typedef typename numerical_traits<R>::field_extension_type F;
      dimension_type n=A.size1();
      
      matrix<R> B(n,n);
      matrix<F> BF(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          B(i,j)=(D(i,j).upper()+D(i,j).lower())/2;
          BF(i,j)=B(i,j);
        }
      }
      
      matrix<F> Binv=inverse(BF);
      matrix< Interval<F> > E(n,n);
      for(size_type i=0; i!=n; ++i) {
        for(size_type j=0; j!=n; ++j) {
          for(size_type k=0; k!=n; ++k) {
            Interval<F> invl(D(i,k).lower(),D(i,k).upper());
            E(i,j)+=invl*Binv(k,j);
          }
        }
      }
      
      std::cerr << "flow derivative centre=" << B << std::endl;
      std::cerr << "flow derivative centre inverse=" << Binv << std::endl;
      std::cerr << "error matrix=" << E << std::endl;
            
      F excess=norm(E);
      uint precision=log_floor(2,F(1/(excess-1)));
      R approx_excess=1+pow(Dyadic(0.5),precision);
    
      return approx_excess*(B*A);

    }
      
    /* Compute an over-approximation to z. */
    template<typename R>
    Geometry::Parallelotope<R>
    over_approximating_parallelotope(const Geometry::Zonotope<R>& z)
    {
      std::cerr << "over_approximating_parallelotope(const Geometry::Zonotope<R>& z)" << std::endl;
      dimension_type n=z.dimension();
      LinearAlgebra::matrix<R> A(n,n);
      for(dimension_type i=0; i!=n; ++i) {
        for(dimension_type j=0; j!=n; ++j) {
          A(i,j)=z.generators()(i,j);
        }
      }
      LinearAlgebra::matrix<R> B=LinearAlgebra::inverse(A);
      const Geometry::Point<R>& c=z.centre();
      
      Geometry::Parallelotope<R> p(c,A);
      while(!subset(z,p)) {
        std::cerr << "A=" << A << std::endl;
        A*=2;
        p=Geometry::Parallelotope<R>(c,A);
      }
      std::cerr << "A=" << A << std::endl;
      return p;
    }
      
    /* Compute an over-approximation to c+DAe. */
    template<typename R>
    Geometry::Parallelotope<R>
    over_approximating_parallelotope(const Geometry::Rectangle<R>& c,
                                     const LinearAlgebra::matrix< Interval<R> >& D,
                                     const LinearAlgebra::matrix<R>& A)
    {
      LinearAlgebra::matrix<R> DA=over_approximating_matrix(D,A);
      Geometry::Zonotope<R> z(c);
      DA=LinearAlgebra::concatenate_columns(DA,z.generators());
      z=Geometry::Zonotope<R>(c.centre(),DA);
      return over_approximating_parallelotope(z);
    }
    
    template<typename R>
    bool
    check_flow_bounds(const Evaluation::VectorField<R>& vf,
                      const Geometry::Rectangle<R>& r,
                      const Geometry::Rectangle<R>& b,
                      const R& h)
    {
      using namespace Geometry;
      return subset(r+Interval<R>(0,h)*vf.apply(b),b);
    }
    
    template<typename R>
    Geometry::Rectangle<R>
    compute_flow_bounds(const Evaluation::VectorField<R>& vf,
                        const Geometry::Rectangle<R>& r,
                        const R& h)
    {
      std::cerr << "\n\n\ncompute_flow_bounds\n";
      using namespace Geometry;
      typedef typename numerical_traits<R>::field_extension_type F;
      uint max_iterations=16;
      uint iteration=0;
      R multiplier=1.125;
      F t=h;
      Rectangle<R> reach=(vf.dimension());
      Rectangle<R> bounds(vf.dimension());
      std::cerr << h << " " << r << std::endl;
      reach=r;
      std::cerr << t << " " << reach << std::endl;
      while(t>0) {
        bounds=reach+Interval<R>(0,multiplier*h)*vf.apply(reach);
        LinearAlgebra::vector< Interval<R> > df=vf.apply(bounds);
        
        F dt=t;
        for(dimension_type i=0; i!=vf.dimension(); ++i) {
          if(df[i].upper()>0) {
            dt=std::min(dt,F((bounds[i].upper()-reach[i].upper())/df[i].upper()));
          }
          if(df[i].lower()<0) {
            dt=std::min(dt,F((bounds[i].lower()-reach[i].lower())/df[i].lower()));
          }
        }
        reach=bounds;
        t-=dt;
        
        std::cerr << t << " " << reach << std::endl;
        ++iteration;
        if(iteration==max_iterations) {
          throw std::runtime_error("Cannot find bounding box for flow");
        }
      }
      return reach;
    }
    
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
      }
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
        
        p=over_approximating_parallelotope(phic,dphi,p.generators());
        std::cerr << "new approximation=" << p << std::endl;
       
        t=t-h;
        h*=2;
        h=std::min(h,t);
      }

      return p;
    }

    template<typename R>
    Geometry::Parallelotope<R> 
    integrate_to(const VectorField<R>& vector_field, 
                 const Geometry::Parallelotope<R>& initial_set, 
                 const R& time)
    {
      typedef typename numerical_traits<R>::field_extension_type F;
      std::cerr << "integrate_to(VectorField<R>, Parallelotope<R>, R)\n";

      using namespace LinearAlgebra;
      using namespace Geometry;

      assert(vector_field.dimension()==initial_set.dimension());

      const VectorField<R>& vf(vector_field);
      Parallelotope<R> p=initial_set;
      const size_type n=p.dimension();
      R t=time;
      R h=time;
      const matrix<R> id=identity_matrix<R>(n);
      
      std::cerr << "time left=" << t << std::endl;
      /* Throws exception if we can't find flow bounds for given stepsize. */
      Rectangle<R> b=compute_flow_bounds(vf,p.bounding_box(),h);
        
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << b << std::endl;
      
      vector< Interval<R> > f=vf.apply(b);
      std::cerr << "flow=" << f << std::endl;
      matrix< Interval<R> > df=vf.derivative(b);
      std::cerr << "jacobian=" << df << std::endl;
      matrix< Interval<R> > dphi=id+Interval<R>(0,h)*df;
      std::cerr << "flow derivative=" << dphi << std::endl;
        
      Point<R> c=p.centre();
      std::cerr << "centre=" << c << std::endl;
      Rectangle<R> phic=refine_flow_bounds(vf,c,b,R(h/2));
      std::cerr << "bounds on centre=" << phic << std::endl;
      
      //interval_vector<R> fh=(R(h/2)*f);
      vector< Interval<R> > fh=(R(h/2)*f);
      std::cerr << "flow times stepsize=" << fh << std::endl;
      zonotopic_vector<R> zfh=symmetrise(fh);
      std::cerr << "symmetrised flow=" << fh << std::endl;
      matrix<R> mdf=over_approximating_matrix(dphi,p.generators());
      std::cerr << "over approximating matrix=" << fh << std::endl;
      zonotopic_vector<R> zv=zfh+zonotopic_vector<R>(vector<R>(n),mdf);
      
      Zonotope<R> z=phic+zv;
      std::cerr << "approximating zonotope " << z;
      
      p=over_approximating_parallelotope(z);
      assert(subset(z,p));
      std::cerr << "new approximation=" << p << std::endl;

      return p;
    }

    
    template<typename R>
    Geometry::Parallelotope<R> 
    integrate(const VectorField<R>& vector_field, 
              const Geometry::Parallelotope<R>& initial_set, 
              const Interval<R>& time) 
    {
      typedef typename numerical_traits<R>::field_extension_type F;
      std::cerr << "integrate(VectorField<R>, Parallelotope<R>, Interval<R>)\n";

      using namespace LinearAlgebra;
      using namespace Geometry;

      assert(vector_field.dimension()==initial_set.dimension());

      const VectorField<R>& vf(vector_field);
      Parallelotope<R> p=initial_set;
      R t1=time.lower();
      R t2=time.upper();
      
      p=integrate(vf,p,t1);
      p=integrate_to(vf,p,R(t2-t1));
      
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
