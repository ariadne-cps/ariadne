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
#include "../geometry/zonotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_set.h"

#include "../evaluation/evaluation_declarations.h"
#include "../evaluation/vector_field.h"
#include "../evaluation/affine_vector_field.h"

namespace Ariadne {
  namespace Evaluation {
    
    template<typename R>
    class IntervalParallelotope
    {
     public:
      IntervalParallelotope(const LinearAlgebra::interval_vector<R>& c, const LinearAlgebra::interval_matrix<R>& A) : _c(c), _A(A) { }
      
      size_type size() const { return _c.size(); }
      
      const LinearAlgebra::interval_vector<R>& centre() const { return _c; }
      LinearAlgebra::interval_vector<R>& centre() { return _c; }

      const LinearAlgebra::interval_matrix<R>& generators() const { return _A; }
      LinearAlgebra::interval_matrix<R>& generators() { return _A; }

      Geometry::Parallelotope<R> over_approximating_parallelotope() const;
     private:
      LinearAlgebra::interval_vector<R> _c;
      LinearAlgebra::interval_matrix<R> _A;
    };
   
    template<typename R>
    Geometry::Parallelotope<R> 
    IntervalParallelotope<R>::over_approximating_parallelotope() const 
    {
      std::cerr << "IntervalParallelotope<R>::over_approximating_parallelotope() const" << std::endl;
      typedef typename numerical_traits<R>::field_extension_type F;
      
      size_type n=this->size();
      
      LinearAlgebra::interval_matrix<R> A=this->generators();
      
      LinearAlgebra::vector<R> cmid=LinearAlgebra::centre(this->centre());
      LinearAlgebra::matrix<R> Amid=LinearAlgebra::centre(this->generators());
      
      LinearAlgebra::matrix<R> D(n,n);
      for(size_type i=0; i!=n; ++i) {
        D(i,i)=(this->centre()(i).upper()-this->centre()(i).lower())/2;
      }
      
      /* FIXME: Don't hard-code this error! */
      LinearAlgebra::matrix<F> AinvF=LinearAlgebra::inverse(Amid);
      LinearAlgebra::interval_matrix<R> Ainv=LinearAlgebra::approximate(AinvF,Dyadic(1)/65536);
      
      R err = upper_norm(Ainv*LinearAlgebra::interval_matrix<R>(D+A));
      std::cerr << "error=" << err << std::endl;
      return Geometry::Parallelotope<R>(cmid,err*Amid);
    }
    
    template<typename R>
    std::ostream& 
    operator<<(std::ostream& os, const IntervalParallelotope<R>& ip)
    {
      return os << "IntervalParallelotope(\n  centre=" << ip.centre() << "\n  generators=" << ip.generators() << "\n)\n";
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
        LinearAlgebra::interval_vector<R> df=vf.apply(bounds);
        
        F dt=t;
        for(dimension_type i=0; i!=vf.dimension(); ++i) {
          if(df(i).upper()>0) {
            dt=min(dt,F((bounds[i].upper()-reach[i].upper())/df(i).upper()));
          }
          if(df(i).lower()<0) {
            dt=min(dt,F((bounds[i].lower()-reach[i].lower())/df(i).lower()));
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
      std::cerr << "estimate_flow_bounds" << std::endl;
      using namespace Geometry;
      Rectangle<R> estimate(vf.dimension());
      Rectangle<R> bounds(vf.dimension());
      std::cerr << h/2 << " " << r << std::endl;
      std::cerr << LinearAlgebra::interval_vector<R>(Interval<R>(0,2*h)*vf.apply(r)) << std::endl;
      estimate=r+Interval<R>(0,2*h)*vf.apply(r);
      std::cerr << "estimate=" << estimate << std::endl;
      bounds=r+Interval<R>(0,h)*vf.apply(estimate);
      std::cerr << "bounds=" << bounds << std::endl;
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
    refine_flow_bounds(const Evaluation::VectorField<R>& vf,
                       const Geometry::Point<R>& x,
                       const Geometry::Rectangle<R>& b,
                       const R& h)
    {
      using namespace Geometry;
      using namespace LinearAlgebra;
      Rectangle<R> rx(x,x);
      Rectangle<R> xb=rx+Interval<R>(0,h)*vf.apply(b);
      Rectangle<R> xxb=rx+Interval<R>(0,h)*vf.apply(xb);
      std::cerr << "new bounds " << xxb << "," << xb << " vs old bounds " << b << "  " << subset(xb,b) << std::endl;
      interval_vector<R> ddphi=vf.derivative(xb)*vf.apply(xb);
      interval_vector<R> dfx=vf.apply(x);
      interval_vector<R> hdfx=(h*dfx);
      interval_vector<R> hhddphi=(R(h*h/2)*ddphi);
      interval_vector<R> dx=hdfx+hhddphi;
      return rx+dx;
    }
    
    
    template<typename R>
    Geometry::Rectangle<R> 
    integration_step(const VectorField<R>& vector_field, 
                     const Geometry::Rectangle<R>& initial_set, 
                     R& step_size) 
    {
      using namespace Geometry;
      using namespace LinearAlgebra;
      
      std::cerr << "integrate(const VectorField<R>& vf, const Geometry::Rectangle<R>& r, const R& t)" << std::endl;
      assert(vector_field.dimension()==initial_set.dimension());
      
      const VectorField<R>& vf(vector_field);
      Rectangle<R> r=initial_set;
      R& h=step_size;

      std::cerr << "suggested stepsize=" << step_size << std::endl;
      Geometry::Rectangle<R> q=estimate_flow_bounds(vf,r,h);
                
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << q << std::endl;

      LinearAlgebra::interval_vector<R> fq=vf.apply(q);
      std::cerr << "derivative=" << fq << std::endl;

      r=r+h*fq;
      std::cerr << "position=" << r << std::endl;

      return r;
    }

    template<typename R>
    Geometry::Parallelotope<R> 
    integration_step(const VectorField<R>& vector_field, 
                     const Geometry::Parallelotope<R>& initial_set, 
                     R& step_size) 
    {
      std::cerr << "integration_step(VectorField<R>, Parallelotope<R>, R)\n";
      const VectorField<R>* cvf_ptr=&vector_field;
      const AffineVectorField<R>* cavf_ptr=dynamic_cast< const AffineVectorField<R>* >(cvf_ptr);
      std::cerr << typeid(*cvf_ptr).name() << "  " << std::cerr << ((cavf_ptr) ? typeid(*cavf_ptr).name() : "void") << std::endl;
      if(!cavf_ptr) { cavf_ptr=static_cast< const AffineVectorField<R>* >(cvf_ptr); }
      if(cavf_ptr) {
        return integration_step(*cavf_ptr,initial_set,step_size);
      }
      
      typedef typename numerical_traits<R>::field_extension_type F;
      throw std::runtime_error("integration_step(VectorField<R>, Parallelotope<R>, R) contains major bugs");
      
      using namespace LinearAlgebra;
      using namespace Geometry;

      assert(vector_field.dimension()==initial_set.dimension());

      const VectorField<R>& vf(vector_field);
      Parallelotope<R> p=initial_set;
      const size_type n=p.dimension();
      R& h=step_size;
      const matrix<R> id=identity_matrix<R>(n);
      
      std::cerr << "suggested stepsize=" << step_size << std::endl;
      Rectangle<R> b=estimate_flow_bounds(vf,p.bounding_box(),h);
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << b << std::endl;
      
      interval_vector<R> f=vf.apply(b);
      std::cerr << "flow=" << f << std::endl;
      interval_matrix<R> df=vf.derivative(b);
      std::cerr << "jacobian=" << df << std::endl;
      R l=upper_log_norm(df);
      std::cerr << "logarithmic_norm=" << l << std::endl;
      l=max(l,R(1));
      interval_matrix<R> bdphi=id+Interval<R>(0,R(l*h))*df;
      std::cerr << "flow derivative_bounds=" << bdphi << std::endl;
        
      Point<R> c=p.centre();
      std::cerr << "centre=" << c << std::endl;
      Rectangle<R> phic=refine_flow_bounds(vf,c,b,h);
      std::cerr << "bounds on centre=" << phic << std::endl;
        
      interval_matrix<R> dphi=bdphi;
      matrix<R> approx_dphi=over_approximation(dphi);
      std::cerr << "over_approximation of flow derivative=" << approx_dphi << std::endl;
      zonotopic_vector<R> zv(approx_dphi*p.generators());
      std::cerr << "zonotopic vector to add to image of centre=" << zv << std::endl;
      Zonotope<R> z=phic+zv;
      std::cerr << "approximating zonotope=" << z << std::endl;
      p=over_approximating_parallelotope(z);
      assert(subset(z,p));
      std::cerr << "new approximation=" << p << std::endl;

      return p;
    }

    template<typename R>
    Geometry::Parallelotope<R> 
    integration_step(const AffineVectorField<R>& vector_field, 
                     const Geometry::Parallelotope<R>& initial_set, 
                     R& step_size) 
    {
      std::cerr << "integration_step(AffineVectorField<R>, Parallelotope<R>, R)\n";
      const AffineVectorField<R>& vf=vector_field;
      Geometry::Parallelotope<R> p=initial_set;
      R& h=step_size;
      
      std::cerr << "parallelotope generators=" << p.generators() << std::endl;
      R max_error=LinearAlgebra::norm(p.generators())/256;
      assert(max_error>0);
      std::cerr << "maximum allowed error=" << max_error << std::endl;
      
      /* Write phi(x)=D x0 + P b */
      std::cerr << "jacobian=" << vf.A() << std::endl;
      std::cerr << "step size=" << h << std::endl;
      LinearAlgebra::matrix<R> D=LinearAlgebra::exp_Ah_approx(vf.A(),h,max_error);
      std::cerr << "approximate derivative=" << D << std::endl;
      LinearAlgebra::matrix<R> P=LinearAlgebra::exp_Ah_sub_id_div_A_approx(vf.A(),h,max_error);
      std::cerr << "twist=" << P << std::endl;
      
      LinearAlgebra::interval_matrix<R> iD(LinearAlgebra::interval_matrix<R>(D,max_error));
      LinearAlgebra::interval_matrix<R> iP(LinearAlgebra::interval_matrix<R>(P,max_error));
      std::cerr << "approximating derivative=" << iD << std::endl;
      std::cerr << "approximating twist=" << iP << std::endl;
      
      IntervalParallelotope<R> img(iD*p.centre().position_vector()+iP*vf.b(),iD*p.generators());
      std::cerr << "interval parallelotope=" << img << std::endl;
      p=img.over_approximating_parallelotope();      
      std::cerr << "parallelotope=" << p << std::endl;
      return p;      
    }
    
    template<typename R>
    Geometry::Rectangle<R> 
    reach_step(const VectorField<R>& vector_field, 
               const Geometry::Rectangle<R>& initial_set, 
               R& step_size) 
    {
      using namespace Geometry;
      using namespace LinearAlgebra;
      
      std::cerr << "integrate(const VectorField<R>& vf, const Geometry::Rectangle<R>& r, const R& t)" << std::endl;
      assert(vector_field.dimension()==initial_set.dimension());
      
      const VectorField<R>& vf(vector_field);
      Rectangle<R> r=initial_set;
      R& h=step_size;

      std::cerr << "suggested stepsize=" << step_size << std::endl;
      Geometry::Rectangle<R> q=estimate_flow_bounds(vf,r,h);
                
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << q << std::endl;

      LinearAlgebra::interval_vector<R> fq=vf.apply(q);

      std::cerr << "derivative=" << fq << std::endl;

      r=r+(Interval<R>(R(0),h)*fq);
      std::cerr << "position=" << r << std::endl;

      return r;
    }
    
    
    template<typename R>
    Geometry::Parallelotope<R> 
    reach_step(const VectorField<R>& vector_field, 
               const Geometry::Parallelotope<R>& initial_set, 
               R& step_size)
    {
      return Geometry::over_approximating_parallelotope(reach_step(vector_field,Geometry::Zonotope<R>(initial_set),step_size));
    }

    template<typename R>
    Geometry::Zonotope<R> 
    reach_step(const VectorField<R>& vector_field, 
               const Geometry::Zonotope<R>& initial_set, 
               R& step_size)
    {
      typedef typename numerical_traits<R>::field_extension_type F;
      std::cerr << "integration_step_to(VectorField<R>, Parallelotope<R>, R)\n";
      throw std::runtime_error("reach_step(VectorField<R>, Parallelotope<R>, R) contains major bugs");

      using namespace LinearAlgebra;
      using namespace Geometry;

      assert(vector_field.dimension()==initial_set.dimension());

      const VectorField<R>& vf(vector_field);
      Zonotope<R> z=initial_set;
      const size_type n=z.dimension();
      R h=step_size;
      const matrix<R> id=identity_matrix<R>(n);
      
      std::cerr << "suggested stepsize=" << step_size << std::endl;
      /* Throws exception if we can't find flow bounds for given stepsize. */
      Rectangle<R> b=estimate_flow_bounds(vf,z.bounding_box(),h);
        
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << b << std::endl;
      
      interval_vector<R> f=vf.apply(z.bounding_box());
      std::cerr << "flow=" << f << std::endl;
      interval_matrix<R> df=vf.derivative(b);
      std::cerr << "jacobian=" << df << std::endl;
      interval_matrix<R> dphi=id+Interval<R>(0,h)*df;
      std::cerr << "flow derivative=" << dphi << std::endl;
        
      Point<R> c=z.centre();
      std::cerr << "centre=" << c << std::endl;
      Rectangle<R> phic=refine_flow_bounds(vf,c,b,R(h/2));
      std::cerr << "bounds on centre=" << phic << std::endl;
      
      //interval_vector<R> fh=(R(h/2)*f);
      interval_vector<R> fh=(R(h/2)*f);
      std::cerr << "flow times stepsize=" << fh << std::endl;
      zonotopic_vector<R> zfh=symmetrise(fh);
      std::cerr << "symmetrised flow=" << fh << std::endl;
      matrix<R> mdf=over_approximation(dphi)*z.generators();
      std::cerr << "over approximating matrix=" << fh << std::endl;
      zonotopic_vector<R> zv=zfh+zonotopic_vector<R>(vector<R>(n),mdf);
      
      z=phic+zv;
      std::cerr << "approximating zonotope " << z;
      return z;
    }

    

    
    template<typename R, template<typename> class BS>
    Geometry::ListSet<R,BS> 
    integrate(const VectorField<R>& vf, const Geometry::ListSet<R,BS>& ds, const Interval<R>& time, const R& step_size)
    {
      R& h=step_size;
      Geometry::ListSet<R,BS> result(vf.dimension());
      for(typename Geometry::ListSet<R,BS>::const_iterator iter=ds.begin(); iter!=ds.end(); ++iter) {
        BS<R> bs=*iter;
        
        R t=time.lower();
        while(t>0) {
          h=max(t,h);
          bs=integration_step(vf,bs,h);
          t=t-h;
          h=max(2*h,step_size);
        }
        
        t=time.upper()-time.lower();
        while(t>0) {
          h=max(t,h);
          result.push_back(integration__reach_step(vf,bs,h));
          bs=integration_step(vf,bs,h);
          t=t-h;
          h=max(2*h,step_size);
        }
      }
    }
    
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
        h=max(t,h);
        r=integration_step(vf,r,h);
        t=t-h;
        h=max(R(2*h),step_size);
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
        h=max(t,h);
        p=integration_step(vf,p,h);
        t=t-h;
        h=max(R(2*h),step_size);
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

      Geometry::ListSet<R,BS> result(vector_field.dimension());

      const VectorField<R>& vf=vector_field;
      R h=step_size;
      for(typename Geometry::ListSet<R,BS>::const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
        BS<R> bs=*iter;
        
        R t=time;
        while(t>0) {
          h=max(t,h);
          bs=integration_step(vf,bs,h);
          t=t-h;
          h=max(R(2*h),step_size);
        }
        result.adjoin(bs);
      }
      return result;
    }
    

    template<typename R>
    Geometry::GridMaskSet<R> 
    integrate(const VectorField<R>& vector_field, 
              const Geometry::GridMaskSet<R>& initial_set, 
              const R& time,
              const R& step_size) 
    {
      if(time==0) { 
        return initial_set;
      }
      
      const Geometry::Grid<R>& g(initial_set.grid());
      const Geometry::Rectangle<R> bb(Geometry::GridRectangle<R>(initial_set.grid(),initial_set.bounds()));

      Geometry::GridMaskSet<R> result(vector_field.dimension());
      for(typename Geometry::GridMaskSet<R>::const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
        const VectorField<R>& vf=vector_field;
        Geometry::Parallelotope<R> p(*iter);
        R t=time;
        R h=step_size;
        while(t>0) {
          h=max(t,h);
          p=integration_step(vf,p,h);
          t=t-h;
          h=max(R(2*h),step_size);
        }
        result.adjoin(over_approximation_of_intersection(p,bb,g));
      }
      return result;
    }
    
    
    template<typename R>
    Geometry::GridMaskSet<R> 
    reach(const VectorField<R>& vector_field, 
              const Geometry::GridMaskSet<R>& initial_set, 
              const R& time,
              const R& step_size) 
    {
      const Geometry::Grid<R>& g(initial_set.grid());
      const Geometry::Rectangle<R> bb(Geometry::GridRectangle<R>(initial_set.grid(),initial_set.bounds()));
      
      Geometry::GridMaskSet<R> result(vector_field.dimension());
      for(typename Geometry::GridMaskSet<R>::const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
        const VectorField<R>& vf=vector_field;
        Geometry::Parallelotope<R> p(*iter);
        R t=time;
        R h=step_size;
        while(t>0) {
          h=max(t,h);
          result.adjoin(over_approximation_of_intersection(reach_step(vf,p,h),bb,g));
          p=integration_step(vf,p,h);
          t=t-h;
          h=max(R(2*h),step_size);
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
