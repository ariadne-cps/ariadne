/***************************************************************************
 *            lohner_integrator.tpl
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

#include "lohner_integrator.h"

#include "../base/array.h"
#include "../numeric/arithmetic.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"

#include "../linear_algebra/matrix.h"
#include "../linear_algebra/matrix_function.h"

#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/zonotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid.h"
#include "../geometry/grid_set.h"

#include "../system/vector_field.h"
#include "../system/affine_vector_field.h"

#include "../evaluation/integrator.h"

namespace Ariadne {
  namespace Evaluation {
    
#ifdef DEBUG
    static const int debug_level=1;
#else
    static const int debug_level=0;
#endif
    
    template<typename R>
    C1LohnerIntegrator<R>::C1LohnerIntegrator(const R& maximum_step_size, const R& lock_to_grid_time, const R& maximum_basic_set_radius)
      : C1Integrator<R>(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
    {
    }

    template<typename R>
    Geometry::Rectangle<R> 
    C0LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                            const Geometry::Rectangle<R>& initial_set, 
                                            R& step_size) const
    {
#ifdef DEBUG      
      std::cerr << "integrate(const VectorField<R>& vf, const Geometry::Rectangle<R>& r, const R& t)" << std::endl;
#endif

      assert(vector_field.dimension()==initial_set.dimension());
      
      const System::VectorField<R>& vf(vector_field);
      Geometry::Rectangle<R> r=initial_set;
      R& h=step_size;
      
      Geometry::Rectangle<R> q=estimate_flow_bounds(vf,r,h);
      
      LinearAlgebra::Vector< Interval<R> > fq=vf(q);
      r=r+(h*fq);
#ifdef DEBUG
      std::cerr << "suggested stepsize=" << step_size << std::endl;
                
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << q << std::endl;

      std::cerr << "derivative=" << fq << std::endl;

      std::cerr << "position=" << r << std::endl;
#endif
      return r;
    }

    template<typename R>
    Geometry::Rectangle<R> 
    C0LohnerIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                             const Geometry::Rectangle<R>& initial_set, 
                                             R& step_size) const
    {
#ifdef DEBUG
      std::cerr << "integrate(const VectorField<R>& vf, const Geometry::Rectangle<R>& r, const R& t)" << std::endl;
#endif
      assert(vector_field.dimension()==initial_set.dimension());
      
      const System::VectorField<R>& vf(vector_field);
      Geometry::Rectangle<R> r=initial_set;
      R& h=step_size;

      Geometry::Rectangle<R> q=estimate_flow_bounds(vf,r,h);
      
      LinearAlgebra::Vector< Interval<R> > fq=vf(q);
      
      r=r+LinearAlgebra::Vector< Interval<R> >(Interval<R>(R(0),h)*fq);

#ifdef DEBUG
      std::cerr << "suggested stepsize=" << step_size << std::endl;
                
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << q << std::endl;

      std::cerr << "derivative=" << fq << std::endl;
      std::cerr << "position=" << r << std::endl;
#endif

      return r;
    }
    
    template<typename R>
    Geometry::Parallelotope<R> 
    C1LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                            const Geometry::Parallelotope<R>& initial_set, 
                                            R& step_size) const
    {
      const System::VectorField<R>* cvf_ptr=&vector_field;
      const System::AffineVectorField<R>* cavf_ptr=dynamic_cast< const System::AffineVectorField<R>* >(cvf_ptr);

#ifdef DEBUG
      if(cavf_ptr) { std::cerr << typeid(*cavf_ptr).name(); } else { std::cerr << "void"; } std::cerr << std::endl;
      std::cerr << cvf_ptr << "  " << cavf_ptr << "\n" << std::endl;
#endif
      if(cavf_ptr != NULL) {
        return integration_step(*cavf_ptr,initial_set,step_size);
      }

#ifdef DEBUG
      std::cerr << "integration_step(VectorField<R>, Parallelotope<R>, R)\n";
#endif

      typedef typename numerical_traits<R>::field_extension_type F;
       
      assert(vector_field.dimension()==initial_set.dimension());

      const System::VectorField<R>& vf(vector_field);
      Geometry::Parallelotope<R> p=initial_set;
      const size_type n=p.dimension();
      R& h=step_size;
      const LinearAlgebra::Matrix<R> id=LinearAlgebra::identity_matrix<R>(n);
      
      R err=norm(p.generators())/65536;
      
      Geometry::Rectangle<R> b=this->estimate_flow_bounds(vf,p.bounding_box(),h);
#ifdef DEBUG
      std::cerr << "suggested stepsize=" << step_size << std::endl;
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << b << std::endl;
#endif
      
      LinearAlgebra::Vector< Interval<R> > f=vf(b);
      LinearAlgebra::Matrix< Interval<R> > df=vf.derivative(b);
      
      R l=df.log_norm().upper();

#ifdef DEBUG
      std::cerr << "flow=" << f << std::endl;
      std::cerr << "jacobian=" << df << std::endl;
      std::cerr << "jacobian centre=" << df.centre() << std::endl;
      std::cerr << "logarithmic_norm=" << l << std::endl;
#endif

      l=max(l,R(1));
      
      //integration_step(): the varaible l is useless!
      
      LinearAlgebra::Matrix< Interval<R> > hdf=h*df;
      LinearAlgebra::Matrix< Interval<R> > dphi=exp(hdf);
      
      Geometry::Point<R> c=p.centre();
      Geometry::Rectangle<R> rc=Geometry::Rectangle<R>(c,c);
      Geometry::Rectangle<R> phic=refine_flow_bounds(vf,rc,b,h);
      LinearAlgebra::Matrix< Interval<R> > zv(dphi*p.generators());

      p=Geometry::Parallelotope<R>::over_approximation(phic,zv);
      
#ifdef DEBUG
      std::cerr << "stepsize*jacobian=" << hdf << std::endl;
      std::cerr << "flow derivative=" << dphi << std::endl;

      std::cerr << "centre=" << c << std::endl;
      std::cerr << "rectangular centre=" << rc << std::endl;
      std::cerr << "flow bounds=" << phic << std::endl;
      std::cerr << "zonotopic <vector> to add to image of centre=" << zv << std::endl;
      std::cerr << "new approximation=" << p << std::endl;
      std::cerr << "Done integration step\n\n\n" << std::endl;
#endif  
      return p;
    }

    template<typename R>
    Geometry::Parallelotope<R> 
    C1LohnerIntegrator<R>::integration_step(const System::AffineVectorField<R>& vector_field, 
                                            const Geometry::Parallelotope<R>& initial_set, 
                                            R& step_size) const
    {
#ifdef DEBUG
      std::cerr << "integration_step(AffineVectorField<R>, Parallelotope<R>, R)\n";
#endif
      const System::AffineVectorField<R>& vf=vector_field;
      Geometry::Parallelotope<R> p=initial_set;
      R& h=step_size;
      
      R max_error=LinearAlgebra::norm(p.generators())/16;
      assert(max_error>0);
      
      if(debug_level>0) { 
        std::cerr << "parallelotope generators=" << p.generators() << std::endl;
        std::cerr << "maximum allowed error=" << max_error << std::endl;
      
        std::cerr << "jacobian=" << vf.A() << std::endl;
        std::cerr << "step size=" << h << std::endl;
      }
      

      LinearAlgebra::Matrix<R> D=LinearAlgebra::exp_Ah_approx(vf.A(),h,max_error);
      if(debug_level>0) { std::cerr << "approximate derivative=" << D << std::endl; }
      LinearAlgebra::Matrix<R> P=LinearAlgebra::exp_Ah_sub_id_div_A_approx(vf.A(),h,max_error);
      if(debug_level>0) { std::cerr << "twist=" << P << std::endl; }
      
      LinearAlgebra::Vector< Interval<R> > ib=vf.b();
       LinearAlgebra::Matrix< Interval<R> > iD=D;
      if(debug_level>0) { std::cerr << "approximating derivative=" << iD << std::endl; }
      LinearAlgebra::Matrix< Interval<R> > iP=P;
      if(debug_level>0) { std::cerr << "approximating twist=" << iP << std::endl; }
      //LinearAlgebra::Vector< Interval<R> > iC=iD*p.centre().position_vector()+iP*vf.b();
      LinearAlgebra::Vector< Interval<R> > iv1=iD*p.centre().position_vector();
      if(debug_level>0) { std::cerr << "iv1=" << iv1 << std::endl; }
       LinearAlgebra::Vector< Interval<R> > iv2=iP*ib;
      if(debug_level>0) { std::cerr << "iv2=" << iv2 << std::endl; }
      LinearAlgebra::Vector< Interval<R> > iC=iv1+iv2;
      
      if(debug_level>0) { std::cerr << "interval centre=" << iC << std::endl; }
      
      p=Geometry::Parallelotope<R>::over_approximation(Geometry::Rectangle<R>(iC),iD*p.generators());
      if(debug_level>0) { std::cerr << "parallelotope=" << p << std::endl; }
      //IntervalParallelotope<R> img(iD*p.centre().position_vector()+iP*vf.b(),iD*p.generators());
      
      return p;      
    }
    

template<typename R>
    Geometry::Zonotope<R> 
    C1LohnerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                            const Geometry::Zonotope<R>& initial_set, 
                                            R& step_size) const
    {
      const System::VectorField<R>* cvf_ptr=&vector_field;
      const System::AffineVectorField<R>* cavf_ptr=dynamic_cast< const System::AffineVectorField<R>* >(cvf_ptr);

#ifdef DEBUG
      if(cavf_ptr) { std::cerr << typeid(*cavf_ptr).name(); } else { std::cerr << "void"; } std::cerr << std::endl;
      std::cerr << cvf_ptr << "  " << cavf_ptr << "\n" << std::endl;
#endif
      if(cavf_ptr != NULL) {
        return integration_step(*cavf_ptr,initial_set,step_size);
      }

#ifdef DEBUG
      std::cerr << "integration_step(VectorField<R>, Zonotope<R>, R)" << std::endl<<std::flush;
#endif

      typedef typename numerical_traits<R>::field_extension_type F;
       
      using namespace LinearAlgebra;
      using namespace Geometry;
      using namespace System;

      assert(vector_field.dimension()==initial_set.dimension());

      const VectorField<R>& vf(vector_field);
      Zonotope<R> p=initial_set;
      const size_type n=p.dimension();
      R& h=step_size;
      const Matrix<R> id=identity_matrix<R>(n);
      
      R err=norm(p.generators())/65536;
      
      Rectangle<R> b=this->estimate_flow_bounds(vf,p.bounding_box(),h);
#ifdef DEBUG
      std::cerr << "suggested stepsize=" << step_size << std::endl;
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << b << std::endl;
#endif
      
      Vector< Interval<R> > f=vf(b);
      Matrix< Interval<R> > df=vf.derivative(b);
      
      R l=df.log_norm().upper();

#ifdef DEBUG
      std::cerr << "flow=" << f << std::endl;
      std::cerr << "jacobian=" << df << std::endl;
      std::cerr << "jacobian centre=" << df.centre() << std::endl;
      std::cerr << "logarithmic_norm=" << l << std::endl;
#endif

      l=max(l,R(1));
      
      //integration_step(): the varaible l is useless!
      
      Matrix< Interval<R> > hdf=h*df;
      Matrix< Interval<R> > dphi=exp(hdf);
      
      Point<R> c=p.centre();
      Rectangle<R> rc=Geometry::Rectangle<R>(c,c);
      Rectangle<R> phic=refine_flow_bounds(vf,rc,b,h);
      
      Matrix< Interval<R> > zv(dphi*p.generators());

      p=Geometry::Zonotope<R>::over_approximation(phic,zv);
      
#ifdef DEBUG
      std::cerr << "stepsize*jacobian=" << hdf << std::endl;
      std::cerr << "flow derivative=" << dphi << std::endl;

      std::cerr << "centre=" << c << std::endl;
      std::cerr << "bounds on centre=" << phic << std::endl;
      std::cerr << "zonotopic <vector> to add to image of centre=" << zv << std::endl;
      std::cerr << "new approximation=" << p << std::endl;
#endif  

      return p;
    }

    template<typename R>
    Geometry::Zonotope<R> 
    C1LohnerIntegrator<R>::integration_step(const System::AffineVectorField<R>& vector_field, 
                                            const Geometry::Zonotope<R>& initial_set, 
                                            R& step_size) const
    {
#ifdef DEBUG
      std::cerr << "integration_step(AffineVectorField<R>, Parallelotope<R>, R)\n";
#endif
      const System::AffineVectorField<R>& vf=vector_field;
      Geometry::Zonotope<R> p=initial_set;
      R& h=step_size;
      
      R max_error=LinearAlgebra::norm(p.generators())/16;
      assert(max_error>0);
      
#ifdef DEBUG
      std::cerr << "parallelotope generators=" << p.generators() << std::endl;
      std::cerr << "maximum allowed error=" << max_error << std::endl;
      
      std::cerr << "jacobian=" << vf.A() << std::endl;
      std::cerr << "step size=" << h << std::endl;
#endif

      LinearAlgebra::Matrix<R> D=LinearAlgebra::exp_Ah_approx(vf.A(),h,max_error);
      /* Write phi(x)=D x0 + P b */
      LinearAlgebra::Matrix<R> P=LinearAlgebra::exp_Ah_sub_id_div_A_approx(vf.A(),h,max_error);
      
      LinearAlgebra::Matrix< Interval<R> > iD=D;
      LinearAlgebra::Matrix< Interval<R> > iP=P;
      LinearAlgebra::Vector< Interval<R> > ib=vf.b();
      LinearAlgebra::Vector< Interval<R> > iC=iD*p.centre().position_vector()+iP*ib;
      p=Geometry::Zonotope<R>::over_approximation(Geometry::Rectangle<R>(iC),iD*p.generators());
      
#ifdef DEBUG
      std::cerr << "twist=" << P << std::endl;
      std::cerr << "approximate derivative=" << D << std::endl;
      
      std::cerr << "approximating derivative=" << iD << std::endl;
      std::cerr << "approximating twist=" << iP << std::endl;
      std::cerr << "bounds on centre=" << iC << std::endl;
      
      std::cerr << "parallelotope=" << p << std::endl;
#endif
      return p;      
    }

    
    template<typename R>
    Geometry::Zonotope<R> 
    C1LohnerIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                             const Geometry::Zonotope<R>& initial_set, 
                                             R& step_size) const
    {

#ifdef DEBUG
      std::cerr << "reachability_step(VectorField<R>, Zonotope<R>, R)\n";
#endif

      typedef typename numerical_traits<R>::field_extension_type F;

      using namespace LinearAlgebra;
      using namespace Geometry;
      using namespace System;

      assert(vector_field.dimension()==initial_set.dimension());

      const VectorField<R>& vf(vector_field);
      Zonotope<R> z=initial_set;
      const size_type n=z.dimension();
      R h=step_size;
      const Matrix<R> id=identity_matrix<R>(n);
      
      /* Throws exception if we can't find flow bounds for given stepsize. */
      Rectangle<R> b=estimate_flow_bounds(vf,z.bounding_box(),h,256);
      
      Vector< Interval<R> > f=vf(b);
      Matrix< Interval<R> > df=vf.derivative(b);
      
      Matrix< Interval<R> > dphi=id+Interval<R>(0,h)*df;
      
      Point<R> c=z.centre();
      Rectangle<R> rc=Geometry::Rectangle<R>(c,c);
      Rectangle<R> phic=refine_flow_bounds(vf,rc,b,R(h/2));
      
      Vector< Interval<R> > fh=(R(h/2)*f);
      
      Matrix<R> zfh=symmetrize(fh);
      Matrix<R> mdf=over_approximation(dphi)*z.generators();
      
      z=Zonotope<R>(phic,zfh)+mdf;
#ifdef DEBUG
      std::cerr << "suggested stepsize=" << step_size << std::endl;
        
      std::cerr << "stepsize=" << h << std::endl;
      std::cerr << "bound=" << b << std::endl;
      
      std::cerr << "flow=" << f << "=" << std::endl;
      std::cerr << "jacobian=" << df << std::endl;
      std::cerr << "flow derivative=" << dphi << std::endl;
        
      std::cerr << "centre=" << c << std::endl;
      std::cerr << "bounds on centre=" << phic << std::endl;
      
      std::cerr << "flow times stepsize=" << fh << std::endl;
      std::cerr << "symmetrised flow=" << zfh << std::endl;
      std::cerr << "over approximating Matrix=" << mdf << std::endl;
      
      std::cerr << "approximating zonotope " << z;
#endif

      return z;
    }

  }
}
