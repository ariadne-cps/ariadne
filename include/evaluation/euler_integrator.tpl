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

#include "euler_integrator.h"

#include "../base/array.h"
#include "../numeric/arithmetic.h"
#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../geometry/rectangle.h"

#include "../system/vector_field.h"

#include "../evaluation/integrator.h"

namespace Ariadne {
  namespace Evaluation {
    
#ifdef DEBUG
    static const int verbosity=1;
#else
    static const int verbosity=0;
#endif
    
    template<class R>
    EulerIntegrator<R>::EulerIntegrator(const time_type& maximum_step_size, const time_type& lock_to_grid_time, const R& maximum_basic_set_radius)
      : Integrator<R>(maximum_step_size,lock_to_grid_time,maximum_basic_set_radius)
    {
    }

    template<class R>
    Geometry::Rectangle<R> 
    EulerIntegrator<R>::integration_step(const System::VectorField<R>& vector_field, 
                                            const Geometry::Rectangle<R>& initial_set, 
                                            time_type& step_size) const
    {
      if(verbosity>0) { 
        std::cerr << "integrate(const VectorField<R>& vf, const Geometry::Rectangle<R>& r, const time_type& t)" << std::endl;
      }

      check_dimension(vector_field,initial_set,__PRETTY_FUNCTION__);
      
      const System::VectorField<R>& vf(vector_field);
      Geometry::Rectangle<R> r=initial_set;
      Geometry::Rectangle<R> q=estimate_flow_bounds(vf,r,step_size);
      
      Interval<R> h=step_size;      
      LinearAlgebra::Vector< Interval<R> > fq=vf(q);
      r=r+(h*fq);

      if(verbosity>0) {
        std::cerr << "suggested stepsize=" << step_size << std::endl;
                
        std::cerr << "stepsize=" << h << std::endl;
        std::cerr << "bound=" << q << std::endl;

        std::cerr << "derivative=" << fq << std::endl;

        std::cerr << "position=" << r << std::endl;
      }
      return r;
    }



    template<class R>
    Geometry::Rectangle<R> 
    EulerIntegrator<R>::reachability_step(const System::VectorField<R>& vector_field, 
                                           const Geometry::Rectangle<R>& initial_set, 
                                           time_type& step_size) const
    {
#ifdef DEBUG
      std::cerr << "integrate(const VectorField<R>& vf, const Geometry::Rectangle<R>& r, const time_type& t)" << std::endl;
#endif
      check_dimension(vector_field,initial_set(),__PRETTY_FUNCTION__);
      
      const System::VectorField<R>& vf(vector_field);
      Geometry::Rectangle<R> r=initial_set;
      time_type& h=step_size;

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
    
  }
}
