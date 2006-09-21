/***************************************************************************
 *            apply.tpl
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
 
#include "apply.h"

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../combinatoric/lattice_set.h"

#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_set.h"


#include "../system/map.h"


namespace Ariadne {
  namespace Evaluation {

   
    template<typename R>
    Geometry::Rectangle<R> 
    C0Applicator<R>::apply(const System::Map<R>& f, const Geometry::Rectangle<R>& r) const
    {
      return f(r);
    }
    
    
    template<typename R>
    Geometry::Parallelotope<R> 
    C1Applicator<R>::apply(const System::Map<R>& f, const Geometry::Parallelotope<R>& p) const 
    {
      typedef typename numerical_traits<R>::field_extension_type F;

      const size_type m=p.dimension();
      const size_type n=p.dimension();
      
      LinearAlgebra::Vector< Interval<R> > cuboid_vector(m);
      const Interval<R> unit_interval(-1,1);
      for(size_type i=0; i!=cuboid_vector.size(); ++i) {
        cuboid_vector(i)=Interval<R>(-1,1);
      }
            
      const Geometry::Point<R>& c=p.centre();
      const LinearAlgebra::Matrix<R>& g=p.generators();
      
      Geometry::Point<R> img_centre=f(c);
      LinearAlgebra::Matrix< Interval<R> > df_on_set = f.derivative(p.bounding_box());
      LinearAlgebra::Matrix<R> df_at_centre = f.derivative(c);
      
      LinearAlgebra::Matrix<R> img_generators = df_at_centre*g;
      
      LinearAlgebra::Matrix< Interval<R> > img_generators_inverse = LinearAlgebra::inverse(LinearAlgebra::Matrix< Interval<R> >(img_generators));
      
      LinearAlgebra::Matrix< Interval<R> > img_generators_on_set = df_on_set * g;
      LinearAlgebra::Matrix< Interval<R> > cuboid_transform = img_generators_inverse * img_generators_on_set;
      
      LinearAlgebra::Vector< Interval<R> > new_cuboid = cuboid_transform * cuboid_vector;
      
      R new_cuboid_sup(0);
      for(size_type j=0; j!=n; ++j) {
        new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid(j).lower())) );
        new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid(j).upper())) );
      }
      
      Geometry::Parallelotope<R> result(img_centre,img_generators);
      return result;
    }

    template<typename R, template<typename> class BS>
    Applicator<R,BS>::~Applicator() 
    {
    }
    
    template<typename R, template<typename> class BS>
    Geometry::ListSet<R,BS> 
    Applicator<R,BS>::apply(const System::Map<R>& f, const Geometry::ListSet<R,BS>& ds) const 
    {
      Geometry::ListSet<R,BS> result(f.result_dimension());
      for(typename Geometry::ListSet<R,BS>::const_iterator iter=ds.begin(); iter!=ds.end(); ++iter) {
        result.push_back(this->apply(f,*iter));
      }
      return result;
    }
     
    
    template<typename R, template<typename> class BS>
    Geometry::GridMaskSet<R> 
    Applicator<R,BS>::apply(const System::Map<R>& f, 
                            const Geometry::GridMaskSet<R>& is, 
                            const Geometry::GridMaskSet<R>& bs) const 
    {
      typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
      
      const Geometry::Grid<R>& g=is.grid();
      Combinatoric::LatticeBlock bd=is.block();
      Geometry::GridMaskSet<R> image(g,bd);
      Geometry::Rectangle<R> bb=is.bounding_box();
      
      for(gms_const_iterator iter=is.begin(); iter!=is.end(); ++iter) {
        Geometry::Rectangle<R> r=*iter;
        BS<R> p=convert_to< BS<R> >(r);
        BS<R> fp=this->apply(f,p);
        image.adjoin(over_approximation_of_intersection(fp,bb,g));
      }
      return regular_intersection(image,bs);
    }
    
    
    template<typename R, template<typename> class BS>
    Geometry::GridMaskSet<R> 
    Applicator<R,BS>::chainreach(const System::Map<R>& f, 
                                 const Geometry::GridMaskSet<R>& is, 
                                 const Geometry::GridMaskSet<R>& bs) const
    {
      typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
      
      const Geometry::Grid<R>& g=is.grid();
      Combinatoric::LatticeBlock bd=is.block();
      Geometry::Rectangle<R> bb=is.bounding_box();
      Geometry::GridMaskSet<R> result(g,bd);
      Geometry::GridMaskSet<R> found(g,bd);
      Geometry::GridMaskSet<R> image(g,bd);
      
      found=is;
      while(!subset(found,result)) {
        found=difference(found,result);
        result.adjoin(found);
        image=Geometry::GridMaskSet<R>(g,bd);
        uint size=0;
        for(gms_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
          ++size;
          Geometry::Rectangle<R> r=*iter;
          BS<R> p=convert_to< BS<R> >(r);
          BS<R> fp=this->apply(f,p);
          image.adjoin(over_approximation_of_intersection(fp,bb,g));
        }
        found=regular_intersection(image,bs);
      }
      return result;
    }

  }
}
