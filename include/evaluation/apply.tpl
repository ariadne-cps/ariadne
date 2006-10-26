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

   
    template<class R>
    Geometry::Rectangle<R> 
    C0Applicator<R>::apply(const System::Map<R>& f, const Geometry::Rectangle<R>& r) const
    {
      return f.image(r);
    }
    
    
    template<class R>
    Geometry::Parallelotope<R> 
    C1Applicator<R>::apply(const System::Map<R>& f, const Geometry::Parallelotope<R>& p) const 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;

      const size_type m=p.dimension();
      const size_type n=p.dimension();
      
      LinearAlgebra::Vector< Interval<R> > cuboid_vector(m);
      const Interval<R> unit_interval(-1,1);
      for(size_type i=0; i!=cuboid_vector.size(); ++i) {
        cuboid_vector(i)=Interval<R>(-1,1);
      }
            
      const Geometry::Point<R>& c=p.centre();
      const LinearAlgebra::Matrix<R>& g=p.generators();
      
      Geometry::Point< Interval<R> > img_centre=f(c);
      LinearAlgebra::Matrix< Interval<R> > df_on_set = f.jacobian(p.bounding_box());
      LinearAlgebra::Matrix< Interval<R> > df_at_centre = f.jacobian(c);
      
      LinearAlgebra::Matrix< Interval<R> > img_generators = df_at_centre*g;
      
      LinearAlgebra::Matrix< Interval<R> > img_generators_inverse = LinearAlgebra::inverse(LinearAlgebra::Matrix< Interval<R> >(img_generators));
      
      LinearAlgebra::Matrix< Interval<R> > img_generators_on_set = df_on_set * g;
      LinearAlgebra::Matrix< Interval<R> > cuboid_transform = img_generators_inverse * img_generators_on_set;
      
      LinearAlgebra::Vector< Interval<R> > new_cuboid = cuboid_transform * cuboid_vector;
      
      R new_cuboid_sup(0);
      for(size_type j=0; j!=n; ++j) {
        new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid(j).lower())) );
        new_cuboid_sup=std::max( new_cuboid_sup, R(abs(new_cuboid(j).upper())) );
      }
      
      // FIXME: This is incorrect; need over-approximations
      Geometry::Point<R> nc=Geometry::approximate_value(img_centre);
      LinearAlgebra::Matrix<R> ng=approximate_value(img_generators);
      
      Geometry::Parallelotope<R> result(nc,ng);
      return result;
    }

    template<class R, template<class> class BS>
    Applicator<R,BS>::~Applicator() 
    {
    }
    
    template<class R, template<class> class BS>
    Geometry::ListSet<R,BS> 
    Applicator<R,BS>::apply(const System::Map<R>& f, const Geometry::ListSet<R,BS>& ds) const 
    {
      Geometry::ListSet<R,BS> result(f.result_dimension());
      for(typename Geometry::ListSet<R,BS>::const_iterator iter=ds.begin(); iter!=ds.end(); ++iter) {
        result.push_back(this->apply(f,*iter));
      }
      return result;
    }
     
    
    template<class R, template<class> class BS>
    Geometry::GridMaskSet<R> 
    Applicator<R,BS>::apply(const System::Map<R>& f, 
                            const Geometry::GridMaskSet<R>& initial_set, 
                            const Geometry::GridMaskSet<R>& bounding_set) const 
    {
      typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
      assert(initial_set.bounded() && bounding_set.bounded());
      
      const Geometry::Grid<R>& g=initial_set.grid();
      Combinatoric::LatticeBlock bd=initial_set.block();
      Geometry::GridMaskSet<R> image(g,bd);
      Geometry::Rectangle<R> bb=initial_set.bounding_box();
      
      for(gms_const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
        Geometry::Rectangle<R> r=*iter;
        BS<R> p=convert_to< BS<R> >(r);
        BS<R> fp=this->apply(f,p);
        image.adjoin(over_approximation(fp,g));
      }
      return regular_intersection(image,bounding_set);
    }
    
    
    template<class R, template<class> class BS>
    Geometry::GridMaskSet<R> 
    Applicator<R,BS>::chainreach(const System::Map<R>& f, 
                                 const Geometry::GridMaskSet<R>& initial_set, 
                                 const Geometry::GridMaskSet<R>& bounding_set) const
    {
      typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
      assert(initial_set.bounded() && bounding_set.bounded());
      
      const Geometry::Grid<R>& g=bounding_set.grid();
      Combinatoric::LatticeBlock bd=bounding_set.block();
      Geometry::GridBlock<R> bb(g,bd);
      Geometry::GridMaskSet<R> result(g,bd);
      Geometry::GridCellListSet<R> found(g);
      Geometry::GridCellListSet<R> image(g);
      
      //std::cerr << "result.size() = " << result.size() << "\n";
      found=initial_set;
      //std::cerr << "found.size() = " << found.size() << "\n";
      while(!subset(found,result)) {
        found=difference(found,result);
        //std::cerr << "new.size() = " << found.size() << "\n";
        result.adjoin(found);
        //std::cerr << "result.size() = " << result.size() << "\n";
        image.clear(); 
        uint size=0;
        for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
          ++size;
          Geometry::Rectangle<R> r=*iter;
          BS<R> p=convert_to< BS<R> >(r);
          BS<R> fp=this->apply(f,p);
          if(!disjoint(fp.bounding_box(),bounding_set)) {
            image.adjoin(over_approximation(fp,g));
          }
          //std::cerr << size << "\n";
        }
        //std::cerr << "image.size() = " << image.size() << "\n";
        found=regular_intersection(image,bounding_set);
        //std::cerr << "found.size() = " << found.size() << "\n";
      }
      //std::cerr << "result.size() = " << result.size() << "\n";
      //std::cerr << "result.bounded() = " << result.bounded() << "\n";
      return result;
    }
  
    template<class R, template<class> class BS>
    bool
    Applicator<R,BS>::verify(const System::Map<R>& f, 
                             const Geometry::GridMaskSet<R>& initial_set, 
                             const Geometry::GridMaskSet<R>& safe_set) const
    {
      typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
      assert(initial_set.bounded() && safe_set.bounded());
      
      const Geometry::Grid<R>& g=initial_set.grid();
      Combinatoric::LatticeBlock bd=safe_set.block();
      Geometry::Rectangle<R> bb=safe_set.bounding_box();
      Geometry::GridMaskSet<R> reach(g,bd);
      Geometry::GridCellListSet<R> found(g);
      Geometry::GridCellListSet<R> cell_image(g);
      Geometry::GridCellListSet<R> image(g);
      
      found=initial_set;
      while(!subset(found,reach)) {
        found=difference(found,reach);
        reach.adjoin(found);
        image.clear();
        uint size=0;
        for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
          ++size;
          Geometry::Rectangle<R> r=*iter;
          BS<R> p=convert_to< BS<R> >(r);
          BS<R> fp=this->apply(f,p);
          cell_image.adjoin(over_approximation(fp,g));
          if(!subset(cell_image,safe_set)) {
            return false;
          }
          image.adjoin(cell_image);
          cell_image.clear();
        }
        found=image;
      }
      return true;
    }

  }
}
