/***************************************************************************
 *            applicator.code.h
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
 
#include "applicator.h"

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "../exceptions.h"

#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../combinatoric/lattice_set.h"

#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_set.h"
#include "../system/grid_multimap.h"


#include "../system/map.h"


namespace Ariadne {
  namespace Evaluation {

   
    template<class R>
    Applicator<R>::Applicator() 
      : _maximum_basic_set_radius(0.1)
    {
    }
    


    template<class R>
    Applicator<R>::~Applicator() 
    {
    }
    



    template<class R>
    R
    Applicator<R>::maximum_basic_set_radius() const
    {
      return this->_maximum_basic_set_radius;
    }
  

    template<class R>
    void
    Applicator<R>::set_maximum_basic_set_radius(const R& mbsr) 
    {
      if(mbsr<=0) {
        throw std::runtime_error("maximum basic set radius must be positive");
      }
      this->_maximum_basic_set_radius=mbsr;
    }
  

    template<class R>
    Geometry::Rectangle<R> 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::Rectangle<R>& r) const
    {
      return Geometry::Rectangle<R>(f.image(Geometry::Point< Interval<R> >(r)));
    }
    
    
    


    template<class R>
    Geometry::Zonotope<R> 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::Zonotope<R>& z) const 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;

      const size_type m=z.dimension();
      const size_type n=z.dimension();
      
      LinearAlgebra::Vector< Interval<R> > cuboid_vector(m);
      const Interval<R> unit_interval(-1,1);
      for(size_type i=0; i!=cuboid_vector.size(); ++i) {
        cuboid_vector(i)=Interval<R>(-1,1);
      }
            
      const Geometry::Point<R>& c=z.centre();
      const LinearAlgebra::Matrix<R>& g=z.generators();
      
      Geometry::Point< Interval<R> > img_centre=f(c);
      LinearAlgebra::Matrix< Interval<R> > df_on_set = f.jacobian(z.bounding_box());
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
      
      Geometry::Zonotope<R> result(nc,ng);
      return result;
    }


    
    
    
    template<class R>
    Geometry::Zonotope< Interval<R> > 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::Zonotope< Interval<R> >& z) const 
    {
      typedef Interval<R> I;

      Geometry::Point<I> img_centre=f(z.centre());
      LinearAlgebra::Matrix<I> df_on_set = f.jacobian(over_approximation(z.bounding_box()));
      LinearAlgebra::Matrix<I> img_generators = df_on_set*z.generators();

      Geometry::Zonotope<I> result(img_centre,img_generators);
      return result;
    }

    template<class R>
    template<class BS>
    Geometry::ListSet<BS> 
    Applicator<R>::image_list_set(const System::Map<R>& f, const Geometry::ListSet<BS>& ds) const 
    {
      Geometry::ListSet<BS> result(f.result_dimension());
      for(typename Geometry::ListSet<BS>::const_iterator iter=ds.begin(); iter!=ds.end(); ++iter) {
        result.push_back(this->image(f,*iter));
      }
      return result;
    }
     


    template<class R>
    Geometry::ListSet< Geometry::Rectangle<R> > 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::ListSet< Geometry::Rectangle<R> >& ds) const 
    {
      return this->image_list_set(f,ds);
    }
     
    
    
    template<class R>
    Geometry::ListSet< Geometry::Zonotope<R> > 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::ListSet< Geometry::Zonotope<R> >& ds) const 
    {
      return this->image_list_set(f,ds);
    }
     

    template<class R>
    Geometry::ListSet< Geometry::Zonotope< Interval<R> > > 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::ListSet< Geometry::Zonotope< Interval<R> > >& ds) const 
    {
      return this->image_list_set(f,ds);
    }
     
    
    template<class R>
    Geometry::GridCellListSet<R> 
    Applicator<R>::image(const System::Map<R>& f, 
                         const Geometry::GridCellListSet<R>& initial_set, 
                         const Geometry::Grid<R>& image_grid) const 
    {
      typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
      
      Geometry::GridCellListSet<R> image(image_grid);

      Geometry::Rectangle<R> r(f.argument_dimension());
      Geometry::Rectangle<R> fr(f.result_dimension());
      Geometry::Parallelotope<R> p(f.argument_dimension());
      Geometry::Parallelotope<R> fp(f.result_dimension());
      
      for(gcls_const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
        Geometry::Rectangle<R> r=*iter;
        if(f.smoothness()>=1) {
          p=r;
          fp=this->image(f,p);
          image.adjoin(over_approximation(fp,image_grid));
        } else {
          fr=this->image(f,r);
          image.adjoin(over_approximation(fr,image_grid));
        }
      }
      return image;
    }
    
    
    
    template<class R>
    Geometry::GridMaskSet<R> 
    Applicator<R>::image(const System::Map<R>& f, 
                         const Geometry::GridMaskSet<R>& initial_set, 
                         const Geometry::FiniteGrid<R>& image_grid) const 
    {
      typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
      check_bounded(initial_set,"Applicator<R>::image(...)");
      
      Geometry::GridMaskSet<R> image(image_grid);
      const Geometry::Grid<R>& g=image_grid.grid();
      
      Geometry::Rectangle<R> r(f.argument_dimension());
      Geometry::Rectangle<R> fr(f.result_dimension());
      Geometry::Parallelotope<R> p(f.argument_dimension());
      Geometry::Parallelotope<R> fp(f.result_dimension());
      
      for(gms_const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
        Geometry::Rectangle<R> r=*iter;
        if(f.smoothness()>=1) {
          p=r;
          fp=this->image(f,p);
          image.adjoin(over_approximation(fp,g));
        } else {
          fr=this->image(f,r);
          image.adjoin(over_approximation(fr,g));
        }
      }
      return image;
    }
    
    
    
    template<class R>
    Geometry::GridMaskSet<R> 
    Applicator<R>::image(const System::Map<R>& f, 
                         const Geometry::GridMaskSet<R>& initial_set, 
                         const Geometry::GridMaskSet<R>& bounding_set) const 
    {
      typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
      check_bounded(initial_set,"Applicator<R>::image(...)");
      check_bounded(bounding_set,"Applicator<R>::image(...)");
      
      const Geometry::Grid<R>& g=initial_set.grid();
      Combinatoric::LatticeBlock bd=initial_set.block();
      Geometry::GridMaskSet<R> image(g,bd);
      Geometry::Rectangle<R> bb=initial_set.bounding_box();
      
      Geometry::Rectangle<R> r(g.dimension());
      Geometry::Rectangle<R> fr(g.dimension());
      Geometry::Parallelotope<R> p(g.dimension());
      Geometry::Parallelotope<R> fp(g.dimension());
      
      for(gms_const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
        Geometry::Rectangle<R> r=*iter;
        if(f.smoothness()>=1) {
          p=r;
          fp=this->image(f,p);
          image.adjoin(over_approximation(fp,g));
        } else {
          fr=this->image(f,r);
          image.adjoin(over_approximation(fr,g));
        }
      }
      return regular_intersection(image,bounding_set);
    }
    
    
    
    template<class R>
    Geometry::GridMaskSet<R> 
    Applicator<R>::reach(const System::Map<R>& f, 
                         const Geometry::GridMaskSet<R>& initial_set) const 
    {
      Geometry::GridMaskSet<R> result(initial_set.grid(),initial_set.block());
      typedef typename Geometry::GridMaskSet<R>::const_iterator basic_set_iterator;
      Geometry::Rectangle<R> rectangle(initial_set.dimension());
      Geometry::Zonotope< Interval<R> > zonotope(result.dimension());
      for(basic_set_iterator bs_iter=initial_set.begin(); 
          bs_iter!=initial_set.end(); ++bs_iter)
      {
        rectangle=*bs_iter;
        zonotope=rectangle;
        while(zonotope.radius() < this->maximum_basic_set_radius()) {
          result.adjoin_over_approximation(zonotope);
          zonotope=this->image(f,zonotope);
        }
      }
      return result;
    }
    
    template<class R>
    Geometry::GridMaskSet<R> 
    Applicator<R>::chainreach(const System::Map<R>& f, 
                              const Geometry::GridMaskSet<R>& initial_set, 
                              const Geometry::GridMaskSet<R>& bounding_set) const
    {
      typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
      check_bounded(initial_set,"Applicator<R>::chainreach(...)");
      check_bounded(bounding_set,"Applicator<R>::chainreach(...)");
      
      const Geometry::Grid<R>& g=bounding_set.grid();
      Combinatoric::LatticeBlock bd=bounding_set.block();
      Geometry::GridBlock<R> bb(g,bd);
      Geometry::GridMaskSet<R> result(g,bd);
      Geometry::GridCellListSet<R> found(g);
      Geometry::GridCellListSet<R> image(g);
      
      Geometry::Rectangle<R> r(g.dimension());
      Geometry::Rectangle<R> fr(g.dimension());
      Geometry::Parallelotope<R> p(g.dimension());
      Geometry::Parallelotope<R> fp(g.dimension());
      
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
          if(f.smoothness()>0) {
            p=r;
            fp=this->image(f,p);
            if(!disjoint(fp.bounding_box(),bounding_set)) {
              image.adjoin(over_approximation(fp,g));
            }
          } else {
            fr=this->image(f,r);
            if(!disjoint(fr.bounding_box(),bounding_set)) {
              image.adjoin(over_approximation(fr,g));
            }
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
  
    
    
    template<class R>
    tribool
    Applicator<R>::verify(const System::Map<R>& f, 
                          const Geometry::GridMaskSet<R>& initial_set, 
                          const Geometry::GridMaskSet<R>& safe_set) const
    {
      typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
      check_bounded(initial_set,"Applicator<R>::verify(...)");
      check_bounded(safe_set,"Applicator<R>::verify(...)");
      
      const Geometry::Grid<R>& g=initial_set.grid();
      Combinatoric::LatticeBlock bd=safe_set.block();
      Geometry::Rectangle<R> bb=safe_set.bounding_box();
      Geometry::GridMaskSet<R> reach(g,bd);
      Geometry::GridCellListSet<R> found(g);
      Geometry::GridCellListSet<R> cell_image(g);
      Geometry::GridCellListSet<R> image(g);
      
      Geometry::Rectangle<R> r(g.dimension());
      Geometry::Rectangle<R> fr(g.dimension());
      Geometry::Parallelotope<R> p(g.dimension());
      Geometry::Parallelotope<R> fp(g.dimension());
      
      found=initial_set;
      while(!subset(found,reach)) {
        found=difference(found,reach);
        reach.adjoin(found);
        image.clear();
        uint size=0;
        for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
          ++size;
          Geometry::Rectangle<R> r=*iter;
          if(f.smoothness()>=1) {
            p=r;
            fp=this->image(f,p);
            cell_image.adjoin(over_approximation(fp,g));
            if(!subset(cell_image,safe_set)) {
              return false;
            }
          } else {
            fr=this->image(f,r);
            cell_image.adjoin(over_approximation(fr,g));
            if(!subset(cell_image,safe_set)) {
              return false;
            }
          }
          image.adjoin(cell_image);
          cell_image.clear();
        }
        found=image;
      }
      return true;
    }

    template<class R>
    System::GridMultiMap<R> 
    Applicator<R>::discretize(const System::Map<R>& f, 
                              const Geometry::GridMaskSet<R>& domain,
                              const Geometry::Grid<R>& range_grid) const
    {
      System::GridMultiMap<R> result(domain.grid(),range_grid);
      Geometry::Zonotope<R> basic_set;
      Geometry::Zonotope<R> image_set;
      for(typename Geometry::GridMaskSet<R>::const_iterator dom_iter=domain.begin();
          dom_iter!=domain.end(); ++dom_iter)
      {
        const Geometry::GridCell<R>& cell=*dom_iter;
        basic_set=cell;
        image_set=this->image(f,basic_set);
        result.adjoin_to_image(cell,over_approximation(image_set,range_grid));
      }
      return result;
    }

  }
}
