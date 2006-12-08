/***************************************************************************
 *            applicator.tpl
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


#include "../system/map.h"


namespace Ariadne {
  namespace Evaluation {

   
    template<class R>
    Applicator<R>::Applicator() 
    {
    }
    


    template<class R>
    Applicator<R>::~Applicator() 
    {
    }
    



    template<class R>
    Geometry::Rectangle<R> 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::Rectangle<R>& r) const
    {
      return f.image(r);
    }
    
    
    
    template<class R>
    Geometry::Parallelotope<R> 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::Parallelotope<R>& p) const 
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
    Geometry::Parallelotope< Interval<R> > 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::Parallelotope< Interval<R> >& p) const 
    {
      typedef Interval<R> I;
      const Geometry::Point<I>& c=p.centre();
      //Geometry::Point<I> img_centre=f(c);
      
      Geometry::Point<Interval<R> > img_centre(f(c));
      LinearAlgebra::Matrix<I> df_on_set = f.jacobian(p.bounding_box());
      LinearAlgebra::Matrix<I> img_generators = df_on_set*p.generators();

      Geometry::Parallelotope<I> result(img_centre,img_generators);
      return result;
    }

    
    
    template<class R>
    Geometry::Zonotope< Interval<R> > 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::Zonotope< Interval<R> >& z) const 
    {
      typedef Interval<R> I;

      Geometry::Point<I> img_centre=f(z.centre());
      LinearAlgebra::Matrix<I> df_on_set = f.jacobian(z.bounding_box());
      LinearAlgebra::Matrix<I> img_generators = df_on_set*z.generators();

      Geometry::Zonotope<I> result(img_centre,img_generators);
      return result;
    }

    template<class R>
    template<class Rl,template<class> class BS>
    Geometry::ListSet<Rl,BS> 
    Applicator<R>::image_list_set(const System::Map<R>& f, const Geometry::ListSet<Rl,BS>& ds) const 
    {
      Geometry::ListSet<Rl,BS> result(f.result_dimension());
      for(typename Geometry::ListSet<Rl,BS>::const_iterator iter=ds.begin(); iter!=ds.end(); ++iter) {
        result.push_back(this->image(f,*iter));
      }
      return result;
    }
     


    template<class R>
    Geometry::ListSet<R,Geometry::Rectangle> 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::ListSet<R,Geometry::Rectangle>& ds) const 
    {
      return this->image_list_set(f,ds);
    }
     
    
    template<class R>
    Geometry::ListSet<R,Geometry::Parallelotope> 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::ListSet<R,Geometry::Parallelotope>& ds) const 
    {
      return this->image_list_set(f,ds);
    }
     
    
    template<class R>
    Geometry::ListSet<R,Geometry::Zonotope> 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::ListSet<R,Geometry::Zonotope>& ds) const 
    {
      return this->image_list_set(f,ds);
    }
     
    template<class R>
    Geometry::ListSet<Interval<R>,Geometry::Parallelotope> 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::ListSet<Interval<R>,Geometry::Parallelotope>& ds) const 
    {
      return this->image_list_set(f,ds);
    }
     
    template<class R>
    Geometry::ListSet<Interval<R>,Geometry::Zonotope> 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::ListSet<Interval<R>,Geometry::Zonotope>& ds) const 
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
      throw NotImplemented(__PRETTY_FUNCTION__);
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
    bool
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

  }
}
