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

#include "../numeric/interval.h"

#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

#include "../combinatoric/lattice_set.h"

#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid.h"
#include "../geometry/grid_set.h"

#include "../system/grid_multimap.h"


#include "../system/map.h"
#include "../system/discrete_time_system.h"


namespace Ariadne {
  namespace Evaluation {
  

    template<class R>
    Applicator<R>::Applicator() 
      : _maximum_basic_set_radius(0.25),
        _grid_size(0.125)
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
    R
    Applicator<R>::grid_size() const
    {
      return this->_grid_size;
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
    void
    Applicator<R>::set_grid_size(const R& mgs) 
    {
      if(mgs<=0) {
        throw std::runtime_error("maximum basic set radius must be positive");
      }
      this->_grid_size=mgs;
    }
  

    template<class R>
    Geometry::Rectangle<R> 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::Rectangle<R>& r) const
    {
      return Geometry::Rectangle<R>(f.image(Geometry::Point< Numeric::Interval<R> >(r)));
    }
    
    
    


    template<class R>
    Geometry::Zonotope<R> 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::Zonotope<R>& z) const 
    {
      typedef typename Numeric::traits<R>::arithmetic_type F;

      const size_type m=z.dimension();
      const size_type n=z.dimension();
      
      LinearAlgebra::Vector< Numeric::Interval<R> > cuboid_vector(m);
      const Numeric::Interval<R> unit_interval(-1,1);
      for(size_type i=0; i!=cuboid_vector.size(); ++i) {
        cuboid_vector(i)=Numeric::Interval<R>(-1,1);
      }
            
      const Geometry::Point<R>& c=z.centre();
      const LinearAlgebra::Matrix<R>& g=z.generators();
      
      Geometry::Point< Numeric::Interval<R> > img_centre=f(c);
      LinearAlgebra::Matrix< Numeric::Interval<R> > df_on_set = f.jacobian(z.bounding_box());
      LinearAlgebra::Matrix< Numeric::Interval<R> > df_at_centre = f.jacobian(c);
      
      LinearAlgebra::Matrix< Numeric::Interval<R> > img_generators = df_at_centre*g;
      
      LinearAlgebra::Matrix< Numeric::Interval<R> > img_generators_inverse = LinearAlgebra::inverse(LinearAlgebra::Matrix< Numeric::Interval<R> >(img_generators));
      
      LinearAlgebra::Matrix< Numeric::Interval<R> > img_generators_on_set = df_on_set * g;
      LinearAlgebra::Matrix< Numeric::Interval<R> > cuboid_transform = img_generators_inverse * img_generators_on_set;
      
      LinearAlgebra::Vector< Numeric::Interval<R> > new_cuboid = cuboid_transform * cuboid_vector;
      
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
    Geometry::Zonotope< Numeric::Interval<R> > 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::Zonotope< Numeric::Interval<R> >& z) const 
    {
      typedef Numeric::Interval<R> I;

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
    Geometry::ListSet< Geometry::Zonotope< Numeric::Interval<R> > > 
    Applicator<R>::image(const System::Map<R>& f, const Geometry::ListSet< Geometry::Zonotope< Numeric::Interval<R> > >& ds) const 
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
      ARIADNE_CHECK_BOUNDED(initial_set,"Applicator<R>::image(...)");
      
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
      ARIADNE_CHECK_BOUNDED(initial_set,"Applicator<R>::image(...)");
      ARIADNE_CHECK_BOUNDED(bounding_set,"Applicator<R>::image(...)");
      
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
      Geometry::Zonotope< Numeric::Interval<R> > zonotope(result.dimension());
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
      ARIADNE_CHECK_BOUNDED(initial_set,"Applicator<R>::chainreach(...)");
      ARIADNE_CHECK_BOUNDED(bounding_set,"Applicator<R>::chainreach(...)");
      
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
    Geometry::GridMaskSet<R>
    Applicator<R>::viable(const System::Map<R>& f, 
                          const Geometry::GridMaskSet<R>& bounding_set) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    

    template<class R>
    tribool
    Applicator<R>::verify(const System::Map<R>& f, 
                          const Geometry::GridMaskSet<R>& initial_set, 
                          const Geometry::GridMaskSet<R>& safe_set) const
    {
      typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
      ARIADNE_CHECK_BOUNDED(initial_set,"Applicator<R>::verify(...)");
      ARIADNE_CHECK_BOUNDED(safe_set,"Applicator<R>::verify(...)");
      
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
    Geometry::SetInterface<R>*
    Applicator<R>::image(const System::Map<R>& f, 
                         const Geometry::SetInterface<R>& set) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
 

    template<class R>
    Geometry::SetInterface<R>*
    Applicator<R>::preimage(const System::Map<R>& f, 
                         const Geometry::SetInterface<R>& set) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    
    template<class R>
    Geometry::SetInterface<R>*
    Applicator<R>::reach(const System::Map<R>& f, 
                         const Geometry::SetInterface<R>& initial_set) const
    {
      using namespace Geometry;
      Rectangle<R> r(f.argument_dimension());
      Rectangle<R> fr(f.result_dimension());
      Parallelotope<R> p(f.argument_dimension());
      Parallelotope<R> fp(f.result_dimension());
      
      Rectangle<R> bb;
      try {
        bb=initial_set.bounding_box();
      }
      catch(UnboundedSet) {
        throw std::runtime_error("Applicator::reach(Map,Set) only implemented for bounded initial sets");
      }
      
      ListSet< Parallelotope<R> > is;
      if(dynamic_cast<const Rectangle<R>*>(&initial_set)) {
        is.adjoin(static_cast< Parallelotope<R> >(dynamic_cast<const Rectangle<R>&>(initial_set)));
      } else if(dynamic_cast<const ListSet< Parallelotope<R> >*>(&initial_set)) {
        is=dynamic_cast<const ListSet< Parallelotope<R> >&>(initial_set);
      } else { 
        GridMaskSet<R> gms(Grid<R>(initial_set.dimension(),this->grid_size()),initial_set.bounding_box());
        gms.adjoin_under_approximation(initial_set);
        is=ListSet< Parallelotope<R> >(gms);
      }
      
      
      
      if(f.smoothness()>=1) {
        ListSet< Parallelotope<R> >* reach=new ListSet< Parallelotope<R> >;
        for(typename ListSet< Parallelotope<R> >::const_iterator iter=is.begin(); iter!=is.end(); ++iter) {
          p=*iter;
          //std::cerr << p << std::endl << p.radius() << " of " << this->maximum_basic_set_radius() << std::endl;
          while(p.radius()<this->maximum_basic_set_radius()) {
            reach->adjoin(p);
            p=this->image(f,p);
            //std::cerr << p << std::endl << p.radius() << " of " << this->maximum_basic_set_radius() << std::endl;
          }
        }
        return reach;
      } else {
        ListSet< Rectangle<R> >* reach=new ListSet< Rectangle<R> >;
        for(typename ListSet< Parallelotope<R> >::const_iterator iter=is.begin(); iter!=is.end(); ++iter) {
          r=iter->bounding_box();
          while(r.radius()<this->maximum_basic_set_radius()) {
            reach->adjoin(r);
            r=this->image(f,r);
          }
        }
        return reach;
      }
    }
    


    template<class R>
    Geometry::SetInterface<R>*
    Applicator<R>::chainreach(const System::Map<R>& f, 
                              const Geometry::SetInterface<R>& initial_set, 
                              const Geometry::SetInterface<R>& bounding_set) const
    {
      using namespace Geometry;
      
      Rectangle<R> bb;
      try {
        bb=bounding_set.bounding_box();
      }
      catch(UnboundedSet&) {
        throw UnboundedSet("chainreach(Map,Set,Set): bounding_set unbounded");
      }
      Grid<R> g(bounding_set.dimension(),this->grid_size());
      FiniteGrid<R> fg(g,bb);
      GridMaskSet<R> gbs(fg);
      gbs.adjoin_over_approximation(bounding_set);
      GridMaskSet<R> gis(fg);
      gis.adjoin_over_approximation(initial_set);
      
      GridMaskSet<R> gcrs=this->chainreach(f,gis,gbs);
      return new GridMaskSet<R>(gcrs);
    }
      
      


    template<class R>
    Geometry::SetInterface<R>*
    Applicator<R>::viable(const System::Map<R>& f, 
                          const Geometry::SetInterface<R>& bounding_set) const
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
    }
    


    template<class R>
    tribool
    Applicator<R>::verify(const System::Map<R>& f, 
                          const Geometry::SetInterface<R>& initial_set, 
                          const Geometry::SetInterface<R>& safe_set) const
    {
      using namespace Geometry;
      Rectangle<R> bb;
      try {
        bb=safe_set.bounding_box();
      }
      catch(UnboundedSet&) {
        throw UnboundedSet("verify(Map,Set,Set): safe_set unbounded");
      }
      if(!initial_set.subset(bb)) {
        return false;
      }
      Grid<R> g(safe_set.dimension(),this->grid_size());
      FiniteGrid<R> fg(g,bb);
      GridMaskSet<R> gss(fg);
      gss.adjoin_over_approximation(safe_set);
      GridMaskSet<R> gis(fg);
      gis.adjoin_over_approximation(initial_set);
      
      return this->verify(f,gis,gss);
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



    template<class R>
    System::GridMultiMap<R> 
    Applicator<R>::control_synthesis(const System::DiscreteTimeSystem<R>& f, 
                                     const Geometry::SetInterface<R>& initial_set,
                                     const Geometry::SetInterface<R>& target_set,
                                     const Geometry::GridMaskSet<R>& state_bounding_set,
                                     const Geometry::GridMaskSet<R>& input_bounding_set,
                                     const Geometry::GridMaskSet<R>& noise_bounding_set) const
    {
      // TODO: Use on-the-fly discretization
      // TODO: Use adaptive grid refinement
      
      typedef typename Numeric::traits<R>::interval_type I;

      using namespace Combinatoric;
      using namespace Geometry;
      using namespace System;

      const Grid<R>& state_grid = state_bounding_set.grid();
      const Grid<R>& input_grid = input_bounding_set.grid();
      dimension_type state_space_dimension = f.state_space_dimension();
      dimension_type input_space_dimension = f.control_space_dimension();

      // Discretize function
      std::map< std::pair<LatticeCell,LatticeCell>, LatticeCellListSet > discretization;
      for(typename GridMaskSet<R>::const_iterator state_iter=state_bounding_set.begin();
          state_iter!=state_bounding_set.end(); ++state_iter)
      {
        Point<I> state=static_cast< Rectangle<R> >(*state_iter);

        for(typename GridMaskSet<R>::const_iterator input_iter=input_bounding_set.begin();
            input_iter!=input_bounding_set.end(); ++input_iter)
        {
          std::pair<LatticeCell,LatticeCell> control = std::make_pair(state_iter->lattice_set(),input_iter->lattice_set());
          discretization.insert(std::make_pair(control,LatticeCellListSet(state_space_dimension)));

          Point<I> input = static_cast< Rectangle<R> >(*input_iter);

          for(typename GridMaskSet<R>::const_iterator noise_iter=noise_bounding_set.begin();
              noise_iter!=noise_bounding_set.end(); ++noise_iter)
          {
            Point<I> noise = static_cast< Rectangle<R> >(*noise_iter);
            
            Point<I> image = f.image(state,input,noise);
            
            GridBlock<R> image_set = over_approximation(image,state_grid);
            discretization.find(control)->second.adjoin(image_set.lattice_set());
          }
        }
      }
            
      // Discretize target set
      GridMaskSet<R> target_approximation(state_bounding_set.grid(),state_bounding_set.block());
      target_approximation.adjoin_under_approximation(target_set);
      Combinatoric::LatticeMaskSet target_lattice_set = target_approximation.lattice_set();

      // Discretize initial set
      GridMaskSet<R> initial_approximation(state_bounding_set.grid(),state_bounding_set.block());
      initial_approximation.adjoin_under_approximation(initial_set);
      Combinatoric::LatticeMaskSet initial_lattice_set = initial_approximation.lattice_set();

      Combinatoric::LatticeMaskSet bounding_lattice_set = state_bounding_set.lattice_set();
      Combinatoric::LatticeMaskSet input_lattice_set = input_bounding_set.lattice_set();

      // Solve control problem
      Combinatoric::LatticeMultiMap lattice_control(state_space_dimension,input_space_dimension);
      Combinatoric::LatticeMaskSet controllable_lattice_set(state_bounding_set.block());
      Combinatoric::LatticeMaskSet new_controllable_lattice_set(state_bounding_set.block());;
      
      new_controllable_lattice_set.adjoin(target_lattice_set);
      controllable_lattice_set.adjoin(new_controllable_lattice_set);
      while(!new_controllable_lattice_set.empty() && !subset(initial_lattice_set,controllable_lattice_set)) {
        new_controllable_lattice_set.clear();
        for(LatticeMaskSet::const_iterator state_iter = bounding_lattice_set.begin();
            state_iter!=bounding_lattice_set.end(); ++state_iter)
        {
          if(!subset(*state_iter,controllable_lattice_set)) {
            for(LatticeMaskSet::const_iterator input_iter = input_lattice_set.begin();
                input_iter!=input_lattice_set.end(); ++input_iter)
              
            {
              if(Combinatoric::subset(discretization.find(std::make_pair(*state_iter,*input_iter))->second,controllable_lattice_set)) {
                new_controllable_lattice_set.adjoin(*state_iter);
                lattice_control.adjoin_to_image(*state_iter,*input_iter);
              }
            }
          }
        }
        controllable_lattice_set.adjoin(new_controllable_lattice_set);
      }
      
      System::GridMultiMap<R> result(state_grid,input_grid,lattice_control);
      
      return result;      
    }


  }
}
