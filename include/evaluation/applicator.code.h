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
#include "../geometry/zonotope.h"
#include "../geometry/list_set.h"
#include "../geometry/grid.h"
#include "../geometry/grid_set.h"
#include "../geometry/rectangular_set.h"

#include "../system/grid_multimap.h"


#include "../system/map.h"
#include "../system/discrete_time_system.h"

#include "../output/logging.h"

namespace Ariadne {
  
namespace Evaluation { 
static int& verbosity = applicator_verbosity; 
static const double DEFAULT_MAXIMUM_BASIC_SET_RADIUS=0.25;
static const double DEFAULT_GRID_SIZE=0.125;
}



template<class R>
Evaluation::Applicator<R>::Applicator() 
  : _maximum_basic_set_radius(DEFAULT_MAXIMUM_BASIC_SET_RADIUS),
    _grid_size(DEFAULT_GRID_SIZE)
{
}

template<class R>
Evaluation::Applicator<R>::Applicator(const R& mbsr, const R& gs) 
  : _maximum_basic_set_radius(mbsr),
    _grid_size(gs)
{
}



template<class R>
Evaluation::Applicator<R>::~Applicator() 
{
}


template<class R>
Evaluation::Applicator<R>*
Evaluation::Applicator<R>::clone() const 
{
  return new Applicator<R>(this->_maximum_basic_set_radius,this->_grid_size);
}




template<class R>
R
Evaluation::Applicator<R>::maximum_basic_set_radius() const
{
  return this->_maximum_basic_set_radius;
}


template<class R>
void
Evaluation::Applicator<R>::set_maximum_basic_set_radius(const R& mbsr) 
{
  if(mbsr<=0) {
    throw std::runtime_error("maximum basic set radius must be positive");
  }
  this->_maximum_basic_set_radius=mbsr;
}


template<class R>
R
Evaluation::Applicator<R>::grid_size() const
{
  return this->_grid_size;
}


template<class R>
void
Evaluation::Applicator<R>::set_grid_size(const R& mgs) 
{
  if(mgs<=0) {
    throw std::runtime_error("maximum basic set radius must be positive");
  }
  this->_grid_size=mgs;
}


template<class R>
R
Evaluation::Applicator<R>::default_bound() const
{
  return this->_default_bound;
}


template<class R>
void
Evaluation::Applicator<R>::set_default_bound(const R& db) 
{
  if(db<=0) {
    throw std::runtime_error("default bound must be positive");
  }
  this->_default_bound=db;
}



template<class R>
Geometry::ListSet< Geometry::Rectangle<R> >
Evaluation::Applicator<R>::subdivide(const Geometry::Rectangle<R>& s) const
{
  return s.subdivide();
}


template<class R>
Geometry::ListSet< Geometry::Zonotope<R,R> >
Evaluation::Applicator<R>::subdivide(const Geometry::Zonotope<R,R>& s) const
{
  return s.subdivide();
}


template<class R>
Geometry::ListSet< Geometry::Zonotope<Numeric::Interval<R>,R> >
Evaluation::Applicator<R>::subdivide(const Geometry::Zonotope<I,R>& s) const
{
  return s.subdivide();
}


template<class R>
Geometry::ListSet< Geometry::Zonotope< Numeric::Interval<R>, Numeric::Interval<R> > >
Evaluation::Applicator<R>::subdivide(const Geometry::Zonotope<I,I>& s) const
{
  return s.subdivide();
}



template<class R>
Geometry::Rectangle<R> 
Evaluation::Applicator<R>::evaluate(const System::MapInterface<R>& f, const Geometry::Rectangle<R>& r) const
{
  ARIADNE_LOG(6,"GridMaskSetEvaluation Applicator::evaluate(MapInterface f, Rectangle<Float> r)\n");
  ARIADNE_LOG(7,"  r="<<r<<"\n");
  ARIADNE_LOG(8,"  f(r)="<<Geometry::Rectangle<R>(f.image(Geometry::Point< Numeric::Interval<R> >(r)))<<"\n");
  return Geometry::Rectangle<R>(f.image(Geometry::Point< Numeric::Interval<R> >(r)));
}





template<class R>
Geometry::Zonotope<R> 
Evaluation::Applicator<R>::evaluate(const System::MapInterface<R>& f, const Geometry::Zonotope<R>& z) const 
{
  ARIADNE_LOG(6,"GridMaskSetEvaluation Applicator::evaluate(MapInterface f, Zonotope<Float,Float> z)\n");
  ARIADNE_LOG(7,"  z="<<z<<"\n");
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
  Geometry::Point<R> nc=Geometry::midpoint(img_centre);
  LinearAlgebra::Matrix<R> ng=midpoint(img_generators);
  
  Geometry::Zonotope<R> result(nc,ng);
  ARIADNE_LOG(8,"  f(z)="<<result<<"\n");
  return result;
}





template<class R>
Geometry::Zonotope<Numeric::Interval<R>,R> 
Evaluation::Applicator<R>::evaluate(const System::MapInterface<R>& f, const Geometry::Zonotope<Numeric::Interval<R>,R>& z) const 
{
  ARIADNE_LOG(6,"GridMaskSetEvaluation Applicator::evaluate(MapInterface f, Zonotope<Interval,Float> z)\n");
  ARIADNE_LOG(7,"  z="<<z<<"\n");
  typedef Numeric::Interval<R> I;
  
  Geometry::Point<I> img_centre=f(z.centre());
  LinearAlgebra::Matrix<I> df_on_set = f.jacobian(over_approximation(z.bounding_box()));
  LinearAlgebra::Matrix<I> img_generators = df_on_set*z.generators();
  Geometry::Zonotope<I,I> interval_zonotope(img_centre,img_generators);
  Geometry::Zonotope<I,R> result(over_approximation(interval_zonotope));
  ARIADNE_LOG(8,"  f(z)="<<result<<"\n");
  return result;
}



template<class R>
Geometry::Zonotope< Numeric::Interval<R> > 
Evaluation::Applicator<R>::evaluate(const System::MapInterface<R>& f, const Geometry::Zonotope< Numeric::Interval<R> >& z) const 
{
  ARIADNE_LOG(6,"GridMaskSet Applicator::evaluate(MapInterface f, Zonotope<Interval,Interval> z)\n");
  ARIADNE_LOG(7,"  z="<<z<<"\n");
  typedef Numeric::Interval<R> I;
  
  Geometry::Point<I> img_centre=f(z.centre());
  LinearAlgebra::Matrix<I> df_on_set = f.jacobian(over_approximation(z.bounding_box()));
  LinearAlgebra::Matrix<I> img_generators = df_on_set*z.generators();
  
  Geometry::Zonotope<I> result(img_centre,img_generators);
  ARIADNE_LOG(8,"  f(z)="<<result<<"\n");
  return result;
}






template<class R>
template<class BS>
Geometry::ListSet<BS> 
Evaluation::Applicator<R>::image_list_set(const System::MapInterface<R>& f, const Geometry::ListSet<BS>& ds) const 
{
  Geometry::ListSet<BS> result(f.result_dimension());
  for(typename Geometry::ListSet<BS>::const_iterator iter=ds.begin(); iter!=ds.end(); ++iter) {
    result.push_back(this->evaluate(f,*iter));
  }
  return result;
}



template<class R>
Geometry::ListSet< Geometry::Rectangle<R> > 
Evaluation::Applicator<R>::image(const System::MapInterface<R>& f, const Geometry::ListSet< Geometry::Rectangle<R> >& ds) const 
{
  ARIADNE_LOG(2,"GridMaskSet Applicator::image(MapInterface map, ListSet< Rectangle<Float> > initial_set)\n");
  ARIADNE_LOG(3,"initial_set="<<ds<<"\n");
  return this->image_list_set(f,ds);
}



template<class R>
Geometry::ListSet< Geometry::Zonotope<R> > 
Evaluation::Applicator<R>::image(const System::MapInterface<R>& f, const Geometry::ListSet< Geometry::Zonotope<R> >& ds) const 
{
  ARIADNE_LOG(2,"GridMaskSet Applicator::image(MapInterface map, ListSet< Zonotope<Float,Float> > initial_set)\n")
    ARIADNE_LOG(3,"initial_set="<<ds<<"\n");
  return this->image_list_set(f,ds);
}



template<class R>
Geometry::ListSet< Geometry::Zonotope<Numeric::Interval<R>,R> > 
Evaluation::Applicator<R>::image(const System::MapInterface<R>& f, const Geometry::ListSet< Geometry::Zonotope<Numeric::Interval<R>,R> >& ds) const 
{
  ARIADNE_LOG(2,"GridMaskSet Applicator::image(MapInterface map, ListSet< Zonotope<Interval,Float> > initial_set)\n");
  ARIADNE_LOG(3,"initial_set="<<ds<<"\n");
  return this->image_list_set(f,ds);
}


template<class R>
Geometry::ListSet< Geometry::Zonotope< Numeric::Interval<R> > > 
Evaluation::Applicator<R>::image(const System::MapInterface<R>& f, const Geometry::ListSet< Geometry::Zonotope< Numeric::Interval<R> > >& ds) const 
{
  ARIADNE_LOG(2,"GridMaskSet Applicator::image(MapInterface map, ListSet< Zonotope<Interval,Interval> > initial_set)\n");
  ARIADNE_LOG(3,"initial_set="<<ds<<"\n");
  return this->image_list_set(f,ds);
}


template<class R>
Geometry::GridCellListSet<R> 
Evaluation::Applicator<R>::image(const System::MapInterface<R>& f, 
                                 const Geometry::GridCellListSet<R>& initial_set, 
                                 const Geometry::Grid<R>& image_grid) const 
{
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  typedef typename GridCellListSet<R>::const_iterator gcls_const_iterator;
  ARIADNE_LOG(2,"GridMaskSet Applicator::image(MapInterface map, GridCellListSet initial_set, Grid image_grid)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\nimage_grid="<<image_grid<<"\n");
  
  GridCellListSet<R> image(image_grid);
  
  Rectangle<R> r(f.argument_dimension());
  Rectangle<R> fr(f.result_dimension());
  Zonotope<I,R> z(f.argument_dimension());
  Zonotope<I,R> fz(f.result_dimension());
  
  for(gcls_const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
    Rectangle<R> r=*iter;
    if(f.smoothness()>=1) {
      z=r;
      fz=this->evaluate(f,z);
      image.adjoin_outer_approximation(fz);
    } else {
      fr=this->evaluate(f,r);
      image.adjoin_outer_approximation(fr);
    }
  }
  return image;
}




template<class R>
Geometry::GridMaskSet<R> 
Evaluation::Applicator<R>::image(const System::MapInterface<R>& f, 
                                 const Geometry::GridMaskSet<R>& initial_set, 
                                 const Geometry::GridMaskSet<R>& bounding_set) const 
{
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  typedef typename GridMaskSet<R>::const_iterator gms_const_iterator;
  ARIADNE_LOG(2,"GridMaskSet Applicator::image(MapInterface f, GridMaskSet initial_set, GridMaskSet bounding_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\nbounding_set="<<bounding_set);
  ARIADNE_CHECK_BOUNDED(initial_set,"GridMaskSet Applicator<R>::image(MapInterface,GridMaskSet,GridMaskSet)");
  ARIADNE_CHECK_BOUNDED(bounding_set,"Applicator<R>::image(MapInterface,GridMaskSet,GridMaskSet)");
  
  const Grid<R>& g=initial_set.grid();
  Combinatoric::LatticeBlock bd=bounding_set.block();
  GridMaskSet<R> image(g,bd);
  Rectangle<R> bb=initial_set.bounding_box();
  
  Rectangle<R> r(g.dimension());
  Rectangle<R> fr(g.dimension());
  Zonotope<I,R> z(g.dimension());
  Zonotope<I,R> fz(g.dimension());
  
  for(gms_const_iterator iter=initial_set.begin(); iter!=initial_set.end(); ++iter) {
    Rectangle<R> r=*iter;
    if(f.smoothness()>=1) {
      z=r;
      fz=this->evaluate(f,z);
      ARIADNE_LOG(7,"z="<<z<<", fz="<<fz<<"\n");
      GridCellListSet<R> gz=outer_approximation(fz,g);
      ARIADNE_LOG(7,"gz="<<gz);
      image.adjoin(gz);
      ARIADNE_LOG(9,"image.size()="<<image.size()<<"\n");
      assert((bool)(subset(gz,image)));
    } else {
      fr=this->evaluate(f,r);
      image.adjoin_outer_approximation(fr);
    }
  }
  return regular_intersection(image,bounding_set);
}





template<class R>
Geometry::GridMaskSet<R> 
Evaluation::Applicator<R>::preimage(const System::MapInterface<R>& f, 
                                    const Geometry::GridMaskSet<R>& set, 
                                    const Geometry::GridMaskSet<R>& bounding_set) const 
{
  ARIADNE_LOG(2,"GridMaskSet Applicator::preimage(MapInterface,GridMaskSet,GridMaskSet)\n");
  ARIADNE_LOG(3,"set="<<set<<"\nbounding_set="<<bounding_set);
  using namespace Numeric;
  using namespace Geometry;
  typedef typename GridMaskSet<R>::const_iterator basic_set_iterator;
  GridMaskSet<R> result(bounding_set);
  result.clear();
  Rectangle<R> r;
  Zonotope<Interval<R>,R> z;
  Zonotope<Interval<R>,R> fz;
  GridCellListSet<R> fgcls(set.grid());
  ARIADNE_LOG(7,"Preimage testing "<<bounding_set.size()<<" cells\n");
  size_type tested=0;
  for(typename GridMaskSet<R>::const_iterator bnd_iter=bounding_set.begin(); 
      bnd_iter!=bounding_set.end(); ++bnd_iter)
    {
      if(tested%256==0 && tested!=0) {
        ARIADNE_LOG(7,"Preimage tested "<<tested<<" cells; found "<<result.size()<<" cells in preimage\n");
      }
      ++tested;
      fgcls.clear();
      r=*bnd_iter;
      z=r;
      fz=this->evaluate(f,z);
      fgcls.adjoin_outer_approximation(fz);
      if(subset(fgcls,set)) {
        result.adjoin(*bnd_iter);
      }
    }
  return result;
}



template<class R>
Geometry::ListSet< Geometry::Zonotope<Numeric::Interval<R>,R> > 
Evaluation::Applicator<R>::reach(const System::MapInterface<R>& f, 
                                 const Geometry::ListSet< Geometry::Zonotope<Numeric::Interval<R>,R> >& initial_set) const 
{
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  ARIADNE_LOG(2,"GridMaskSet Applicator::reach(MapInterface,ListSet< Zonotope<I,F> >\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  ListSet< Zonotope<I,R> > result; 
  Zonotope<I,R> bs;
  for(typename ListSet< Zonotope<I,R> >::const_iterator bs_iter=initial_set.begin();
      bs_iter!=initial_set.end(); ++bs_iter)
    {
      ARIADNE_LOG(6,"  computing reach for bs="<<*bs_iter);
      size_type steps=0;
      bs=*bs_iter;
      while(bs.radius()<this->maximum_basic_set_radius()) {
        result.adjoin(bs);
        bs=this->evaluate(f,bs);
        ++steps;
      }
      ARIADNE_LOG(6,"  iterated "<<steps<<" time steps");
    }
  return result;
}


template<class R>
Geometry::GridMaskSet<R> 
Evaluation::Applicator<R>::reach(const System::MapInterface<R>& f, 
                                 const Geometry::GridMaskSet<R>& initial_set) const 
{
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  ARIADNE_LOG(2,"GridMaskSet Applicator::reach(MapInterface,GridMaskSet)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  GridMaskSet<R> result(initial_set.grid(),initial_set.block());
  typedef typename GridMaskSet<R>::const_iterator basic_set_iterator;
  Rectangle<R> rectangle(initial_set.dimension());
  Zonotope< Numeric::Interval<R> > zonotope(result.dimension());
  for(basic_set_iterator bs_iter=initial_set.begin(); 
      bs_iter!=initial_set.end(); ++bs_iter)
    {
      ARIADNE_LOG(6,"  computing reach for cell="<<*bs_iter);
      size_type steps=0;
      rectangle=*bs_iter;
      zonotope=rectangle;
      while(zonotope.radius() < this->maximum_basic_set_radius()) {
        result.adjoin_outer_approximation(zonotope);
        zonotope=this->evaluate(f,zonotope);
        ++steps;
      }
      ARIADNE_LOG(6,"  iterated "<<steps<<" time steps");
    }
  return result;
}


template<class R>
Geometry::GridMaskSet<R> 
Evaluation::Applicator<R>::chainreach(const System::MapInterface<R>& f, 
                                      const Geometry::GridMaskSet<R>& initial_set, 
                                      const Geometry::GridMaskSet<R>& bounding_set) const
{
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
  ARIADNE_LOG(2,"GridMaskSet Applicator::chainreach(MapInterface map, GridMaskSet initial_set, GridMaskSet bounding_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\nbounding_set="<<bounding_set);
  ARIADNE_CHECK_BOUNDED(initial_set,"GridMaskSet Applicator<R>::chainreach(MapInterface,GridMaskSet,GridMaskSet)");
  ARIADNE_CHECK_BOUNDED(bounding_set,"GridMaskSet Applicator<R>::chainreach(MapInterface,GridMaskSet,GridMaskSet)");
  
  const Grid<R>& g=bounding_set.grid();
  Combinatoric::LatticeBlock bd=bounding_set.block();
  GridBlock<R> bb(g,bd);
  GridMaskSet<R> result(g,bd);
  GridCellListSet<R> found(g);
  GridCellListSet<R> image(g);
  
  Rectangle<R> r(g.dimension());
  Rectangle<R> fr(g.dimension());
  Zonotope<I,R> z(g.dimension());
  Zonotope<I,R> fz(g.dimension());
  
  uint step=0;
  found=initial_set;
  while(!subset(found,result)) {
    ARIADNE_LOG(3,"Chainreach step "<<step<<": found "<<found.size()<<" cells, ");
    found=difference(found,result);
    ARIADNE_LOG(3,found.size()<<" of which are new.\n");
    ARIADNE_LOG(3,"reached "<<result.size()<<" cells in total.\n");
    result.adjoin(found);
    image.clear(); 
    uint size=0;
    for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
      ++size;
      Rectangle<R> r=*iter;
      if(f.smoothness()>0) {
        z=r;
        fz=this->evaluate(f,z);
        if(!disjoint(fz.bounding_box(),bounding_set)) {
          image.adjoin(outer_approximation(fz,g));
        }
      } else {
        fr=this->evaluate(f,r);
        if(!disjoint(fr.bounding_box(),bounding_set)) {
          image.adjoin(outer_approximation(fr,g));
        }
      }
    }
    image.unique_sort();
    found=regular_intersection(image,bounding_set);
    ++step;
  }
  return result;
}



template<class R>
Geometry::GridMaskSet<R>
Evaluation::Applicator<R>::viable(const System::MapInterface<R>& f, 
                                  const Geometry::GridMaskSet<R>& bounding_set) const
{
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  typedef typename Geometry::GridMaskSet<R>::const_iterator gms_const_iterator;
  ARIADNE_LOG(2,"GridMaskSet Applicator::viable(MapInterface map, GridMaskSet bounding_set)\n");
  ARIADNE_LOG(3,"bounding_set="<<bounding_set<<"\n");
  ARIADNE_CHECK_BOUNDED(bounding_set,"Applicator<R>::viable(MapInterface,GridMaskSet)");
  
  const Grid<R>& g=bounding_set.grid();
  Combinatoric::LatticeBlock bd=bounding_set.block();
  GridBlock<R> bb(g,bd);
  GridMaskSet<R> result(g,bd);
  GridCellListSet<R> unsafe(g);
  
  Rectangle<R> r(g.dimension());
  Rectangle<R> fr(g.dimension());
  GridBlock<R> fgb(g);
  Zonotope<R> z(g.dimension());
  Zonotope<R> fz(g.dimension());
  GridCellListSet<R> fgcls(g);
  
  result=bounding_set;
  size_type step=0;

  ARIADNE_LOG(3,"Computing discretization...\n");
  System::GridMultiMap<R> discretization=this->discretize(f,bounding_set,bounding_set.grid());
  ARIADNE_LOG(3,"   Done computing discretization.\n");
  do {
    ARIADNE_LOG(3,"Viability step "<<step+1<< ": testing "<<result.size()<<" cells.\n");
    ++step;
    unsafe.clear();
    for(gms_const_iterator iter=result.begin(); iter!=result.end(); ++iter) {
      fgcls=discretization.image(*iter);
      ARIADNE_LOG(7,"cell="<<*iter<<", image.size()="<<fgcls.size()<<"\n");
      if(!overlap(result,fgcls)) {
        unsafe.adjoin(*iter);
      }
    }
    result.remove(unsafe);
  } while(!unsafe.empty());
  return result;
}


template<class R>
tribool
Evaluation::Applicator<R>::verify(const System::MapInterface<R>& f, 
                                  const Geometry::GridMaskSet<R>& initial_set, 
                                  const Geometry::GridMaskSet<R>& safe_set) const
{
  ARIADNE_LOG(2,"triboolEvaluation::Applicator::verify(MapInterface map, GridMaskSet initial_set,GridMaskSet  safe_set)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n"<<"safe_set="<<safe_set);
  typedef Numeric::Interval<R> I;
  using namespace Geometry;
  typedef typename Geometry::GridCellListSet<R>::const_iterator gcls_const_iterator;
  ARIADNE_CHECK_BOUNDED(initial_set,"Applicator<R>::verify(...)");
  ARIADNE_CHECK_BOUNDED(safe_set,"Applicator<R>::verify(...)");
  
  const Grid<R>& g=initial_set.grid();
  Combinatoric::LatticeBlock bd=safe_set.block();
  Rectangle<R> bb=safe_set.bounding_box();
  GridMaskSet<R> reach(g,bd);
  GridCellListSet<R> found(g);
  GridCellListSet<R> cell_image(g);
  GridCellListSet<R> image(g);
  
  Rectangle<R> r(g.dimension());
  Rectangle<R> fr(g.dimension());
  Zonotope<I,R> z(g.dimension());
  Zonotope<I,R> fz(g.dimension());
  
  found=initial_set;
  while(!subset(found,reach)) {
    found=difference(found,reach);
    reach.adjoin(found);
    image.clear();
    uint size=0;
    for(gcls_const_iterator iter=found.begin(); iter!=found.end(); ++iter) {
      ++size;
      Rectangle<R> r=*iter;
      if(f.smoothness()>=1) {
        z=r;
        fz=this->evaluate(f,z);
        cell_image.adjoin(outer_approximation(fz,g));
        if(!subset(cell_image,safe_set)) {
          return false;
        }
      } else {
        fr=this->evaluate(f,r);
        cell_image.adjoin(outer_approximation(fr,g));
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
Evaluation::Applicator<R>::image(const System::MapInterface<R>& f, 
                                 const Geometry::SetInterface<R>& set) const
{
  // FIXME: Only computes an over-approximation to the image.
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  ARIADNE_LOG(2,"SetInterface* Applicator::image(MapInterface,SetInterface)\n");
  ARIADNE_LOG(3,"set="<<set<<"\n");
  ListSet< Zonotope<I,R> >* result=new ListSet< Zonotope<I,R> >;
  Rectangle<R> bb=set.bounding_box();
  Grid<R> arg_grid(f.argument_dimension(),this->grid_size());
  GridMaskSet<R> grid_arg_set(arg_grid,bb);
  grid_arg_set.adjoin_outer_approximation(set);
  Rectangle<R> r(f.argument_dimension());
  Zonotope<I,R> z(f.argument_dimension());
  for(typename GridMaskSet<R>::const_iterator cell_iter=grid_arg_set.begin();
      cell_iter!=grid_arg_set.end(); ++cell_iter)
    {
      z=(r=*cell_iter);
      z=this->evaluate(f,z);
      result->adjoin(z);
    }
  return result;
}



template<class R>
Geometry::SetInterface<R>*
Evaluation::Applicator<R>::preimage(const System::MapInterface<R>& map, 
                                    const Geometry::SetInterface<R>& set,
                                    const Geometry::SetInterface<R>& bound) const
{
  typedef Numeric::Interval<R> I;
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* Applicator::preimage(MapInterface map, SetInterface set)\n");
  ARIADNE_LOG(3,"set="<<set<<"\n");
  Rectangle<R> preimage_bounding_box=bound.bounding_box();
  Grid<R> preimage_grid(map.argument_dimension(),this->grid_size());
  Grid<R> image_grid(map.result_dimension(),this->grid_size());
  GridMaskSet<R> grid_image_set(image_grid,set.bounding_box());
  grid_image_set.adjoin_inner_approximation(set);
  GridMaskSet<R> grid_preimage_bounding_set(preimage_grid,preimage_bounding_box);
  grid_preimage_bounding_set.adjoin_outer_approximation(bound);
  return new GridMaskSet<R>(this->preimage(map,grid_image_set,grid_preimage_bounding_set));
}



template<class R>
Geometry::SetInterface<R>*
Evaluation::Applicator<R>::reach(const System::MapInterface<R>& f, 
                                 const Geometry::SetInterface<R>& initial_set) const
{
  typedef Numeric::Interval<R> I;
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* Applicator::reach(MapInterface,SetInterface)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n");
  const ListSet< Zonotope<I,R> >* list_initial_set_ptr=dynamic_cast<const ListSet< Zonotope<I,R> >*>(&initial_set);
  ARIADNE_LOG(4,"list_initial_set_ptr="<<list_initial_set_ptr<<"\n");
  if(list_initial_set_ptr) {
    return new ListSet< Zonotope<I,R> >(this->reach(f,*list_initial_set_ptr));
  }
  
  ARIADNE_CHECK_BOUNDED(initial_set,"SetInterface* Applicator::reach(MapInterface,SetInterface)");
  
  ListSet< Zonotope<I,R> > list_initial_set;
  if(dynamic_cast<const RectangularSet<R>*>(&initial_set)) {
    ARIADNE_LOG(4,"Cast to RectangularSet<R>\n");
  }
  
  Grid<R> grid(initial_set.dimension(),this->grid_size());
  Rectangle<R> bounding_box=initial_set.bounding_box();
  GridMaskSet<R> gms(grid,bounding_box);
  gms.adjoin_inner_approximation(initial_set);
  ARIADNE_LOG(4,"gms="<<gms);
  list_initial_set=ListSet< Zonotope<I,R> >(gms);
  
  return new ListSet< Zonotope<I,R> >(this->reach(f,list_initial_set));      
}



template<class R>
Geometry::SetInterface<R>*
Evaluation::Applicator<R>::chainreach(const System::MapInterface<R>& map, 
                                      const Geometry::SetInterface<R>& initial_set, 
                                      const Geometry::SetInterface<R>& bounding_set) const
{
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* Applicator::chainreach(MapInterface,SetInterface,SetInterface)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n"<<"bounding_set="<<bounding_set<<"\n");
  
  ARIADNE_CHECK_BOUNDED(bounding_set,"SetInterface* Applicator::chainreach(MapInterface map, SetInterface initial_set, SetInterface bounding_set)");
  Rectangle<R> bounding_box=bounding_set.bounding_box();
  Grid<R> grid(bounding_set.dimension(),this->grid_size());
  FiniteGrid<R> finite_grid(grid,bounding_box);
  GridMaskSet<R> grid_bounding_set(finite_grid);
  grid_bounding_set.adjoin_outer_approximation(bounding_set);
  GridMaskSet<R> grid_initial_set(finite_grid);
  grid_initial_set.adjoin_outer_approximation(initial_set);
  
  return new GridMaskSet<R>(this->chainreach(map,grid_initial_set,grid_bounding_set));
}




template<class R>
Geometry::SetInterface<R>*
Evaluation::Applicator<R>::viable(const System::MapInterface<R>& map, 
                                  const Geometry::SetInterface<R>& bounding_set) const
{
  using namespace Geometry;
  ARIADNE_LOG(2,"SetInterface* Applicator::viable(MapInterface,SetInterface)\n");
  ARIADNE_LOG(3,"bounding_set="<<bounding_set<<"\n");
  ARIADNE_CHECK_BOUNDED(bounding_set,"SetInterface* Applicator::viable(MapInterface map, SetInterface bounding_set)");
  Rectangle<R> bounding_box=bounding_set.bounding_box();
  Grid<R> grid(bounding_set.dimension(),this->grid_size());
  GridMaskSet<R> grid_bounding_set(grid,bounding_box);
  grid_bounding_set.adjoin_outer_approximation(bounding_set);
  return new GridMaskSet<R>(this->viable(map,grid_bounding_set));
}



template<class R>
tribool
Evaluation::Applicator<R>::verify(const System::MapInterface<R>& f, 
                                  const Geometry::SetInterface<R>& initial_set, 
                                  const Geometry::SetInterface<R>& safe_set) const
{
  using namespace Geometry;
  ARIADNE_LOG(2,"triboolEvaluation::Applicator::verify(MapInterface,SetInterface,SetInterface)\n");
  ARIADNE_LOG(3,"initial_set="<<initial_set<<"\n"<<"safe_set="<<safe_set<<"\n");
  ARIADNE_CHECK_BOUNDED(safe_set,"SetInterface* Applicator::verify(MapInterface map, SetInterface initial_set, SetInterface safe_set)");
  Rectangle<R> bounding_box=safe_set.bounding_box();
  Grid<R> grid(safe_set.dimension(),this->grid_size());
  FiniteGrid<R> finite_grid(grid,bounding_box);
  GridMaskSet<R> grid_inner_safe_set(finite_grid);
  grid_inner_safe_set.adjoin_inner_approximation(safe_set);
  GridMaskSet<R> grid_initial_set(finite_grid);
  grid_initial_set.adjoin_outer_approximation(initial_set);
  
  return this->verify(f,grid_initial_set,grid_inner_safe_set);
}




template<class R>
System::GridMultiMap<R> 
Evaluation::Applicator<R>::discretize(const System::MapInterface<R>& f, 
                                      const Geometry::GridMaskSet<R>& domain,
                                      const Geometry::Grid<R>& range_grid) const
{
  ARIADNE_LOG(2,"GridMultiMap*Evaluation::Applicator::discretize(MapInterface map, GridMaskSet domain, Grid range_grid)\n");
  ARIADNE_LOG(3,"domain="<<domain<<"\n"<<"range_grid="<<range_grid);
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  System::GridMultiMap<R> result(domain.grid(),range_grid);
  Zonotope<R> basic_set;
  Zonotope<R> image_set;
  for(typename GridMaskSet<R>::const_iterator dom_iter=domain.begin();
      dom_iter!=domain.end(); ++dom_iter)
    {
      const GridCell<R>& cell=*dom_iter;
      basic_set=cell;
      image_set=this->evaluate(f,basic_set);
      result.adjoin_to_image(cell,outer_approximation(image_set,range_grid));
    }
  return result;
}



template<class R>
System::GridMultiMap<R> 
Evaluation::Applicator<R>::control_synthesis(const System::DiscreteTimeSystem<R>& f, 
                                             const Geometry::SetInterface<R>& initial_set,
                                             const Geometry::SetInterface<R>& target_set,
                                             const Geometry::GridMaskSet<R>& state_bounding_set,
                                             const Geometry::GridMaskSet<R>& input_bounding_set,
                                             const Geometry::GridMaskSet<R>& noise_bounding_set) const
{
  // TODO: Use on-the-fly discretization
  // TODO: Use adaptive grid refinement
  
  ARIADNE_LOG(2,"GridMultiMap* Applicator::control_synthesis(...)\n");
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
              
              GridBlock<R> image_set = outer_approximation(image,state_grid);
              discretization.find(control)->second.adjoin(image_set.lattice_set());
            }
        }
    }
  
  // Discretize target set
  GridMaskSet<R> target_approximation(state_bounding_set.grid(),state_bounding_set.block());
  target_approximation.adjoin_inner_approximation(target_set);
  Combinatoric::LatticeMaskSet target_lattice_set = target_approximation.lattice_set();
  
  // Discretize initial set
  GridMaskSet<R> initial_approximation(state_bounding_set.grid(),state_bounding_set.block());
  initial_approximation.adjoin_inner_approximation(initial_set);
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
