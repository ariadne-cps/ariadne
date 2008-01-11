/***************************************************************************
 *            grid_approximation.code.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "grid.h"
#include "grid_cell.h"
#include "grid_block.h"
#include "grid_cell_list_set.h"
#include "grid_mask_set.h"

#include "point.h"
#include "box.h"
#include "rectangle.h"
#include "zonotope.h"
#include "polytope.h"
#include "polyhedron.h"
#include "box_list_set.h"
#include "list_set.h"

#include "set_interface.h"

#include "grid_approximation.h"

#include "output/logging.h"



namespace {

using namespace Ariadne;


template<class R, class BS>
inline
Geometry::GridCellListSet<R>
outer_approximation_of_basic_set(const BS& bs, const Geometry::Grid<R>& g) 
{
  //std::cout << __PRETTY_FUNCTION__ << std::endl;
  Geometry::GridCellListSet<R> gcls(g);
  Geometry::GridBlock<R> gbb=outer_approximation(Geometry::bounding_box(bs),g);
  Geometry::Box<R> r(bs.dimension());
  //std::cerr << "bs="<<bs<<" bb="<<bs.bounding_box()<<" gbb="<<gbb<<std::endl;
  for(typename Geometry::GridBlock<R>::const_iterator iter=gbb.begin(); iter!=gbb.end(); ++iter) {
    r=*iter;
    if(disjoint(bs,r)) {
      //std::cerr << "disjoint("<<r<<","<<bs<<")\n";
    } else {
      //std::cerr << "adjoining "<<*iter<<"\n";
      gcls.adjoin(*iter);
    }
  }
  return gcls;
}


template<class R, class BS>
inline
Geometry::GridCellListSet<R>
inner_approximation_of_basic_set(const BS& bs, const Geometry::Grid<R>& g) 
{
  Geometry::GridCellListSet<R> gcls(g);
  Geometry::GridBlock<R> gbb=outer_approximation(Geometry::bounding_box(bs),g);
  Geometry::Box<R> r(bs.dimension());
  
  for(typename Geometry::GridBlock<R>::const_iterator iter=gbb.begin(); iter!=gbb.end(); ++iter) {
    r=*iter;
    if(subset(r,bs)) {
      gcls.adjoin(*iter);
    }
  }
  return gcls;
}


template<class R, class DS>
inline
Geometry::GridCellListSet<R>
outer_approximation_of_denotable_set(const DS& ds, const Geometry::Grid<R>& g) 
{
  //std::cout << __PRETTY_FUNCTION__ << std::endl;
  Geometry::GridCellListSet<R> gcls(g);
  for(typename DS::const_iterator iter=ds.begin(); iter!=ds.end(); ++iter) {
    gcls.adjoin(outer_approximation(*iter,g));
  }
  gcls.unique_sort();
  return gcls;
}

template<class R, class DS>
inline
Geometry::GridMaskSet<R>
outer_approximation_of_denotable_set(const DS& ds, const Geometry::FiniteGrid<R>& fg) 
{
  //std::cout << __PRETTY_FUNCTION__ << std::endl;
  Geometry::GridMaskSet<R> gms(fg);
  Geometry::Box<R> r(ds.dimension());
  for(typename DS::const_iterator bs_iter=ds.begin(); bs_iter!=ds.end(); ++bs_iter) {
    Geometry::GridBlock<R> gbb=outer_approximation(bs_iter->bounding_box(),fg.grid());
    for(typename Geometry::GridBlock<R>::const_iterator gc_iter=gbb.begin(); gc_iter!=gbb.end(); ++gc_iter) {
      r=*gc_iter;
      if(not subset(*gc_iter,gms) && possibly(not disjoint(*bs_iter,r))) {
        gms.adjoin(*gc_iter);
      }
    }
  }
  return gms;
}

template<class R>
inline
Geometry::GridCellListSet<R>
fuzzy_outer_approximation_of_zonotope(const Geometry::Zonotope<R>& z, const Geometry::Grid<R>& g) 
{
  typedef Numeric::Interval<R> I;
  Geometry::GridCellListSet<R> gcls(g);

  LinearAlgebra::Vector<I> ze;
  Geometry::Point<I> pt;
  Geometry::Box<R> r;

  dimension_type d=z.dimension();

  if(z.dimension()==z.number_of_generators()) {
    LinearAlgebra::Matrix<I> Ginv=LinearAlgebra::inverse(z.generators());
    const Geometry::Point<I> c=z.centre()+I(-1,1)*z.error();
    LinearAlgebra::Vector<I> e(z.dimension(),I(-1,1));
    Geometry::GridBlock<R> gbb=outer_approximation(Geometry::bounding_box(z),g);
    for(typename Geometry::GridBlock<R>::const_iterator iter=gbb.begin(); iter!=gbb.end(); ++iter) {
      r=*iter;
      pt=r;
      ze=Ginv*(pt-c);
      if(possibly(ze==e)) {
        gcls.adjoin(*iter);
      }
    }
    return gcls;
  } else if(z.dimension()+1u==z.number_of_generators()) {
    LinearAlgebra::Vector<I> e(d+1,I(-1,1));
    LinearAlgebra::Vector<I> x(d);
    const Geometry::Point<I> c=z.centre();
    const LinearAlgebra::Matrix<I> G=z.generators();
    LinearAlgebra::Matrix<I> A=LinearAlgebra::MatrixSlice<const I>(d,d,G.begin(),G.row_increment(),G.column_increment());
    LinearAlgebra::Vector<I> b=G.column(d);
    LinearAlgebra::Matrix<I> Ainv=LinearAlgebra::inverse(A);
    LinearAlgebra::Vector<I> y=Ainv*b;
    
    Geometry::GridBlock<R> gbb=outer_approximation(Geometry::bounding_box(z),g);
    for(typename Geometry::GridBlock<R>::const_iterator iter=gbb.begin(); iter!=gbb.end(); ++iter) {
      r=*iter;
      pt=r;
      x=Ainv*(pt-c);
      I err(-1,1);
      for(size_type i=0; i!=d; ++i) {
        I diff=x[i]-e[i];
        if(Numeric::encloses(y[i],static_cast<R>(0))) {
          if(!Numeric::encloses(I(x[i]-y[i]*err),0)) {
            err=I(1,0);
          }
        } else {
          err=intersection(err,I(diff/y[i]));
        }
        if(Numeric::empty(err)) {
          break;
        }
      }
      if(!Numeric::empty(err)) {
        gcls.adjoin(*iter);
      }
    }
    return gcls;
  } else {
    return outer_approximation_of_basic_set(z,g);
  }
}




}





namespace Ariadne {

template<class R>
void
Geometry::instantiate_grid_approximation()
{
  typedef Numeric::Interval<R> I;
  tribool tb;
  Point<I>* ipt=0;
  Box<R>* bx=0;
  Rectangle<R>* r=0;
  Zonotope<R>* z=0;
  Polytope<R>* pltp=0;
  Polyhedron<R>* plhd=0;
  SetInterface<R>* set=0;
  
  Grid<R>* g=0;
  FiniteGrid<R>* fg=0;
  GridBlock<R>* gb=0;
  GridCellListSet<R>* gcls=0;
  GridMaskSet<R>* gms=0;
  BoxListSet<R>* bxls=0;

  ListSet< Rectangle<R> >* rls=0;
  ListSet< Zonotope<R> >* zls=0;
  
  *gb=outer_approximation(*ipt,*g);
  
  *gb=over_approximation(*bx,*g);
  *gb=under_approximation(*bx,*g);
  *gb=outer_approximation(*bx,*g);
  *gb=inner_approximation(*bx,*g);
  
  *gcls=outer_approximation(*r,*g);
  *gcls=inner_approximation(*r,*g);

  *gcls=outer_approximation(*pltp,*g);
  *gcls=outer_approximation(*plhd,*g);
  *gcls=outer_approximation(*z,*g);
  *gcls=outer_approximation(*set,*g);
  
  //  *gcls=fuzzy_outer_approximation(*pltp,*g);
  //  *gcls=fuzzy_outer_approximation(*plhd,*g);
  *gcls=fuzzy_outer_approximation(*z,*g);
  
  *gcls=inner_approximation(*pltp,*g);
  *gcls=inner_approximation(*plhd,*g);
  *gcls=inner_approximation(*z,*g);
  *gcls=inner_approximation(*set,*g);

  *gcls=outer_approximation(*bxls,*g);
  *gms=outer_approximation(*bxls,*fg);
  
  *gms=outer_approximation(*rls,*fg);
  *gms=outer_approximation(*zls,*fg);
  *gcls=outer_approximation(*zls,*g);
  
  *gms=outer_approximation(*set,*fg);
  *gms=inner_approximation(*set,*fg);

  *bxls=lower_approximation(*set,*g);
  *bxls=lower_approximation(*set,*fg);
  
  *bxls=point_approximation(*set,*g);
  *bxls=point_approximation(*set,*fg);
  
}


// Approximations of rectangles ------------------------------------------




template<class R>
Geometry::GridBlock<R>
Geometry::over_approximation(const Box<R>& bx, const Grid<R>& g) 
{
  IndexArray lower(bx.dimension());
  IndexArray upper(bx.dimension());
  for(size_type i=0; i!=bx.dimension(); ++i) {
    if(bx.lower_bound(i)==bx.upper_bound(i)) {
      ARIADNE_THROW(EmptyInterior,"GridBlock over_approximation(Box bx, Grid g)"," with bx="<<bx<<" (use outer_approximation(r,g) instead)");
    }
    lower[i]=g.subdivision_lower_index(i,bx.lower_bound(i));
    upper[i]=g.subdivision_upper_index(i,bx.upper_bound(i));
  }
  return GridBlock<R>(g,lower,upper);
}


template<class R>
Geometry::GridBlock<R>
Geometry::under_approximation(const Box<R>& bx, const Grid<R>& g) 
{
  IndexArray lower(bx.dimension());
  IndexArray upper(bx.dimension());
  for(size_type i=0; i!=bx.dimension(); ++i) {
    if(bx.lower_bound(i)==bx.upper_bound(i)) {
      ARIADNE_THROW(EmptyInterior,"GridBlock under_approximation(Box r, Grid g)"," with bx="<<bx<<" (use outer_approximation(r,g) instead)");
    }
    lower[i]=g.subdivision_upper_index(i,bx.lower_bound(i));
    upper[i]=g.subdivision_lower_index(i,bx.upper_bound(i));
  }
  return GridBlock<R>(g,lower,upper);
}


template<class R>
Geometry::GridBlock<R>
Geometry::outer_approximation(const Box<R>& bx, const Grid<R>& g) 
{
  if(bx.empty()) {
    return GridBlock<R>(g);
  }
  IndexArray lower(bx.dimension());
  IndexArray upper(bx.dimension());
  for(size_type i=0; i!=bx.dimension(); ++i) {
    lower[i]=g.subdivision_upper_index(i,bx.lower_bound(i))-1;
    upper[i]=g.subdivision_lower_index(i,bx.upper_bound(i))+1;
  }
  return GridBlock<R>(g,lower,upper);
}

template<class R>
Geometry::GridBlock<R>
Geometry::inner_approximation(const Box<R>& bx, const Grid<R>& g) 
{
  if(bx.empty()) {
    return GridBlock<R>(g);
  }
  IndexArray lower(bx.dimension());
  IndexArray upper(bx.dimension());
  for(size_type i=0; i!=bx.dimension(); ++i) {
    lower[i]=g.subdivision_lower_index(i,bx.lower_bound(i))+1;
    upper[i]=g.subdivision_upper_index(i,bx.upper_bound(i))-1;
  }
  return GridBlock<R>(g,lower,upper);
}



template<class R>
Geometry::GridBlock<R>
Geometry::outer_approximation(const Point< Numeric::Interval<R> >& ipt, const Grid<R>& g) 
{
  return outer_approximation(Box<R>(ipt),g);
}

template<class R>
Geometry::GridBlock<R>
Geometry::inner_approximation(const Point< Numeric::Interval<R> >& ipt, const Grid<R>& g) 
{
  return inner_approximation(Box<R>(ipt),g);
}

template<class R>
Geometry::GridCellListSet<R>
Geometry::outer_approximation(const Rectangle<R>& r, const Grid<R>& g) 
{
  GridCellListSet<R> gcls(g);
  gcls.adjoin_outer_approximation(Box<R>(r));
  return gcls;
}


template<class R>
Geometry::GridCellListSet<R>
Geometry::inner_approximation(const Rectangle<R>& r, const Grid<R>& g) 
{
  GridCellListSet<R> gcls(g);
  gcls.adjoin_inner_approximation(Box<R>(r));
  return gcls;
}

template<class R>
Geometry::GridCellListSet<R>
Geometry::outer_approximation(const Polyhedron<R>& ph, const Grid<R>& g) 
{
  return ::outer_approximation_of_basic_set(ph,g);
}

template<class R>
Geometry::GridCellListSet<R>
Geometry::inner_approximation(const Polyhedron<R>& ph, const Grid<R>& g) 
{
  return ::inner_approximation_of_basic_set(ph,g);
}


template<class R>
Geometry::GridCellListSet<R>
Geometry::outer_approximation(const Polytope<R>& p, const Grid<R>& g) 
{
  return ::outer_approximation_of_basic_set(p,g);
}

template<class R>
Geometry::GridCellListSet<R>
Geometry::inner_approximation(const Polytope<R>& p, const Grid<R>& g) 
{
  return ::inner_approximation_of_basic_set(p,g);
}


template<class R>
Geometry::GridCellListSet<R>
Geometry::outer_approximation(const Zonotope<R>& z, const Grid<R>& g) 
{
  return ::outer_approximation_of_basic_set(z,g);
}

template<class R>
Geometry::GridCellListSet<R>
Geometry::inner_approximation(const Zonotope<R>& z, const Grid<R>& g) 
{
  return ::inner_approximation_of_basic_set(z,g);
}


template<class R>
Geometry::GridCellListSet<R>
Geometry::outer_approximation(const BoxListSet<R>& bxls, const Grid<R>& g) 
{
  return ::outer_approximation_of_denotable_set(bxls,g);
}

template<class R>
Geometry::GridMaskSet<R>
Geometry::outer_approximation(const BoxListSet<R>& bxls, const FiniteGrid<R>& fg) 
{
  return ::outer_approximation_of_denotable_set(bxls,fg);
}

template<class R>
Geometry::GridCellListSet<R>
Geometry::outer_approximation(const SetInterface<R>& set, const Grid<R>& g) 
{
  ARIADNE_LOG(4,"GridCellListSet outer_approximation(SetInterface, Grid)\n");
  GridCellListSet<R> result(g);
  ARIADNE_CHECK_EQUAL_DIMENSIONS(set,g,"outer_approximation(SetInterface<R>,Grid<R>)");
  ARIADNE_ASSERT(set.bounding_box().bounded());
  ARIADNE_LOG(4,"bb="<<set.bounding_box()<<"\n");

  const GridBlock<R> gb=outer_approximation(set.bounding_box(),g);
  ARIADNE_LOG(4,"gb="<<gb<<"\n");
  Box<R> r(g.dimension());
  uint i=0;
  for(typename GridBlock<R>::const_iterator iter=gb.begin(); iter!=gb.end(); ++iter) {
    ARIADNE_LOG(4,*iter<<"\n");
    r=*iter;
    if(!bool(set.disjoint(r))) {
      result.adjoin(*iter);
    }
  }
  return result;
}


template<class R>
Geometry::GridCellListSet<R>
Geometry::inner_approximation(const SetInterface<R>& set, const Grid<R>& g) 
{
  ARIADNE_LOG(4,"GridCellListSet inner_approximation(SetInterface, Grid)\n");
  GridCellListSet<R> result(g);
  ARIADNE_CHECK_EQUAL_DIMENSIONS(set,g,"inner_approximation(SetInterface<R>,Grid<R>)\n");
  
  const GridBlock<R> gb=outer_approximation(set.bounding_box(),g);
  Box<R> r(g.dimension());
  for(typename GridBlock<R>::const_iterator iter=gb.begin(); iter!=gb.end(); ++iter) {
    r=*iter;
    if((set.superset(r))) {
      result.adjoin(*iter);
    }
  }
  return result;
}


template<class R>
Geometry::GridCellListSet<R>
Geometry::fuzzy_outer_approximation(const Zonotope<R>& z, const Grid<R>& g) 
{
  return ::fuzzy_outer_approximation_of_zonotope(z,g);
}

template<class R>
Geometry::BoxListSet<R>
Geometry::lower_approximation(const SetInterface<R>& s, const Grid<R>& g) 
{
  ARIADNE_LOG(4,"ListSet<Rectangle> lower_approximation(SetInterface s, Grid fg)\n");
  FiniteGrid<R> fg(g,s.bounding_box());
  return lower_approximation(s,fg);
}

template<class R>
Geometry::BoxListSet<R>
Geometry::point_approximation(const SetInterface<R>& s, const Grid<R>& g) 
{
  ARIADNE_LOG(4,"ListSet<Rectangle> point_approximation(SetInterface s, Grid fg)\n");
  FiniteGrid<R> fg(g,s.bounding_box());
  return point_approximation(s,fg);
}


template<class R, class BS>
Geometry::GridMaskSet<R>
Geometry::outer_approximation(const ListSet<BS>& ls, const FiniteGrid<R>& fg) 
{
  GridMaskSet<R> result(fg);
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls,fg,"outer_approximation(ListSet<BS> ls, FiniteGrid<R> fg)\n");
  
  for(typename ListSet<BS>::const_iterator bs_iter=ls.begin(); bs_iter!=ls.end(); ++bs_iter) {
    result.adjoin_outer_approximation(*bs_iter);
  }
  return result;
}

template<class R, class BS>
Geometry::GridCellListSet<R>
Geometry::outer_approximation(const ListSet<BS>& ls, const Grid<R>& g) 
{
  GridCellListSet<R> result(g);
  ARIADNE_CHECK_EQUAL_DIMENSIONS(ls,g,"outer_approximation(ListSet<BS> ls, Grid<R> g)\n");
  
  for(typename ListSet<BS>::const_iterator bs_iter=ls.begin(); bs_iter!=ls.end(); ++bs_iter) {
    result.adjoin_outer_approximation(*bs_iter);
  }
  return result;
}


template<class R>
Geometry::GridMaskSet<R>
Geometry::outer_approximation(const SetInterface<R>& set, const FiniteGrid<R>& fg) 
{
  ARIADNE_LOG(4,"GridMaskSet outer_approximation(SetInterface, FiniteGrid)\n");
  GridMaskSet<R> result(fg);
  ARIADNE_CHECK_EQUAL_DIMENSIONS(set,fg,"outer_approximation(PartitionTreeSet<R>,FiniteGrid<R>)");
  
  const Grid<R>& g=fg.grid();
  const Combinatoric::LatticeBlock& lb=fg.lattice_block();
  
  for(typename Combinatoric::LatticeBlock::const_iterator iter=lb.begin(); iter!=lb.end(); ++iter) {
    GridCell<R> gc(g,*iter);
    if(!bool(set.disjoint(Box<R>(gc)))) {
      result.adjoin(gc);
    }
  }
  return result;
}




template<class R>
Geometry::GridMaskSet<R>
Geometry::inner_approximation(const SetInterface<R>& s, const FiniteGrid<R>& fg) 
{
  ARIADNE_LOG(4,"GridMaskSet inner_approximation(SetInterface s, FiniteGrid fg)\n");
  GridMaskSet<R> result(fg);
  GridBlock<R> gb(fg.grid(),fg.lattice_block());
  Box<R> r;
  ARIADNE_LOG(5,"  testing cells in "<<gb<<"\n");
  for(typename GridBlock<R>::const_iterator gb_iter=gb.begin();
      gb_iter!=gb.end(); ++gb_iter)
  {
    r=*gb_iter;
    if(s.superset(r)) {
      result.adjoin(*gb_iter);
    }
    ARIADNE_LOG(6,"s.superset("<<r<<")="<<s.superset(r)<<"\n");
  }
  return result;
}





template<class R>  
Geometry::BoxListSet<R>
Geometry::lower_approximation(const SetInterface<R>& s, const FiniteGrid<R>& fg) 
{
  ARIADNE_LOG(4,"ListSet<Box> lower_approximation(SetInterface s, FiniteGrid fg)\n"); 
  BoxListSet<R> result;
  GridBlock<R> gb(fg.grid(),fg.lattice_block());
  Box<R> r;
  Box<R> nr;
  ARIADNE_LOG(5,"  testing cells in "<<gb<<"\n");
  for(typename GridBlock<R>::const_iterator gb_iter=gb.begin();
      gb_iter!=gb.end(); ++gb_iter)
  {
    r=*gb_iter;
    if(!bool(s.disjoint(r))) {
      nr=gb_iter->neighbourhood();
      if(s.intersects(nr)) {
        result.adjoin(nr);
      }
    }
    ARIADNE_LOG(6,"s.interscts("<<r<<")="<<s.intersects(r)<<"\n");
  }
  return result;
}


template<class R>  
Geometry::BoxListSet<R>
Geometry::point_approximation(const SetInterface<R>& s, const FiniteGrid<R>& fg) 
{
  typedef Numeric::Interval<R> I;
  ARIADNE_LOG(4,"ListSet<Box> point_approximation(SetInterface s, FiniteGrid fg)\n"); 
  BoxListSet<R> result;
  GridBlock<R> gb(fg.grid(),fg.lattice_block());
  Box<R> r;
  Box<R> cr;
  Point<I> cpt;
  ARIADNE_LOG(5,"  testing cells in "<<gb<<"\n");
  for(typename GridBlock<R>::const_iterator gb_iter=gb.begin();
      gb_iter!=gb.end(); ++gb_iter)
  {
    r=*gb_iter;
    if(s.superset(r)) {
      cpt=r.centre();
      cr=static_cast< Box<R> >(cpt);
      result.adjoin(cr);
    }
    ARIADNE_LOG(6,"s.superset("<<r<<")="<<s.superset(r)<<"\n");
  }
  return result;
}














}
