/****************************************************************************
 *            epsstream.template.h
 *
 *  Copyright  2005-8  Alberto Casagrande, Pieter Collins
 *
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

namespace Ariadne {



template<class R> 
epsstream&
operator<<(epsstream& eps, const Point<R>& pt) 
{
  Point2d dpt=eps.projection_map()(pt);
  eps.draw(dpt);
  return eps;
}

template<class R> 
epsstream&
operator<<(epsstream& eps, const Box<R>& r) 
{
  Rectangle2d dr=eps.projection_map()(r);
  eps.draw(dr);
  return eps;
}

template<class R> 
epsstream&
operator<<(epsstream& eps, const Rectangle<R>& r) 
{
  Rectangle2d dr=eps.projection_map()(r);
  eps.draw(dr);
  return eps;
}

template<class R> 
epsstream&
operator<<(epsstream& eps, const RectangularSet<R>& rs)
{
  return eps << Box<R>(rs);
}




template<class R> 
epsstream&
operator<<(epsstream& eps, const Zonotope<R>& z)
{ 
  Zonotope<Rational> qz(error_free_over_approximation(z));
  Polygon2d vertices=eps.projection_map()(qz);      
  eps.draw(vertices);
  return eps;
}


template<class R> 
epsstream&
operator<<(epsstream& eps, const Polytope<R>& p)
{
  if(p.bounded()) {
		eps.draw(eps.projection_map()(p));
	}
  return eps;
}

template<class R> 
epsstream&
operator<<(epsstream& eps, const Polyhedron<R>& p)
{
	if(p.bounded()) {
		return eps << Polytope<Rational>(p);
	} else {
		return eps;
	}
}

template<class R> 
epsstream&
operator<<(epsstream& eps, const PolyhedralSet<R>& ps)
{
  return eps << Polyhedron<R>(ps);
}


template<class R> 
epsstream&
operator<<(epsstream& eps, const BoxListSet<R>& ds)
{
  typedef typename BoxListSet<R>::const_iterator const_iterator;
  if(eps.fill_style) {
    // draw without lines
    bool line_style=eps.line_style; 
    eps.line_style=false;
    for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
      eps << *set_iter;
    }
    eps.line_style=line_style;
  }
  if(eps.line_style) {
    bool fill_style=eps.fill_style; 
    eps.fill_style=false;
    for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
      eps << *set_iter;
    }
    eps.fill_style=fill_style;
  }
  return eps;
}

template<class BS> 
epsstream&
operator<<(epsstream& eps, const ListSet<BS>& ds)
{
  typedef typename ListSet<BS>::const_iterator const_iterator;
  if(eps.fill_style) {
    // draw without lines
    bool line_style=eps.line_style; 
    eps.line_style=false;
    for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
      eps << *set_iter;
    }
    eps.line_style=line_style;
  }
  if(eps.line_style) {
    bool fill_style=eps.fill_style; 
    eps.fill_style=false;
    for(const_iterator set_iter=ds.begin(); set_iter!=ds.end(); ++set_iter) {
      eps << *set_iter;
    }
    eps.fill_style=fill_style;
  }
  return eps;
}



template<class R> 
epsstream&
operator<<(epsstream& eps, const GridCell<R>& bs)
{
  return eps << Box<R>(bs);
}


template<class R> 
epsstream&
operator<<(epsstream& eps, const GridBlock<R>& bs)
{
  return eps << Box<R>(bs);
}


template<class R> 
epsstream&
operator<<(epsstream& eps, const GridCellListSet<R>& ds)
{
  std::vector<Rectangle2d> rl(ds.size());
  std::transform(ds.begin(),ds.end(),rl.begin(),eps.projection_map());
  eps.draw(rl);
  return eps;
}


template<class R> 
epsstream&
operator<<(epsstream& eps, const GridMaskSet<R>& ds)
{
  std::vector<Rectangle2d> rl(ds.size());
  std::transform(ds.begin(),ds.end(),rl.begin(),eps.projection_map());
  eps.draw(rl);
  return eps;
}



template<class R> 
epsstream&
operator<<(epsstream& eps, const PartitionTreeSet<R>& ds)
{
  std::vector<Rectangle2d> rl(ds.size());
  std::transform(ds.begin(),ds.end(),rl.begin(),eps.projection_map());
  eps.draw(rl);
  return eps;
}

template<class R> 
epsstream&
operator<<(epsstream& eps, const SetInterface< Box<R> >& set)
{
  
  typedef Interval<R> I;
  
  if(dynamic_cast<const RectangularSet<R>*>(&set)) {
    return eps << dynamic_cast<const RectangularSet<R>&>(set);
  } else if(dynamic_cast<const PolyhedralSet<R>*>(&set)) {
    return eps << dynamic_cast<const PolyhedralSet<R>&>(set);
  } else if(dynamic_cast<const ListSet< Zonotope<R> >*>(&set)) {
    return eps << dynamic_cast<const ListSet< Zonotope<R> >&>(set);
  } else if(dynamic_cast<const GridCellListSet<R>*>(&set)) {
    return eps << dynamic_cast<const GridCellListSet<R>&>(set);
  } else if(dynamic_cast<const GridMaskSet<R>*>(&set)) {
    return eps << dynamic_cast<const GridMaskSet<R>&>(set);
  } else if(dynamic_cast<const PartitionTreeSet<R>*>(&set)) {
    return eps << dynamic_cast<const PartitionTreeSet<R>&>(set);
  }  else {
    Box<R> bb;
    try {
      bb=set.bounding_box();
    } 
    catch(UnboundedSet& e) {
      if(set.space()==EuclideanSpace(2)) {
        Rectangle2d bbox=eps.bounding_box();
        bb=Box<R>(2);
        bb.set_lower_bound(0,bbox.lower_bound(0));
        bb.set_upper_bound(0,bbox.upper_bound(0));
        bb.set_lower_bound(1,bbox.lower_bound(1));
        bb.set_upper_bound(1,bbox.upper_bound(1));
      } else {
        throw e;
      }
    }
    PartitionScheme<R> ps(bb);
    int depth=16;
    PartitionTreeSet<R> pts=outer_approximation(set,ps,depth);
    return eps << pts;
  }
}

template<class R> 
epsstream&
operator<<(epsstream& eps, const FiniteGrid<R>& fg)
{
  bool fill_style=eps.fill_style;
  if(fill_style) { eps.fill_style=false; }
  GridCellListSet<R> gcls(fg.grid());
  gcls.adjoin(GridBlock<R>(fg.grid(),fg.lattice_block()));
  eps << gcls;
  if(fill_style) { eps.fill_style=true; }
  return eps;
}


template<class R> 
epsstream&
operator<<(epsstream& eps, const PartitionTree<R>& pt)
{
  bool fill_style=eps.fill_style;
  if(fill_style) { eps.fill_style=false; }
  for(typename PartitionTree<R>::const_iterator iter = pt.begin(); iter!=pt.end(); ++iter) {
    eps << Box<R>(*iter);
  }
  if(fill_style) { eps.fill_style=true; }
  return eps;
}

}
