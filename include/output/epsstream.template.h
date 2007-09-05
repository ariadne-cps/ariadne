/****************************************************************************
 *            epsstream.template.h
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
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

namespace Ariadne {



template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::Point<R>& pt) 
{
  Point2d dpt=eps.projection_map()(pt);
  eps.draw(dpt);
  return eps;
}

template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::Rectangle<R>& r) 
{
  Rectangle2d dr=eps.projection_map()(r);
  eps.draw(dr);
  return eps;
}

template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::RectangularSet<R>& rs)
{
  return eps << Geometry::Rectangle<R>(rs);
}


template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::Zonotope< Numeric::Interval<R>,Numeric::Interval<R> >& iz)
{ 
  Geometry::Zonotope<Numeric::Interval<R>,R> ez=Geometry::over_approximation(iz);
  Geometry::Zonotope<R> z=Geometry::over_approximation(ez);
  Geometry::Zonotope<Numeric::Rational> qz(z);
  Polygon2d vertices=eps.projection_map()(qz);      
  eps.draw(vertices);
  return eps;
}

template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::Zonotope<Numeric::Interval<R>,R>& ez)
{ 
  Geometry::Zonotope<R> z=Geometry::over_approximation(ez);
  Geometry::Zonotope<Numeric::Rational> qz(z);
  Polygon2d vertices=eps.projection_map()(qz);      
  eps.draw(vertices);
  return eps;
}

template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::Zonotope<R>& z)
{ 
  Geometry::Zonotope<Numeric::Rational> qz(z);
  Polygon2d vertices=eps.projection_map()(qz);      
  eps.draw(vertices);
  return eps;
}

template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::Parallelotope<R>& p)
{
  const Geometry::Zonotope<R>& z=p;
  return eps << z;
}

template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::Polytope<R>& p)
{
  eps.draw(eps.projection_map()(p.vertices()));
  return eps;
}

template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::Polyhedron<R>& p)
{
  return eps << Geometry::Polytope<Numeric::Rational>(p);
}

template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::PolyhedralSet<R>& ps)
{
  return eps << Geometry::Polyhedron<R>(ps);
}


template<class BS> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::ListSet<BS>& ds)
{
  typedef typename Geometry::ListSet<BS>::const_iterator const_iterator;
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
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::GridCell<R>& bs)
{
  return eps << Geometry::Rectangle<R>(bs);
}


template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::GridBlock<R>& bs)
{
  return eps << Geometry::Rectangle<R>(bs);
}


template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::GridCellListSet<R>& ds)
{
  std::vector<Rectangle2d> rl(ds.size());
  std::transform(ds.begin(),ds.end(),rl.begin(),eps.projection_map());
  eps.draw(rl);
  return eps;
}


template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::GridMaskSet<R>& ds)
{
  std::vector<Rectangle2d> rl(ds.size());
  std::transform(ds.begin(),ds.end(),rl.begin(),eps.projection_map());
  eps.draw(rl);
  return eps;
}



template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::PartitionTreeSet<R>& ds)
{
  std::vector<Rectangle2d> rl(ds.size());
  std::transform(ds.begin(),ds.end(),rl.begin(),eps.projection_map());
  eps.draw(rl);
  return eps;
}

template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::SetInterface<R>& set)
{
  using namespace Geometry;
  typedef Numeric::Interval<R> I;
  
  if(dynamic_cast<const RectangularSet<R>*>(&set)) {
    return eps << dynamic_cast<const RectangularSet<R>&>(set);
  } else if(dynamic_cast<const PolyhedralSet<R>*>(&set)) {
    return eps << dynamic_cast<const PolyhedralSet<R>&>(set);
  } else if(dynamic_cast<const ListSet< Rectangle<R> >*>(&set)) {
    return eps << dynamic_cast<const ListSet< Rectangle<R> >&>(set);
  } else if(dynamic_cast<const ListSet< Zonotope<R,R> >*>(&set)) {
    return eps << dynamic_cast<const ListSet< Zonotope<R,R> >&>(set);
  } else if(dynamic_cast<const ListSet< Zonotope<I,R> >*>(&set)) {
    return eps << dynamic_cast<const ListSet< Zonotope<I,R> >&>(set);
  } else if(dynamic_cast<const ListSet< Zonotope<I,I> >*>(&set)) {
    return eps << dynamic_cast<const ListSet< Zonotope<I,I> >&>(set);
  } else if(dynamic_cast<const GridCellListSet<R>*>(&set)) {
    return eps << dynamic_cast<const GridCellListSet<R>&>(set);
  } else if(dynamic_cast<const GridMaskSet<R>*>(&set)) {
    return eps << dynamic_cast<const GridMaskSet<R>&>(set);
  } else if(dynamic_cast<const PartitionTreeSet<R>*>(&set)) {
    return eps << dynamic_cast<const PartitionTreeSet<R>&>(set);
  }  else {
    Rectangle<R> bb;
    try {
      bb=set.bounding_box();
    } 
    catch(Geometry::UnboundedSet& e) {
      if(set.dimension()==2) {
        Rectangle2d bbox=eps.bounding_box();
        bb=Geometry::Rectangle<R>(2);
        bb.set_lower_bound(0,bbox.lower_bound(0));
        bb.set_upper_bound(0,bbox.upper_bound(0));
        bb.set_lower_bound(1,bbox.lower_bound(1));
        bb.set_upper_bound(1,bbox.upper_bound(1));
      } else {
        throw e;
      }
    }
    Geometry::PartitionScheme<R> ps(bb);
    int depth=16;
    Geometry::PartitionTreeSet<R> pts=Geometry::outer_approximation(set,ps,depth);
    return eps << pts;
  }
}

template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::FiniteGrid<R>& fg)
{
  bool fill_style=eps.fill_style;
  if(fill_style) { eps.fill_style=false; }
  Geometry::GridCellListSet<R> gcls(fg.grid());
  gcls.adjoin(Geometry::GridBlock<R>(fg.grid(),fg.lattice_block()));
  eps << gcls;
  if(fill_style) { eps.fill_style=true; }
  return eps;
}


template<class R> 
Output::epsstream&
Output::operator<<(epsstream& eps, const Geometry::PartitionTree<R>& pt)
{
  bool fill_style=eps.fill_style;
  if(fill_style) { eps.fill_style=false; }
  for(typename Geometry::PartitionTree<R>::const_iterator iter = pt.begin(); iter!=pt.end(); ++iter) {
    eps << Geometry::Rectangle<R>(*iter);
  }
  if(fill_style) { eps.fill_style=true; }
  return eps;
}

}
