/***************************************************************************
 *            grid_tree_set.code.h
 *
 *
 *  Copyright  2005-7  Pieter Collins
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


#include "grid_tree_set.h"

#include "base/stlio.h"

#include "geometry/box.h"
#include "geometry/list_set.h"
#include "geometry/grid_mask_set.h"
#include "output/logging.h"

#include <vector>

namespace Ariadne {


    
    


template<class R>
Geometry::GridTreeSet<R>::GridTreeSet(const GridMaskSet<R>& gms) 
  : _unit_box(gms.bounding_box()),
    _subdivision_set(gms._lattice_set)
{ }



template<class R>
Geometry::GridTreeSet<R>*
Geometry::GridTreeSet<R>::clone() const
{
  return new GridTreeSet<R>(*this);
}


template<class R>
tribool
Geometry::GridTreeSet<R>::contains(const Point<R>& pt) const
{
  return Geometry::contains(*this,pt);
}


template<class R>
tribool
Geometry::GridTreeSet<R>::intersects(const Box<R>& r) const
{
  return !Geometry::disjoint(*this,r);
}


template<class R>
tribool
Geometry::GridTreeSet<R>::disjoint(const Box<R>& r) const
{
  return Geometry::disjoint(*this,r);
}


template<class R>
tribool
Geometry::GridTreeSet<R>::superset(const Box<R>& r) const
{
  return Geometry::subset(r,*this);
}


template<class R>
tribool
Geometry::GridTreeSet<R>::subset(const Box<R>& r) const
{
  return Geometry::subset(*this,r);
}







template<class R>
Geometry::GridTreeSet<R>::operator ListSet< Box<R> >() const 
{
  ListSet< Box<R> > res(this->dimension());
  for(const_iterator iter=begin(); iter!=end(); ++iter) {
    res.push_back(Box<R>(*iter));
  }
  return res;
}



template<class R>
void
Geometry::GridTreeSet<R>::_instantiate_geometry_operators()
{
  uint d=0;
  Box<R>* r=0;
  Zonotope<R>* z=0;
  SetInterface<R>* s=0;
  GridMaskSet<R>* gms=0;
  GridScheme<R>* ps=0;
  GridTreeCell<R>* ptc=0;
  GridTreeSet<R>* pts=0;
  
  *ptc=over_approximation(*r,*ps);
  
  *pts=outer_approximation(*r,*ps,d);
  *pts=inner_approximation(*r,*ps,d);
  
  *pts=outer_approximation(*z,*ps,d);
  *pts=inner_approximation(*z,*ps,d);
  
  *pts=outer_approximation(*gms,*ps,d);
  *pts=inner_approximation(*gms,*ps,d);
  
  *pts=outer_approximation(*s,*ps,d);
  *pts=inner_approximation(*s,*ps,d);
}


template<class R>
tribool
Geometry::contains(const GridTreeSet<R>& pts, const Point<R>& pt)
{
  return subset(Box<R>(pt),pts);
}


template<class R>
tribool
Geometry::disjoint(const GridTreeSet<R>& pts, const Box<R>& r)
{
  tribool result=true;
  Box<R> cell(pts.dimension());
  for(typename GridTreeSet<R>::const_iterator pts_iter=pts.begin();
      pts_iter!=pts.end(); ++pts_iter)
    {
      cell=*pts_iter;
      result = result && Geometry::disjoint(cell,r);
      if(result==false) {
        return result;
      }
    }
  return result;
}


template<class R>
tribool
Geometry::subset(const GridTreeSet<R>& pts, const Box<R>& r)
{
  tribool result=true;
  Box<R> cell(pts.dimension());
  for(typename GridTreeSet<R>::const_iterator pts_iter=pts.begin();
      pts_iter!=pts.end(); ++pts_iter)
    {
      cell=*pts_iter;
      result = result && Geometry::subset(cell,r);
      if(result==false) {
        return result;
      }
    }
  return result;
}


template<class R>
tribool
Geometry::subset(const Box<R>& r, const GridTreeSet<R>& pts)
{
  return subset(r,ListSet< Box<R> >(pts));
}




template<class R>
Geometry::GridTreeCell<R>
Geometry::over_approximation(const Box<R>& r, const GridScheme<R>& ps)
{
  if(possibly(r.empty())) {
    ARIADNE_THROW(EmptyInterior,"GridTreeCell<R> over_approximation(const Box<R>& r, const GridScheme<R>& ps)"," with r="<<r);
  }
  
  if(!subset(r,ps.unit_box())) {
    ARIADNE_THROW(std::runtime_error,"ParttitionTreeCell over_approximation(Box r, GridScheme ps, uint d)"," with r="<<r<<", ps="<<ps<<": r is not a subset of ps.unit_box()");
  }
  
  const Box<R>& unit_box=ps.unit_box();
  const Combinatoric::SubdivisionSequence& subdivisions=ps.subdivisions();
  
  Combinatoric::BinaryWord word;
  Combinatoric::BinaryWord new_word;
  
  while(true) {
    new_word=word;
    new_word.push_back(Combinatoric::BinaryTree::leaf);
    if(!subset(r,Box<R>(GridTreeCell<R>(unit_box,subdivisions,new_word)))) {
      new_word.pop_back();
      word.push_back(Combinatoric::BinaryTree::right);
      if(!subset(r,Box<R>(GridTreeCell<R>(unit_box,subdivisions,word)))) {
        break;
      }
    }
  }
  return GridTreeCell<R>(unit_box,subdivisions,word);
}


template<class R, class S>
Geometry::GridTreeSet<R>
Geometry::outer_approximation(const S& s, const Grid<R>& g, const uint depth)
{
  ARIADNE_LOG(4,"outer_approximation(S set, GridScheme ps, uint depth");
  const Box<R>& bounding_box=ps.unit_box();
  const Combinatoric::SubdivisionSequence& subdivisions=ps.subdivisions();
  std::vector<bool> tree;
  std::vector<bool> mask;
  
  Combinatoric::BinaryWord word;
  
  do {
    Box<R> cell=Box<R>(GridTreeCell<R>(bounding_box,subdivisions,word));
    if(word.size()==depth+1 || subset(cell,s)) {
      tree.push_back(Combinatoric::BinaryTree::leaf);
      mask.push_back(true);
      Combinatoric::BinaryTree::advance(word);
    }  
    else if(disjoint(cell,s)) {
      tree.push_back(Combinatoric::BinaryTree::leaf);
      mask.push_back(false);
      Combinatoric::BinaryTree::advance(word);
    }
    else {
      tree.push_back(Combinatoric::BinaryTree::branch);
      word.push_back(Combinatoric::BinaryTree::left);
    }
  } while(!word.empty());
  
  return GridTreeSet<R>(bounding_box,subdivisions,Combinatoric::BinaryTree(tree),BooleanArray(mask));
}

template<class R, class S>
Geometry::GridTreeSet<R>
Geometry::inner_approximation(const S& s, const GridScheme<R>& ps, const uint depth)
{
  const Box<R>& unit_box=ps.unit_box();
  const Combinatoric::SubdivisionSequence& subdivisions=ps.subdivisions();
  std::vector<bool> tree;
  std::vector<bool> mask;
  
  Combinatoric::BinaryWord word;
  
  do {
    Box<R> cell=Box<R>(GridTreeCell<R>(unit_box,subdivisions,word));
    if(word.size()==depth+1 || disjoint(cell,s)) {
      tree.push_back(Combinatoric::BinaryTree::leaf);
      mask.push_back(false);
      Combinatoric::BinaryTree::advance(word);
    }  
    else if(subset(cell,s)) {
      tree.push_back(Combinatoric::BinaryTree::leaf);
      mask.push_back(true);
      Combinatoric::BinaryTree::advance(word);
    }
    else {
      tree.push_back(Combinatoric::BinaryTree::branch);
      word.push_back(Combinatoric::BinaryTree::left);
    }
  } while(!word.empty());
  
  return GridTreeSet<R>(unit_box,subdivisions,Combinatoric::BinaryTree(tree),BooleanArray(mask));
}


template<class R, class S>
Geometry::GridTreeSet<R>
Geometry::under_approximation(const S& s, const GridScheme<R>& ps, const uint depth)
{
  const Box<R>& unit_box=ps.unit_box();
  const Combinatoric::SubdivisionSequence& subdivisions=ps.subdivisions();
  std::vector<bool> tree;
  std::vector<bool> mask;
  
  Combinatoric::BinaryWord word;
  
  do {
    Box<R> cell=Box<R>(GridTreeCell<R>(unit_box,subdivisions,word));
    if(word.size()==depth+1 || disjoint(cell,s)) {
      tree.push_back(Combinatoric::BinaryTree::leaf);
      mask.push_back(false);
      Combinatoric::BinaryTree::advance(word);
    }  
    else if(subset(cell,s)) {
      tree.push_back(Combinatoric::BinaryTree::leaf);
      mask.push_back(true);
      Combinatoric::BinaryTree::advance(word);
    }
    else {
      tree.push_back(Combinatoric::BinaryTree::branch);
      word.push_back(Combinatoric::BinaryTree::left);
    }
  } while(!word.empty());
  
  return GridTreeSet<R>(unit_box,subdivisions,Combinatoric::BinaryTree(tree),BooleanArray(mask));
}



template<class R>
std::ostream&
Geometry::GridScheme<R>::write(std::ostream& os) const
{
  os << "GridScheme<" << Numeric::name<R>() << ">(\n";
  os << "  unit_box=" << this->unit_box() << ",\n";
  os << "  subdivision_coordinates=" << this->subdivisions() << "\n";
  os << ")\n";
  return os;
}

template<class R>
std::ostream&
Geometry::GridTreeCell<R>::write(std::ostream& os) const
{
  os << "GridTreeCell<" << Numeric::name<R>() << ">(\n";
  os << "  bounds=" << this->subdivision_cell() << ",\n";
  os << "  rectangle=" << Box<R>(*this) << "\n";
  os << ")\n";
  return os;
}

template<class R>
std::ostream&
Geometry::GridTree<R>::write(std::ostream& os) const
{
  os << "GridTree<" << Numeric::name<R>() << ">(\n";
  os << "  unit_box=" << this->unit_box() << ",\n";
  os << "  subdivisions=" << this->subdivisions() << "\n";
  os << "  words="; Base::write_sequence(os, this->binary_tree().begin(), this->binary_tree().end()); os << ",\n";
  os << "  blocks=["; Base::write_sequence(os,  this->subdivision_tree().begin(), this->subdivision_tree().end()); os << ",\n";
  os << "  cells=["; Base::write_sequence(os,  this->begin(), this->end()); os << ",\n";
  os << ")\n";
  return os;
}

template<class R>
std::ostream&
Geometry::GridTreeSet<R>::write(std::ostream& os) const
{
  os << "GridTreeSet<" << Numeric::name<R>() << ">(\n";
  os << "  unit_box=" << this->unit_box() << ",\n";
  os << "  subdivisions=" << this->subdivisions() << ",\n";
  os << "  tree=" << this->binary_tree() << ",\n";
  os << "  mask=" << this->mask() << ",\n";
  os << "  words="; Base::write_sequence(os, this->subdivision_set().words().begin(), this->subdivision_set().words().end()); os << ",\n";
  os << "  blocks=["; Base::write_sequence(os,  this->subdivision_set().begin(), this->subdivision_set().end()); os << ",\n";
  os << "  cells=["; Base::write_sequence(os,  this->begin(), this->end()); os << ",\n";
  os << ")\n";
  return os;
}

template<class R>
std::ostream&
Geometry::GridTreeSet<R>::summarize(std::ostream& os) const
{
  os << "GridTreeSet<" << Numeric::name<R>() << ">(";
  os << " unit_box=" << this->unit_box() << ",";
  os << " size=" << this->size() << " )";
  return os;
}


} // namespace Ariadne
