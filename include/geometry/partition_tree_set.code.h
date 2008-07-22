/***************************************************************************
 *            partition_tree_set.code.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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


#include "partition_tree_set.h"

#include "base/stlio.h"

#include "geometry/euclidean_space.h"
#include "geometry/box.h"
#include "geometry/box_list_set.h"
#include "geometry/grid_cell_list_set.h"
#include "geometry/grid_mask_set.h"
#include "geometry/list_set.h"
#include "output/logging.h"

#include <vector>

namespace Ariadne {


    
    
template<class R>
PartitionScheme<R>::PartitionScheme(const Box<R>& bb) 
  : _unit_box(bb), 
    _subdivisions(bb.dimension()) 
{ }


// FIXME: Should use add_approx etc here.
// We use approximate values since the cell boundaries are defined
// by the lower bound an upper bound formulae (which must be the same given 
// the cell bound.
template<class R>
R
PartitionTreeCell<R>::lower_bound(dimension_type i) const 
{
  return add_approx(_unit_box.lower_bound(i),
                    mul_approx(_subdivision_cell.lower_bound(i),
                               sub_approx(_unit_box.upper_bound(i),_unit_box.lower_bound(i))));
}

template<class R>
R
PartitionTreeCell<R>::upper_bound(dimension_type i) const 
{
  return add_approx(_unit_box.lower_bound(i),
                    mul_approx(_subdivision_cell.upper_bound(i),
                               sub_approx(_unit_box.upper_bound(i),_unit_box.lower_bound(i))));
}



template<class R>
PartitionTreeSet<R>::PartitionTreeSet(const GridMaskSet<R>& gms) 
  : _unit_box(gms.bounding_box()),
    _subdivision_set(gms._lattice_set)
{ }



template<class R>
PartitionTreeSet<R>*
PartitionTreeSet<R>::clone() const
{
  return new PartitionTreeSet<R>(*this);
}


template<class R> inline
EuclideanSpace
PartitionTreeSet<R>::space() const 
{
  return EuclideanSpace(this->dimension()); 
}


template<class R>
tribool
PartitionTreeSet<R>::contains(const Point<R>& pt) const
{
  return Ariadne::contains(*this,pt);
}


template<class R>
tribool
PartitionTreeSet<R>::intersects(const Box<R>& r) const
{
  return !Ariadne::disjoint(*this,r);
}


template<class R>
tribool
PartitionTreeSet<R>::disjoint(const Box<R>& r) const
{
  return Ariadne::disjoint(*this,r);
}


template<class R>
tribool
PartitionTreeSet<R>::superset(const Box<R>& r) const
{
  return Ariadne::subset(r,*this);
}


template<class R>
tribool
PartitionTreeSet<R>::subset(const Box<R>& r) const
{
  return Ariadne::subset(*this,r);
}







template<class R>
PartitionTreeSet<R>::operator BoxListSet<R>() const 
{
  BoxListSet<R> res(this->dimension());
  for(const_iterator iter=begin(); iter!=end(); ++iter) {
    res.push_back(Box<R>(*iter));
  }
  return res;
}



template<class R>
void
PartitionTreeSet<R>::_instantiate_geometry_operators()
{
  uint d=0;
  Box<R>* bx=0;
  Box<R>* r=0;
  Zonotope<R>* z=0;
  SetInterface< Box<R> >* s=0;
  GridMaskSet<R>* gms=0;
  PartitionScheme<R>* ps=0;
  PartitionTreeCell<R>* ptc=0;
  PartitionTreeSet<R>* pts=0;
  
  *ptc=Ariadne::over_approximation(*bx,*ps);
  
  *pts=Ariadne::outer_approximation(*r,*ps,d);
  *pts=Ariadne::inner_approximation(*r,*ps,d);
  
  *pts=Ariadne::outer_approximation(*z,*ps,d);
  *pts=Ariadne::boundary_approximation(*z,*ps,d);
  *pts=Ariadne::inner_approximation(*z,*ps,d);
  
  *pts=Ariadne::outer_approximation(*gms,*ps,d);
  *pts=Ariadne::inner_approximation(*gms,*ps,d);
  
  *pts=Ariadne::outer_approximation(*s,*ps,d);
  *pts=Ariadne::inner_approximation(*s,*ps,d);
}


template<class R>
tribool
contains(const PartitionTreeSet<R>& pts, const Point<R>& pt)
{
  return subset(Box<R>(pt),pts);
}


template<class R>
tribool
disjoint(const PartitionTreeSet<R>& pts, const Box<R>& r)
{
  tribool result=true;
  Box<R> cell(pts.dimension());
  for(typename PartitionTreeSet<R>::const_iterator pts_iter=pts.begin();
      pts_iter!=pts.end(); ++pts_iter)
    {
      cell=*pts_iter;
      result = result && disjoint(r,cell);
      if(result==false) {
        return result;
      }
    }
  return result;
}


template<class R>
tribool
subset(const PartitionTreeSet<R>& pts, const Box<R>& r)
{
  tribool result=true;
  Box<R> cell(pts.dimension());
  for(typename PartitionTreeSet<R>::const_iterator pts_iter=pts.begin();
      pts_iter!=pts.end(); ++pts_iter)
    {
      cell=*pts_iter;
      result = result && subset(cell,r);
      if(result==false) {
        return result;
      }
    }
  return result;
}


template<class R>
tribool
subset(const Box<R>& r, const PartitionTreeSet<R>& pts)
{
  return subset(r,BoxListSet<R>(pts));
}




template<class R>
PartitionTreeCell<R>
over_approximation(const Box<R>& r, const PartitionScheme<R>& ps)
{
  if(possibly(r.empty())) {
    ARIADNE_THROW(EmptyInterior,"PartitionTreeCell<R> over_approximation(const Box<R>& r, const PartitionScheme<R>& ps)"," with r="<<r);
  }
  
  if(!subset(r,ps.unit_box())) {
    ARIADNE_THROW(std::runtime_error,"ParttitionTreeCell over_approximation(Box r, PartitionScheme ps, uint d)"," with r="<<r<<", ps="<<ps<<": r is not a subset of ps.unit_box()");
  }
  
  const Box<R>& unit_box=ps.unit_box();
  const SubdivisionSequence& subdivisions=ps.subdivisions();
  
  BinaryWord word;
  BinaryWord new_word;
  
  while(true) {
    new_word=word;
    new_word.push_back(BinaryTree::leaf);
    if(!subset(r,Box<R>(PartitionTreeCell<R>(unit_box,subdivisions,new_word)))) {
      new_word.pop_back();
      word.push_back(BinaryTree::right);
      if(!subset(r,Box<R>(PartitionTreeCell<R>(unit_box,subdivisions,word)))) {
        break;
      }
    }
  }
  return PartitionTreeCell<R>(unit_box,subdivisions,word);
}


template<class R, class S>
PartitionTreeSet<R>
outer_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
{
  ARIADNE_LOG(4,"outer_approximation(S set, PartitionScheme ps, uint depth");
  const Box<R>& bounding_box=ps.unit_box();
  const SubdivisionSequence& subdivisions=ps.subdivisions();
  std::vector<bool> tree;
  std::vector<bool> mask;
  
  BinaryWord word;
  
  do {
    Box<R> cell=Box<R>(PartitionTreeCell<R>(bounding_box,subdivisions,word));
    if(word.size()==depth+1 || superset(s,cell)) {
      tree.push_back(BinaryTree::leaf);
      mask.push_back(true);
      BinaryTree::advance(word);
    }  
    else if(disjoint(s,cell)) {
      tree.push_back(BinaryTree::leaf);
      mask.push_back(false);
      BinaryTree::advance(word);
    }
    else {
      tree.push_back(BinaryTree::branch);
      word.push_back(BinaryTree::left);
    }
  } while(!word.empty());
  
  return PartitionTreeSet<R>(bounding_box,subdivisions,BinaryTree(tree),BooleanArray(mask));
}

template<class R, class S>
PartitionTreeSet<R>
 boundary_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
{
  ARIADNE_LOG(4,"boundary_approximation(S set, PartitionScheme ps, uint depth");
  const Box<R>& bounding_box=ps.unit_box();
  const SubdivisionSequence& subdivisions=ps.subdivisions();
  std::vector<bool> tree;
  std::vector<bool> mask;
  
  BinaryWord word;
  
  do {
    Box<R> cell=Box<R>(PartitionTreeCell<R>(bounding_box,subdivisions,word));
    if(word.size()==depth+1) {
      tree.push_back(BinaryTree::leaf);
      mask.push_back(true);
      BinaryTree::advance(word);
    }  
    else if(disjoint(s,cell) || superset(s,cell)) {
      tree.push_back(BinaryTree::leaf);
      mask.push_back(false);
      BinaryTree::advance(word);
    }
    else {
      tree.push_back(BinaryTree::branch);
      word.push_back(BinaryTree::left);
    }
  } while(!word.empty());
  
  return PartitionTreeSet<R>(bounding_box,subdivisions,BinaryTree(tree),BooleanArray(mask));
}

template<class R, class S>
PartitionTreeSet<R>
inner_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
{
  const Box<R>& unit_box=ps.unit_box();
  const SubdivisionSequence& subdivisions=ps.subdivisions();
  std::vector<bool> tree;
  std::vector<bool> mask;
  
  BinaryWord word;
  
  do {
    Box<R> cell=Box<R>(PartitionTreeCell<R>(unit_box,subdivisions,word));
    if(word.size()==depth+1 || disjoint(s,cell)) {
      tree.push_back(BinaryTree::leaf);
      mask.push_back(false);
      BinaryTree::advance(word);
    }  
    else if(superset(s,cell)) {
      tree.push_back(BinaryTree::leaf);
      mask.push_back(true);
      BinaryTree::advance(word);
    }
    else {
      tree.push_back(BinaryTree::branch);
      word.push_back(BinaryTree::left);
    }
  } while(!word.empty());
  
  return PartitionTreeSet<R>(unit_box,subdivisions,BinaryTree(tree),BooleanArray(mask));
}


template<class R, class S>
PartitionTreeSet<R>
under_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
{
  const Box<R>& unit_box=ps.unit_box();
  const SubdivisionSequence& subdivisions=ps.subdivisions();
  std::vector<bool> tree;
  std::vector<bool> mask;
  
  BinaryWord word;
  
  do {
    Box<R> cell=Box<R>(PartitionTreeCell<R>(unit_box,subdivisions,word));
    if(word.size()==depth+1 || disjoint(s,cell)) {
      tree.push_back(BinaryTree::leaf);
      mask.push_back(false);
      BinaryTree::advance(word);
    }  
    else if(superset(s,cell)) {
      tree.push_back(BinaryTree::leaf);
      mask.push_back(true);
      BinaryTree::advance(word);
    }
    else {
      tree.push_back(BinaryTree::branch);
      word.push_back(BinaryTree::left);
    }
  } while(!word.empty());
  
  return PartitionTreeSet<R>(unit_box,subdivisions,BinaryTree(tree),BooleanArray(mask));
}



template<class R>
std::ostream&
PartitionScheme<R>::write(std::ostream& os) const
{
  os << "PartitionScheme<" << name<R>() << ">(\n";
  os << "  unit_box=" << this->unit_box() << ",\n";
  os << "  subdivision_coordinates=" << this->subdivisions() << "\n";
  os << ")\n";
  return os;
}

template<class R>
std::ostream&
PartitionTreeCell<R>::write(std::ostream& os) const
{
  os << "PartitionTreeCell<" << name<R>() << ">(\n";
  os << "  bounds=" << this->subdivision_cell() << ",\n";
  os << "  rectangle=" << Box<R>(*this) << "\n";
  os << ")\n";
  return os;
}

template<class R>
std::ostream&
PartitionTree<R>::write(std::ostream& os) const
{
  os << "PartitionTree<" << name<R>() << ">(\n";
  os << "  unit_box=" << this->unit_box() << ",\n";
  os << "  subdivisions=" << this->subdivisions() << "\n";
  os << "  words="; write_sequence(os, this->binary_tree().begin(), this->binary_tree().end()); os << ",\n";
  os << "  blocks=["; write_sequence(os,  this->subdivision_tree().begin(), this->subdivision_tree().end()); os << ",\n";
  os << "  cells=["; write_sequence(os,  this->begin(), this->end()); os << ",\n";
  os << ")\n";
  return os;
}

template<class R>
std::ostream&
PartitionTreeSet<R>::write(std::ostream& os) const
{
  os << "PartitionTreeSet<" << name<R>() << ">(\n";
  os << "  unit_box=" << this->unit_box() << ",\n";
  os << "  subdivisions=" << this->subdivisions() << ",\n";
  os << "  tree=" << this->binary_tree() << ",\n";
  os << "  mask=" << this->mask() << ",\n";
  os << "  words="; write_sequence(os, this->subdivision_set().words().begin(), this->subdivision_set().words().end()); os << ",\n";
  os << "  blocks=["; write_sequence(os,  this->subdivision_set().begin(), this->subdivision_set().end()); os << ",\n";
  os << "  cells=["; write_sequence(os,  this->begin(), this->end()); os << ",\n";
  os << ")\n";
  return os;
}

template<class R>
std::ostream&
PartitionTreeSet<R>::summarize(std::ostream& os) const
{
  os << "PartitionTreeSet<" << name<R>() << ">(";
  os << " unit_box=" << this->unit_box() << ",";
  os << " size=" << this->size() << " )";
  return os;
}


} // namespace Ariadne
