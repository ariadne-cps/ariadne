/***************************************************************************
 *            partition_tree_set.code.h
 *
 *  1 July 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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


#include "partition_tree_set.h"

#include "../base/stlio.h"

#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_mask_set.h"
#include "../output/logging.h"

#include <vector>

namespace Ariadne {


    
    
template<class R>
Geometry::PartitionScheme<R>::PartitionScheme(const Rectangle<R>& bb) 
  : _unit_box(bb), 
    _subdivisions(bb.dimension()) 
{ }


// FIXME: Should use add_approx etc here.
// We use approximate values since the cell boundaries are defined
// by the lower bound an upper bound formulae (which must be the same given 
// the cell bound.
template<class R>
R
Geometry::PartitionTreeCell<R>::lower_bound(dimension_type i) const 
{
  return add_approx(_unit_box.lower_bound(i),
                    mul_approx(_subdivision_cell.lower_bound(i),
                               sub_approx(_unit_box.upper_bound(i),_unit_box.lower_bound(i))));
}

template<class R>
R
Geometry::PartitionTreeCell<R>::upper_bound(dimension_type i) const 
{
  return add_approx(_unit_box.lower_bound(i),
                    mul_approx(_subdivision_cell.upper_bound(i),
                               sub_approx(_unit_box.upper_bound(i),_unit_box.lower_bound(i))));
}



template<class R>
Geometry::PartitionTreeSet<R>::PartitionTreeSet(const GridMaskSet<R>& gms) 
  : _unit_box(gms.bounding_box()),
    _subdivision_set(gms._lattice_set)
{ }



template<class R>
Geometry::PartitionTreeSet<R>*
Geometry::PartitionTreeSet<R>::clone() const
{
  return new PartitionTreeSet<R>(*this);
}


template<class R>
tribool
Geometry::PartitionTreeSet<R>::contains(const Point<R>& pt) const
{
  return Geometry::contains(*this,pt);
}


template<class R>
tribool
Geometry::PartitionTreeSet<R>::intersects(const Rectangle<R>& r) const
{
  return !Geometry::disjoint(*this,r);
}


template<class R>
tribool
Geometry::PartitionTreeSet<R>::disjoint(const Rectangle<R>& r) const
{
  return Geometry::disjoint(*this,r);
}


template<class R>
tribool
Geometry::PartitionTreeSet<R>::superset(const Rectangle<R>& r) const
{
  return Geometry::subset(r,*this);
}


template<class R>
tribool
Geometry::PartitionTreeSet<R>::subset(const Rectangle<R>& r) const
{
  return Geometry::subset(*this,r);
}







template<class R>
Geometry::PartitionTreeSet<R>::operator ListSet< Rectangle<R> >() const 
{
  ListSet< Rectangle<R> > res(this->dimension());
  for(const_iterator iter=begin(); iter!=end(); ++iter) {
    res.push_back(Rectangle<R>(*iter));
  }
  return res;
}



template<class R>
void
Geometry::PartitionTreeSet<R>::_instantiate_geometry_operators()
{
  uint d=0;
  Rectangle<R>* r=0;
  Zonotope<R>* z=0;
  SetInterface<R>* s=0;
  GridMaskSet<R>* gms=0;
  PartitionScheme<R>* ps=0;
  PartitionTreeCell<R>* ptc=0;
  PartitionTreeSet<R>* pts=0;
  
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
Geometry::contains(const PartitionTreeSet<R>& pts, const Point<R>& pt)
{
  return subset(Rectangle<R>(pt),pts);
}


template<class R>
tribool
Geometry::disjoint(const PartitionTreeSet<R>& pts, const Rectangle<R>& r)
{
  tribool result=true;
  Rectangle<R> cell(pts.dimension());
  for(typename PartitionTreeSet<R>::const_iterator pts_iter=pts.begin();
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
Geometry::subset(const PartitionTreeSet<R>& pts, const Rectangle<R>& r)
{
  tribool result=true;
  Rectangle<R> cell(pts.dimension());
  for(typename PartitionTreeSet<R>::const_iterator pts_iter=pts.begin();
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
Geometry::subset(const Rectangle<R>& r, const PartitionTreeSet<R>& pts)
{
  return subset(r,ListSet< Rectangle<R> >(pts));
}




template<class R>
Geometry::PartitionTreeCell<R>
Geometry::over_approximation(const Rectangle<R>& r, const PartitionScheme<R>& ps)
{
  if(possibly(r.empty())) {
    ARIADNE_THROW(EmptyInterior,"PartitionTreeCell<R> over_approximation(const Rectangle<R>& r, const PartitionScheme<R>& ps)"," with r="<<r);
  }
  
  if(!subset(r,ps.unit_box())) {
    ARIADNE_THROW(std::runtime_error,"ParttitionTreeCell over_approximation(Rectangle r, PartitionScheme ps, uint d)"," with r="<<r<<", ps="<<ps<<": r is not a subset of ps.unit_box()");
  }
  
  const Rectangle<R>& unit_box=ps.unit_box();
  const Combinatoric::SubdivisionSequence& subdivisions=ps.subdivisions();
  
  Combinatoric::BinaryWord word;
  Combinatoric::BinaryWord new_word;
  
  while(true) {
    new_word=word;
    new_word.push_back(Combinatoric::BinaryTree::leaf);
    if(!subset(r,Rectangle<R>(PartitionTreeCell<R>(unit_box,subdivisions,new_word)))) {
      new_word.pop_back();
      word.push_back(Combinatoric::BinaryTree::right);
      if(!subset(r,Rectangle<R>(PartitionTreeCell<R>(unit_box,subdivisions,word)))) {
        break;
      }
    }
  }
  return PartitionTreeCell<R>(unit_box,subdivisions,word);
}


template<class R, class S>
Geometry::PartitionTreeSet<R>
Geometry::outer_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
{
  ARIADNE_LOG(4,"outer_approximation(S set, PartitionScheme ps, uint depth");
  const Rectangle<R>& bounding_box=ps.unit_box();
  const Combinatoric::SubdivisionSequence& subdivisions=ps.subdivisions();
  std::vector<bool> tree;
  std::vector<bool> mask;
  
  Combinatoric::BinaryWord word;
  
  do {
    Rectangle<R> cell=Rectangle<R>(PartitionTreeCell<R>(bounding_box,subdivisions,word));
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
  
  return PartitionTreeSet<R>(bounding_box,subdivisions,Combinatoric::BinaryTree(tree),BooleanArray(mask));
}

template<class R, class S>
Geometry::PartitionTreeSet<R>
Geometry::inner_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
{
  const Rectangle<R>& unit_box=ps.unit_box();
  const Combinatoric::SubdivisionSequence& subdivisions=ps.subdivisions();
  std::vector<bool> tree;
  std::vector<bool> mask;
  
  Combinatoric::BinaryWord word;
  
  do {
    Rectangle<R> cell=Rectangle<R>(PartitionTreeCell<R>(unit_box,subdivisions,word));
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
  
  return PartitionTreeSet<R>(unit_box,subdivisions,Combinatoric::BinaryTree(tree),BooleanArray(mask));
}


template<class R, class S>
Geometry::PartitionTreeSet<R>
Geometry::under_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
{
  const Rectangle<R>& unit_box=ps.unit_box();
  const Combinatoric::SubdivisionSequence& subdivisions=ps.subdivisions();
  std::vector<bool> tree;
  std::vector<bool> mask;
  
  Combinatoric::BinaryWord word;
  
  do {
    Rectangle<R> cell=Rectangle<R>(PartitionTreeCell<R>(unit_box,subdivisions,word));
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
  
  return PartitionTreeSet<R>(unit_box,subdivisions,Combinatoric::BinaryTree(tree),BooleanArray(mask));
}



template<class R>
std::ostream&
Geometry::PartitionScheme<R>::write(std::ostream& os) const
{
  os << "PartitionScheme<" << Numeric::name<R>() << ">(\n";
  os << "  unit_box=" << this->unit_box() << ",\n";
  os << "  subdivision_coordinates=" << this->subdivisions() << "\n";
  os << ")\n";
  return os;
}

template<class R>
std::ostream&
Geometry::PartitionTreeCell<R>::write(std::ostream& os) const
{
  os << "PartitionTreeCell<" << Numeric::name<R>() << ">(\n";
  os << "  bounds=" << this->subdivision_cell() << ",\n";
  os << "  rectangle=" << Rectangle<R>(*this) << "\n";
  os << ")\n";
  return os;
}

template<class R>
std::ostream&
Geometry::PartitionTree<R>::write(std::ostream& os) const
{
  os << "PartitionTree<" << Numeric::name<R>() << ">(\n";
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
Geometry::PartitionTreeSet<R>::write(std::ostream& os) const
{
  os << "PartitionTreeSet<" << Numeric::name<R>() << ">(\n";
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
Geometry::PartitionTreeSet<R>::summarize(std::ostream& os) const
{
  os << "PartitionTreeSet<" << Numeric::name<R>() << ">(";
  os << " unit_box=" << this->unit_box() << ",";
  os << " size=" << this->size() << " )";
  return os;
}


} // namespace Ariadne
