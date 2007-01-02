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
#include "../geometry/grid_set.h"

#include <vector>

namespace Ariadne {
  namespace Geometry {

    
    
    template<class R>
    PartitionScheme<R>::PartitionScheme(const Rectangle<R>& bb) 
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
    PartitionTreeSet<R>::operator ListSet<R,Rectangle>() const 
    {
      ListSet<R,Rectangle> res(this->dimension());
      for(const_iterator iter=begin(); iter!=end(); ++iter) {
        res.push_back(Rectangle<R>(*iter));
      }
      return res;
    }
    
    
    
    template<class R>
    void
    PartitionTreeSet<R>::_instantiate_geometry_operators()
    {
      uint d=0;
      Rectangle<R>* r=0;
      Parallelotope<R>* pl=0;
      GridMaskSet<R>* gms=0;
      PartitionScheme<R>* ps=0;
      PartitionTreeSet<R>* pts=0;
      
      *pts=outer_approximation(*r,*ps,d);
      *pts=inner_approximation(*r,*ps,d);
      *pts=over_approximation(*r,*ps,d);
      *pts=under_approximation(*r,*ps,d);
    
      *pts=outer_approximation(*pl,*ps,d);
      *pts=inner_approximation(*pl,*ps,d);
      *pts=over_approximation(*pl,*ps,d);
      *pts=under_approximation(*pl,*ps,d);
    
      *pts=outer_approximation(*gms,*ps,d);
      *pts=inner_approximation(*gms,*ps,d);
      *pts=over_approximation(*gms,*ps,d);
      *pts=under_approximation(*gms,*ps,d);
    }

    
    template<class R, class S>
    PartitionTreeSet<R>
    outer_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
    {
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
    PartitionTreeSet<R>
    inner_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
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
    PartitionTreeSet<R>
    over_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
    {
      const Rectangle<R>& unit_box=ps.unit_box();
      const Combinatoric::SubdivisionSequence& subdivisions=ps.subdivisions();
      std::vector<bool> tree;
      std::vector<bool> mask;
      
      Combinatoric::BinaryWord word;
      
      do {
        Rectangle<R> cell=Rectangle<R>(PartitionTreeCell<R>(unit_box,subdivisions,word));
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
      
      return PartitionTreeSet<R>(unit_box,subdivisions,Combinatoric::BinaryTree(tree),BooleanArray(mask));
    }
    

    template<class R, class S>
    PartitionTreeSet<R>
    under_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
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
      os << "  rectangle=" << Rectangle<R>(*this) << "\n";
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
      os << "  words="; Base::write_sequence(os, this->binary_tree().begin(), this->binary_tree().end()); os << ",\n";
      os << "  blocks=["; Base::write_sequence(os,  this->subdivision_tree().begin(), this->subdivision_tree().end()); os << ",\n";
      os << "  cells=["; Base::write_sequence(os,  this->begin(), this->end()); os << ",\n";
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
      os << "  words="; Base::write_sequence(os, this->subdivision_set().words().begin(), this->subdivision_set().words().end()); os << ",\n";
      os << "  blocks=["; Base::write_sequence(os,  this->subdivision_set().begin(), this->subdivision_set().end()); os << ",\n";
      os << "  cells=["; Base::write_sequence(os,  this->begin(), this->end()); os << ",\n";
      os << ")\n";
      return os;
    }

  }
}
