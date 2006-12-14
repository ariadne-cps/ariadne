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
    PartitionTreeSet<R>::operator GridBlockListSet<R>() const 
    {
      throw NotImplemented(__PRETTY_FUNCTION__);
      //FiniteGrid<R>* grid_ptr=new FiniteGrid<R>(this->bounding_box(),this->subdivisions());
      //GridBlockListSet<R> res(*grid_ptr);
      //for(const_iterator iter=begin(); iter!=end(); ++iter) {
      //  res.push_back(GridBlock<R>(res.grid(),Rectangle<R>(*iter)));
      //}
      //return res;
    }



    template<class R>
    PartitionTreeSet<R>::operator ListSet<R,Rectangle>() const 
    {
      ListSet<R,Rectangle> res(this->dimension());
      for(const_iterator iter=begin(); iter!=end(); ++iter) {
        res.push_back(Rectangle<R>(*iter));
      }
      return res;
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
    operator<<(std::ostream& os, const PartitionScheme<R>& g)
    {
      os << "PartitionScheme<" << name<R>() << ">(\n";
      os << "  unit_box=" << g.unit_box() << ",\n";
      os << "  subdivision_coordinates=" << g.subdivisions() << "\n";
      os << ")\n";
      return os;
    }

    template<class R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTreeCell<R>& c)
    {
      os << "PartitionTreeCell<" << name<R>() << ">(\n";
      os << "  bounds=" << c.subdivision_cell() << ",\n";
      os << "  rectangle=" << Rectangle<R>(c) << "\n";
      os << ")\n";
      return os;
    }

    template<class R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTree<R>& pt)
    {
      os << "PartitionTree<" << name<R>() << ">(\n";
      os << "  unit_box=" << pt.unit_box() << ",\n";
      os << "  subdivisions=" << pt.subdivisions() << "\n";
      os << "  words="; Base::write_sequence(os, pt.binary_tree().begin(), pt.binary_tree().end()); os << ",\n";
      os << "  blocks=["; Base::write_sequence(os,  pt.subdivision_tree().begin(), pt.subdivision_tree().end()); os << ",\n";
      os << "  cells=["; Base::write_sequence(os,  pt.begin(), pt.end()); os << ",\n";
      os << ")\n";
      return os;
    }

    template<class R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTreeSet<R>& pts)
    {
      os << "PartitionTreeSet<" << name<R>() << ">(\n";
      os << "  unit_box=" << pts.unit_box() << ",\n";
      os << "  subdivisions=" << pts.subdivisions() << ",\n";
      os << "  tree=" << pts.binary_tree() << ",\n";
      os << "  mask=" << pts.mask() << ",\n";
      os << "  words="; Base::write_sequence(os, pts.subdivision_set().words().begin(), pts.subdivision_set().words().end()); os << ",\n";
      os << "  blocks=["; Base::write_sequence(os,  pts.subdivision_set().begin(), pts.subdivision_set().end()); os << ",\n";
      os << "  cells=["; Base::write_sequence(os,  pts.begin(), pts.end()); os << ",\n";
      os << ")\n";
      return os;
    }

  }
}
