/***************************************************************************
 *            partition_tree_set.tpl
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

#include "../utility/stlio.h"

#include "../geometry/rectangle.h"
#include "../geometry/list_set.h"
#include "../geometry/grid_set.h"

#include <vector>

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    PartitionScheme<R>::PartitionScheme(const Rectangle<R>& bb) 
      : _bounding_box(bb), 
        _subdivisions(bb.dimension()) 
      { }

    template<typename R>
    PartitionTreeCell<R>::operator Rectangle<R>() const {
      Rectangle<R> result(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result[i]=Interval<R>(
            _bounding_box.lower_bound(i)+_unit_cell.lower_bound(i)
                *(_bounding_box.upper_bound(i)-_bounding_box.lower_bound(i)),
            _bounding_box.lower_bound(i)+_unit_cell.upper_bound(i)
                *(_bounding_box.upper_bound(i)-_bounding_box.lower_bound(i)) 
          ); 
      }
      return result;
    }

    template<typename R>
    PartitionTreeSet<R>::PartitionTreeSet(const GridMaskSet<R>& gms) 
      : _bounding_box(gms.bounding_box()),
        _unit_set(gms._lattice_set)
    { }

    template<typename R>
    PartitionTreeSet<R>::operator ListSet<R,Rectangle>() const 
    {
      ListSet<R,Rectangle> res(this->dimension());
      for(const_iterator iter=begin(); iter!=end(); ++iter) {
        res.push_back(Rectangle<R>(*iter));
      }
      return res;
    }

    void advance(BinaryWord& word) {
      while(!word.empty() && word.back()==BinaryTree::right) {
        word.pop_back();
      }
      if(!word.empty()) {
        word.set_back(BinaryTree::right);
      }
    }
    
    template<typename R, class S>
    PartitionTreeSet<R>
    outer_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
    {
      const Rectangle<R>& bounding_box=ps.bounding_box();
      const SubdivisionSequence& subdivisions=ps.subdivisions();
      std::vector<bool> tree;
      std::vector<bool> mask;
      
      BinaryWord word;
      
      do {
        Rectangle<R> cell=Rectangle<R>(PartitionTreeCell<R>(bounding_box,subdivisions,word));
        if(word.size()==depth+1 || subset(cell,s)) {
          tree.push_back(BinaryTree::leaf);
          mask.push_back(true);
          advance(word);
        }  
        else if(disjoint(cell,s)) {
          tree.push_back(BinaryTree::leaf);
          mask.push_back(false);
          advance(word);
        }
        else {
          tree.push_back(BinaryTree::branch);
          word.push_back(BinaryTree::left);
        }
      } while(!word.empty());
      
      return PartitionTreeSet<R>(bounding_box,subdivisions,BinaryTree(tree),BooleanArray(mask));
    }
    
    template<typename R, class S>
    PartitionTreeSet<R>
    inner_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
    {
      const Rectangle<R>& bounding_box=ps.bounding_box();
      const SubdivisionSequence& subdivisions=ps.subdivisions();
      std::vector<bool> tree;
      std::vector<bool> mask;
      
      BinaryWord word;
      
      do {
        Rectangle<R> cell=Rectangle<R>(PartitionTreeCell<R>(bounding_box,subdivisions,word));
        if(word.size()==depth+1 || !interiors_intersect(cell,s)) {
          tree.push_back(BinaryTree::leaf);
          mask.push_back(false);
          advance(word);
        }  
        else if(inner_subset(cell,s)) {
          tree.push_back(BinaryTree::leaf);
          mask.push_back(true);
          advance(word);
        }
        else {
          tree.push_back(BinaryTree::branch);
          word.push_back(BinaryTree::left);
         }
      } while(!word.empty());
      
      return PartitionTreeSet<R>(bounding_box,subdivisions,BinaryTree(tree),BooleanArray(mask));
    }
    

    template<typename R, class S>
    PartitionTreeSet<R>
    over_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
    {
      const Rectangle<R>& bounding_box=ps.bounding_box();
      const SubdivisionSequence& subdivisions=ps.subdivisions();
      std::vector<bool> tree;
      std::vector<bool> mask;
      
      BinaryWord word;
      
      do {
        Rectangle<R> cell=Rectangle<R>(PartitionTreeCell<R>(bounding_box,subdivisions,word));
        if(word.size()==depth+1 || subset(cell,s)) {
          tree.push_back(BinaryTree::leaf);
          mask.push_back(true);
          advance(word);
        }  
        else if(!interiors_intersect(cell,s)) {
          tree.push_back(BinaryTree::leaf);
          mask.push_back(false);
          advance(word);
        }
        else {
          tree.push_back(BinaryTree::branch);
          word.push_back(BinaryTree::left);
        }
      } while(!word.empty());
      
      return PartitionTreeSet<R>(bounding_box,subdivisions,BinaryTree(tree),BooleanArray(mask));
    }
    

    template<typename R, class S>
    PartitionTreeSet<R>
    under_approximation(const S& s, const PartitionScheme<R>& ps, const uint depth)
    {
      const Rectangle<R>& bounding_box=ps.bounding_box();
      const SubdivisionSequence& subdivisions=ps.subdivisions();
      std::vector<bool> tree;
      std::vector<bool> mask;
      
      BinaryWord word;
      
      do {
        Rectangle<R> cell=Rectangle<R>(PartitionTreeCell<R>(bounding_box,subdivisions,word));
        if(word.size()==depth+1 || !interiors_intersect(cell,s)) {
          tree.push_back(BinaryTree::leaf);
          mask.push_back(false);
          advance(word);
        }  
        else if(subset(cell,s)) {
          tree.push_back(BinaryTree::leaf);
          mask.push_back(true);
          advance(word);
        }
        else {
          tree.push_back(BinaryTree::branch);
          word.push_back(BinaryTree::left);
         }
      } while(!word.empty());
      
      return PartitionTreeSet<R>(bounding_box,subdivisions,BinaryTree(tree),BooleanArray(mask));
    }
    

    
/*
    template<typename R>
    PartitionTreeSet<R>::operator GridRectangleListSet<R>() const 
    {
      std::cerr << "PartitionTreeSet<R>::operator GridRectangleListSet<R>() const" << std::endl;
      FiniteGrid<R>* grid_ptr=new FiniteGrid<R>(this->bounding_box(),this->subdivisions());
      GridRectangleListSet<R> res(*grid_ptr);
      for(const_iterator iter=begin(); iter!=end(); ++iter) {
        res.push_back(GridRectangle<R>(res.grid(),Rectangle<R>(*iter)));
      }
      return res;
    }
*/

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionScheme<R>& g)
    {
      os << "PartitionScheme<" << name<R>() << ">(\n";
      os << "  bounding_box=" << g.bounding_box() << ",\n";
      os << "  subdivision_coordinates=" << g.subdivisions() << "\n";
      os << ")\n";
      return os;
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTreeCell<R>& c)
    {
      os << "PartitionTreeCell<" << name<R>() << ">(\n";
      os << "  bounds=" << c.unit_cell() << ",\n";
      os << "  rectangle=" << Rectangle<R>(c) << "\n";
      os << ")\n";
      return os;
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTree<R>& pt)
    {
      os << "PartitionTree<" << name<R>() << ">(\n";
      os << "  bounding_box=" << pt.bounding_box() << ",\n";
      os << "  subdivisions=" << pt.subdivisions() << "\n";
      os << "  words="; Utility::write_sequence(os, pt.binary_tree().begin(), pt.binary_tree().end()); os << ",\n";
      os << "  blocks=["; Utility::write_sequence(os,  pt.unit_tree().begin(), pt.unit_tree().end()); os << ",\n";
      os << "  cells=["; Utility::write_sequence(os,  pt.begin(), pt.end()); os << ",\n";
      os << ")\n";
      return os;
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTreeSet<R>& pts)
    {
      os << "PartitionTreeSet<" << name<R>() << ">(\n";
      os << "  bounding_box=" << pts.bounding_box() << ",\n";
      os << "  subdivisions=" << pts.subdivisions() << ",\n";
      os << "  tree=" << pts.binary_tree() << ",\n";
      os << "  mask=" << pts.mask() << ",\n";
      os << "  words="; Utility::write_sequence(os, pts.unit_set().words().begin(), pts.unit_set().words().end()); os << ",\n";
      os << "  blocks=["; Utility::write_sequence(os,  pts.unit_set().begin(), pts.unit_set().end()); os << ",\n";
      os << "  cells=["; Utility::write_sequence(os,  pts.begin(), pts.end()); os << ",\n";
      os << ")\n";
      return os;
    }

  }
}
