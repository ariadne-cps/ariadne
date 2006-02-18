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


#include "geometry/partition_tree_set.h"

#include "geometry/rectangle.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"

#include <vector>

namespace Ariadne {
  namespace Geometry {

    template<typename R>
    PartitionTreeCell<R>::operator Rectangle<R>() const {
      Rectangle<R> res(_scheme.bounding_box());
      SubdivisionSequence::const_iterator coord_iter=_scheme.subdivision_coordinates().begin();
      BinaryWord::const_iterator word_iter=_word.begin();

      while(word_iter!=_word.end()) {
        size_type i=(*coord_iter);
        R centre = ( res.lower_bound(i) + res.upper_bound(i) ) / 2;
        if( (*word_iter)==left ) {
          res.set_upper_bound(i,centre); }
        else {
          res.set_lower_bound(i,centre);
        }
        ++word_iter;
        ++coord_iter;
      }

      return res;
    }

    template<typename R>
    PartitionTreeSet<R>::PartitionTreeSet(const GridMaskSet<R>& gms) 
      : _ptree(gms.bounding_box()), _mask()
    {
      dimension_type n=this->dimension();
      SizeArray grid_subdivisions=gms.bounds().sizes();
      SizeArray subdivisions(this->dimension());
      SizeArray depths(this->dimension());
      size_type depth=0;

      /* Compute grid sizes as powers of two */
      for(dimension_type i=0; i!=dimension(); ++i) {
        subdivisions[i]=1;
        depths[i]=1;
        size_type sd=grid_subdivisions[i];
        while(sd!=1) {
          if(sd%2==1) {
            sd+=1;
          }
          sd/=2;
          subdivisions[i]*=2;
          depths[i]+=1;
          depth+=1;
       }
      }
      
      /* Compute subdivision coordinates */
      std::vector<dimension_type> sc;
      for(size_type i=0; i!=depth; ++i) {
        dimension_type coordinate=0;
        for(dimension_type j=1; j!=n; ++j) {
          if(depths[j]>depths[coordinate]) {
             // || depths[j]=depth[coordinate] && gms.grid().upper_bound(j) > gms.grid().upper_bound(i)) {
            coordinate=j;
          }
        }
        sc.push_back(coordinate);
        depths[coordinate]-=1;
      }
      
      SubdivisionSequence subdivision_coordinates(sc.begin(),sc.end()-n,sc.end());
      Rectangle<R> bounding_box=GridRectangle<R>(gms.grid(),gms.bounds());
      std::vector<bool> tree;
      std::vector<bool> mask;
      BinaryWord word;
      
      do {
        //TODO: Write non-templated subset routine direct from grid
        PartitionTreeCell<R> c(bounding_box,subdivision_coordinates,word);
        Rectangle<R> r(c);
        GridRectangle<R> gr=over_approximation_of_intersection(r,gms.bounding_box(),gms.grid());
        if(subset(r,gms)) {
          tree.push_back(leaf);
          mask.push_back(true);
        }
        else if(!interiors_intersect(r,gms)) {
          tree.push_back(leaf);
          mask.push_back(false);
        }
        else {
          tree.push_back(branch);
        }          
        if(tree.back()==leaf) {
          while(!word.empty() && word.back()==right) {
            word.pop_back();
          }
          if(!word.empty()) {
            word.set_back(right);
          }
        }
        else {
          word.push_back(left);
        }
      } 
      while(!word.empty());
      
      this->_ptree=PartitionTree<R>(bounding_box,subdivision_coordinates,BinaryTree(tree));
      this->_mask=BooleanArray(mask.begin(),mask.end());
    }

    


    template<typename R>
    PartitionTreeSet<R>::operator ListSet<R,Rectangle>() const 
    {
      ListSet<R,Rectangle> res(this->dimension());
      for(const_iterator iter=begin(); iter!=end(); ++iter) {
        res.push_back(Rectangle<R>(*iter));
      }
      return res;
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
      os << "  subdivision_coordinates=" << g.subdivision_coordinates() << "\n";
      os << ")\n";
      return os;
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTree<R>& t)
    {
      os << "PartitionTree<" << name<R>() << ">(\n";
      os << "  bounding_box=" << t.bounding_box() << ",\n";
      os << "  subdivision_coordinates=" << t.subdivision_coordinates() << "\n";
      os << "  words=";
      write_sequence(os, t.tree().begin(), t.tree().end());
      os << "\n";
      os << "  rectangles=[ " << Rectangle<R>(*t.begin());
      for(typename PartitionTree<R>::const_iterator ptree_iter=++t.begin() ; ptree_iter!=t.end(); ++ptree_iter) {
        os << ", " << Rectangle<R>(*ptree_iter);
      }
      os << " ]\n";
      os << ")\n";
      return os;
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTreeCell<R>& c)
    {
      os << "PartitionTreeCell<" << name<R>() << ">(\n";
      os << "  word=" << c.index() << ",\n";
      os << "  rectangle=" << Rectangle<R>(c) << "\n";
      os << ")\n";
      return os;
    }

    template<typename R>
    std::ostream&
    operator<<(std::ostream& os, const PartitionTreeSet<R>& set)
    {
      os << "PartitionTreeSet<" << name<R>() << ">(\n";
      os << "  bounding_box=" << set.bounding_box() << ",\n";
      os << "  subdivision_coordinates=" << set.subdivision_coordinates() << ",\n";
      os << "  tree=" << set.tree() << ",\n";
      os << "  mask=" << set.mask() << ",\n";
      os << "  cells="; write_sequence(os, set.tree().begin(), set.tree().end()); os << ",\n";
      BinarySubtreeIterator bstb(set.tree().begin(),set.mask().begin());
      BinarySubtreeIterator bste(set.tree().end(),set.mask().end());
      os << "  words="; write_sequence(os, bstb, bste); os << ",\n";
      os << ")\n";
      return os;
    }

  }
}
