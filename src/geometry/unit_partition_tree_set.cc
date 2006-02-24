/***************************************************************************
 *            unit_partition_tree_set.cc
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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

#include "base/arithmetic.h"

#include "geometry/rectangle.h"
#include "geometry/unit_grid_set.h"

#include "geometry/unit_partition_tree_set.h"
#include "geometry/partition_tree_operations.h"

namespace Ariadne {
  namespace Geometry {
    
    UnitPartitionTreeCell::UnitPartitionTreeCell(const dimension_type& d,
                                                 const SubdivisionSequence& ss,
                                                 const BinaryWord& w)
      : _subdivisions(ss), _word(w), _bounds(2*d)
    {
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        _bounds[2*i]=0;
        _bounds[2*i+1]=1;
      }
      for(size_type j=0; j!=_word.size(); ++j) {
        dimension_type i=_subdivisions[j];
        dyadic_type c=(_bounds[2*i]+_bounds[2*i+1])/2;
        if(w[j]==left) {
          _bounds[2*i+1]=c;
        }
        else {
          _bounds[2*i]=c;
        }
      }
    }

    Rectangle<UnitPartitionTreeCell::dyadic_type>
    UnitPartitionTreeCell::bounds() const
    {
      Rectangle<dyadic_type> result(this->dimension());
      for(dimension_type i=0; i!=this->dimension(); ++i) {
        result.set_lower_bound(i,this->lower_bound(i));
        result.set_upper_bound(i,this->upper_bound(i));
      }
      return result;
    }

    UnitPartitionTreeSet::UnitPartitionTreeSet(const SubdivisionSequence& ss)
      : _dimension(compute_dimension(ss)), _subdivisions(ss), 
        _tree(), _mask(1,false)
    {
    }
    
    UnitPartitionTreeSet::UnitPartitionTreeSet(const SubdivisionSequence& ss, 
                                               const BinaryTree& bt,
                                               const BooleanArray& ba)
      : _dimension(compute_dimension(ss)), _subdivisions(ss), 
        _tree(bt), _mask(ba)
    {
      assert(bt.size()==ba.size());
      this->reduce();
    }
    
    UnitPartitionTreeSet::UnitPartitionTreeSet(const UnitGridMaskSet& ms) 
      : _dimension(ms.dimension()),
        _subdivisions(default_subdivision_coordinates(ms.dimension())),
        _tree(),
        _mask()
    {
      dimension_type n=ms.dimension();
      dimension_type dimension=ms.dimension();
      UnitGridRectangle bounds=ms.bounds();
      SizeArray grid_sizes=bounds.sizes();
      SizeArray new_sizes(this->dimension());
      SizeArray depths(this->dimension());
      size_type depth=0;

      /* Compute grid sizes as powers of two */
      for(dimension_type i=0; i!=n; ++i) {
        depths[i]=log_ceil(2,grid_sizes[i]);
        new_sizes[i]=pow(2,depths[i]);
        depth+=depths[i];
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
      for(dimension_type j=0; j!=n; ++j) {
        sc.push_back(j);
      }
      while(sc.size()>n && (sc[sc.size()-1]==sc[sc.size()-n-1])) {
        sc.pop_back();
      }
      SubdivisionSequence subdivisions(sc.begin(),sc.end()-n,sc.end());
      
      std::vector<bool> tree;
      std::vector<bool> mask;
      BinaryWord word;
      
      do {
        //TODO: Construct full tree, then reduce
        UnitPartitionTreeCell c(dimension,subdivisions,word);
        UnitGridRectangle r=compute_block(c,bounds);
        if(subset(r,ms)) {
          //std::cerr <<  c.bounds() << " " << r << " subset" << std::endl;
          tree.push_back(leaf);
          mask.push_back(true);
        }
        else if(!interiors_intersect(r,ms)) {
          //std::cerr <<  c.bounds() << " " << r << " disjoint" << std::endl;
          tree.push_back(leaf);
          mask.push_back(false);
        }
        else {
          //std::cerr <<  c.bounds() << " " << r << " neither" << std::endl;
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
      
      this->_dimension=dimension;
      this->_subdivisions=subdivisions;
      this->_tree=BinaryTree(tree);
      this->_mask=mask;
    
      this->reduce();
    }

    SizeArray
    UnitPartitionTreeSet::depths() const
    {
      SizeArray result(this->dimension(),0);
      for(size_type j=0; j!=this->depth(); ++j) {
        result[_subdivisions[j]]+=1;
      }
      return result;
    }
    
    void
    UnitPartitionTreeSet::reduce()
    {
      const BooleanArray& tree(this->_tree.array());
      const BooleanArray& mask(this->_mask);

      std::vector<bool> new_tree;
      std::vector<bool> new_mask;
      
      if(tree.size()==1) {
        return; 
      }
      
      BooleanArray::const_iterator tree_iterator=tree.begin();
      BooleanArray::const_iterator mask_iterator=mask.begin();
      
      while(tree_iterator!=tree.end()) {
        new_tree.push_back(*tree_iterator);
        if(*tree_iterator==leaf) {
          new_mask.push_back(*mask_iterator);
          ++mask_iterator;
      
          size_type n=new_tree.size();
          size_type m=new_mask.size();
          while(new_tree.size()>=3 
                  && new_tree[n-1]==leaf && new_tree[n-2]==leaf && new_tree[n-3]==branch 
                  && new_mask[m-1]==new_mask[m-2]) 
          {
            new_tree.pop_back();
            new_tree.pop_back();
            new_tree.pop_back();
            new_tree.push_back(leaf);
            new_mask.pop_back();
            n=new_tree.size();
            m=new_mask.size();
          }
        }
        ++tree_iterator;
      }
      //std::cerr << tree << " " << mask << std::endl;
      //std::cerr << new_tree << " " << new_mask << std::endl;
      this->_tree=BinaryTree(new_tree);
      this->_mask=BooleanArray(new_mask);
    }

    std::ostream&
    operator<<(std::ostream& os, const UnitPartitionTreeCell& uptc) 
    {
      os << "UnitPartitionTreeCell(" 
         << "subdivisions=" << uptc.subdivisions() << ", "
         << "word=" << uptc.word() << ", "
         << "block=" << "[" << uptc.lower_bound(0) << "," << uptc.upper_bound(0) << "]";
      for(dimension_type i=1; i!=uptc.dimension(); ++i) {
        os << "x" << "[" << uptc.lower_bound(i) << "," << uptc.upper_bound(i) << "]";
      }
      os << ")";
      return os;
    }
    
    
    index_type compute_index(const SubdivisionSequence& ss, const BinaryWord& bw, const UnitGridRectangle& r)
    {
      IndexArray lower=r.lower();
      IndexArray upper=r.upper();
      for(size_type j=0; j!=bw.size(); ++j) {
        dimension_type i=ss[j];
        assert((lower[i]+upper[i])%2==0);
        index_type c=(lower[i]+upper[i])/2;
        if(bw[j]==left) {
          upper[i]=c;
        }
        else {
          lower[i]=c;
        }
      }
      for(dimension_type i=0; i!=r.dimension(); ++i) {
        assert(lower[i]+1 == upper[i]);
      }
      return compute_index(lower,r.lower(),r.strides());
    }
    
      
    UnitGridRectangle 
    compute_block(const UnitPartitionTreeCell& c, 
                  const UnitGridRectangle& r)
    {
      IndexArray lower(r.dimension());
      IndexArray upper(r.dimension());
      SizeArray sizes=r.sizes();
      for(dimension_type i=0; i!=r.dimension(); ++i) {
        lower[i]=r.lower_bound(i)+index_type(c.lower_bound(i)*sizes[i]);
        upper[i]=r.lower_bound(i)+index_type(c.upper_bound(i)*sizes[i]);
      }
      return UnitGridRectangle(lower,upper);
    }


  }
}
