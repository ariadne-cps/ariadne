/***************************************************************************
 *            subdivision_tree_set.cc
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

#include <iostream>
#include <string>
#include <sstream>
#include <exception>
#include <cassert>


#include "combinatoric/lattice_set.h"
#include "combinatoric/subdivision_tree_set.h"
#include "combinatoric/exceptions.h"

namespace Ariadne {
namespace Combinatoric {

const bool branch=BinaryTree::branch;
const bool leaf=BinaryTree::leaf;
const bool left=BinaryTree::left;
const bool right=BinaryTree::right;
    
 
sequence<dimension_type> 
SubdivisionSequence::_default(dimension_type n) 
{
  dimension_type coords[n];
  for(dimension_type i=0; i!=n; ++i) {
    coords[i]=i;
  }
  return sequence<dimension_type>(coords,coords,coords+n);
}

dimension_type 
SubdivisionSequence::_compute_dimension() 
{
  dimension_type result=0;
  const sequence<dimension_type>& ss(this->_sequence);
  for(size_type i=0; i!=ss.body_size()+ss.tail_size(); ++i) {
    if(ss[i]>=result) {
      result=ss[i]+1;
    }
  }
  return result;
}
    
void
MaskedBinaryTree::reduce()
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
  this->_tree=BinaryTree(new_tree);
  this->_mask=BooleanArray(new_mask);
}

    
SubdivisionCell::SubdivisionCell(const SubdivisionSequence& ss,
                                 const BinaryWord& w)
  : _subdivisions(w), _bounds(ss.dimension())
{
  _compute_bounds(ss,w);
}
    

SubdivisionCell::SubdivisionCell(const dimension_type& dim, 
                                 const depth_type& dpth, 
                                 const size_type& pos)
  : _subdivisions(), _bounds(dim)
{
  ARIADNE_PRECONDITION(pos<(1u<<dpth));
  _compute_word(dpth,pos);
  _compute_bounds(dim,this->_subdivisions);
}
    
SubdivisionCell::SubdivisionCell(const dimension_type& dim, 
                                 const word_type& subd)
  : _subdivisions(subd), _bounds(dim)
{
  _compute_bounds(dim,subd);
}
    


void
SubdivisionCell::_compute_word(const depth_type& dpth, 
                               const size_type& pos)
{
  while(!_subdivisions.size()==0) { _subdivisions.pop_back(); }
  for(depth_type i=0; i!=dpth;++i) {
    _subdivisions.push_back( (pos & (1<<i)) ? 1 : 0 );
  }
}

void
SubdivisionCell::_compute_bounds(const dimension_type& d,
                                 const word_type& w)
{
  _bounds.resize(d);
  for(dimension_type i=0; i!=d; ++i) {
    _bounds[i].lower=0.0;
    _bounds[i].upper=1.0;
  }
  for(size_type j=0; j!=w.size(); ++j) {
    dimension_type i=j%d;
    dyadic_type c=(_bounds[i].lower+_bounds[i].upper)/2;
    if(w[j]==left) {
      _bounds[i].upper=c;
    } else {
      _bounds[i].lower=c;
    }
  }
}

void
SubdivisionCell::_compute_bounds(const SubdivisionSequence& ss,
                                 const BinaryWord& w)
{
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    _bounds[i].lower=0.0;
    _bounds[i].upper=1.0;
  }
  for(size_type j=0; j!=w.size(); ++j) {
    dimension_type i=ss[j];
    dyadic_type c=(_bounds[i].lower+_bounds[i].upper)/2;
    if(w[j]==left) {
      _bounds[i].upper=c;
    }
    else {
      _bounds[i].lower=c;
    }
  }
}

void
SubdivisionCell::_split(bool lr) 
{
  uint i=this->depth()%this->dimension();
  dyadic_type c=(this->_bounds[i].lower+this->_bounds[i].upper)/2;
  if(lr) { this->_bounds[i].lower=c; } else { this->_bounds[i].upper=c; }
  this->_subdivisions.push_back(lr); 
}

dyadic_type
SubdivisionBox::volume() const
{
  dyadic_type result=1;
  for(dimension_type i=0; i!=this->dimension(); ++i) {
    result*=(this->upper_bound(i)-this->lower_bound(i));
  }
  return result;
}

SubdivisionBox 
hull(const SubdivisionBox& sbx1, const SubdivisionBox& sbx2)
{
  ARIADNE_PRECONDITION(sbx1.dimension()==sbx2.dimension());
  dimension_type dim=sbx1.dimension();
  SubdivisionBox result(dim);
  result._depth=std::max(sbx1.depth(),sbx2.depth());
  for(dimension_type i=0; i!=dim; ++i) {
    result._bounds[i].lower=std::min(sbx1.lower_bound(i),sbx2.lower_bound(i));
    result._bounds[i].upper=std::max(sbx1.upper_bound(i),sbx2.upper_bound(i));
  }
  return result;
}


void
SubdivisionCellListSet::unique_sort()
{
  std::sort(this->_words.begin(),this->_words.end());
  std::vector<word_type>::iterator new_end=std::unique(this->_words.begin(),this->_words.end());
  this->_words.erase(new_end,this->_words.end());
}

SubdivisionCellListSetIterator 
SubdivisionCellListSet::begin() const 
{
  return SubdivisionCellListSetIterator(this->_dimension,this->_words.begin());
}

SubdivisionCellListSetIterator 
SubdivisionCellListSet::end() const 
{
  return SubdivisionCellListSetIterator(this->_dimension,this->_words.end());
}






SubdivisionTree::SubdivisionTree(const SubdivisionSequence& ss,
                                 const BinaryTree& t)
  : _subdivisions(ss), _tree(t)
{
}




SubdivisionTreeSet::SubdivisionTreeSet(const SubdivisionSequence& ss)
  :  _subdivisions(ss), _words()
{
}

SubdivisionTreeSet::SubdivisionTreeSet(const SubdivisionSequence& ss, 
                                       const BinaryTree& bt,
                                       const BooleanArray& ba)
  : _subdivisions(ss), _words(bt,ba)
{
  if(bt.size()!=ba.size()) {
    ARIADNE_THROW(LengthError,
                  "SubdivisionTreeSet::SubdivisionTreeSet(SubdivisionSequence ss, BinaryTree bt, BooleanArray ba)",
                  " bt.size()="<<bt.size()<<", ba.size()="<<ba.size());
  }
  this->reduce();
}

SubdivisionTreeSet::SubdivisionTreeSet(const SubdivisionTreeSet& sts)
  : _subdivisions(sts._subdivisions), _words(sts._words)
{
}

SubdivisionTreeSet::SubdivisionTreeSet(const LatticeMaskSet& ms) 
  : _subdivisions(ms.dimension()), _words()
{      
  dimension_type n=ms.dimension();
  
  LatticeBlock block=ms.block();
  SizeArray grid_sizes=block.sizes();
  SizeArray new_sizes(this->dimension());
  SizeArray depths(this->dimension());
  size_type depth=0;
  
  /* Compute grid sizes as powers of two */
  for(dimension_type i=0; i!=n; ++i) {
    depths[i]=Numeric::log2_ceil(grid_sizes[i]);
    new_sizes[i]=Numeric::exp2(depths[i]);
    depth+=depths[i];
  }
  
  if(grid_sizes!=new_sizes) {
    throw std::runtime_error("Can only convert LatticeMaskSet to SubdivisionTreeSet "
                             "if all dimensions are subdivided as powers of two.");
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
    SubdivisionCell c(subdivisions,word);
    LatticeBlock r=compute_block(c,block);
    if(subset(r,ms)) {
      tree.push_back(leaf);
      mask.push_back(true);
    }
    else if(!overlap(r,ms)) {
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
  
  this->_subdivisions=subdivisions;
  this->_words=MaskedBinaryTree(BinaryTree(tree),mask);
  
  this->reduce();  
}

SizeArray
SubdivisionTreeSet::depths() const
{
  SizeArray result(this->dimension(),0);
  for(size_type j=0; j!=this->depth(); ++j) {
    result[_subdivisions[j]]+=1;
  }
  return result;
}




SubdivisionMaskSet::SubdivisionMaskSet(const dimension_type& dim, const depth_type& dpth) 
  : _dimension(dim), _depth(dpth), _mask(1<<dpth)
{
}

size_type
SubdivisionMaskSet::_position(const SubdivisionCell& c) const
{
  assert(c.depth()==this->_depth);
  size_type pos=0;
  for(size_type i=0; i!=c.depth(); ++i) {
    pos+= c.subdivisions()[i] * (1<<i);
  }
  return pos;
}

void
SubdivisionMaskSet::_adjoin_block(const SubdivisionCell& c) 
{
  ARIADNE_ASSERT(false);
}

void
SubdivisionMaskSet::_remove_block(const SubdivisionCell& c) 
{
  ARIADNE_ASSERT(false);
}

void
SubdivisionMaskSet::_restrict_block(const SubdivisionCell& c) 
{
  ARIADNE_ASSERT(false);
}

SubdivisionMaskSetIterator 
SubdivisionMaskSet::begin() const 
{
  return SubdivisionMaskSetIterator(this->_dimension,this->_mask.begin(),this->_mask.begin(),this->_mask.end());
}

SubdivisionMaskSetIterator 
SubdivisionMaskSet::end() const 
{
  return SubdivisionMaskSetIterator(this->_dimension,this->_mask.end(),this->_mask.begin(),this->_mask.end());
}






IndexArray 
compute_position(const SubdivisionSequence& ss, const BinaryWord& bw, const LatticeBlock& r)
{
  IndexArray lower=r.lower_corner();
  IndexArray upper=r.upper_corner();
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
  return lower;
}


LatticeBlock 
compute_block(const SubdivisionCell& c, 
              const LatticeBlock& r)
{
  IndexArray lower(r.dimension());
  IndexArray upper(r.dimension());
  SizeArray sizes=r.sizes();
  for(dimension_type i=0; i!=r.dimension(); ++i) {
    lower[i]=r.lower_bound(i)+index_type(c.lower_bound(i)*sizes[i]);
    upper[i]=r.lower_bound(i)+index_type(c.upper_bound(i)*sizes[i]);
  }
  return LatticeBlock(lower,upper);
}

std::ostream&
operator<<(std::ostream& os, const SubdivisionSequence& ss)
{
  return os << ss._sequence;
}

std::ostream&
operator<<(std::ostream& os, const SubdivisionCell& sc) 
{
  os << "SubdivisionCell" 
     << "( subdivisions=" << sc.subdivisions()
     << ", box=" << sc.box()
     << " )";
  return os;
}

std::ostream&
operator<<(std::ostream& os, const SubdivisionBox& sbx) 
{
  for(dimension_type i=0; i!=sbx.dimension(); ++i) {
    if(i!=0) { os << 'x'; }
    os << '[' << sbx.lower_bound(i) << ':' << sbx.upper_bound(i) << ']';
  }
  return os;
}  

std::ostream&
operator<<(std::ostream& os, const SubdivisionCellListSet& scls) 
{
  os << "SubdivisionCellListSet" 
     << "(\n  dimension=" << scls.dimension()
     << ",\n  size=" << scls.size()
     << ",\n  words=" << scls._words
     << "\n)\n";
  return os;
}


std::ostream&
operator<<(std::ostream& os, const SubdivisionMaskSet& sms) 
{
  os << "\nSubdivisionMaskSet" 
     << "(\n  dimension=" << sms.dimension()
     << ",\n  depth=" << sms.depth()
     << ",\n  size=" << sms.size()
     << ",\n  capacity=" << sms.capacity()
     << ",\n  mask=" << sms.mask()
     << "\n)\n";
  return os;
} 



} // namespace Combinatoric
} // namespace Ariadne
