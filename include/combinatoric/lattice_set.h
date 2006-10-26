/***************************************************************************
 *            lattice_set.h
 *
 *  Copyright  2005,6  Alberto Casagrande, Pieter Collins
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

/*! \file lattice_set.h
 *  \brief Sets on an integer grid.
 */

#ifndef _ARIADNE_LATTICE_SET_H
#define _ARIADNE_LATTICE_SET_H

#include <vector>
#include <iosfwd>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include "../base/array.h"
#include "../base/iterator.h"

#include "../numeric/interval.h"
#include "../combinatoric/array_operations.h"

namespace Ariadne {
  namespace Combinatoric {
    class LatticeCell;

    class LatticePoint;
      
    class LatticeBlock;
    class LatticeBlockIterator;

    class LatticeCellListSet;
    class LatticeCellListSetIterator;

    class LatticeBlockListSet;
    class LatticeBlockListSetIterator;

    class LatticeMaskSet;
    class LatticeMaskSetIterator;

    bool disjoint(const LatticeBlock&, const LatticeBlock&);
    bool disjoint(const LatticeBlock&, const LatticeMaskSet&);
    bool disjoint(const LatticeMaskSet&, const LatticeMaskSet&);
    bool overlap(const LatticeBlock&, const LatticeBlock&);
    bool overlap(const LatticeBlock&, const LatticeMaskSet&);
    bool overlap(const LatticeMaskSet&, const LatticeMaskSet&);
    bool subset(const LatticeCell&, const LatticeBlock&);
    bool subset(const LatticeBlock&, const LatticeBlock&);
    bool subset(const LatticeCellListSet&, const LatticeBlock&);
    bool subset(const LatticeMaskSet&, const LatticeBlock&);
    bool subset(const LatticeCell&, const LatticeMaskSet&);
    bool subset(const LatticeBlock&, const LatticeMaskSet&);
    bool subset(const LatticeCellListSet&, const LatticeMaskSet&);
    bool subset(const LatticeMaskSet&, const LatticeMaskSet&);
    
    LatticeBlock regular_intersection(const LatticeBlock&, const LatticeBlock&);
    LatticeCellListSet regular_intersection(const LatticeCellListSet&, const LatticeMaskSet&);
    LatticeCellListSet regular_intersection(const LatticeMaskSet&, const LatticeCellListSet&);
    LatticeMaskSet regular_intersection(const LatticeBlock&, const LatticeMaskSet&);
    LatticeMaskSet regular_intersection(const LatticeMaskSet&, const LatticeBlock&);
    LatticeMaskSet regular_intersection(const LatticeMaskSet&, const LatticeMaskSet&);
    LatticeMaskSet join(const LatticeMaskSet&, const LatticeMaskSet&);
    LatticeCellListSet difference(const LatticeCellListSet&, const LatticeMaskSet&);
    LatticeMaskSet difference(const LatticeMaskSet&, const LatticeMaskSet&);

    bool operator<(const LatticeCell& lc1, const LatticeCell& lc2);
    
    std::istream& operator>>(std::istream& os, LatticeCell&);
    std::istream& operator>>(std::istream& os, LatticeBlock&);
    
    std::ostream& operator<<(std::ostream& os, const LatticeCell&);
    std::ostream& operator<<(std::ostream& os, const LatticeBlock&);
    std::ostream& operator<<(std::ostream& os, const LatticeMaskSet&);
    std::ostream& operator<<(std::ostream& os, const LatticeCellListSet&);
    std::ostream& operator<<(std::ostream& os, const LatticeBlockListSet&);


    /*!\ingroup Lattice
     * \brief A cell in an integer lattice. 
     */
    class LatticeCell {
      friend class LatticeBlockIterator;
      friend class LatticeMaskSet;
     public:
      /*!\brief Construct a cell of dimension \a n with lower corner at the origin. */
      explicit LatticeCell(dimension_type n=0) : _lower(n) { }
      /*!\brief Construct a cell with lower corner at \a l. */
      explicit LatticeCell(const IndexArray& l) : _lower(l) { }
      /*!\brief Construct a cell with lower corner at \a l. */
      explicit LatticeCell(const std::string& str);
  
      /*!\brief Copy constructor. */
      LatticeCell(const LatticeCell& c) : _lower(c._lower) { }
      /*!\brief Assignment operator. */
      LatticeCell& operator=(const LatticeCell& c) {
        if(this!=&c) { _lower=c._lower; } return *this; }

      /*!\brief Equality operator. */
      bool operator==(const LatticeCell& c) const { 
        return this->_lower==c._lower; }
      /*!\brief Inequality operator. */
      bool operator!=(const LatticeCell& c) const { 
        return !(*this==c); }
         
      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { return this->_lower.size(); }
      /*!\brief The \a i th interval. */
      Interval<index_type> operator[](dimension_type i) const {
        return Interval<index_type>(this->_lower[i],this->_lower[i]+1); } 
      /*!\brief The lower bound in the \a i th dimension. */
      index_type lower_bound(dimension_type i) const { return this->_lower[i]; }
      /*!\brief The upper bound in the \a i th dimension. */
      index_type upper_bound(dimension_type i) const { return this->_lower[i]+1; }

      /*!\brief The position of the lower corner in the latiice. */
      const IndexArray& position() const { return this->_lower; }
      /*!\brief The position of the lower corner in the latiice. */
      IndexArray lower() const { return this->_lower; }
      /*!\brief The position of the upper corner in the latiice. */
      IndexArray upper() const { 
        IndexArray result(this->dimension());
        for(dimension_type i=0; i!=this->dimension(); ++i) {
          result[i]=this->upper_bound(i); }
        return result;
      }
     private:
      IndexArray _lower;
    };
    

    /*!\ingroup Lattice
     * \brief A block of indices in a grid. */
    class LatticeBlock {
      friend class LatticeBlockIterator;
      friend class LatticeMaskSet;
     public:
      typedef LatticeBlockIterator const_iterator;
      
      /*!\brief Construct an empty lattice rectangle of dimension \a n. */
      explicit LatticeBlock(dimension_type n=0) : _lower(n), _upper(n) { }
      /*!\brief Construct a lattice rectangle specified by lower and upper corners. */
      explicit LatticeBlock(const IndexArray& l, const IndexArray& u)
        : _lower(l), _upper(u) { assert(l.size()==u.size()); }
      /*!\brief Construct a lattice rectangle defined by a string literal. */
      explicit LatticeBlock(const std::string& s);

      /*!\brief Convert from a lattice cell. */
      LatticeBlock(const LatticeCell& c)
        : _lower(c.lower()), _upper(c.upper()) { }
      /*!\brief Copy constructor. */
      LatticeBlock(const LatticeBlock& r)
        : _lower(r._lower), _upper(r._upper) { }
      
      /*!\brief Equality operator. */
      bool operator==(const LatticeBlock& other) const { 
        return this->_lower==other._lower && this->_upper==other._upper; }
      /*!\brief Inequality operator. */
      bool operator!=(const LatticeBlock& other) const { 
        return !(*this==other); }
        
      /*!\brief The dimension of the lattice rectangle. */
      dimension_type dimension() const { return this->_lower.size(); }
      /*!\brief Returns true if the lattice rectangle is empty. */
      bool empty() const;
      /*!\brief Returns true if the lattice rectangle has empty interior. */
      bool empty_interior() const;
      

      /*!\brief The \a i th interval. */
      Interval<index_type> operator[](dimension_type i) const {
        return Interval<index_type>(this->_lower[i],this->_upper[i]); } 
      /*!\brief The lower bound in the \a i th dimension. */
      index_type lower_bound(dimension_type i) const { return this->_lower[i]; }
      /*!\brief The upper bound in the \a i th dimension. */
      index_type upper_bound(dimension_type i) const { return this->_upper[i]; }

      /*!\brief Set the lower bound in the \a i th dimension to \a n. */
      void set_lower_bound(dimension_type i, index_type n) { this->_lower[i]=n; }
      /*!\brief Set the upper bound in the \a i th dimension to \a n. */
      void set_upper_bound(dimension_type i, index_type n) { this->_upper[i]=n; }

      /*!\brief The position of the lower corner in the lattice. */
      const IndexArray& lower() const { return this->_lower; }
      /*!\brief The position of the upper corner in the lattice. */
      const IndexArray& upper() const { return this->_upper; }

      /*!\brief The number of cells in each dimension. */
      SizeArray sizes() const;
      /*!\brief The product of the number of cells in each lower dimension. */
      SizeArray strides() const;
      /*!\brief The totel number of cells. */
      size_type size() const;
      
      /*!\brief A constant iterator to the lower cell in the lattice rectangle. */
      const_iterator begin() const;
      /*!\brief A constant iterator to the past-the-end cell of the lattice rectangle. */
      const_iterator end() const;
      
      /*!\brief The one-box neighbourhood of the rectangle. */
      LatticeBlock neighbourhood() const;

      /*! \brief Write to an output stream */
      std::ostream& write(std::ostream& os) const;
     private:
      IndexArray _lower;
      IndexArray _upper;
    };
    
    
    class LatticeBlockIterator 
      : public boost::iterator_facade<LatticeBlockIterator,
                                      LatticeCell,
                                      boost::forward_traversal_tag,
                                      LatticeCell>
    {
     public:
      LatticeBlockIterator(const IndexArray& l, const IndexArray& u)
        : _lower(l), _upper(u), _position(l) { }
      LatticeBlockIterator(const IndexArray& l, const IndexArray& u, const IndexArray& p)
        : _lower(l), _upper(u), _position(p) { }
     private:
      bool equal(const LatticeBlockIterator& other) const {
        return (this->_position==other._position) 
          && (this->_lower==other._lower) && (this->_upper==other._upper);
      }
      void increment() {
        dimension_type d=0;
        _position[d]+=1;
        while(_position[d]>=_upper[d] && (d+1u)!=_position.size() ) {
          _position[d]=_lower[d];
          d+=1;
          _position[d]+=1;
        }
      }
      LatticeCell dereference() const { return LatticeCell(_position); }
     private:
      friend class boost::iterator_core_access;
      dimension_type dimension() const { return _position.size(); }
      const IndexArray& _lower;
      const IndexArray& _upper;
      IndexArray _position;
    };

    inline
    LatticeBlock::const_iterator 
    LatticeBlock::begin() const 
    { 
      return const_iterator(this->_lower, this->_upper, this->_lower);
    }
    
    inline
    LatticeBlock::const_iterator 
    LatticeBlock::end() const 
    { 
      IndexArray end_position=this->_lower;
      end_position[this->dimension()-1]=_upper[this->dimension()-1];
      return const_iterator(this->_lower, this->_upper, end_position);
    }
    


    /*!\ingroup Lattice
     * \brief A list of cells in a lattice. */
    class LatticeCellListSet {
     public:
      typedef conversion_iterator<array_vector<index_type>::const_iterator, LatticeCell> iterator;
      typedef conversion_iterator<array_vector<index_type>::const_iterator, LatticeCell> const_iterator;
     public:
      /*!\brief Construct an empty list, to hold cells of dimension \a n. */
     LatticeCellListSet(dimension_type n=0) 
        : _list(n) { }

      /*!\brief Convert a single cell \a c into a one-element list of cells. */
      LatticeCellListSet(const LatticeCell& c);
      /*!\brief Convert a lattice rectangle \a r into a list of cells. */
      LatticeCellListSet(const LatticeBlock& r);
      /*!\brief Convert from a list of rectangles. */
      LatticeCellListSet(const LatticeBlockListSet& rls);
      /*!\brief Convert from a lattice set defined by a mask on a block of cells. */
      LatticeCellListSet(const LatticeMaskSet& ms);
      
      /*!\brief Copy constructor. */
      LatticeCellListSet(const LatticeCellListSet& cls) 
        : _list(cls._list) { }
      /*!\brief Assignment operator. */
      LatticeCellListSet& operator=(const LatticeCellListSet& cls) {
        if(this!=&cls) { this->_list=cls._list; } return *this; }
      
      /*!\brief The dimension of the cells in the list. */
      dimension_type dimension() const { return _list.array_size(); }
      /*!\brief True if the list is empty. */
      size_type empty() const { return _list.empty(); }
      /*!\brief The number of cells in the list. */
      size_type size() const { return _list.size(); }
      /*!\brief The \a i th cell in the list. */
      LatticeCell operator[] (size_type i) const { return LatticeCell(_list[i]); }
      
      /*!\brief A rectangular block containing all cells in the list. */
      LatticeBlock bounding_block() const;

      /*! \brief Adjoins a LatticeCell to the set. */
      void adjoin(const LatticeCell& c) { 
        assert(this->dimension() == c.dimension());
        this->_list.push_back(c.position()); 
      }
      /*! \brief Adjoins all cells in a LatticeBlock to the set. */
      void adjoin(const LatticeBlock& r);
      /*! \brief Adjoins a LatticeCellListSet to the set. */
      void adjoin(const LatticeCellListSet& cl);
      /*! \brief Adjoins a LatticeBlockListSet to the set. */
      void adjoin(const LatticeBlockListSet& rl);
      /*! \brief Adjoins a LatticeMaskSet to the set. */
      void adjoin(const LatticeMaskSet& ms);
      
      /*! \brief Empties the set. */
      void clear();

      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const { return const_iterator(_list.begin()); }
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const { return const_iterator(_list.end()); }
      
      /*! \brief Sorts the cells lexicographically, removing duplicates. */
      void unique_sort(); 

      /*! \brief Write to an output stream */
      std::ostream& write(std::ostream& os) const;
     private:
      array_vector<index_type> _list;
    };
      


    /*!\ingroup Lattice
     * \brief A list of rectangles in a lattice. */
    class LatticeBlockListSet {
     public:
      //      typedef LatticeBlockListSetIterator iterator;
      //typedef LatticeBlockListSetIterator const_iterator;
      typedef pair_constructor_iterator<array_vector<index_type>::const_iterator, LatticeBlock> iterator;
      typedef pair_constructor_iterator<array_vector<index_type>::const_iterator, LatticeBlock> const_iterator;
     public:
      /*!\brief Construct an empty list, to hold rectangular blocks of dimension \a n. */
      LatticeBlockListSet(dimension_type n) 
        : _list(n) { }

      /*!\brief Copy constructor. */
      LatticeBlockListSet(const LatticeBlockListSet& rls) 
        : _list(rls._list) { }
      /*!\brief Assignment operator. */
      LatticeBlockListSet& operator=(const LatticeBlockListSet& rls) {
        if(this!=&rls) { this->_list=rls._list; } return *this; }
        
      /*!\brief The dimension of the cells in the list. */
      dimension_type dimension() const { return _list.array_size(); }
      /*!\brief True if the list is empty. */
      size_type empty() const { return _list.empty(); }
      /*!\brief The number of rectangles in the list. */
      size_type size() const { return _list.size()/2; }
      /*!\brief The \a i th rectangle in the list. */
      LatticeBlock operator[] (size_type i) const { 
        return LatticeBlock(_list[2*i],_list[2*i+1]); 
      }

      /*!\brief A rectangular block containing all cells in the list. */
      LatticeBlock bounding_block() const;

      /*! \brief Adjoins a LatticeCell to the set. */
      void adjoin(const LatticeCell& c) { 
        assert(this->dimension() == c.dimension());
        this->_list.push_back(c.lower()); 
        this->_list.push_back(c.upper()); 
      }
      /*! \brief Adjoins all cells in a LatticeBlock to the set. */
      void adjoin(const LatticeBlock& r) { 
        assert(this->dimension() == r.dimension());
        this->_list.push_back(r.lower()); 
        this->_list.push_back(r.upper()); 
      }
      /*! \brief Adjoins a LatticeCellListSet to the set. */
      void adjoin(const LatticeCellListSet& cl);
      /*! \brief Adjoins a LatticeBlockListSet to the set. */
      void adjoin(const LatticeBlockListSet& rl);
      /*! \brief Adjoins a LatticeMaskSet to the set. */
      void adjoin(const LatticeMaskSet& ms);
      
      /*! \brief Empties the set. */
      void clear();

      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const  { return const_iterator(_list.begin()); }
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const { return const_iterator(_list.end()); }

      /*! \brief Write to an output stream */
      std::ostream& write(std::ostream& os) const;
     private:
      array_vector<index_type> _list;
    };
       
  


    /*!\ingroup Lattice
     * \brief A list of arrays of integers of the same size, representing rectangles in a grid. */
    class LatticeMaskSet {
     public:
      typedef mask_iterator<LatticeBlock::const_iterator,BooleanArray::const_iterator> iterator;
      typedef iterator const_iterator;
     public:
      /*! \brief Construct an lattice mask set on an empty block set of dimension \a n. */
      LatticeMaskSet(const size_type& n) 
        : _block(n), _mask(), _unbounded(false) { this->_compute_cached_attributes(); }
      /*! \brief Construct an empty lattice mask set on the block \a bb. */
      LatticeMaskSet(const LatticeBlock& bb) 
        : _block(bb), _mask(bb.size(),false), _unbounded(false) { this->_compute_cached_attributes(); }
      /*! \brief Construct a lattice mask set n the block \a bb with cells given by the mast \a ma. */
      LatticeMaskSet(const LatticeBlock& bb, const BooleanArray& ma) 
        : _block(bb), _mask(ma), _unbounded(false) { this->_compute_cached_attributes(); }
      /*! \brief Construct a lattice mask set n the block \a bb with cells given by the mast \a ma. */
      LatticeMaskSet(const LatticeBlock& bb, const BooleanArray& ma, const bool& ub) 
        : _block(bb), _mask(ma), _unbounded(ub) { this->_compute_cached_attributes(); }
      /*! \brief Construct a lattice mask oset n the block \a bb with cells given by the mast \a ma. */
      LatticeMaskSet(const LatticeBlock& bb, const LatticeCellListSet& cls);
      /*! \brief Construct a lattice mask oset n the block \a bb with cells given by the mast \a ma. */
      LatticeMaskSet(const LatticeBlock& bb, const LatticeBlockListSet& rls);

      /*! \brief Convert from a %LatticeCellListSet. */
      LatticeMaskSet(const LatticeCellListSet& cls);
      /*! \brief Convert from a %LatticeBlockListSet. */
      LatticeMaskSet(const LatticeBlockListSet& rls);
      /*! \brief Copy constructor. */
      LatticeMaskSet(const LatticeMaskSet& ms);
      
      /*!\brief The rectangular block of cells covered by the mask. */
      const LatticeBlock& block() const { return _block; }
      /*!\brief The number of cells in each dimension. */
      const SizeArray& sizes() const { return _sizes; }
      /*!\brief An array of boolean values determining whether the ith cell in the block is in the set. */
      const BooleanArray& mask() const { return _mask; };
            
      /*!\brief The dimension of the set. */
      dimension_type dimension() const { return _block.dimension(); }
      /*!\brief The number of cells in the block. */
      size_type capacity() const { return _mask.size(); }
      /*!\brief True if the set is empty. */
      size_type empty() const { return this->size()==0; }
      /*!\brief True if the set is bounded. */
      size_type bounded() const { return !this->_unbounded; }
      /*!\brief The number of cells in the set. */
      size_type size() const { return std::count(_mask.begin(),_mask.end(),true); }
      /*!\brief The \a i th cell in the set.<br>This method is expensive, since
       * it must search sequentially to find the next nonempty cell. */
      LatticeCell operator[](size_type i) const; 

      /*! \brief Compute the index of a cell in a grid. */
      size_type index(const LatticeCell& pos) const;

      /*! \brief Compute the cell at an index in a grid. */
      LatticeCell cell(size_type index) const;

      /*! \brief Empties the set. */
      void clear();
    
      /*! \brief Adjoins a LatticeCell to the set. */
      void adjoin(const LatticeCell& c) { 
        if (subset(c,this->block())) {
           this->_mask[this->index(c)]=true;
        }
      }
    
      /*! \brief Adjoins a LatticeBlock to the set. */
      void adjoin(const LatticeBlock& r);

      /*! \brief Adjoins a LatticeCellListSet to the set. */
      void adjoin(const LatticeCellListSet& cl);
      /*! \brief Adjoins a LatticeBlockListSet to the set. */
      void adjoin(const LatticeBlockListSet& rl);
      /*! \brief Adjoins a LatticeMaskSet to the set. */
      void adjoin(const LatticeMaskSet& ms);
        
      /*! \brief The one-box neighbourhood on the same grid. */
      LatticeMaskSet neighbourhood() const;
      /*! \brief The adjoining elements on the same grid. */
      LatticeMaskSet adjoining() const;
    
      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const { return const_iterator(_block.begin(),_mask.begin(),_mask.end()); }
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const { return const_iterator(_block.end(),_mask.end(),_mask.end()); }

      /*! \brief Write to an output stream */
      std::ostream& write(std::ostream& os) const;
     private:
      void _compute_cached_attributes();
      size_type _index(const IndexArray& cl) const;
     private:
      LatticeBlock _block;
      IndexArray _lower;
      IndexArray _upper;
      SizeArray _sizes;
      SizeArray _strides;
      BooleanArray _mask;
      bool _unbounded;
    };

    
    
    /*!\ingroup Lattice
     * \brief A transformation rule for lattice points, which maps rectangles to rectangles of the same dimension. 
     *
     * The transformation is described by an array of sequences \f$T_{i,j}\f$, 
     * \f$i\in[0,d)\f$ and \f$j\in\mathbb{Z}\f$, such that the image of vertex
     * \f$v\f$ is given by \f$(Tv)_i=T_{i,v[i]}\f$.
     */
    class LatticeTransformation {
      /*!\brief Construct from an array of sequences, one for each dimension. */
      LatticeTransformation(const array< std::vector<index_type> >& tr)
        : _transformation(tr) { }

      /*!\brief The dimension of lattice the transformation rule works on. */
      dimension_type dimension() const { return _transformation.size(); }
      /*!\brief The transformation applied to a point. */
      IndexArray operator() (const IndexArray& c) const;
      /*!\brief The transformation applied to a cell. */
      LatticeBlock operator() (const LatticeCell& c) const;
      /*!\brief The dimension of lattice the transformation rule works on. */
      LatticeBlock operator() (const LatticeBlock& r) const;
     private:
      array< std::vector<index_type> > _transformation;
    };
    
  
  }
}

#endif /* _ARIADNE_LATTICE_SET_H */
