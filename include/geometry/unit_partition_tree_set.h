/***************************************************************************
 *            unit_partition_tree_set.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

/*! \file unit_partition_tree_set.h
 *  \brief Cuboidal partition trees on a unit cuboid.
 */

#ifndef _ARIADNE_UNIT_PARTITION_TREE_SET_H
#define _ARIADNE_UNIT_PARTITION_TREE_SET_H

#include <iosfwd>
#include <algorithm>
#include <vector>
#include <iosfwd>

#include "base/basic_type.h"
#include "base/array.h"
#include "base/interval.h"
#include "base/binary_word.h"
#include "base/binary_tree.h"
#include "base/utility.h"

#include "geometry/geometry_declarations.h"

namespace Ariadne {
  namespace Geometry {

    class UnitPartitionTreeCell;
    class UnitPartitionTreeSet;
    class UnitPartitionTreeSetIterator;
    
    std::ostream& operator<<(std::ostream&, const UnitPartitionTreeCell&);
    
    class UnitGridMaskSet;
    class UnitGridRectangle;
      
    //TODO: Replace by a general-purpose mask iterator
    class BinarySubtreeIterator {
      typedef BinaryWord value_type;
      typedef BinaryWord reference;
     public:
      BinarySubtreeIterator(BinaryTree::const_iterator ti, 
                            BooleanArray::const_iterator mi, 
                            BooleanArray::const_iterator me) 
        : _word_iter(ti), _mask_iter(mi), _mask_end(me) 
      { 
        while(_mask_iter!=_mask_end && !*_mask_iter) { 
          ++_mask_iter; ++_word_iter; 
        }
      }
      
      bool operator==(const BinarySubtreeIterator& other) const { 
        return this->equal(other); }
      bool operator!=(const BinarySubtreeIterator& other) const { 
        return !this->equal(other); }
      BinaryWord operator*() const { 
        return this->dereference(); }
      BinarySubtreeIterator& operator++() { 
        this->increment(); return *this; }
     private:
      bool equal(const BinarySubtreeIterator& other) const {
        return this->_word_iter==other._word_iter && this->_mask_iter==other._mask_iter;
      }
      BinaryWord dereference() const { 
        return _word_iter.operator*(); 
      }
      void increment() { 
        do { 
          ++_mask_iter; ++_word_iter; } 
        while(_mask_iter!=_mask_end && !*_mask_iter);
      }
     private:
      BinaryTree::const_iterator _word_iter;
      BooleanArray::const_iterator _mask_iter;
      BooleanArray::const_iterator _mask_end;
    };

    
    
    class UnitPartitionTreeCell {
      typedef double dyadic_type;
     public:
      UnitPartitionTreeCell(const dimension_type& d, const SubdivisionSequence& ss, 
                            const BinaryWord& bw);

      bool operator==(const UnitPartitionTreeCell& other) const {
        return this->_bounds==other._bounds; }

      const SubdivisionSequence& subdivisions() const { 
        return _subdivisions; }
      BinaryWord word() const { 
        return _word; }
      dimension_type dimension() const { 
        return _bounds.size()/2; }
      Interval<dyadic_type> operator[](dimension_type i) const { 
        return Interval<dyadic_type>(_bounds[2*i],_bounds[2*i+1]); }
      dyadic_type lower_bound(dimension_type i) const { return _bounds[2*i]; }
      dyadic_type upper_bound(dimension_type i) const { return _bounds[2*i+1]; }
      Rectangle<dyadic_type> bounds() const;
     private:
      const SubdivisionSequence _subdivisions;
      BinaryWord _word;
      array<dyadic_type> _bounds;
    };
    
    class UnitPartitionTreeSet {
     public:
      typedef UnitPartitionTreeSetIterator iterator;
      typedef UnitPartitionTreeSetIterator const_iterator;
     
      UnitPartitionTreeSet(const SubdivisionSequence& ss);
      UnitPartitionTreeSet(const SubdivisionSequence& ss, 
                           const BinaryTree& bt,
                           const BooleanArray& ba); 
      UnitPartitionTreeSet(const UnitGridMaskSet& ms); 

     /*! \brief The space dimension of the tree. */
      dimension_type dimension() const { return _dimension; }

      /*! \brief The sequence describing the order of subdivisions. */
      const SubdivisionSequence& subdivisions() const { return _subdivisions; }

      /*! \brief The array describing the tree. */
      const BinaryTree& tree() const { return _tree; }

      /*! \brief The array describing the tree. */
      const BooleanArray& mask() const { return _mask; }

      /*! \brief The number of cells in the tree. */
      size_type capacity() const { return _mask.size(); }

      /*! \brief The number of cells in the UnitPartitionTreeSet. */
      size_type size() const { return std::count(_mask.begin(),_mask.end(), true); }
      
      /*! \brief The depth of the smallest cell in the set. */
      size_type depth() const { return tree().depth(); }
      /*! \brief The maximum number of subdivisions in each dimension. */
      SizeArray depths() const;      

      /*! \brief Constant iterator to the beginning of the cells in the tree. */
      const_iterator begin() const;
      /*! \brief Constant iterator to the end of the cells in the tree. */
      const_iterator end() const;
     private:
      void reduce();
     private:
      dimension_type _dimension;
      SubdivisionSequence _subdivisions;
      BinaryTree _tree;
      BooleanArray _mask;
    };
    


    class UnitPartitionTreeSetIterator {
      typedef UnitPartitionTreeCell value_type;
      typedef UnitPartitionTreeCell reference;
     public:
      UnitPartitionTreeSetIterator(const dimension_type& d,
                                const SubdivisionSequence& ss,
                                BinaryTree::const_iterator ti, 
                                BooleanArray::const_iterator mi, 
                                BooleanArray::const_iterator me) 
        : _dimension(d), _subdivisions(ss), _word_iter(ti), _mask_iter(mi), _mask_end(me) 
      { 
        while(_mask_iter!=_mask_end && !*_mask_iter) { 
          ++_mask_iter; ++_word_iter; 
        }
      }
      
      bool operator==(const UnitPartitionTreeSetIterator& other) const { 
        return this->equal(other); }
      bool operator!=(const UnitPartitionTreeSetIterator& other) const { 
        return !this->equal(other); }
      UnitPartitionTreeCell operator*() const { 
        return this->dereference(); }
      UnitPartitionTreeSetIterator& operator++() { 
        this->increment(); return *this; }
     private:
      bool equal(const UnitPartitionTreeSetIterator& other) const {
        return this->_word_iter==other._word_iter && this->_mask_iter==other._mask_iter;
      }
      UnitPartitionTreeCell dereference() const { 
        return UnitPartitionTreeCell(_dimension,_subdivisions,*_word_iter); 
      }
      void increment() { 
        do { 
          ++_mask_iter; ++_word_iter; } 
        while(_mask_iter!=_mask_end && !*_mask_iter);
      }
     private:
      dimension_type _dimension;
      SubdivisionSequence _subdivisions;
      BinaryTree::const_iterator _word_iter;
      BooleanArray::const_iterator _mask_iter;
      BooleanArray::const_iterator _mask_end;
    };
  
    inline 
    UnitPartitionTreeSet::const_iterator 
    UnitPartitionTreeSet::begin() const 
    { 
      return const_iterator(_dimension,_subdivisions,_tree.begin(),_mask.begin(),_mask.end()); 
    }

    inline 
    UnitPartitionTreeSet::const_iterator 
    UnitPartitionTreeSet::end() const 
    { 
      return const_iterator(_dimension,_subdivisions,_tree.end(),_mask.end(),_mask.end()); 
    }

    index_type compute_index(const SubdivisionSequence& ss, const BinaryWord& bw,const UnitGridRectangle& r);
      
    UnitGridRectangle 
    compute_block(const UnitPartitionTreeCell& c, 
                  const UnitGridRectangle& r);

  }
  
  
}

#endif /* _ARIADNE_UNIT_PARTITION_TREE_SET_H */
