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
#include "base/sequence.h"
#include "base/interval.h"
#include "base/binary_word.h"
#include "base/binary_tree.h"
#include "base/iterator.h"
#include "base/utility.h"

#include "geometry/geometry_declarations.h"
#include "geometry/rectangle.h"

namespace Ariadne {
  namespace Geometry {
    class UnitPartitionTreeCell;
    class UnitPartitionTree;
    class UnitPartitionTreeSet;
    
    std::ostream& operator<<(std::ostream&, const UnitPartitionTreeCell&);
    
    class LatticeMaskSet;
    class LatticeRectangle;
      
    class SubdivisionSequence
      : public sequence<dimension_type> 
    {
     public:
      SubdivisionSequence(const dimension_type& n)
        : _sequence(_default(n)), _dimension(n) { }
      
      template<typename FwdIter> 
      SubdivisionSequence(FwdIter b, FwdIter tb, FwdIter te)
        : _sequence(b,tb,te), _dimension(_compute_dimension())
      { }

      dimension_type dimension() const { 
        return _dimension; }
      dimension_type operator[](const size_type& i) const { 
        return _sequence[i]; }
     private:
      dimension_type _compute_dimension();
      sequence<dimension_type> _default(dimension_type n);
     private:
      sequence<dimension_type> _sequence;
      dimension_type _dimension;
    };

    class MaskedBinaryTree {
     public:
      typedef mask_iterator<BinaryTree::const_iterator, 
                            BooleanArray::const_iterator> const_iterator;
      typedef const_iterator iterator;

      MaskedBinaryTree() 
        : _tree(), _mask(1,false) { }
      MaskedBinaryTree(const BooleanArray& t, const BooleanArray& m) 
        : _tree(t), _mask(m) { }
      MaskedBinaryTree(const BinaryTree& t, const BooleanArray& m) 
        : _tree(t), _mask(m) { }

      const BinaryTree& tree() const { return _tree; }
      const BooleanArray& mask() const { return _mask; }
      
      size_type capacity() const { return mask().size(); }
      size_type size() const { return std::count(_mask.begin(),_mask.end(),true); }

      void reduce();

      const_iterator begin() const { 
        return const_iterator(_tree.begin(),_mask.begin(),_mask.end()); }
      const_iterator end() const { 
        return const_iterator(_tree.end(),_mask.end(),_mask.end()); }
     private:
      BinaryTree _tree;
      BooleanArray _mask;
    };

    class UnitPartitionTreeCell {
     public:
      typedef double dyadic_type;

      UnitPartitionTreeCell(const SubdivisionSequence& ss, 
                            const BinaryWord& bw);

      bool operator==(const UnitPartitionTreeCell& other) const {
        return this->_bounds==other._bounds; }

      dimension_type dimension() const { 
        return _bounds.dimension(); }
      const dyadic_type& lower_bound(dimension_type i) const {
        return _bounds.lower_bound(i); }
      const dyadic_type& upper_bound(dimension_type i) const {
        return _bounds.upper_bound(i); }
      const Rectangle<dyadic_type>& bounds() const { return _bounds; };
     private:
      void _compute_bounds(const SubdivisionSequence& ss, const BinaryWord& bw);
     private:
      Rectangle<dyadic_type> _bounds;
    };
    


    class UnitPartitionTree {
     public:
      typedef binary_constructor_iterator<BinaryTree::const_iterator, 
                                          UnitPartitionTreeCell, 
                                          SubdivisionSequence> const_iterator;
      typedef const_iterator iterator;
     
      UnitPartitionTree(const SubdivisionSequence& ss, 
                        const BinaryTree& bt);

      /*! \brief The space dimension of the tree. */
      dimension_type dimension() const { return _subdivisions.dimension(); }

      /*! \brief The sequence describing the order of subdivisions. */
      const SubdivisionSequence& subdivisions() const { return _subdivisions; }

      /*! \brief The array describing the tree. */
      const BinaryTree& binary_tree() const { return _tree; }

      /*! \brief The number of cells in the tree. */
      size_type size() const { return _tree.size(); }

      /*! \brief The depth of the smallest cell in the set. */
      size_type depth() const { return _tree.depth(); }
      /*! \brief The maximum number of subdivisions in each dimension. */
      SizeArray depths() const;      

      /*! \brief Constant iterator to the beginning of the cells in the tree. */
      const_iterator begin() const { return const_iterator(_subdivisions,_tree.begin()); }
      /*! \brief Constant iterator to the end of the cells in the tree. */
      const_iterator end() const { return const_iterator(_subdivisions,_tree.end()); }
     private:
      void reduce();
     private:
      SubdivisionSequence _subdivisions;
      BinaryTree _tree;
    };
    
    class UnitPartitionTreeSetIterator 
      : public boost::iterator_adaptor<UnitPartitionTreeSetIterator,
                                       mask_iterator<BinaryTree::const_iterator, BooleanArray::const_iterator>,
                                       UnitPartitionTreeCell,
                                       boost::use_default,
                                       UnitPartitionTreeCell>
      
    {
      typedef mask_iterator<BinaryTree::const_iterator, BooleanArray::const_iterator> Base;
     public:
      UnitPartitionTreeSetIterator(const SubdivisionSequence& ss,
                                   BinaryTree::const_iterator ti, 
                                   BooleanArray::const_iterator mi, 
                                   BooleanArray::const_iterator me) 
      //        : UnitPartitionTreeSetIterator::iterator_adaptor_(mask_iterator<BinaryTree::const_iterator, BooleanArray::const_iterator>(ti,mi,me)),
        : UnitPartitionTreeSetIterator::iterator_adaptor_(Base(ti,mi,me)),
          _subdivisions(ss)
      { }
     private:
      friend class boost::iterator_core_access;
      UnitPartitionTreeCell dereference() const { 
        return UnitPartitionTreeCell(_subdivisions,*this->base_reference()); 
      }
     private:
      SubdivisionSequence _subdivisions;
    };
  
    class UnitPartitionTreeSet {
     public:
      typedef binary_constructor_iterator<MaskedBinaryTree::const_iterator,
                                          UnitPartitionTreeCell,
                                          SubdivisionSequence> const_iterator;
      typedef const_iterator iterator;
     
      UnitPartitionTreeSet(const SubdivisionSequence& ss);
      UnitPartitionTreeSet(const SubdivisionSequence& ss, 
                           const BinaryTree& bt,
                           const BooleanArray& ba); 
      UnitPartitionTreeSet(const LatticeMaskSet& ms); 

     /*! \brief The space dimension of the tree. */
      dimension_type dimension() const { return _subdivisions.dimension(); }

      /*! \brief The sequence describing the order of subdivisions. */
      const SubdivisionSequence& subdivisions() const { return _subdivisions; }

      /*! \brief The array describing the tree. */
      const MaskedBinaryTree& words() const { return _words; }

      /*! \brief The array describing the tree. */
      const BinaryTree& binary_tree() const { return _words.tree(); }

      /*! \brief The array describing the tree. */
      const BooleanArray& mask() const { return _words.mask(); }

      /*! \brief The number of cells in the tree. */
      size_type capacity() const { return _words.capacity(); }

      /*! \brief The number of cells in the UnitPartitionTreeSet. */
      size_type size() const { return _words.size(); }
      
      /*! \brief The depth of the smallest cell in the set. */
      size_type depth() const { return binary_tree().depth(); }
      /*! \brief The maximum number of subdivisions in each dimension. */
      SizeArray depths() const;      

      /*! \brief Constant iterator to the beginning of the cells in the tree. */
      const_iterator begin() const { return const_iterator(_subdivisions,_words.begin()); }
      /*! \brief Constant iterator to the end of the cells in the tree. */
      const_iterator end() const { return const_iterator(_subdivisions,_words.end()); };
     private:
      void reduce() { this->_words.reduce(); }
     private:
      SubdivisionSequence _subdivisions;
      MaskedBinaryTree _words;
    };
    

    index_type compute_index(const SubdivisionSequence& ss, const BinaryWord& bw,const LatticeRectangle& r);
      
    LatticeRectangle 
    compute_block(const UnitPartitionTreeCell& c, 
                  const LatticeRectangle& r);

  }
  
  
}

#endif /* _ARIADNE_UNIT_PARTITION_TREE_SET_H */
