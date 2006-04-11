/***************************************************************************
 *            subdivision_tree_set.h
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

/*! \file subdivision_tree_set.h
 *  \brief Cuboidal partition trees on a unit cuboid.
 */

#ifndef _ARIADNE_SUBDIVISION_TREE_SET_H
#define _ARIADNE_SUBDIVISION_TREE_SET_H

#include <iosfwd>

#include "../declarations.h"

#include "../base/array.h"
#include "../base/sequence.h"
#include "../base/binary_word.h"
#include "../base/binary_tree.h"
#include "../base/iterator.h"

#include "../geometry/rectangle.h"

namespace Ariadne {
  namespace Geometry {
    class SubdivisionTreeCell;
    class SubdivisionTree;
    class SubdivisionTreeSet;
    
    std::ostream& operator<<(std::ostream&, const SubdivisionTreeCell&);
    
    class LatticeMaskSet;
    class LatticeRectangle;
      
    /*! \brief A sequence of coordinate giving axes of subdivision for a subdivision tree. */
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

    /*! \brief A binary tree with a boolean labelling on the leaves. */
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

    /*! \brief A cell in a subdivision tree. */
    class SubdivisionTreeCell {
     public:
      typedef double dyadic_type;

      SubdivisionTreeCell(const SubdivisionSequence& ss, 
                            const BinaryWord& bw);

      bool operator==(const SubdivisionTreeCell& other) const {
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
    


    /*! \brief A subdivision structure on the unit hypercube determined by a sequence of subdivision coordinates and a binary tree. */
    class SubdivisionTree {
     public:
      typedef binary_constructor_iterator<BinaryTree::const_iterator, 
                                          SubdivisionTreeCell, 
                                          SubdivisionSequence> const_iterator;
      typedef const_iterator iterator;
     
      SubdivisionTree(const SubdivisionSequence& ss, 
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
    
    class SubdivisionTreeSetIterator 
      : public boost::iterator_adaptor<SubdivisionTreeSetIterator,
                                       mask_iterator<BinaryTree::const_iterator, BooleanArray::const_iterator>,
                                       SubdivisionTreeCell,
                                       boost::use_default,
                                       SubdivisionTreeCell>
      
    {
      typedef mask_iterator<BinaryTree::const_iterator, BooleanArray::const_iterator> Base;
     public:
      SubdivisionTreeSetIterator(const SubdivisionSequence& ss,
                                   BinaryTree::const_iterator ti, 
                                   BooleanArray::const_iterator mi, 
                                   BooleanArray::const_iterator me) 
      //        : SubdivisionTreeSetIterator::iterator_adaptor_(mask_iterator<BinaryTree::const_iterator, BooleanArray::const_iterator>(ti,mi,me)),
        : SubdivisionTreeSetIterator::iterator_adaptor_(Base(ti,mi,me)),
          _subdivisions(ss)
      { }
     private:
      friend class boost::iterator_core_access;
      SubdivisionTreeCell dereference() const { 
        return SubdivisionTreeCell(_subdivisions,*this->base_reference()); 
      }
     private:
      SubdivisionSequence _subdivisions;
    };
  
    /*! \brief A subset of the unit hypercube described by a subdivision structure. */
    class SubdivisionTreeSet {
     public:
      typedef binary_constructor_iterator<MaskedBinaryTree::const_iterator,
                                          SubdivisionTreeCell,
                                          SubdivisionSequence> const_iterator;
      typedef const_iterator iterator;
     
      SubdivisionTreeSet(const SubdivisionSequence& ss);
      SubdivisionTreeSet(const SubdivisionSequence& ss, 
                           const BinaryTree& bt,
                           const BooleanArray& ba); 
      SubdivisionTreeSet(const LatticeMaskSet& ms); 

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

      /*! \brief The number of cells in the SubdivisionTreeSet. */
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
    

    IndexArray 
    compute_position(const SubdivisionSequence& ss, 
                     const BinaryWord& bw,
                     const LatticeRectangle& r);
      
    LatticeRectangle 
    compute_block(const SubdivisionTreeCell& c, 
                  const LatticeRectangle& r);

  }
  
  
}

#endif /* _ARIADNE_SUBDIVISION_TREE_SET_H */
