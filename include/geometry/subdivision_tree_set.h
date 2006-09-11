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
    class SubdivisionSequence;
    class SubdivisionTreeCell;
    class SubdivisionTree;
    class SubdivisionTreeSet;
    
    std::ostream& operator<<(std::ostream&, const SubdivisionTreeCell&);
    std::ostream& operator<<(std::ostream&, const SubdivisionSequence&);
    std::istream& operator>>(std::istream&, SubdivisionSequence&);
    
    class LatticeMaskSet;
    class LatticeRectangle;
      
    
    /*!\ingroup SubdivisionTree
     * \brief A binary tree with a boolean labelling on the leaves. */
    class MaskedBinaryTree {
     public:
      typedef mask_iterator<BinaryTree::const_iterator, 
                            BooleanArray::const_iterator> const_iterator;
      typedef const_iterator iterator;

      /*!\brief Construct a tree with a single leaf marked false. */
      MaskedBinaryTree()
        : _tree(), _mask(1,false) { }
      /*!\brief Construct a from a binary tree, labelling all leaves as false. */
      MaskedBinaryTree(const BinaryTree& t) 
        : _tree(t), _mask(t.size(),false) { }
      /*!\brief Construct a from the binary tree \a t, labelling leaves as given by the values in \a m. */
      MaskedBinaryTree(const BinaryTree& t, const BooleanArray& m) 
        : _tree(t), _mask(m) { }

      /*!\brief A constant reference to the binary tree. */
      const BinaryTree& tree() const { return _tree; }
      /*!\brief A constant reference to the mask on the leaves. */
      const BooleanArray& mask() const { return _mask; }
      
      /*!\brief The number of leaves of the tree. */
      size_type capacity() const { return mask().size(); }
      /*!\brief The number of leaves of the tree marked as true.  */
      size_type size() const { return std::count(_mask.begin(),_mask.end(),true); }

      /*!\brief Reduce the representation of the tree by combining elements
       * the same marking. 
       */
      void reduce();

      /*!\brief A constant iterator to the beginning of the marked leaves of the tree. */
      const_iterator begin() const { 
        return const_iterator(_tree.begin(),_mask.begin(),_mask.end()); }
      /*!\brief A constant iterator to the end of the marked leaves of the tree. */
      const_iterator end() const { 
        return const_iterator(_tree.end(),_mask.end(),_mask.end()); }
     private:
      BinaryTree _tree;
      BooleanArray _mask;
    };

    /*!\ingroup SubdivisionTree
     * \brief A sequence of coordinate giving axes of subdivision for a subdivision tree. 
     */
    class SubdivisionSequence
    {
     public:
      /*!\brief Construct the default sequence in dimension \a n, which consists of 
       * the sequence \f$0,1,2,\ldots,n-1\f$ repeated.
       */
      SubdivisionSequence(const dimension_type& n)
        : _sequence(_default(n)), _dimension(n) { }
      
      /*!\brief Construct from a sequence starting at \a b, where \a tb and \a te
       * describe the periodic tail of the sequence.
       */
      template<typename FwdIter> 
      SubdivisionSequence(FwdIter b, FwdIter tb, FwdIter te)
        : _sequence(b,tb,te), _dimension(_compute_dimension())
      { }

      /*!\brief Construct from a string literal of the for,
       * "\f$[a_1,a_2,\ldots,a_{k-1};b_1,\ldots,b_{l-1}]\f$", where
       * the values before the ';' denote the body of the sequence and the values
       * after denote the periodic tail.
       */
      SubdivisionSequence(const std::string& str);
       
      /*!\brief Equality operator. */
      bool operator==(const SubdivisionSequence& ss) const {
        return this->_sequence==ss._sequence && this->_dimension==ss._dimension; }

      /*!\brief Inequality operator. */
      bool operator!=(const SubdivisionSequence& ss) const {
        return !(*this==ss); }

      /*!\brief The dimension of the space the sequence describes subdivisions of. */
      dimension_type dimension() const { return _dimension; }
      /*!\brief The number of elements in the aperiodic body. */
      size_type body_size() const { return _sequence.body_size(); }
      /*!\brief The number of elements in the periodic tail. */
      size_type tail_size() const { return _sequence.tail_size(); }
      /*!\brief The \a i th element. */
      dimension_type operator[](const size_type& i) const { 
        return _sequence[i]; }
     private:
      dimension_type _compute_dimension();
      sequence<dimension_type> _default(dimension_type n);
      
      friend std::ostream& operator<<(std::ostream&, const SubdivisionSequence&);
      friend std::istream& operator>>(std::istream&, SubdivisionSequence&);
     private:
      sequence<dimension_type> _sequence;
      dimension_type _dimension;
    };

    /*!\ingroup SubdivisionTree
     * \brief A cell in a subdivision tree. */
    class SubdivisionTreeCell {
     public:
      /*!\brief The type used to represent the dyadic numbers giving the upper
       * and lower bounds of the cell. */
      typedef double dyadic_type;

      /*!\brief Construct from a sequence giving the subdivision dimensions, 
       * and a binary word giving the cell to be chosen at each subdivision. */
      SubdivisionTreeCell(const SubdivisionSequence& ss, 
                          const BinaryWord& bw);

      /*!\brief Equality operator. */
      bool operator==(const SubdivisionTreeCell& other) const {
        return this->_bounds==other._bounds; }

      /*!\brief Inequality operator. */
      bool operator!=(const SubdivisionTreeCell& other) const {
        return !(*this==other); }
        
      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { 
        return _bounds.dimension(); }
      /*!\brief A rectangle giving the cell. (Deprecated)
       *
       * \deprecated This class should be independent of the Geometry library.
       */
      const Rectangle<dyadic_type>& bounds() const { return _bounds; };
      /*!\brief The lower bound in the \a i th dimension. */
      const dyadic_type& lower_bound(dimension_type i) const {
        return _bounds.lower_bound(i); }
      /*!\brief The upper bound in the \a i th dimension. */
      const dyadic_type& upper_bound(dimension_type i) const {
        return _bounds.upper_bound(i); }
     private:
      void _compute_bounds(const SubdivisionSequence& ss, const BinaryWord& bw);
     private:
      Rectangle<dyadic_type> _bounds;
    };
    


    /*!\ingroup SubdivisionTree
     * \brief A subdivision structure on the unit hypercube determined by a sequence of subdivision coordinates and a binary tree. */
    class SubdivisionTree {
     public:
      typedef binary_constructor_iterator<BinaryTree::const_iterator, 
                                          SubdivisionTreeCell, 
                                          SubdivisionSequence> const_iterator;
      typedef const_iterator iterator;
     
      /*! Construct from a subdivision sequence and a binary tree. */
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
      /* Reduce the tree by combining cells where possible. */
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
  
    /*!\ingroup SubdivisionTree
     * \brief A subset of the unit hypercube described by a subdivision structure. */
    class SubdivisionTreeSet {
     public:
      typedef binary_constructor_iterator<MaskedBinaryTree::const_iterator,
                                          SubdivisionTreeCell,
                                          SubdivisionSequence> const_iterator;
      typedef const_iterator iterator;
     
      /*!\brief Construct an empty set with a single cell. */
      SubdivisionTreeSet(const SubdivisionSequence& ss);
      /*!\brief Construct an empty set with a single cell. */
      SubdivisionTreeSet(const SubdivisionSequence& ss, 
                         const BinaryTree& bt,
                         const BooleanArray& ba); 

      /*!\brief Copy constructor. */
      SubdivisionTreeSet(const SubdivisionTreeSet& sts);

      /*!\brief Convert from a lattice mask set \a ms. 
       *
       * The supporting block of \a ms must have sides which are a power of 2. */
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
