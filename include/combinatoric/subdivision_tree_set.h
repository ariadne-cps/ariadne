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

#ifndef ARIADNE_SUBDIVISION_TREE_SET_H
#define ARIADNE_SUBDIVISION_TREE_SET_H

#include <iosfwd>

#include "macros/precondition.h"

#include "base/types.h"
#include "base/array.h"
#include "base/sequence.h"
#include "base/iterator.h"

#include "combinatoric/binary_word.h"
#include "combinatoric/binary_tree.h"

namespace Ariadne {


    typedef BinaryWord word_type;

    class SubdivisionSequence;
    class SubdivisionCell;
    class SubdivisionBox;
    class SubdivisionTree;
    class SubdivisionCellListSet;
    class SubdivisionMaskSet;
    class SubdivisionTreeSet;
 
    class SubdivisionCellListSetIterator;
    class SubdivisionMaskSetIterator;
    class SubdivisionTreeSetIterator;
   
    std::ostream& operator<<(std::ostream&, const SubdivisionCell&);
    std::ostream& operator<<(std::ostream&, const SubdivisionBox&);
    std::ostream& operator<<(std::ostream&, const SubdivisionCellListSet&);
    std::ostream& operator<<(std::ostream&, const SubdivisionMaskSet&);
    std::ostream& operator<<(std::ostream&, const SubdivisionTreeSet&);
    std::ostream& operator<<(std::ostream&, const SubdivisionSequence&);
    
    class LatticeMaskSet;
    class LatticeBlock;
      
    struct dyadic_interval { double lower, upper; };
    inline bool operator==(const dyadic_interval& ivl1, const dyadic_interval& ivl2) {
      return ivl1.lower==ivl2.lower && ivl1.upper==ivl2.upper; }
    inline bool operator!=(const dyadic_interval& ivl1, const dyadic_interval& ivl2) {
      return ivl1.lower!=ivl2.lower || ivl1.upper!=ivl2.upper; }
 
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
     * \brief A sequence of coordinates giving axes of subdivision for a subdivision tree. 
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
      template<class FwdIter> 
      SubdivisionSequence(FwdIter b, FwdIter tb, FwdIter te)
        : _sequence(b,tb,te), _dimension(_compute_dimension())
      { }

      /* \brief Construct from a string literal of the for,
       * "\f$[a_1,a_2,\ldots,a_{k-1};b_1,\ldots,b_{l-1}]\f$", where
       * the values before the ';' denote the body of the sequence and the values
       * after denote the periodic tail.
       */
       
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
     private:
      sequence<dimension_type> _sequence;
      dimension_type _dimension;
    };




    /*!\ingroup SubdivisionTree
     * \brief A cell in a subdivision tree. 
     *
     * A %SubdivisionCell is defined by a BinaryWord \a w and a SubdivisionSequence \a ss.
     * Start with a unit cell \f$[0,1]^d\f$ in \f$\mathbb{R}^d\f$ and successively partition.
     * At the \a i th step, we subdivide in coordinate \a k=ss[i]. 
     * If \f$bw[i]=0\f$, then the \a k th coordinate \f$[a_k,b_k]\f$ becomes \f$[a_k,c_k]\f$, 
     * where \f$c_k=(a_k+b_k)/2\f$.
     * If \f$bw[i]=1\f$, then the \a k th coordinate becomes \f$[c_k,b_k]\f$.
     *
     * Since the coordinates are described by dyadic numbers in the unit interval,
     * they can be exactly computed and represented by double-precision floating point
     * numbers up to about 50 subdivisions.         
     */
    class SubdivisionCell {
      friend class SubdivisionBox;
      friend class SubdivisionTreeSet;
     public:
      /*!\brief The type used to represent the dyadic numbers giving the upper
       * and lower bounds of the cell. */

      /*!\brief Construct from a sequence giving the subdivision dimensions, 
       * and a binary word giving the cell to be chosen at each subdivision. */
      SubdivisionCell(const SubdivisionSequence& ss, 
                          const BinaryWord& bw);

      /*!\brief Construct from a dimension, a depth and an integer giving the position among cells of that depth. */
      SubdivisionCell(const dimension_type& dim, const depth_type& dpth, const size_type& pos); 

      /*!\brief Construct from a dimension, and a word giving the subdivisions. */
      SubdivisionCell(const dimension_type& dim, const word_type& word); 

      /*!\brief Equality operator. */
      bool operator==(const SubdivisionCell& other) const {
        return this->_bounds==other._bounds; }

      /*!\brief Inequality operator. */
      bool operator!=(const SubdivisionCell& other) const {
        return !(*this==other); }
        
      /*!\brief The dimension of the cell. */
      dimension_type dimension() const { 
        return this->_bounds.size(); }

      /*!\brief The depth of the cell in the tree. */
      size_type depth() const { 
        return this->_subdivisions.size(); }

      /*!\brief The subdivisions needed to get to the cell. */
      BinaryWord subdivisions() const { 
        return this->_subdivisions; }

      /*!\brief The box represented by the cell. */
      SubdivisionBox box() const;

      /*!\brief Split the cell. */
      SubdivisionCell split(bool lr) const { 
        SubdivisionCell result=*this; result._split(lr); return result; }

      /*!\brief The lower bound in the \a i th dimension. */
      const dyadic_type& lower_bound(dimension_type i) const {
        return this->_bounds[i].lower; }

      /*!\brief The upper bound in the \a i th dimension. */
      const dyadic_type& upper_bound(dimension_type i) const {
        return this->_bounds[i].upper; }
        
      /*! The volume of the cell as a fraction of the unit cell. */
      dyadic_type volume() const { return 1<<this->depth(); }
     private:
      void _split(bool lr);
     private:
      void _compute_word(const depth_type& dpth, const size_type& pos);
      void _compute_bounds(const SubdivisionSequence& ss, const BinaryWord& bw);
      void _compute_bounds(const dimension_type& dim, const word_type& bw);
     private:
      word_type _subdivisions;
      array<dyadic_interval> _bounds;
    };
    


    /*!\ingroup SubdivisionTree
     * \brief A box in a subdivision tree. 
     *
     * A %SubdivisionBox is defined by an array of intervals with dyadic coefficients. 
     */
    class SubdivisionBox {
      friend SubdivisionBox hull(const SubdivisionBox& sbx1, const SubdivisionBox& sbx2);
     public:
      SubdivisionBox(const dimension_type& dim) : _depth(0), _bounds(dim) { 
        for(dimension_type i=0; i!=dim; ++i) { this->_bounds[i].lower=0; this->_bounds[i].upper=1; } }
      SubdivisionBox(const SubdivisionCell& c) : _depth(c.depth()), _bounds(c._bounds) { }
      dimension_type dimension() const { return this->_bounds.size(); }
      depth_type depth() const { return this->_depth; }
      dyadic_type lower_bound(dimension_type i) const { return this->_bounds[i].lower; }
      dyadic_type upper_bound(dimension_type i) const { return this->_bounds[i].upper; }
      dyadic_type volume() const;
     private:
      depth_type _depth;
      array<dyadic_interval> _bounds;
    };

    SubdivisionBox hull(const SubdivisionBox& sbx1, const SubdivisionBox& sbx2);


    /*!\ingroup SubdivisionCellListSet
     * \brief A list of subdivision cells.
     *
     * A %SubdivisionCellListSet represents a list of subdivision cells of the 
     * same dimension using efficient arrayed storage. 
     */
    class SubdivisionCellListSet
    {
     public:
      typedef SubdivisionCellListSetIterator const_iterator;
      typedef SubdivisionCellListSetIterator iterator;
      SubdivisionCellListSet() : _dimension() { }
      SubdivisionCellListSet(dimension_type d) : _dimension(d) { }
      dimension_type dimension() const { return this->_dimension; }
      size_type size() const { return this->_words.size(); }
      SubdivisionCell pop() { 
        const word_type& word=this->_words.back();
        this->_words.pop_back(); 
        return SubdivisionCell(this->_dimension,word); }
      void adjoin(const SubdivisionCell& c) {
        if(this->_words.empty()) { this->_dimension=c.dimension(); }
        ARIADNE_PRECONDITION(this->_dimension==c.dimension());
        _words.push_back(c.subdivisions()); }
      void unique_sort();
      const_iterator begin() const;
      const_iterator end() const;
     private:
      friend std::ostream& operator<<(std::ostream&, const SubdivisionCellListSet&);
     private:
      dimension_type _dimension;
      std::vector<word_type> _words;
    };


    /*!\ingroup SubdivisionTree
     * \brief A subdivision structure on the unit hypercube determined by a sequence of subdivision coordinates and a binary tree. 
     *
     * A %SubdivisionTree represents a partition of the unit hypercube into hypercubes described by successive subdivisions 
     * parallel to the coordinate axes.
     *
     * A subdivision into \a n cells uses 2n+1 bits of data, independently of the dimension of the set. This means that
     * a %SubdivisionTree provides a highly memory-efficient adaptive partitioning scheme.
     */
    class SubdivisionTree {
     public:
      typedef binary_constructor_iterator<BinaryTree::const_iterator, 
                                          SubdivisionCell, 
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
    

    /*!\ingroup SubdivisionMaskSet
     * \brief A subset of the unit hypercube described by a subdivision structure with constant depth. 
     *
     * The \em depth of the set is the maximum number of subdivisions needed to specify the smallest cell.
     *
     * The operations of intersection, union and set difference can all be performed in
     * linear time on the capacity (NOT the size) of the two operands. This makes a the set slow for describing
     * very sparse sets. However, adjoining and removing cells and testing for inclusing only require constant time. 
     * 
     * A %SubdivisionMaskSet defined on a depth-n partition uses \f$2^n\f$ bits of data.
     */
    class SubdivisionMaskSet {
     public:
      typedef SubdivisionMaskSetIterator iterator;
      typedef SubdivisionMaskSetIterator const_iterator;

      /*!\brief Construct an empty set with the given dimension and depth. */
      SubdivisionMaskSet(const dimension_type& dim, const depth_type& depth);

      /*! \brief The space dimension of the tree. */
      dimension_type dimension() const { return this->_dimension; }

      /*! \brief The depth of the smallest cell in the set. */
      size_type depth() const { return this->_depth; }

      /*! \brief The array describing the tree. */
      const array<bool>& mask() const { return _mask; }

      /*! \brief The number of cells in the tree. */
      size_type capacity() const { return _mask.size(); }

      /*! \brief The number of cells in the SubdivisionTreeSet. */
      size_type size() const { return std::count(_mask.begin(),_mask.end(),true); }
      
      /*! \brief The volume of the set. */
      dyadic_type volume() const { return this->size() / this->capacity(); }

      /*! \brief Adjoin a cell. */
      void adjoin(const SubdivisionCell& c) {
        if(c.depth()==this->depth()) { this->_mask[this->_position(c)]=true; }
        else { this->_adjoin_block(c); } }
      
      /*! \brief Remove a cell. */
      void remove(const SubdivisionCell& c) {
        if(c.depth()==this->depth()) { this->_mask[this->_position(c)]=false; }
        else { this->_remove_block(c); } }
      
      /*! \brief Constant iterator to the beginning of the cells in the tree. */
      const_iterator begin() const;
      /*! \brief Constant iterator to the end of the cells in the tree. */
      const_iterator end() const;
     private:
      size_type _position(const SubdivisionCell& c) const;
      void _adjoin_block(const SubdivisionCell& c);
      void _remove_block(const SubdivisionCell& c);
      void _restrict_block(const SubdivisionCell& c);
     private:
      uint _dimension;
      uint _depth; // store depth for quick access
      BooleanArray _mask;
    };
    
    /*!\ingroup SubdivisionTree
     * \brief A subset of the unit hypercube described by a subdivision structure. 
     *
     * The \em depth of the set is the maximum number of subdivisions needed to specify the smallest cell.
     *
     * The operations of intersection, union and set difference can all be performed in
     * linear time on the two operands. However, adjoining a new cell also requires
     * linear time, which means that it is often preferable to store cells to be
     * added in a LatticeCellListSet and adjoin them all at once.
     *
     * A %SubdivisionTreeSet defined on an n-cell partition uses 3n+1 bits of data.
     * This means that a %SubdivisionTreeSet provides a highly efficient representation of structured data sets.
     * In the worst case where every cell of a depth-\a d subdivision needs to be
     * explicitly specified, a %SubdivisionTreeSet uses three times as much memory as a LatticeMaskSet, but in 
     * typical cases a %SubdivisionTreeSet performs much better.
     */
    class SubdivisionTreeSet {
     public:
      typedef binary_constructor_iterator<MaskedBinaryTree::const_iterator,
                                          SubdivisionCell,
                                          SubdivisionSequence> const_iterator;
      typedef const_iterator iterator;
     
      typedef double dyadic_type;
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
      const_iterator end() const { return const_iterator(_subdivisions,_words.end()); }
      
      /*! \brief The volume of the set. */
      dyadic_type volume() const {
        dyadic_type result=0;
        for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
          result+=iter->volume();
        }
        return result;
      }
     private:
      void reduce() { this->_words.reduce(); }
     private:
      SubdivisionSequence _subdivisions;
      MaskedBinaryTree _words;
    };
    

    IndexArray 
    compute_position(const SubdivisionSequence& ss, 
                     const BinaryWord& bw,
                     const LatticeBlock& r);
      
    LatticeBlock 
    compute_block(const SubdivisionCell& c, 
                  const LatticeBlock& r);


    inline 
    SubdivisionBox 
    SubdivisionCell::box() const {
      return SubdivisionBox(*this); 
    }


  
} // namespace Ariadne



#include "subdivision_tree_set_iterators.h"

#endif /* ARIADNE_SUBDIVISION_TREE_SET_H */
