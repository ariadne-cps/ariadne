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
#include "../base/array_operations.h"

namespace Ariadne {
  namespace Geometry {
    class LatticeCell;

    class LatticeRectangle;
    class LatticeRectangleIterator;

    class LatticeCellListSet;
    class LatticeCellListSetIterator;

    class LatticeRectangleListSet;
    class LatticeRectangleListSetIterator;

    class LatticeMaskSet;
    class LatticeMaskSetIterator;

    bool interiors_intersect(const LatticeRectangle&, const LatticeRectangle&);
    bool interiors_intersect(const LatticeRectangle&, const LatticeMaskSet&);
    bool interiors_intersect(const LatticeMaskSet&, const LatticeMaskSet&);
    bool subset(const LatticeCell&, const LatticeRectangle&);
    bool subset(const LatticeRectangle&, const LatticeRectangle&);
    bool subset(const LatticeRectangle&, const LatticeMaskSet&);
    bool subset(const LatticeMaskSet&, const LatticeMaskSet&);
    
    LatticeRectangle regular_intersection(const LatticeRectangle&, const LatticeRectangle&);
    LatticeMaskSet regular_intersection(const LatticeMaskSet&, const LatticeMaskSet&);
    LatticeMaskSet join(const LatticeMaskSet&, const LatticeMaskSet&);
    LatticeMaskSet difference(const LatticeMaskSet&, const LatticeMaskSet&);

    bool operator<(const LatticeCell& lc1, const LatticeCell& lc2);
    
    std::istream& operator>>(std::istream& os, LatticeRectangle&);
    
    std::ostream& operator<<(std::ostream& os, const LatticeCell&);
    std::ostream& operator<<(std::ostream& os, const LatticeRectangle&);
    std::ostream& operator<<(std::ostream& os, const LatticeMaskSet&);
    std::ostream& operator<<(std::ostream& os, const LatticeCellListSet&);
    std::ostream& operator<<(std::ostream& os, const LatticeRectangleListSet&);


    /*! \brief A transformation rule for lattice points. */
    class LatticeTransformation {
      LatticeTransformation(const array< std::vector<index_type> >& tr)
        : _transformation(tr) { }

      dimension_type dimension() const { return _transformation.size(); }
      LatticeRectangle operator() (const LatticeCell& c) const;
      LatticeRectangle operator() (const LatticeRectangle& r) const;
     private:
      array< std::vector<index_type> > _transformation;
    };

 
      
 


    /*! \brief A cell in a unit grid. */
    class LatticeCell {
      friend class LatticeRectangleIterator;
      friend class LatticeMaskSet;
     public:
      LatticeCell(dimension_type n) : _lower(n) { }
      LatticeCell(const IndexArray& l) : _lower(l) { }
      LatticeCell(const LatticeCell& c) : _lower(c._lower) { }
      
      /*!\brief Addignment operator. */
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
    
    /*! \brief A block of indices in a grid. */
    class LatticeRectangle {
      friend class LatticeRectangleIterator;
      friend class LatticeMaskSet;
     public:
      typedef LatticeRectangleIterator const_iterator;
      
      /*!\brief An empty lattice rectangle. */
      explicit LatticeRectangle() : _lower(), _upper() { }
      /*!\brief A lattice rectangle of dimension \a n. */
      explicit LatticeRectangle(dimension_type n) : _lower(n), _upper(n) { }
      /*!\brief A lattice rectangle specified by lower and upper corners. */
      explicit LatticeRectangle(const IndexArray& l, const IndexArray& u)
        : _lower(l), _upper(u) { }
      /*!\brief A lattice rectangle defined by a string literal. */
      explicit LatticeRectangle(const std::string& s);
      /*!\brief Convert from a lattice cell. */
      LatticeRectangle(const LatticeCell& c)
        : _lower(c.lower()), _upper(c.upper()) { }
      /*!\brief Copy constructor. */
      LatticeRectangle(const LatticeRectangle& r)
        : _lower(r._lower), _upper(r._upper) { }
      
      /*!\brief Equality operator. */
      bool operator==(const LatticeRectangle& other) const { 
        return this->_lower==other._lower && this->_upper==other._upper; }
        
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

      /*!\brief The position of the lower corner in the latiice. */
      const IndexArray& lower() const { return this->_lower; }
      /*!\brief The position of the upper corner in the latiice. */
      const IndexArray& upper() const { return this->_upper; }

      /*! \brief The number of cells in each dimension. */
      SizeArray sizes() const;
      /*! \brief The product of the number of cells in each lower dimension. */
      SizeArray strides() const;
      /*! \brief The totel number of cells. */
      size_type size() const;
      
      /*!\brief A constant iterator to the lower cell in the lattice rectangle. */
      const_iterator begin() const;
      /*!\brief A constant iterator to the past-the-end cell of the lattice rectangle. */
      const_iterator end() const;
     private:
      IndexArray _lower;
      IndexArray _upper;
    };
    
    /*! \brief An iterator for positions in rectangular piece of a grid. */
    class LatticeRectangleIterator 
      : public boost::iterator_facade<LatticeRectangleIterator,
                                      LatticeCell,
                                      boost::forward_traversal_tag,
                                      LatticeCell>
    {
     public:
      LatticeRectangleIterator(const IndexArray& l, const IndexArray& u)
        : _lower(l), _upper(u), _position(l) { }
      LatticeRectangleIterator(const IndexArray& l, const IndexArray& u, const IndexArray& p)
        : _lower(l), _upper(u), _position(p) { }
     private:
      bool equal(const LatticeRectangleIterator& other) const {
        return (this->_position==other._position) 
          && (this->_lower==other._lower) && (this->_upper==other._upper);
      }
      void increment() {
        dimension_type d=0;
        _position[d]+=1;
        while(_position[d]==_upper[d] && (d+1u)!=_position.size() ) {
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
    LatticeRectangle::const_iterator 
    LatticeRectangle::begin() const 
    { 
      return const_iterator(this->_lower, this->_upper, this->_lower);
    }
    
    inline
    LatticeRectangle::const_iterator 
    LatticeRectangle::end() const 
    { 
      IndexArray end_position=this->_lower;
      end_position[this->dimension()-1]=_upper[this->dimension()-1];
      return const_iterator(this->_lower, this->_upper, end_position);
    }
    


    /*!\brief A list of cells in a lattice. */
    class LatticeCellListSet {
     public:
      typedef conversion_iterator<array_vector<index_type>::const_iterator, LatticeCell> iterator;
      typedef conversion_iterator<array_vector<index_type>::const_iterator, LatticeCell> const_iterator;
     public:
      LatticeCellListSet(dimension_type n) 
        : _list(n) { }

      LatticeCellListSet(const LatticeCell& c);
      LatticeCellListSet(const LatticeRectangle& r);
      LatticeCellListSet(const LatticeMaskSet& ms);
      LatticeCellListSet(const LatticeRectangleListSet& rls);
      
      LatticeCellListSet(const LatticeCellListSet& cls) 
        : _list(cls._list) { }
      
      LatticeCellListSet& operator=(const LatticeCellListSet& cls) {
        if(this!=&cls) { this->_list=cls._list; } return *this; }
      
      dimension_type dimension() const { return _list.array_size(); }
      size_type empty() const { return _list.empty(); }
      size_type size() const { return _list.size(); }
      LatticeCell operator[] (size_type i) const { return LatticeCell(_list[i]); }
      
      LatticeRectangle bounds() const;

      /*! \brief Adjoins a LatticeCell to the set. */
      void adjoin(const LatticeCell& c) { 
        assert(this->dimension() == c.dimension());
        this->_list.push_back(c.position()); 
      }
      /*! \brief Adjoins all cells in a LatticeRectangle to the set. */
      void adjoin(const LatticeRectangle& r);
      /*! \brief Adjoins a LatticeCellListSet to the set. */
      void adjoin(const LatticeCellListSet& cl);
      /*! \brief Adjoins a LatticeRectangleListSet to the set. */
      void adjoin(const LatticeRectangleListSet& rl);
      /*! \brief Adjoins a LatticeMaskSet to the set. */
      void adjoin(const LatticeMaskSet& ms);
      
      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const { return const_iterator(_list.begin()); }
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const { return const_iterator(_list.end()); }
      
      /*! \brief Sorts the cells lexicographically, removing duplicates. */
      void unique_sort(); 
     private:
      array_vector<index_type> _list;
    };
      


    /*!\brief A list of rectangles in a lattice. */
    class LatticeRectangleListSet {
     public:
      //      typedef LatticeRectangleListSetIterator iterator;
      //typedef LatticeRectangleListSetIterator const_iterator;
      typedef pair_constructor_iterator<array_vector<index_type>::const_iterator, LatticeRectangle> iterator;
      typedef pair_constructor_iterator<array_vector<index_type>::const_iterator, LatticeRectangle> const_iterator;
     public:
      LatticeRectangleListSet(dimension_type n) 
        : _list(n) { }

      LatticeRectangleListSet(const LatticeRectangleListSet& rls) 
        : _list(rls._list) { }
      
      dimension_type dimension() const { return _list.array_size(); }
      size_type empty() const { return _list.empty(); }
      size_type size() const { return _list.size()/2; }
      LatticeRectangle operator[] (size_type i) const { 
        return LatticeRectangle(_list[2*i],_list[2*i+1]); 
      }
      LatticeRectangle bounds() const;

      /*! \brief Adjoins a LatticeCell to the set. */
      void adjoin(const LatticeCell& c) { 
        assert(this->dimension() == c.dimension());
        this->_list.push_back(c.lower()); 
        this->_list.push_back(c.upper()); 
      }
      /*! \brief Adjoins all cells in a LatticeRectangle to the set. */
      void adjoin(const LatticeRectangle& r) { 
        assert(this->dimension() == r.dimension());
        this->_list.push_back(r.lower()); 
        this->_list.push_back(r.upper()); 
      }
      /*! \brief Adjoins a LatticeCellListSet to the set. */
      void adjoin(const LatticeCellListSet& cl);
      /*! \brief Adjoins a LatticeRectangleListSet to the set. */
      void adjoin(const LatticeRectangleListSet& rl);
      /*! \brief Adjoins a LatticeMaskSet to the set. */
      void adjoin(const LatticeMaskSet& ms);
      
      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const  { return const_iterator(_list.begin()); }
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const { return const_iterator(_list.end()); }
     private:
      array_vector<index_type> _list;
    };
       



    /*!\brief A list of arrays of integers of the same size, representing rectangles in a grid. */
    class LatticeMaskSet {
     public:
      typedef mask_iterator<LatticeRectangle::const_iterator,BooleanArray::const_iterator> iterator;
      typedef iterator const_iterator;
     public:
      LatticeMaskSet(const size_type& n) 
        : _bounds(n), _mask() { this->_compute_cached_attributes(); }
      LatticeMaskSet(const LatticeRectangle& bb) 
        : _bounds(bb), _mask(bb.size(),false) { this->_compute_cached_attributes(); }
      LatticeMaskSet(const LatticeRectangle& bb, const BooleanArray& ma) 
        : _bounds(bb), _mask(ma) { this->_compute_cached_attributes(); }
      LatticeMaskSet(const LatticeRectangle& bb, const LatticeCellListSet& cls);
      LatticeMaskSet(const LatticeRectangle& bb, const LatticeRectangleListSet& rls);

      LatticeMaskSet(const LatticeMaskSet& ms);
      
      const LatticeRectangle& bounds() const { return _bounds; }
      const BooleanArray& mask() const { return _mask; };
      dimension_type dimension() const { return _bounds.dimension(); }
      size_type capacity() const { return _mask.size(); }
      size_type empty() const { return this->size()==0; }
      size_type size() const { return std::count(_mask.begin(),_mask.end(),true); }
      LatticeCell operator[](size_type i) const; 

      const IndexArray& lower() const { return _lower; }
      const IndexArray& upper() const { return _upper; }
      const SizeArray& sizes() const { return _sizes; }
      const SizeArray& strides() const { return _strides; }
      
      /*! Compute the index of a position in a grid. */
      size_type index(const IndexArray& pos) const;

      /*! Compute the position of an index in a grid. */
      IndexArray position(size_type index) const;

      /*! \brief Empties the set. */
      void clear();

      /*! \brief Adjoins a LatticeCell to the set. */
      void adjoin(const LatticeCell& c) { 
        assert(subset(c,this->bounds()));
        this->_mask[this->index(c.position())]=true;
      }
      
      /*! \brief Adjoins a LatticeRectangle to the set. */
      void adjoin(const LatticeRectangle& r);

      /*! \brief Adjoins a GridCellListSet to the set. */
      void adjoin(const LatticeCellListSet& cl);
      /*! \brief Adjoins a GridRectangleListSet to the set. */
      void adjoin(const LatticeRectangleListSet& rl);
      /*! \brief Adjoins a LatticeMaskSet to the set. */
      void adjoin(const LatticeMaskSet& ms);
        
      /*! \brief The one-box neighbourhood on the same grid. */
      LatticeMaskSet neighbourhood() const;
      /*! \brief The adjoining elements on the same grid. */
      LatticeMaskSet adjoining() const;

      /*! \brief Constant iterator to the beginning of the cells in the set. */
      const_iterator begin() const { return const_iterator(_bounds.begin(),_mask.begin(),_mask.end()); }
      /*! \brief Constant iterator to the end of the cells in the set. */
      const_iterator end() const { return const_iterator(_bounds.end(),_mask.end(),_mask.end()); }
     private:
      void _compute_cached_attributes();

     private:
      LatticeRectangle _bounds;
      IndexArray _lower;
      IndexArray _upper;
      SizeArray _sizes;
      SizeArray _strides;
      BooleanArray _mask;
    };

  }
}

#endif /* _ARIADNE_LATTICE_SET_H */
