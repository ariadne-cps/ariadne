/***************************************************************************
 *            iterator.h
 *
 *  Copyright  2006  Pieter Collins, Alberto Casagrande
 *  Email: Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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

/*! \file iterator.h
 *  \brief General-purpose iterators
 */

#ifndef ARIADNE_ITERATOR_H
#define ARIADNE_ITERATOR_H

#include <boost/iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include "base/array.h"
#include "base/types.h"

namespace Ariadne {
  namespace Base {

    //! \ingroup Traversal
    /*! \brief An iterator traversing a subset of values defined by a mask */
    template<class Base, class Mask>
    class mask_iterator
      : public boost::iterator_facade<mask_iterator<Base,Mask>,
                                      typename Base::value_type,
                                      boost::forward_traversal_tag,
                                      typename Base::reference,
                                      typename Base::pointer>
    {
     public:
      mask_iterator(Base i, Mask mi, Mask me)
        : _base(i), _mask_iter(mi), _mask_end(me) 
      { 
        while(_mask_iter!=_mask_end && !*_mask_iter) { 
          ++_mask_iter; ++_base; 
        }
      }
     private:
      bool equal(const mask_iterator& other) const {
        return this->_base==other._base && this->_mask_iter==other._mask_iter;
      }
      typename Base::reference dereference() const { 
        return _base.operator*(); 
      }
      void increment() { 
        do { 
          ++_mask_iter; ++_base; } 
        while(_mask_iter!=_mask_end && !*_mask_iter);
      }
     private:
      friend class boost::iterator_core_access;
      Base _base;
      Mask _mask_iter;
      Mask _mask_end;
    };
    
    
    
    //! \ingroup Traversal
    /*! \brief An iterator adaptor which uses a conversion operator to convert its argument into another type. */
    template<class Base, class Value>
    class conversion_iterator
      : public boost::iterator_adaptor<conversion_iterator<Base,Value>,
                                       Base,
                                       Value,
                                       boost::use_default,
                                       Value>
    {
     public:
      conversion_iterator(Base i)
        : conversion_iterator::iterator_adaptor_(i) { }
     private:
      friend class boost::iterator_core_access;
      Value dereference() const { 
        return Value(*this->base_reference());
      }
    };
    
    
    
    //! \ingroup Traversal
    /*! \brief An iterator which constructs its results from a piece of data and an iterator. */
    template<class Base, class Value, class Data>
    class binary_constructor_iterator 
      : public boost::iterator_adaptor<binary_constructor_iterator<Base,Value,Data>,
                                       Base,
                                       Value,
                                       boost::use_default,
                                       Value>
    {
     public:
      binary_constructor_iterator(const Data& d, Base i)
        : binary_constructor_iterator::iterator_adaptor_(i),
          _data(d) { }
     private:
      friend class boost::iterator_core_access;
      Value dereference() const { return Value(_data,*this->base()); }
      Data _data;
    };
    
    
    
    //! \ingroup Traversal
    /*! \brief An iterator which constructs its results from a piece of data and an iterator. */
    template<class Base, class Value, class Data1, class Data2>
    class ternary_constructor_iterator 
      : public boost::iterator_adaptor<ternary_constructor_iterator<Base, Value,
                                                                    Data1,Data2>,
                                       Base,
                                       Value,
                                       boost::use_default,
                                       Value>
    {
     public:
      ternary_constructor_iterator(const Data1& d1, const Data2& d2, Base i)
        : ternary_constructor_iterator::iterator_adaptor_(i),
          _data1(d1), _data2(d2) { }
     private:
      friend class boost::iterator_core_access;
      Value dereference() const { return Value(_data1,_data2,*this->base()); }
      Data1 _data1;
      Data2 _data2;
    };
    
    
    
    //! \ingroup Traversal
    /*! \brief An iterator which constructs its results from a pair or adjacent values of another iterator. */
    template<class Base, class Value>
    class pair_constructor_iterator
      : public boost::iterator_adaptor<pair_constructor_iterator<Base,Value>,
                                       Base,
                                       Value,
                                       boost::use_default,
                                       Value>
    {
     public:
      pair_constructor_iterator(Base i)
        : pair_constructor_iterator::iterator_adaptor_(i) { }
     private:
      void increment() { 
        ++(++this->base_reference()); }
      void advance(typename Base::difference_type n) { 
        this->base_reference()+=(2*n); }
      Value dereference() const { 
        Base tmp=this->base(); ++tmp; return Value(*this->base_reference(),*tmp); }
      friend class boost::iterator_core_access;
    };
    
    
    
    //! \ingroup Traversal
    /*! \brief An iterator for positions in rectangular piece of a grid. */
    class lattice_iterator 
      : public boost::iterator_facade<lattice_iterator,
                                      array<index_type>,
                                      boost::forward_traversal_tag,
                                      const array<index_type>&>
    {
     public:
      lattice_iterator(const array<index_type>& p, 
                       const array<index_type>& l, 
                       const array<index_type>& u)
        : _position(p), _lower(l), _upper(u) { }
     private:
      bool equal(const lattice_iterator& other) const {
        return (this->_position==other._position) 
          && (this->_lower==other._lower) && (this->_upper==other._upper);
      }
      void increment();
      const array<index_type>& dereference() const { return _position; }
      friend class boost::iterator_core_access;
     private:
      array<index_type> _position;
      array<index_type> _lower;
      array<index_type> _upper;
    };
    
    
    inline
    void 
    lattice_iterator::increment()
    {
      dimension_type d=0;
      _position[d]+=1;
      while(_position[d]==_upper[d] && (d+1u)!=_position.size() ) {
        _position[d]=_lower[d];
        d+=1;
        _position[d]+=1;
      }
    }
    
    
  } // namespace Base
    
} // namespace Ariadne
  
  
#endif /* ARIADNE_ITERATOR_H */
