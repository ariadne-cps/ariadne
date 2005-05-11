/***************************************************************************
 *            sequence.h
 *
 *  9 May 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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
 
/*! \file sequence.h
 *  \brief An eventually-periodic sequence.
 */

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <cstddef>
#include <cassert>

#include <exception>
#include <stdexcept>
#include <iterator>


namespace Ariadne {
  
  template<typename T> class sequence;
  template<typename T> class _sequence_const_iterator;

  /*! \brief An eventually-periodic sequence with an STL style interface.
   *
   * A sequence<T> is an eventually-periodic sequence of values with STL style forward iterators.
   */
  template<typename T> class sequence {
    friend class _sequence_const_iterator<T>;
   public:
    typedef _sequence_const_iterator<T> iterator;
    typedef _sequence_const_iterator<T> const_iterator;
    
    typedef T value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef pointer iterator;
    typedef const_pointer const_iterator;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    
    ~sequence() { delete[] _ptr; }
    sequence() : _body_size(0), _tail_size(1), _ptr(new value_type[1]) { _ptr[0]=value_type(); }
    template<typename ForwardIterator> sequence(ForwardIterator& b, ForwardIterator tb, ForwardIterator te) 
      : _body_size(std::difference(b,tb)), _tail_size(std::difference(tb,te)), _ptr(new value_type[std::difference(b,te)]) 
    { fill(b); }
    
    sequence(const sequence& a) : _body_size(a._body_size), _tail_size(a._tail_size), _ptr(new value_type[_body_size+_tail_size]) { fill(a.begin()); }
    sequence& operator=(const sequence& a) { resize(a.body_size(),a.tail_size()); fill(a.begin()); return *this; }
    
    size_type body_size() const { _body_size; }
    size_type tail_size() const { _tail_size; }
    void resize(size_type bs, size_type ts) { delete[] _ptr; _body_size=bs; _tail_size=ts;  _ptr=new value_type[_body_size+_tail_size]; } 
  
    const_reference operator[](size_type i) const {  if(i>=_body_size) { i = ( (i-_body_size) % _tail_size ) + _body_size; } return _ptr[i]; }
    const_reference at(size_type i) const { return _ptr[i]; }
  
    const_iterator begin() const;
     
    /*
      bool operator==(const sequence& other) const {
      if(_body_size!=other._body_size) return false; 
      const_iterator first=begin(); const_iterator last=end(); const_iterator curr=other.begin(); 
      while(first!=last) { if((*first)!=(*curr)) { return false; } ++first; ++curr; } return true; }
      bool operator!=(const array& other) const { return !((*this)==other); }
    */
    
    template<typename InputIterator> void fill(InputIterator iter) { 
      pointer curr=_ptr; pointer end=curr+_size(); while(curr!=end) { *curr=*iter; ++curr; ++iter; } }
    template<typename ForwardIterator> void assign(ForwardIterator first, ForwardIterator tail, ForwardIterator last) { 
      resize(std::distance(first,tail),std::distance(tail,last)); fill(first); }
   private:
    size_type _size() const { return _body_size + _tail_size; }
   private:
    size_type _body_size; 
    size_type _tail_size; 
    pointer _ptr; 
  };

  template<class T>
  class _sequence_const_iterator : public std::iterator<std::forward_iterator_tag,T> {
    typedef typename sequence<T>::const_pointer const_pointer;
    typedef typename sequence<T>::const_reference reference;
    typedef typename sequence<T>::const_pointer pointer;
   public:
    _sequence_const_iterator (const sequence<T>& seq, size_t n=0)
      : _tail_begin(seq.pointer()+seq.body_size()), _tail_end(seq.pointer() + seq.body_size()),
	_current(seq.pointer() + n) { }
    _sequence_const_iterator (const_pointer tb, const_pointer te, const_pointer c) 
      : _tail_begin(tb), _tail_end(te), _current(c) { }
    _sequence_const_iterator (const _sequence_const_iterator<T>& iter) 
      : _tail_begin(iter._tail_begin), _tail_end(iter._tail_end), _current(iter._current) { }
    _sequence_const_iterator<T>& operator++() { ++_current; if(_current==_tail_end) { _current=_tail_begin; } return *this; }
    reference operator*() const { return *_current; }
    pointer operator->() const { return _current; }
    bool operator==(const _sequence_const_iterator<T>& other) const { return ( this->_current == other._current); }
   private:
    const_pointer _tail_begin;
    const_pointer _tail_end;
    const_pointer _current;
  };

  template<typename T>
  inline
  typename sequence<T>::const_iterator
  sequence<T>::begin() const
  {
    return sequence<T>::const_iterator(_ptr,_ptr+_body_size,_ptr+_body_size+_tail_size);
  }
 
} // namespace Ariadne

#endif
