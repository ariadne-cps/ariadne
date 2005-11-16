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
 *  \brief Infinite sequences.
 */

#ifndef _ARIADNE_SEQUENCE_H
#define _ARIADNE_SEQUENCE_H

#include <cstddef>
#include <cassert>

#include <exception>
#include <stdexcept>
#include <iterator>
#include <iostream>


namespace Ariadne {
  
  template<typename T> class sequence;
  template<typename T> class compact_sequence;
  template<typename T> class _sequence_const_iterator;
  template<typename T> class _compact_sequence_const_iterator;
  
  template<typename T> compact_sequence<T> convolution(const compact_sequence<T>& u, const compact_sequence<T>& v);

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
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    
    ~sequence() { delete[] _ptr; }
    sequence() : _body_size(0), _tail_size(1), _ptr(new value_type[1]) { _ptr[0]=value_type(); }

    template<typename ForwardIterator> sequence(ForwardIterator b, ForwardIterator tb, ForwardIterator te)
      : _body_size(std::distance(b,tb)), _tail_size(std::distance(tb,te)), _ptr(new value_type[std::distance(b,te)])
    { fill(b); }

    sequence(const sequence<T>& a)
      : _body_size(a._body_size), _tail_size(a._tail_size), _ptr(new value_type[_body_size+_tail_size])
    { fill(a._ptr); }

    sequence<T>& operator=(const sequence<T>& a) {
      if(&a!=this) { resize(a.body_size(),a.tail_size()); fill(a._ptr); } return *this; }

    bool operator==(const sequence& a) const {
      /* FIXME: How far do we need to go? */
      for(size_type i=0; i!=a.body_size()+a.tail_size()+body_size()+tail_size(); ++i) {
        if((*this)[i]!=a[i]) { return false; } }
      return true; }

    size_type body_size() const { return _body_size; }
    size_type tail_size() const { return _tail_size; }
    void resize(size_type bs, size_type ts) {
      delete[] _ptr; _body_size=bs; _tail_size=ts;  _ptr=new value_type[_body_size+_tail_size]; }
  
    const_reference operator[](size_type i) const {
      if(i>=_body_size) { i = ( (i-_body_size) % _tail_size ) + _body_size; } return _ptr[i]; }
    value_type get(size_type i) const {
      return _ptr[i]; }
    const_reference at(size_type i) const {
      return _ptr[i]; }
  
    const_iterator begin() const;
     
    /*
      bool operator==(const sequence& other) const {
      if(_body_size!=other._body_size) return false; 
      const_iterator first=begin(); const_iterator last=end(); const_iterator curr=other.begin(); 
      while(first!=last) { if((*first)!=(*curr)) { return false; } ++first; ++curr; } return true; }
      bool operator!=(const array& other) const { return !((*this)==other); }
    */

    template<typename InputIterator> void fill(InputIterator iter) { 
      pointer curr=_ptr; pointer end=curr+size(); while(curr!=end) { *curr=*iter; ++curr; ++iter; } }
    template<typename ForwardIterator> void assign(ForwardIterator first, ForwardIterator tail, ForwardIterator last) {
      resize(std::distance(first,tail),std::distance(tail,last)); fill(first); }
   private:
    size_type size() const { return _body_size + _tail_size; }
    const_pointer pointer_begin() const { return _ptr; }
   private:
    size_type _body_size;
    size_type _tail_size;
    pointer _ptr;
  };

  template<class C>
  class _container_reference {
    typedef typename C::value_type value_type;
    typedef typename C::size_type index_type;
   public:
    _container_reference(C& container, index_type& index) : _container(container), _index(index) { }
    operator value_type() const { return _container.get(_index); }
    void operator=(const value_type& x) { _container.set(_index,x); }
   private:
    C& _container;
    index_type& _index;
  };

  template<class T>
  class _sequence_const_iterator : public std::iterator<std::forward_iterator_tag,T> {
    typedef typename sequence<T>::const_pointer const_pointer;
    typedef typename sequence<T>::const_reference reference;
    typedef typename sequence<T>::const_pointer pointer;
   public:
    _sequence_const_iterator (const sequence<T>& seq, size_t n=0)
      : _current(seq._ptr + n), _tail_begin(seq._ptr+seq._body_size), _tail_end(seq._ptr+seq._body_size+seq._tail_size) { }
    _sequence_const_iterator (const_pointer c, const_pointer tb, const_pointer te)
      : _current(c), _tail_begin(tb), _tail_end(te) { }
    _sequence_const_iterator (const _sequence_const_iterator<T>& iter) 
      : _current(iter._current), _tail_begin(iter._tail_begin), _tail_end(iter._tail_end) { }
    _sequence_const_iterator<T>& operator++() { ++_current; if(_current==_tail_end) { _current=_tail_begin; } return *this; }
    reference operator*() const { return *_current; }
    pointer operator->() const { return _current; }
    bool operator==(const _sequence_const_iterator<T>& other) const { return ( this->_current == other._current); }
   private:
    const_pointer _current;
    const_pointer _tail_begin;
    const_pointer _tail_end;
  };

  template<typename T>
  inline
  typename sequence<T>::const_iterator
  sequence<T>::begin() const
  {
    return typename sequence<T>::const_iterator(_ptr,_ptr+_body_size,_ptr+_body_size+_tail_size);
  }

  /*! \brief An eventually-zero biinfinite sequence.
   */
  template<typename T> class compact_sequence {
    friend class _compact_sequence_const_iterator<T>;
//    friend template<> compact_sequence<T> convolution(const compact_sequence<T>& u, const compact_sequence<T>& v);
   public:
    typedef _compact_sequence_const_iterator<T> iterator;
    typedef _compact_sequence_const_iterator<T> const_iterator;
    
    typedef T value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef ptrdiff_t size_type;
    typedef ptrdiff_t difference_type;
    
    ~compact_sequence() { delete[] (_ptr+_start); }
    compact_sequence() : _start(0), _finish(0), _ptr(new value_type[0]) { }
    compact_sequence(size_type start, size_type finish, value_type x) : _start(start), _finish(finish), _ptr(new value_type[finish-start]) { 
      for(size_type i=_start; i!=_finish; ++i) { _ptr[i]=x; }
    }
    template<typename ForwardIterator> compact_sequence(const ForwardIterator& b, ForwardIterator c, ForwardIterator e)
      : _start(std::distance(c,b)), _finish(std::distance(c,e)), _ptr(new value_type[std::distance(b,c)]-_start)
    { fill(c); }

    compact_sequence(const compact_sequence<T>& s)
      : _start(s._start), _finish(s._finish), _ptr(new value_type[-_start+_finish]+_start)
    { fill(s._ptr); }

    sequence<T>& operator=(const sequence<T>& s) {
      if(&s!=this) { resize(s._start,s._finish); fill(s._ptr); } return *this; }

    bool operator==(const sequence<T>& s) const {
      /* FIXME: How far do we need to go? */
      for(size_type i=_start; i!=_finish; ++i) {
        if((*this)[i]!=s[i]) { return false; } }
      return true; }

    size_type start() const { return _start; }
    size_type finish() const { return _finish; }

    void resize(size_type s, size_type f) {
      pointer _old_pointer=_ptr;
      size_type _old_start=_start;
      size_type _old_finish=_finish;
      _start=s; 
      _finish=f;  
      _ptr=new value_type[-_start+_finish]-_start;
      for(size_type i=_start; i!=_finish; ++i) {
        _ptr[i]=0;
      }
      for(size_type i=_old_start; i!=_old_finish; ++i) {
        _ptr[i]=_old_pointer[i];
      }
      delete[] (_old_pointer+_old_start); 
    }

    const_reference get(difference_type i) const {
      if(i<_start || i>=_finish) { return 0; } return _ptr[i]; 
    }
      
    void set(difference_type i, reference x) {
      if(i<_start) {
        resize(i,_finish);
      } 
      else if(i>=_finish) {
        resize(_start,i+1);
      }
      _ptr[i]=x;
    }

    const_iterator begin() const;
     
    template<typename InputIterator> void fill(InputIterator iter) { 
      pointer curr=_ptr+_start; iter=iter+_start; pointer end=curr+_finish; while(curr!=end) { *curr=*iter; ++curr; ++iter; } }
    template<typename ForwardIterator> void assign(ForwardIterator first, ForwardIterator center, ForwardIterator last) {
      resize(std::distance(center,first),std::distance(center,last)); fill(center); }
   private:
    size_type size() const { return _finish-_start; }
    const_pointer pointer_begin() const { return _ptr; }
   private:
    size_type _start;
    size_type _finish;
    pointer _ptr;
  };

  template<typename T>
  compact_sequence<T>
  convolution(const compact_sequence<T>& u, const compact_sequence<T>& v) {
    typedef typename compact_sequence<T>::size_type size_type;
    compact_sequence<T> r(u.start()-v.finish(),u.finish()-v.start(),0);
    for(size_type j=u._start; j!=u.finish; ++j) {
      for(size_type k=v._start; k!=v.finish; ++k) {
       r._ptr[j-k] += u._ptr[j] * v._ptr[k];
      }
    }
    return r;
  }

  template<typename T>
  std::ostream&
  operator<<(std::ostream& os, const sequence<T>& s)
  {
    typedef typename sequence<T>::size_type size_type;
    os << "[";
    for(size_type i=0; i!=s.body_size(); ++i) {
      if(i!=0) { os << ","; }
      os << s[i];
    }
    os << ";" << s[s.body_size()];
    for(size_type i=s.body_size()+1; i!=s.body_size()+s.tail_size(); ++i) {
      os << "," << s[i];
    }
    os << ",...]";
    return os;
  }

} // namespace Ariadne

#endif /* _ARIADNE_SEQUENCE_H */
