/***************************************************************************
 *            array.h
 *
 *  4 October 2004
 *  Copyright  2004  Alberto Casagrande, Pieter Collins
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
 
/*! \file array.h
 *  \brief STL style arrays.
 */

#ifndef _ARRAY_H
#define _ARRAY_H

#include <cstddef>
#include <exception>

namespace Ariadne {

template<typename T, size_t N=0> class array;

/*! \brief STL style interface to dynaically-sized arrays. 
 *
 * An array<T> is a variable-size array which can be resized and is allocated 
 * on the heap. Arrays provide checked access using at and unchecked access using operator[].
 *
 * FIXME: Exception safety in constructors!
 */
template<typename T> class array<T> {
public:
  typedef T* iterator;
  typedef const T* const_iterator;
  
  typedef T value_type;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef pointer iterator;
  typedef const_pointer const_iterator;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

  ~array() { delete[] _ptr; }
  array() : _size(0), _ptr(new value_type[_size]) { }
  explicit array(const size_type __n) : _size(__n), _ptr(new value_type[_size]) { }
  array(const size_type __n, const value_type& __val) : _size(__n), _ptr(new value_type[_size]) { fill(__val); }
  template<class In> array(In __first, In __last) : _size(std::distance(__first,__last)), _ptr(new value_type[_size]) { 
    fill(__first); }
	
  array(const array& a) : _size(a.size()), _ptr(new value_type[_size]) { fill(a.begin()); }
  array& operator=(const array& a) { resize(a.size()); fill(a.begin()); return *this; }
  template<size_t N> array(const array<T,N>& a); /* inline conversion */ 
  template<size_t N> array& operator=(const array<T,N>& a); /* inline conversion */ 

  size_type empty() const { return _size==0u; }
  size_type size() const { return _size; }
  size_type max_size() const { return (size_type) (-1); }
  void resize(size_type __n) { if(size()!=__n) { delete[] _ptr; _size=__n; _ptr=new value_type[_size]; } }

  reference operator[](size_type i) { return _ptr[i]; } 
  const_reference operator[](size_type i) const { return _ptr[i]; }
  reference at(size_type i) { if(i<_size) { return _ptr[i]; } else { throw std::out_of_range(); } } 
  const_reference at(size_type i) const { if(i<_size) { return _ptr[i]; } else { throw std::out_of_range(); } }

  iterator begin() { return iterator(_ptr); }
  const_iterator begin() const { return const_iterator(_ptr); }
  iterator end() { return iterator(_ptr+_size); }
  const_iterator end() const { return const_iterator(_ptr+_size); }

  bool operator==(const array& other) const {
    if(size()!=other.size()) return false; 
    const_iterator __first=begin(); const_iterator __last=end(); const_iterator __other=other.begin(); 
    while(__first!=__last) { if((*__first)!=(*__other)) { return false; } ++__first; ++__other; } return true; }
  bool operator!=(const array& other) const { return !((*this)==other); }

  void fill(value_type __val) { 
    pointer __curr=_ptr; pointer __end=__curr+size(); while(__curr!=__end) { *__curr=__val; ++__curr; } }
  template<class Inputiterator> void fill(Inputiterator __iter) { 
    pointer __curr=_ptr; pointer __end=__curr+size(); while(__curr!=__end) { *__curr=*__iter; ++__curr; ++__iter; } }
  template<class ForwardIterator> void assign(ForwardIterator __first, ForwardIterator __last) { 
    resize(std::distance(__first,__last)); fill(__first); }
 private:
  size_type _size; 
  pointer _ptr; 
};


/*! \brief STL style interface to fixed-size arrays. 
 *
 * An array<T,N> is a fixed-size array of size \a N which is allocated 
 * on the stack. Arrays provide checked access using at and unchecked access using operator[].
 *
 * FIXME: Constructors for small arrays.
 */
template<typename T, size_t N> class array {
public:
  typedef T value_type;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef pointer iterator;
  typedef const_pointer const_iterator;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;

  ~array() { }
  array() { }
  explicit array(const value_type& __val) { fill(__val); }
  template<class In> array(In __first, In __last) { assert(distance(__first,__last==N)); fill(__first); }
  array(const array& a) { fill(a.begin()); }
  array& operator=(const array& a) { fill(a.begin()); return *this; }
	
  array(const value_type& __x, const value_type& __y) { 
    assert(N==2); _ptr[0]=__x; _ptr[1]=__y; }
  array(const value_type& __x, const value_type& __y, const value_type& __z) { 
    assert(N==3); _ptr[0]=__x; _ptr[1]=__y; _ptr[2]=__z; }
  array(const value_type& __w, const value_type& __x, const value_type& __y, const value_type& __z) { 
    assert(N==4); _ptr[0]=__w; _ptr[1]=__x; _ptr[2]=__y; _ptr[3]=__z; }
  
  size_type empty() const { return N==0u; }
  size_type size() const { return N; }
  size_type max_size() const { return N; }

  reference operator[](size_type i) { return _ptr[i]; }
  const_reference operator[](size_type i) const { return _ptr[i]; }
  reference at(size_type i) { assert(i<N); return _ptr[i]; }
  const_reference at(size_type i) const { assert(i<N); return _ptr[i]; }
  
  iterator begin() { return _ptr; }
  iterator end() { return _ptr+N; }
  const_iterator begin() const { return _ptr; }
  const_iterator end() const { return _ptr+N; }
  
  bool operator==(const array& other) const {
    const_iterator __first=begin(); const_iterator __last=end(); const_iterator __other=other.begin(); 
    while(__first!=__last) { if((*__first)!=(*__other)) { return false; } ++__first; ++__other; } return true; }
  bool operator!=(const array& other) const { return !((*this)==other); }
  void fill(value_type __val) { 
    for(size_type i=0; i!=N; ++i) { _ptr[i]=__val; } }
  template<class In> void fill(In __iter) { 
    _assign_iter(__first); }
  template<class In> void assign(In __first, In __last) { 
    assert(distance(__first,__last)==size()); fill(__first); }
 private:
  template<class Iter> void _assign_iter(Iter __iter); 
 private:
  value_type _ptr[N];
};
	
	
template<typename T, size_t N> template<class Iter> 
inline 
void array<T,N>::_assign_iter(Iter __iter)
{ 
  if(N==1) { 
    *_ptr=*__iter; }
  else if(N==2) { 
    value_type* __curr=_ptr; *__curr=*__iter; ++__curr; ++__iter;  *__curr=*__iter; }
  else if(N==3) { 
    value_type* __curr=_ptr; *__curr=*__iter; 
    ++__curr; ++__iter;  *__curr=*__iter; 
    ++__curr; ++__iter;  *__curr=*__iter; }
  else if(N==4) { 
    value_type* __curr=_ptr; *__curr=*__iter; 
    ++__curr; ++__iter; *__curr=*__iter; 
    ++__curr; ++__iter;  *__curr=*__iter;  
    ++__curr; ++__iter;  *__curr=*__iter; }
  else {
    value_type* __curr=_ptr; value_type* __end=_ptr+N; *__curr=*__iter;
    while(__curr!=__end) { ++__curr; ++__iter;  *__curr=*__iter;  } }
}
 
template<typename T> template<size_t N> 
inline
array<T,0>::array(const array<T,N>& a) : _size(a.size()), _ptr(new value_type[_size])
{ 
  fill(a.begin()); 
}

template<typename T> template<size_t N> 
inline
array<T,0>& array<T,0>::operator=(const array<T,N>& a) 
{ 
  resize(a.size()); fill(a.begin()); return *this; 
}

} // namespace Ariadne

#endif
