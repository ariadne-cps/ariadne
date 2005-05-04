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
#include <stdexcept>
#include <cassert>

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
  explicit array(const size_type n) : _size(n), _ptr(new value_type[_size]) { }
  array(const size_type n, const value_type& val) : _size(n), _ptr(new value_type[_size]) { fill(val); }
  template<class In> array(In first, In last) : _size(std::distance(first,last)), _ptr(new value_type[_size]) { 
    fill(first); }
	
  array(const array& a) : _size(a.size()), _ptr(new value_type[_size]) { fill(a.begin()); }
  array& operator=(const array& a) { resize(a.size()); fill(a.begin()); return *this; }
  template<size_t N> array(const array<T,N>& a); /* inline conversion */ 
  template<size_t N> array& operator=(const array<T,N>& a); /* inline conversion */ 

  size_type empty() const { return _size==0u; }
  size_type size() const { return _size; }
  size_type max_size() const { return (size_type) (-1); }
  void resize(size_type n) { if(size()!=n) { delete[] _ptr; _size=n; _ptr=new value_type[_size]; } }

  reference operator[](size_type i) { return _ptr[i]; } 
  const_reference operator[](size_type i) const { return _ptr[i]; }
  reference at(size_type i) { if(i<_size) { return _ptr[i]; } else { throw std::out_of_range("array index out-of-range"); } } 
  const_reference at(size_type i) const { if(i<_size) { return _ptr[i]; } else { throw std::out_of_range("array index out-of-range"); } }

  iterator begin() { return iterator(_ptr); }
  const_iterator begin() const { return const_iterator(_ptr); }
  iterator end() { return iterator(_ptr+_size); }
  const_iterator end() const { return const_iterator(_ptr+_size); }

  bool operator==(const array& other) const {
    if(size()!=other.size()) return false; 
    const_iterator first=begin(); const_iterator last=end(); const_iterator curr=other.begin(); 
    while(first!=last) { if((*first)!=(*curr)) { return false; } ++first; ++curr; } return true; }
  bool operator!=(const array& other) const { return !((*this)==other); }

  void fill(value_type val) { 
    pointer curr=_ptr; pointer end=curr+size(); while(curr!=end) { *curr=val; ++curr; } }
  template<class Inputiterator> void fill(Inputiterator iter) { 
    pointer curr=_ptr; pointer end=curr+size(); while(curr!=end) { *curr=*iter; ++curr; ++iter; } }
  template<class ForwardIterator> void assign(ForwardIterator first, ForwardIterator last) { 
    resize(std::distance(first,last)); fill(first); }
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
  explicit array(const value_type& val) { fill(val); }
  template<class In> array(In first, In last) { assert(distance(first,last==N)); fill(first); }
  array(const array& a) { fill(a.begin()); }
  array& operator=(const array& a) { fill(a.begin()); return *this; }
	
  array(const value_type& x, const value_type& y) { 
    assert(N==2); _ptr[0]=x; _ptr[1]=y; }
  array(const value_type& x, const value_type& y, const value_type& z) { 
    assert(N==3); _ptr[0]=x; _ptr[1]=y; _ptr[2]=z; }
  array(const value_type& w, const value_type& x, const value_type& y, const value_type& z) { 
    assert(N==4); _ptr[0]=w; _ptr[1]=x; _ptr[2]=y; _ptr[3]=z; }
  
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
    const_iterator first=begin(); const_iterator last=end(); const_iterator curr=other.begin(); 
    while(first!=last) { if((*first)!=(*curr)) { return false; } ++first; ++curr; } return true; }
  bool operator!=(const array& other) const { return !((*this)==other); }
  void fill(value_type val) { 
    for(size_type i=0; i!=N; ++i) { _ptr[i]=val; } }
  template<class In> void fill(In iter) { 
    _assign_iter(iter); }
  template<class In> void assign(In first, In last) { 
    assert(distance(first,last)==size()); fill(first); }
 private:
  template<class Iter> void _assign_iter(Iter iter); 
 private:
  value_type _ptr[N];
};
	
	
template<typename T, size_t N> template<class Iter> 
inline 
void array<T,N>::_assign_iter(Iter iter)
{ 
  if(N==1) { 
    *_ptr=*iter; }
  else if(N==2) { 
    value_type* curr=_ptr; *curr=*iter; ++curr; ++iter;  *curr=*iter; }
  else if(N==3) { 
    value_type* curr=_ptr; *curr=*iter; 
    ++curr; ++iter;  *curr=*iter; 
    ++curr; ++iter;  *curr=*iter; }
  else if(N==4) { 
    value_type* curr=_ptr; *curr=*iter; 
    ++curr; ++iter; *curr=*iter; 
    ++curr; ++iter;  *curr=*iter;  
    ++curr; ++iter;  *curr=*iter; }
  else {
    value_type* curr=_ptr; value_type* end=_ptr+N; *curr=*iter;
    while(curr!=end) { ++curr; ++iter;  *curr=*iter;  } }
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
