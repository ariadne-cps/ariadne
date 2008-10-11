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
 *  the Free Software Foundation; either version 2 of the License, or5
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

#ifndef ARIADNE_ARRAY_H
#define ARIADNE_ARRAY_H

#include <cstddef>

#include <exception>
#include <stdexcept>
#include <vector>
#include <iterator>

#include <iostream>

#include "array.decl.h"

namespace Ariadne {
  
  namespace Base {
    
    /*! \ingroup Storage
     *  \brief STL style interface to dynamically-sizable arrays. 
     *
     * An array<T> is a variable-size array which can be resized and is allocated 
     * on the heap. Arrays provide checked access using at and unchecked access using operator[].
     */
    template<class T> class array<T> {
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
      
      /*! \brief Destructor */
      ~array() { delete[] _ptr; }
      
      /*! \brief Default constructor. Constructs an empty array. */
      array() : _size(0), _ptr(new value_type[_size]) { }
      /*! \brief Constructs an array of size \a n with uninitialised elements. */
      explicit array(const size_type n) : _size(n), _ptr(new value_type[_size]) { }
      /*! \brief Constructs an array of size \a n with elements initialised to \a x. */
      array(const size_type n, const value_type& x) : _size(n), _ptr(new value_type[_size]) { fill(x); }
      /*! \brief Constructs an array from the range \a first to \a last. */
      template<class ForwardIterator> array(ForwardIterator first, ForwardIterator last)
        : _size(std::distance(first,last)), _ptr(new value_type[_size]) { fill(first); }
      
      /*! \brief Conversion constructor. */
      template<class T1> array(const array<T1>& a) : _size(a.size()), _ptr(new value_type[_size]) { fill(a.begin()); }
      /*! \brief Copy constructor. */
      array(const array<T>& a) : _size(a.size()), _ptr(new value_type[_size]) { fill(a.begin()); }
      /*! \brief Copy assignment. */ 
      array<T>& operator=(const array<T>& a) { resize(a.size()); fill(a.begin()); return *this; }
      /*! \brief Conversion from an array of fixed size \p N. */ 
      template<unsigned short int N> array(const array<T,N>& a); /* inline conversion */ 
      /*! \brief Assignment from an array of fixed size \p N. */ 
      template<unsigned short int N> array<T>& operator=(const array<T,N>& a); /* inline conversion */ 
      
      /*! \brief True if the array's size is 0. */
      bool empty() const { return _size==0u; }
      /*! \brief The size of the array. */
      size_type size() const { return _size; }
      /*! \brief The maximum possible size of the array. */
      size_type max_size() const { return (size_type) (-1); }
      /*! \brief Resizes the array to hold \a n elements. If \a n is larger than the current size, the extra elements are initialised. */
      void resize(size_type n) { if(size()!=n) { pointer _new_ptr=new value_type[n]; 
        for(size_type i=0; i!=n; ++i) { if(i<_size) { _new_ptr[i]=_ptr[i]; } else { _new_ptr[i]=value_type(); } }
        delete[] _ptr; _size=n; _ptr=_new_ptr; } }
      /*! \brief Reallocates the array to hold \a n elements. The new elements are uninitialised. */
      void reallocate(size_type n) { if(size()!=n) { delete[] _ptr; _size=n; _ptr=new value_type[_size]; } }
      /*! \brief Efficiently swap two arrays.  */
      void swap(array<T>& a) { std::swap(_size,a._size); std::swap(_ptr,a._ptr); }
      
      /*! \brief The \a n th element. */
      reference operator[](size_type i) { return _ptr[i]; } 
      /*! \brief The \a n th element. */
      const_reference operator[](size_type i) const { return _ptr[i]; }
      /*! \brief Checked access to the \a n th element. */
      reference at(size_type i) { if(i<_size) { return _ptr[i]; } else { throw std::out_of_range("array: index out-of-range"); } } 
      /*! \brief Checked access to the \a n th element. */
      const_reference at(size_type i) const { if(i<_size) { return _ptr[i]; } else { throw std::out_of_range("array: index out-of-range"); } }
      
      /*! \brief A reference to the first element of the array. */
      reference front() { return _ptr[0]; }
      /*! \brief A constant reference to the first element of the array. */
      const_reference front() const { return _ptr[0]; }
      /*! \brief A reference to the last element of the array. */
      reference back() { return _ptr[_size-1]; }
      /*! \brief A constantreference  to the last element of the array. */
      const_reference back() const { return _ptr[_size-1]; }
      
      /*! \brief An iterator pointing to the beginning of the array. */
      iterator begin() { return _ptr; }
      /*! \brief A constant iterator pointing to the beginning of the array. */
      const_iterator begin() const { return _ptr; }
      /*! \brief An iterator pointing to the end of the array. */
      iterator end() { return _ptr+_size; }
      /*! \brief A constant iterator pointing to the end of the array. */
      const_iterator end() const { return _ptr+_size; }
      
      /*! \brief Tests two arrays for equality */
      bool operator==(const array& other) const {
        if(size()!=other.size()) return false; 
        const_iterator first=begin(); const_iterator last=end(); const_iterator curr=other.begin(); 
        while(first!=last) { if((*first)!=(*curr)) { return false; } ++first; ++curr; } return true; }
      /*! \brief Tests two arrays for inequality */
      bool operator!=(const array& other) const { return !((*this)==other); }
      
      /*! \brief Fills the array with copies of \a x. */
      void fill(const value_type& x) { 
        iterator curr=_ptr; iterator last=end(); while(curr!=last) { *curr=x; ++curr; } }
      /*! \brief Fills the array from the sequence starting at \a first. */
      template<class InputIterator> void fill(InputIterator first) { 
        iterator curr=begin(); iterator this_end=end(); while(curr!=this_end) { *curr=*first; ++curr; ++first; } }
      /*! \brief Assigns the sequence from \a first to \a last. */
      template<class ForwardIterator> void assign(ForwardIterator first, ForwardIterator last) { 
        resize(std::distance(first,last)); fill(first); }
     private:
      size_type _size; 
      pointer _ptr; 
    };
    
    
    
    
  /*! \brief Specialisation of array template to hold bit values. */
    template<> class array<bool> {
     public:
      typedef std::vector<bool>::value_type value_type;
      typedef std::vector<bool>::reference reference;
      typedef std::vector<bool>::const_reference const_reference;
      typedef std::vector<bool>::pointer pointer;
      typedef std::vector<bool>::const_pointer const_pointer;
      typedef std::vector<bool>::iterator iterator;
      typedef std::vector<bool>::const_iterator const_iterator;
      typedef std::vector<bool>::size_type size_type;
      typedef std::vector<bool>::difference_type difference_type;
      
      ~array() {  }
      array() : _vector() { }
      explicit array(const size_type n) : _vector(n) { }
      array(const size_type n, const value_type& val) : _vector(n,val) { }
      template<class In> array(In first, In last) : _vector(first,last) { }
      
      array(const array& a) : _vector(a._vector) { }
      array(const std::vector<bool>& v) : _vector(v) { }
      array& operator=(const array& a) { if(this!=&a) { _vector=a._vector; } return *this; }
      
      size_type empty() const { return _vector.empty(); }
      size_type size() const { return  _vector.size(); }
      size_type max_size() const { return  _vector.max_size(); }
      void resize(size_type n) { _vector.resize(n); }
      void reallocate(size_type n) { _vector.resize(n); }
      void swap(array<bool>& a) { _vector.swap(a._vector); }
 
      reference operator[](size_type i) { return _vector[i]; }
      const_reference operator[](size_type i) const { return _vector[i]; }
      reference at(size_type i) { if(i<size()) { return _vector[i]; } else { throw std::out_of_range("array::index out-of-range"); } } 
      const_reference at(size_type i) const { if(i<size()) { return _vector[i]; } else { throw std::out_of_range("array::index out-of-range"); } }
      
      reference front() { return _vector.front(); }
      const_reference front() const { return _vector.front(); }
      reference back() { return _vector.back(); }
      const_reference back() const { return _vector.back(); }
      
      iterator begin() { return _vector.begin(); }
      const_iterator begin() const { return _vector.begin(); }
      iterator end() { return _vector.end(); }
      const_iterator end() const { return _vector.end(); }
      
      bool operator==(const array& other) const { return _vector==other._vector; }
      bool operator!=(const array& other) const { return !((*this)==other); }
      
      void fill(value_type val) { 
        iterator curr=begin(); while(curr!=end()) { *curr=val; ++curr; } }
      template<class Inputiterator> void fill(Inputiterator iter) { 
        iterator curr=begin(); while(curr!=end()) { *curr=*iter; ++curr; ++iter; } }
      template<class ForwardIterator> void assign(ForwardIterator first, ForwardIterator last) { 
        resize(std::distance(first,last)); fill(first); }
     private:
      std::vector<bool> _vector;
    };
    
    
    
    
    
    //! \ingroup Storage
    /*! \brief STL style interface to sized arrays. 
     *
     * An array<T,N> is a fixed-size array of size \a N which is allocated 
     * on the stack. Arrays provide checked access using at and unchecked access using operator[].
     */
    // FIXME: Exception safety in constructors
    template<class T, unsigned short int N> class array {
     public:
      typedef T value_type;
      typedef value_type& reference;
      typedef const value_type& const_reference;
      typedef value_type* pointer;
      typedef const value_type* const_pointer;
      typedef pointer iterator;
      typedef const_pointer const_iterator;
      typedef unsigned short int size_type;
      typedef ptrdiff_t difference_type;
      
      /*! \brief Destructor */
      ~array() { }
      /*! \brief Default constructor. Constructs an empty array. */
      array() : _size(0) { }
      /*! \brief Constructs an array of size \a n with uninitialised elements. \a n must be less than or equal to \a N. */
      explicit array(const unsigned short int n) : _size(n) { 
        if(n>N) { throw std::length_error("array<T,N>::array(size_type)"); } }
      /*! \brief Constructs an array of size \a n with elements initialised to \a x. \a n must be less than or equal to \a N. */
      array(const unsigned short int n, const value_type& val) : _size(n) { 
        if(n>N) { throw std::length_error("array<T,N>::array(size_type)"); } fill(val); }
      /*! \brief Constructs an array from the range \a first to \a last. */
      template<class In> array(In first, In last) : _size(distance(first,last)) { 
        if(_size>N) { _size=N; throw std::length_error("array<T,N>::array(size_type)"); } _fill_iter(first); }
      /*! \brief Copy constructor. */
      array(const array<T,N>& a) : _size(a.size()) { _fill_iter(a.begin()); }
      /*! \brief Copy assignment. */ 
      array<T,N>& operator=(const array<T,N>& a) { _size=a.size(); _fill_iter(a.begin()); return *this; }
      
      /*! \brief Construct an array of size 2. */ 
      array(const value_type& x, const value_type& y) : _size(2) { 
        if(N<2) { throw std::length_error("array<T,N>::array(value_type,value_type)"); }
        _ptr[0]=x; _ptr[1]=y; }
      /*! \brief Construct an array of size 3. */ 
      array(const value_type& x, const value_type& y, const value_type& z) : _size(3) { 
        if(N<3) { throw std::length_error("array<T,N>::array(value_type,value_type,value_type)"); }
        _ptr[0]=x; _ptr[1]=y; _ptr[2]=z; }
      /*! \brief Construct an array of size 4. */ 
      array(const value_type& w, const value_type& x, const value_type& y, const value_type& z) : _size(4) { 
        if(N<4) { throw std::length_error("array<T,N>::array(value_type,value_type,value_type,value_type)"); }
        _ptr[0]=w; _ptr[1]=x; _ptr[2]=y; _ptr[3]=z; }
      
      /*! \brief True if the array's size is 0. */
      size_type empty() const { return _size==0u; }
      /*! \brief The size of the array. */
      size_type size() const { return _size; }
      /*! \brief The maximum possible size of the array. */
      size_type max_size() const { return N; }
      /*! \brief Resizes the array to hold \a n elements. \a n must be less than or equal to \a N. If \a n is larger than the current size, the extra elements are initialised. */
      void resize(size_type n) { 
        if(n>N) { throw std::length_error("array<T,N>::resize(size_type)"); }
        _size=n; }
      
      /*! \brief A reference to the \a n th element. */
      reference operator[](size_type i) { return _ptr[i]; }
      /*! \brief A constant reference to the \a n th element. */
      const_reference operator[](size_type i) const { return _ptr[i]; }
      /*! \brief Checked access to the \a n th element. */
      reference at(size_type i) { 
        if(i>=N) { throw std::out_of_range("array<T,N>::at(size_type)"); }
        return _ptr[i]; }
      /*! \brief Checked access to the \a n th element. */
      const_reference at(size_type i) const { 
        if(i>=N) { throw std::out_of_range("array<T,N>::at(size_type)"); }
        return _ptr[i]; }
      
      /*! \brief A reference to the first element of the array. */
      reference front() { return _ptr[0]; }
      /*! \brief A constant reference to the first element of the array. */
      const_reference front() const { return _ptr[0]; }
      /*! \brief A reference to the last element of the array. */
      reference back() { return _ptr[_size-1]; }
      /*! \brief A constantreference  to the last element of the array. */
      const_reference back() const { return _ptr[_size-1]; }
      
      /*! \brief An iterator pointing to the beginning of the array. */
      iterator begin() { return _ptr; }
      /*! \brief An iterator pointing to the beginning of the array. */
      const_iterator begin() const { return _ptr; }
      /*! \brief A constant iterator pointing to the end of the array. */
      iterator end() { return _ptr+_size; }
      /*! \brief A constant iterator pointing to the end of the array. */
      const_iterator end() const { return _ptr+_size; }
      
      /*! \brief Tests two arrays for equality */
      bool operator==(const array& other) const {
        const_iterator first=begin(); const_iterator last=end(); const_iterator curr=other.begin(); 
        while(first!=last) { if((*first)!=(*curr)) { return false; } ++first; ++curr; } return true; }
      /*! \brief Tests two arrays for inequality */
      bool operator!=(const array& other) const { return !((*this)==other); }

      /*! \brief Fills the array with copies of \a x. */
      void fill(value_type val) { 
        for(size_type i=0; i!=_size; ++i) { _ptr[i]=val; } }
      /*! \brief Fills the array from the sequence starting at \a first. */
      template<class In> void fill(In iter) { 
        _fill_iter(iter); }
      /*! \brief Assigns the sequence from \a first to \a last. */
      template<class In> void assign(In first, In last) { 
        if(distance(first,last)!=_size) { throw std::length_error("array<T,N>::assign(In,In)"); }
        _fill_iter(first); }
     private:
      template<class InputIterator> inline void _fill_iter(InputIterator);
     private:
      size_type _size;
      value_type _ptr[N];
    };
    
    template<class T, unsigned short int N> template<class InputIterator> inline
    void array<T,N>::_fill_iter(InputIterator iter)
    { 
      if(_size==1) { 
        *_ptr=*iter; }
      else if(_size==2) { 
        value_type* curr=_ptr; *curr=*iter; ++curr; ++iter;  *curr=*iter; }
      else if(_size==3) { 
        value_type* curr=_ptr; *curr=*iter; 
        ++curr; ++iter;  *curr=*iter; 
        ++curr; ++iter;  *curr=*iter; }
      else if(_size==4) { 
        value_type* curr=_ptr; *curr=*iter; 
        ++curr; ++iter; *curr=*iter; 
        ++curr; ++iter;  *curr=*iter;  
        ++curr; ++iter;  *curr=*iter; }
      else {
        iterator curr=_ptr; iterator this_end=this->end(); *curr=*iter;
        while(curr!=this_end) { ++curr; ++iter;  *curr=*iter;  } }
    }
    
    template<class T> template<unsigned short int N> 
    inline
    array<T>::array(const array<T,N>& a)
      : _size(a.size()), _ptr(new value_type[_size])
    { 
      fill(a.begin()); 
    }
    
    template<class T> template<unsigned short int N> 
    inline
    array<T>& array<T>::operator=(const array<T,N>& a) 
    { 
      resize(a.size()); fill(a.begin()); return *this; 
    }
    
    
    
    
    
    //! \ingroup Storage
    /*! \brief A range of values defined by a beginning and end iterator, which can be used as a reference to an array. */
    template<class Iter>
    class range {
     public:
      typedef Iter iterator;
      typedef typename std::iterator_traits<Iter>::value_type value_type;
      typedef typename std::iterator_traits<Iter>::reference reference;
      typedef typename std::iterator_traits<Iter>::pointer pointer;
      typedef typename std::iterator_traits<Iter>::difference_type difference_type;
      typedef typename array<value_type>::size_type size_type;
      
      /*! \brief Destructor. */
      ~range() { }
      /*! \brief Construct a range giving a reference to an array. */
      range(const array<value_type>& a) : _first(a.begin()), _last(a.end()) { }
      
      /*! \brief Construct a range from a beginning and end iterator. */
      template<class ForwardIterator>
      range(ForwardIterator first, ForwardIterator last)
        : _first(first), _last(last) { }
      
      /*! \brief Convert to an array by copying elements. */
      operator array<value_type>() const { return array<value_type>(begin(),end()); }
      
      /*! \brief Checks if the range is empty. */
      size_type empty() const { return (_first==_last); }
      /*! \brief The number of elements between the beginning and end iterators. */
      size_type size() const { return std::distance(_first,_last); }
      
      /*! \brief A reference to the \a i th element. */
      reference operator[](size_type i) const { return _first[i]; }
      /*! \brief A checked reference to the \a i th element. */
      reference at(size_type i) const { if(i<size()) { return _first[i]; } else { throw std::out_of_range("array index out-of-range"); } }
      
      /*! \brief A reference to the first element. */
      reference front() const { return _first[0]; }
      /*! \brief A reference to the last element. */
      reference back() const { return _last[-1]; }
      
      /*! \brief The beginning of the range. */
      iterator begin() const { return _first; }
      /*! \brief The end of the range. */
      iterator end() const { return _last; }
      
      /*! \brief Equality operator (compares the referenced arrays). */
      bool operator==(const range<iterator>& other) const { return equals(other.begin(),other.end()); }
      /*! \brief Inequality operator. */
      bool operator!=(const range<iterator>& other) const { return !((*this)==other); }
      /*! \brief Equality operator (compares the referenced arrays). */
      bool operator==(const array<value_type>& other) const { return equals(other.begin(),other.end()); }
      /*! \brief Inequality operator. */
      bool operator!=(const array<value_type>& other) const { return !((*this)==other); }
     private:
      template<class ForwardIterator>
      bool equals(ForwardIterator first, ForwardIterator last) const {
        if( size() != size_type(std::distance(first,last)) ) { return false; }
        iterator curr=_first;
        while(first!=last) { if((*first)!=(*curr)) { return false; } ++first; ++curr; }
        return true;
      }
     private:
      iterator _first, _last;
    };
    
    
    
    
    
    //! \ingroup Storage
    /*! \brief A reference to an array of values, defined by a size and a beginning iterator. */
    template<class Iter>
    class array_reference {
      typedef array_reference<Iter> Self;
     public:
      typedef typename std::iterator_traits<Iter>::value_type value_type;
      typedef typename std::iterator_traits<Iter>::reference reference;
      typedef typename std::iterator_traits<Iter>::pointer pointer;
      typedef typename array<value_type>::size_type size_type;
      typedef Iter iterator;
      
      /*! Destructor. */
      ~array_reference() { }
      /*! A refernce to an array. */
      array_reference(array<value_type>& a) : _size(a.size()), _begin(a.begin()) { }
      /*! A refernce to an array of size \a n beginning at \a b. */
      array_reference(const size_type& n, iterator b) : _size(n), _begin(b) { }
      /*! Copy constructor, allowing change of underlying iterator type. */
      template<class I> array_reference(const array_reference<I>& a) : _size(a.size()), _begin(a.begin()) { }
      
      /*! Assign the referenced array. */
      Self& operator=(const array<value_type>& a) { 
        if(this->size()!=a.size()) { throw std::length_error("array_reference<T>::operator=(array<T> const&)"); }
        _assign(a.begin(),a.end()); return *this; }
      /*! Assign the referenced array. */
      template<class Iter2> Self& operator=(const array_reference<Iter2>& a) { 
        if(this->size()!=a.size()) { throw std::length_error("array_reference<T>::operator=(array<X> const&)"); }
        _assign(a.begin(),a.end()); return *this; }
      /*! Assign the referenced array. */
      template<class RanIter> Self& operator=(const range<RanIter>& r) { 
        if(this->size()!=r.size()) { throw std::length_error("array_reference<T>::operator=(range<X> const&)"); }
        _assign(r.begin(),r.end()); return *this; }
      
      /*! Equality operator. */
      bool operator==(const array<value_type>& a) { return this->_equals(a.begin(),a.end()); }  
      /*! Inequality operator. */
      bool operator!=(const array<value_type>& a) { return !((*this)==a); }
      
      /*! Equality operator. */
      template<class I> bool operator==(const array_reference<I>& a) { return this->_equals(a.begin(),a.end()); }  
      /*! Inequality operator. */
      template<class I> bool operator!=(const array_reference<I>& a) { return !((*this)==a); }
      
      /*! Convert to an array with own storage. */
      operator array<value_type>() const { return array<value_type>(begin(),end()); }
      
      /*! \brief True if the referenced array's size is 0. */
      size_type empty() const { return _size==0u; }
      /*! \brief The size of the referenced array. */
      size_type size() const { return _size; }
      
      /*! \brief A reference to the \a i th element. */
      reference operator[](size_type i) const { return _begin[i]; }
      /*! \brief A checked reference to the \a i th element. */
      reference at(size_type i) const { if(i<size()) { return _begin[i]; } else { throw std::out_of_range("array index out-of-range"); } }
      
      /*! \brief A reference to the first element. */
      reference front() const { return (*this)[0]; }
      /*! \brief A reference to the last element. */
      reference back() const { return (*this)[size()-1]; }
      
      /*! \brief An iterator to the beginning of the referenced array. */
      iterator begin() const { return _begin; }
      /*! \brief An iterator to the end of the referenced array. */
      iterator end() const { return _begin+_size; }
     private:
      template<class FwdIter> void _assign(FwdIter b, FwdIter e) { 
        iterator curr=_begin; while(b!=e) { *curr=*b; ++curr; ++b; } }
      template<class FwdIter> bool _equals(FwdIter b, FwdIter e) { 
        if(std::distance(b,e)!=int(this->size())) { return false; }
        iterator curr=_begin; while(b!=e) { if(*curr!=*b) { return false; } ++curr; ++b; }
        return true;
      }
     private:
      size_type _size;
      iterator _begin;
    };
    
    
    
    
    
    template<class BaseIter>
    class _vector_array_iterator
      : public std::iterator<std::random_access_iterator_tag,
                             array_reference<BaseIter>,
                             array_reference<BaseIter>
                            >
    {
     private:
      typedef typename BaseIter::value_type element_type;
      typedef _vector_array_iterator<BaseIter> Self;
      typedef std::vector<BaseIter> position_iterator;
      typedef BaseIter element_iterator;
     public:
      typedef typename array<element_type>::size_type size_type;
      typedef typename BaseIter::difference_type difference_type;
      typedef array_reference<BaseIter> reference;
     private:
      position_iterator _position;
      element_iterator _element;
     private:
      size_type array_size() const { return (*(_position+1)-(*_position)); }
     public:
      _vector_array_iterator(position_iterator pos, element_iterator ptr)
        : _position(pos), _element(ptr) { }
      template<class B> _vector_array_iterator(const _vector_array_iterator<B>& iter) 
        : _position(iter._position), _element(iter._element) { }
      template<class B> bool operator==(const _vector_array_iterator<B>& other) const { 
        return _element==other._element && _position==other._position; }
      template<class B> bool operator!=(const _vector_array_iterator<B>& other) const { 
        return !(*this==other); }
      reference operator*() const { return reference(array_size(),_element); }
      Self& operator++() { _element+=array_size(); ++_position; return *this; }
      Self& operator--() { --_position; _element-=array_size(); return *this; }
      Self& operator+=(difference_type i) { _element+=(_position[i]-_position[0]); _position+=i; return *this; }
      Self operator+(difference_type i) const { Self tmp=*this; tmp+=i; return tmp; }
      Self operator-(difference_type i) const { Self tmp=*this; tmp+=(-i); return tmp; }
      difference_type operator-(const Self& other) const { return (this->_position-other._position); }
    };  
    
    
    //! \ingroup Storage
    /*! \brief A vector of arrays of different sizes. */
    template<class T>
    class vector_array {
     private:
      typedef typename std::vector<T>::iterator vector_iterator;
      typedef typename std::vector<T>::const_iterator vector_const_iterator;
     public:
      typedef T array_value_type;
      typedef typename std::vector<T>::size_type array_size_type;
      typedef typename std::vector<T>::size_type size_type;
      typedef array_reference<vector_iterator> reference;
      typedef array_reference<vector_const_iterator> const_reference;
      typedef _vector_array_iterator<vector_iterator> iterator;
      typedef _vector_array_iterator<vector_const_iterator> const_iterator;
      typedef vector_iterator element_iterator;
      typedef vector_const_iterator element_const_iterator;
     public:
      /*!\brief Construct a vector to hold arrays of size as. */
      vector_array() : _positions(), _elements() { _positions.push_back(0); }
      
      /*!\brief Copy constructor. */
      vector_array(const vector_array<T>& av)
        : _positions(av._positions), _elements(av._elements) { }
      
      /*!\brief Assignment operator. */
      vector_array<T>& operator=(const vector_array<T>& av) {
        if(this!=&av) { _positions=av._positions; _elements=av._elements; } return *this; }
        
      /*!\brief The total number of array elements. */
      size_type length() const { return _positions.back(); }
      
      /*!\brief The number of arrays in the vector.  */
      size_type size() const { return _positions.size()-1; }
      
      /*!\brief True if the vector is empty.  */
      bool empty() const { return _elements.empty(); }
      
      /*!\brief The size of the \a i th array. */
      array_size_type array_size(size_type i) const { return _positions[i+1]-_positions[i]; }
      
      /*!\brief The size of each array word in the list. */
      void clear() { _positions.clear(); _elements.clear(); _positions.push_back(0); }
      
      /*!\brief Insert an element at the back of the vector. */
      void push_back(const array<T>& a) { 
        _positions.push_back(_positions.back()+a.size());
        for(size_type i=0; i!=a.size(); ++i) { _elements.push_back(a[i]); }
      }
      
      /*!\brief Removes the element at the back of the vector. */
      void pop_back() { _positions.pop_back(); _elements.resize(_positions.back()); }
      
      /*!\brief A reference to the array at the front of the vector. */
      reference front() { return operator[](0); }
      /*!\brief A reference to the array at the front of the vector. */
      const_reference front() const { return operator[](0); }
      
      /*!\brief A reference to the array at the back of the vector. */
      reference back() { return operator[](size()-1); }
      /*!\brief A reference to the array at the back of the vector. */
      const_reference back() const { operator[](size()-1); }

      /*!\brief Returns the nth array of the vector. */
      reference operator[] (size_type i) {
        return reference(_positions[i+1]-_positions[i],_elements.begin()+_positions[i]);
      }
      
      /*!\brief Returns the nth array of the vector. */
      const_reference operator[] (size_type i) const {
        return const_reference(_positions[i+1]-_positions[i],_elements.begin()+_positions[i]);
      }
      
      /*!\brief Checked access to the nth array of the vector. */
      reference at(size_type i) { 
        if(i<size()) { return this->operator[](i); } else { throw std::out_of_range("array index out-of-range"); } 
      }
      
      /*!\brief Checked access to the nth array of the vector. */
      const_reference at(size_type i) const { 
        if(i<size()) { return this->operator[](i); } else { throw std::out_of_range("array index out-of-range"); }
      }
      
      /*!\brief Efficiently swap two array vectors.  */
      void swap(vector_array<T>& other) {
        _positions.swap(other._positions); _elements.swap(other._elements); }
        
      /*!\brief A random-access constant iterator pointing to the beginning of the vector of arrays.  */
      const_iterator begin() const { return const_iterator(_positions.begin(),_elements.begin()); }
      
      /*!\brief A random-access constant iterator pointing to the end of the vector of arrays. */
      const_iterator end() const { return const_iterator(_positions.end(),_elements.end()); }

      /*!\brief A random-access constant iterator pointing to the beginning of the vector of arrays.  */
      iterator begin() { return iterator(_positions.begin(),_elements.begin()); }
      
      /*!\brief A random-access constant iterator pointing to the end of the vector of arrays. */
      iterator end() { return iterator(_positions.end(),_elements.end()); }
     private:
      std::vector<element_iterator> _positions;
      std::vector<T> _elements;
    };
    
    
    
    
    
    
    
    template<class BaseIter>
    class _array_vector_iterator
      : public std::iterator<std::random_access_iterator_tag,
                             array_reference<BaseIter>,
                             array_reference<BaseIter>
                            >
    {
     private:
      typedef typename BaseIter::value_type element_type;
      typedef _array_vector_iterator<BaseIter> Self;
      typedef BaseIter element_iterator;
     public:
      typedef typename array<element_type>::size_type size_type;
      typedef typename BaseIter::difference_type difference_type;
      typedef array_reference<BaseIter> reference;
     private:
      size_type _array_size;
      element_iterator _curr;
     public:
      _array_vector_iterator() { }
      _array_vector_iterator(size_type n, element_iterator ptr)
        : _array_size(n), _curr(ptr) { }
      template<class B> _array_vector_iterator(const _array_vector_iterator<B>& iter) 
        : _array_size((*iter).size()), _curr((*iter).begin()) { }
      template<class B> bool operator==(const _array_vector_iterator<B>& other) const { 
        return _curr==(*other).begin(); }
      template<class B> bool operator!=(const _array_vector_iterator<B>& other) const { 
        return !(*this==other); }
      reference operator*() const { return reference(_array_size,_curr); }
      Self& operator++() { _curr+=_array_size; return *this; }
      Self& operator--() { _curr-=_array_size; return *this; }
      Self& operator+=(difference_type i) { _curr+=(i*_array_size); return *this; }
      Self operator+(difference_type i) const { Self tmp=*this; tmp+=i; return tmp; }
      Self operator-(difference_type i) const { Self tmp=*this; tmp+=(-i); return tmp; }
      difference_type operator-(const Self& other) const { return (this->_curr-other._curr)/_array_size; }
    };
    
    
    //! \ingroup Storage
    /*! \brief A vector of arrays of a fixed size. */
    template<class T>
    class array_vector {
     private:
      typedef typename std::vector<T>::iterator vector_iterator;
      typedef typename std::vector<T>::const_iterator vector_const_iterator;
     public:
      typedef T array_value_type;
      typedef typename std::vector<T>::size_type array_size_type;
      typedef typename std::vector<T>::size_type size_type;
      typedef array_reference<vector_iterator> reference;
      typedef array_reference<vector_const_iterator> const_reference;
      typedef _array_vector_iterator<vector_iterator> iterator;
      typedef _array_vector_iterator<vector_const_iterator> const_iterator;
      typedef vector_iterator element_iterator;
      typedef vector_const_iterator element_const_iterator;
     public:
      /*!\brief Construct a vector to hold arrays of size as. */
      array_vector(array_size_type as) : _array_size(as), _elements() { }
      
      /*!\brief Copy constructor. */
      array_vector(const array_vector<T>& av)
        : _array_size(av._array_size), _elements(av._elements) { }
      
      /*!\brief Assignment operator. */
      array_vector<T>& operator=(const array_vector<T>& av) {
        if(this!=&av) { _array_size=av._array_size; _elements=av._elements; } return *this; }
        
      /*!\brief The total number of array elements. */
      size_type length() const { return _elements.size(); }
      
      /*!\brief The number of arrays in the vector.  */
      size_type size() const { return _elements.size()/array_size(); }
      
      /*!\brief True if the vector is empty.  */
      bool empty() const { return _elements.empty(); }
      
      /*!\brief The size of each array word in the list. */
      array_size_type array_size() const { return _array_size; }
      
      /*!\brief Resizes the array to hold \a n elements. */
      void resize(size_type n) { _elements.resize(_array_size*n); }
      
      /*!\brief The size of each array word in the list. */
      void clear() { _elements.clear(); }
      
      /*!\brief Efficiently swap two array vectors.  */
      void swap(array_vector<T>& other) {
        std::swap(_array_size,other._array_size); _elements.swap(other._elements); }
        
      /*!\brief Insert an element at the back of the vector. */
      void push_back(const array<T>& a) { 
        if(a.size()!=array_size()) { throw std::length_error("array_vector<T>::push_back(array<T> const&)"); }
        for(size_type i=0; i!=array_size(); ++i) { _elements.push_back(a[i]); }
      }
      
      /*!\brief Removes the element at the back of the vector. */
      void pop_back() { _elements.resize(_elements.size()-array_size()); }
      
      /*!\brief A reference to the array at the front of the vector. */
      reference front() { return this->operator[](0); }
      /*!\brief A reference to the array at the front of the vector. */
      const_reference front() const { return this->operator[](0); }
      
      /*!\brief A reference to the array at the back of the vector. */
      reference back() { return this->operator[](this->size()-1); }
      /*!\brief A reference to the array at the back of the vector. */
      const_reference back() const { return this->operator[](this->size()-1); }
      
      /*!\brief Returns the nth array of the vector. */
      reference operator[] (size_type i) {
        vector_iterator first=_elements.begin()+(i*array_size()); 
        return reference(array_size(),first);
      }
      
      /*!\brief Returns the nth array of the vector. */
      const_reference operator[] (size_type i) const {
        vector_const_iterator first=_elements.begin()+(i*array_size()); 
        return const_reference(array_size(),first); 
      }
      
      /*!\brief Checked access to the nth array of the vector. */
      reference at(size_type i) { 
        if(i<size()) { return this->operator[](i); } else { throw std::out_of_range("array index out-of-range"); } 
      }
      
      /*!\brief Checked access to the nth array of the vector. */
      const_reference at(size_type i) const { 
        if(i<size()) { return this->operator[](i); } else { throw std::out_of_range("array index out-of-range"); }
      }
      
      /*!\brief A random-access constant iterator pointing to the beginning of the vector of arrays.  */
      const_iterator begin() const { return const_iterator(array_size(),_elements.begin()); }
      
      /*!\brief A random-access constant iterator pointing to the end of the vector of arrays. */
      const_iterator end() const { return const_iterator(array_size(),_elements.end()); }

      /*!\brief A random-access constant iterator pointing to the beginning of the vector of arrays.  */
      iterator begin() { return iterator(array_size(),_elements.begin()); }
      
      /*!\brief A random-access constant iterator pointing to the end of the vector of arrays. */
      iterator end() { return iterator(array_size(),_elements.end()); }
     private:
      array_size_type _array_size;
      std::vector<T> _elements;
    };
    
    
    
  } // namespace Base
} // namespace Ariadne
  
#endif /* ARIADNE_ARRAY_H */
