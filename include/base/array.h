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

#ifndef _ARIADNE_ARRAY_H
#define _ARIADNE_ARRAY_H

#include <cstddef>

#include <exception>
#include <stdexcept>
#include <vector>
#include <iterator>

#include "../declarations.h"


namespace Ariadne {
  
  namespace Base {
    
    /*! \brief STL style interface to sized arrays. 
     *
     * An array<T> is a variable-size array which can be resized and is allocated 
     * on the heap. Arrays provide checked access using at and unchecked access using operator[].
     *
     * An array<T,N> is a fixed-size array of size \a N which is allocated 
     * on the stack. Arrays provide checked access using at and unchecked access using operator[].
     *
     * FIXME: Constructors for small arrays.
     * FIXME: Exception safety in constructors!
     */
    template<typename T, size_t N> class array;
    
    // Need to give default size in first declaration.
    //template<typename T, size_t N=0> class array;
    
    template<typename T> class array<T> {
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
      
      /*! \brief Copy constructor. */
      array(const array<T>& a) : _size(a.size()), _ptr(new value_type[_size]) { fill(a.begin()); }
      /*! \brief Copy assignment. */ 
      array<T>& operator=(const array<T>& a) { resize(a.size()); fill(a.begin()); return *this; }
      /*! \brief Conversion from an array of fixed size \p N. */ 
      template<size_type N> array(const array<T,N>& a); /* inline conversion */ 
      /*! \brief Assignment from an array of fixed size \p N. */ 
      template<size_type N> array<T>& operator=(const array<T,N>& a); /* inline conversion */ 
      
      /*! \brief True if the array's size is 0. */
      bool empty() const { return _size==0u; }
      /*! \brief The size of the array. */
      size_type size() const { return _size; }
      /*! \brief The maximum possible size of the array. */
      size_type max_size() const { return (size_type) (-1); }
      /*! \brief Resizes the array to hold \a n elements. If \a n is larger than the current size, the extra elements are uninitialised. */
      void resize(size_type n) { if(size()!=n) { delete[] _ptr; _size=n; _ptr=new value_type[_size]; } }
      
      /*! \brief The \a n'th element. */
      reference operator[](size_type i) { return _ptr[i]; } 
      const_reference operator[](size_type i) const { return _ptr[i]; }
      /*! \brief Checked access to the \a n'th element. */
      reference at(size_type i) { if(i<_size) { return _ptr[i]; } else { throw std::out_of_range("array index out-of-range"); } } 
      const_reference at(size_type i) const { if(i<_size) { return _ptr[i]; } else { throw std::out_of_range("array index out-of-range"); } }
      
      /*! \brief An iterator pointing to the beginning of the array. */
      iterator begin() { return _ptr; }
      const_iterator begin() const { return _ptr; }
      /*! \brief An iterator pointing to the end of the array. */
      iterator end() { return _ptr+_size; }
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
      
      reference operator[](size_type i) { return _vector[i]; }
      const_reference operator[](size_type i) const { return _vector[i]; }
      reference at(size_type i) { if(i<size()) { return _vector[i]; } else { throw std::out_of_range("array index out-of-range"); } } 
      const_reference at(size_type i) const { if(i<size()) { return _vector[i]; } else { throw std::out_of_range("array index out-of-range"); } }
      
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
    
    template<typename T, size_t N> template<class InputIterator> inline
    void array<T,N>::_assign_iter(InputIterator iter)
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
        iterator curr=_ptr; iterator this_end=this->end(); *curr=*iter;
        while(curr!=this_end) { ++curr; ++iter;  *curr=*iter;  } }
    }
    
    
    template<typename T> template<size_t N> 
    inline
    array<T,0>::array(const array<T,N>& a)
      : _size(a.size()), _ptr(new value_type[_size])
    { 
      fill(a.begin()); 
    }
    
    template<typename T> template<size_t N> 
    inline
    array<T,0>& array<T,0>::operator=(const array<T,N>& a) 
    { 
      resize(a.size()); fill(a.begin()); return *this; 
    }
    
    
    
    
    
    
    
    /*! \brief A range of values, which can be used as a reference to an array. */
    template<class RandomAccessIterator>
    class range {
     public:
      typedef RandomAccessIterator iterator;
      typedef typename std::iterator_traits<RandomAccessIterator>::value_type value_type;
      typedef typename std::iterator_traits<RandomAccessIterator>::reference reference;
      typedef typename std::iterator_traits<RandomAccessIterator>::pointer pointer;
      typedef typename std::iterator_traits<RandomAccessIterator>::difference_type difference_type;
      typedef typename array<value_type>::size_type size_type;
      
      ~range() { }
      range(const array<value_type>& a) : _first(a.begin()), _last(a.end()) { }
      
      template<typename ForwardIterator>
      range(ForwardIterator first, ForwardIterator last)
        : _first(first), _last(last) { }
      
      operator array<value_type>() const { return array<value_type>(begin(),end()); }
      
      size_type empty() const { return (_first==_last); }
      size_type size() const { return std::distance(_first,_last); }
      
      reference operator[](size_type i) const { return _first[i]; }
      reference at(size_type i) const { if(i<size()) { return _first[i]; } else { throw std::out_of_range("array index out-of-range"); } }
      
      iterator begin() const { return _first; }
      iterator end() const { return _last; }
      
      bool operator==(const range<iterator>& other) const { return equals(other.begin(),other.end()); }
      bool operator==(const array<value_type>& other) const { return equals(other.begin(),other.end()); }
      bool operator!=(const range<iterator>& other) const { return !((*this)==other); }
      bool operator!=(const array<value_type>& other) const { return !((*this)==other); }
     private:
      template<typename ForwardIterator>
      bool equals(ForwardIterator first, ForwardIterator last) const {
        if( size() != size_type(std::distance(first,last)) ) { return false; }
        iterator curr=_first;
        while(first!=last) { if((*first)!=(*curr)) { return false; } ++first; ++curr; }
        return true;
      }
     private:
      iterator _first, _last;
    };
    
    
    
    
    
    template<typename T>
    class _array_vector_iterator
      : public std::iterator<std::random_access_iterator_tag, T>
    {
     private:
      typedef typename std::vector<T>::const_iterator element_iterator;
      typedef _array_vector_iterator<T> Self;
     public:
      typedef typename std::vector<T>::size_type size_type;
      typedef typename std::vector<T>::difference_type difference_type;
      typedef range<element_iterator> reference;
     private:
      size_type _array_size;
      element_iterator _curr;
     public:
      _array_vector_iterator(size_type n, element_iterator ptr)
        : _array_size(n), _curr(ptr) { }
      bool operator==(const Self& other) const { return _curr==other._curr; }
      bool operator!=(const Self& other) const { return !(*this==other); }
      reference operator*() const { return reference(_curr,_curr+_array_size); }
      Self& operator++() { _curr+=_array_size; return *this; }
      Self& operator--() { _curr-=_array_size; return *this; }
      Self& operator+=(difference_type i) { _curr+=(i*_array_size); return *this; }
      Self operator+(difference_type i) const { Self tmp=*this; tmp+=i; return tmp; }
    };
    
    template<class T> class array_vector;
    template<typename T> std::ostream& operator<<(std::ostream&, const array_vector<T>&);
    
    
    /*! \brief A vector of arrays of a fixed size. */
    template<class T>
    class array_vector {
      friend std::ostream& operator<< < >(std::ostream&, const array_vector<T>&);
     private:
      typedef typename std::vector<T>::iterator Vector_iterator;
      typedef typename std::vector<T>::const_iterator Vector_const_iterator;
     public:
      typedef T array_value_type;
      typedef typename std::vector<T>::size_type array_size_type;
      typedef typename std::vector<T>::size_type size_type;
      typedef range<Vector_iterator> reference;
      typedef range<Vector_const_iterator> const_reference;
      typedef _array_vector_iterator<T> const_iterator;
     public:
      /*!\brief Construct a vector to hold arrays of size as. */
      array_vector(array_size_type as) : _array_size(as), _elements() { }
      
      /*!\brief The total number of array elements. */
      size_type length() const { return _elements.size(); }
      
      /*!\brief The number of arrays in the vector.  */
      size_type size() const { return _elements.size()/array_size(); }
      
      /*!\brief True if the vector is empty.  */
      bool empty() const { return _elements.empty(); }
      
      /*!\brief The size of each array word in the list. */
      array_size_type array_size() const { return _array_size; }
      
      /*!\brief The size of each array word in the list. */
      void resize(size_type n) { _elements.resize(_array_size*n); }
      
      /*!\brief Insert an element at the back of the vector. */
      void push_back(const array<T>& a) { 
        assert(a.size()==array_size()); 
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
        Vector_iterator first=_elements.begin()+(i*array_size()); 
        return reference(first,first+array_size()); 
      }
      
      /*!\brief Returns the nth array of the vector. */
      const_reference operator[] (size_type i) const {
        Vector_const_iterator first=_elements.begin()+(i*array_size()); 
        return const_reference(first,first+array_size()); 
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
     private:
      array_size_type _array_size;
      std::vector<T> _elements;
    };
    
  } // namespace Base
} // namespace Ariadne
  
#endif /* _ARIADNE_ARRAY_H */
