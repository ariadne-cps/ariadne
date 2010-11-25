/***************************************************************************
 *            array.h
 *
 *  Copyright 2008  Pieter Collins
 *
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
 *  \brief STL-style fixed-size arrays.
 */

#ifndef ARIADNE_ARRAY_H
#define ARIADNE_ARRAY_H

#include <cstddef>
#include <cstdarg>
#include <stdexcept>
#include <new>

namespace Ariadne {

template<class T>
class array {
  private:
    static T* uninitialized_new(size_t n) { return static_cast<T*>(::operator new(n*sizeof(T))); }
    static void uninitialized_delete(T* p) { ::operator delete(p); }
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
    ~array() { this->_destroy_elements(); uninitialized_delete(_ptr); }

    /*! \brief Default constructor. Constructs an empty array. */
    array() : _size(0), _ptr(0) { }
    /*! \brief Constructs an array of size \a n with uninitialised elements. */
    explicit array(const size_type n) : _size(n), _ptr(uninitialized_new(n)) { for(size_type i=0; i!=n; ++i) { new (_ptr+i) T(); } }
    /*! \brief Constructs an array of size \a n with elements initialised to \a x. */
    array(const size_type n, const value_type& x) : _size(n), _ptr(uninitialized_new(n)) { this->_uninitialized_fill(x); }
    /*! \brief Constructs an array of size \a n with elements initialised by the variable argument list x0,x1,... . */
    array(const size_type n, const value_type& x0, const value_type& x1, ...)
        : _size(n), _ptr(uninitialized_new(n))
    {
        if(n<2) { throw std::out_of_range("array: size is less than number of arguments"); }
        va_list args; va_start(args,x1); new (_ptr) T(x0); new (_ptr+1u) T(x1);
        for(size_t i=2; i!=n; ++i) { new (_ptr+i) T(va_arg(args,T)); } va_end(args); }
    /*! \brief Constructs an array from the range \a first to \a last. */
    template<class ForwardIterator> array(ForwardIterator first, ForwardIterator last)
        : _size(std::distance(first,last)), _ptr(uninitialized_new(_size))
    {
        this->_uninitialized_fill(first);
    }
    /*! \brief Conversion constructor. */
    template<class T1> array(const array<T1>& a) : _size(a.size()), _ptr(uninitialized_new(_size)) {
        this->_uninitialized_fill(a.begin()); }
    /*! \brief Copy constructor. */
    array(const array<T>& a) : _size(a.size()), _ptr(uninitialized_new(_size)) {
        this->_uninitialized_fill(a.begin()); }
    /*! \brief Copy assignment. */
    array<T>& operator=(const array<T>& a) {
        if(this->size()==a.size()) { fill(a.begin()); }
        else { this->_destroy_elements(); uninitialized_delete(_ptr); _size=a.size(); _ptr=uninitialized_new(_size); this->_uninitialized_fill(a.begin()); }
        return *this; }

    /*! \brief True if the array's size is 0. */
    bool empty() const { return _size==0u; }
    /*! \brief The size of the array. */
    size_type size() const { return _size; }
    /*! \brief The maximum possible size of the array. */
    size_type max_size() const { return (size_type) (-1); }
    /*! \brief Resizes the array to hold \a n elements. If \a n is larger than the current size, the extra elements are default initialised. */
    void resize(size_type n) {
        if(size()!=n) {
            pointer _new_ptr=uninitialized_new(n);
            for(size_type i=0; i!=n; ++i) { if(i<_size) { new (_new_ptr+i) T(_ptr[i]); } else { new (_new_ptr+i) T(); } }
            this->_destroy_elements(); uninitialized_delete(_ptr); _size=n; _ptr=_new_ptr; } }
    /*! \brief Reallocates the array to hold \a n elements. The new elements are default-constructed. */
    void reallocate(size_type n) { if(size()!=n) { this->_destroy_elements(); uninitialized_delete(_ptr);
        _size=n; _ptr=uninitialized_new(_size); for(size_type i=0; i!=_size; ++i) { new (_ptr+i) T(); } } }
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
        iterator curr=begin(); iterator end=this->end(); while(curr!=end) { *curr=x; ++curr; } }
    /*! \brief Fills the array from the sequence starting at \a first. */
    template<class InputIterator> void fill(InputIterator first) {
        iterator curr=begin(); iterator end=this->end(); while(curr!=end) { *curr=*first; ++curr; ++first; } }
    /*! \brief Assigns the sequence from \a first to \a last. */
    template<class ForwardIterator> void assign(ForwardIterator first, ForwardIterator last) {
        resize(std::distance(first,last)); fill(first); }
  private:
    void _destroy_elements() { pointer curr=_ptr+_size; while(curr!=_ptr) { --curr; curr->~T(); } }
    void _uninitialized_fill(const value_type& x) {
        pointer curr=_ptr; pointer end=_ptr+_size; while(curr!=end) { new (curr) T(x); ++curr; } }
     template<class InputIterator> void _uninitialized_fill(InputIterator first) {
        pointer curr=_ptr; pointer end=_ptr+_size;
        while(curr!=end) { new (curr) T(*first); ++curr; ++first; } }
  private:
    size_type _size;
    pointer _ptr;
};


template<class A, class X> void serialize(A& a, array<X>& ary, const unsigned int v) {
    // We can't use separate save/load unless serialize is a member (I think).
    // We therefore need the same code to read and write.
    // We use a trick by which we first store the current value in a temporary
    // variable, read/write the temporary, and set it back to the archive.
    size_t m=ary.size(); a & m; if(m!=ary.size()) { ary.resize(m); }
    for(size_t i=0; i!=m; ++i) { X& k=ary[i]; a & k; }
}



} // namespace Ariadne;

#endif /* ARIADNE_ARRAY_H */
