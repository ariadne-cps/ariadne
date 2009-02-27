/***************************************************************************
 *            multi_index.h
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
 
/*! \file multi_index.h
 *  \brief An index specifying the degree of differentiation.
 */
#ifndef ARIADNE_MULTI_INDEX_H
#define ARIADNE_MULTI_INDEX_H

#include <cassert>
#include <boost/iterator.hpp>
#include <boost/iterator_adaptors.hpp>
#include "macros.h"
#include "array.h"

namespace Ariadne {

uint32_t fac(uint8_t);
uint32_t bin(uint8_t,uint8_t);

typedef unsigned char uchar;

class MultiIndex {
  public:
    typedef unsigned int size_type;
    typedef unsigned char value_type;
    typedef unsigned char word_type;
    //typedef unsigned int word_type;
    //typedef unsigned long long int word_type;
  public:
    /*! Destructor. */
    ~MultiIndex() { _deallocate(this->_p); }
    /*! Construct a multi index of degree \a 0 with \a nv variables. */
    explicit MultiIndex(size_type nv);
    /*! Construct a multi index with \a nv variables from the array \a ary. */
    explicit MultiIndex(size_type nv, const uint* ary);
    /*! Construct a multi index with \a nv variables from variable arguments. */
    explicit MultiIndex(size_type nv, int a1, ...);

    /*! Copy constructor. */
    MultiIndex(const MultiIndex& a);
    /*! Copy assignment operator. */
    MultiIndex& operator=(const MultiIndex& a);
  
    /*! Construct the zero multi index with \a nv variables. */
    static MultiIndex zero(size_type nv);
    /*! Construct the unit multi index in variable \a j with \a nv variables. */
    static MultiIndex unit(size_type nv, size_type j);
    /*! Construct the first multi index of degree \a d with \a nv variables. */
    static MultiIndex first(size_type nv, value_type d);

    /*! Resize to hold n variables. */
    void resize(size_type n);
    /*! Set all values to zero. */
    void clear();
    /*! The number of variables. */
    size_type size() const;
    /*! The degree of the multi-index, equal to the sum of the number of occurrences of the variables. */
    value_type degree() const;
     /*! The number of variables. */
    size_type number_of_variables() const;
    /*! The number of occurrences of the \a i th variable. */
    value_type const& operator[](size_type i) const;
  
    /*! Equality operator. */
    bool operator==(const MultiIndex& a) const;
    /*! Inequality operator. */
    bool operator!=(const MultiIndex& a) const;
    /*! Comparison operator. */
    bool operator<(const MultiIndex& a) const; // inline
  
    /*! Increment. No post-increment operator as we sometimes pass MultiIndex by reference. */
    MultiIndex& operator++(); // inline
    // No post-increment operator as we sometimes pass MultiIndex by reference.Post increment. 
    // MultiIndex operator++(int); 
    /*! Inplace sum. */
    MultiIndex& operator+=(const MultiIndex& a); // inline
    /*! Inplace difference. */
    MultiIndex& operator-=(const MultiIndex& a); // inline
    /*! Inplace scalar product. */
    MultiIndex& operator*=(const value_type& a); // inline
    /*! Sum. */
    friend MultiIndex operator+(const MultiIndex& a1, const MultiIndex& a2); // inline
    /*! Difference. */
    friend MultiIndex operator-(const MultiIndex& a1, const MultiIndex& a2); // inline
    /*! Scalar product. */
    friend MultiIndex operator*(const MultiIndex& a, value_type s); // inline
    friend MultiIndex operator*(value_type s, const MultiIndex& a); // inline
  
    /*! The position of the element in the array of tensor values. */
    uint position() const;
  
    /*! The product of the factorials of the indices. */
    uint factorial() const;
    /*! The number of ordered index arrays with each element occurring the number of times specified by the multi index. */
    uint number() const;
  
    /*! Set the value of the \a i th index to \a j. */
    void set_index(size_type i, value_type j);
    /*! Increment the the \a i th index, thereby increasing the degree. */
    void increment_index(size_type i);
    /*! Decrement the the \a i th index, thereby decreasing the degree. */
    void decrement_index(size_type i);
  
    /*! Set the value of the \a i th index to \a j. */
    void set(size_type i, value_type j);
    friend std::ostream& operator<<(std::ostream&, const MultiIndex&);
  private:
    size_type word_size() const { return _word_size(_n); }
    word_type& word_at(size_type j) { return reinterpret_cast<word_type*>(_p)[j]; }
    const word_type& word_at(size_type j) const { return reinterpret_cast<const word_type*>(_p)[j]; }
    const value_type& at(size_type i) const { return _p[i]; }
  private:
    static size_type _word_size(size_type n) { return ((n*sizeof(value_type))/sizeof(word_type)+1); }
    static value_type* _allocate(size_type n) { return reinterpret_cast<value_type*>(new word_type[_word_size(n)]); }
    static void _deallocate(value_type* p) { delete[] reinterpret_cast<word_type*>(p); }
  private:
    size_type _n;
    value_type* _p;
};
      


    
inline MultiIndex::MultiIndex(size_type n)
    : _n(n), _p(_allocate(_n)) 
{
    for(uint j=0; j!=word_size(); ++j) { word_at(j)=0; }
}

inline MultiIndex::MultiIndex(size_type n, const uint* ary)
    : _n(n), _p(_allocate(_n))
{
    for(uint j=0; j!=word_size(); ++j) { word_at(j)=0; }
    _p[n]=0; for(uint i=0; i!=n; ++i) { _p[i]=ary[i]; _p[n]+=ary[i]; } 
}

inline MultiIndex::MultiIndex(const MultiIndex& a)
    : _n(a.size()), _p(_allocate(_n))
{
    for(uint j=0; j!=word_size(); ++j) { this->word_at(j)=a.word_at(j); }
}

inline MultiIndex& MultiIndex::operator=(const MultiIndex& a) { 
    if(this!=&a) { this->resize(a.size());
        for(uint j=0; j!=word_size(); ++j) { this->word_at(j)=a.word_at(j); } }
    return *this;
}

inline MultiIndex MultiIndex::zero(size_type n)
{
    return MultiIndex(n);
}

inline MultiIndex MultiIndex::unit(size_type n, size_type j)
{
    MultiIndex result(n);
    result._p[j]=1u;
    result._p[n]=1u;
    return result;
}

inline MultiIndex MultiIndex::first(size_type n, value_type d)
{
    MultiIndex result(n);
    result._p[0]=d;
    result._p[n]=d;
    return result;
}

inline void MultiIndex::clear()
{
    for(uint j=0; j!=word_size(); ++j) { this->word_at(j)=0; }
}

inline void MultiIndex::resize(size_type n) { 
    if(this->_n!=n) { _deallocate(this->_p); this->_n=n; this->_p=_allocate(this->_n); }
}

inline MultiIndex::size_type MultiIndex::size() const { 
    return this->_n; 
}

inline MultiIndex::value_type MultiIndex::degree() const { 
    return this->_p[this->_n]; 
}

inline MultiIndex::size_type MultiIndex::number_of_variables() const { 
    return this->_n; 
}

inline MultiIndex::value_type const& MultiIndex::operator[](size_type i) const {
    assert(i<this->size()); return this->_p[i];
}


inline bool MultiIndex::operator==(const MultiIndex& a) const { 
    if(this->_n!=a._n) { return false; }
    for(size_type i=0; i!=this->size(); ++i) { 
        if(this->at(i)!=a.at(i)) { return false; } }
    return true;
    for(size_type j=0; j!=this->word_size(); ++j) { 
        if(this->word_at(j)!=a.word_at(j)) { return false; } }
    return true;
}

inline bool MultiIndex::operator!=(const MultiIndex& a) const { 
    return !(*this==a); 
}

inline void MultiIndex::set_index(size_type i, const value_type j) { 
    this->_p[this->_n]+=(j-this->_p[i]); this->_p[i]=j; 
}

inline void MultiIndex::increment_index(size_type i) { 
    ++this->_p[i]; ++this->_p[this->_n]; 
}

inline void MultiIndex::decrement_index(size_type i) { 
    if(!(this->_p[i]>0)) { 
        throw std::runtime_error("MultiIndex::decrement_index: the number of occurence of the index must be positive"); 
    }
    --this->_p[i]; 
    --this->_p[this->_n]; 
}

inline void MultiIndex::set(size_type i, value_type j) { 
    this->set_index(i,j); 
}





inline
uint MultiIndex::number() const
{
    uint result=fac(this->degree());
    for(uint k=0; k!=this->size(); ++k) {
        result/=fac((*this)[k]);
    }
    return result;
}

inline
uint MultiIndex::factorial() const
{
    uint result=1;
    for(uint k=0; k!=this->size(); ++k) {
        result*=fac((*this)[k]);
    }
    return result;
}

inline
uint MultiIndex::position() const
{
    uint deg=this->degree()-1;
    uint nvar=this->size();
    uint result=bin(deg+nvar,nvar);
    for(uint k=0; k!=this->size()-1; ++k) {
        --nvar;
        deg-=(*this)[k];
        result+=bin(deg+nvar,nvar);
    }
    return result;
}

inline
bool MultiIndex::operator<(const MultiIndex& a2) const {
    const MultiIndex& a1=*this;
  
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("operator<(MultiIndex,MultiIndex): number of variables must match");
    }
  
    if(a1.degree()!=a2.degree()) {
        return a1.degree()<a2.degree();
    } else {
        //for(size_type j=0; j!=a1.word_size(); ++j) {
        for(int j=a1.word_size()-1; j!=-1; --j) {
            if(a1.word_at(j)!=a2.word_at(j)) { 
                return a1.word_at(j)>a2.word_at(j);
            }
        }
        return false;
        for(size_type i=0; i!=a1.size(); ++i) {
            if(a1[i]!=a2[i]) { 
                return a1[i]>a2[i];
            }
        }
        return false;
    }
}


inline
MultiIndex& MultiIndex::operator++() 
{
    MultiIndex& a=*this;
    const uint& nv=a.size();
    assert(nv>0);
    if(nv==1) {
        a.increment_index(0);
        return a;
    }
    if(a[nv-2]!=0) {
        a.decrement_index(nv-2);
        a.increment_index(nv-1);
        return a;
    } else {
        uint li=a[nv-1];
        a.set_index(nv-1,0);
        for(uint k=nv-1; k!=0; --k) {
            if(a[k-1]!=0) {
                a.decrement_index(k-1);
                a.set_index(k,li+1);
                return a;
            }
        }
        a.set_index(0,li+1);
    }
    return a;
}

inline
MultiIndex& MultiIndex::operator+=(const MultiIndex& a) {
    for(size_type j=0; j!=this->word_size(); ++j) {
        this->word_at(j)+=a.word_at(j);
    }
    return *this;
    for(size_type i=0; i!=this->size(); ++i) {
        this->_p[i]+=a._p[i];
    }
}

inline
MultiIndex& MultiIndex::operator-=(const MultiIndex& a) {
    for(size_type i=0; i!=this->size()+1; ++i) {
        ARIADNE_ASSERT(this->_p[i]>=a._p[i]);
        this->_p[i]-=a._p[i];
    }
    return *this;
}

inline
MultiIndex& MultiIndex::operator*=(const value_type& s) {
    for(size_type j=0; j!=this->word_size(); ++j) {
        this->word_at(j)*=s;
    }
    return *this;
    for(size_type i=0; i!=this->size()+1; ++i) {
        this->_p[i]*=s;
    }
}

inline
MultiIndex operator+(const MultiIndex& a1, const MultiIndex& a2) {
    MultiIndex r(a1); r+=a2; return r;
}

inline
MultiIndex operator-(const MultiIndex& a1, const MultiIndex& a2) {
    MultiIndex r(a1); r-=a2; return r;
}

inline
MultiIndex operator*(const MultiIndex& a, MultiIndex::value_type s) {
    MultiIndex r(a); r*=s; return r;
}

inline
MultiIndex operator*(MultiIndex::value_type s, const MultiIndex& a) {
    MultiIndex r(a); r*=s; return r;
}






inline
uint 
number(const MultiIndex& i)
{
    uint result=fac(i.degree());
    for(uint k=0; k!=i.size(); ++k) {
        result/=fac(i[k]);
    }
    return result;
}

inline
uint 
fac(const MultiIndex& i)
{
    uint result=1;
    for(uint k=0; k!=i.size(); ++k) {
        result*=fac(i[k]);
    }
    return result;
}

inline
uint 
bin(const MultiIndex& n, const MultiIndex& k)
{
    assert(n.size()==k.size());
    uint result=1;
    for(uint i=0; i!=n.size(); ++i) {
        result*=bin(n[i],k[i]);
    }
    return result;
}



inline 
std::ostream& operator<<(std::ostream& os, const MultiIndex& a) {
    for(uint i=0; i!=a.size(); ++i) { os << (i==0?'(':',') << uint(a[i]); } 
    return os << ')';
}

}

#endif /* ARIADNE_MULTI_INDEX_H */
