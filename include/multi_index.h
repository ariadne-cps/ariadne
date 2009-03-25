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
#include <iostream>

#include "macros.h"
#include "array.h"

namespace Ariadne {

uint32_t fac(uint8_t);
uint32_t bin(uint8_t,uint8_t);

typedef unsigned char uchar;
typedef double Float;

class MultiIndex;
class MultiIndexValueReference;

class MultiIndexBound;

template<class T> class Reference;
template<> class Reference<MultiIndex>;
template<> class Reference<const MultiIndex>;

//! \brief An array of non-negative integers, suitable for storing the
//! powers of a term in some polynomial expansion, ordered by degree.
class MultiIndex {
  public:
    typedef unsigned int size_type;
    typedef unsigned char byte_type;
    typedef unsigned int word_type;

    typedef byte_type value_type;
    typedef MultiIndexValueReference reference;
    typedef const value_type& const_reference;
    //typedef unsigned int word_type;
    //typedef unsigned long long int word_type;
  public:
    /*! Destructor. */
    ~MultiIndex() { _deallocate(this->_p); }
    /*! Construct a multi index with no coefficients. */
    explicit MultiIndex();
    /*! Construct a multi index of degree \a 0 with \a nv variables. */
    explicit MultiIndex(size_type nv);
    /*! Construct a multi index with \a nv variables from the array \a ary. */
    explicit MultiIndex(size_type nv, const int* ary);
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
    const_reference operator[](size_type i) const;
    /*! The number of occurrences of the \a i th variable. */
    reference operator[](size_type i);

    /*! Equality operator. */
    friend bool operator==(const MultiIndex& a1, const MultiIndex& a2); // inline
    /*! Inequality operator. */
    friend bool operator!=(const MultiIndex& a1, const MultiIndex& a2); // inline
    /*! Comparison operator. */
    friend bool operator<(const MultiIndex& a1, const MultiIndex& a2); // inline

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
    /*! Scalar product. */
    friend MultiIndex operator*(value_type s, const MultiIndex& a); // inline

    /*! The position of the element in the array of tensor values. */
    unsigned int position() const;

    /*! The product of the factorials of the indices. */
    unsigned int factorial() const;
    /*! The number of ordered index arrays with each element occurring the number of times specified by the multi index. */
    unsigned int number() const;

    /*! Write to an output stream. */
    friend std::ostream& operator<<(std::ostream&, const MultiIndex&);
  public:
    value_type& at(size_type i) { return _p[i]; }
    const value_type& at(size_type i) const { return _p[i]; }
    value_type* begin() { return _p; }
    const value_type* begin() const { return _p; }
  public:
    size_type word_size() const { return _word_size(_n); }
    word_type& word_at(size_type j) { return reinterpret_cast<word_type*>(_p)[j]; }
    const word_type& word_at(size_type j) const { return reinterpret_cast<const word_type*>(_p)[j]; }
    word_type* word_begin() { return reinterpret_cast<word_type*>(_p); }
    const word_type* word_begin() const { return reinterpret_cast<const word_type*>(_p); }
  public:
    static size_type _word_size(size_type n) {
        return ((n*sizeof(byte_type))/sizeof(word_type)+1); }
    static size_type _element_size(size_type n) {
        return ((n*sizeof(byte_type))/sizeof(word_type)+1+sizeof(double)/sizeof(word_type)); }
    static value_type* _allocate(size_type n) {
        return reinterpret_cast<value_type*>(new word_type[_word_size(n)]); }
    static void _deallocate(value_type* p) {
        delete[] reinterpret_cast<word_type*>(p); }
  private:
    size_type _n;
    value_type* _p;
};


class MultiIndexValueReference {
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::value_type value_type;
    typedef MultiIndex::word_type word_type;
  private:
    size_type _n; value_type* _p; size_type _i; 
  public:
    MultiIndexValueReference(size_type n, value_type* p, size_type i) : _n(n), _p(p), _i(i) { }
    operator const value_type& () { return _p[_i]; }
    MultiIndexValueReference& operator=(const value_type& d) { _p[_n]+=(d-_p[_i]); _p[_i]=d; return *this; }
    MultiIndexValueReference& operator++() { ++_p[_n]; ++_p[_i]; return *this; }
    MultiIndexValueReference& operator--();
};

inline MultiIndexValueReference& MultiIndexValueReference::operator--() {
    if(_p[_i]==0) { 
        ARIADNE_THROW(std::runtime_error,"--MultiIndex[i]"," decrementing zero value at "<<_i<<" in "<<reinterpret_cast<const MultiIndex&>(*this)); }
    --_p[_n]; --_p[_i]; return *this;
}


inline MultiIndex::MultiIndex()
    : _n(0), _p(_allocate(_n))
{
    for(size_type j=0; j!=word_size(); ++j) { word_at(j)=0; }
}

inline MultiIndex::MultiIndex(size_type n)
    : _n(n), _p(_allocate(_n))
{
    for(size_type j=0; j!=word_size(); ++j) { word_at(j)=0; }
}

inline MultiIndex::MultiIndex(size_type n, const int* ary)
    : _n(n), _p(_allocate(_n))
{
    for(size_type j=0; j!=word_size(); ++j) { word_at(j)=0; }
    _p[n]=0; for(size_type i=0; i!=n; ++i) { _p[i]=ary[i]; _p[n]+=ary[i]; }
}

inline MultiIndex::MultiIndex(const MultiIndex& a)
    : _n(a.size()), _p(_allocate(_n))
{
    for(size_type j=0; j!=word_size(); ++j) { this->word_at(j)=a.word_at(j); }
}

inline MultiIndex& MultiIndex::operator=(const MultiIndex& a) {
    if(this!=&a) { this->resize(a.size());
        for(size_type j=0; j!=word_size(); ++j) { this->word_at(j)=a.word_at(j); } }
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
    for(size_type j=0; j!=word_size(); ++j) { this->word_at(j)=0; }
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

inline MultiIndexValueReference MultiIndex::operator[](size_type i) {
    return MultiIndexValueReference(this->_n,this->_p,i);
}


inline bool operator==(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) { return false; }
    for(MultiIndex::size_type i=0; i!=a1.size(); ++i) {
        if(a1.at(i)!=a2.at(i)) { return false; } }
    return true;
    for(MultiIndex::size_type j=0; j!=a1.word_size(); ++j) {
        if(a1.word_at(j)!=a2.word_at(j)) { return false; } }
    return true;
}

inline bool operator!=(const MultiIndex& a1, const MultiIndex& a2) {
    return !(a1==a2);
}

inline bool operator<(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("operator<(MultiIndex,MultiIndex): number of variables must match");
    }

    if(a1.degree()!=a2.degree()) {
        return a1.degree()<a2.degree();
    } else {
        for(MultiIndex::size_type i=0; i!=a1.size(); ++i) {
            if(a1[i]!=a2[i]) {
                return a1[i]>a2[i];
            }
        }
        return false;
        //for(size_type j=0; j!=a1.word_size(); ++j) {
        for(int j=a1.word_size()-1; j!=-1; --j) {
            if(a1.word_at(j)!=a2.word_at(j)) {
                return a1.word_at(j)<a2.word_at(j);
            }
        }
        return false;
    }
}


inline
unsigned int MultiIndex::number() const
{
    unsigned int result=fac(this->degree());
    for(unsigned int k=0; k!=this->size(); ++k) {
        result/=fac((*this)[k]);
    }
    return result;
}

inline
unsigned int MultiIndex::factorial() const
{
    unsigned int result=1;
    for(unsigned int k=0; k!=this->size(); ++k) {
        result*=fac((*this)[k]);
    }
    return result;
}

inline
unsigned int MultiIndex::position() const
{
    unsigned int deg=this->degree()-1;
    unsigned int nvar=this->size();
    unsigned int result=bin(deg+nvar,nvar);
    for(unsigned int k=0; k!=this->size()-1; ++k) {
        --nvar;
        deg-=(*this)[k];
        result+=bin(deg+nvar,nvar);
    }
    return result;
}


inline
MultiIndex& MultiIndex::operator++()
{
    //std::cerr<<"MultiIndex::operator++() with *this="<<*this<<" "<<std::flush;
    assert(_n>0);
    if(_n==1) {
        ++_p[0];
        ++_p[_n];
        return *this;
    }
    if(_p[_n-2]!=0) {
        --_p[_n-2];
        ++_p[_n-1];
        return *this;
    } else {
        value_type li=_p[_n-1];
        _p[_n-1]=0;
        for(size_type k=_n-1; k!=0; --k) {
            if(_p[k-1]!=0) {
                --_p[k-1];
                _p[k]=li+1;
                return *this;
            }
        }
        _p[0]=li+1;
        ++_p[_n];
    }
    return *this;
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
unsigned int
number(const MultiIndex& i)
{
    unsigned int result=fac(i.degree());
    for(unsigned int k=0; k!=i.size(); ++k) {
        result/=fac(i[k]);
    }
    return result;
}

inline
unsigned int
fac(const MultiIndex& i)
{
    unsigned int result=1;
    for(unsigned int k=0; k!=i.size(); ++k) {
        result*=fac(i[k]);
    }
    return result;
}

inline
unsigned int
bin(const MultiIndex& n, const MultiIndex& k)
{
    assert(n.size()==k.size());
    unsigned int result=1;
    for(unsigned int i=0; i!=n.size(); ++i) {
        result*=bin(n[i],k[i]);
    }
    return result;
}


inline
std::ostream& operator<<(std::ostream& os, const MultiIndex& a) {
    //os << "("<<int(a.degree());
    //for(MultiIndex::size_type i=0; i!=a.size(); ++i) { os << (i==0?';':',') << int(a[i]); }
    for(MultiIndex::size_type i=0; i!=a.size(); ++i) { os << (i==0?'(':',') << int(a[i]); }
    return os << ')';
}






/*! \brief A bound on a MultiIndex object, allowing different groups of
 *  variables to have different maximum degrees.
 */
class MultiIndexBound {
  public:
    typedef MultiIndex::size_type size_type;
    MultiIndexBound(size_type as, size_type d);
    MultiIndexBound(const MultiIndex& a);
    MultiIndexBound(size_type as, size_type ng, size_type g1s, size_type g1d, ...);
    size_type size() const { return _groups.size(); }
    friend bool operator<=(const MultiIndex& a, const MultiIndexBound& b);
  private:
    array<size_type> _groups;
    array<size_type> _max_degrees;
};

inline MultiIndexBound::MultiIndexBound(size_type as, size_type d)
    : _groups(as), _max_degrees(1u)
{
    for(size_type i=0; i!=as; ++i) {
        _groups[i]=0u;
    }
    _max_degrees[0]=d;
}


inline MultiIndexBound::MultiIndexBound(const MultiIndex& a)
    : _groups(a.size()), _max_degrees(a.size())
{
    for(size_type i=0; i!=a.size(); ++i) {
        _groups[i]=i;
        _max_degrees[i]=a[i];
    }
}


inline MultiIndexBound::MultiIndexBound(size_type as, size_type ng, size_type g1s, size_type g1d, ...)
    : _groups(as), _max_degrees(ng)
{
    va_list args; va_start(args,g1d);
    size_type k=0;
    for( ; k!=g1s; ++k) {
        _groups[k]=0;
    }
    _max_degrees[0]=g1d;
    for(size_type i=1; i!=ng; ++i) {
        size_type nk=k+va_arg(args,int);
        for( ; k!=nk; ++k) {
            _groups[k]=i;
        }
        _max_degrees[i]=va_arg(args,int);
    }
    va_end(args);
    ARIADNE_ASSERT(k==as);
}

inline bool operator<=(const MultiIndex& a, const MultiIndexBound& b) {
    typedef MultiIndex::size_type size_type;
    array<size_type> degrees(b._max_degrees.size());
    for(size_type j=0; j!=a.size(); ++j) {
        degrees[b._groups[j]]+=a[j];
    }
    for(size_type k=0; k!=degrees.size(); ++k) {
        if(degrees[k]>b._max_degrees[k]) {
            return false;
        }
    }
    return true;
}


} // namespace Ariadne

#endif /* ARIADNE_MULTI_INDEX_H */
