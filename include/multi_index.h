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

/*! \brief \file multi_index.h
 *  \brief An index specifying the degree of differentiation.
 */

#ifndef ARIADNE_MULTI_INDEX_H
#define ARIADNE_MULTI_INDEX_H

#include <cassert>
#include <cstdarg>
#include <iostream>

#include "macros.h"
#include "array.h"
#include "numeric.h"

namespace Ariadne {

uint32_t fac(uint8_t);
uint32_t bin(uint8_t,uint8_t);

typedef unsigned char uchar;

class MultiIndex;
class MultiIndexValueReference;

class MultiIndexBound;

template<class T> class Reference;
template<> class Reference<MultiIndex>;
template<> class Reference<const MultiIndex>;

//struct MultiIndexData { unsigned int _n; unsigned int* _p; };

//! \brief An Array of non-negative integers, suitable for storing the
//! powers of a term in some polynomial expansion, ordered by degree.
//!
//! \par Python interface
//! In the Python interface, multi-indices can be constructed and automatically converted from Python Tuple literals \c (a1,...,am)
//!
//! \b Rationale: The reason why tuples are used for multi-index literals is that they can be used as keys in Python \c dict objects.
class MultiIndex {
  public:
    typedef unsigned int size_type;
    typedef unsigned char byte_type;
    typedef unsigned int word_type;

    typedef byte_type index_type;
    typedef MultiIndexValueReference reference;
    typedef const index_type& const_reference;
    //typedef unsigned int word_type;
    //typedef unsigned long long int word_type;
  public:
    //! \brief Destructor.
    ~MultiIndex();
    //! \brief Construct a multi index with no coefficients.
    explicit MultiIndex();
    //! \brief Construct a multi index of degree \a 0 with \a nv variables.
    explicit MultiIndex(size_type nv);
    //! \brief Construct a multi index with \a nv variables from the Array \a ary.
    explicit MultiIndex(size_type nv, const int* ary);
    //! \brief Construct a multi index with \a nv variables from the Array \a ary.
    explicit MultiIndex(size_type nv, const unsigned char* ary);
    //! \brief Construct a multi index with \a nv variables from variable arguments.
    explicit MultiIndex(size_type nv, int a1, ...);

    //! \brief Copy constructor.
    MultiIndex(const MultiIndex& a);
    //! \brief Copy assignment operator.
    MultiIndex& operator=(const MultiIndex& a);

    //! \brief Construct the zero multi index with \a nv variables.
    static MultiIndex zero(size_type nv);
    //! \brief Construct the unit multi index in variable \a j with \a nv variables.
    static MultiIndex unit(size_type nv, size_type j);
    //! \brief Construct the first multi index of degree \a d with \a nv variables.
    static MultiIndex first(size_type nv, index_type d);

    //! \brief Resize to hold n variables.
    void resize(size_type n);
    //! \brief Assigns values from another index. Precondition: the size of \a a must equal the current size.
    void assign(const MultiIndex& a);
    //! \brief Set all values to zero.
    void clear();
    //! \brief The number of variables.
    size_type size() const;
    //! \brief The degree of the multi-index, equal to the sum of the number of occurrences of the variables.
    index_type degree() const;
     //! \brief The number of variables.
    size_type number_of_variables() const;
    //! \brief The number of occurrences of the \a i th variable.
    index_type get(size_type i) const;
    //! \brief Set the number of occurrences of the \a i th variable to \a n.
    void set(size_type i, index_type n);
    //! \brief The number of occurrences of the \a i th variable.
    const_reference operator[](size_type i) const;
    //! \brief The number of occurrences of the \a i th variable.
    reference operator[](size_type i);
    //! \brief Increment the value of the \a ith element
    void increment(size_type i);
    //! \brief Decrement the value of the \a ith element
    void decrement(size_type i);

    //! \brief Equality operator.
    friend bool operator==(const MultiIndex& a1, const MultiIndex& a2); // inline
    //! \brief Inequality operator.
    friend bool operator!=(const MultiIndex& a1, const MultiIndex& a2); // inline

    //! \brief Increment. No post-increment operator as we sometimes pass MultiIndex by reference.
    MultiIndex& operator++(); // inline
    // No post-increment operator as we sometimes pass MultiIndex by reference.Post increment.
    // MultiIndex operator++(int);
    //! \brief Inplace sum.
    MultiIndex& operator+=(const MultiIndex& a); // inline
    //! \brief Inplace difference.
    MultiIndex& operator-=(const MultiIndex& a); // inline
    //! \brief Inplace scalar product.
    MultiIndex& operator*=(const index_type& a); // inline
    //! \brief Sum.
    friend MultiIndex operator+(const MultiIndex& a1, const MultiIndex& a2); // inline
    //! \brief Difference.
    friend MultiIndex operator-(const MultiIndex& a1, const MultiIndex& a2); // inline
    //! \brief Scalar product.
    friend MultiIndex operator*(const MultiIndex& a, index_type s); // inline
    //! \brief Scalar product.
    friend MultiIndex operator*(index_type s, const MultiIndex& a); // inline

    //! \brief The position of the element in the Array of tensor values.
    unsigned int position() const;

    //! \brief The product of the factorials of the indices.
    unsigned int factorial() const;
    //! \brief The number of ordered index arrays with each element occurring the number of times specified by the multi index.
    unsigned int number() const;

    //! \brief Write to an output stream.
    friend std::ostream& operator<<(std::ostream&, const MultiIndex&);
  public:
    //index_type& at(size_type i) { return _p[i]; }
    const index_type& at(size_type i) const { return reinterpret_cast<const index_type*>(_p)[i]; }
    index_type& at(size_type i) { return reinterpret_cast<index_type*>(_p)[i]; }
    index_type* begin() { return reinterpret_cast<index_type*>(_p); }
    const index_type* begin() const { return reinterpret_cast<const index_type*>(_p); }
  public:
    size_type word_size() const { return _nw; }
    word_type& word_at(size_type j) { return _p[j]; }
    const word_type& word_at(size_type j) const { return _p[j]; }
    word_type* word_begin() { return _p; }
    const word_type* word_begin() const { return _p; }
    word_type* word_end() { return _p+word_size(); }
    const word_type* word_end() const { return _p+word_size(); }
  public:
    static size_type _word_size(size_type n) {
        return ((n*sizeof(byte_type))/sizeof(word_type)+1); }
    static word_type* _allocate_words(size_type nw) {
        word_type* p=new word_type[nw]; return p; }
    static void _deallocate_words(word_type* p) {
        delete[] p; }
  private:
    size_type _n;
    size_type _nw;
    word_type* _p;
};


class MultiIndexValueReference {
    typedef MultiIndex::size_type size_type;
    typedef MultiIndex::index_type index_type;
    typedef MultiIndex::word_type word_type;
  private:
    size_type _n; index_type* _p; size_type _i;
  public:
    MultiIndexValueReference(size_type n, index_type* p, size_type i) : _n(n), _p(p), _i(i) { }
    operator const index_type& () { return _p[_i]; }
    MultiIndexValueReference& operator=(const index_type& d) { _p[_n]+=(d-_p[_i]); _p[_i]=d; return *this; }
    MultiIndexValueReference& operator++() { ++_p[_n]; ++_p[_i]; return *this; }
    MultiIndexValueReference& operator--();
    MultiIndexValueReference& operator+=(int k) { _p[_n]+=k; _p[_i]+=k; return *this; }
    MultiIndexValueReference& operator-=(int k) { _p[_n]-=k; _p[_i]-=k; return *this; }
};

inline MultiIndexValueReference& MultiIndexValueReference::operator--() {
    if(_p[_i]==0) {
        ARIADNE_THROW(std::runtime_error,"--MultiIndex[i]"," decrementing zero value at "<<_i<<" in "<<reinterpret_cast<const MultiIndex&>(*this)); }
    --_p[_n]; --_p[_i]; return *this;
}


inline MultiIndex::~MultiIndex()
{
    _deallocate_words(_p);
}

inline MultiIndex::MultiIndex()
    : _n(0), _nw(_word_size(_n)), _p(_allocate_words(_nw))
{
    std::fill(this->word_begin(),this->word_end(),0);
}

inline MultiIndex::MultiIndex(size_type n)
    : _n(n), _nw(_word_size(_n)), _p(_allocate_words(_nw))
{
    std::fill(this->word_begin(),this->word_end(),0);
}

inline MultiIndex::MultiIndex(size_type n, const unsigned char* ary)
    : _n(n), _nw(_word_size(_n)), _p(_allocate_words(_nw))
{
    for(size_type j=0; j!=word_size(); ++j) { word_at(j)=0; }
    index_type* p=reinterpret_cast<index_type*>(_p);
    p[n]=0; for(size_type i=0; i!=n; ++i) { p[i]=ary[i]; p[n]+=ary[i]; }
}

inline MultiIndex::MultiIndex(size_type n, const int* ary)
    : _n(n), _nw(_word_size(_n)), _p(_allocate_words(_nw))
{
    for(size_type j=0; j!=word_size(); ++j) { word_at(j)=0; }
    index_type* p=reinterpret_cast<index_type*>(_p);
    p[n]=0; for(size_type i=0; i!=n; ++i) { p[i]=ary[i]; p[n]+=ary[i]; }
}

inline MultiIndex::MultiIndex(size_type n, int a0, ...)
    : _n(n), _nw(_word_size(_n)), _p(_allocate_words(_nw))
{
    ARIADNE_ASSERT(n>0);
    index_type* p=reinterpret_cast<index_type*>(_p);
    p[0]=a0; p[_n]=a0;
    va_list args; va_start(args,a0);
    for(size_type i=1; i!=n; ++i) {
        p[i]=va_arg(args,int);
        p[_n]+=p[i];
    }
    va_end(args);
}

inline MultiIndex::MultiIndex(const MultiIndex& a)
    : _n(a.size()), _nw(a._nw), _p(_allocate_words(_nw))
{
    this->assign(a);
}

inline MultiIndex& MultiIndex::operator=(const MultiIndex& a) {
    if(this!=&a) { this->resize(a.size()); this->assign(a); }
    return *this;
}

inline void MultiIndex::assign(const MultiIndex& a) {
    for(size_type j=0; j!=this->word_size(); ++j) { this->word_at(j)=a.word_at(j); }
    //std::copy(a.word_begin(),a.word_end(),this->word_begin());
}

inline MultiIndex MultiIndex::zero(size_type n)
{
    return MultiIndex(n);
}

inline MultiIndex MultiIndex::unit(size_type n, size_type i)
{
    MultiIndex result(n);
    reinterpret_cast<index_type*>(result._p)[i]=1u;
    reinterpret_cast<index_type*>(result._p)[n]=1u;
    return result;
}

inline MultiIndex MultiIndex::first(size_type n, index_type d)
{
    MultiIndex result(n);
    reinterpret_cast<index_type*>(result._p)[0]=d;
    reinterpret_cast<index_type*>(result._p)[n]=d;
    return result;
}

inline void MultiIndex::clear()
{
    std::fill(this->word_begin(),this->word_end(),0);
}

inline void MultiIndex::resize(size_type n) {
    if(this->_n!=n) { _deallocate_words(this->_p); this->_n=n; this->_nw=_word_size(_n); this->_p=_allocate_words(this->_nw); }
}

inline MultiIndex::size_type MultiIndex::size() const {
    return this->_n;
}

inline MultiIndex::index_type MultiIndex::degree() const {
    return reinterpret_cast<const index_type*>(this->_p)[this->_n];
}

inline MultiIndex::size_type MultiIndex::number_of_variables() const {
    return this->_n;
}

inline MultiIndex::index_type MultiIndex::get(size_type i) const {
    assert(i<this->size()); return reinterpret_cast<const index_type*>(this->_p)[i];
}

inline void MultiIndex::set(size_type i, index_type k) {
    assert(i<this->size());
    index_type& ai=reinterpret_cast<index_type*>(this->_p)[i];
    index_type& d=reinterpret_cast<index_type*>(this->_p)[this->_n];
    d+=k; d-=ai; ai=k;
}

inline MultiIndex::index_type const& MultiIndex::operator[](size_type i) const {
    assert(i<this->size()); return reinterpret_cast<const index_type*>(this->_p)[i];
}

inline MultiIndexValueReference MultiIndex::operator[](size_type i) {
    return MultiIndexValueReference(this->_n,reinterpret_cast<index_type*>(this->_p),i);
}

inline void MultiIndex::increment(size_type i) {
    ++reinterpret_cast<index_type*>(this->_p)[i]; ++reinterpret_cast<index_type*>(this->_p)[this->_n];
}

inline void MultiIndex::decrement(size_type i) {
    ARIADNE_ASSERT(reinterpret_cast<index_type*>(this->_p)[i]>0u);
    --reinterpret_cast<index_type*>(this->_p)[i]; --reinterpret_cast<index_type*>(this->_p)[this->_n];
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

inline bool graded_less(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("graded_less(MultiIndex,MultiIndex): number of variables must match");
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

inline bool lexicographic_less(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("lexicographic_less(MultiIndex,MultiIndex): number of variables must match");
    }

    for(MultiIndex::size_type i=0; i!=a1.size(); ++i) {
        if(a1[i]!=a2[i]) {
            return a1[i]<a2[i];
        }
    }
    return false;
}

inline bool reverse_lexicographic_less(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("reverse_lexicographic_less(MultiIndex,MultiIndex): number of variables must match");
    }

    MultiIndex::size_type i=a1.size();
    while(i!=0) {
        --i;
        if(a1[i]!=a2[i]) {
            return a1[i]>a2[i];
        }
    }
    return false;
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

    size_type const n=this->_n;
    index_type* const p=reinterpret_cast<index_type*>(this->_p);

    if(n==1) {
        ++p[0];
        ++p[n];
        return *this;
    }
    if(p[n-2]!=0) {
        --p[n-2];
        ++p[n-1];
        return *this;
    } else {
        index_type li=p[n-1];
        p[n-1]=0;
        for(size_type k=n-1; k!=0; --k) {
            if(p[k-1]!=0) {
                --p[k-1];
                p[k]=li+1;
                return *this;
            }
        }
        p[0]=li+1;
        ++p[n];
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
MultiIndex& MultiIndex::operator*=(const index_type& s) {
    for(size_type j=0; j!=this->word_size(); ++j) {
        this->word_at(j)*=s;
    }
    return *this;
    for(size_type i=0; i!=this->size()+1; ++i) {
        this->_p[i]*=s;
    }
}

inline
void iadd(MultiIndex& r, const MultiIndex& a1, const MultiIndex& a2) {
    for(MultiIndex::size_type j=0; j!=r.word_size(); ++j) {
        r.word_at(j)=a1.word_at(j)+a2.word_at(j);
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
MultiIndex operator*(const MultiIndex& a, MultiIndex::index_type s) {
    MultiIndex r(a); r*=s; return r;
}

inline
MultiIndex operator*(MultiIndex::index_type s, const MultiIndex& a) {
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
    if(a.size()==0) { os << '('; }
    for(MultiIndex::size_type i=0; i!=a.size(); ++i) { os << (i==0?'(':',') << int(a[i]); }
    os<<";"<<int(a.degree());
    return os << ')';
}






//! \brief \brief A bound on a MultiIndex object, allowing different groups of
//!variables to have different maximum degrees.

class MultiIndexBound {
  public:
    typedef MultiIndex::size_type size_type;
    MultiIndexBound(size_type as, size_type d);
    MultiIndexBound(const MultiIndex& a);
    MultiIndexBound(size_type as, size_type ng, size_type g1s, size_type g1d, ...);
    size_type size() const { return _groups.size(); }
    friend bool operator<=(const MultiIndex& a, const MultiIndexBound& b);
  private:
    Array<size_type> _groups;
    Array<size_type> _max_degrees;
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
    Array<size_type> degrees(b._max_degrees.size());
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
