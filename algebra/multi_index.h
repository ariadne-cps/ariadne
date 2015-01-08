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
#include <initializer_list>
#include <iostream>

#include "utility/macros.h"
#include "utility/array.h"
#include "numeric/numeric.h"

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
    typedef unsigned int SizeType;
    typedef unsigned char ByteType;
    typedef unsigned int WordType;

    typedef ByteType IndexType;
    typedef MultiIndexValueReference Reference;
    typedef const IndexType& ConstReference;
    //typedef unsigned int WordType;
    //typedef unsigned long long int WordType;
  public:
    //! \brief Destructor.
    ~MultiIndex();
    //! \brief Construct a multi index with no coefficients.
    explicit MultiIndex();
    //! \brief Construct a multi index of degree \a 0 with \a nv variables.
    explicit MultiIndex(SizeType nv);
    //! \brief Construct a multi index with \a nv variables from the Array \a ary.
    explicit MultiIndex(SizeType nv, const Int* ary);
    //! \brief Construct a multi index with \a nv variables from the Array \a ary.
    explicit MultiIndex(SizeType nv, const unsigned char* ary);
    //! \brief Construct a multi index with from an initializer list.
    MultiIndex(InitializerList<Int> lst);

    //! \brief Copy constructor.
    MultiIndex(const MultiIndex& a);
    //! \brief Copy assignment operator.
    MultiIndex& operator=(const MultiIndex& a);

    //! \brief Construct the zero multi index with \a nv variables.
    static MultiIndex zero(SizeType nv);
    //! \brief Construct the unit multi index in variable \a j with \a nv variables.
    static MultiIndex unit(SizeType nv, SizeType j);
    //! \brief Construct the first multi index of degree \a d with \a nv variables.
    static MultiIndex first(SizeType nv, IndexType d);

    //! \brief Resize to hold n variables.
    Void resize(SizeType n);
    //! \brief Assigns values from another index. Precondition: the size of \a a must equal the current size.
    Void assign(const MultiIndex& a);
    //! \brief Set all values to zero.
    Void clear();
    //! \brief The number of variables.
    SizeType size() const;
    //! \brief The degree of the multi-index, equal to the sum of the number of occurrences of the variables.
    IndexType degree() const;
     //! \brief The number of variables.
    SizeType number_of_variables() const;
    //! \brief The number of occurrences of the \a i th variable.
    IndexType get(SizeType i) const;
    //! \brief Set the number of occurrences of the \a i th variable to \a n.
    Void set(SizeType i, IndexType n);
    //! \brief The number of occurrences of the \a i th variable.
    const IndexType& operator[](SizeType i) const;
    //! \brief The number of occurrences of the \a i th variable.
    MultiIndexValueReference operator[](SizeType i);
    //! \brief Increment the value of the \a ith element
    Void increment(SizeType i);
    //! \brief Decrement the value of the \a ith element
    Void decrement(SizeType i);

    //! \brief Equality operator.
    friend Bool operator==(const MultiIndex& a1, const MultiIndex& a2); // inline
    //! \brief Inequality operator.
    friend Bool operator!=(const MultiIndex& a1, const MultiIndex& a2); // inline

    //! \brief Increment. No post-increment operator as we sometimes pass MultiIndex by reference.
    MultiIndex& operator++(); // inline
    // No post-increment operator as we sometimes pass MultiIndex by reference.Post increment.
    // MultiIndex operator++(Int);
    //! \brief Inplace sum.
    MultiIndex& operator+=(const MultiIndex& a); // inline
    //! \brief Inplace difference.
    MultiIndex& operator-=(const MultiIndex& a); // inline
    //! \brief Inplace scalar product.
    MultiIndex& operator*=(const IndexType& a); // inline
    //! \brief Sum.
    friend MultiIndex operator+(const MultiIndex& a1, const MultiIndex& a2); // inline
    //! \brief Difference.
    friend MultiIndex operator-(const MultiIndex& a1, const MultiIndex& a2); // inline
    //! \brief Scalar product.
    friend MultiIndex operator*(const MultiIndex& a, IndexType s); // inline
    //! \brief Scalar product.
    friend MultiIndex operator*(IndexType s, const MultiIndex& a); // inline

    //! \brief The position of the element in the Array of tensor values.
    unsigned int position() const;

    //! \brief The product of the factorials of the indices.
    unsigned int factorial() const;
    //! \brief The number of ordered index arrays with each element occurring the number of times specified by the multi index.
    unsigned int number() const;

    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream&, const MultiIndex&);
  public:
    //IndexType& at(SizeType i) { return _p[i]; }
    const IndexType& at(SizeType i) const { return reinterpret_cast<const IndexType*>(_p)[i]; }
    IndexType& at(SizeType i) { return reinterpret_cast<IndexType*>(_p)[i]; }
    IndexType* begin() { return reinterpret_cast<IndexType*>(_p); }
    const IndexType* begin() const { return reinterpret_cast<const IndexType*>(_p); }
  public:
    SizeType word_size() const { return _nw; }
    WordType& word_at(SizeType j) { return _p[j]; }
    const WordType& word_at(SizeType j) const { return _p[j]; }
    WordType* word_begin() { return _p; }
    const WordType* word_begin() const { return _p; }
    WordType* word_end() { return _p+word_size(); }
    const WordType* word_end() const { return _p+word_size(); }
  public:
    static SizeType _word_size(SizeType n) {
        return ((n*sizeof(ByteType))/sizeof(WordType)+1); }
    static WordType* _allocate_words(SizeType nw) {
        WordType* p=new WordType[nw]; return p; }
    static Void _deallocate_words(WordType* p) {
        delete[] p; }
  private:
    SizeType _n;
    SizeType _nw;
    WordType* _p;
};


class MultiIndexValueReference {
    typedef MultiIndex::SizeType SizeType;
    typedef MultiIndex::IndexType IndexType;
    typedef MultiIndex::WordType WordType;
  private:
    SizeType _n; IndexType* _p; SizeType _i;
  public:
    MultiIndexValueReference(SizeType n, IndexType* p, SizeType i) : _n(n), _p(p), _i(i) { }
    operator const IndexType& () { return _p[_i]; }
    MultiIndexValueReference& operator=(const IndexType& d) { _p[_n]+=(d-_p[_i]); _p[_i]=d; return *this; }
    MultiIndexValueReference& operator++() { ++_p[_n]; ++_p[_i]; return *this; }
    MultiIndexValueReference& operator--();
    MultiIndexValueReference& operator+=(Int k) { _p[_n]+=k; _p[_i]+=k; return *this; }
    MultiIndexValueReference& operator-=(Int k) { _p[_n]-=k; _p[_i]-=k; return *this; }
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

inline MultiIndex::MultiIndex(SizeType n)
    : _n(n), _nw(_word_size(_n)), _p(_allocate_words(_nw))
{
    std::fill(this->word_begin(),this->word_end(),0);
}

inline MultiIndex::MultiIndex(SizeType n, const unsigned char* ary)
    : _n(n), _nw(_word_size(_n)), _p(_allocate_words(_nw))
{
    for(SizeType j=0; j!=word_size(); ++j) { word_at(j)=0; }
    IndexType* p=reinterpret_cast<IndexType*>(_p);
    p[n]=0; for(SizeType i=0; i!=n; ++i) { p[i]=ary[i]; p[n]+=ary[i]; }
}

inline MultiIndex::MultiIndex(SizeType n, const Int* ary)
    : _n(n), _nw(_word_size(_n)), _p(_allocate_words(_nw))
{
    for(SizeType j=0; j!=word_size(); ++j) { word_at(j)=0; }
    IndexType* p=reinterpret_cast<IndexType*>(_p);
    p[n]=0; for(SizeType i=0; i!=n; ++i) { p[i]=ary[i]; p[n]+=ary[i]; }
}

inline MultiIndex::MultiIndex(InitializerList<Int> lst)
    : _n(lst.size()), _nw(_word_size(_n)), _p(_allocate_words(_nw))
{
    IndexType* ptr=reinterpret_cast<IndexType*>(_p);
    InitializerList<Int>::const_iterator iter=lst.begin();
    IndexType* deg_ptr=ptr+_n;
    *deg_ptr=0;
    while(iter!=lst.end()) {
        *ptr=*iter;
        *deg_ptr+=*ptr;
        ++ptr; ++iter;
    }
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

inline Void MultiIndex::assign(const MultiIndex& a) {
    for(SizeType j=0; j!=this->word_size(); ++j) { this->word_at(j)=a.word_at(j); }
    //std::copy(a.word_begin(),a.word_end(),this->word_begin());
}

inline MultiIndex MultiIndex::zero(SizeType n)
{
    return MultiIndex(n);
}

inline MultiIndex MultiIndex::unit(SizeType n, SizeType i)
{
    MultiIndex result(n);
    reinterpret_cast<IndexType*>(result._p)[i]=1u;
    reinterpret_cast<IndexType*>(result._p)[n]=1u;
    return result;
}

inline MultiIndex MultiIndex::first(SizeType n, IndexType d)
{
    MultiIndex result(n);
    reinterpret_cast<IndexType*>(result._p)[0]=d;
    reinterpret_cast<IndexType*>(result._p)[n]=d;
    return result;
}

inline Void MultiIndex::clear()
{
    std::fill(this->word_begin(),this->word_end(),0);
}

inline Void MultiIndex::resize(SizeType n) {
    if(this->_n!=n) { _deallocate_words(this->_p); this->_n=n; this->_nw=_word_size(_n); this->_p=_allocate_words(this->_nw); }
}

inline MultiIndex::SizeType MultiIndex::size() const {
    return this->_n;
}

inline MultiIndex::IndexType MultiIndex::degree() const {
    return reinterpret_cast<const IndexType*>(this->_p)[this->_n];
}

inline MultiIndex::SizeType MultiIndex::number_of_variables() const {
    return this->_n;
}

inline MultiIndex::IndexType MultiIndex::get(SizeType i) const {
    assert(i<this->size()); return reinterpret_cast<const IndexType*>(this->_p)[i];
}

inline Void MultiIndex::set(SizeType i, IndexType k) {
    assert(i<this->size());
    IndexType& ai=reinterpret_cast<IndexType*>(this->_p)[i];
    IndexType& d=reinterpret_cast<IndexType*>(this->_p)[this->_n];
    d+=k; d-=ai; ai=k;
}

inline MultiIndex::IndexType const& MultiIndex::operator[](SizeType i) const {
    assert(i<this->size()); return reinterpret_cast<const IndexType*>(this->_p)[i];
}

inline MultiIndexValueReference MultiIndex::operator[](SizeType i) {
    return MultiIndexValueReference(this->_n,reinterpret_cast<IndexType*>(this->_p),i);
}

inline Void MultiIndex::increment(SizeType i) {
    ++reinterpret_cast<IndexType*>(this->_p)[i]; ++reinterpret_cast<IndexType*>(this->_p)[this->_n];
}

inline Void MultiIndex::decrement(SizeType i) {
    ARIADNE_ASSERT(reinterpret_cast<IndexType*>(this->_p)[i]>0u);
    --reinterpret_cast<IndexType*>(this->_p)[i]; --reinterpret_cast<IndexType*>(this->_p)[this->_n];
}


inline Bool operator==(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) { return false; }
    for(MultiIndex::SizeType i=0; i!=a1.size(); ++i) {
        if(a1.at(i)!=a2.at(i)) { return false; } }
    return true;
    for(MultiIndex::SizeType j=0; j!=a1.word_size(); ++j) {
        if(a1.word_at(j)!=a2.word_at(j)) { return false; } }
    return true;
}

inline Bool operator!=(const MultiIndex& a1, const MultiIndex& a2) {
    return !(a1==a2);
}

inline Bool graded_less(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("graded_less(MultiIndex,MultiIndex): number of variables must match");
    }

    if(a1.degree()!=a2.degree()) {
        return a1.degree()<a2.degree();
    } else {
        for(MultiIndex::SizeType i=0; i!=a1.size(); ++i) {
            if(a1[i]!=a2[i]) {
                return a1[i]>a2[i];
            }
        }
        return false;
        //for(SizeType j=0; j!=a1.word_size(); ++j) {
        for(Int j=a1.word_size()-1; j!=-1; --j) {
            if(a1.word_at(j)!=a2.word_at(j)) {
                return a1.word_at(j)<a2.word_at(j);
            }
        }
        return false;
    }
}

inline Bool lexicographic_less(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("lexicographic_less(MultiIndex,MultiIndex): number of variables must match");
    }

    for(MultiIndex::SizeType i=0; i!=a1.size(); ++i) {
        if(a1[i]!=a2[i]) {
            return a1[i]<a2[i];
        }
    }
    return false;
}

inline Bool reverse_lexicographic_less(const MultiIndex& a1, const MultiIndex& a2) {
    if(a1.size()!=a2.size()) {
        throw std::runtime_error("reverse_lexicographic_less(MultiIndex,MultiIndex): number of variables must match");
    }

    MultiIndex::SizeType i=a1.size();
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

    SizeType const n=this->_n;
    IndexType* const p=reinterpret_cast<IndexType*>(this->_p);

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
        IndexType li=p[n-1];
        p[n-1]=0;
        for(SizeType k=n-1; k!=0; --k) {
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
    for(SizeType j=0; j!=this->word_size(); ++j) {
        this->word_at(j)+=a.word_at(j);
    }
    return *this;
    for(SizeType i=0; i!=this->size(); ++i) {
        this->_p[i]+=a._p[i];
    }
}

inline
MultiIndex& MultiIndex::operator-=(const MultiIndex& a) {
    for(SizeType i=0; i!=this->size()+1; ++i) {
        ARIADNE_ASSERT(this->_p[i]>=a._p[i]);
        this->_p[i]-=a._p[i];
    }
    return *this;
}

inline
MultiIndex& MultiIndex::operator*=(const IndexType& s) {
    for(SizeType j=0; j!=this->word_size(); ++j) {
        this->word_at(j)*=s;
    }
    return *this;
    for(SizeType i=0; i!=this->size()+1; ++i) {
        this->_p[i]*=s;
    }
}

inline
Void iadd(MultiIndex& r, const MultiIndex& a1, const MultiIndex& a2) {
    for(MultiIndex::SizeType j=0; j!=r.word_size(); ++j) {
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
MultiIndex operator*(const MultiIndex& a, MultiIndex::IndexType s) {
    MultiIndex r(a); r*=s; return r;
}

inline
MultiIndex operator*(MultiIndex::IndexType s, const MultiIndex& a) {
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
OutputStream& operator<<(OutputStream& os, const MultiIndex& a) {
    //os << "("<<Int(a.degree());
    //for(MultiIndex::SizeType i=0; i!=a.size(); ++i) { os << (i==0?';':',') << Int(a[i]); }
    if(a.size()==0) { os << '('; }
    for(MultiIndex::SizeType i=0; i!=a.size(); ++i) { os << (i==0?'(':',') << Int(a[i]); }
    os<<";"<<Int(a.degree());
    return os << ')';
}






//! \brief \brief A bound on a MultiIndex object, allowing different groups of
//!variables to have different maximum degrees.

class MultiIndexBound {
  public:
    typedef MultiIndex::SizeType SizeType;
    MultiIndexBound(SizeType as, SizeType d);
    MultiIndexBound(const MultiIndex& a);
    SizeType size() const { return _groups.size(); }
    friend Bool operator<=(const MultiIndex& a, const MultiIndexBound& b);
  private:
    Array<SizeType> _groups;
    Array<SizeType> _max_degrees;
};

inline MultiIndexBound::MultiIndexBound(SizeType as, SizeType d)
    : _groups(as), _max_degrees(1u)
{
    for(SizeType i=0; i!=as; ++i) {
        _groups[i]=0u;
    }
    _max_degrees[0]=d;
}


inline MultiIndexBound::MultiIndexBound(const MultiIndex& a)
    : _groups(a.size()), _max_degrees(a.size())
{
    for(SizeType i=0; i!=a.size(); ++i) {
        _groups[i]=i;
        _max_degrees[i]=a[i];
    }
}


inline Bool operator<=(const MultiIndex& a, const MultiIndexBound& b) {
    typedef MultiIndex::SizeType SizeType;
    Array<SizeType> degrees(b._max_degrees.size());
    for(SizeType j=0; j!=a.size(); ++j) {
        degrees[b._groups[j]]+=a[j];
    }
    for(SizeType k=0; k!=degrees.size(); ++k) {
        if(degrees[k]>b._max_degrees[k]) {
            return false;
        }
    }
    return true;
}


} // namespace Ariadne

#endif /* ARIADNE_MULTI_INDEX_H */
