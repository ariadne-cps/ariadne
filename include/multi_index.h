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
#include "array.h"

namespace Ariadne {

uint fac(uint);
uint bin(uint,uint);

typedef unsigned char uchar;

class MultiIndex {
 public:
  /*! Construct a multi index of degree \a 0 with \a nv variables. */
  explicit MultiIndex(uint nv);
  /*! Construct a multi index with \a nv variables from the array \a ary. */
  explicit MultiIndex(uint nv, const uint* ary);
  /*! Copy constructor. */
  MultiIndex(const MultiIndex& a);
  /*! Copy assignment operator. */
  MultiIndex& operator=(const MultiIndex& a);
  
  /*! Construct the zero multi index with \a nv variables. */
  static MultiIndex zero(uint nv);
  /*! Construct the unit multi index in variable \a j with \a nv variables. */
  static MultiIndex unit(uint nv, uint j);
  /*! Construct the first multi index of degree \a d with \a nv variables. */
  static MultiIndex first(uint nv, uint d);

  /*! The number of variables. */
  const uint size() const;
  /*! The degree of the multi-index, equal to the sum of the number of occurrences of the variables. */
  const uint degree() const;
  /*! The number of variables. */
  const uint number_of_variables() const;
  /*! The number of occurrences of the \a i th variable. */
  const uint& operator[](const uint& i) const;
  
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
  /*! Sum. */
  MultiIndex operator+(const MultiIndex& a) const; // inline
  /*! Product. */
  MultiIndex& operator*=(const uint& a); // inline
  MultiIndex operator*(const uint& a) const; // inline
  
  /*! The position of the element in the array of tensor values. */
  uint position() const;
  
  /*! The product of the factorials of the indices. */
  uint factorial() const;
  /*! The number of ordered index arrays with each element occurring the number of times specified by the multi index. */
  uint number() const;
  
  /*! Set the value of the \a i th index to \a j. */
  void set_index(const uint& i, const uint j);
  /*! Increment the the \a i th index, thereby increasing the degree. */
  void increment_index(const uint& i);
  /*! Decrement the the \a i th index, thereby decreasing the degree. */
  void decrement_index(const uint& i);
  
  /*! Set the value of the \a i th index to \a j. */
  void set(const uint& i, const uint j);
  friend std::ostream& operator<<(std::ostream&, const MultiIndex&);
 private:
  uint _degree;
  array<uint> _occurrences;
};
      


    
inline MultiIndex::MultiIndex(uint nv)
  : _degree(0), _occurrences(nv,0) 
{
}

inline MultiIndex::MultiIndex(uint nv, const uint* ary)
  : _degree(0), _occurrences(nv) 
{
  for(uint i=0; i!=nv; ++i) { this->_occurrences[i]=ary[i]; this->_degree+=ary[i]; }
}

inline MultiIndex::MultiIndex(const MultiIndex& a)
  : _degree(a._degree), _occurrences(a._occurrences) 
{
}

inline MultiIndex& MultiIndex::operator=(const MultiIndex& a) { 
  this->_degree=a._degree; this->_occurrences=a._occurrences; return *this; }


inline MultiIndex MultiIndex::zero(uint nv)
{
  return MultiIndex(nv);
}

inline MultiIndex MultiIndex::unit(uint nv, uint j)
{
  MultiIndex result(nv);
  result._occurrences[j]=1u;
  result._degree=1u;
  return result;
}

inline MultiIndex MultiIndex::first(uint nv, uint d)
{
  MultiIndex result(nv);
  result._occurrences[0]=d;
  result._degree=d;
  return result;
}

inline const uint MultiIndex::size() const { 
  return this->_occurrences.size(); 
}

inline const uint MultiIndex::degree() const { 
  return this->_degree; 
}

inline const uint MultiIndex::number_of_variables() const { 
  return this->_occurrences.size(); 
}

inline const uint& MultiIndex::operator[](const uint& i) const {
  assert(i<this->number_of_variables()); return this->_occurrences[i]; 
}


inline bool MultiIndex::operator==(const MultiIndex& a) const { 
  return this->_occurrences==a._occurrences; 
}

inline bool MultiIndex::operator!=(const MultiIndex& a) const { 
  return !(*this==a); 
}

inline void MultiIndex::set_index(const uint& i, const uint j) { 
  this->_degree+=(j-this->_occurrences[i]); this->_occurrences[i]=j; 
}

inline void MultiIndex::increment_index(const uint& i) { 
  ++this->_occurrences[i]; ++this->_degree; 
}

inline void MultiIndex::decrement_index(const uint& i) { 
  if(!(this->_occurrences[i]>0)) { 
    throw std::runtime_error("MultiIndex::decrement_index: the number of occurence of the index must be positive"); 
  }
  --this->_occurrences[i]; 
  --this->_degree; 
}

inline void MultiIndex::set(const uint& i, const uint j) { 
  this->set_index(i,j); 
}





inline
uint MultiIndex::number() const
{
  uint result=fac(this->degree());
  for(uint k=0; k!=this->number_of_variables(); ++k) {
    result/=fac((*this)[k]);
  }
  return result;
}

inline
uint MultiIndex::factorial() const
{
  uint result=1;
  for(uint k=0; k!=this->number_of_variables(); ++k) {
    result*=fac((*this)[k]);
  }
  return result;
}

inline
uint MultiIndex::position() const
{
  uint deg=this->degree()-1;
  uint nvar=this->number_of_variables();
  uint result=bin(deg+nvar,nvar);
  for(uint k=0; k!=this->number_of_variables()-1; ++k) {
    --nvar;
    deg-=(*this)[k];
    result+=bin(deg+nvar,nvar);
  }
  return result;
}

inline
bool MultiIndex::operator<(const MultiIndex& a2) const {
  const MultiIndex& a1=*this;
  
  if(a1.number_of_variables()!=a2.number_of_variables()) {
    throw std::runtime_error("operator<(SortedIndex,SortedIndex): number of variables must match");
  }
  
  if(a1.degree()!=a2.degree()) {
    return a1.degree()<a2.degree();
  } else {
    for(uint i=0; i!=a1.number_of_variables(); ++i) {
      if(a1[i]!=a2[i]) { 
        return a1[i]<a2[i];
      }
    }
    return false;
  }
}


inline
MultiIndex& MultiIndex::operator++() 
{
  MultiIndex& a=*this;
  const uint& nv=a.number_of_variables();
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
MultiIndex& MultiIndex::operator+=(const MultiIndex& a2) {
  for(uint i=0; i!=this->number_of_variables(); ++i) {
    this->set_index(i,(*this)[i]+a2[i]);
  }
  return *this;
}

inline
MultiIndex& MultiIndex::operator-=(const MultiIndex& a2) {
  for(uint i=0; i!=this->number_of_variables(); ++i) {
    this->set_index(i,(*this)[i]-a2[i]);
  }
  return *this;
}

inline
MultiIndex& MultiIndex::operator*=(const uint& n) {
  for(uint i=0; i!=this->number_of_variables(); ++i) {
    this->set_index(i,(*this)[i]*n);
  }
  return *this;
}

inline
MultiIndex MultiIndex::operator+(const MultiIndex& a2) const {
  return MultiIndex(*this)+=a2;
}

inline
MultiIndex operator-(const MultiIndex& a, const MultiIndex& b) {
  MultiIndex c(a); c-=b; return c;
}

inline
MultiIndex operator*(const MultiIndex& a, const uint& n) {
  MultiIndex c(a); c*=n; return c;
}

inline
MultiIndex operator*(const MultiIndex& a, const int& n) {
  MultiIndex c(a); c*=n; return c;
}


inline
uint 
number(const MultiIndex& i)
{
  uint result=fac(i.degree());
  for(uint k=0; k!=i.number_of_variables(); ++k) {
    result/=fac(i[k]);
  }
  return result;
}

inline
uint 
fac(const MultiIndex& i)
{
  uint result=1;
  for(uint k=0; k!=i.number_of_variables(); ++k) {
    result*=fac(i[k]);
  }
  return result;
}

inline
uint 
bin(const MultiIndex& n, const MultiIndex& k)
{
  assert(n.number_of_variables()==k.number_of_variables());
  uint result=1;
  for(uint i=0; i!=n.number_of_variables(); ++i) {
    result*=bin(n[i],k[i]);
  }
  return result;
}




// New code



/*
MultiIndex& operator+=(MultiIndex& a, const MultiIndex& b) {
  assert(a._as==b._as);
  a._d+=b._d;
  for(uint i=0; i!=a._as; ++i) {
    a._o[i]+=b._o[i];
  }
  return a;
}

MultiIndex& operator-=(MultiIndex& a, const MultiIndex& b) {
  assert(a._as==b._as);
  a._d-=b._d;
  for(uint i=0; i!=a._as; ++i) {
    assert(a._o[i]>=b._o[i]);
    a._o[i]-=b._o[i];
  }
  return a;
}

MultiIndex& operator*=(MultiIndex& a, const uint& n) {
  a._d*=n;
  for(uint i=0; i!=a._as; ++i) {
    a._o[i]*=n;
  }
  return a;
}

MultiIndex operator+(const MultiIndex& a, const MultiIndex& b) {
  MultiIndex c(a); c+=b; return c;
}

MultiIndex operator*(const MultiIndex& a, const uint& n) {
  MultiIndex c(a); c*=n; return c;
}

bool operator==(const MultiIndex& a, const MultiIndex& b) {
  if(a.argument_size()!=b.argument_size()) { return false; }
  for(uint i=0; i!=a.argument_size(); ++i) {
    if(a[i]!=b[i]) { return false; } }
  return true;
}
  
bool operator<(const MultiIndex& a, const MultiIndex& b) {
  if(a.degree()==b.degree()) {
    for(uint i=0; i!=a.argument_size(); ++i) {
      if(a[i]!=b[i]) { return a[i]>b[i]; } 
    }
    return false;
  } else {
    return a.degree()<b.degree();
  }
}



class MultiIndexIterator
  : public boost::iterator_facade<MultiIndexIterator,
                                  MultiIndex,
                                  boost::forward_traversal_tag,
                                  const MultiIndex&,
                                  const MultiIndex*>
{
  MultiIndex _index;
 public:
  MultiIndexIterator(const MultiIndex& i) : _index(i) { }
  MultiIndexIterator(const uint& nv, const uint& d)
    : _index(MultiIndex(nv)) { this->_index.set_index(0,d); }
  
  void increment() {
    const uint nv=this->_index.number_of_variables();
    if(this->_index[nv-2]!=0) {
      this->_index.decrement_index(nv-2);
      this->_index.increment_index(nv-1);
      return;
    } else {
      uint li=this->_index[nv-1];
      this->_index.set_index(nv-1,0);
      for(uint k=nv-1; k!=0; --k) {
        if(this->_index[k-1]!=0) {
          this->_index.decrement_index(k-1);
          this->_index.set_index(k,li+1);
          return;
        }
      }
      this->_index.set_index(0,li+1);
    }
  }
  
  const MultiIndex& dereference() const { return this->_index; }
  
  bool equal(const MultiIndexIterator& other) const { return this->_index==other._index; }
};



*/

inline 
std::ostream& operator<<(std::ostream& os, const MultiIndex& a) {
  return os << a._occurrences;
}

}

#endif /* ARIADNE_MULTI_INDEX_H */
