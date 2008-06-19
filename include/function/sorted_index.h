/***************************************************************************
 *            sorted_index.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
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
 
/*! \file sorted_index.h
 *  \brief Sorted indices for polynomials
 */

#ifndef ARIADNE_SORTED_INDEX_H
#define ARIADNE_SORTED_INDEX_H

#include <algorithm>

#include "base/iterator.h"
#include "base/array.h"
#include "base/stlio.h"

#include "numeric/integer.h"

namespace Ariadne {
    
    class MultiIndex;

    class SortedIndex;
    std::ostream& operator<<(std::ostream&, const SortedIndex&);

    /*! \ingroup LinearAlgebra
     *  \brief An index of a tensor object. 
     */
    class SortedIndex {
      friend class MultiIndex;
     public:
      explicit SortedIndex(uint nv);
      explicit SortedIndex(uint nv, uint d);
      SortedIndex(uint nv, uint d, const uint* ary);
      SortedIndex(const SortedIndex& a);
      SortedIndex& operator=(const SortedIndex& a);

      const uint degree() const;
      const uint number_of_variables() const;
      const uint& operator[](const uint& i) const;
     
      bool operator==(const SortedIndex& a) const;
      bool operator!=(const SortedIndex& a) const;
      bool operator<(const SortedIndex& a) const;
     
      SortedIndex& operator++();
      SortedIndex& operator+=(const SortedIndex& a);
      SortedIndex operator+(const SortedIndex& a) const;

      uint position() const;
      
       /*! Increment the the \a i th index, thereby increasing the degree. */
      void increment_index(const uint& i);
      /*! Decrement the the \a i th index, thereby decreasing the degree. */
      void decrement_index(const uint& i);

      void sort();
      friend std::ostream& operator<<(std::ostream&, const SortedIndex&);
     private:
      uint _number_of_variables;
      std::vector<uint> _entries;
    };
    
  
} // namespace Ariadne

namespace Ariadne {

inline 
SortedIndex::SortedIndex(uint nv) : _number_of_variables(nv), _entries(0u) { }

inline 
SortedIndex::SortedIndex(uint nv, uint d) : _number_of_variables(nv), _entries(d) { }

inline 
SortedIndex::SortedIndex(uint nv, uint d, const uint* ary) : _number_of_variables(nv), _entries(ary,ary+d) { this->sort(); }

inline 
SortedIndex::SortedIndex(const SortedIndex& a)
  : _number_of_variables(a._number_of_variables), _entries(a._entries) { }

inline    
SortedIndex& 
SortedIndex::operator=(const SortedIndex& a) 
{ 
  this->_number_of_variables=a._number_of_variables; this->_entries=a._entries; return *this; 
}

inline
const uint 
SortedIndex::degree() const { 
  return this->_entries.size(); 
}

inline
const uint 
SortedIndex::number_of_variables() const { 
  return this->_number_of_variables; 
}
 
inline     
const uint& 
SortedIndex::operator[](const uint& i) const { 
  return this->_entries[i]; 
}
     
inline   
bool 
SortedIndex::operator==(const SortedIndex& a) const { 
  if(this->_number_of_variables!=a._number_of_variables) {
    throw std::runtime_error("operator==(SortedIndex,SortedIndex): number of variables must match");
  }
  return this->_entries==a._entries; 
}

inline  
bool 
SortedIndex::operator!=(const SortedIndex& a) const { 
  return !(*this==a); 
}

inline
SortedIndex &
SortedIndex::operator++()  
{
  uint i=this->degree();
  if(i==0) {
    this->_entries.push_back(0);
    return *this;
  }
  
  while(i!=0) {
    --i;
    ++this->_entries[i];
    if(this->_entries[i]<this->_number_of_variables) {
      break;
    }
  }
  if(this->_entries[i]==this->_number_of_variables) {
    this->_entries.push_back(this->_number_of_variables);
    this->_entries[i]=0;
  }
  uint n=this->_entries[i];
  ++i;
  while(i!=this->degree()) {
    this->_entries[i]=n;
    ++i;
  }
  return *this;
}


inline
SortedIndex&
SortedIndex::operator+=(const SortedIndex& a2) 
{
  SortedIndex& a1=*this;
  if(a1._number_of_variables!=a2._number_of_variables) {
    throw std::runtime_error("operator<(SortedIndex,SortedIndex): number of variables must match");
  }
  uint a1d=a1.degree();
  a1._entries.resize(a1.degree()+a2.degree());
  std::copy(a2._entries.begin(),a2._entries.end(),a1._entries.begin()+a1d);
  std::inplace_merge(a1._entries.begin(),a1._entries.begin()+a1d,a1._entries.end());
  return a1;
}


inline
SortedIndex 
SortedIndex::operator+(const SortedIndex& a2) const 
{
  const SortedIndex& a1=*this;
  if(a1._number_of_variables!=a2._number_of_variables) {
    throw std::runtime_error("operator<(SortedIndex,SortedIndex): number of variables must match");
  }

  SortedIndex r(a1.number_of_variables(),a1.degree()+a2.degree());
  std::merge(a1._entries.begin(),a1._entries.end(),a2._entries.begin(),a2._entries.end(),r._entries.begin());
  return r;
}

inline
bool 
SortedIndex::operator<(const SortedIndex& a) const {
  if(this->_number_of_variables!=a._number_of_variables) {
    throw std::runtime_error("operator<(SortedIndex,SortedIndex): number of variables must match");
  }

  if(this->_entries.size()==a._entries.size()) {
    for(uint i=0; i!=this->_entries.size(); ++i) {
      if(this->_entries[i]!=a._entries[i]) {
        return this->_entries[i]<a._entries[i];
      }
    }
    return false;
  } else {
    return a._entries.size()<a._entries.size();
  }
}

     
inline 
void SortedIndex::increment_index(const uint& i) { 
  this->_entries.push_back(i); 
  std::inplace_merge(this->_entries.begin(),this->_entries.end()-1,this->_entries.end());
}

inline
void SortedIndex::decrement_index(const uint& i) 
{ 
  std::vector<uint>::iterator pos=std::lower_bound(this->_entries.begin(),this->_entries.end(),i);
  if(pos==this->_entries.end()) {
    throw std::runtime_error("SortedIndex::decrement_index(uint i): the number of occurence of the index must be positive"); 
  } else {
    this->_entries.erase(pos);
  }
}

inline 
void 
SortedIndex::sort() { 
  std::sort(this->_entries.begin(),this->_entries.end()); 
}
    
inline
uint 
SortedIndex::position() const
{
  const uint& nv=this->number_of_variables();
  const uint d=this->degree();
  const uint* p=&this->_entries[0];
  long result=bin(nv+d,d);
  //std::cerr << "\n"<<*this<<"\n";
  //std::cerr << "k=-" << " (" << (nv+d) << "," << nv << ")=" << bin(nv+d,nv) << " r=" << result << std::endl;
  for(uint k=0; k!=d; ++k) {
    int m=d-k;
    int n=nv-p[k]-2;
    //std::cerr << "k=" << k << " p[k]="<<p[k]<<" ("<< m+n << "," << m << ")=" << bin(m+n,m) << " r=" << result << std::endl;
    result-=bin(m+n,m);
  }
  result-=1;
  //std::cerr << "position(" << (*this) << ")=" << result << " ";
  return result;
}
      

      
inline 
std::ostream& 
operator<<(std::ostream& os, const SortedIndex& a) {
  return os << a._entries;
}
  
    
} // namespace Ariadne

#endif /* ARIADNE_SORTED_INDEX_H */
