/***************************************************************************
 *            multi_index.h
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
 
/*! \file multi_index.h
 *  \brief Multi-indices of symmetric tensors.
 */

#ifndef ARIADNE_MULTI_INDEX_H
#define ARIADNE_MULTI_INDEX_H

#include "../base/iterator.h"
#include "../base/array.h"
#include "../base/stlio.h"

namespace Ariadne {
  namespace LinearAlgebra {
    
    template<class C>
    class ContainerElementReference {
     public:
      ContainerElementReference(C& c, const typename C::size_type& i) : _c(c), _i(i) { }
      operator const typename C::value_type& () const { return _c[_i]; }
      void operator=(const typename C::value_type& x) { _c.set(_i,x); }
     private:
      C& _c; const typename C::size_type& _i;
    };
    
   /*! \ingroup LinearAlgebra
     *  \brief An index of a tensor object. 
     */
    class Index {
     public:
      explicit Index(size_type nv) : _number_of_variables(nv), _entries(0u) { }
      Index(const Index& a) : _number_of_variables(a._number_of_variables), _entries(a._entries) { }
      Index& operator=(const Index& a) { 
        this->_number_of_variables=a._number_of_variables; this->_entries=a._entries; return *this; }

      const size_type degree() const { return this->_entries.size(); }
      const size_type number_of_variables() const { return this->_number_of_variables; };
      const size_type& operator[](const size_type& i) const { return this->_entries[i]; }
     
      bool operator==(const Index& a) const { 
        return this->_entries==a._entries && this->_number_of_variables==a._number_of_variables; }
      bool operator!=(const Index& a) const { return !(*this==a); }
      bool operator<(const Index& a) const; 
     
      size_type position() const;
      
      void push_back(const size_type& i) { 
        if(!(i<this->_number_of_variables)) { throw InvalidIndex(__PRETTY_FUNCTION__); }
        this->_entries.push_back(i); }
      void pop_back() { this->_entries.pop_back(); }
      
      friend std::ostream& operator<<(std::ostream&, const Index&);
     private:
      size_type _number_of_variables;
      std::vector<size_type> _entries;
    };
    
    inline
    size_type Index::position() const
    {
      size_type result=this->_entries[0];
      for(size_type k=1; k!=this->_entries.size(); ++k) {
        result*=this->_number_of_variables;
        result+=this->_entries[k];
      }
      //std::cerr << "position(" << (*this) << ")=" << result << " ";
      return result;
    }
      
   /*! \ingroup LinearAlgebra
     *  \brief An index of a symmetric object. 
     */
    class MultiIndex {
     public:
      typedef Ariadne::Base::size_type size_type;
      typedef size_type value_type;
     
      /*! Construct a multi index of degree \a 0 with \a nv variables. */
      explicit MultiIndex(size_type nv) : _degree(0u), _occurrences(nv,0u) { }
      /*! Construct a multi index from the tensor index \a a. */
      explicit MultiIndex(const Index& a) : _degree(0u), _occurrences(a.number_of_variables(),0u) { 
        for(size_type i=0; i!=a.degree(); ++i) { this->increment_index(a[i]); } }
      /*! Construct a multi index with \a nv variables from the ordered list of indices \a a. */
      explicit MultiIndex(const size_type nv, const IndexArray& a) : _degree(0u), _occurrences(nv,0u) { 
        for(size_type i=0; i!=a.size(); ++i) { this->increment_index(a[i]); } 
      }
      /*! Copy constructor. */
      MultiIndex(const MultiIndex& a) : _degree(a._degree), _occurrences(a._occurrences) { }
      /*! Copy assignment operator. */
      MultiIndex& operator=(const MultiIndex& a) { 
        this->_degree=a._degree; this->_occurrences=a._occurrences; return *this; }

      /*! The degree of the multi-index, equal to the sum of the number of occurrences of the variables. */
      const size_type degree() const { return this->_degree; }
      /*! The number of variables. */
      const size_type number_of_variables() const { return this->_occurrences.size(); };
      /*! The number of occurrences of the \a i th variable. */
      const size_type& operator[](const size_type& i) const { return this->_occurrences[i]; }
      /*! The number of occurrences of the \a i th variable. */
      //ContainerElementReference<MultiIndex> operator[](const size_type& i) { 
      //  return ContainerElementReference<MultiIndex>(*this,i); }
     
      /*! Equality operator. */
      bool operator==(const MultiIndex& a) const { return this->_occurrences==a._occurrences; }
      /*! Inequality operator. */
      bool operator!=(const MultiIndex& a) const { return !(*this==a); }
      /*! Comparison operator. */
      bool operator<(const MultiIndex& a) const; // inline
      /*! Sum. */
      MultiIndex operator+(const MultiIndex& a) const; // inline
     
      /*! Convert to a normal tensor index, with elements ordered lowest to highest. */
      operator IndexArray () const;
     
      /*! The position of the element in the array of tensor values. */
      size_type position() const;
      
      /*! The number of ordered index arrays with each element occurring the number of times specified by the multi index. */
      size_type number() const;
      
      /*! Set the value of the \a i th index to \a j. */
      void set_index(const size_type& i, const size_type j) { this->_degree+=(j-this->_occurrences[i]); this->_occurrences[i]=j; }
      /*! Increment the the \a i th index, thereby increasing the degree. */
      void increment_index(const size_type& i) { ++this->_occurrences[i]; ++this->_degree; }
      /*! Decrement the the \a i th index, thereby decreasing the degree. */
      void decrement_index(const size_type& i) { 
        if(!(this->_occurrences[i]>0)) { 
          //throw std::runtime_error(__PRETTY_FUNCTION__ ": the number of occurence of the index must be positive"); }
          throw std::runtime_error(__PRETTY_FUNCTION__ ); }
        --this->_occurrences[i]; --this->_degree; }
      
      /*! Set the value of the \a i th index to \a j. */
      void set(const size_type& i, const size_type j) { this->set_index(i,j); }
      friend std::ostream& operator<<(std::ostream&, const MultiIndex&);
     private:
      size_type _degree;
      array<size_type> _occurrences;
    };
      
    inline
    size_type MultiIndex::number() const
    {
      size_type result=Numeric::fac(this->degree());
      for(size_type k=0; k!=this->number_of_variables(); ++k) {
        result/=Numeric::fac((*this)[k]);
      }
      //std::cerr << "number(" << (*this) << ")=" << result << " " << std::flush;
      return result;
    }
      
    inline
    size_type MultiIndex::position() const
    {
      size_type result=0;
      size_type deg=this->degree();
      size_type nvar=this->number_of_variables();
      for(size_type k=0; k!=this->number_of_variables()-1; ++k) {
        --nvar;
        deg-=(*this)[k];
        result+=Numeric::bin(deg+nvar-1,nvar);
      }
      //std::cerr << "position(" << (*this) << ")=" << result << " " << std::flush;
      return result;
    }
      
    inline
    MultiIndex::operator IndexArray () const 
    {
      IndexArray result(this->degree()); 
      size_type k=0;
      for(size_type i=0; i!=this->number_of_variables(); ++i) {
        for(size_type j=k; j!=k+(*this)[i]; ++j) {
          result[j]=i;
        }
        k=k+(*this)[i];
      }
      return result;
    }
      
    inline
    bool MultiIndex::operator<(const MultiIndex& a2) const {
      const MultiIndex& a1=*this;
      if(a1.degree()!=a2.degree()) {
        return a1.degree()<a2.degree();
      } else {
        for(size_type i=0; i!=a1.number_of_variables(); ++i) {
          if(a1[i]!=a2[i]) { 
            return a1[i]<a2[i];
          }
        }
        return false;
      }
    }

    inline
    MultiIndex MultiIndex::operator+(const MultiIndex& a2) const {
      MultiIndex result(a2.number_of_variables());
      const MultiIndex& a1=*this;
      for(size_type i=0; i!=a1.number_of_variables(); ++i) {
        result.set_index(i,a1[i]+a2[i]);
      }
      return result;
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
      MultiIndexIterator(const size_type& nv, const size_type& d)
        : _index(MultiIndex(nv)) { this->_index.set_index(0,d); }
       
      void increment() {
        const size_type nv=this->_index.number_of_variables();
        if(this->_index[nv-2]!=0) {
          this->_index.decrement_index(nv-2);
          this->_index.increment_index(nv-1);
          return;
        } else {
          size_type li=this->_index[nv-1];
          this->_index.set_index(nv-1,0);
          for(size_type k=nv-1; k!=0; --k) {
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
              
      
            
      
      
    inline 
    std::ostream& operator<<(std::ostream& os, const Index& a) {
      return os << a._entries;
    }
    
    inline 
    std::ostream& operator<<(std::ostream& os, const MultiIndex& a) {
      return os << a._occurrences;
    }
    
  }
}
#endif /* ARIADNE_MULTI_INDEX_H */
