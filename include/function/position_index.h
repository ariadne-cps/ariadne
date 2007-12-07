/***************************************************************************
 *            position_index.h
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
 
/*! \file symmetric_index.h
 *  \brief Position-based index for dense symmetric tensors
 */

#ifndef ARIADNE_SYMMETRIC_POSITION_H
#define ARIADNE_SYMMETRIC_POSITION_H

#include <cassert>

#include "../base/iterator.h"
#include "../base/array.h"
#include "../base/stlio.h"

#include "../numeric/integer.h"

#include "sorted_index.h"
#include "multi_index.h"

namespace Ariadne {
  namespace Function {
    
    class SortedIndex;
    class MultiIndex;

    /*! \ingroup Function
     *  \brief An index of a dense symmetric object. 
     */
    class PositionIndex {
     public:
      typedef Ariadne::Base::size_type size_type;
      typedef size_type value_type;
     
      /*! Construct a multi index of degree \a 0 with \a nv variables. */
      explicit PositionIndex(size_type nv);
      /*! Construct a multi index with \a nv variables from the array \a ary. */
      explicit PositionIndex(size_type nv, const size_type pos);
      /*! Construct a multi index from the tensor index \a a. */
      explicit PositionIndex(const SortedIndex& a);
      /*! Construct a multi index from the tensor index \a a. */
      explicit PositionIndex(const MultiIndex& a);
      /*! Copy constructor. */
      PositionIndex(const PositionIndex& a);
      /*! Copy assignment operator. */
      PositionIndex& operator=(const PositionIndex& a);

      /*! The degree of the multi-index, equal to the sum of the number of occurrences of the variables. */
      const size_type degree() const;
      /*! The number of variables. */
      const size_type& number_of_variables() const;
      /*! The number of occurrences of the \a i th variable. */
      const size_type& position() const;
     
      /*! Equality operator. */
      bool operator==(const PositionIndex& a) const;
      /*! Inequality operator. */
      bool operator!=(const PositionIndex& a) const;
      /*! Comparison operator. */
      bool operator<(const PositionIndex& a) const; // inline
      /*! Increment. */
      PositionIndex& operator++(); // inline
     
      /*! Convert to a normal tensor index, with elements ordered lowest to highest. */
      operator SortedIndex () const;
     
      /*! Convert to a normal tensor index, with elements ordered lowest to highest. */
      operator MultiIndex () const;
     
      friend std::ostream& operator<<(std::ostream&, const PositionIndex&);
     private:
      size_type _number_of_variables;
      size_type _position;
    };
      


    
    inline PositionIndex::PositionIndex(size_type nv)
      : _number_of_variables(nv), _position(0)
    {
    }

    inline PositionIndex::PositionIndex(size_type nv, size_type pos)
      : _number_of_variables(nv), _position(pos)
    {
    }

    inline PositionIndex::PositionIndex(const PositionIndex& a)
      : _number_of_variables(a._number_of_variables), _position(a._position) 
    {
    }

    inline PositionIndex& PositionIndex::operator=(const PositionIndex& a) { 
      this->_number_of_variables=a._number_of_variables; this->_position=a._position; return *this; }

    inline PositionIndex::PositionIndex(const SortedIndex& a)
      : _number_of_variables(a.number_of_variables()), _position(a.position())
    { 
    }

    inline PositionIndex::PositionIndex(const MultiIndex& a)
      : _number_of_variables(a.number_of_variables()), _position(a.position())
    { 
    }


    inline const size_type PositionIndex::degree() const { 
      return SortedIndex(*this).degree();
    }

    inline const size_type& PositionIndex::number_of_variables() const { 
      return this->_number_of_variables;
    }

    inline
    const size_type& PositionIndex::position() const
    {
      return this->_position;
    }
      
    inline bool PositionIndex::operator==(const PositionIndex& a) const { 
      assert(this->_number_of_variables==a._number_of_variables);
      return this->_position==a._position;
    }
    
    inline bool PositionIndex::operator!=(const PositionIndex& a) const { 
      return !(*this==a); 
    }

    inline PositionIndex& PositionIndex::operator++() { 
      ++this->_position; return *this;
    }

    inline
    PositionIndex::operator SortedIndex () const 
    {
      const size_type& nv=this->_number_of_variables;
      const size_type& p=this->_position;
      
      size_type d=0;
      while(Numeric::bin(nv+d,nv)<p) {
        ++d;
      }
      std::cout << "nv=" << nv << " p=" << p << std::endl;
      std::cout << "d=" << d << std::endl;
      
      array<size_type> a(d,nv);
      size_type cp=Numeric::bin(nv+d,nv);
      for(uint k=0; k!=d; ++k) {
        int m=d-k;
        int l=m+(nv-1);
        int j=nv-1;
        long c=Numeric::bin(l-j,m);
        //std::cout << "k=" << k <<" j="<<j<<" cp="<<cp<<" ("<<(l-j)<<","<<(m)<<")="<<c<<"\n";
        while(cp-c>p) {
          --j;
          c=Numeric::bin(l-j,m);
          //std::cout << "k=" << k <<" j="<<j<<" cp="<<cp<<" ("<<(l-j)<<","<<(m)<<")="<<c<<"\n";
          assert(j>-2);
        }
        a[k]=j;
        cp-=Numeric::bin(l-j-1,m);
      }
      SortedIndex r(nv,d,a.begin());
      //std::cout << "r=" << r << " r.position()=" << r.position() << "\n";
      return r;
    }
      
    inline
    PositionIndex::operator MultiIndex () const 
    {
      //std::cerr << "PositionIndex::operator MultiIndex() const\n";
      //std::cerr << "nv=" << this->number_of_variables() << " p="<<this->position() << "\n";
      MultiIndex r(this->number_of_variables());
      const size_type& p=this->position();
      while(r.position()<=p) {
        //std::cerr << r << " " << r.position() << std::endl;
        r.increment_index(0);
        assert(r[0]<20);
      }
      r.decrement_index(0);

      for(uint k=1; k!=r.number_of_variables(); ++k) {
        //std::cerr << r << " " << r.position() << std::endl;
        while(r.position()<=p && r[k-1]>0) {
          r.decrement_index(k-1);
          r.increment_index(k);
          //std::cerr << r << " " << r.position() << std::endl;
        }
        if(r.position()>p) {
          r.increment_index(k-1);
          r.decrement_index(k);
        }
      }
  
      //std::cerr << r << " " << r.position() << std::endl << std::endl;
      assert(r.position()==p);
      return r;
    }

      
    inline 
    std::ostream& operator<<(std::ostream& os, const PositionIndex& a) {
      return os << a._position;
    }
    
    inline 
    MultiIndex::MultiIndex(const PositionIndex& p) {
      *this=p.operator MultiIndex();
    }

  }
}
#endif /* ARIADNE_MULTI_INDEX_H */
