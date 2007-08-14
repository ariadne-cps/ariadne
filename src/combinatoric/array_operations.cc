/***************************************************************************
 *            array_operations.cc
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "base/stlio.h"
#include "base/array.h"
#include "base/exceptions.h"

#include "combinatoric/exceptions.h"
#include "combinatoric/array_operations.h"

namespace Ariadne {
  namespace Combinatoric {

    size_type
    inner_product(const SizeArray& a1, const SizeArray& a2)
    {
      size_type result=0;
      for(size_type i=0; i!=a1.size(); ++i) {
        result += (a1[i] * a2[i]);
      }
      return result;
    }

    index_type
    inner_product(const IndexArray& a1, const IndexArray& a2)
    {
      index_type result=0;
      for(size_type i=0; i!=a1.size(); ++i) {
        result += (a1[i] * a2[i]);
      }
      return result;
    }

    BooleanArray&
    operator&=(BooleanArray& a1, const BooleanArray& a2)
    {
      ARIADNE_CHECK_EQUAL_ARRAY_SIZES(a1,a2,"BooleanArray& operator&=(BooleanArray a1, BooleanArray a2)");
      typedef BooleanArray::const_iterator const_iterator;
      typedef BooleanArray::iterator iterator;

      iterator a1_iter=a1.begin();
      const_iterator a2_iter=a2.begin();
      iterator a1_end=a1.end();

      while(a1_iter!=a1_end) {
        (*a1_iter) = ( (*a1_iter) & (*a2_iter) );
        ++a1_iter;
        ++a2_iter;
      }
      return a1;
    }

    BooleanArray&
    operator|=(BooleanArray& a1, const BooleanArray& a2)
    {
      ARIADNE_CHECK_EQUAL_ARRAY_SIZES(a1,a2,"BooleanArray& operator|=(BooleanArray a1, BooleanArray a2)");
      typedef BooleanArray::const_iterator const_iterator;
      typedef BooleanArray::iterator iterator;

      iterator a1_iter=a1.begin();
      const_iterator a2_iter=a2.begin();
      iterator a1_end=a1.end();

      while(a1_iter!=a1_end) {
        (*a1_iter) = ( (*a1_iter) | (*a2_iter) );
        ++a1_iter;
        ++a2_iter;
      }
      return a1;
    }

    BooleanArray&
    operator-=(BooleanArray& a1, const BooleanArray& a2)
    {
      ARIADNE_CHECK_EQUAL_ARRAY_SIZES(a1,a2,"BooleanArray& operator-=(BooleanArray a1, BooleanArray a2)");
      typedef BooleanArray::const_iterator const_iterator;
      typedef BooleanArray::iterator iterator;

      iterator a1_iter=a1.begin();
      const_iterator a2_iter=a2.begin();
      iterator a1_end=a1.end();

      while(a1_iter!=a1_end) {
        (*a1_iter) = ( (*a1_iter) & (!(*a2_iter)) );
        ++a1_iter;
        ++a2_iter;
      }
      return a1;
    }

    BooleanArray
    operator&(const BooleanArray& a1, const BooleanArray& a2)
    {
      BooleanArray result(a1);
      return result&=a2;
    }

    BooleanArray
    operator|(const BooleanArray& a1, const BooleanArray& a2)
    {
      BooleanArray result(a1);
      return result|=a2;
    }

    BooleanArray
    operator-(const BooleanArray& a1, const BooleanArray& a2)
    {
      BooleanArray result(a1);
      return result-=a2;
    }

    bool
    operator<=(const BooleanArray& a1, const BooleanArray& a2)
    {
      ARIADNE_CHECK_EQUAL_ARRAY_SIZES(a1,a2,"bool operator<=(BooleanArray a1, BooleanArray a2)");
      typedef BooleanArray::const_iterator const_iterator;
      const_iterator a1_iter=a1.begin();
      const_iterator a2_iter=a2.begin();
      const_iterator a1_end=a1.end();
      while(a1_iter!=a1_end) {
        if((*a1_iter) && (!*a2_iter)) {
          return false;
        }
        ++a1_iter;
        ++a2_iter;
      }
      return true;
    }

    bool
    lexicographic_less(const IndexArray& a1, const IndexArray& a2)
    {
      ARIADNE_CHECK_EQUAL_ARRAY_SIZES(a1,a2,"bool lexicographic_less(BooleanArray a1, BooleanArray a2)");
      for(dimension_type i=0; i!=a1.size(); ++i) {
        if(a1[i]<a2[i]) {
          return true;
        }
        if(a1[i]>a2[i]) {
          return false;
        }
      }
      return false;
    }

    bool
    lexicographic_less_equal(const IndexArray& a1, const IndexArray& a2)
    {
      ARIADNE_CHECK_EQUAL_ARRAY_SIZES(a1,a2,"bool lexicographic_less_equal(BooleanArray a1, BooleanArray a2)");
      for(dimension_type i=0; i!=a1.size(); ++i) {
        if(a1[i]<a2[i]) {
          return true;
        }
        if(a1[i]>a2[i]) {
          return false;
        }
      }
      return true;
    }

    bool
    coordinate_less_equal(const IndexArray& a1, const IndexArray& a2)
    {
      ARIADNE_CHECK_EQUAL_ARRAY_SIZES(a1,a2,"bool coordinate_less_equal(BooleanArray a1, BooleanArray a2)");
      for(dimension_type i=0; i!=a1.size(); ++i) {
        if(a1[i]>a2[i]) {
          return false;
        }
      }
      return true;
    }

    IndexArray
    operator+(const IndexArray& l, const SizeArray& s)
    {
      IndexArray result(l.size());
      for(dimension_type i=0; i!=result.size(); ++i) {
         result[i]=l[i]+s[i];
      }
      return result;
    }

    IndexArray
    operator-(const IndexArray& u, const IndexArray& l)
    {
      IndexArray result(l.size());
      for(dimension_type i=0; i!=result.size(); ++i) {
        result[i]=size_type(u[i]-l[i]);
      }
      return result;
    }

    void
    assign_max(IndexArray& a, const IndexArray& l) {
      for(dimension_type i=0; i!=a.size(); ++i) {
        if(l[i]>a[i]) {
          a[i]=l[i];
        }
      }
    }

    void
    assign_min(IndexArray& a, const IndexArray& u) {
      for(dimension_type i=0; i!=a.size(); ++i) {
        if(u[i]<a[i]) {
          a[i]=u[i];
        }
      }
    }



  } 
}
