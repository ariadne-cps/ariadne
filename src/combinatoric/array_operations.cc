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

#include "base/array.h"
#include "base/exceptions.h"

#include "combinatoric/exceptions.h"
#include "combinatoric/array_operations.h"

namespace Ariadne {
  namespace Base {

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
    operator&=(BooleanArray& v1, const BooleanArray& v2)
    {
      check_equal_array_sizes(v1,v2,__PRETTY_FUNCTION__);
      typedef BooleanArray::const_iterator const_iterator;
      typedef BooleanArray::iterator iterator;

      iterator v1_iter=v1.begin();
      const_iterator v2_iter=v2.begin();
      iterator v1_end=v1.end();

      while(v1_iter!=v1_end) {
        (*v1_iter) = ( (*v1_iter) & (*v2_iter) );
        ++v1_iter;
        ++v2_iter;
      }
      return v1;
    }

    BooleanArray&
    operator|=(BooleanArray& v1, const BooleanArray& v2)
    {
      check_equal_array_sizes(v1,v2,__PRETTY_FUNCTION__);
      typedef BooleanArray::const_iterator const_iterator;
      typedef BooleanArray::iterator iterator;

      iterator v1_iter=v1.begin();
      const_iterator v2_iter=v2.begin();
      iterator v1_end=v1.end();

      while(v1_iter!=v1_end) {
        (*v1_iter) = ( (*v1_iter) | (*v2_iter) );
        ++v1_iter;
        ++v2_iter;
      }
      return v1;
    }

    BooleanArray&
    operator-=(BooleanArray& v1, const BooleanArray& v2)
    {
      check_equal_array_sizes(v1,v2,__PRETTY_FUNCTION__);
      typedef BooleanArray::const_iterator const_iterator;
      typedef BooleanArray::iterator iterator;

      iterator v1_iter=v1.begin();
      const_iterator v2_iter=v2.begin();
      iterator v1_end=v1.end();

      while(v1_iter!=v1_end) {
        (*v1_iter) = ( (*v1_iter) & (!(*v2_iter)) );
        ++v1_iter;
        ++v2_iter;
      }
      return v1;
    }

    BooleanArray
    operator&(const BooleanArray& v1, const BooleanArray& v2)
    {
      BooleanArray result(v1);
      return result&=v2;
    }

    BooleanArray
    operator|(const BooleanArray& v1, const BooleanArray& v2)
    {
      BooleanArray result(v1);
      return result|=v2;
    }

    BooleanArray
    operator-(const BooleanArray& v1, const BooleanArray& v2)
    {
      BooleanArray result(v1);
      return result-=v2;
    }

    bool
    operator<=(const BooleanArray& v1, const BooleanArray& v2)
    {
      check_equal_array_sizes(v1,v2,__PRETTY_FUNCTION__);
      typedef BooleanArray::const_iterator const_iterator;
      const_iterator v1_iter=v1.begin();
      const_iterator v2_iter=v2.begin();
      const_iterator v1_end=v1.end();
      while(v1_iter!=v1_end) {
        if((*v1_iter) && (!*v2_iter)) {
          return false;
        }
        ++v1_iter;
        ++v2_iter;
      }
      return true;
    }

    bool
    lexicographic_less(const IndexArray& s1, const IndexArray& s2)
    {
      check_equal_array_sizes(s1,s2,__PRETTY_FUNCTION__);
      for(dimension_type i=0; i!=s1.size(); ++i) {
        if(s1[i]<s2[i]) {
          return true;
        }
        if(s1[i]>s2[i]) {
          return false;
        }
      }
      return false;
    }

    bool
    lexicographic_less_equal(const IndexArray& s1, const IndexArray& s2)
    {
      check_equal_array_sizes(s1,s2,__PRETTY_FUNCTION__);
      for(dimension_type i=0; i!=s1.size(); ++i) {
        if(s1[i]<s2[i]) {
          return true;
        }
        if(s1[i]>s2[i]) {
          return false;
        }
      }
      return true;
    }

    bool
    coordinate_less_equal(const IndexArray& s1, const IndexArray& s2)
    {
      check_equal_array_sizes(s1,s2,__PRETTY_FUNCTION__);
      for(dimension_type i=0; i!=s1.size(); ++i) {
        if(s1[i]>s2[i]) {
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
