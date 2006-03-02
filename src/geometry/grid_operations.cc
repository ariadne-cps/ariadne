/***************************************************************************
 *            grid_operations.cc
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

#include <iostream>

#include "geometry/grid_operations.h"
#include "geometry/unit_grid_set.h"

namespace Ariadne {
  namespace Geometry {

    /* Inner product of two positive arrays */
    size_type
    inner_product(const array<size_type>& a1, const array<size_type>& a2)
    {
      size_type result=0;
      for(array<size_type>::size_type i=0; i!=a1.size(); ++i) {
        result += (a1[i] * a2[i]);
      }
      return result;
    }

    /* Compute the sum of an index array and a size. */
    IndexArray
    operator+(const IndexArray& l, const SizeArray& s)
    {
      IndexArray result(l.size());
      for(dimension_type i=0; i!=result.size(); ++i) {
         result[i]=l[i]+s[i];
      }
      return result;
    }

    /*! Compute a positive offset from two index sets */
    SizeArray
    operator-(const IndexArray& u, const IndexArray& l)
    {
      SizeArray result(l.size());
      for(dimension_type i=0; i!=result.size(); ++i) {
        assert(l[i]<=u[i]);
        result[i]=size_type(u[i]-l[i]);
      }
      return result;
    }

    BooleanArray&
    operator&=(BooleanArray& v1, const BooleanArray& v2)
    {
      assert(v1.size()==v2.size());
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
      assert(v1.size()==v2.size());
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
      assert(v1.size()==v2.size());
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
      assert(v1.size()==v2.size());
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
    lexicographic_order(const IndexArray& s1, const IndexArray& s2)
    {
      assert(s1.size() == s2.size());
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
    coordinate_order(const IndexArray& s1, const IndexArray& s2)
    {
      assert(s1.size() == s2.size());
      for(dimension_type i=0; i!=s1.size(); ++i) {
        if(s1[i]>s2[i]) {
          return false;
        }
      }
      return true;
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


    /* Compute strides from a list of sizes. */
    SizeArray
    compute_strides(const SizeArray& s) {
      SizeArray result(s.size()+1);
      result[0]=1;
      for(dimension_type i=0; i!=s.size(); ++i) {
        result[i+1]=s[i]*result[i];
      }
      return result;
    }


    /* Compute the index of a position in a grid. */
    size_type
    compute_index(const IndexArray& pos, const IndexArray& lower, const SizeArray& strides)
    {
      dimension_type dim=pos.size();
      size_type result=0;
      for(dimension_type i=0; i!=dim; ++i) {
        result += size_type(pos[i]-lower[i])*strides[i];
      }
      return result;
    }

    /* Compute the index of a position in a grid. */
    size_type
    compute_index(const IndexArray& pos, const SizeArray& strides, const index_type offset)
    {
      dimension_type dim=pos.size();
      index_type result=offset;
      for(dimension_type i=0; i!=dim; ++i) {
        result += pos[i]*strides[i];
      }
      return size_type(result);
    }

    /* Compute the position of an index in a grid. */
    IndexArray
    compute_position(size_type index, const IndexArray& lower, const SizeArray& strides)
    {
      dimension_type dim=lower.size();
      IndexArray result(dim);
      for(dimension_type i=dim-1; i!=0; --i) {
        result[i] = index/strides[i]+lower[i];
        index = index%strides[i];
      }
      result[0]=index;
      return result;
    }

  } 
}
