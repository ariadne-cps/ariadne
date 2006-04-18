/***************************************************************************
 *            list_set.tpl
 *
 *  23 June 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
#include "list_set.h"
#include <iostream>
#include "../geometry/rectangle.h"

namespace Ariadne {
  namespace Geometry {

    template <typename R, template<class> class BS>
    bool
    disjoint (const ListSet<R,BS>& A,
              const ListSet<R,BS>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable sets have different space dimensions.");
      }

      size_type i,j;
      for (i=0; i<A.size() ; i++) {
        for (j=0; j<B.size() ; j++) {
          if (!disjoint(A[i],B[j])) { return false; }
        }
      }

      return true;
    }

    template <typename R, template <typename> class BS>
    bool
    interiors_intersect(const ListSet<R,BS>& A,
                        const ListSet<R,BS>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable sets have different space dimensions.");
      }

      size_type i,j;

      for (i=0; i<A.size() ; i++) {
        for (j=0; j<B.size() ; j++) {
          if (interiors_intersect(A[i],B[j])) { return true; }
        }
      }

      return false;
    }

    template <typename R, template <typename> class BS>
    bool
    inner_subset(const ListSet<R,BS>& A,
                 const ListSet<R,BS>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable sets have different space dimensions.");
      }

      for (size_type i=0; i<A.size() ; i++) {
        if (!inner_subset(A[i], B)) { return false; }
      }

      return true;
    }

    template <typename R, template <typename> class BS>
    bool
    subset(const ListSet<R,BS>& A,
           const ListSet<R,BS>& B)
    {
      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable sets have different space dimensions.");
      }

      throw std::runtime_error("Not implemented");
    }

    
    template <typename R, template <typename> class BS>
    ListSet<R,BS>
    join(const ListSet<R,BS>& A,
         const ListSet<R,BS>& B)
    {
      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable sets have different space dimensions.");
      }

      ListSet<R,BS> ds_union(A);

      ds_union.inplace_union(B);

      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      return ds_union;
    }

    template <typename R, template <typename> class BS>
    ListSet<R,BS>
    regular_intersection(const ListSet<R,BS>& A,
                         const ListSet<R,BS>& B)
    {
      ListSet<R,BS> ds_inter;
      std::vector< BS<R> >& vec=ds_inter._vector;

      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      if (A.dimension()!=B.dimension()) {
        throw std::invalid_argument("The two denotable sets have different space dimensions.");
      }

      for (size_type i=0; i<A.size(); i++) {
        for (size_type j=0; j<B.size(); j++) {
          if (interiors_intersect(A[i],B[j])) {
              vec.push_back(regular_intersection(A[i],B[j]));
          }
        }
      }

      ds_inter._dimension=A._dimension;

      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif

      return ds_inter;
    }
    
    
    
    template <typename R, template<typename> class BS>
    std::ostream& operator<<(std::ostream& os,
                             const ListSet<R,BS>& A)
    {
      os << "ListSet<" << name<R>() << ",BS>(\n  ";
      os << "[";
      if (A.size() >0 ) {
        os << A[0];
      }
      for (size_type i=1; i<A.size(); i++) {
        os << ",\n    " << A[i];
      }
      os << "\n  ]\n}" << std::endl;
      return os;
    }

    template <typename R, template<typename> class BS>
    std::istream& operator>>(std::istream& is,
                             ListSet<R,BS>& A)
    {
      std::vector< BS<R> >& vec(A._vector);
      is >> vec;

      if(vec.size()==0) {
        A._dimension = 0;
      }
      else {
        A._dimension=vec[0].dimension();
      }
      return is;
    }


  }
}
