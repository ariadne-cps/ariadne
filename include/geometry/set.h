/***************************************************************************
 *            set.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file set.h
 *  \brief General sets.
 */

#ifndef _ARIADNE_SET_H
#define _ARIADNE_SET_H

#include <iosfwd>

#include "../declarations.h"

#include "../base/tribool.h"

namespace Ariadne {
  namespace Geometry {

    //! \ingroup ExactSet
    /*! An abstract base class for general sets. */
    template<class R>
    class Set {
     public:
      typedef R real_type;
      typedef Point<R> state_type;
     
      virtual ~Set();
      virtual Set<R>* clone() const = 0;
     
      virtual dimension_type dimension() const = 0;
      virtual tribool contains(const Point<R>&) const = 0;
     
      virtual Rectangle<R> bounding_box() const = 0;
      virtual tribool disjoint(const Rectangle<R>&) const = 0;
      virtual tribool superset(const Rectangle<R>&) const = 0;
      virtual tribool subset(const Rectangle<R>&) const = 0;
      
      virtual std::ostream& write(std::ostream& os) const = 0;
    };
     
    template<class R> 
    Set<R>::~Set() 
    {
    }


    template<class R> inline tribool disjoint(const Set<R>& A, const Rectangle<R>& B) {
      return A.disjoint(B);
    }
    
    template<class R> inline tribool disjoint(const Rectangle<R>& A, const Set<R>& B) {
      return B.disjoint(A);
    }
    
    template<class R> inline tribool subset(const Rectangle<R>& A, const Set<R>& B) {
      return B.superset(A);
    }
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Set<R>& set) {
      return set.write(os);
    }
    
  }
}


#endif /* _ARIADNE_SET_H */
