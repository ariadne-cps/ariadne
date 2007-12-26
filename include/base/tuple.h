/***************************************************************************
 *            tuple.h
 *
 *  Copyright  2007  Pieter Collins
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
 
/*! \file tuple.h
 *  \brief Pair and tuple types, and types to be used as lvalues in assignments.
 */

#ifndef ARIADNE_TUPLE_H
#define ARIADNE_TUPLE_H

#include <boost/tuple/tuple.hpp>

namespace Ariadne { 

  namespace Base {
    using std::pair;
    using boost::tuple;

    /*! \brief A pair of references, suitable as use as an lvalue for a function returning a pair. */
    template<class T1, class T2>
    struct lpair 
    {
      inline lpair(T1& t1, T2& t2) : first(t1), second(t2) { }
      inline lpair<T1,T2> operator=(const std::pair<T1,T2>& rv) { 
        this->first=rv.first; this->second=rv.second; return *this; }
      T1& first; T2& second;
    };
  
    template<class T1,class T2> inline
    lpair<T1,T2> make_lpair(T1& t1, T2& t2) {
      return lpair<T1,T2>(t1,t2);
    }
  
  }
}

#endif /* ARIADNE_TUPLE_H */
