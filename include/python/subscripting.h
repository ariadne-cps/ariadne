/***************************************************************************
 *            python/subscripting.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
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

/*! \file python/subscripting.h
 *  Commonly used inline methods for the Python interface.
 */
 
#ifndef ARIADNE_PYTHON_SUBSCRIPTING_H
#define ARIADNE_PYTHON_SUBSCRIPTING_H

#include <cstring>
#include <functional>

#include "base/types.h"

#include "numeric/rational.h"
#include "python/float.h"

namespace Ariadne {
namespace Python {

  template<class C> 
  typename C::value_type 
  __getitem__(const C& c, int n) {
    if(n<0) {
      n+=c.size();
    }
    if(n<0) { throw std::out_of_range("Index out-of-range"); }
    size_t m=size_t(n);
    if(c.size()<=m) { throw std::out_of_range("Index out-of-range"); }
    return c[m];
  }


  template<class C, class T> 
  void
  __setitem__(C& c, int n, const T& x) {
    if(n<0) {
      n+=c.size();
    }
    if(n<0) { throw std::out_of_range("Index out-of-range"); }
    size_t m=size_t(n);
    if(c.size()<=m) { throw std::out_of_range("Index out-of-range"); }
    c[n]=x;
  }


/*
  template<class C> 
  inline
  void
  __setitem__(C& c, int n, const typename C::value_type& x) {
    if(n<0) {
      n+=c.size();
    }
    if(n<0) { throw std::out_of_range("Index out-of-range"); }
    size_t m=size_t(n);
    if(c.size()<=m) { throw std::out_of_range("Index out-of-range"); }
    c[n]=x;
  }
*/


}}

#endif /* ARIADNE_PYTHON_SUBSCRIPTING_H */
