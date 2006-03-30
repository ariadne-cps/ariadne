/***************************************************************************
 *            python/python_utilities.h
 *
 *  16 November 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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

/*! \file python_utilities.h
 *  Commonly used inline methods for the Python interface.
 */
 
#ifndef _ARIADNE_PYTHON_UTILITIES_H
#define _ARIADNE_PYTHON_UTILITIES_H

template<class C> 
inline
typename C::value_type 
get(const C& c, const typename C::size_type n) {
  assert(n<c.size());
  return c[n];
}

template<class C> 
inline
void
set(C& c, const typename C::size_type n, const typename C::value_type x) {
  assert(n<c.size());
  c[n]=x;
}

#endif /* _ARIADNE_PYTHON_UTILITIES_H */
