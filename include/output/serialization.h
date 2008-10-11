/****************************************************************************
 *            serialization.h
 *
 *  Copyright  2008  Pieter Collins
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

#ifndef ARIADNE_SERIALIZATION_H
#define ARIADNE_SERIALIZATION_H

/*! \file serialization.h
 *  \brief Reading and writing to a boost archive.
 */
 
#include <iostream>
#include <iomanip>
#include <fstream>

#include <boost/serialization/map.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "base/array.h"
#include "numeric/float.h"
#include "numeric/interval.h"

namespace Ariadne {
  
template<class A, class X> void serialize(A& a, array<X>& ary, const uint v) {
  // We can't use separate save/load unless serialize is a member (I think).
  // We therefore need the same code to read and write. 
  // We use a trick by which we first store the current value in a temporary
  // variable, read/write the temporary, and set it back to the archive.
  uint m=ary.size(); a & m; if(m!=ary.size()) { ary.resize(m); }
  for(uint i=0; i!=m; ++i) { X k=ary[i]; a & k; ary[i]=k; }
}

template<class A, class T> void serialize(A& a, Float<T>& x, const uint v) {
  a & x._value;
}

template<class A, class R> void serialize(A& a, Interval<R>& ivl, const uint v) {
  a & ivl._lower;
  a & ivl._upper;
}

template<class A> void serialize(A& a, DiscreteState& ds, const uint v) {
  id_type& id = const_cast<id_type&>(ds.id());
  a & id;
}



}


#endif /* ARIADNE_SERIALIZATION_H */
