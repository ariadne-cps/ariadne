/***************************************************************************
 *            utility/functional.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file utility/functional.hpp
 *  \brief Classes for functional-style operations.
 */

#ifndef ARIADNE_FUNCTIONAL_PROGRAMMING_HPP
#define ARIADNE_FUNCTIONAL_PROGRAMMING_HPP

#include <cassert>

#include "metaprogramming.hpp"

namespace Ariadne {

template<class T> class List;

template<class F, class T> List<ResultOf<F(T)>> elementwise(F const& f, List<T> const& l) {
    typedef ResultOf<F(T)> R;
    List<R> r; r.reserve(l.size());
    for(SizeType i=0; i!=l.size(); ++i) { r.append(f(l[i])); }
    return r;
}

template<class F, class T1, class T2> List<ResultOf<F(T1,T2)>> elementwise(F const& f, List<T1> const& l1, List<T2> const& l2) {
    typedef ResultOf<F(T1,T2)> R;
    assert(l1.size()==l2.size());
    List<R> r; r.reserve(l1.size());
    for(SizeType i=0; i!=l1.size(); ++i) { r.append(f(l1[i],l2[i])); }
    return r;
}

} // namespace Ariadne

#endif
