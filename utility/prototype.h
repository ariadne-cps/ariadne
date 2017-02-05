/***************************************************************************
 *            utility/prototype.h
 *
 *  Copyright 2013-17  Pieter Collins
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

/*! \file utility/prototype.h
 *  \brief
 */

#ifndef ARIADNE_PROTOTYPE_H
#define ARIADNE_PROTOTYPE_H

#include <utility>

namespace Ariadne {

template<class T> class Prototype {
    T _prototype;
  public:
    Prototype(T t) : _prototype(t) { }
    template<class X> X convert(X const& x) { return _prototype.convert(x); }
};
template<class T> inline Prototype<T> make_prototype(T const& t) {
    return Prototype<T>(t); }

} // namespace Ariadne

#endif
