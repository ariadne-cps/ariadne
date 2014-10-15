/***************************************************************************
 *            utility/writable.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file utility/writable.h
 *  \brief 
 */



#ifndef ARIADNE_WRITABLE_H
#define ARIADNE_WRITABLE_H

#include "typedefs.h"
#include "metaprogramming.h"

namespace Ariadne {

/************ WritableInterface **********************************************/

class WritableInterface {
  public:
    virtual ~WritableInterface() = default;
    friend OutputStream& operator<<(OutputStream& os, const WritableInterface& w);
  protected:
    virtual OutputStream& _write(OutputStream&) const = 0;
};
inline OutputStream& operator<<(OutputStream& os, const WritableInterface& w) { w._write(os); return os; }


template<class T, class = Fallback> struct IsWritable : False { };
template<class T> struct IsWritable<T, EnableIf<IsDefined<decltype(declval<T>()._write(declval<OutputStream>()))>,Fallback>> : True { };

template<class T> EnableIf<IsWritable<T>,OutputStream&> operator<<(OutputStream& os, const T& t) {
    return t._write(os);
}

} // namespace Ariadne

#endif
