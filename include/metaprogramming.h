/***************************************************************************
 *            metaprogramming.h
 *
 *  Copyright 2008-11  Pieter Collins
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

/*! \file metaprogramming.h
 *  \brief Classes for template metaprogramming.
 */
#ifndef ARIADNE_METAPROGRAMMING_H
#define ARIADNE_METAPROGRAMMING_H


namespace Ariadne {

typedef void Void;

struct True { static const bool value = true; };
struct False { static const bool value = false; };

template<class P1, class P2, class P3=void> struct And { static const bool value = P1::value && P2::value && P3::value; };
template<class P1, class P2> struct And<P1,P2> { static const bool value = P1::value && P2::value; };
template<class P1, class P2> struct Or { static const bool value = P1::value || P2::value; };
template<class P> struct Not { static const bool value = !P::value; };

template<bool,class> struct EnableIfBool;
template<class T> struct EnableIfBool<false,T> { };
template<class T> struct EnableIfBool<true,T> { typedef T Type; };

template<class P, class T> struct EnableIfClass : public EnableIfBool<P::value,T> { };

//template<bool B, class T> struct EnableIf : public EnableIfBool<B,T> { };
template<class P, class T=Void> struct EnableIf : public EnableIfClass<P,T> { };

template<class T1, class T2> struct IsSame : False { };
template<class T> struct IsSame<T,T> : True { };

template<class X> struct IsDefined : public True { };


} // namespace Ariadne

#endif
