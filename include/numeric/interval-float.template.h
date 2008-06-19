/***************************************************************************
 *            interval-float.template.h
 *
 *  Copyright  2007-8  Pieter Collins
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
 

namespace Ariadne {


template<class T> void sqrt_(Interval< Float<T> >& r, const Interval< Float<T> >& x);
template<class T> void hypot_(Interval< Float<T> >& r, const Interval< Float<T> >& x, const Interval< Float<T> >& y);

template<class T> void exp_(Interval< Float<T> >& r, const Interval< Float<T> >& x);
template<class T> void log_(Interval< Float<T> >& r, const Interval< Float<T> >& x);
template<class T> void pi_(Interval< Float<T> >& r);
template<class T> void sin_(Interval< Float<T> >& r, const Interval< Float<T> >& x);
template<class T> void cos_(Interval< Float<T> >& r, const Interval< Float<T> >& x);
template<class T> void tan_(Interval< Float<T> >& r, const Interval< Float<T> >& x);
template<class T> void asin_(Interval< Float<T> >& r, const Interval< Float<T> >& x);
template<class T> void acos_(Interval< Float<T> >& r, const Interval< Float<T> >& x);
template<class T> void atan_(Interval< Float<T> >& r, const Interval< Float<T> >& x);



} // namespace Ariadne

