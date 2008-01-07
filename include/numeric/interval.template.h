/***************************************************************************
 *            interval.template.h
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
 

namespace Ariadne {
namespace Numeric {

template<class R, class X> void sqrt_(Interval<R>& r, const Interval<X>& x);
template<class R, class X, class Y> void hypot_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y);

template<class R, class X> void exp_(Interval<R>& r, const Interval<X>& x);
template<class R, class X> void log_(Interval<R>& r, const Interval<X>& x);
template<class R> void pi_(Interval<R>& r);
template<class R, class X> void sin_(Interval<R>& r, const Interval<X>& x);
template<class R, class X> void cos_(Interval<R>& r, const Interval<X>& x);
template<class R, class X> void tan_(Interval<R>& r, const Interval<X>& x);
template<class R, class X> void asin_(Interval<R>& r, const Interval<X>& x);
template<class R, class X> void acos_(Interval<R>& r, const Interval<X>& x);
template<class R, class X> void atan_(Interval<R>& r, const Interval<X>& x);

template<class R> std::ostream& operator<<(std::ostream& os, const Interval<R>& x);
template<class R> std::istream& operator>>(std::istream& is, Interval<R>& x);


 
template<class R> 
std::ostream& operator<<(std::ostream& os, const Interval<R>& ivl) 
{
  if(ivl.lower()>ivl.upper()&&false) {
    return os << "[1:0]";
  }
  else {
    return os << "[" << ivl.lower() << ":" << ivl.upper() << "]";
  }
}


template<class R> 
std::istream& 
operator>>(std::istream& is, Interval<R>& ivl)
{
  char c;
  R l;
  R u;
  is >> c;
  if(c=='[') {
    is >> l >> c;
    if(c!=',' && c!=':') {
      is.setstate(std::ios_base::failbit);
    }
    is >> u >> c;
    if(c!=']') {
      is.setstate(std::ios_base::failbit);
    }
    ivl=Interval<R>(l,u);
  } else {
    is.putback(c);
    is >> l;
    ivl=Interval<R>(l);
  }
  return is;
}


} // namespace Numeric
} // namespace Ariadne

