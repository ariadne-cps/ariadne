/***************************************************************************
 *            taylor_series.template.h
 *
 *  Copyright 2007  Pieter Collins
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


template<class X> template<class XX>
TaylorSeries<X>
TaylorSeries<X>::constant(smoothness_type d, const XX& c)
{
  TaylorSeries<X> result(d);
  result[0]=c;
  return result;
}


template<class X> template<class XX>
TaylorSeries<X>
TaylorSeries<X>::variable(smoothness_type d, const XX& c)
{
  TaylorSeries<X> result(d);
  result[0]=c;
  if(d>=1) {
    result[1]=1;
  }
  return result;
}


}
