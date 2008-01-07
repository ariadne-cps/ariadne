/***************************************************************************
 *            affine_variable.template.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
#include "affine_variable.h"

namespace Ariadne {


template<class X> inline
void
Function::neg(AffineVariable<X>& rv, const AffineVariable<X>& av) 
{
  ARIADNE_LOG(3,"neg(AffineVariable av)\n");
  rv._x=-av._x;
  rv._dx=-av._dx;
}

template<class X> inline
void
Function::rec(AffineVariable<X>& rv, const AffineVariable<X>& av) 
{
  ARIADNE_LOG(3,"rec(AffineVariable av)\n");
  rv._x=1/av._x;
  //rv._dx=-av1._dx;
}

template<class X> inline
void
Function::add(AffineVariable<X>& rv, const AffineVariable<X>& av1, const AffineVariable<X>& av2) 
{
  ARIADNE_LOG(3,"add(AffineVariable av1, AffineVariable av2)\n");
  rv._x=av1._x+av2._x;
  rv._dx=av1._dx+av2._dx;
}

template<class X> inline
void
Function::sub(AffineVariable<X>& rv, const AffineVariable<X>& av1, const AffineVariable<X>& av2) 
{
  ARIADNE_LOG(3,"sub(AffineVariable av1, AffineVariable av2)\n");
  rv._x=av1._x-av2._x;
  rv._dx=av1._dx-av2._dx;
}

template<class X> inline
void
Function::mul(AffineVariable<X>& rv, const AffineVariable<X>& av1, const AffineVariable<X>& av2) 
{
  ARIADNE_LOG(3,"mul(AffineVariable av1, AffineVariable av2)\n");
  //  rv._x=av._x-av._dx;
  //  rv._dx=av1._dx-av2._dx;
}

template<class X> inline
void
Function::div(AffineVariable<X>& rv, const AffineVariable<X>& av1, const AffineVariable<X>& av2) 
{
  ARIADNE_LOG(3,"div(AffineVariable av1, AffineVariable av2)\n");
  AffineVariable<X> tv(av2.argument_dimension());
  rec(tv,av2);
  mul(rv,av1,tv);
}


template<class X> inline
void
Function::compose(AffineVariable<X>& rv, const AffineVariable<X>& av1, const AffineVariable<X>& av2) 
{
  ARIADNE_LOG(3,"AffineVariable compose(AffineVariable av1, AffineVariable av2)\n");
  //FIXME: Use slices
  typedef typename AffineVariable<X>::I I;
  using namespace LinearAlgebra;
  assert(av1._dx.size()==1u);
  rv._x=av1._x+av1._dx[0]*av2._a;
  rv._dx=av1._dx[0]*av2._dx;
}

template<class X> inline
void
Function::reduce(AffineVariable<X>& rv, const AffineVariable<X>& av, smoothness_type s) 
{
  assert(s<=av._s);
  if(s==av._s) { return av;  }
  
  const dimension_type ad=av._ad;
  
  AffineVariable<X> res(ad);
  typename AffineVariable<X>::I tmp;
  res._x=av._x;
  for(size_type j=0; j!=ad; ++j) {
    rv._dx[j]=av._dx[j].midpoint();
    rv._x[j]+=(av._dx[j]-rv._dx[j]);
  }
  return res;
}



template<class X>
std::ostream& 
Function::operator<<(std::ostream& os, const AffineVariable<X>& av)
{
  return os << "AffineVariable( a0=" << av._x << ", a1=" << av._dx << ")";
}



}
