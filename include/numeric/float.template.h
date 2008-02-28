/***************************************************************************
 *            numeric/float.template.h
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
namespace Numeric {

template<class T> class Float;

template<class T, class Rnd>
void pow_(Float<T>& r, const Float<T>& x, const unsigned int& n, Rnd rnd) 
{
  Float<T> p(x); 
  unsigned int m=n;
  r=1; 
  while(m) { 
    if(m%2) { 
      mul_(r,r,p,rnd); 
    } 
    mul_(p,p,p,rnd); 
    m/=2; 
  }
}      


template<class T, class Rnd> 
inline void pow_(Float<T>& r, const Float<T>& x, const int& n, Rnd rnd) 
{
  if(n>=0) { 
    pow_(r,x,(unsigned int)(n),rnd); 
  }
  else { 
    Float<T> t(x); 
    div_(t,1,t,rnd); 
    pow_(r,t,(unsigned int)(-n),rnd); 
  }
}


template<class T, class Rnd>  
void hypot_(Float<T>& r, const Float<T>& x, const Float<T>& y, Rnd rnd) {
  if(abs_(x)>abs_(y)) { 
    hypot_(r,y,x,Rnd()); 
  }
  Float<T> t; 
  div_(t,x,y,rnd); 
  mul_(t,t,t,rnd); 
  add_(t,1,t,rnd); 
  sqrt_(t,t,rnd); 
  mul_(r,t,y,rnd);
}

} // namespace Ariadne
} // namespace Numeric
