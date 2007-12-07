/***************************************************************************
 *            integer.template.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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
	

template<class R, class A> 
void fac_(R& r, const A& n) {
  if(n<=1) { r=1; return; }
  if(n<=12) { int m(n); r=factorials[m]; return; }
  A i=n; r=479001600;
  while(i!=12) { r*=i; --i; }	
}


template<class R, class A1, class A2>  
void bin_(R& r, const A1& n, const A2& k) 
{
  //std::cerr << "bin(" << n << "," << k << ")=" << std::flush;
  if(k==0 || k==n) { r=1; return; }
  if(k<0 || k>n) { r=0; return; }
  A2 m=(n-k < k) ? k : static_cast<A2>(n-k);
  R result=1;
  for(A1 i=n; i!=n-m; --i) { result*=i; }
  for(A1 i=m; i!=1; --i) { result/=i; }
  //std::cerr << result << std::endl;
  r=result;
}

   
   

template<class R, class A1, class A2>  
void gcd_(R& r, const A1& a, const A2& b) {
  R aa=a; R bb=b; R cc=aa%bb;
  while(cc!=0) { aa=bb; bb=cc; cc=aa%bb; }
  r=bb;
}

template<class R, class A1, class A2>  
void lcm_(R& r, const A1& a, const A2& b) {
  R res; quot(res,a*b,gcd(a,b)); return res;
}



template<class R, class N> 
void log2_floor_(R& r, const N& n) {
  if(n<1) { throw std::invalid_argument(__PRETTY_FUNCTION__); }
  N y=n;
  r=0;
  while(y>=n) {
    y/=2;
    r+=1;
  }
}

template<class R, class N> 
void log2_ceil_(R& r, const N& n) {
  if(n<1) { throw std::invalid_argument(__PRETTY_FUNCTION__); }
  N y=n;
  r=0;
  while(y>1) {
    y/=2;
    r+=1;
  }
}





}}
