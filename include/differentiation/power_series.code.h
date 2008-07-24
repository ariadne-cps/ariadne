/***************************************************************************
 *            power_series.code.h
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
 
#include "differential.h"


namespace Ariadne {



namespace {

template<class X> 
void 
compute_product(PowerSeries<X>& z, const PowerSeries<X>& y, const PowerSeries<X>& x)
{
  for(size_type n=0; n<=z.degree(); ++n) {
    z[n]*=0;
    for(size_type i=0; i<=n; ++i) {
      z[n] += x[i]*y[n-i];
    }
  }
}


template<class X> 
void 
compute_composition(PowerSeries<X>& y, const PowerSeries<X>& x)
{
  

  int d=std::min(x.degree(),y.degree());
  using namespace std;
  //cerr<<"y="<<y<<"\nx="<<x<<endl;
  PowerSeries<X> w=x;
  w.value()=0;
  //cerr<<"w="<<w<<endl<<endl;
  PowerSeries<X> t(d);
  PowerSeries<X> u(d);
  t[0]=y[d];
  //cerr<<"t="<<t<<endl;
  for(int n=1; n<=d; ++n) {
    u=t*w;
    u.value()+=y[d-n];
    //cerr<<"u="<<u<<endl;
    t=u;
    //cerr<<"t="<<t<<endl;
  };
  y=t;
  //cerr<<endl;
  return;

  for(int n=0; n<d; ++n) {
    for(int i=0; i<=n; ++i) {
      //std::cout<<"y["<<d-i<<"] = 1*y["<<d-i<<"]*x[1]"<<std::flush;
      y[d-i]*=x[1];
      for(int j=1; j<=n-i; ++j) {
        //std::cout<<" + "<<bin(n-i,j)<<"*y["<<d-i-j<<"]*x["<<j+1<<"]"<<std::flush;
        y[d-i] += y[d-i-j] * x[j+1];
      }
      //std::cout << std::endl;
    }
    //std::cout << y << std::endl;
  }

#ifdef DEBUG
  std::cerr<<"compose(y,x)\nx="<<x<<"\ny="<<t<<"\nr="<<y<<"\nz="<<z<<std::endl;
#endif
  return;
}

// Compute inverse algebraically using the exact formula. 
// Only inplemented for degree up to 5.
// Useful for testing.
template<class X>
PowerSeries<X>
algebraic_inverse(const PowerSeries<X>& x, const X& c)
{
  smoothness_type d=x.degree();
  PowerSeries<X> y(d);

  assert(d<=5);

  if(d>=0) { y[0]=c; }
  if(d>=1) { y[1]=1/x[1]; }
  if(d>=2) { y[2]=-( (x[2])*y[1] ) / pow(x[1],2); }
  if(d>=3) { y[3]=-( (x[3])*y[1] + (2*x[1]*x[2])*y[2] ) / pow(x[1],3); }
  if(d>=4) { y[4]=-( (x[4])*y[1] + (2*x[1]*x[3]+x[2]*x[2])*y[2] + (3*x[1]*x[1]*x[2])*y[3] ) / pow(x[1],4); }
  if(d>=5) { y[5]=-( (x[5])*y[1] + (2*x[1]*x[4]+2*x[2]*x[3])*y[2] + (3*x[1]*x[1]*x[3]+3*x[1]*x[2]*x[2])*y[3] + (4*x[1]*x[1]*x[1]*x[2])*y[4] ) / pow(x[1],5); }
  std::cerr << "algebraic_inverse(x,c)" << std::endl;
  std::cerr << "x=" << x << std::endl;
  std::cerr << "c=" << c << std::endl;
  std::cerr << "y=" << y << std::endl;
  std::cerr << "compose(x,y)=" << compose(x,y) << std::endl;
  std::cerr << "compose(y,x)=" << compose(y,x) << std::endl;
  std::cerr << std::endl;

  return y;
}

// Compute inverse iteratively using a contractor. 
// Useful for interval methods.
// This method is explicit; see recursive inverse for details
template<class X>
PowerSeries<X>
iterative_inverse(const PowerSeries<X>& x, const X& c)
{
  assert(bool(x[1]!=0));
  
  smoothness_type d=x.degree();
  X one=1;
  X r=1/x[1];

  PowerSeries<X> id(d,x.value(),one);
  PowerSeries<X> y(d,c,one);

  std::cerr << "iterative_inverse(x,c)" << std::endl;
  std::cerr << "x=" << x << std::endl;
  std::cerr << "c=" << c << std::endl;
  std::cerr << "y=" << y << std::endl;
  for(smoothness_type i=0; i!=10; ++i) {
    y=y-r*(compose(x,y)-id);
    std::cerr << "y=" << y << std::endl;
  }  
  std::cerr << "compose(x,y)=" << compose(x,y) << std::endl;
  std::cerr << "compose(y,x)=" << compose(y,x) << std::endl;
  std::cerr << std::endl;

  return y;
}

// Compute inverse iteratively using a contractor. 
// This should be used as the default method
// TODO: Only compute to the degree necessary at each iteration.
template<class X>
PowerSeries<X>
recursive_inverse(const PowerSeries<X>& x, const X& c)
{
  assert(bool(x[1]!=0));
  
  smoothness_type d=x.degree();
  X one=1;
  X r=one/x[1];

  PowerSeries<X> id=PowerSeries<X>::variable(d,1);
  PowerSeries<X> y(d);
  
#ifdef DEBUG
  std::cerr << "recursive_inverse(x,c)" << std::endl;
  std::cerr << "x=" << x << std::endl;
  std::cerr << "c=" << c << std::endl;
  std::cerr << "y=" << y << std::endl;
#endif
  y[0]=c;
  y[1]=r;
  for(smoothness_type i=2; i<=d; ++i) {
    y=y-r*(compose(x,y)-id);
#ifdef DEBUG
    std::cerr << "y=" << y << std::endl;
#endif
  }  
#ifdef DEBUG
  std::cerr << "compose(x,y)=" << compose(x,y) << std::endl;
  std::cerr << "compose(y,x)=" << compose(y,x) << std::endl;
  std::cerr << std::endl;
#endif

  return y;
}


} // namespace









template<class X> 
PowerSeries<X> 
mul(const PowerSeries<X>& x, const PowerSeries<X>& y)
{
  PowerSeries<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) {
    result[n]=x[0]*y[n];
    for(size_type i=1; i<=n; ++i) {
      result[n] += x[i]*y[n-i];
    }
  }
  return result;
}


template<class X> 
PowerSeries<X> 
derivative(const PowerSeries<X>& x)
{
  PowerSeries<X> result(x.degree()-1);
  for(uint n=1; n<=x.degree(); ++n) { result[n-1]=X(n)*x[n]; }
  return result;
}

template<class X> 
PowerSeries<X> 
antiderivative(const PowerSeries<X>& x)
{
  PowerSeries<X> result(x.degree()+1);
  for(uint n=1; n<=x.degree()+1u; ++n) { result[n]=x[n-1]/n; }
  result[0]=x[0]*0;
  return result;
}

template<class X> 
PowerSeries<X> 
antiderivative(const PowerSeries<X>& x, const X& c)
{
  PowerSeries<X> result(x.degree()+1);
  for(uint n=1; n<=x.degree()+1u; ++n) { result[n]=x[n-1]/n; }
  result[0]=c;
  return result;
}

template<class X> 
PowerSeries<X> 
compose(const PowerSeries<X>& y, const PowerSeries<X>& x)
{
  using namespace std;

  size_type d=std::min(x.degree(),y.degree());
  PowerSeries<X> r(d);

  //cerr<<"y="<<y<<"\nx="<<x<<endl;
  PowerSeries<X> w=x;
  w.value()=0;
  //cerr<<"w="<<w<<endl<<endl;
  PowerSeries<X> t(d);
  r[0]=y[d];
  //cerr<<"t="<<t<<endl;
  for(uint n=1; n<=d; ++n) {
    t=r*w;
    t.value()+=y[d-n];
    //cerr<<"u="<<u<<endl;
    r=t;
    //cerr<<"t="<<t<<endl;
  };
  return r;
}

template<class X>
PowerSeries<X>
inverse(const PowerSeries<X>& x, const X& c)
{
  return recursive_inverse(x,c);
}


template<class X> 
std::ostream& 
operator<<(std::ostream& os, const PowerSeries<X>& x) {
  os << "S";
  for(size_type i=0; i<=x.degree(); ++i) {
    os << (i==0 ? '(' : ',') << x[i]; 
  }
  os << ")";
  return os;
}


template<class X> 
void
PowerSeries<X>::instantiate() 
{
  X* x=0;
  PowerSeries<X>* ts=0;
  std::ostream* os = 0;

  mul(*ts,*ts);
  compose(*ts,*ts);
  inverse(*ts,*x);
  derivative(*ts);
  antiderivative(*ts);
  antiderivative(*ts,*x);

  operator<<(*os,*ts);
}






} // namespace Ariadne

