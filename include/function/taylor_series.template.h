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
Function::TaylorSeries<X>
Function::TaylorSeries<X>::constant(smoothness_type d, const XX& c)
{
  TaylorSeries<X> result(d);
  result[0]=c;
  return result;
}


template<class X> template<class XX>
Function::TaylorSeries<X>
Function::TaylorSeries<X>::variable(smoothness_type d, const XX& c)
{
  TaylorSeries<X> result(d);
  result[0]=c;
  result[1]=1;
  return result;
}



namespace {

template<class X> 
void 
compute_product(Function::TaylorSeries<X>& z, const Function::TaylorSeries<X>& y, const Function::TaylorSeries<X>& x)
{
  for(size_type n=0; n<=z.degree(); ++n) {
    z[n]=0;
    for(size_type i=0; i<=n; ++i) {
      z[n] += Numeric::bin<int>(n,i)*x[i]*y[n-i];
    }
  }
}


template<class X> 
void 
compute_composition(Function::TaylorSeries<X>& y, const Function::TaylorSeries<X>& x)
{
  using namespace Function;

  assert(y.degree()==x.degree());
  int d=y.degree();

#ifdef DEBUG
  // Compute using explicit formula
  TaylorSeries<X> t(y);
  TaylorSeries<X> z(std::min(d,5));
  if(d>=0) { z[0] = y[0]; }
  if(d>=1) { z[1] = y[1]*x[1]; }
  if(d>=2) { z[2] = y[2]*x[1]*x[1] + y[1]*x[2]; }
  if(d>=3) { z[3] = y[3]*x[1]*x[1]*x[1] + 3*y[2]*x[1]*x[2] + y[1]*x[3]; }
  if(d>=4) { z[4] = y[4]*x[1]*x[1]*x[1]*x[1] + 6*y[3]*x[1]*x[1]*x[2] + 4*y[2]*x[1]*x[3] + 3*y[2]*x[2]*x[2] + y[1]*x[4]; }
  if(d>=5) { z[5] = y[5]*x[1]*x[1]*x[1]*x[1]*x[1] + 10*y[4]*x[1]*x[1]*x[1]*x[2] + 15*y[3]*x[1]*x[2]*x[2] 
      + 10*y[3]*x[1]*x[1]*x[3] + 10*y[2]*x[2]*x[3] + 5*y[2]*x[1]*x[4] + y[1]*x[5]; }
#endif
  
  using namespace std;
  //cerr<<"y="<<y<<"\nx="<<x<<endl;
  TaylorSeries<X> w=x;
  w.value()=0;
  //cerr<<"w="<<w<<endl<<endl;
  TaylorSeries<X> t(d);
  TaylorSeries<X> u(d);
  t[0]=y[d]/Numeric::fac<Numeric::Integer>(d);;
  //cerr<<"t="<<t<<endl;
  for(int n=1; n<=d; ++n) {
    u=t*w;
    u.value()+=y[d-n]/Numeric::fac<Numeric::Integer>(d-n);
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
        //std::cout<<" + "<<Numeric::bin(n-i,j)<<"*y["<<d-i-j<<"]*x["<<j+1<<"]"<<std::flush;
        y[d-i] += Numeric::bin<int>(n-i,j) * y[d-i-j] * x[j+1];
      }
      //std::cout << std::endl;
    }
    cout << y << endl;
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
Function::TaylorSeries<X>
algebraic_inverse(const Function::TaylorSeries<X>& x, const X& c)
{
  smoothness_type d=x.degree();
  Function::TaylorSeries<X> y(d);

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
Function::TaylorSeries<X>
iterative_inverse(const Function::TaylorSeries<X>& x, const X& c)
{
  assert(bool(x[1]!=0));
  
  smoothness_type d=x.degree();
  X one=1;
  X r=1/x[1];

  Function::TaylorSeries<X> id(d,x.value(),one);
  Function::TaylorSeries<X> y(d,c,one);

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
Function::TaylorSeries<X>
recursive_inverse(const Function::TaylorSeries<X>& x, const X& c)
{
  assert(bool(x[1]!=0));
  
  smoothness_type d=x.degree();
  X one=1;
  X r=one/x[1];

  Function::TaylorSeries<X> id=Function::TaylorSeries<X>::variable(d,1);
  Function::TaylorSeries<X> y(d);
  
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
std::ostream& 
Function::operator<<(std::ostream& os, const TaylorSeries<X>& x) {
  os << "TaylorSeries";
  for(size_type i=0; i<=x.degree(); ++i) {
    os << (i==0 ? '(' : ',') << x[i]; 
  }
  os << ")";
  return os;
}


template<class X> 
Function::TaylorSeries<X> 
Function::antiderivative(const TaylorSeries<X>& x, const X& c)
{
  TaylorSeries<X> result(x.degree());
  for(size_type n=0; n<x.degree(); ++n) { result[n+1]=x[n]/(n+1); }
  result[0]=c;
  return result;
}

template<class X> 
Function::TaylorSeries<X> 
Function::compose(const TaylorSeries<X>& y, const TaylorSeries<X>& x)
{
  TaylorSeries<X> result(std::min(x.degree(),y.degree()));
  for(size_type n=0; n<=result.degree(); ++n) { result[n]=y[n]; }
  compute_composition(result,x);
  return result;
}

template<class X>
Function::TaylorSeries<X>
Function::inverse(const TaylorSeries<X>& x, const X& c)
{
  return recursive_inverse(x,c);
}









}
