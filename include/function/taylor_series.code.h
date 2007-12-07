/***************************************************************************
 *            taylor_series.code.h
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

namespace {

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
Function::TaylorSeries<X>
Function::TaylorSeries<X>::constant(smoothness_type d, const X& c)
{
  TaylorSeries<X> result(d);
  result[0]=c;
  return result;
}

template<class X> 
Function::TaylorSeries<X>
Function::TaylorSeries<X>::variable(smoothness_type d, const X& c)
{
  TaylorSeries<X> result(d);
  result[0]=c;
  result[1]=1;
  return result;
}



template<class X> 
Function::TaylorSeries<X> 
Function::TaylorSeries<X>::rec(smoothness_type d, const X& c) 
{
  TaylorSeries<X> y(d);
  X mr = X(-1)/c;
  for(size_type i=0; i<=y.degree(); ++i) {
    y[i]=(-Numeric::fac<int>(i))*Numeric::pow(mr,i+1u);
  }
  return y;
}


template<class X> 
Function::TaylorSeries<X> 
Function::TaylorSeries<X>::pow(smoothness_type d, const X& c, const uint& k)
{
  size_type n=k;
  TaylorSeries<X> y(d);
  for(size_type i=0; i<=std::min(size_type(d),n); ++i) {
    int j=n-i;
    y[i]=(Numeric::fac<int>(n)/Numeric::fac<int>(j))*Numeric::pow(c,j);
  }
  return y;
}

template<class X> 
Function::TaylorSeries<X> 
Function::TaylorSeries<X>::sqrt(smoothness_type d, const X& c)
{
  TaylorSeries<X> y(d);
  y[0]=Numeric::sqrt(c);
  X mhr=(-0.5)/c;
  for(size_type i=1; i<=y.degree(); ++i) {
    y[i]=(2*int(i)-3)*mhr*y[i-1];
  }
  return y;
}

template<class X> 
Function::TaylorSeries<X> 
Function::TaylorSeries<X>::exp(smoothness_type d, const X& c)
{
  TaylorSeries<X> y(d);
  y[0]=Numeric::exp(c);
  for(size_type i=1; i<=y.degree(); ++i) {
    y[i]=y[0];
  }
  return y;
}

template<class X>  
Function::TaylorSeries<X> 
Function::TaylorSeries<X>::log(smoothness_type d, const X& c)
{
  TaylorSeries<X> y(d);
  y[0]=Numeric::log(c);
  X mr=(-1)/c;
  for(size_type i=1; i!=y.degree();++i) {
    y[i]=(-Numeric::fac<int>(i-1))*Numeric::pow(c,i);
  }
  return y;
}

template<class X> 
Function::TaylorSeries<X> 
Function::TaylorSeries<X>::sin(smoothness_type d, const X& c)
{
  TaylorSeries<X> y(d);
  y[0]=Numeric::sin(c);
  y[1]=Numeric::cos(c);
  for(size_type i=2; i!=d; ++i) {
    y[i]=-y[i-2];
  }
  return y;
}

template<class X> 
Function::TaylorSeries<X> 
Function::TaylorSeries<X>::cos(smoothness_type d, const X& c)
{
  TaylorSeries<X> y(d);
  y[0]=Numeric::cos(c);
  y[1]=-Numeric::sin(c);
  for(size_type i=2; i!=d; ++i) {
    y[i]=-y[i-2];
  }
  return y;
}

template<class X> 
Function::TaylorSeries<X> 
Function::TaylorSeries<X>::tan(smoothness_type d, const X& c)
{
  return sin(d,c)/cos(d,c);
}

template<class X>  
Function::TaylorSeries<X> 
Function::TaylorSeries<X>::asin(smoothness_type d, const X& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X>  
Function::TaylorSeries<X> 
Function::TaylorSeries<X>::acos(smoothness_type d, const X& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class X>  
Function::TaylorSeries<X> 
Function::TaylorSeries<X>::atan(smoothness_type d, const X& c)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}





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








Function::TaylorSeries<Numeric::Rational> 
Function::TaylorSeries<Numeric::Rational>::rec(smoothness_type d, const Numeric::Rational& c) 
{
  TaylorSeries<Numeric::Rational> y(d);
  Numeric::Rational mr = Numeric::Rational(-1)/c;
  for(size_type i=0; i<=y.degree(); ++i) {
    y[i]=(-Numeric::fac<int>(i))*Numeric::pow(mr,i+1u);
  }
  return y;
}

Function::TaylorSeries<Numeric::Rational> 
Function::TaylorSeries<Numeric::Rational>::pow(smoothness_type d, const Numeric::Rational& c, const uint& k)
{
  size_type n=k;
  TaylorSeries<Numeric::Rational> y(d);
  for(size_type i=0; i<=std::min(size_type(d),n); ++i) {
    int j=n-i;
    y[i]=(Numeric::fac<int>(n)/Numeric::fac<int>(j))*Numeric::pow(c,j);
  }
  return y;
}




template<class X>
void
Function::TaylorSeries<X>::instantiate()
{
  TaylorSeries<X>* ts=0;
  std::ostream* os=0;
  inverse(*ts,X(0));
  compose(*ts,*ts);
  *os << *ts;
}

void
Function::TaylorSeries<Numeric::Rational>::instantiate()
{
  typedef Numeric::Rational X;
  TaylorSeries<X>* ts=0;
  std::ostream* os=0;
  inverse(*ts,X(0));
  compose(*ts,*ts);
  *os << *ts;
}




}
