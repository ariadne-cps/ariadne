#include <iostream>
#include "macros.h"
#include "numeric.h"

namespace Ariadne {

uint 
fac(uint n) 
{ 
  ARIADNE_ASSERT(n<13); 
  uint r=1; 
  for(uint i=1; i<=n; ++i) { 
    r*=i; 
  } 
  return r; 
}


uint 
bin(uint n, uint k) 
{ 
  ARIADNE_ASSERT(n<13); 
  if(k>n+1) { std::cerr << "bin("<<n<<","<<k<<")\n"; }
  if(k==n+1) { return 0; }
  ARIADNE_ASSERT(k>=0 && k<=n);
  uint r=1; 
  for(uint i=1; i<=k; ++i) { 
    r*=(n+1-i); 
    r/=i; 
  } 
  return r; 
}


Interval rec(Interval i) 
{
  ARIADNE_ASSERT(i.u<0 || i.l>0);
  return Interval(down(1/i.u),up(1/i.l)); 
}

Interval add(Interval i1, Interval i2) 
{
  return Interval(down(i1.l+i2.l),up(i1.u+i2.u));
}

Interval sub(Interval i1, Interval i2) 
{
  return Interval(down(i1.l-i2.u),up(i1.u-i2.l));
}

Interval mul(Interval i1, Interval i2) 
{
  if(i1.l>=0) {
    if(i2.l>=0) {
      return Interval(down(i1.l*i2.l),up(i1.u*i2.u));
    } else if(i2.u<=0) {
      return Interval(down(i1.l*i2.u),up(i1.u*i2.l));
    } else {
      return Interval(down(i1.u*i2.l),up(i1.u*i2.u));
    }
  }
  else if(i1.u<=0) {
    if(i2.l>=0) {
      return Interval(down(i1.l*i2.u),up(i1.u*i2.l));
    } else if(i2.u<=0) {
      return Interval(down(i1.u*i2.u),up(i1.l*i2.l));
    } else {
      return Interval(down(i1.l*i2.u),up(i1.l*i2.l));
    }
  } else {
    if(i2.l>=0) {
      return Interval(down(i1.l*i2.u),up(i1.u*i2.u));
    } else if(i2.u<=0) {
      return Interval(down(i1.u*i2.l),up(i1.l*i2.l));;
    } else {
      return Interval(down(min(i1.u*i2.l,i1.l*i2.u)),up(max(i1.l*i2.l,i1.u*i2.u)));
    }
  }
}


Interval mul(Float x, Interval i) 
{
  if(x>=0) {
    return Interval(down(x*i.l),up(x*i.u)); 
  } else {
    return Interval(down(x*i.u),up(x*i.l));
  }
}


Interval trunc(Interval x, uint n) 
{
  Interval e=Interval(pow(2.0,52-n));
  Interval y=x+e;
  return y-e;
}

Interval abs(Interval i) 
{
  if(i.lower()>=0) {
    return i;
  } else if(i.upper()<=0) {
    return -i;
  } else {
    return Interval(0,max(-i.lower(),i.upper()));
  }
}


Interval sqr(Interval i) 
{
  if(i.l >=0) {
    return Interval(down(i.l*i.l),up(i.u*i.u));
  } else if(i.u<=0) {
    return Interval(down(i.u*i.u),up(i.l*i.l));
  } else {
    return Interval(0.0,up(max(i.l*i.l,i.u*i.u)));
  }
}

Interval pow(Interval i, int n) 
{
  if(n<0) { return pow(rec(i),-n); }
  Interval r=1; Interval p=i;
  while(n>0) { if(n%2==1) { r*=p; } p=sqr(p); n/=2; }
  return r;
}

Interval pow(Interval i, uint m) 
{
  Interval r=1; Interval p=i;
  while(m>0) { if(m%2==1) { r*=p; } p=sqr(p); m/=2; }
  return r;
}

#warning Interval transcendental functions not correct

Interval sqrt(Interval i)
{
  return Interval(down(sqrt(i.l)),up(sqrt(i.u)));
}

Interval exp(Interval i)
{
  return Interval(down(exp(i.l)),up(exp(i.u)));
}

Interval log(Interval i)
{
  return Interval(down(log(i.l)),up(log(i.u)));
}


std::ostream& 
operator<<(std::ostream& os, const Interval& i)
{
  return os << '[' << i.l << ':' << i.u << ']';
}

} // namespace Ariadne

