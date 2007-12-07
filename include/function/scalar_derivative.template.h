/***************************************************************************
 *            scalar_derivative.template.h
 *
 *  Copyright 2007  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

// Compute inverse algebraically using the exact formula. 
// Only inplemented for degree up to 5.
// Useful for testing.
template<class X>
Function::ScalarDerivative<X>
algebraic_inverse(const Function::ScalarDerivative<X>& x, const X& c)
{
  smoothness_type d=x.degree();
  Function::ScalarDerivative<X> y(d);

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
Function::ScalarDerivative<X>
iterative_inverse(const Function::ScalarDerivative<X>& x, const X& c)
{
  assert(bool(x[1]!=0));
  
  smoothness_type d=x.degree();
  X one=1;
  X r=1/x[1];

  Function::ScalarDerivative<X> id(d,x.value(),one);
  Function::ScalarDerivative<X> y(d,c,one);

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
Function::ScalarDerivative<X>
recursive_inverse(const Function::ScalarDerivative<X>& x, const X& c)
{
  assert(bool(x[1]!=0));
  
  smoothness_type d=x.degree();
  X one=1;
  X r=one/x[1];

  Function::ScalarDerivative<X> id(d,x.value(),one);
  Function::ScalarDerivative<X> y(d,c,r);

#ifdef DEBUG
  std::cerr << "recursive_inverse(x,c)" << std::endl;
  std::cerr << "x=" << x << std::endl;
  std::cerr << "c=" << c << std::endl;
  std::cerr << "y=" << y << std::endl;
#endif
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


// The composition of two derivatives, computed in-place. 
//
//  The composition inductively by
//  \f[ y^\[n\] = \sum_{i=0}^{n-1} \choose{n}{i} {\dot{y}}^{[i]} (x-c)^{(n-i)}
//
template<class X> 
void 
Function::compute_composition(ScalarDerivative<X>& y, const ScalarDerivative<X>& x)
{
  assert(y.degree()==x.degree());
  int d=y.degree();

#ifdef DEBUG
  // Compute using explicit formula
  ScalarDerivative<X> t(y);
  ScalarDerivative<X> z(std::min(d,5));
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

template<class X>
Function::ScalarDerivative<X>
Function::inverse(const ScalarDerivative<X>& x, const X& c)
{
  return recursive_inverse(x,c);
}






} //namespace Ariadne
