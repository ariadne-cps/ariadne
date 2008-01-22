/***************************************************************************
 *            taylor_derivative.code.h
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
 
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

namespace Ariadne {

template<class X>
const array< Function::TaylorVariable<X> >&
Function::TaylorDerivative<X>::variables() const
{
  return this->_variables;
}


template<class X>
LinearAlgebra::Vector<X>
Function::TaylorDerivative<X>::value() const
{
  size_type rs=this->result_size();
  LinearAlgebra::Vector<X> result(rs);
  for(uint i=0; i!=rs; ++i) {
    result[i]=this->_variables[i].value();
  }
  return result;
}


template<class X>
LinearAlgebra::Matrix<X>
Function::TaylorDerivative<X>::jacobian() const
{
  size_type rs=this->result_size();
  size_type as=this->argument_size();
  LinearAlgebra::Matrix<X> result(rs,as);
  for(uint i=0; i!=rs; ++i) {
    for(uint j=0; j!=as; ++j) {
      result[i][j]=this->_variables[i]._data[j+1u];
    }
  }
  return result;
}


template<class X>
void
Function::TaylorDerivative<X>::set_value(const LinearAlgebra::Vector<X>& v)
{
  ARIADNE_ASSERT(this->result_size()==v.size());
  size_type rs=this->result_size();
  for(uint i=0; i!=rs; ++i) {
    this->_variables[i]._data[0]=v[i];
  }
}


template<class X>
void
Function::TaylorDerivative<X>::set_jacobian(const LinearAlgebra::Matrix<X>& A)
{
  ARIADNE_ASSERT(this->result_size()==A.number_of_rows());
  ARIADNE_ASSERT(this->argument_size()==A.number_of_columns());
  size_type rs=this->result_size();
  size_type as=this->argument_size();
  for(uint i=0; i!=rs; ++i) {
    for(uint j=0; j!=as; ++j) {
      this->_variables[i]._data[j+1u]=A[i][j];
    }
  }
}


template<class X>
Function::TaylorDerivative<X>&
Function::operator-=(TaylorDerivative<X>& x, const LinearAlgebra::Vector<X>& v)
{
  ARIADNE_ASSERT(x.result_size()==v.size());
  for(size_type i=0; i!=v.size(); ++i) {
    x[i]-=v[i];
  }
  return x;
}

template<class X>
Function::TaylorDerivative<X>
Function::operator-(const TaylorDerivative<X>& x, const LinearAlgebra::Vector<X>& v)
{
  Function::TaylorDerivative<X> r(x);
  return r-=v;
}

template<class X>
Function::TaylorDerivative<X>
Function::operator*(const LinearAlgebra::Matrix<X>& A, const Function::TaylorDerivative<X>& x)
{
  ARIADNE_ASSERT(A.number_of_columns()==x.result_size());
  TaylorDerivative<X> y(A.number_of_rows(),x.argument_size(),x.degree());
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      y[i]+=A(i,j)*x[j];
    }
  }
  return y;
}


template<class X>
Function::TaylorDerivative<X>&
Function::operator*=(Function::TaylorDerivative<X>& x, const X& c)
{
  for(size_type i=0; i!=x.result_size(); ++i) {
    x[i]*=c;
  }
  return x;
}

template<class X>
Function::TaylorDerivative<X>
Function::operator*(const X& c, const Function::TaylorDerivative<X>& x)
{
  TaylorDerivative<X> y(x);
  return y*=c;
}

template<class X>
Function::TaylorDerivative<X>
Function::operator*(const Function::TaylorDerivative<X>& x, const X& c)
{
  TaylorDerivative<X> y(x);
  return y*=c;
}


template<class X> 
X
Function::evaluate(const TaylorVariable<X>& y, const LinearAlgebra::Vector<X>& x)
{
  return y.evaluate(x);
}

template<class X> 
Function::TaylorVariable<X>
Function::evaluate(const TaylorVariable<X>& y, const TaylorDerivative<X>& x)
{
  return Function::evaluate(y,x.variables());
}

template<class X> 
LinearAlgebra::Vector<X>
Function::evaluate(const TaylorDerivative<X>& y, const LinearAlgebra::Vector<X>& x)
{
  LinearAlgebra::Vector<X> r;
  r.data()=evaluate(y,x.data());
  return r;
}



template<class X> 
Function::TaylorVariable<X>
Function::compose(const TaylorVariable<X>& y, const TaylorDerivative<X>& x)
{
  TaylorDerivative<X> z(1u,y.argument_size(),y.degree());
  z[0]=y;
  return compose(z,x)[0];
}





template<class X> 
Function::TaylorDerivative<X>
Function::evaluate(const TaylorDerivative<X>& y, const TaylorDerivative<X>& x)
{
  using namespace std;
  ARIADNE_ASSERT(y.argument_size()==x.result_size());
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;

  size_type d=std::min(x.degree(),y.degree());
  size_type rs=y.result_size();
  size_type ms=x.result_size();
  size_type as=x.argument_size();
  
  TaylorDerivative<X> r(rs,as,d);
  TaylorVariable<X> t(as,d);

  // Use inefficient brute-force approach with lots of storage...
  array< array< TaylorVariable<X> > > val(ms, array< TaylorVariable<X> >(d+1));
  for(uint j=0; j!=ms; ++j) {
    val[j][0]=TaylorVariable<X>::constant(as,d,1.0);
    for(uint k=1; k<=d; ++k) {
      val[j][k]=val[j][k-1]*x[j];
    }
  }
  for(MultiIndex j(ms); j.degree()<=d; ++j) {
    t=TaylorVariable<X>::constant(as,d,1.0);
    for(uint k=0; k!=ms; ++k) {
      t=t*val[k][j[k]];
    }
    for(uint i=0; i!=rs; ++i) {
      r[i]+=y[i][j]*t;
    }
  }
  return r;
}


template<class X> 
Function::TaylorDerivative<X>
Function::compose(const TaylorDerivative<X>& y, const TaylorDerivative<X>& x)
{
  using namespace std;
  ARIADNE_ASSERT(y.argument_size()==x.result_size());
  size_type ms=x.result_size();
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;
  
  TaylorDerivative<X> w=x;
  for(uint i=0; i!=ms; ++i) {
    w[i].value()=0;
  }

  return evaluate(y,w);
}



template<class X> 
Function::TaylorDerivative<X> 
Function::concatenate(const TaylorDerivative<X>& x, const TaylorDerivative<X>& y)
{
  ARIADNE_ASSERT(x.argument_size()==y.argument_size());
  ARIADNE_ASSERT(x.degree()==y.degree());
  Function::TaylorDerivative<X> r(x.result_size()+y.result_size(),x.argument_size(),x.degree());
  for(size_type i=0; i!=x.result_size(); ++i) {
    r[i]=x[i];
  }
  for(size_type i=0; i!=y.result_size(); ++i) {
    r[i+x.result_size()]=y[i];
  }
  return r;
}

template<class X> 
Function::TaylorDerivative<X> 
Function::reduce(const TaylorDerivative<X>& x, const size_type& d)
{
  ARIADNE_ASSERT(x.degree()>=d);
  Function::TaylorDerivative<X> r(x.argument_size(),d);
  for(MultiIndex i(x.argument_size()); i.degree() <= x.degree(); ++i) {
    r[i]=x[i];
  }
}


template<class X> 
Function::TaylorDerivative<X> 
Function::derivative(const TaylorDerivative<X>& x, const size_type& k)
{
  TaylorDerivative<X> r( x.result_size(), x.argument_size(), ( x.degree()==0 ? 0u : x.degree()-1 ) );
  for(uint i=0; i!=x.result_size(); ++i) {
    r[i]=derivative(x[i],k);
  }
  return r;
}



template<class X> 
Function::TaylorDerivative<X> 
Function::inverse(const TaylorDerivative<X>& x, const LinearAlgebra::Vector<X>& c)
{
  using namespace std;
  ARIADNE_ASSERT(x.result_size()==x.argument_size());
  //std::cerr << "x=" << x << std::endl;
  size_type n=x.result_size();
  smoothness_type d=x.degree();
  LinearAlgebra::Vector<X> z(n,0);
  LinearAlgebra::Matrix<X> J=inverse(x.jacobian());

  TaylorDerivative<X> y(n,n,d);
  TaylorDerivative<X> id=TaylorDerivative<X>::variable(n,n,d,z);
  
  y.set_value(c);
  y.set_jacobian(J);
  for(smoothness_type i=2; i<=d; ++i) {
    y=y-J*(compose(x,y)-id);
  }
  return y;
}

template<class X> 
Function::TaylorDerivative<X> 
Function::implicit(const TaylorDerivative<X>& x, const LinearAlgebra::Vector<X>& c)
{
  using namespace std;
  using namespace LinearAlgebra;
  ARIADNE_ASSERT(x.result_size()<=x.argument_size());
  ARIADNE_ASSERT(c.size()==x.result_size());
  //std::cerr << "x=" << x << std::endl;
  
  size_type rs=x.result_size();
  size_type as=x.argument_size()-x.result_size();
  smoothness_type d=x.degree();

  LinearAlgebra::Matrix<X> A1(as,rs);
  for(size_type i=0; i!=as; ++i) {
    for(size_type j=0; j!=rs; ++j) {
      A1(i,j)=x[i].data()[1u+j];
    }
  }
  
  LinearAlgebra::Matrix<X> A2(rs,rs);
  for(size_type i=0; i!=rs; ++i) {
    for(size_type j=0; j!=rs; ++j) {
      A2(i,j)=x[i].data()[1u+as+j];
    }
  }
  
  Matrix<X> J(as+rs,rs);
  J(slice(as,rs),slice(0,rs))=inverse(A2);

  TaylorDerivative<X> y(as+rs,as,d);
  for(size_type i=0; i!=as; ++i) {
    y[i]=TaylorVariable<X>::variable(as,d,1.0,i);
  }
  for(size_type i=0; i!=rs; ++i) {
    // y[as+i]=TaylorVariable<X>::constant(as,d,0.0);
  }

  for(smoothness_type i=0; i!=d; ++i) {
    TaylorDerivative<X> z=compose(x,y);
    y=y-J*z;
  }

  TaylorDerivative<X> r(rs,as,d);
  for(size_type i=0; i!=rs; ++i) {
    r[i]=y[as+i];
  }
  return r;
}



template<class X> 
array< Function::TaylorSeries< Function::TaylorVariable<X> > >
Function::integrate(const array< TaylorVariable<X> >& f, const array< TaylorVariable<X> >& x)
{
  size_type n=x.size();
  smoothness_type d=std::max(f[0].degree(),x[0].degree());

  array< TaylorSeries< TaylorVariable<X> > > y = x;
  array< TaylorSeries< TaylorVariable<X> > > yp(n);
  for(uint j=0; j<d; ++j) {
    yp=compose(f,y);
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
  } 
  return y;
}



template<class X> 
std::ostream& 
Function::operator<<(std::ostream& os, const TaylorDerivative<X>& x) {
  //  return os << "TaylorDerivative( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
  if(x.result_size()==1) { return os << "[" << x[0] << "]"; }
  for(size_type i=0; i!=x.result_size(); ++i) {
    if(i==0) { os << "\n["; } else { os << ",\n "; }
    size_type degree=0;
    for(MultiIndex j(x.argument_size()); j.degree()<=x.degree(); ++j) {
      if(j.degree()==0) {
        os << "[ ";
      } else if(j.degree()==degree) {
        os << ",";
      } else {
        degree=j.degree();
        os << "; ";
      }
      os << x.get(i,j);
    }
    os << " ]";
  }
  os << "]\n";
  return os;

//  return os << "TaylorDerivative( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
}


template<class X> 
void
Function::TaylorDerivative<X>::instantiate() 
{
  X* c=0;
  LinearAlgebra::Vector<X>* v=0;
  TaylorVariable<X>* tv=0;
  TaylorDerivative<X>* td=0;
  std::ostream* os = 0;

  operator-=(*td,*v);
  operator-(*td,*v);

  operator*=(*td,*c);
  operator*(*td,*c);
  operator*(*c,*td);

  evaluate(*tv,*td);
  compose(*tv,*td);

  evaluate(*td,*v);
  evaluate(*td,*td);
  compose(*td,*td);
  inverse(*td,*v);
  implicit(*td,*v);
  derivative(*td,0u);
  concatenate(*td,*td);

  operator<<(*os,*td);
}


}
