/***************************************************************************
 *            differential_vector.code.h
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
 
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

namespace Ariadne {

template<class X>
const array< Differential<X> >&
DifferentialVector<X>::variables() const
{
  return this->_variables;
}


template<class X>
Vector<X>
DifferentialVector<X>::value() const
{
  size_type rs=this->result_size();
  Vector<X> result(rs);
  for(uint i=0; i!=rs; ++i) {
    result[i]=this->_variables[i].value();
  }
  return result;
}


template<class X>
Matrix<X>
DifferentialVector<X>::jacobian() const
{
  size_type rs=this->result_size();
  size_type as=this->argument_size();
  Matrix<X> result(rs,as);
  for(uint i=0; i!=rs; ++i) {
    for(uint j=0; j!=as; ++j) {
      result[i][j]=this->_variables[i]._data[j+1u];
    }
  }
  return result;
}


template<class X>
void
DifferentialVector<X>::set_value(const Vector<X>& v)
{
  ARIADNE_ASSERT(this->result_size()==v.size());
  size_type rs=this->result_size();
  for(uint i=0; i!=rs; ++i) {
    this->_variables[i]._data[0]=v[i];
  }
}


template<class X>
void
DifferentialVector<X>::set_jacobian(const Matrix<X>& A)
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
DifferentialVector<X>&
operator-=(DifferentialVector<X>& x, const Vector<X>& v)
{
  ARIADNE_ASSERT(x.result_size()==v.size());
  for(size_type i=0; i!=v.size(); ++i) {
    x[i]-=v[i];
  }
  return x;
}

template<class X>
DifferentialVector<X>
operator-(const DifferentialVector<X>& x, const Vector<X>& v)
{
  DifferentialVector<X> r(x);
  return r-=v;
}

template<class X>
DifferentialVector<X>
operator*(const Matrix<X>& A, const DifferentialVector<X>& x)
{
  ARIADNE_ASSERT(A.number_of_columns()==x.result_size());
  DifferentialVector<X> y(A.number_of_rows(),x.argument_size(),x.degree());
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    for(size_type j=0; j!=A.number_of_columns(); ++j) {
      y[i]+=A(i,j)*x[j];
    }
  }
  return y;
}


template<class X>
DifferentialVector<X>&
operator*=(DifferentialVector<X>& x, const double& c)
{
  for(size_type i=0; i!=x.result_size(); ++i) {
    x[i]*=c;
  }
  return x;
}

template<class X>
DifferentialVector<X>&
operator*=(DifferentialVector<X>& x, const X& c)
{
  for(size_type i=0; i!=x.result_size(); ++i) {
    x[i]*=c;
  }
  return x;
}

template<class X>
DifferentialVector<X>
operator*(const X& c, const DifferentialVector<X>& x)
{
  DifferentialVector<X> y(x);
  return y*=c;
}

template<class X>
DifferentialVector<X>
operator*(const DifferentialVector<X>& x, const X& c)
{
  DifferentialVector<X> y(x);
  return y*=c;
}


template<class X> 
X
evaluate(const Differential<X>& y, const Vector<X>& x)
{
  return y.evaluate(x);
}

template<class X> 
Differential<X>
evaluate(const Differential<X>& y, const DifferentialVector<X>& x)
{
  return evaluate(y,x.variables());
}

template<class X> 
Vector<X>
evaluate(const DifferentialVector<X>& y, const Vector<X>& x)
{
  Vector<X> r;
  r.data()=evaluate(y,x.data());
  return r;
}



template<class X>
DifferentialVector<X> 
translate(const DifferentialVector<X>& x, const Vector<X>& c)
{
  assert(x.argument_size()==c.size());
  DifferentialVector<X> r(x);
  for(uint k=0; k!=x.argument_size(); ++k) {
    if(c[k]!=0) { 
      DifferentialVector<X> y(r);
      r*=0.0;
      MultiIndex e(y.argument_size());
      e.set(k,1);
      for(MultiIndex j(y.argument_size()); j.degree()<=y.degree(); ++j) {
        MultiIndex tj=j;
        for(uint l=0; l<=y.degree()-j.degree(); ++l) {
          assert(tj.degree()<=y.degree());
          for(uint i=0; i!=x.result_size(); ++i) {
            r[i][j]+=y[i][tj]*pow(c[k],l)*bin(tj[k],l);
          }
          tj+=e;
        }
      }
    }
  }
  return r;
}
 

template<class X> 
Differential<X>
compose(const Differential<X>& y, const DifferentialVector<X>& x)
{
  DifferentialVector<X> z(1u,y.argument_size(),y.degree());
  z[0]=y;
  return compose(z,x)[0];
}





template<class X> 
DifferentialVector<X>
evaluate(const DifferentialVector<X>& y, const DifferentialVector<X>& x)
{
  using namespace std;
  ARIADNE_ASSERT(y.argument_size()==x.result_size());
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;

  size_type d=std::min(x.degree(),y.degree());
  size_type rs=y.result_size();
  size_type ms=x.result_size();
  size_type as=x.argument_size();
  
  DifferentialVector<X> r(rs,as,d);
  Differential<X> t(as,d);

  // Use inefficient brute-force approach with lots of storage...
  array< array< Differential<X> > > val(ms, array< Differential<X> >(d+1));
  for(uint j=0; j!=ms; ++j) {
    val[j][0]=Differential<X>::constant(as,d,1.0);
    for(uint k=1; k<=d; ++k) {
      val[j][k]=val[j][k-1]*x[j];
    }
  }
  for(MultiIndex j(ms); j.degree()<=d; ++j) {
    t=Differential<X>::constant(as,d,1.0);
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
DifferentialVector<X>
compose(const DifferentialVector<X>& y, const DifferentialVector<X>& x)
{
  using namespace std;
  ARIADNE_ASSERT(y.argument_size()==x.result_size());
  size_type ms=x.result_size();
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;
  
  DifferentialVector<X> w=x;
  for(uint i=0; i!=ms; ++i) {
    w[i].value()=0;
  }

  return evaluate(y,w);
}



template<class X> 
DifferentialVector<X> 
concatenate(const DifferentialVector<X>& x, const DifferentialVector<X>& y)
{
  ARIADNE_ASSERT(x.argument_size()==y.argument_size());
  ARIADNE_ASSERT(x.degree()==y.degree());
  DifferentialVector<X> r(x.result_size()+y.result_size(),x.argument_size(),x.degree());
  for(size_type i=0; i!=x.result_size(); ++i) {
    r[i]=x[i];
  }
  for(size_type i=0; i!=y.result_size(); ++i) {
    r[i+x.result_size()]=y[i];
  }
  return r;
}

template<class X> 
DifferentialVector<X> 
reduce(const DifferentialVector<X>& x, const size_type& d)
{
  ARIADNE_ASSERT(x.degree()>=d);
  DifferentialVector<X> r(x.argument_size(),d);
  for(MultiIndex i(x.argument_size()); i.degree() <= x.degree(); ++i) {
    r[i]=x[i];
  }
}


template<class X> 
DifferentialVector<X> 
derivative(const DifferentialVector<X>& x, const size_type& k)
{
  DifferentialVector<X> r( x.result_size(), x.argument_size(), ( x.degree()==0 ? 0u : x.degree()-1 ) );
  for(uint i=0; i!=x.result_size(); ++i) {
    r[i]=derivative(x[i],k);
  }
  return r;
}


template<class X> 
DifferentialVector<X> 
antiderivative(const DifferentialVector<X>& x, const size_type& k)
{
  DifferentialVector<X> r( x.result_size(), x.argument_size(), x.degree()+1);
  for(uint i=0; i!=x.result_size(); ++i) {
    r[i]=antiderivative(x[i],k);
  }
  return r;
}



template<class X> 
DifferentialVector<X> 
inverse(const DifferentialVector<X>& x, const Vector<X>& c)
{
  using namespace std;
  ARIADNE_ASSERT(x.result_size()==x.argument_size());
  //std::cerr << "x=" << x << std::endl;
  size_type n=x.result_size();
  smoothness_type d=x.degree();
  Vector<X> z(n,0);
  Matrix<X> J=inverse(x.jacobian());

  DifferentialVector<X> y(n,n,d);
  DifferentialVector<X> id=DifferentialVector<X>::variable(n,n,d,z);
  
  y.set_value(c);
  y.set_jacobian(J);
  for(smoothness_type i=2; i<=d; ++i) {
    y=y-J*(compose(x,y)-id);
  }
  return y;
}

template<class X> 
DifferentialVector<X> 
implicit(const DifferentialVector<X>& x, const Vector<X>& c)
{
  using namespace std;
  
  ARIADNE_ASSERT(x.result_size()<=x.argument_size());
  ARIADNE_ASSERT(c.size()==x.result_size());
  //std::cerr << "x=" << x << std::endl;
  
  size_type rs=x.result_size();
  size_type as=x.argument_size()-x.result_size();
  smoothness_type d=x.degree();

  Matrix<X> A1(as,rs);
  for(size_type i=0; i!=as; ++i) {
    for(size_type j=0; j!=rs; ++j) {
      A1(i,j)=x[i].data()[1u+j];
    }
  }
  
  Matrix<X> A2(rs,rs);
  for(size_type i=0; i!=rs; ++i) {
    for(size_type j=0; j!=rs; ++j) {
      A2(i,j)=x[i].data()[1u+as+j];
    }
  }
  
  Matrix<X> J(as+rs,rs);
  J(slice(as,rs),slice(0,rs))=inverse(A2);

  DifferentialVector<X> y(as+rs,as,d);
  for(size_type i=0; i!=as; ++i) {
    y[i]=Differential<X>::variable(as,d,1.0,i);
  }
  for(size_type i=0; i!=rs; ++i) {
    // y[as+i]=Differential<X>::constant(as,d,0.0);
  }

  for(smoothness_type i=0; i!=d; ++i) {
    DifferentialVector<X> z=compose(x,y);
    y=y-J*z;
  }

  DifferentialVector<X> r(rs,as,d);
  for(size_type i=0; i!=rs; ++i) {
    r[i]=y[as+i];
  }
  return r;
}



template<class X> 
array< PowerSeries< Differential<X> > >
integrate(const array< Differential<X> >& f, const array< Differential<X> >& x)
{
  size_type n=x.size();
  smoothness_type d=std::max(f[0].degree(),x[0].degree());

  array< PowerSeries< Differential<X> > > y = x;
  array< PowerSeries< Differential<X> > > yp(n);
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
operator<<(std::ostream& os, const DifferentialVector<X>& x) {
  //  return os << "DifferentialVector( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
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

//  return os << "DifferentialVector( argument_size=" << x.argument_size() << ", degree=" << x.degree() << ", data=" << x.data() << ")";
}


template<class X> 
void
DifferentialVector<X>::instantiate() 
{
  X* c=0;
  Vector<X>* v=0;
  Differential<X>* tv=0;
  DifferentialVector<X>* td=0;
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
  translate(*td,*v);
  inverse(*td,*v);
  implicit(*td,*v);
  derivative(*td,0u);
  antiderivative(*td,0u);
  concatenate(*td,*td);

  operator<<(*os,*td);
}


}
