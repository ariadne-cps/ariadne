/***************************************************************************
 *            taylor_model.code.h
 *
 *  Copyright  2007  Alberto Casagrande
 *  pieter.collins@cwi.nl
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

#include "taylor_model.h"
#include "exceptions.h"

#include <string>
#include <sstream>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "base/stlio.h"
#include "base/array.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"

#include "function/position_index.h"
#include "function/sorted_index.h"
#include "function/multi_index.h"

#include "geometry/rectangle.h"

#include "evaluation/newton.h"

#include "output/logging.h"
#include "output/latexstream.h"

namespace Ariadne {

template<class R>
Function::TaylorModel<R>::TaylorModel() 
  : _result_size(0), 
    _argument_size(0), 
    _order(0),
    _smoothness(0),
    _data() 
{
}


template<class R>
Function::TaylorModel<R>::TaylorModel(const std::string& s) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Function::TaylorModel<R>::TaylorModel(const size_type& rs, const size_type& as, const size_type& d)
  : _result_size(rs), 
    _argument_size(as), 
    _order(d),
    _smoothness(d),
    _data(rs*Numeric::bin(d+as,as),static_cast<R>(0)) 
{
}


template<class R>
Function::TaylorModel<R>::TaylorModel(const size_type& rs, const size_type& as, const size_type& d, const size_type& s)
  : _result_size(rs), 
    _argument_size(as), 
    _order(d),
    _smoothness(s),
    _data(rs*Numeric::bin(d+as,as),static_cast<R>(0)) 
{
}


template<class R>
void
Function::TaylorModel<R>::resize(const size_type& rs, const size_type& as, const size_type& d, const size_type& s)
{
  this->_result_size=rs; 
  this->_argument_size=as; 
  this->_order=d;
  this->_smoothness=s;
  this->_data.resize(rs*Numeric::bin(d+as,as));
}


template<class R>
Function::TaylorModel<R>
Function::TaylorModel<R>::zero(const size_type& rs, const size_type& as)  
{
  return TaylorModel<R>(rs,as,0,0);
}


template<class R>
Function::TaylorModel<R>
Function::TaylorModel<R>::one(const size_type& as)  
{
  R o=1;
  return TaylorModel<R>(1u,as,0,0,&o);
}


template<class R>
Function::TaylorModel<R>
Function::TaylorModel<R>::constant(const size_type& as, const R& c)  
{
  return TaylorModel<R>(1u,as,0,0,&c);
}


template<class R>
bool
Function::TaylorModel<R>::operator==(const TaylorModel<R>& p2) const
{
  const TaylorModel<R>& p1=*this;
  if(p1._result_size!=p2._result_size || p1._argument_size!=p2._argument_size) {
    return false;
  }
  if(p1._order==p2._order) {
    return p1._data==p2._data;
  } else {
    size_type smin=std::min(p1._data.size(),p2._data.size());
    size_type smax=std::max(p1._data.size(),p2._data.size());
    for(size_type i=0; i!=smin; ++i) {
      if(p1._data[i]!=p2._data[i]) {
        return false;
      }
    } 
    if(p1._data.size()>p2._data.size()) {
      for(size_type i=smin; i!=smax; ++i) {
        if(p1._data[i]!=0) {
          return false;
        }
      } 
    } else {
      for(size_type i=smin; i!=smax; ++i) {
        if(p2._data[i]!=0) {
          return false;
        }
      } 
    }
    return true;
  }
}


template<class R>
bool
Function::TaylorModel<R>::operator!=(const TaylorModel<R>& p2) const
{
  return !(*this==p2);
}


template<class R>
size_type 
Function::TaylorModel<R>::argument_size() const
{ 
  return this->_argument_size; 
}

template<class R>
size_type 
Function::TaylorModel<R>::result_size() const 
{ 
  return this->_result_size;
}

template<class R>
size_type 
Function::TaylorModel<R>::order() const 
{
  return this->_order; 
}
      
template<class R>
size_type 
Function::TaylorModel<R>::smoothness() const 
{ 
  return (size_type) -1; 
}
      
template<class R>
const array<R>&
Function::TaylorModel<R>::data() const 
{ 
  return this->_data;
}
      

template<class R>
void
Function::TaylorModel<R>::set(const size_type& i, const MultiIndex& j, const R& x) 
{ 
  assert(i<this->result_size());
  assert(j.degree()<=this->order());
  this->_data[this->result_size()*j.position()+i]=x;
}
      
template<class R>
R&
Function::TaylorModel<R>::at(const size_type& i, const MultiIndex& j)
{ 
  assert(i<this->result_size());
  assert(j.degree()<=this->order());
  return this->_data[this->result_size()*j.position()+i];
}
      
template<class R>
const R&
Function::TaylorModel<R>::get(const size_type& i, const MultiIndex& j) const 
{ 
  assert(i<this->result_size());
  assert(j.degree()<=this->order());
  return this->_data[this->result_size()*j.position()+i];
}
      


template<class R>
Function::TaylorModel<R> 
Function::TaylorModel<R>::component(const size_type& i) const
{
  TaylorModel<R> result(1u,this->argument_size(),this->order());
  R* rptr=result._data.begin();
  R* eptr=result._data.end();
  const R* aptr=this->_data.begin()+i;
  const size_type& inc=this->result_size();
  while(rptr!=eptr) {
    *rptr=*aptr;
    rptr+=1;
    aptr+=inc;
  }
  return result;
}



template<class R>
Function::TaylorModel<typename Function::TaylorModel<R>::I> 
Function::TaylorModel<R>::truncate(const size_type& order, const size_type& smoothness, const Geometry::Rectangle<R>& domain) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
LinearAlgebra::Vector<typename Function::TaylorModel<R>::F> 
Function::TaylorModel<R>::evaluate(const LinearAlgebra::Vector<F>& x) const
{
  if(this->argument_size()!=x.size()) {
    ARIADNE_THROW(IncompatibleSizes,"TaylorModel::evaluate(Vector)","Incompatible argument size");
  }

  // TODO: Make this more efficient
  LinearAlgebra::Vector<F> result(this->result_size());
  for(MultiIndex j(this->argument_size()); j.degree()<=this->order(); ++j) {
    F xa=1;
    for(size_type k=0; k!=j.number_of_variables(); ++k) {
      xa*=Numeric::pow(x[k],int(j[k]));
    }
    for(size_type i=0; i!=this->result_size(); ++i) {
      result[i]+=this->get(i,j)*xa;
    }
  }
  return result;
}



template<class R0, class R1, class R2>
void
Function::add(TaylorModel<R0>& p0, const TaylorModel<R1>& p1, const TaylorModel<R2>& p2)
{
  if(p1.result_size()!=p2.result_size()) {
    ARIADNE_THROW(IncompatibleSizes,"add(TaylorModel,TaylorModel)","Incompatible result sizes");
  }
  if(p1.argument_size()!=p2.argument_size()) {
    ARIADNE_THROW(IncompatibleSizes,"add(TaylorModel,TaylorModel)","Incompatible argument sizes");
  }
  size_type rs=p1.result_size();
  size_type as=p1.argument_size();
  size_type d1=p1.order();
  size_type d2=p2.order();
  size_type d=std::max(d1,d2);
  size_type s1=p1.smoothness();
  size_type s2=p2.smoothness();
  size_type s=std::max(s1,s2);

  p0.resize(rs,as,d,s);
  size_type kmin=std::min(p1._data.size(),p2._data.size());
  size_type kmax=p0._data.size();
  for(size_type i=0; i!=kmin; ++i) {
    p0._data[i]=p1._data[i]+p2._data[i];
  }
  if(d1>d2) {
    for(size_type i=kmin; i!=kmax; ++i) {
      p0._data[i]=p1._data[i];
    }
  } else {
    for(size_type i=kmin; i!=kmax; ++i) {
      p0._data[i]=p2._data[i];
    }
  }
}


template<class R0, class R1, class R2>
void
Function::sub(TaylorModel<R0>& p0, const TaylorModel<R1>& p1, const TaylorModel<R2>& p2)
{
  if(p1.result_size()!=p2.result_size()) {
    ARIADNE_THROW(IncompatibleSizes,"add(TaylorModel,TaylorModel)","Incompatible result sizes");
  }
  if(p1.argument_size()!=p2.argument_size()) {
    ARIADNE_THROW(IncompatibleSizes,"add(TaylorModel,TaylorModel)","Incompatible argument sizes");
  }
  size_type rs=p1.result_size();
  size_type as=p1.argument_size();
  size_type d1=p1.order();
  size_type d2=p2.order();
  size_type d=std::max(d1,d2);
  size_type s1=p1.smoothness();
  size_type s2=p2.smoothness();
  size_type s=std::max(s1,s2);

  p0.resize(rs,as,d,s);
  size_type kmin=std::min(p1._data.size(),p2._data.size());
  size_type kmax=p0._data.size();
  for(size_type i=0; i!=kmin; ++i) {
    p0._data[i]=p1._data[i]-p2._data[i];
  }
  if(d1>d2) {
    for(size_type i=kmin; i!=kmax; ++i) {
      p0._data[i]=p1._data[i];
    }
  } else {
    for(size_type i=kmin; i!=kmax; ++i) {
      p0._data[i]=-p2._data[i];
    }
  }
}


template<class R0, class R1, class R2>
void
Function::mul(TaylorModel<R0>& p0, const TaylorModel<R1>& p1, const TaylorModel<R2>& p2)
{
  if(p1.result_size()!=1u) {
    ARIADNE_THROW(IncompatibleSizes,"mul(TaylorModel,TaylorModel)","p1.result_size()="<<p1.result_size());
  }
  if(p2.result_size()!=1u) {
    ARIADNE_THROW(IncompatibleSizes,"mul(TaylorModel,TaylorModel)","p2.result_size()="<<p2.result_size());
  }
  if(p1.argument_size()!=p2.argument_size()) {
    ARIADNE_THROW(IncompatibleSizes,"add(TaylorModel p1,TaylorModel p2)","p1.result_size()="<<p1.result_size()<<", p2.result_size()="<<p2.result_size());
  }

  size_type rs=1u;
  size_type as=p1.argument_size();
  size_type d1=p1.order();
  size_type d2=p2.order();
  size_type d=d1+d2;
  size_type s1=p1.smoothness();
  size_type s2=p2.smoothness();
  size_type s=std::max(s1,s2);

  p0.resize(rs,as,d,s);
  MultiIndex j0(as);
  for(MultiIndex j1(as); j1.degree()<=d1; ++j1) {
    const R1& x1=p1.get(0u,j1);
    for(MultiIndex j2(as); j2.degree()<=d2; ++j2) {
      const R2& x2=p2.get(0u,j2);
      j0=j1+j2;
      p0.at(0u,j0)+=x1*x2;
    }
  }
}



template<class R0,class R1>
void
Function::pow(TaylorModel<R0>& p0, const TaylorModel<R1>& p1, const unsigned int& n)
{
  assert(p1.result_size()==1);

  if(n==1) {
    p0=p1;
    return;
  }

  R0 one=1;
  p0=TaylorModel<R0>::one(p1.argument_size());
  if(n==0) {
    return;
  }

  TaylorModel<R0> tmp(p1);
  for(uint i=1; i<=n; i*=2) {
    if(i&n) {
      p0=tmp*p0;
    }
    tmp=tmp*tmp;
  }
}



template<class R0,class R1>
void
Function::scale(TaylorModel<R0>& p0, const R1& x1)
{
  for(size_type i=0; i!=p0._data.size(); ++i) {
    p0._data[i]*=x1;
  }
}


template<class R0, class R1, class R2>
void
Function::compose(TaylorModel<R0>& p0, const TaylorModel<R1>& p1, const TaylorModel<R2>& p2)
{
  // TODO: Improve this algorithm as it's critical!!
  if(p1.argument_size()!=p2.result_size()) {
    ARIADNE_THROW(IncompatibleSizes,"compose(TaylorModel p1,TaylorModel p2)","p1.argument_size()="<<p1.argument_size()<<", p2.result_size()="<<p2.result_size());
  }

  if(p1.order()==0) { 
    p0=static_cast< TaylorModel<R0> >(p1); 
    return;
  }

  p0.resize(p1.result_size(),p2.argument_size(),p1.order()*p2.order(),std::max(p1.smoothness(),p2.smoothness()));
  
  TaylorModel<R0>* all_powers=new TaylorModel<R0>[p2.result_size()*(p1.order()+1)];
  TaylorModel<R0>* powers[p2.result_size()];
  for(size_type i=0; i!=p2.result_size(); ++i) {
    powers[i]=all_powers+i*(p1.order()+1);
  }

  for(size_type i=0; i!=p2.result_size(); ++i) {
    powers[i][0]=TaylorModel<R0>::one(p2.argument_size());
    powers[i][1]=p2.component(i);
    if(p1.order()>=2) {
      powers[i][2]=Function::pow(powers[i][1],2);
    }
    for(size_type j=3; j<=p1.order(); ++j) {
      powers[i][j]=powers[i][2]*powers[i][j-2];
    }
  }
  
  TaylorModel<R0>* results=new TaylorModel<R0>[p1.result_size()];
  for(size_type i=0; i!=p1.result_size(); ++i) {
    results[i]=TaylorModel<R0>::zero(1u,p2.argument_size());
  }

  for(size_type i=0; i!=p1.result_size(); ++i) {
    for(MultiIndex j(p1.argument_size()); j.degree()<=p1.order(); ++j) {
      TaylorModel<R0> t=TaylorModel<R0>::constant(p2.argument_size(),p1.get(i,j));
      for(size_type k=0; k!=p1.argument_size(); ++k) {
        t=t*powers[k][j[k]];
      }
      results[i]=results[i]+t;
    }
  }
  
  for(size_type i=0; i!=p0.result_size(); ++i) {
    for(size_type j=0; j!=p0.data().size()/p0.result_size(); ++j) {
      p0._data[i+j*p0.result_size()]=results[i].data()[j];
    }
  }
  
  delete[] results;
  delete[] all_powers;
}


template<class R0,class R1>
void
Function::derivative(TaylorModel<R0>& p0, const TaylorModel<R1>& p1, const size_type& k)
{
  if(p1.smoothness()==0) {
    ARIADNE_THROW(std::runtime_error,"derivative(TaylorModel,uint)"," model has smoothness 0");
  }

  p0.resize(p1.result_size(),p1.argument_size(),p1.order()-1,p1.smoothness()-1);
  
  MultiIndex dj(p1.argument_size());

  for(size_type i=0; i!=p0.result_size(); ++i) {
    for(MultiIndex j(p1.argument_size()); j.degree()<=p1.order(); ++j) {
      if(j[k]!=0) {
        dj=j;
        dj.decrement_index(k);
        p0.at(i,dj)+=static_cast<int>(j[k])*p1.get(i,j);
      }
    }
  }
}

template<class R>
LinearAlgebra::Matrix<typename Function::TaylorModel<R>::F> 
Function::TaylorModel<R>::jacobian(const LinearAlgebra::Vector<F>& s) const
{
  typedef typename Function::TaylorModel<R>::F R0;
  typedef R R1;

  LinearAlgebra::Matrix<R0> J(this->result_size(),this->argument_size());
  array< array<R0> > powers=this->_powers(s);
  
  for(size_type j=0; j!=this->argument_size(); ++j) {
    for(MultiIndex m(this->argument_size()); m.degree()<=this->order(); ++m) {
      MultiIndex n=m;
      int c=n[j];
      if(c!=0) {
        n.decrement_index(j);
        R0 a=c;
        for(size_type k=0; k!=this->argument_size(); ++k) {
          a*=powers[k][n[k]];
        }
        for(size_type i=0; i!=this->result_size(); ++i) {
          J(i,j)+=a*this->get(i,m);
        }
      }
    }
  }
  
  return J;
}



template<class R1, class R2> 
Function::TaylorModel<typename Numeric::traits<R1>::arithmetic_type> 
Function::inverse(const TaylorModel<R1>& p, const LinearAlgebra::Vector<R2>& v)
{
  assert(p.result_size()==p.argument_size());
  assert(p.argument_size()==v.size());

  // The following are only to simplfy testing.
  assert(v==LinearAlgebra::Vector<R2>(v.size(),0));
  assert(p.evaluate(v)==LinearAlgebra::Vector<R2>(p.result_size(),0));
  typedef typename Numeric::traits<R1>::arithmetic_type F0;
  typedef typename Numeric::traits<F0>::number_type R0;
  typedef typename Numeric::traits<F0>::interval_type I0;

  ARIADNE_LOG(2,"inverse(TaylorModel p, Vector v)\n");
  ARIADNE_LOG(3,"  p="<<p<<"\n  v="<<v<<"\n");
  LinearAlgebra::Vector<R0> c=midpoint(LinearAlgebra::Vector<I0>(p.evaluate(v)));
  LinearAlgebra::Matrix<F0> J=p.jacobian(v);
  
  LinearAlgebra::Matrix<F0> invJ=inverse(J);

  // FIXME: Need to re-solve for image of centre. What should initial set be? Different code needed for Rational?
  LinearAlgebra::Vector<F0> invf=v;

  TaylorModel<F0> result(p.argument_size(),p.result_size(),p.order(),p.smoothness());

  for(MultiIndex m(p.result_size()); m.degree()<=p.order(); ++m) {
    if(m.degree()==0) {
      for(size_type i=0; i!=p.argument_size(); ++i) {
        result.at(i,m)=v(i);
      }
    } else if(m.degree()==1) {
      for(size_type i=0; i!=p.argument_size(); ++i) {
        result.at(i,m)=invJ(i,m.position()-1);
      }
    } else {
      // FIXME: Add code for higher indices
    }
  }
  return result;
}


template<class R> template<class RR> 
Base::array< Base::array<typename Function::TaylorModel<R>::F> >
Function::TaylorModel<R>::_powers(const LinearAlgebra::Vector<RR>& v) const
{
  array< array<F> > powers(this->argument_size(), array<F>(this->order()+1));
  for(size_type i=0; i!=this->argument_size(); ++i) {
    powers[i][0]=static_cast<F>(1);
    if(this->order()>=1) {
      powers[i][1]=v(i);
      if(this->order()>=2) {
        powers[i][2]=Numeric::pow(v(i),2);
        for(size_type j=3; j<=this->order(); ++j) {
          powers[i][j]=powers[i][2]*powers[i][j-2];
        }
      }
    }
  }
  return powers;
}


template<class R>
std::ostream&
Function::TaylorModel<R>::write(std::ostream& os) const 
{
  os << "TaylorModel(\n";
  for(uint i=0; i!=this->result_size(); ++i) {
    os << "  " << i << ": " << std::flush;
    bool first=true;
    for(MultiIndex j(this->argument_size()); j.degree()<=this->order(); ++j) {
      if(first) { first=false; } else { os << ", "; } 
      os << j << ":" << j.position() << ":"  << std::flush; os << this->get(i,j) << std::flush;
    }
    os << ";\n";
  }
  os << ")";
  return os;
}


template<class R>
std::ostream&
Function::operator<<(std::ostream& os, const TaylorModel<R>& p)
{
  return p.write(os);
}


template<class R>
Output::latexstream&
Output::operator<<(Output::latexstream& texs, const Function::TaylorModel<R>& p)
{
  using namespace Function;
  texs << "%TaylorModel\n";
  texs << "\\ensuremath{\n";
  texs << "\\left( \\begin{array}{c}\n";
  char var='x';
  for(size_type i=0; i!=p.result_size(); ++i) {
    bool first = true;
    if(i!=0) { texs << "\\\\"; }
    for(MultiIndex j(p.argument_size()); j.degree()<=p.order(); ++j) {
      const R& a=p.get(i,j);
      if(a!=0) {
        if(first) { first=false; }
        else { if(a>0) { texs << '+'; } }
        if(a==1) { if(j.degree()==0) { texs << a; } }
        else if(a==-1) { if(j.degree()==0) { texs << a; } else { texs << '-'; } }
        else { texs << a << ' '; }
        for(size_type k=0; k!=p.argument_size(); ++k) {
          if(j[k]!=0) {
            texs << var << "_{ " << k+1 << "}";
            if(j[k]!=1) {
              texs << "^{" << j[k] << "}";
            }
          }
          texs << " ";
        }
      } 
    }
    texs << "\n"; 
  }
  texs << "\\end{array}\\right)\n}\n";
  return texs;
}


template<class R>
void
Function::TaylorModel<R>::instantiate()
{
  typedef typename Numeric::traits<R>::arithmetic_type I;
  size_type* k=0;
  R* x=0;
  LinearAlgebra::Vector<R>* v=0;
  TaylorModel<R>* p=0;
  std::ostream* os = 0;
  Output::latexstream* texs = 0;

  Function::operator+(*p,*p);
  Function::operator-(*p,*p);
  Function::operator*(*p,*p);
  Function::operator*(*x,*p);
  Function::operator*(*p,*x);
  Function::operator/(*p,*x);
  Function::pow(*p,0u);
  Function::compose(*p,*p);
  Function::derivative(*p,*k);
  Function::inverse(*p,*v);
  *os << *p;
  *texs << *p;
}




}
