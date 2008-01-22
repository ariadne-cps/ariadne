/***************************************************************************
 *            taylor_model.code.h
 *
 *  Copyright  2007  Pieter Collins
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

#include "geometry/box.h"

#include "evaluation/newton.h"

#include "output/logging.h"
#include "output/latexstream.h"

namespace Ariadne {

template<class R>
Function::TaylorModel<R>::TaylorModel() 
  : _domain(), 
    _centre(),
    _centre_derivatives(),
    _domain_derivatives() 
{
}

template<class R>
Function::TaylorModel<R>::TaylorModel(size_type rs, size_type as, smoothness_type o, smoothness_type s) 
  : _domain(Geometry::Box<R>::entire_space(as)),
    _centre(LinearAlgebra::Vector<R>(as)),
    _centre_derivatives(rs,as,o),
    _domain_derivatives(rs,as,o) 
{
}

template<class R>
Function::TaylorModel<R>::TaylorModel(const Geometry::Box<R>& d, const LinearAlgebra::Vector<R>& c, 
                                      const TaylorDerivative<I>& cd, const TaylorDerivative<I>& dd)
  : _domain(d),
    _centre(c),
    _centre_derivatives(cd),
    _domain_derivatives(dd)
{
}

template<class R>
Function::TaylorModel<R>::TaylorModel(const Geometry::Box<R>& d, const Geometry::Point<R>& c, 
                                      const TaylorDerivative<I>& cd, const TaylorDerivative<I>& dd)
  : _domain(d),
    _centre(c.position_vector()),
    _centre_derivatives(cd),
    _domain_derivatives(dd)
{
}

template<class R>
Function::TaylorModel<R>::TaylorModel(const LinearAlgebra::Vector<I>& d, const LinearAlgebra::Vector<R>& c,
                                      const TaylorDerivative<I>& cd, const TaylorDerivative<I>& dd)
  : _domain(d),
    _centre(c),
    _centre_derivatives(cd),
    _domain_derivatives(dd)
{
}

template<class R>
Function::TaylorModel<R>::TaylorModel(const Geometry::Box<R>& d, const Geometry::Point<R>& c, 
                                      smoothness_type o, smoothness_type s,
                                      const Function::FunctionInterface<R>& f)
  : _domain(d),
    _centre(c.position_vector()),
    _centre_derivatives(f.derivative(LinearAlgebra::Vector<I>(c.position_vector()),o)),
    _domain_derivatives(f.derivative(d.position_vectors(),o))
{
}

template<class R>
Function::TaylorModel<R>::TaylorModel(const LinearAlgebra::Vector<I>& d, const LinearAlgebra::Vector<R>& c, 
                                      smoothness_type o, smoothness_type s, 
                                      const Function::FunctionInterface<R>& f)
  : _domain(d),
    _centre(c),
    _centre_derivatives(f.derivative(LinearAlgebra::Vector<I>(c),o)),
    _domain_derivatives(f.derivative(d,o))
{
}

/*
template<class R>
void
Function::TaylorModel<R>::resize(size_type rs, size_type as, size_type d, size_type s)
{
  this->_smoothness=s;
  //this->_derivatives.resize(rs,as,d);
}
*/

template<class R>
Function::TaylorModel<R>
Function::TaylorModel<R>::zero(size_type rs, size_type as)  
{
  return TaylorModel<R>(rs,as,0,0);
}


template<class R>
Function::TaylorModel<R>
Function::TaylorModel<R>::one(size_type as)  
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Function::TaylorModel<R>
Function::TaylorModel<R>::constant(size_type as, const R& c)  
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
bool
Function::TaylorModel<R>::operator==(const TaylorModel<R>& tm) const
{
  return this->_centre==tm._centre
    && this->_centre_derivatives==tm._centre_derivatives;
}


template<class R>
bool
Function::TaylorModel<R>::operator!=(const TaylorModel<R>& p2) const
{
  return !(*this==p2);
}


template<class R>
Geometry::Box<R>
Function::TaylorModel<R>::domain() const
{ 
  return this->_domain; 
}

template<class R>
LinearAlgebra::Vector<R>
Function::TaylorModel<R>::centre() const
{ 
  return this->_centre; 
}

template<class R>
Geometry::Box<R>
Function::TaylorModel<R>::range() const
{ 
  return Geometry::Box<R>(this->_domain_derivatives.value());
}


template<class R>
const Function::TaylorDerivative<typename Function::TaylorModel<R>::I>&
Function::TaylorModel<R>::centre_derivatives() const
{ 
  return this->_centre_derivatives; 
}

template<class R>
const Function::TaylorDerivative<typename Function::TaylorModel<R>::I>&
Function::TaylorModel<R>::domain_derivatives() const
{ 
  return this->_domain_derivatives; 
}

template<class R>
size_type 
Function::TaylorModel<R>::argument_size() const
{ 
  return this->_centre_derivatives.argument_size(); 
}

template<class R>
size_type 
Function::TaylorModel<R>::result_size() const 
{ 
  return this->_centre_derivatives.result_size();
}

template<class R>
smoothness_type 
Function::TaylorModel<R>::order() const 
{
  return this->_centre_derivatives.degree();
}
      
template<class R>
smoothness_type 
Function::TaylorModel<R>::smoothness() const 
{ 
  return this->_domain_derivatives.degree();
}
      

















template<class R>
Function::TaylorModel<R> 
Function::TaylorModel<R>::truncate(const Geometry::Box<R>& domain, const LinearAlgebra::Vector<R>& centre,
                                   smoothness_type order, smoothness_type smoothness) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R>
LinearAlgebra::Vector<typename Function::TaylorModel<R>::I> 
Function::TaylorModel<R>::evaluate(const LinearAlgebra::Vector<R>& x) const
{
  return this->evaluate(LinearAlgebra::Vector<I>(x));
}


template<class R>
LinearAlgebra::Vector<typename Function::TaylorModel<R>::I> 
Function::TaylorModel<R>::evaluate(const LinearAlgebra::Vector<I>& x) const
{
  if(this->argument_size()!=x.size()) {
    ARIADNE_THROW(IncompatibleSizes,"TaylorModel::evaluate(Vector)","Incompatible argument size");
  }

  // TODO: Make this more efficient
  LinearAlgebra::Vector<I> w=x-this->centre();
  LinearAlgebra::Vector<I> result(this->result_size());
  for(MultiIndex j(this->argument_size()); j.degree()<this->order(); ++j) {
    I wa=1;
    for(size_type k=0; k!=j.number_of_variables(); ++k) {
      wa*=Numeric::pow(w[k],int(j[k]));
    }
    for(size_type i=0; i!=this->result_size(); ++i) {
      result[i]+=this->_centre_derivatives[i][j]*wa;
    }
  }
  for(MultiIndex j(this->argument_size(),this->order()); j.degree()<=this->order(); ++j) {
    I wa=1;
    for(size_type k=0; k!=j.number_of_variables(); ++k) {
      wa*=Numeric::pow(w[k],int(j[k]));
    }
    for(size_type i=0; i!=this->result_size(); ++i) {
      result[i]+=this->_domain_derivatives[i][j]*wa;
    }
  }
  return result;
}


template<class R>
Function::TaylorModel<R>
Function::recentre(const TaylorModel<R>& tm, const Geometry::Box<R>& bx, const LinearAlgebra::Vector<R>& c)
{
  ARIADNE_ASSERT(bx.subset(tm.domain()));
  ARIADNE_ASSERT(bx.contains(Geometry::Point<R>(c)));
  typedef Numeric::Interval<R> I;

  TaylorDerivative<I> translation=TaylorDerivative<I>::variable(tm.argument_size(),tm.argument_size(),tm.order(),c);
  
  //FIXME: This is incorrect...
  TaylorDerivative<I> new_centre_derivatives=compose(tm.centre_derivatives(),translation);
  TaylorDerivative<I> new_domain_derivatives=compose(tm.domain_derivatives(),translation);

  return TaylorModel<R>(bx,c,new_centre_derivatives,new_domain_derivatives);
}


template<class R>
Function::TaylorModel<R>
Function::add(const TaylorModel<R>& p1, const TaylorModel<R>& p2)
{
  ARIADNE_ASSERT(!open_intersection(p1.domain(),p2.domain()).empty());
  if(p1.centre()==p2.centre()) {
    return TaylorModel<R>(closed_intersection(p1.domain(),p2.domain()),
                          p1.centre(),
                          p1.centre_derivatives()+p2.centre_derivatives(),
                          p1.domain_derivatives()+p2.domain_derivatives());
  } else {
    Geometry::Box<R> new_domain=closed_intersection(p1.domain(),p2.domain());
    LinearAlgebra::Vector<R> new_centre=new_domain.centre().position_vector();
    return add(recentre(p1,new_domain,new_centre),recentre(p2,new_domain,new_centre));
  }
}

template<class R>
Function::TaylorModel<R>
Function::sub(const TaylorModel<R>& p1, const TaylorModel<R>& p2)
{
  ARIADNE_ASSERT(!open_intersection(p1.domain(),p2.domain()).empty());
  if(p1.centre()==p2.centre()) {
    return TaylorModel<R>(closed_intersection(p1.domain(),p2.domain()),
                          p1.centre(),
                          p1.centre_derivatives()-p2.centre_derivatives(),
                          p1.domain_derivatives()-p2.domain_derivatives());
  } else {
    Geometry::Box<R> new_domain=closed_intersection(p1.domain(),p2.domain());
    //Geometry::Point<R> new_centre=new_domain.centre();
    LinearAlgebra::Vector<R> new_centre=new_domain.centre().position_vector();
    return add(recentre(p1,new_domain,new_centre),recentre(p2,new_domain,new_centre));
  }
}



template<class R>
Function::TaylorModel<R>
Function::mul(const TaylorModel<R>& tm, const R& x)
{
  return TaylorModel<R>(tm.domain(),tm.centre(),tm.smoothness(),tm.derivatives()*x);
}


template<class R>
void
Function::mul(TaylorModel<R>& p0, const TaylorModel<R>& p1, const TaylorModel<R>& p2)
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
    const R& x1=p1.get(0u,j1);
    for(MultiIndex j2(as); j2.degree()<=d2; ++j2) {
      const R& x2=p2.get(0u,j2);
      j0=j1+j2;
      p0.at(0u,j0)+=x1*x2;
    }
  }
}



template<class R>
void
Function::pow(TaylorModel<R>& p0, const TaylorModel<R>& p1, const unsigned int& n)
{
  assert(p1.result_size()==1);

  if(n==1) {
    p0=p1;
    return;
  }

  R one=1;
  p0=TaylorModel<R>::one(p1.argument_size());
  if(n==0) {
    return;
  }

  TaylorModel<R> tmp(p1);
  for(uint i=1; i<=n; i*=2) {
    if(i&n) {
      p0=tmp*p0;
    }
    tmp=tmp*tmp;
  }
}



template<class R>
Function::TaylorModel<R>&
Function::operator*=(TaylorModel<R>& p0, const R& x1)
{
  for(size_type i=0; i!=p0._data.size(); ++i) {
    p0._data[i]*=x1;
  }
  return p0;
}


template<class R>
Function::TaylorModel<R>
Function::compose(const TaylorModel<R>& p2, const TaylorModel<R>& p1)
{
  typedef Numeric::Interval<R> I;
  ARIADNE_ASSERT(p2.centre()==p1.centre_derivatives().value());
  ARIADNE_ASSERT(subset(p1.range(),p2.domain()));
  Geometry::Box<R> new_domain=p1.domain();
  LinearAlgebra::Vector<R> new_centre=p1.centre();
  TaylorDerivative<I> new_centre_derivatives=compose(p1.centre_derivatives(),p2.centre_derivatives());
  TaylorDerivative<I> new_domain_derivatives=compose(p1.domain_derivatives(),p2.domain_derivatives());

  return TaylorModel<R>(new_domain,new_centre,new_centre_derivatives,new_domain_derivatives);
}


template<class R>
Function::TaylorModel<R>
Function::derivative(const TaylorModel<R>& tm, size_type k) 
{
  return TaylorModel<R>(tm.domain(),
                        tm.centre(),
                        derivative(tm.centre_derivatives(),k),
                        derivative(tm.domain_derivatives(),k));
}


template<class R>
void
Function::compose(TaylorModel<R>& p0, const TaylorModel<R>& p1, const TaylorModel<R>& p2)
{
  // TODO: Improve this algorithm as it's critical!!
  if(p1.argument_size()!=p2.result_size()) {
    ARIADNE_THROW(IncompatibleSizes,"compose(TaylorModel p1,TaylorModel p2)","p1.argument_size()="<<p1.argument_size()<<", p2.result_size()="<<p2.result_size());
  }

  if(p1.order()==0) { 
    p0=static_cast< TaylorModel<R> >(p1); 
    return;
  }

  p0.resize(p1.result_size(),p2.argument_size(),p1.order()*p2.order(),std::max(p1.smoothness(),p2.smoothness()));
  
  TaylorModel<R>* all_powers=new TaylorModel<R>[p2.result_size()*(p1.order()+1)];
  TaylorModel<R>* powers[p2.result_size()];
  for(size_type i=0; i!=p2.result_size(); ++i) {
    powers[i]=all_powers+i*(p1.order()+1);
  }

  for(size_type i=0; i!=p2.result_size(); ++i) {
    powers[i][0]=TaylorModel<R>::one(p2.argument_size());
    powers[i][1]=p2.component(i);
    if(p1.order()>=2) {
      powers[i][2]=Function::pow(powers[i][1],2);
    }
    for(size_type j=3; j<=p1.order(); ++j) {
      powers[i][j]=powers[i][2]*powers[i][j-2];
    }
  }
  
  TaylorModel<R>* results=new TaylorModel<R>[p1.result_size()];
  for(size_type i=0; i!=p1.result_size(); ++i) {
    results[i]=TaylorModel<R>::zero(1u,p2.argument_size());
  }

  for(size_type i=0; i!=p1.result_size(); ++i) {
    for(MultiIndex j(p1.argument_size()); j.degree()<=p1.order(); ++j) {
      TaylorModel<R> t=TaylorModel<R>::constant(p2.argument_size(),p1.get(i,j));
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


template<class R>
void
Function::derivative(TaylorModel<R>& p0, const TaylorModel<R>& p1, size_type k)
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
LinearAlgebra::Matrix<typename Function::TaylorModel<R>::I> 
Function::TaylorModel<R>::jacobian(const LinearAlgebra::Vector<R>& s) const
{
  return this->jacobian(LinearAlgebra::Vector<I>(s));
}

template<class R>
LinearAlgebra::Matrix<typename Function::TaylorModel<R>::I> 
Function::TaylorModel<R>::jacobian(const LinearAlgebra::Vector<I>& x) const
{
  LinearAlgebra::Matrix<I> J(this->result_size(),this->argument_size());
  LinearAlgebra::Vector<I> w=x-this->centre();
  array< array<I> > powers=this->_powers(w);

  for(size_type j=0; j!=this->argument_size(); ++j) {
    for(MultiIndex m(this->argument_size()); m.degree()<this->order(); ++m) {
      MultiIndex n=m;
      int c=n[j];
      if(c!=0) {
        n.decrement_index(j);
        I a=c;
        for(size_type k=0; k!=this->argument_size(); ++k) {
          a*=powers[k][n[k]];
        }
        for(size_type i=0; i!=this->result_size(); ++i) {
          J[i][j]+=a*this->_centre_derivatives[i][m];
        }
      }
    }
    for(MultiIndex m(this->argument_size(),this->order()); m.degree()<=this->order(); ++m) {
      MultiIndex n=m;
      int c=n[j];
      if(c!=0) {
        n.decrement_index(j);
        I a=c;
        for(size_type k=0; k!=this->argument_size(); ++k) {
          a*=powers[k][n[k]];
        }
        for(size_type i=0; i!=this->result_size(); ++i) {
          J[i][j]+=a*this->_domain_derivatives[i][m];
        }
      }
    }
  }
  
  return J;
}



template<class R> 
Function::TaylorModel<R> 
Function::inverse(const TaylorModel<R>& p, const LinearAlgebra::Vector<R>& v)
{
  assert(p.result_size()==p.argument_size());
  assert(p.argument_size()==v.size());

  // The following are only to simplfy testing.
  assert(v==LinearAlgebra::Vector<R>(v.size(),0));
  assert(p.evaluate(v)==LinearAlgebra::Vector<R>(p.result_size(),0));
  typedef typename Numeric::traits<R>::interval_type I;

  ARIADNE_LOG(2,"inverse(TaylorModel p, Vector v)\n");
  ARIADNE_LOG(3,"  p="<<p<<"\n  v="<<v<<"\n");
  LinearAlgebra::Vector<R> c=midpoint(LinearAlgebra::Vector<I>(p.evaluate(v)));
  LinearAlgebra::Matrix<I> J=p.jacobian(v);
  
  LinearAlgebra::Matrix<I> invJ=inverse(J);

  // FIXME: Need to re-solve for image of centre. What should initial set be? Different code needed for Rational?
  LinearAlgebra::Vector<I> invf=v;

  TaylorModel<R> result(p.argument_size(),p.result_size(),p.order(),p.smoothness());

  for(MultiIndex m(p.result_size()); m.degree()<=p.order(); ++m) {
    if(m.degree()==0) {
      for(size_type i=0; i!=p.argument_size(); ++i) {
        result._centre_derivatives[i][m]=v[i];
      }
    } else if(m.degree()==1) {
      for(size_type i=0; i!=p.argument_size(); ++i) {
        result._centre_derivatives[i][m]=invJ[i][m.position()-1];
      }
    } else {
      // FIXME: Add code for higher indices
    }
  }
  return result;
}


template<class R> 
Function::TaylorModel<R> 
Function::implicit(const TaylorModel<R>& p, const LinearAlgebra::Vector<R>& v)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R>
Base::array< Base::array<typename Function::TaylorModel<R>::I> >
Function::TaylorModel<R>::_powers(const LinearAlgebra::Vector<I>& v) const
{
  array< array<I> > powers(this->argument_size(), array<I>(this->order()+1));
  for(size_type i=0; i!=this->argument_size(); ++i) {
    powers[i][0]=1;
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
    os << "  domain=" << this->domain() << ",\n" << std::flush;
    os << "  centre=" << this->centre() << ",\n" << std::flush;
    os << "  range=" << this->range() << ",\n" << std::flush;
    os << "  expansion=" << this->_centre_derivatives << ",\n" << std::flush;
    os << "  bounds=" << this->_domain_derivatives << ",\n" << std::flush;
  }
  os << ")\n";
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
      const Numeric::Interval<R>& a=p.centre_derivatives()[i][j];
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
  //R* x=0;
  LinearAlgebra::Vector<R>* v=0;
  //Geometry::Point<R>* pt=0;
  Geometry::Box<R>* bx=0;
  TaylorModel<R>* tm=0;
  std::ostream* os = 0;
  Output::latexstream* texs = 0;

  Function::recentre(*tm,*bx,*v);
  Function::add(*tm,*tm);
  Function::sub(*tm,*tm);
  Function::compose(*tm,*tm);
  Function::derivative(*tm,*k);
  Function::inverse(*tm,*v);
  Function::implicit(*tm,*v);
  *os << *tm;
  *texs << *tm;
}




}
