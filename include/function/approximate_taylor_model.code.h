/***************************************************************************
 *            approximate_taylor_model.code.h
 *
 *  Copyright 2008  Pieter Collins
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

#include "base/exceptions.h"
#include "numeric/approximate_float.h"
#include "numeric/interval.h"
#include "linear_algebra/slice.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/multi_index.h"
#include "differentiation/sparse_differential.h"
#include "function/exceptions.h"
#include "function/approximate_taylor_model.h"
#include "function/function_interface.h"
#include "function/function_model_concept.h"
#include "output/latexstream.h"
#include "output/logging.h"


namespace Ariadne {

template<class R>
struct ApproximateTaylorModel<R>::Data {
  Data() : _domain(), _centre(), _expansion() { }

  Data(const Vector<I>& domain, const Vector<R>& centre,
       const SparseDifferentialVector<A>& expansion) 
    : _domain(domain), _centre(centre), _expansion(expansion) { }

  // Domain of definition.
  Ariadne::Vector<I> _domain;
  // The centre of the derivative expansion.
  Ariadne::Vector<A> _centre;
  // The derivatives of the model
  SparseDifferentialVector<A> _expansion;
};


template<class R>
ApproximateTaylorModel<R>::~ApproximateTaylorModel() 
{
  delete this->_data;
}


template<class R>
ApproximateTaylorModel<R>::ApproximateTaylorModel() 
  : _data(new Data)
{
}


template<class R>
ApproximateTaylorModel<R>::ApproximateTaylorModel(uint rs, uint as, ushort o, ushort s) 
  : _data(new Data(Vector<I>(as,I(-1,1)),Vector<A>(as),SparseDifferentialVector<A>(rs,as,o)))
{
}


template<class R>
ApproximateTaylorModel<R>::ApproximateTaylorModel(const Vector<I>& d, 
                                                  const Vector<R>& c, 
                                                  const SparseDifferentialVector<A>& e)
  : _data(new Data(d,c,e)) 
{
  ARIADNE_ASSERT(d.size()==e.argument_size());
  ARIADNE_ASSERT(c.size()==e.argument_size());
  for(uint i=0; i!=e.size(); ++i) {
    ARIADNE_ASSERT(e[i].argument_size()==e[0].argument_size());
  }
}




template<class R>
ApproximateTaylorModel<R>::ApproximateTaylorModel(const Vector<I>& d, const Vector<R>& c,
                                                  const FunctionInterface<R>& f, 
                                                  ushort o, ushort s)
  : _data(new Data(d,c,f.expansion(Vector<A>(c),o)))
{
  ARIADNE_ASSERT(d.size()==f.argument_size());
  ARIADNE_ASSERT(c.size()==f.argument_size());
}


template<class R>
ApproximateTaylorModel<R>::ApproximateTaylorModel(const Vector<I>& d, 
                                                  const FunctionInterface<R>& f, 
                                                  ushort o, ushort s)
  : _data(new Data(d,midpoint(d),f.expansion(Vector<A>(midpoint(d)),o)))
{
  ARIADNE_ASSERT(d.size()==f.argument_size());
}

template<class R>
ApproximateTaylorModel<R>::ApproximateTaylorModel(const ApproximateTaylorModel<R>& atm)
  : _data(new Data(*atm._data))
{
}


template<class R>
ApproximateTaylorModel<R>&
ApproximateTaylorModel<R>::operator=(const ApproximateTaylorModel<R>& atm)
{
  if(this!=&atm) {
    *this->_data=*atm._data;
  }
  return *this;
}




template<class R>
ApproximateTaylorModel<R>
ApproximateTaylorModel<R>::identity(const Vector<I>& domain, const Vector<R>& centre,
                                    ushort order, ushort smoothness)
{
  SparseDifferentialVector<A> expansion=SparseDifferentialVector<A>::variable(centre.size(),centre.size(),order,centre);
  return ApproximateTaylorModel<R>(domain,centre,expansion);
}


template<class R>
ApproximateTaylorModel<R>
ApproximateTaylorModel<R>::constant(const Vector<I>& domain, const Vector<R>& centre,
                                    const Vector<A>& value,
                                    ushort order, ushort smoothness)
{
  SparseDifferentialVector<A> expansion=SparseDifferentialVector<A>::constant(value.size(),domain.size(),order,value);
  return ApproximateTaylorModel<R>(domain,centre,expansion);
}


template<class R>
ApproximateTaylorModel<R>
ApproximateTaylorModel<R>::affine(const Vector<I>& domain, const Vector<R>& centre,
                                  const Vector<A>& value, const Matrix<A>& jacobian,
                                  ushort order, ushort smoothness)
{
  SparseDifferentialVector<A> expansion=SparseDifferentialVector<A>::affine(value.size(),domain.size(),order,value,jacobian);
  return ApproximateTaylorModel<R>(domain,centre,expansion);
}


template<class R>
ApproximateTaylorModel<R>
ApproximateTaylorModel<R>::affine(const I& domain, const R& centre,
                                  const A& value, const A& derivative,
                                  ushort order, ushort smoothness)
{
  SparseDifferentialVector<A> expansion(1u,1u,order);
  expansion[0].set_value(value);
  expansion[0].set_gradient(0,derivative);
  return ApproximateTaylorModel<R>(Vector<I>(1u,domain),Vector<R>(1u,centre),expansion);
}




template<class R>
Vector< Interval<R> >
ApproximateTaylorModel<R>::domain() const
{ 
  return this->_data->_domain; 
}


template<class R>
//Vector<typename ApproximateTaylorModel<R>::A>
Vector<R>
ApproximateTaylorModel<R>::centre() const
{ 
  Vector<R> r(this->_data->_centre.size()); 
  for(uint i=0; i!=this->_data->_centre.size(); ++i) {
    r[i]=this->_data->_centre[i]._value; 
  }
  return r;
}


template<class R>
Vector< Interval<R> >
ApproximateTaylorModel<R>::range() const
{ 
  return Ariadne::evaluate(this->expansion(),this->domain());
}



template<class R>
const SparseDifferentialVector<typename ApproximateTaylorModel<R>::A>&
ApproximateTaylorModel<R>::expansion() const
{ 
  return this->_data->_expansion; 
  //return *this->_expansion_ptr; 
}


template<class R>
std::pair< Vector<Interval<R> >, Matrix<R> >
affine_model(const ApproximateTaylorModel<R>& f)
{ 
  typedef typename traits<R>::interval_type I;
  typedef typename traits<R>::approximate_arithmetic_type A;
  

  uint rs=f.result_size();
  uint as=f.argument_size();
  const SparseDifferentialVector<A>& fe=f.expansion();

  // FIXME: Currently only assume domain is [-1,1], centre is 0
  //std::cerr<<"domain="<<f.domain()<<" centre="<<f.centre()<<std::endl;
  ARIADNE_ASSERT(f.domain()==Vector<I>(as,I(-1,1)));
  ARIADNE_ASSERT(f.centre()==Vector<R>(as,R(0)));

  Vector<I> value_range(rs);
  Matrix<R> jacobian(rs,as);
  for(uint i=0; i!=rs; ++i) {
    value_range[i]=fe[i].value()._value;
    for(uint j=0; j!=as; ++j) {
      jacobian[i][j]=R(fe[i].gradient(j)._value);
    }
    for(typename SparseDifferential<A>::const_iterator iter=fe[i].begin();
        iter!=fe[i].end(); ++iter)
    {
      if(iter->first.degree()>=2) {
        // FIXME: Change domain
        value_range[i]+=I(-1,1)*iter->second;
      }
    }
  }
  return std::make_pair(value_range,jacobian);
}


template<class R>
uint 
ApproximateTaylorModel<R>::argument_size() const
{ 
  return this->_data->_expansion.argument_size(); 
}


template<class R>
uint 
ApproximateTaylorModel<R>::result_size() const 
{ 
  return this->_data->_expansion.result_size();
}


template<class R>
ushort 
ApproximateTaylorModel<R>::order() const 
{
  return this->_data->_expansion.degree();
}
      

template<class R>
ushort 
ApproximateTaylorModel<R>::smoothness() const 
{ 
  return this->_data->_expansion.degree();
}
      

template<class R> template<class X>
array< array<X> >
ApproximateTaylorModel<R>::_powers(const Vector<X>& v) const
{
  ARIADNE_ASSERT(this->argument_size()==v.size());
  array< array<X> > powers(this->argument_size(), array<X>(this->order()+1));
  for(uint i=0; i!=this->argument_size(); ++i) {
    powers[i][0]=1;
    if(this->order()>=1) {
      powers[i][1]=v(i);
      if(this->order()>=2) {
        powers[i][2]=pow(v(i),2);
        for(uint j=3; j<=this->order(); ++j) {
          powers[i][j]=powers[i][2]*powers[i][j-2];
        }
      }
    }
  }
  return powers;
}


template<class R> 
Vector<typename ApproximateTaylorModel<R>::A> 
ApproximateTaylorModel<R>::evaluate(const Vector<A>& x) const
{
  ARIADNE_ASSERT(this->argument_size()==x.size());
  return Ariadne::evaluate(this->expansion(),x);
}


template<class R> 
Vector< Interval<R> > 
ApproximateTaylorModel<R>::evaluate(const Vector<I>& x) const
{
  ARIADNE_ASSERT(this->argument_size()==x.size());
  return Ariadne::evaluate(this->expansion(),x);
}




template<class R> 
Matrix<typename ApproximateTaylorModel<R>::A> 
ApproximateTaylorModel<R>::jacobian(const Vector<A>& x) const
{
  ARIADNE_ASSERT(this->argument_size()==x.size());
  SparseDifferentialVector<A> affine=SparseDifferentialVector<A>::variable(x.size(),x.size(),1u,x);
  affine=Ariadne::evaluate(this->expansion(),affine);
  return affine.jacobian();
}


template<class R> 
Matrix< Interval<R> > 
ApproximateTaylorModel<R>::jacobian(const Vector<I>& x) const
{
  ARIADNE_ASSERT(this->argument_size()==x.size());
  SparseDifferentialVector<I> affine=SparseDifferentialVector<I>::variable(x.size(),x.size(),1u,x);
  Vector< SparseDifferential<I> >& vec=affine;
  vec=Ariadne::evaluate(this->expansion(),vec);
  return affine.jacobian();
}


template<class R> 
ApproximateTaylorModel<R> 
project(const ApproximateTaylorModel<R>& f, const Slice& slc)
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  ARIADNE_ASSERT(slc.stop()<=f.result_size());
  return ApproximateTaylorModel<R>(f.domain(),f.centre(),SparseDifferentialVector<A>(project(f.expansion(),slc)));
}

template<class R> 
ApproximateTaylorModel<R> 
embed(const ApproximateTaylorModel<R>& f, const Vector< Interval<R> >& new_domain, const Vector<R>& new_centre, uint start)
{
  return ApproximateTaylorModel<R>(new_domain,new_centre,embed(f.expansion(),new_centre.size(),start));
}

template<class R> 
ApproximateTaylorModel<R> 
join(const ApproximateTaylorModel<R>& f, const ApproximateTaylorModel<R>& g)
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  ARIADNE_ASSERT(f.domain()==g.domain());
  ARIADNE_ASSERT(f.centre()==g.centre());

  return ApproximateTaylorModel<R> (f.domain(),f.centre(),join(f.expansion(),g.expansion()));
}

template<class R> 
ApproximateTaylorModel<R> 
combine(const ApproximateTaylorModel<R>& f, const ApproximateTaylorModel<R>& g)
{
  typedef typename traits<R>::interval_type I;
  typedef typename traits<R>::approximate_arithmetic_type A;
  uint fas=f.argument_size();
  uint gas=g.argument_size();
  Vector<I> hd=join(f.domain(),g.domain());
  Vector<R> hc=join(f.centre(),g.centre());
  SparseDifferentialVector<A> fex=embed(f.expansion(),fas+gas,0u);
  SparseDifferentialVector<A> gex=embed(g.expansion(),fas+gas,fas);
  SparseDifferentialVector<A> he=join(fex,gex);
  return ApproximateTaylorModel<R>(hd,hc,he);
}


/*

ApproximateTaylorModel<R> 
compose(const ApproximateTaylorModel<R>& f,
        const AffineTransformation<R>& g)
{
  ARIADNE_ASSERT(f.centre()==g.b());
  SparseDifferentialVector<R> const& fe=f.expansion();
  Vector<R> const& b=g.b();
  Matrix<R> const& A=g.A();

  SparseDifferentialVector<R> he;

  // Check that the matrix has zero or one unit in each row
  for(uint i=0; i!=A.row_size(); ++i) {
    uint ones=0;
    for(uint j=0; i!=A.column_size(); ++i) {
      if(A[i][j]==1) { ++ones; }
      else { ARIADNE_ASSERT(A[i][j]==0); }
    }
    ARIADNE_ASSERT(ones<=1);
  }
  
  ARIADNE_ASSERT(false); // not called
}

ApproximateTaylorModel<R> 
compose(const AffineTransformation<R>& f,
        const ApproximateTaylorModel<R>& g)
{
  ARIADNE_ASSERT(f.centre()==g.expansion().value());
  Vector<R> const& c=f.centre();
  Vector<R> const& b=f.b();
  Matrix<R> const& A=f.A();
  SparseDifferentialVector<R> const& ge=g.expansion();

  SparseDifferentialVector<R> he(f.result_size(),g.argument_size(),g.order());
  for(uint i=0; i!=f.result_size(); ++i) {
    he[i]+=b[i];
    for(uint j=0; j!=f.argument_size(); ++j) {
      he[i]+=A[i][j]*ge[j];
    }
  }
  return ApproximateTaylorModel<R>(g.domain(),g.centre(),he);
}

*/



template<class R> 
ApproximateTaylorModel<R> 
operator+(const ApproximateTaylorModel<R>& f,
          const ApproximateTaylorModel<R>& g)
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  ARIADNE_ASSERT(f.domain()==g.domain());
  ARIADNE_ASSERT(f.centre()==g.centre());
  ARIADNE_ASSERT(f.result_size()==g.result_size());
  SparseDifferentialVector<A> he;
  Vector< SparseDifferential<A> >& hev=he;
  hev=f.expansion(); hev+=g.expansion();
  return ApproximateTaylorModel<R>(f.domain(),f.centre(),he);
}


template<class R> 
ApproximateTaylorModel<R> 
operator-(const ApproximateTaylorModel<R>& f,
          const ApproximateTaylorModel<R>& g)
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  ARIADNE_ASSERT(f.domain()==g.domain());
  ARIADNE_ASSERT(f.centre()==g.centre());
  ARIADNE_ASSERT(f.result_size()==g.result_size());
  SparseDifferentialVector<A> he;
  Vector< SparseDifferential<A> >& hev=he;
  hev=f.expansion(); hev-=g.expansion();
  return ApproximateTaylorModel<R>(f.domain(),f.centre(),he);
}


template<class R> 
ApproximateTaylorModel<R> 
operator*(const ApproximateTaylorModel<R>& f,
          const ApproximateTaylorModel<R>& g)
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  ARIADNE_ASSERT(f.domain()==g.domain());
  ARIADNE_ASSERT(f.centre()==g.centre());
  ARIADNE_ASSERT(f.result_size()==1u || g.result_size()==1u);
  const SparseDifferential<A>* se;
  const Vector< SparseDifferential<A> >* ve;
  if(f.result_size()==1u) {
    se=&(f.expansion()[0]);
    ve=&g.expansion();
  } else {
    se=&(g.expansion()[0]);
    ve=&f.expansion();
  }
  SparseDifferentialVector<A> he(ve->size());
  Vector< SparseDifferential<A> >& hev=he;
  for(uint i=0; i!=ve->size(); ++i) {
    hev[i] = (*se) * ((*ve)[i]);
  }
  return ApproximateTaylorModel<R>(f.domain(),f.centre(),he);
}


template<class R> 
ApproximateTaylorModel<R> 
operator*(const ApproximateTaylorModel<R>& f,
          const R& c)
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  return ApproximateTaylorModel<R>(f.domain(),f.centre(),f.expansion()*A(c));
}


template<class R> 
ApproximateTaylorModel<R> 
compose(const ApproximateTaylorModel<R>& f,
        const ApproximateTaylorModel<R>& g)
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  ARIADNE_ASSERT(f.argument_size()==g.result_size());
  typedef Interval<R> I;
  typedef typename traits<R>::approximate_arithmetic_type A;
  SparseDifferentialVector<A> const& fe=f.expansion();
  //std::cerr<<"fe="<<fe<<std::endl;
  SparseDifferentialVector<A> const& ge=g.expansion();
  //std::cerr<<"ge="<<ge<<std::endl;
  Vector<A> tr=-Vector<A>(f.centre())+ge.value();
  //std::cerr<<"tr="<<tr<<std::endl;
  SparseDifferentialVector<A> fet=translate(fe,tr);
  //std::cerr<<"fet="<<fet<<std::endl;
  Vector<I> const& hd=g.domain();
  Vector<A> hc=g.centre();
  //std::cerr<<"he="<<std::flush;
  SparseDifferentialVector<A> he=compose(fet,ge);
  // FIXME: Change domain
  //std::cerr<<he<<std::endl;
  return ApproximateTaylorModel<R>(hd,hc,he);
}


template<class R> 
ApproximateTaylorModel<R> 
recentre(const ApproximateTaylorModel<R>& f,
         const Vector< Interval<R> >& d,
         const Vector<R>& c)
{
  ARIADNE_ASSERT(f.argument_size()==d.size());
  ARIADNE_ASSERT(f.argument_size()==c.size());
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R> 
ApproximateTaylorModel<R> 
derivative(const ApproximateTaylorModel<R>& f, uint k)
{
  return ApproximateTaylorModel<R>(f.domain(),f.centre(),derivative(f.expansion(),k));
}

template<class R> 
ApproximateTaylorModel<R> 
antiderivative(const ApproximateTaylorModel<R>& f, uint k)
{
  return ApproximateTaylorModel<R>(f.domain(),f.centre(),antiderivative(f.expansion(),k));
}

template<class R> 
ApproximateTaylorModel<R> 
inverse(const ApproximateTaylorModel<R>& f)
{
  return inverse(f,f.centre());
}


template<class R> 
ApproximateTaylorModel<R> 
inverse(const ApproximateTaylorModel<R>& f,
        const Vector<R>& x)
{
  ARIADNE_ASSERT(f.argument_size()==x.size());
  typedef typename traits<R>::interval_type I;
  typedef typename traits<R>::approximate_arithmetic_type A;
  Vector<A> tr=f.centre()-x;
  SparseDifferentialVector<A> fet=translate(f.expansion(),tr);
  Vector<I> gd=f.domain();
  Vector<A> gc=fet.value();
  
  SparseDifferentialVector<A> ge;
  try {
    ge=inverse(fet);
  } 
  catch(const SingularMatrixException&) {
    throw NonInvertibleFunctionException("inverse function model");
  }

  // FIXME: Change domain
  return ApproximateTaylorModel<R>(gd,gc,ge);
}


template<class R> 
ApproximateTaylorModel<R>
implicit(const ApproximateTaylorModel<R>& f)
{
  //std::cerr << __PRETTY_FUNCTION__ << std::endl;
  typedef typename traits<R>::interval_type I;
  typedef typename traits<R>::approximate_arithmetic_type A;
  ARIADNE_ASSERT(f.argument_size()>f.result_size());
  uint m=f.argument_size(); 
  uint n=f.result_size();

  array<uint> p(n);
  for(uint i=0; i!=n; ++i) { p[i]=i+m-n; }
  //std::cerr << "p=" << p << std::endl;

  // Construct the Taylor model for g(y)=f(c,y)
  Vector<I> gd=project(f.domain(),range(m-n,m));
  //std::cerr << "gd=" << gd << std::endl;
  Vector<A> gc=project(f.centre(),range(m-n,m));
  //std::cerr << "gc=" << gc << std::endl;
  SparseDifferentialVector<A> projection(m,n,f.order());
  for(uint i=0; i!=n; ++i) { projection[m-n+i][i]=1.0; }
  SparseDifferentialVector<A> ge=compose(f.expansion(),projection);
  //std::cerr << "ge=" << ge << std::endl;
  ApproximateTaylorModel<R> g(gd,gc,ge);

  Vector<R> z(n);
  Vector<I> iv = solve(g,z);
  Vector<A> v = midpoint(iv);
  //std::cerr<<"iv="<<iv<<std::endl;
  //std::cerr<<"v="<<v<<std::endl;
  Vector<A> t(m);
  project(t,range(m-n,m))=v;
  //std::cerr<<"t="<<t<<std::endl;

  SparseDifferentialVector<A> fe=f.expansion();
  //std::cerr<<"fe="<<fe<<std::endl;
  SparseDifferentialVector<A> fet=translate(fe,t);
  //std::cerr<<"fet="<<fet<<std::endl;
  fet.set_value(Vector<A>(fe.result_size(),0.0));
  //std::cerr<<"fet="<<fet<<std::endl;

  SparseDifferentialVector<A> he;
  try {
    he=implicit(fet);
  } 
  catch(const SingularMatrixException&) {
    ARIADNE_THROW(NonInvertibleFunctionException,__FUNCTION__,fet);
  }

  //std::cerr<<"he="<<he<<std::endl;
  he+=v;
  //std::cerr<<"he="<<he<<std::endl;
  Vector<I> hd=project(f.domain(),range(0,m-n));
  Vector<A> hc=project(f.centre(),range(0,m-n));
  return ApproximateTaylorModel<R>(hd,hc,he);
}

// If \f$f:\R^m\rightarrow \R^n\f$, compute the function
// \f$ h:\R^{m-n}\rightarrow \R^n\f$ satisfying 
// \f$ f(x,h(x)) = f(c) \f$ for all \f$x\f$.
// The function \f$ h\f$ is centred at
// \f$ c=(c_1,\ldots,c_{m-n}) \f$
// and has value \f$h(x)=(c_{m-n+1},\ldots,c_n)\f$.
template<class R> 
ApproximateTaylorModel<R>
implicit1(const ApproximateTaylorModel<R>& f,
          const Vector<R>& c)
{
  typedef Interval<R> I;
  uint m=f.argument_size(); 
  uint n=f.result_size();
  Vector<R> v=f.evaluate(c);
  Vector<R> hc=project(c,range(0,m-n));
  Vector<I> hd=project(f.domain(),range(0,m-n));
  
}

// If \f$f:\R^m\rightarrow \R^n\f$, compute the function
// \f$ h:\R^{m-n}\rightarrow \R^n\f$ satisfying 
// \f$ f(x,h(x)) = z \f$ for all \f$c\f$.
// The function \f$ h\f$ is centred at
// \f$ c=(c_1,\ldots,c_{m-n}) \f$.
//
// Precondition: there is a unique solution 
// of f(c,y)=z$.
template<class R> 
ApproximateTaylorModel<R>
implicit2(const ApproximateTaylorModel<R>& f,
          const Vector<R>& c,
          const Vector<R>& z)
{
  typedef Interval<R> I;
  // Solve for the value of h(c)
  uint m=f.argument_size(); 
  uint n=f.result_size();
  Vector<R> v=f.evaluate(c);
  Matrix<R> A=f.expansion().jacobian();
  //A=A[slice(0,n)][slice(m-n,n)];
  Matrix<R> B=inverse(A);
  Vector<R> y(n);
  Vector<R> tr=y-f.centre();
  
  SparseDifferentialVector<R> fet=translate(f.expansion(),tr)-z;
  SparseDifferentialVector<R> he=implicit(fet)+y;
  
  Vector<I> hd(m-n,I(-inf(),+inf()));
  Vector<R> hc=c;
  return ApproximateTaylorModel<R>(hd,hc,he);
}

template<class R> 
ApproximateTaylorModel<R> 
flow(const ApproximateTaylorModel<R>& vf)
{
  typedef Interval<R> I;
  typedef typename traits<R>::approximate_arithmetic_type A;
  uint ox=vf.order();
  uint ot=4u;

  Vector<A> vfc=vf.centre();
  SparseDifferentialVector<A> x=flow(vf.expansion(),vfc,ot,ox);

  /*
  uint n=vf.result_size();
  SparseDifferentialVector<A> x(n,n+1,ox);
  for(MultiIndex jx(n+1); jx.degree()<=ox; ++jx) {
    MultiIndex jy(n);
    for(uint k=0; k!=n; ++k) { 
      jy.set(k,jx[k]); 
    }
    uint jt=jx[n];
    for(uint i=0; i!=n; ++i) {
      std::cerr<<"i="<<i<<" jx="<<jx<<std::flush; std::cerr<<" x[i][jx]="<<x[i][jx]<<std::endl;
      std::cerr<<"i="<<i<<" jt="<<jt<<" jy="<<jy<<std::endl; 
      std::cerr<<"   y[i]]="<<y[i]<<std::endl;
      std::cerr<<"   y[i][jt]="<<y[i][jt]<<std::endl;
      std::cerr<<"   y[i][jt][jy]="<<y[i][jt][jy]<<std::endl;
      x[i][jx]=y[i][jt][jy];
    }
  }
  std::cerr << "x=" << x << std::endl;
  */

  I h(-0.1,0.1);
  Vector<I> d=join(vf.domain(),h);
  Vector<R> c=join(vf.centre(),R(0.0));
  return ApproximateTaylorModel<R>(d,c,x);
}


template<class R> 
ApproximateTaylorModel<R> 
integrate(const ApproximateTaylorModel<R>& vf, const R& h, uint ox)
{
  typedef Interval<R> I;
  ARIADNE_ASSERT(vf.result_size()==vf.argument_size());
  uint n=vf.result_size();
  SparseDifferentialVector<R> xe=vector_variable(n,ox,Vector<R>(n+1,0.0));
  SparseDifferential<R> te=scalar_constant(n,ox,h);
  ApproximateTaylorModel<R> f=flow(vf,ox);
  Vector<I> xtd=join(vf.domain(),I(-1,1));
  Vector<R> xtc=join(vf.centre(),0.0);
  ApproximateTaylorModel<R> xt(xtd,xtc,join(xe,te));
  return compose(f,xt);
}


// Find the set of values such that \f$f(x)=y\f$
// assuming a unique solution in 
template<class R> 
Vector< Interval<R> >
solve(const ApproximateTaylorModel<R>& f,
      const Vector<R>& y)
{
  ARIADNE_ASSERT(f.argument_size()==y.size());
  typedef typename traits<R>::approximate_arithmetic_type A;
  typedef Interval<R> I;
  Vector<I> x=f.centre();
  //Vector<I> x=f.domain();
  Matrix<A> J,Jinv;
  Vector<I> nx,fx,fxmy;
  Vector<A> m;
  for(uint i=0; i!=6; ++i) {
    //std::cerr << "  x[" << i << "]="<<x <<" y="<<y<<"\n";
    m=midpoint(x);
    J=f.jacobian(m);
    Jinv=inverse(J);
    fx=f.evaluate(x);
    fxmy=fx-y;
    nx=m-Vector<I>(Jinv*Vector<I>(fx-y));
    if(disjoint(x,nx)) {
      x=nx;
    } else {
      x=intersection(x,nx);
    }
    //std::cerr<<"    m="<<m<<" J="<<J<<" Jinv="<<Jinv<<" fx="<<fx<<" fx-y="<<fxmy<<" nx="<<nx<<"\n";
  }
  return x;
}



// Compute the hitting map of the flow of vf(x)
// with the guard condition g(x)=0
template<class R> 
ApproximateTaylorModel<R>
hitting(const ApproximateTaylorModel<R>& vf,
        const ApproximateTaylorModel<R>& g)
{
  typedef typename traits<R>::approximate_arithmetic_type A;
  typedef Interval<R> I;
  ARIADNE_ASSERT(vf.result_size()==vf.argument_size());
  ARIADNE_ASSERT(g.argument_size()==vf.result_size());
  ARIADNE_ASSERT(g.result_size()==1);
  ARIADNE_ASSERT(g.expansion().degree()==vf.expansion().degree());
  uint n=vf.result_size();
  uint d=vf.order();
  //std::cerr<<"vf="<<vf<<std::endl; // n,n
  //std::cerr<<"g="<<g<<std::endl; // 1,n+1
  //Differential<R> t=variable(n+1,d,0.0,n); // 1,n+1
  //std::cerr<<"t="<<t<<std::endl;
  ApproximateTaylorModel<R> f=flow(vf); // n,n+1
  //std::cerr<<"f="<<f<<std::endl;
  ApproximateTaylorModel<R> gf=compose(g,f); // 1,n+1
  //std::cerr<<"gf="<<gf<<std::endl;
  ApproximateTaylorModel<R> ht=implicit(gf); // 1,n
  //std::cerr<<"ht="<<ht<<std::endl;

  SparseDifferentialVector<A> xe=SparseDifferentialVector<A>::variable(n,n,d,vf.centre());
  //std::cerr<<"  xe="<<xe<<std::endl;
  ApproximateTaylorModel<R> hxt(ht.domain(),ht.centre(),join(xe,ht.expansion())); // n+1,n
  //std::cerr<<"hxt="<<hxt<<std::endl;
  
  ApproximateTaylorModel<R> hy=compose(f,hxt); // n,n
  //std::cerr<<"hy="<<hy<<std::endl;
  return join(hy,ht);
}



template<class R> 
ApproximateTaylorModel<R> 
ApproximateTaylorModel<R>::truncate(const Vector<I>& domain, const Vector<R>& centre,
                                    ushort order, ushort smoothness) const
{
  ARIADNE_ASSERT(this->argument_size()==domain.size());
  ARIADNE_ASSERT(this->argument_size()==centre.size());
  throw NotImplemented(__PRETTY_FUNCTION__);
}



template<class R> 
std::ostream&
operator<<(std::ostream& os, const ApproximateTaylorModel<R>& tm)
{
  os << "ApproximateTaylorModel(\n";
  os << "  argument_size=" << tm.expansion().argument_size() << ",\n" << std::flush;
  os << "  result_size=" << tm.expansion().result_size() << ",\n" << std::flush;
  os << "  degree=" << tm.expansion().degree() << ",\n" << std::flush;
  os << "  domain=" << tm.domain() << ",\n" << std::flush;
  os << "  centre=" << tm.centre() << ",\n" << std::flush;
  os << "  expansion=" << tm.expansion() << ",\n" << std::flush;
  os << ")\n";
  return os;
}


template<class R> 
latexstream&
operator<<(latexstream& texs, const ApproximateTaylorModel<R>& p)
{
  texs << "%ApproximateTaylorModel<R>\n";
  texs << "\\ensuremath{\n";
  texs << "\\left( \\begin{array}{c}\n";
  char var='x';
  for(uint i=0; i!=p.result_size(); ++i) {
    bool first = true;
    if(i!=0) { texs << "\\\\"; }
    for(MultiIndex j(p.argument_size()); j.degree()<=p.order(); ++j) {
      const R& a=p.expansion()[i][j];
      if(a!=0) {
        if(first) { first=false; }
        else { if(a>0) { texs << '+'; } }
        if(a==1) { if(j.degree()==0) { texs << a; } }
        else if(a==-1) { if(j.degree()==0) { texs << a; } else { texs << '-'; } }
        else { texs << a << ' '; }
        for(uint k=0; k!=p.argument_size(); ++k) {
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
ApproximateTaylorModel<R>::_instantiate()
{
  Slice* slc=0;
  Vector<R>* v=0;
  Vector< Interval<R> >* iv=0;
  ApproximateTaylorModel<R>* atm=0;
  std::ostream* os = 0;

  operator+(*atm,*atm);
  operator-(*atm,*atm);
  operator*(*atm,*atm);
  
  project(*atm,*slc);
  embed(*atm,*iv,*v,0u);
  combine(*atm,*atm);
  derivative(*atm,0u);
  antiderivative(*atm,0u);
  compose(*atm,*atm);
  inverse(*atm);
  implicit(*atm);
  flow(*atm);
  hitting(*atm,*atm);
  solve(*atm,*v);

  affine_model(*atm);
  
  *os << *atm;
}

} // namespace Ariadne
