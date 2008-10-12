/***************************************************************************
 *            file
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
 
/*! \file orbit.h
 *  \brief Orbits of dynamic systems
 */
#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "sparse_differential.h"
#include "function.h"
#include "approximate_taylor_model.h"


namespace Ariadne {

ApproximateTaylorModel implicit1(const ApproximateTaylorModel&, const Vector<Float>&);
ApproximateTaylorModel implicit2(const ApproximateTaylorModel&, const Vector<Float>&, const Vector<Float>&);

struct ApproximateTaylorModel::Data {
  Data() : _domain(), _centre(), _expansion() { }

  Data(const Vector<Interval>& domain, const Vector<Float>& centre,
       const DifferentialVector< SparseDifferential<Float> >& expansion) 
    : _domain(domain), _centre(centre), _expansion(expansion) { }

  // Domain of definition.
  Ariadne::Vector<Interval> _domain;
  // The centre of the derivative expansion.
  Ariadne::Vector<Float> _centre;
  // The smoothness to which the model can be computed
  uint _smoothness; 
  // The derivatives of the model
  DifferentialVector< SparseDifferential<Float> > _expansion;
};

ApproximateTaylorModel::~ApproximateTaylorModel() 
{
}


ApproximateTaylorModel::ApproximateTaylorModel() 
  : _data(new Data()) 
{
}


ApproximateTaylorModel::ApproximateTaylorModel(uint rs, uint as, ushort o, ushort s) 
  : _data(new Data(Vector<I>(as,I(-1,1)),Vector<R>(as),DifferentialVectorType(rs,as,o)))
{
}


ApproximateTaylorModel::ApproximateTaylorModel(const Vector<I>& d, 
                                               const Vector<A>& c, 
                                               const DifferentialVectorType& e)
  : _data(new Data(d,c,e)) 
{
  ARIADNE_ASSERT(d.size()==e.argument_size());
  ARIADNE_ASSERT(c.size()==e.argument_size());
  for(uint i=0; i!=e.size(); ++i) {
    ARIADNE_ASSERT(e[i].argument_size()==e[0].argument_size());
  }
}




ApproximateTaylorModel::ApproximateTaylorModel(const Vector<I>& d, const Vector<A>& c,
                                               const FunctionInterface& f, 
                                               ushort o, ushort s)
  : _data(new Data(d,c,f.expansion(Vector<A>(c),o)))
{
  ARIADNE_ASSERT(d.size()==f.argument_size());
  ARIADNE_ASSERT(c.size()==f.argument_size());
}




ApproximateTaylorModel::ApproximateTaylorModel(const Vector<I>& d, 
                                               const FunctionInterface& f, 
                                               ushort o, ushort s)
  : _data(new Data(d,midpoint(d),f.expansion(Vector<A>(midpoint(d)),o)))
{
  ARIADNE_ASSERT(d.size()==f.argument_size());
}

ApproximateTaylorModel::ApproximateTaylorModel(const ApproximateTaylorModel& atm)
  : _data(new Data(*atm._data))
{
}


ApproximateTaylorModel&
ApproximateTaylorModel::operator=(const ApproximateTaylorModel& atm)
{
  if(this!=&atm) { *this->_data=*atm._data; }
  return *this;
}




ApproximateTaylorModel
ApproximateTaylorModel::zero(uint rs, uint as)  
{
  return ApproximateTaylorModel(rs,as,0,0);
}



ApproximateTaylorModel
ApproximateTaylorModel::one(uint as)  
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



ApproximateTaylorModel
ApproximateTaylorModel::constant(uint as, const R& c)  
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


ApproximateTaylorModel
ApproximateTaylorModel::constant(const Vector<I>& d, const Vector<R>& c, const Vector<R>& x, ushort o, ushort s)  
{
  ApproximateTaylorModel result(x.size(),d.size(),o,s);
  result._data->_domain=d;
  result._data->_centre=c;
  result._data->_expansion=x;
  return result;
}

ApproximateTaylorModel
ApproximateTaylorModel::affine(const I& d, const R& c, const R& x, const R& g, ushort o, ushort s)  
{
  ApproximateTaylorModel result(1u,1u,o,s);
  result._data->_domain[0]=d;
  result._data->_centre[0]=c;
  result._data->_expansion[0].value()=x;
  result._data->_expansion[0].gradient(0)=g;
  return result;
}


Vector<Interval>
ApproximateTaylorModel::domain() const
{ 
  return this->_data->_domain; 
}


Vector<Float>
ApproximateTaylorModel::centre() const
{ 
  return this->_data->_centre;
}


Vector<Interval>
ApproximateTaylorModel::range() const
{ 
  return Ariadne::evaluate(this->expansion(),this->domain());
}



const ApproximateTaylorModel::DifferentialVectorType&
ApproximateTaylorModel::expansion() const
{ 
  return this->_data->_expansion; 
}


std::pair< Vector<Interval >, Matrix<Float> >
affine(const ApproximateTaylorModel& f)
{ 
  typedef Interval I;
  typedef Float A;

  uint rs=f.result_size();
  uint as=f.argument_size();
  const ApproximateTaylorModel::DifferentialVectorType& fe=f.expansion();
  Vector<I> value_range(rs);
  Matrix<A> jacobian(rs,as);
  for(uint i=0; i!=rs; ++i) {
    value_range[i]=fe[i].value();
    for(uint j=0; j!=as; ++j) {
      jacobian[i][j]=fe[i].gradient(j);
    }
    for(SparseDifferential<A>::const_iterator iter=fe[i].begin();
        iter!=fe[i].end(); ++iter)
    {
      if(iter->first.degree()>=2) {
        // FIXME: Change domain
        const A& value=iter->second;
        value_range[i]+=I(-1,1)*value;
      }
    }
  }
  return std::make_pair(value_range,jacobian);
}


uint 
ApproximateTaylorModel::argument_size() const
{ 
  return this->_data->_expansion.argument_size(); 
}


uint 
ApproximateTaylorModel::result_size() const 
{ 
  return this->_data->_expansion.result_size();
}


ushort 
ApproximateTaylorModel::order() const 
{
  return this->_data->_expansion.degree();
}
      

ushort 
ApproximateTaylorModel::smoothness() const 
{ 
  return this->_data->_expansion.degree();
}
      

template<class X>
array< array<X> >
ApproximateTaylorModel::_powers(const Vector<X>& v) const
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


Vector<Float> 
ApproximateTaylorModel::evaluate(const Vector<Float>& x) const
{
  ARIADNE_ASSERT(this->argument_size()==x.size());
  return Ariadne::evaluate(this->expansion(),x);
}


Vector<Interval> 
ApproximateTaylorModel::evaluate(const Vector<Interval>& x) const
{
  ARIADNE_ASSERT(this->argument_size()==x.size());
  return Ariadne::evaluate(this->expansion(),x);
}




Matrix<Float> 
ApproximateTaylorModel::jacobian(const Vector<Float>& x) const
{
  ARIADNE_ASSERT(this->argument_size()==x.size());
  DifferentialVectorType affine=DifferentialVectorType::variable(x.size(),x.size(),1u,x);
  affine=Ariadne::evaluate(this->expansion(),affine);
  return affine.jacobian();
}


Matrix<Interval> 
ApproximateTaylorModel::jacobian(const Vector<Interval>& x) const
{
  ARIADNE_ASSERT(this->argument_size()==x.size());
  IntervalDifferentialVectorType affine=IntervalDifferentialVectorType::variable(x.size(),x.size(),1u,x);
  Vector< SparseDifferential<I> >& vec=affine;
  vec=Ariadne::evaluate(this->expansion(),vec);
  return affine.jacobian();
}


ApproximateTaylorModel
ApproximateTaylorModel::identity(const Vector<I>& d)
{
  uint n=d.size();
  uint o=255;
  DifferentialVectorType ide(n,n,o);
  for(uint i=0; i!=n; ++i) {
    ide[i]=SparseDifferential<A>::variable(n,o,0.0,i);
  }
  return ApproximateTaylorModel(d,midpoint(d),ide);
}                     


ApproximateTaylorModel 
project_model(const ApproximateTaylorModel& f, const Slice& slc)
{
  //ARIADNE_ASSERT(slc.stop()<=f.result_size());
  return ApproximateTaylorModel(f.domain(),f.centre(),project(f.expansion(),slc));
}

ApproximateTaylorModel 
embed(const ApproximateTaylorModel& f, const Vector<Interval>& new_domain, const Vector<Float>& new_centre, uint start)
{
  return ApproximateTaylorModel(new_domain,new_centre,embed(f.expansion(),new_centre.size(),start));
}

ApproximateTaylorModel 
join(const ApproximateTaylorModel& f, const ApproximateTaylorModel& g)
{
  typedef Float A;
  ARIADNE_ASSERT(f.domain()==g.domain());
  ARIADNE_ASSERT(f.centre()==g.centre());

  ApproximateTaylorModel::DifferentialVectorType he(f.result_size()+g.result_size(),f.argument_size(),f.order());
  project(he,range(0,f.result_size()))=f.expansion();
  project(he,range(f.result_size(),f.result_size()+g.result_size()))=g.expansion();

  return ApproximateTaylorModel (f.domain(),f.centre(),he);
}


/*

ApproximateTaylorModel 
compose(const ApproximateTaylorModel& f,
        const AffineTransformation& g)
{
  ARIADNE_ASSERT(f.centre()==g.b());
  SparseDifferentialVector const& fe=f.expansion();
  Vector const& b=g.b();
  Matrix const& A=g.A();

  SparseDifferentialVector he;

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

ApproximateTaylorModel 
compose(const AffineTransformation& f,
        const ApproximateTaylorModel& g)
{
  ARIADNE_ASSERT(f.centre()==g.expansion().value());
  Vector const& c=f.centre();
  Vector const& b=f.b();
  Matrix const& A=f.A();
  SparseDifferentialVector const& ge=g.expansion();

  SparseDifferentialVector he(f.result_size(),g.argument_size(),g.order());
  for(uint i=0; i!=f.result_size(); ++i) {
    he[i]+=b[i];
    for(uint j=0; j!=f.argument_size(); ++j) {
      he[i]+=A[i][j]*ge[j];
    }
  }
  return ApproximateTaylorModel(g.domain(),g.centre(),he);
}

*/



ApproximateTaylorModel 
operator+(const ApproximateTaylorModel& f,
          const ApproximateTaylorModel& g)
{
  typedef Float A;
  ARIADNE_ASSERT(f.domain()==g.domain());
  ARIADNE_ASSERT(f.centre()==g.centre());
  ARIADNE_ASSERT(f.result_size()==g.result_size());
  ApproximateTaylorModel::DifferentialVectorType he;
  Vector< SparseDifferential<A> >& hev=he;
  hev=f.expansion(); hev+=g.expansion();
  return ApproximateTaylorModel(f.domain(),f.centre(),he);
}


ApproximateTaylorModel 
operator-(const ApproximateTaylorModel& f,
          const ApproximateTaylorModel& g)
{
  typedef Float A;
  ARIADNE_ASSERT(f.domain()==g.domain());
  ARIADNE_ASSERT(f.centre()==g.centre());
  ARIADNE_ASSERT(f.result_size()==g.result_size());
  ApproximateTaylorModel::DifferentialVectorType he;
  Vector< SparseDifferential<Float> >& hev=he;
  hev=f.expansion(); hev-=g.expansion();
  return ApproximateTaylorModel(f.domain(),f.centre(),he);
}



ApproximateTaylorModel
operator*(const ApproximateTaylorModel& f,
          const ApproximateTaylorModel& g)
{
  ARIADNE_ASSERT(f.domain()==g.domain());
  ARIADNE_ASSERT(f.centre()==g.centre());
  ARIADNE_ASSERT(f.result_size()==1u || g.result_size()==1u);
  const SparseDifferential<Float>* se;
  const Vector< SparseDifferential<Float> >* ve;
  if(f.result_size()==1u) {
    se=&(f.expansion()[0]);
    ve=&g.expansion();
  } else {
    se=&(g.expansion()[0]);
    ve=&f.expansion();
  }
  ApproximateTaylorModel::DifferentialVectorType he(ve->size());
  Vector< SparseDifferential<Float> >& hev=he;
  for(uint i=0; i!=ve->size(); ++i) {
    hev[i] = (*se) * ((*ve)[i]);
  }
  return ApproximateTaylorModel(f.domain(),f.centre(),he);
}


ApproximateTaylorModel 
compose(const ApproximateTaylorModel& f,
        const ApproximateTaylorModel& g)
{
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  ARIADNE_ASSERT(f.argument_size()==g.result_size());
  typedef Interval I;
  typedef Float A;
  typedef ApproximateTaylorModel::DifferentialVectorType DifferentialVectorType;

  DifferentialVectorType const& fe=f.expansion();
  //std::cerr<<"fe="<<fe<<std::endl;
  DifferentialVectorType const& ge=g.expansion();
  //std::cerr<<"ge="<<ge<<std::endl;
  Vector<A> tr=-Vector<A>(f.centre())+ge.value();
  //std::cerr<<"tr="<<tr<<std::endl;
  DifferentialVectorType fet=translate(fe,tr);
  //std::cerr<<"fet="<<fet<<std::endl;
  Vector<I> const& hd=g.domain();
  Vector<A> hc=g.centre();
  //std::cerr<<"he="<<std::flush;
  DifferentialVectorType he=compose(fet,ge);
  // FIXME: Change domain
  //std::cerr<<he<<std::endl;
  return ApproximateTaylorModel(hd,hc,he);
}


ApproximateTaylorModel 
recentre(const ApproximateTaylorModel& f,
         const Vector<Interval>& d,
         const Vector<Float>& c)
{
  ARIADNE_ASSERT(f.argument_size()==d.size());
  ARIADNE_ASSERT(f.argument_size()==c.size());
  throw NotImplemented(__PRETTY_FUNCTION__);
}


ApproximateTaylorModel 
inverse(const ApproximateTaylorModel& f,
        const Vector<Float>& x)
{
  ARIADNE_ASSERT(f.argument_size()==x.size());
  typedef Interval I;
  typedef Float A;
  Vector<A> tr=f.centre()-x;
  ApproximateTaylorModel::DifferentialVectorType fet=translate(f.expansion(),tr);
  Vector<I> gd=f.domain();
  Vector<A> gc=fet.value();
  ApproximateTaylorModel::DifferentialVectorType ge=inverse(fet);
  // FIXME: Change domain
  return ApproximateTaylorModel(gd,gc,ge);
}


ApproximateTaylorModel
implicit(const ApproximateTaylorModel& f)
{
  //std::cerr << __PRETTY_FUNCTION__ << std::endl;
  typedef Interval I;
  typedef Float A;
  typedef ApproximateTaylorModel::DifferentialVectorType DifferentialVectorType;
  
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
  DifferentialVectorType projection(m,n,f.order());
  for(uint i=0; i!=n; ++i) { projection[m-n+i][i]=1.0; }
  DifferentialVectorType ge=compose(f.expansion(),projection);
  //std::cerr << "ge=" << ge << std::endl;
  ApproximateTaylorModel g(gd,gc,ge);

  Vector<A> z(n);
  Vector<I> iv = solve(g,z);
  Vector<A> v = midpoint(iv);
  //std::cerr<<"iv="<<iv<<std::endl;
  //std::cerr<<"v="<<v<<std::endl;
  Vector<A> t(m);
  project(t,range(m-n,m))=v;
  //std::cerr<<"t="<<t<<std::endl;

  DifferentialVectorType fe=f.expansion();
  //std::cerr<<"fe="<<fe<<std::endl;
  DifferentialVectorType fet=translate(fe,t);
  //std::cerr<<"fet="<<fet<<std::endl;
  fet.set_value(Vector<A>(fe.result_size(),0.0));
  //std::cerr<<"fet="<<fet<<std::endl;

  ApproximateTaylorModel::DifferentialVectorType he=implicit(fet);
  //std::cerr<<"he="<<he<<std::endl;
  he+=v;
  //std::cerr<<"he="<<he<<std::endl;
  Vector<I> hd=project(f.domain(),range(0,m-n));
  Vector<A> hc=project(f.centre(),range(0,m-n));
  return ApproximateTaylorModel(hd,hc,he);
}

// If \f$f:\R^m\rightarrow \R^n\f$, compute the function
// \f$ h:\R^{m-n}\rightarrow \R^n\f$ satisfying 
// \f$ f(x,h(x)) = f(c) \f$ for all \f$x\f$.
// The function \f$ h\f$ is centred at
// \f$ c=(c_1,\ldots,c_{m-n}) \f$
// and has value \f$h(x)=(c_{m-n+1},\ldots,c_n)\f$.
ApproximateTaylorModel
implicit1(const ApproximateTaylorModel& f,
          const Vector<Float>& c)
{
  typedef Float A;
  typedef Interval I;
  uint m=f.argument_size(); 
  uint n=f.result_size();
  Vector<A> v=f.evaluate(c);
  Vector<A> hc=project(c,range(0,m-n));
  Vector<I> hd=project(f.domain(),range(0,m-n));
  throw NotImplemented(__PRETTY_FUNCTION__);
}

// If \f$f:\R^m\rightarrow \R^n\f$, compute the function
// \f$ h:\R^{m-n}\rightarrow \R^n\f$ satisfying 
// \f$ f(x,h(x)) = z \f$ for all \f$c\f$.
// The function \f$ h\f$ is centred at
// \f$ c=(c_1,\ldots,c_{m-n}) \f$.
//
// Precondition: there is a unique solution 
// of f(c,y)=z$.
ApproximateTaylorModel
implicit2(const ApproximateTaylorModel& f,
          const Vector<Float>& c,
          const Vector<Float>& z)
{
  typedef Float A;
  typedef Interval I;
  // Solve for the value of h(c)
  uint m=f.argument_size(); 
  uint n=f.result_size();
  Vector<A> v=f.evaluate(c);
  Matrix<A> a=f.expansion().jacobian();
  //A=A[slice(0,n)][slice(m-n,n)];
  Matrix<A> b=inverse(a);
  Vector<A> y(n);
  Vector<A> tr=y-f.centre();
  
  ApproximateTaylorModel::DifferentialVectorType fet=translate(f.expansion(),tr)-z;
  ApproximateTaylorModel::DifferentialVectorType he=implicit(fet)+y;
  
  Vector<I> hd(m-n,I(-inf(),+inf()));
  Vector<A> hc=c;
  return ApproximateTaylorModel(hd,hc,he);
}

ApproximateTaylorModel 
flow(const ApproximateTaylorModel& vf)
{
  typedef Interval I;
  typedef Float A;
  uint ox=vf.order();
  uint ot=4u;

  Vector<A> vfc=vf.centre();
  ApproximateTaylorModel::DifferentialVectorType x=flow(vf.expansion(),vfc,ot,ox);

  //std::cerr << "x=" << x << std::endl;

  /*
  uint n=vf.result_size();
  DifferentialVectorType x(n,n+1,ox);
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
  Vector<A> c=join(vf.centre(),A(0.0));
  return ApproximateTaylorModel(d,c,x);
}


ApproximateTaylorModel 
integrate(const ApproximateTaylorModel& vf, const Float& h, uint ox)
{
  typedef Float A;
  typedef Interval I;
  typedef ApproximateTaylorModel::DifferentialVectorType DifferentialVectorType;
  
  ARIADNE_ASSERT(vf.result_size()==vf.argument_size());
  uint n=vf.result_size();
  DifferentialVectorType xe=DifferentialVectorType::variable(n,n,ox,Vector<A>(n+1,0.0));
  SparseDifferential<A> te=SparseDifferential<A>::constant(n,ox,h);
  ApproximateTaylorModel f=flow(vf);
  Vector<I> xtd=join(vf.domain(),I(-1,1));
  Vector<A> xtc=join(vf.centre(),0.0);
  ApproximateTaylorModel xt(xtd,xtc,join(xe,te));
  return compose(f,xt);
}


// Find the set of values such that \f$f(x)=y\f$
// assuming a unique solution in 
Vector<Interval>
solve(const ApproximateTaylorModel& f,
      const Vector<Float>& y)
{
  ARIADNE_ASSERT(f.argument_size()==y.size());
  typedef Float A;
  typedef Interval I;
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
ApproximateTaylorModel
hitting(const ApproximateTaylorModel& vf,
        const ApproximateTaylorModel& g)
{
  typedef Float A;
  typedef Interval I;
  typedef ApproximateTaylorModel::DifferentialVectorType DifferentialVectorType;

  ARIADNE_ASSERT(vf.result_size()==vf.argument_size());
  ARIADNE_ASSERT(g.argument_size()==vf.result_size());
  ARIADNE_ASSERT(g.result_size()==1);
  ARIADNE_ASSERT(g.expansion().degree()==vf.expansion().degree());
  uint n=vf.result_size();
  uint d=vf.order();
  //std::cerr<<"vf="<<vf<<std::endl; // n,n
  //std::cerr<<"g="<<g<<std::endl; // 1,n+1
  //Differential t=variable(n+1,d,0.0,n); // 1,n+1
  //std::cerr<<"t="<<t<<std::endl;
  ApproximateTaylorModel f=flow(vf); // n,n+1
  //std::cerr<<"f="<<f<<std::endl;
  ApproximateTaylorModel gf=compose(g,f); // 1,n+1
  //std::cerr<<"gf="<<gf<<std::endl;
  ApproximateTaylorModel ht=implicit(gf); // 1,n
  //std::cerr<<"ht="<<ht<<std::endl;

  DifferentialVectorType xe=DifferentialVectorType::variable(n,n,d,vf.centre());
  //std::cerr<<"  xe="<<xe<<std::endl;
  ApproximateTaylorModel hxt(ht.domain(),ht.centre(),join(xe,ht.expansion())); // n+1,n
  //std::cerr<<"hxt="<<hxt<<std::endl;
  
  ApproximateTaylorModel hy=compose(f,hxt); // n,n
  //std::cerr<<"hy="<<hy<<std::endl;
  return join(hy,ht);
}



ApproximateTaylorModel 
ApproximateTaylorModel::truncate(const Vector<Interval>& domain, const Vector<Float>& centre,
                                 ushort order, ushort smoothness) const
{
  ARIADNE_ASSERT(this->argument_size()==domain.size());
  ARIADNE_ASSERT(this->argument_size()==centre.size());
  throw NotImplemented(__PRETTY_FUNCTION__);
}



std::ostream&
operator<<(std::ostream& os, const ApproximateTaylorModel& tm)
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


class latexstream : public std::ostream { };

latexstream&
operator<<(latexstream& texs, const ApproximateTaylorModel& p)
{
  typedef Float A;
  texs << "%ApproximateTaylorModel\n";
  texs << "\\ensuremath{\n";
  texs << "\\left( \\begin{array}{c}\n";
  char var='x';
  for(uint i=0; i!=p.result_size(); ++i) {
    bool first = true;
    if(i!=0) { texs << "\\\\"; }
    for(MultiIndex j(p.argument_size()); j.degree()<=p.order(); ++j) {
      const A& a=p.expansion()[i][j];
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



} // namespace Ariadne
