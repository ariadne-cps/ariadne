#ifndef ARIADNE_DIFFERENTIAL_VECTOR_H
#define ARIADNE_DIFFERENTIAL_VECTOR_H

#include <map>

#include "macros.h"
#include "array.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"

namespace Ariadne {

template<class X> class Vector;

/*! \brief A class representing the derivatives of a vector quantity depending on multiple arguments. */
template<class DIFF>
class DifferentialVector
  : public Vector< DIFF >
{
  //BOOST_CONCEPT_ASSERT((DifferentialVectorConcept<SparseDifferentialVector<X> >));
  typedef typename DIFF::ScalarType X;
 public:
  DifferentialVector() 
    : Vector<DIFF>(0,DIFF()) { }
  DifferentialVector(uint rs, uint as, uint d) 
    : Vector<DIFF>(rs,DIFF(as,d)) { }
  DifferentialVector(const Vector<DIFF>& vsd) 
    : Vector<DIFF>(vsd) { }
  template<class XX> DifferentialVector(uint rs, uint as, uint d, const XX* ptr) 
    : Vector<DIFF>(rs,DIFF(as,d)) 
  { 
    for(uint i=0; i!=rs; ++i) { for(MultiIndex j(as); j.degree()<=d; ++j) {
        if(*ptr!=0) { (*this)[i][j]=*ptr; } ++ptr; } } 
  }
  DifferentialVector(uint rs, uint as, uint d, 
                     const Vector<X>& v, const Matrix<X>& A)
    :  Vector<DIFF>(rs,DIFF()) {
    ARIADNE_ASSERT(rs==v.size());
    ARIADNE_ASSERT(rs==A.row_size());
    ARIADNE_ASSERT(as==A.column_size());
    for(uint i=0; i!=this->result_size(); ++i) { 
      (*this)[i]=v[i]; 
      for(uint j=0; j!=this->argument_size(); ++j) { 
        const X& x=A[i][j];
        if(x!=0) { (*this)[i][j]=x; } } } 
  }
  template<class E> DifferentialVector(const ublas::vector_expression<E>& ve) 
    : Vector<DIFF>(ve) { }

  uint result_size() const { return this->Vector<DIFF>::size(); }
  uint argument_size() const { return (*this)[0].argument_size(); }
  uint degree() const { return (*this)[0].degree(); }

  Vector<X> value() const { 
    Vector<X> r(this->result_size()); for(uint i=0; i!=r.size(); ++i) { r[i]=(*this)[i].value(); } return r; }
  Matrix<X> jacobian() const { Matrix<X> r(this->result_size(),this->argument_size()); 
    for(uint i=0; i!=r.row_size(); ++i) { for(uint j=0; j!=r.column_size(); ++j) { r[i][j]=(*this)[i].gradient(j); } } return r; }

  void set_value(const Vector<X>& c) {
    ARIADNE_ASSERT(this->result_size()==c.size());
    for(uint i=0; i!=c.size(); ++i) { (*this)[i].set_value(c[i]); } }

  static DifferentialVector<DIFF> constant(uint rs, uint as, uint d, const Vector<X>& c) {
    ARIADNE_ASSERT(c.size()==rs);
    DifferentialVector<DIFF> result(rs,as,d);
    for(uint i=0; i!=rs; ++i) { result[i]=c[i]; }
    return result;
  }

  static DifferentialVector<DIFF> variable(uint rs, uint as, uint d, const Vector<X>& x) {
    ARIADNE_ASSERT(x.size()==rs);
    DifferentialVector<DIFF> result(rs,as,d);
    for(uint i=0; i!=rs; ++i) { result[i]=x[i]; result[i][i]=X(1.0); }
    return result;
  }

};

template<class DIFF, class Y>
DifferentialVector<DIFF>&
operator+=(DifferentialVector<DIFF>& x, const Vector<Y>& c)
{  
  assert(x.result_size()==c.size());
  for(uint i=0; i!=c.size();++i) {
    x[i]+=c[i];
  }
  return x;
}

  
template<class DIFF>
DifferentialVector<DIFF>&
operator-=(DifferentialVector<DIFF>& x, const Vector<typename DIFF::ScalarType>& c)
{  
  assert(x.result_size()==c.size());
  for(uint i=0; i!=c.size();++i) {
    x[i]-=c[i];
  }
  return x;
}

  
template<class DIFF>
DifferentialVector<DIFF> 
operator+(const DifferentialVector<DIFF>& x, const Vector<typename DIFF::ScalarType>& c)
{  
  DifferentialVector<DIFF> r(x);
  return r+=c;
}


template<class DIFF>
DifferentialVector<DIFF> 
operator-(const DifferentialVector<DIFF>& x, const Vector<typename DIFF::ScalarType>& c)
{  
  DifferentialVector<DIFF> r(x);
  return r-=c;
}


template<class DIFF>
DifferentialVector<DIFF> 
operator*(const Matrix<typename DIFF::ScalarType>& A, const DifferentialVector<DIFF>& x)
{  
  assert(A.column_size()==x.result_size());
  DifferentialVector<DIFF> r(A.row_size(),x.argument_size(),x.degree());
  for(uint i=0; i!=A.row_size();++i) {
    for(uint j=0; j!=A.column_size();++j) {
      r[i]+=A[i][j]*x[j];
    }
  }
  return r;
}


  
template<class DIFF>
DifferentialVector<DIFF> 
operator+(const DifferentialVector<DIFF>& x, const DifferentialVector<DIFF>& y)
{  
  DifferentialVector<DIFF> r(x);
  return r+=y;
}


template<class DIFF>
DifferentialVector<DIFF> 
operator-(const DifferentialVector<DIFF>& x, const DifferentialVector<DIFF>& y)
{  
  DifferentialVector<DIFF> r(x);
  return r-=y;
}


template<class DIFF>
DifferentialVector<DIFF>
join(const DifferentialVector<DIFF>& f, const DifferentialVector<DIFF>& g)
{
  ARIADNE_ASSERT(f.argument_size()==g.argument_size());
  DifferentialVector<DIFF> h(f.result_size()+g.result_size(),f.argument_size(),std::max(f.degree(),g.degree()));
  for(uint i=0; i!=f.result_size(); ++i) {
    h[i]=f[i];
  }
  for(uint i=0; i!=g.result_size(); ++i) {
    h[i+f.result_size()]=g[i];
  }
  return h;
}


template<class DIFF>
DifferentialVector<DIFF>
join(const DifferentialVector<DIFF>& f, const DIFF& g)
{
  ARIADNE_ASSERT(f.argument_size()==g.argument_size());
  DifferentialVector<DIFF> h(f.result_size()+1u,f.argument_size(),std::max(f.degree(),g.degree()));
  for(uint i=0; i!=f.result_size(); ++i) {
    h[i]=f[i];
  }
  h[f.result_size()]=g;
  return h;
}


template<class DIFF, class Y>
Y
evaluate(const DIFF& x, 
         const Vector<Y>& y)
{  
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  using namespace std;
  assert(x.argument_size()==y.size());
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;

  uint d=x.degree();
  uint s=y.size();

  Y zero=y[0]; zero*=0;
  Y one=zero; one+=1;

  Y r=zero;

  //std::cerr << "zero="<<zero<<std::endl;
  //std::cerr << "one="<<one<<std::endl;
  // Use inefficient brute-force approach with lots of storage...
  array< array< Y > > val(s, array< Y >(d+1,zero));
  for(uint j=0; j!=s; ++j) {
    val[j][0]=one;
    for(uint k=1; k<=d; ++k) {
      val[j][k]=val[j][k-1]*y[j];
    }
  }
  //std::cerr << "val="<<val<<std::endl;
  for(MultiIndex j(s); j.degree()<=d; ++j) {
    Y t=one;
    for(uint k=0; k!=s; ++k) {
      t=t*val[k][j[k]];
    }
    r+=x[j]*t;
  }
  
  return r;
}


template<class DIFF, class Y>
Vector<Y>
evaluate(const DifferentialVector<DIFF>& x, 
         const Vector<Y>& y)
{  
  //std::cerr<<__PRETTY_FUNCTION__<<std::endl;
  using namespace std;
  assert(x.argument_size()==y.size());
  typedef typename DIFF::ScalarType X;
  //std::cerr << "y=" << y << std::endl;
  //std::cerr << "x=" << x << std::endl;

  uint d=x.degree();
  uint rs=x.result_size();
  uint as=x.argument_size();

  Y zero=y[0]; zero*=0;
  Y one=zero; one+=1;

  Vector<Y> r(rs,zero);

  // Use inefficient brute-force approach with lots of storage...
  array< array< Y > > val(as, array< Y >(d+1,zero));
  for(uint j=0; j!=as; ++j) {
    val[j][0]=one;
    for(uint k=1; k<=d; ++k) {
      val[j][k]=val[j][k-1]*y[j];
    }
  }
  for(MultiIndex j(as); j.degree()<=d; ++j) {
    Y t=one;
    for(uint k=0; k!=as; ++k) {
      t=t*val[k][j[k]];
    }
    for(uint i=0; i!=rs; ++i) {
      const X& xij=x[i][j];
      Y& ri=r[i];
      Y txij=xij*t;
      ri+=txij;
    }
  }
  
  return r;
}



template<class DIFF>
DifferentialVector<DIFF> 
evaluate(const DifferentialVector<DIFF>& x, 
         const DifferentialVector<DIFF>& y)
{  
  assert(x.argument_size()==y.result_size());
  DifferentialVector<DIFF> r;
  
  static_cast<Vector<DIFF>&>(r) = 
    evaluate(x,static_cast<const Vector<DIFF>&>(y));
  for(uint i=0; i!=r.result_size(); ++i) { r[i].cleanup(); }
  return r;
}


template<class DIFF>
DIFF
embed(const DIFF& x, 
      uint size, uint start)
{  
  assert(start+x.argument_size()<=size);
  DIFF r(size,x.degree());
  MultiIndex jr(size);
  for(typename DIFF::const_iterator iter=x.begin();
      iter!=x.end(); ++iter)
  {
    const MultiIndex& jx=iter->first;
    for(uint k=0; k!=x.argument_size(); ++k) {
      jr.set(start+k,jx[k]);
    }
    r[jr]=iter->second;
  }
  return r;
}

template<class DIFF>
DifferentialVector<DIFF> 
embed(const DifferentialVector<DIFF>& x, 
      uint size, uint start)
{  
  assert(start+x.argument_size()<=size);
  DifferentialVector<DIFF> r(x.result_size(),size,x.degree());
  for(uint i=0; i!=x.result_size(); ++i) { r[i]=embed(x[i],size,start); }
  return r;
}


/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class DIFF>
DIFF 
compose(const DIFF& x, 
        const DifferentialVector<DIFF>& y)
{  
  typedef typename DIFF::ScalarType X;
  Vector<X> yv=y.value();
  DifferentialVector<DIFF>& ync=const_cast<DifferentialVector<DIFF>&>(y); 
  for(uint i=0; i!=ync.result_size(); ++i) { ync[i].value()=0; }
  DIFF r=evaluate(x,ync);
  ync+=yv;
  return r;
}


/*! \brief Compose the series of \a x and the series of \a y, assuming that \a x is centred at the value of \a y. The value of \a y is therefore unused by this method. */
template<class DIFF>
DifferentialVector<DIFF> 
compose(const DifferentialVector<DIFF>& x, 
        const DifferentialVector<DIFF>& y)
{  
  typedef typename DIFF::ScalarType X;
  Vector<X> yv=y.value();
  DifferentialVector<DIFF>& ync=const_cast<DifferentialVector<DIFF>&>(y); 
  for(uint i=0; i!=ync.result_size(); ++i) { ync[i].value()=0; }
  DifferentialVector<DIFF> r=evaluate(x,ync);
  ync+=yv;
  return r;
}

template<class DIFF> 
DifferentialVector<DIFF> 
implicit(const DifferentialVector<DIFF>& x)
{
  typedef typename DIFF::ScalarType X;
  assert(x.result_size()<=x.argument_size());
  //std::cerr << "x=" << x << std::endl;
  
  uint rs=x.result_size();
  uint xas=x.argument_size();
  uint zas=x.argument_size()-x.result_size();
  uint d=x.degree();

  Matrix<X> A1(rs,zas);
  for(uint i=0; i!=rs; ++i) {
    for(uint j=0; j!=zas; ++j) {
      A1(i,j)=x[i].gradient(j);
    }
  }
  
  Matrix<X> A2(rs,rs);
  for(uint i=0; i!=rs; ++i) {
    for(uint j=0; j!=rs; ++j) {
      A2(i,j)=x[i].gradient(zas+j);
    }
  }
  
  Matrix<X> J(xas,rs);
  //J(range(zas,zas+rs),range(0,rs))=inverse(A2);
  project(J,range(zas,zas+rs),range(0,rs)) = inverse(A2);

  DifferentialVector<DIFF> y(xas,zas,d);
  for(uint i=0; i!=zas; ++i) {
    y[i]=DIFF::variable(zas,d,1.0,i);
  }
  for(uint i=0; i!=rs; ++i) {
    // y[as+i]=TaylorVariable<X>::constant(as,d,0.0);
  }

  for(uint i=0; i!=d; ++i) {
    DifferentialVector<DIFF> z=compose(x,y);
    y-=J*z;
  }

  DifferentialVector<DIFF> r(rs,zas,d);
  for(uint i=0; i!=rs; ++i) {
    r[i]=y[zas+i];
  }
  return r;
}



template<class DIFF> 
DifferentialVector<DIFF> 
inverse(const DifferentialVector<DIFF>& x)
{
  typedef typename DIFF::ScalarType X;
  return inverse(x,Vector<X>(x.argument_size()));
}

template<class DIFF> 
DifferentialVector<DIFF> 
inverse(const DifferentialVector<DIFF>& x, const Vector<typename DIFF::ScalarType>& c)
{
  using namespace std;
  typedef typename DIFF::ScalarType X;
  assert(x.result_size()==x.argument_size());
  assert(x.result_size()==c.size());
  //std::cerr << "x=" << x << std::endl;
  uint n=x.result_size();
  uint d=x.degree();
  Vector<X> z(n,0);
  Matrix<X> J=inverse(x.jacobian());

  DifferentialVector<DIFF> y(n,n,d);
  DifferentialVector<DIFF> id(n,n,d);
  for(uint i=0; i!=n; ++i) { id[i][i]=1.0; }
  
  for(uint i=0; i!=n; ++i) { 
    y[i].value()=c[i]; 
    for(uint j=0; j!=n; ++j) { 
      y[i][j]=J[i][j];
    }
  }

  for(uint i=2; i<=d; ++i) {
    DifferentialVector<DIFF> z=compose(x,y);
    z-=id;
    y-=J*z;
  }
  return y;
}



template<class DIFF> 
DifferentialVector<DIFF>
flow1(const DifferentialVector<DIFF>& f, const Vector<typename DIFF::ScalarType>& x, uint to, uint so)
{
  typedef typename DIFF::ScalarType X;
  // f is an untimed vector field
  assert(f.result_size()==f.argument_size());
  assert(f.result_size()==x.size());
  uint n=x.size();
  uint d=f.degree();

  Matrix<X> Q(n+1,n); for(uint i=0; i!=n; ++i) { Q[i][i]=1; }
  
  Vector<X> tx(n+1);
  for(uint i=0; i!=n; ++i) { tx[i]=x[i]; } 
  tx[n]=0.0;

  /*
    DifferentialVector<DIFF> tf(n,n+1,d);
  MultiIndex ta(n+1);
  for(uint i=0; i!=n; ++i) {
    for(typename DIFF::const_iterator iter=f[i].begin(); iter!=f[i].end(); ++iter) {
      const MultiIndex& a=iter->first;
      for(uint k=0; k!=n; ++k) { ta.set(k,a[k]); }
      tf[i][ta]=iter->second;
    }
  }
  */
  //std::cerr << "f=" << f << std::endl;
  //std::cerr << "tf=" << tf << std::endl << std::endl;

  DifferentialVector<DIFF> y(n,n+1,d);
  for(uint i=0; i!=n; ++i) { y[i].gradient(i)=1; }
  //std::cerr << "y=" << y << std::endl;
  DifferentialVector<DIFF> yp(n,n+1,d);
  for(uint j=0; j<d; ++j) {
    yp=compose(f,y);
    //std::cerr << "yp=" << yp << std::endl;
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],n);
      y[i].value()=0;
      y[i].gradient(i)=1;
    }
    //std::cerr << "y=" << y << std::endl << std::endl;
  } 
  for(uint i=0; i!=n; ++i) { y[i].value()=x[i]; }
  return y;
}



template<class DIFF> 
DifferentialVector<DIFF>
flow(const DifferentialVector<DIFF>& f, const Vector<typename DIFF::ScalarType>& x, uint to, uint so)
{
  return flow1(f,x,to,so);
}



//! Compute the flow map to the crossing set g under the vector field \a vf
template<class DIFF> 
DifferentialVector<DIFF>  
hitting(const DifferentialVector<DIFF>& vf, const DIFF& g)
{
}


//! Translate the polynomial given by \a x to one with centre \a v.
template<class DIFF> 
DifferentialVector<DIFF>  
translate(const DifferentialVector<DIFF>& x, const Vector<typename DIFF::ScalarType>& v)
{
  uint as=v.size();
  uint d=x.degree();
  DifferentialVector<DIFF> t=DifferentialVector<DIFF>::variable(as,as,d,v);
  return evaluate(x,t);
}


//! Scale the polynomial given by \a x by the values in the array \a s. 
template<class DIFF> 
DifferentialVector<DIFF>  
scale(const DifferentialVector<DIFF>& x, const Vector<typename DIFF::ScalarType>& s)
{
  uint as=s.size();
  uint d=x.degree();
  DifferentialVector<DIFF> t(as,as,d);
  for(uint i=0; i!=as; ++i) { t[i][i]=s[i]; } 
  return evaluate(x,t);
}










} //namespace Ariadne

#endif /* ARIADNE_SPARSE_DIFFERENTIAL_H */
