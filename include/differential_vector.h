/***************************************************************************
 *            differential_vector.h
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
 
/*! \file differential_vector.h
 *  \brief Operations on elements of a differential algebra.
 */
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




template<class DIFF> DifferentialVector<DIFF> differential_vector(uint rs, uint as, uint d, const typename DIFF::RealType* ptr);

template<class DIFF> DifferentialVector<DIFF> vector_variable(uint n, uint d, const Vector<typename DIFF::RealType>& c);

template<class DIFF> DifferentialVector<DIFF> project(const DifferentialVector<DIFF>& y, Slice rng);
template<class DIFF> DifferentialVector<DIFF> restrict(const DifferentialVector<DIFF>& y, const array<uint>& p);
template<class DIFF> DifferentialVector<DIFF> expand(const DifferentialVector<DIFF>& y, uint as, const array<uint>& p);
template<class DIFF> DifferentialVector<DIFF> join(const DifferentialVector<DIFF>& x1, const DifferentialVector<DIFF>& x2);
template<class DIFF> DifferentialVector<DIFF> join(const DifferentialVector<DIFF>& x1, const DIFF& x2);

template<class DIFF> DifferentialVector<DIFF> operator-(const DifferentialVector<DIFF>& x, const DifferentialVector<DIFF>& y);
template<class DIFF> DifferentialVector<DIFF> operator+(const DifferentialVector<DIFF>& v, const Vector<typename DIFF::RealType>& c);
template<class DIFF> DifferentialVector<DIFF> operator-(const DifferentialVector<DIFF>& v, const Vector<typename DIFF::RealType>& c);

template<class DIFF, class Y> Vector<Y> evaluate(const DifferentialVector<DIFF>& x, const Vector<Y>& y);
template<class DIFF> DifferentialVector<DIFF> evaluate(const DifferentialVector<DIFF>& x, const DifferentialVector<DIFF>& y);

template<class DIFF> DifferentialVector<DIFF> translate(const DifferentialVector<DIFF>& y, const Vector<typename DIFF::RealType>& c);
template<class DIFF> DifferentialVector<DIFF> compose(const DifferentialVector<DIFF>& y, const DifferentialVector<DIFF>& z);
template<class DIFF> DifferentialVector<DIFF> inverse(const DifferentialVector<DIFF>& y);
template<class DIFF> DifferentialVector<DIFF> implicit(const DifferentialVector<DIFF>& y);
template<class DIFF> DifferentialVector<DIFF> derivative(const DifferentialVector<DIFF>& x, uint i);
template<class DIFF> DifferentialVector<DIFF> antiderivative(const DifferentialVector<DIFF>& x, uint j);
//template<class DIFF> DifferentialVector<DIFF> flow(const DifferentialVector<DIFF>& vf, const Vector<typename DIFF::RealType>& x, uint ox);




/*! \brief A class representing the derivatives of a vector quantity depending on multiple arguments. */
template<class DIFF>
class DifferentialVector
    : public Vector< DIFF >
{
    //BOOST_CONCEPT_ASSERT((DifferentialVectorConcept<SparseDifferentialVector<X> >));
    typedef typename DIFF::ScalarType X;
  public:
    // The type used for accessing elements
    typedef uint IndexType;
    // The value stored in the vector.
    typedef DIFF ValueType;
    // The type used for scalars.
    typedef typename DIFF::ScalarType ScalarType;

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

    Vector<X> get_value() const { 
        Vector<X> r(this->result_size()); for(uint i=0; i!=r.size(); ++i) { r[i]=(*this)[i].value(); } return r; }
    Matrix<X> get_jacobian() const { Matrix<X> r(this->result_size(),this->argument_size()); 
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

    static DifferentialVector<DIFF> affine(uint rs, uint as, uint d, const Vector<X>& b, const Matrix<X>& A) {
        ARIADNE_ASSERT(b.size()==rs);
        ARIADNE_ASSERT(A.row_size()==rs);
        ARIADNE_ASSERT(A.column_size()==as);
        DifferentialVector<DIFF> result(rs,as,d);
        for(uint i=0; i!=rs; ++i) { 
            result[i]=b[i]; 
            for(uint j=0; j!=as; ++j) {
                result[i][j]=A[i][j]; 
            }
        }
        return result;
    }

};


template<class DIFF> DifferentialVector<DIFF> 
differential_vector(uint rs, uint as, uint d, const typename DIFF::RealType* ptr)
{
    return DifferentialVector<DIFF>(rs,as,d,ptr);
}


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
    return static_cast<Vector<DIFF>&>(r)-=y;
}


template<class DIFF>
DifferentialVector<DIFF> 
restrict(const DifferentialVector<DIFF>& x, const array<uint>& p)
{
    uint d=x.degree();
    uint rs=x.result_size();
    DifferentialVector<DIFF> r(x.result_size(),p.size(),x.degree());
    MultiIndex rj(r.argument_size());
    MultiIndex xj(x.argument_size());
    for( ; rj.degree()<=d; ++rj) {
        for(uint k=0; k!=p.size(); ++k) {
            xj.set(p[k],rj[k]);
        }
        for(uint i=0; i!=rs; ++i) {
            r[i][rj]=x[i][xj];
        }
    }
    return r;
}


template<class DIFF>
DifferentialVector<DIFF> 
expand(const DifferentialVector<DIFF>& x, uint as, const array<uint>& p)
{
    uint d=x.degree();
    uint rs=x.result_size();
    DifferentialVector<DIFF> r(rs,as,d);
    MultiIndex rj(r.argument_size());
    MultiIndex xj(x.argument_size());
    for( ; xj.degree()<=d; ++xj) {
        for(uint k=0; k!=p.size(); ++k) {
            rj.set(p[k],xj[k]);
        }
        for(uint i=0; i!=rs; ++i) {
            r[i][rj]=x[i][xj];
        }
    }
    return r;
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
    //std::cerr<<" y="<<y<<std::endl;
    //std::cerr<<" x="<<x<<std::endl;

    using namespace std;
    assert(x.argument_size()==y.size());

    uint d=x.degree();
    uint s=y.size();

    Y zero=y[0]; zero*=0;  
    Y one=zero; one+=1;

    Y r=zero;

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
    //std::cerr<<" y="<<y<<std::endl;
    //std::cerr<<" x="<<x<<std::endl;
    using namespace std;
    assert(x.argument_size()==y.size());
    typedef typename DIFF::ScalarType X;

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
DifferentialVector<DIFF> 
project(const DifferentialVector<DIFF>& x, Slice slc)
{
    DifferentialVector<DIFF> r(slc.size(),x.argument_size(),x.degree());
    for(uint i=0; i!=slc.size(); ++i) {
        r[i]=x[i+slc.start()];
    }
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
    Vector<X> yv=y.get_value();
    DifferentialVector<DIFF>& ync=const_cast<DifferentialVector<DIFF>&>(y); 
    for(uint i=0; i!=ync.result_size(); ++i) { ync[i].set_value(0); }
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
    //std::cerr<<"compose(DV x, DV y)\n x="<<x<<"\n y="<<y<<std::endl;
    typedef typename DIFF::ScalarType X;
    Vector<X> yv=y.get_value();
    DifferentialVector<DIFF>& ync=const_cast<DifferentialVector<DIFF>&>(y); 
    for(uint i=0; i!=ync.result_size(); ++i) { ync[i].set_value(0); }
    DifferentialVector<DIFF> r=evaluate(x,ync);
    //std::cerr<<"r="<<r<<"\n"<<std::endl;
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
    // std::cerr<<"A2="<<A2<<std::endl;
    // std::cerr<<"J="<<J<<std::endl;

    DifferentialVector<DIFF> y(xas,zas,d);
    for(uint i=0; i!=zas; ++i) {
        y[i]=DIFF::variable(zas,d,1.0,i);
    }
    for(uint i=0; i!=rs; ++i) {
        y[zas+i]=DIFF::constant(zas,d,0.0);
    }

    // std::cerr<<"\nx="<<x<<std::endl;
    // std::cerr<<"y="<<y<<std::endl;
    for(uint i=0; i!=d; ++i) {
        DifferentialVector<DIFF> z=compose(x,y);
        // std::cerr<<"z="<<z<<std::endl;
        // Need cast to avoid ambiguity
        static_cast<Vector<DIFF>&>(y)-=J*z;
        // std::cerr<<"y="<<y<<std::endl;
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
    //std::cerr << "\ninverse(x)\nx=" << x << "\n" << std::endl;
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
    //std::cerr << "\ninverse(x,c)\nx=" << x << "\nc=" << c << "\n" << std::endl;
    uint n=x.result_size();
    uint d=x.degree();
    Vector<X> z(n,0);
    Matrix<X> J=inverse(x.get_jacobian());

    DifferentialVector<DIFF> y(n,n,d);
    DifferentialVector<DIFF> id(n,n,d);
    for(uint i=0; i!=n; ++i) { id[i][i]=1.0; }
  
    for(uint i=0; i!=n; ++i) { 
        y[i].set_value(c[i]); 
        for(uint j=0; j!=n; ++j) { 
            y[i].set_gradient(j,J[i][j]);
        }
    }
  
    for(uint i=2; i<=d; ++i) {
        DifferentialVector<DIFF> z=compose(x,y);
        static_cast<Vector<DIFF>&>(z)-=id;
        static_cast<Vector<DIFF>&>(y)-=J*z;
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

    // Perform antiderivatives starting from zero
    // The code below will compute all derivates correctly up to order d
    // However, in principle, derivatives with respect to time can be 
    // computed up to order d+1. This has not been done in the code below
    // since in the dense representation, we have a fixed order
    // for all variables. Of course, in the sparse representation, 
    // we could modify the code by setting the gradients at the beginning
  
  
    //std::cerr << "\nf=" << f << "\n" << std::endl;
    DifferentialVector<DIFF> y(n,n+1,0);
    //for(uint i=0; i!=n; ++i) { y[i].gradient(i)=1; }
    //std::cerr << "y[0]=" << y << std::endl;
    DifferentialVector<DIFF> yp(n,n+1,d);
    for(uint j=0; j<d; ++j) {
        yp=compose(f,y);
        //std::cerr << "yp["<<j<<"]=" << yp << std::endl;
        for(uint i=0; i!=n; ++i) {  
            y[i]=antiderivative(yp[i],n);
            y[i].set_value(0);
            y[i].set_gradient(i,1);
        }
        //std::cerr << "y["<<j+1<<"]=" << y << std::endl << std::endl;
    } 
    for(uint i=0; i!=n; ++i) { y[i].set_value(x[i]); }
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
