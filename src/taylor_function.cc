/***************************************************************************
 *            taylor_function.cc
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

#include <iostream>
#include <iomanip>

#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "polynomial.h"
#include "differential.h"
#include "function_interface.h"
#include "taylor_variable.h"
#include "taylor_function.h"

namespace Ariadne {

typedef unsigned int uint;

typedef Vector<Float> Point;
typedef Vector<Interval> Box;

TaylorFunction::TaylorFunction()
    : _domain(),
      _variables()
{
}


TaylorFunction::TaylorFunction(uint rs, uint as)
    : _domain(as,I(-1,1)),
      _variables(rs,as)
{
}


TaylorFunction::TaylorFunction(uint rs, uint as, ushort o, ushort s)
    : _domain(as,I(-1,1)),
      _variables(rs,as)
{
}


TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector< Expansion<Float> >& f)
    : _domain(d),
      _variables(f.size())
{
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
    }
    for(uint i=0; i!=f.size(); ++i) {
        _variables[i]=TaylorVariable(d,f[i]);
    }
}

TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector< Expansion<Float> >& f,
                               const Vector<Float>& e)
    : _domain(d),
      _variables(e.size())
{
    for(uint i=0; i!=f.size(); ++i) {
        ARIADNE_ASSERT(d.size()==f[i].argument_size());
    }
    ARIADNE_ASSERT(f.size()==e.size());
    for(uint i=0; i!=e.size(); ++i) {
        _variables[i]=TaylorVariable(d,f[i],e[i]);
    }
}


TaylorFunction::TaylorFunction(const Vector<TaylorVariable>& v)
    : _domain(v[0].domain()),
      _variables(v)
{
    for(uint i=1; i!=v.size(); ++i) {
        ARIADNE_ASSERT(v[0].domain()==v[i].domain());
    }
}

TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const FunctionInterface& f)
    : _domain(d),
      _variables(f.result_size())
{
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<TaylorVariable> x=TaylorVariable::scaling(d);
    this->_variables=f.evaluate(x);
}


TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                               const Vector< Polynomial<Float> >& p)
    : _domain(d),
      _variables(p.size())
{
    ARIADNE_ASSERT(d.size()==p[0].argument_size());

    Vector<TaylorVariable> x=TaylorVariable::scaling(d);
    this->_variables=Ariadne::evaluate(p,x);
}





TaylorFunction::TaylorFunction(const Vector<Interval>& d,
                         ushort o, ushort s,
                         const FunctionInterface& f)
    : _domain(d),
      _variables(f.result_size())
{
    ARIADNE_NOT_IMPLEMENTED;
}




/*
TaylorFunction
TaylorFunction::zero(uint rs, uint as)
{
    return TaylorFunction(Vector<Interval>(as,Interval(-inf(),+inf()),TaylorVariable::constants(rs,as,Vector<Float>(rs,0.0))));
}
*/




TaylorFunction
TaylorFunction::constant(const Vector<Interval>& d, const Vector<Float>& c)
{
    return TaylorFunction(TaylorVariable::constants(d,c));
}

TaylorFunction
TaylorFunction::constant(const Vector<Interval>& d, const Vector<Interval>& r)
{
    Vector<TaylorVariable> e(r.size(),d);
    for(uint i=0; i!=r.size(); ++i) {
        e[i]=r[i];
    }
    return TaylorFunction(e);
}

TaylorFunction
TaylorFunction::identity(const Vector<Interval>& d)
{
    return TaylorFunction(TaylorVariable::scaling(d));
}



bool
TaylorFunction::operator==(const TaylorFunction& tm) const
{
        return this->_domain==tm._domain
            && this->_variables==tm._variables;
}



bool
TaylorFunction::operator!=(const TaylorFunction& p2) const
{
    return !(*this==p2);
}



const Vector<Interval>&
TaylorFunction::domain() const
{
    return this->_domain;
}


const Vector<Interval>
TaylorFunction::range() const
{
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_variables[i].range();
    }
    return result;
}



const Vector<TaylorVariable>&
TaylorFunction::variables() const
{
    return this->_variables;
}




uint
TaylorFunction::argument_size() const
{
    return this->_variables[0].argument_size();
}


uint
TaylorFunction::result_size() const
{
    return this->_variables.size();
}






















TaylorFunction
TaylorFunction::truncate(ushort degree) const
{
    ARIADNE_NOT_IMPLEMENTED;
}



Vector<Interval>
TaylorFunction::evaluate(const Vector<Float>& x) const
{
    return this->evaluate(Vector<Interval>(x));
}


Vector<Interval>
TaylorFunction::evaluate(const Vector<Interval>& x) const
{
    if(this->argument_size()!=x.size()) {
        ARIADNE_THROW(std::runtime_error,"TaylorFunction::evaluate(Vector)","Incompatible argument size");
    }

    //TODO: Make this MUCH more efficient!

    // Scale x to domain
    Vector<Interval> scaled_x(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        const Interval& d=this->_domain[i];
        Interval dm=add_ivl(d.l/2,d.u/2);
        Interval dr=sub_ivl(d.u/2,d.l/2);
        scaled_x[i]=(x[i]-dm)/dr;
    }

    // Compute result
    Vector<Interval> result(this->result_size());
    for(uint i=0; i!=result.size(); ++i) {
        result[i]=this->_variables[i].evaluate(scaled_x);
    }
    return result;
}


Matrix<Float>
TaylorFunction::jacobian(const Vector<Float>& x) const
{
    Vector< Differential<Float> > s(this->argument_size(),this->argument_size(),1u);
    for(uint j=0; j!=this->argument_size(); ++j) {
        Interval dj=this->domain()[j];
        s[j].set_value((x[j]-dj.midpoint())/dj.radius());
        s[j].set_gradient(j,1/dj.radius());
    }
    Vector< Expansion<Float> > p(this->result_size(),this->argument_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        p[i]=this->variables()[i].expansion();
    }
    Vector< Differential<Float> > d=Ariadne::evaluate(p,s);
    //std::cerr<<"  x="<<x<<"\n  p="<<p<<"\n"<<"  s="<<s<<"\n  p.s="<<d<<"\n  J="<<d.jacobian()<<"\n"<<std::endl;
    return d.jacobian();
}


Matrix<Interval>
TaylorFunction::jacobian(const Vector<Interval>& x) const
{
    Vector< Differential<Interval> > s(this->argument_size(),this->argument_size(),1u);
    for(uint j=0; j!=this->argument_size(); ++j) {
        Interval dj=this->domain()[j];
        s[j].set_value((x[j]-dj.midpoint())/dj.radius());
        s[j].set_gradient(j,1/dj.radius());
    }
    Vector< Expansion<Float> > p(this->result_size(),this->argument_size());
    for(uint i=0; i!=this->result_size(); ++i) {
        p[i]=this->variables()[i].expansion();
    }
    Vector< Differential<Interval> > d=Ariadne::evaluate(p,s);
    //std::cerr<<"  x="<<x<<"\n  p="<<p<<"\n"<<"  s="<<s<<"\n  p.s="<<d<<"\n  J="<<d.jacobian()<<"\n"<<std::endl;
    return d.jacobian();
}


TaylorFunction
join(const TaylorFunction& f, const TaylorFunction& g)
{
    ARIADNE_ASSERT(f.domain()==g.domain());
    return TaylorFunction(join(f.variables(),g.variables()));
}


TaylorFunction
restrict(const TaylorFunction& tm, const Box& nd)
{
    const Vector<Interval>& od=tm.domain();
    ARIADNE_ASSERT(subset(nd,od));
    if(nd==od) { return tm; }
    Vector<TaylorVariable> sm=Vector<TaylorVariable>(nd.size(),nd.size());
    for(uint i=0; i!=nd.size(); ++i) {
        set_rounding_mode(upward);
        const Float& l=nd[i].lower();
        const Float& u=nd[i].upper();
        volatile Float ce,re,mcl,cu,mrl,ru;
        mcl=-l; mcl-=u;
        cu=l+u;
        ce=(mcl+cu)/4;
        mrl=l; mrl-=u;
        ru=u-l;
        re=mrl+ru;
        sm[i].set_error(ce+re);
        set_rounding_mode(to_nearest);
        sm[i].set_value((l+u)/2);
        sm[i].set_gradient(i,(u-l)/2);
    }
    return TaylorFunction(compose(tm._variables,tm._domain,sm));
}



TaylorFunction&
operator+=(TaylorFunction& f, const Vector<Interval>& e)
{
    ARIADNE_ASSERT(f.result_size()==e.size());
    f._variables+=e;
    return f;
}


TaylorFunction
operator+(const TaylorFunction& p1, const TaylorFunction& p2)
{
    ARIADNE_ASSERT(!intersection(p1.domain(),p2.domain()).empty());
    if(p1.domain()==p2.domain()) {
        return TaylorFunction(Vector<TaylorVariable>(p1._variables+p2._variables));
    } else {
        Box new_domain=intersection(p1.domain(),p2.domain());
        return operator+(restrict(p1,new_domain),restrict(p2,new_domain));
    }
}


TaylorFunction
operator-(const TaylorFunction& p1, const TaylorFunction& p2)
{
    ARIADNE_ASSERT(!intersection(p1.domain(),p2.domain()).empty());
    if(p1.domain()==p2.domain()) {
        return TaylorFunction(Vector<TaylorVariable>(p1._variables-p2._variables));
    } else {
        ARIADNE_NOT_IMPLEMENTED;
        Box new_domain=intersection(p1.domain(),p2.domain());
        return operator-(restrict(p1,new_domain),restrict(p2,new_domain));
    }
}



TaylorFunction
operator+(const TaylorFunction& f, const Vector<Float>& c)
{
    return TaylorFunction(Vector<TaylorVariable>(f._variables+c));
}

TaylorFunction
operator+(const TaylorFunction& f, const Vector<Interval>& c)
{
    return TaylorFunction(Vector<TaylorVariable>(f._variables+c));
}

TaylorFunction
operator-(const TaylorFunction& f, const Vector<Float>& c)
{
    return TaylorFunction(Vector<TaylorVariable>(f._variables-c));
}

TaylorFunction
operator-(const TaylorFunction& f, const Vector<Interval>& c)
{
    return TaylorFunction(Vector<TaylorVariable>(f._variables-c));
}

TaylorFunction
operator*(const Matrix<Interval>& A, const TaylorFunction& f)
{
    return TaylorFunction(Vector<TaylorVariable>(prod(A,f._variables)));
}


/*

TaylorFunction
mul(const TaylorFunction& tm, const Float& x)
{
    return TaylorFunction(tm.domain(),tm.centre(),tm.smoothness(),tm.centre_derivatives()*x);
}

*/


/*

void
mul(TaylorFunction& p0, const TaylorFunction& p1, const TaylorFunction& p2)
{
    if(p1.result_size()!=1u) {
        ARIADNE_THROW(IncompatibleSizes,"mul(TaylorFunction,TaylorFunction)","p1.result_size()="<<p1.result_size());
    }
    if(p2.result_size()!=1u) {
        ARIADNE_THROW(IncompatibleSizes,"mul(TaylorFunction,TaylorFunction)","p2.result_size()="<<p2.result_size());
    }
    if(p1.argument_size()!=p2.argument_size()) {
        ARIADNE_THROW(IncompatibleSizes,"add(TaylorFunction p1,TaylorFunction p2)","p1.result_size()="<<p1.result_size()<<", p2.result_size()="<<p2.result_size());
    }

    uint rs=1u;
    uint as=p1.argument_size();
    uint d1=p1.order();
    uint d2=p2.order();
    uint d=d1+d2;
    uint s1=p1.smoothness();
    uint s2=p2.smoothness();
    uint s=std::max(s1,s2);

    p0.resize(rs,as,d,s);
    MultiIndex j0(as);
    for(MultiIndex j1(as); j1.degree()<=d1; ++j1) {
        const Float& x1=p1.get(0u,j1);
        for(MultiIndex j2(as); j2.degree()<=d2; ++j2) {
            const R& x2=p2.get(0u,j2);
            j0=j1+j2;
            p0.at(0u,j0)+=x1*x2;
        }
    }
}


*/





/*
TaylorFunction&
operator*=(TaylorFunction& p0, const Interval& x1)
{
    for(uint i=0; i!=p0._data.size(); ++i) {
        p0._data[i]*=x1;
    }
    return p0;
}
*/

TaylorFunction
combine(const TaylorFunction& f1, const TaylorFunction& f2)
{
    uint as1=f1.argument_size();
    uint as2=f2.argument_size();
    Vector<Interval> d=join(f1.domain(),f2.domain());
    Vector<TaylorVariable> ev1=embed(f1.variables(),as1+as2,0u);
    Vector<TaylorVariable> ev2=embed(f2.variables(),as1+as2,as1);
    Vector<TaylorVariable> ev=join(ev1,ev2);
    return TaylorFunction(ev);
}


TaylorFunction
compose(const TaylorFunction& g, const TaylorFunction& f)
{
    if(!subset(f.range(),g.domain())) {
        std::cerr<<"f.range()="<<f.range()<<" is not a subset of g.domain()="<<g.domain()<<std::endl;
        ARIADNE_ASSERT(subset(f.range(),g.domain()));
    }
    return TaylorFunction(Ariadne::compose(g.variables(),g.domain(),f.variables()));
}



TaylorFunction
antiderivative(const TaylorFunction& tm, uint k)
{
    return TaylorFunction(antiderivative(tm.variables(),tm.domain()[k],k));
}



/*
void
compose(TaylorFunction& p0, const TaylorFunction& p1, const TaylorFunction& p2)
{
    // TODO: Improve this algorithm as it's critical!!
    if(p1.argument_size()!=p2.result_size()) {
        ARIADNE_THROW(IncompatibleSizes,"compose(TaylorFunction p1,TaylorFunction p2)","p1.argument_size()="<<p1.argument_size()<<", p2.result_size()="<<p2.result_size());
    }

    if(p1.order()==0) {
        p0=static_cast< TaylorFunction >(p1);
        return;
    }

    p0.resize(p1.result_size(),p2.argument_size(),p1.order()*p2.order(),std::max(p1.smoothness(),p2.smoothness()));

    TaylorFunction* all_powers=new TaylorFunction[p2.result_size()*(p1.order()+1)];
    TaylorFunction* powers[p2.result_size()];
    for(uint i=0; i!=p2.result_size(); ++i) {
        powers[i]=all_powers+i*(p1.order()+1);
    }

    for(uint i=0; i!=p2.result_size(); ++i) {
        powers[i][0]=TaylorFunction::one(p2.argument_size());
        powers[i][1]=p2.component(i);
        if(p1.order()>=2) {
            powers[i][2]=pow(powers[i][1],2);
        }
        for(uint j=3; j<=p1.order(); ++j) {
            powers[i][j]=powers[i][2]*powers[i][j-2];
        }
    }

    TaylorFunction* results=new TaylorFunction[p1.result_size()];
    for(uint i=0; i!=p1.result_size(); ++i) {
        results[i]=TaylorFunction::zero(1u,p2.argument_size());
    }

    for(uint i=0; i!=p1.result_size(); ++i) {
        for(MultiIndex j(p1.argument_size()); j.degree()<=p1.order(); ++j) {
            TaylorFunction t=TaylorFunction::constant(p2.argument_size(),p1.get(i,j));
            for(uint k=0; k!=p1.argument_size(); ++k) {
                t=t*powers[k][j[k]];
            }
            results[i]=results[i]+t;
        }
    }

    for(uint i=0; i!=p0.result_size(); ++i) {
        for(uint j=0; j!=p0.data().size()/p0.result_size(); ++j) {
            p0._data[i+j*p0.result_size()]=results[i].data()[j];
        }
    }

    delete[] results;
    delete[] all_powers;
}

*/

/*

void
derivative(TaylorFunction& p0, const TaylorFunction& p1, uint k)
{
    if(p1.smoothness()==0) {
        ARIADNE_THROW(std::runtime_error,"derivative(TaylorFunction,uint)"," model has smoothness 0");
    }

    p0.resize(p1.result_size(),p1.argument_size(),p1.order()-1,p1.smoothness()-1);

    MultiIndex dj(p1.argument_size());

    for(uint i=0; i!=p0.result_size(); ++i) {
        for(MultiIndex j(p1.argument_size()); j.degree()<=p1.order(); ++j) {
            if(j[k]!=0) {
                dj=j;
                dj.decrement_index(k);
                p0.at(i,dj)+=static_cast<int>(j[k])*p1.get(i,j);
            }
        }
    }
}

*/

 /*
Matrix<Interval>
TaylorFunction::jacobian(const Vector<Float>& s) const
{
    return this->jacobian(Vector<Interval>(s));
}


Matrix<Interval>
TaylorFunction::jacobian(const Vector<Interval>& x) const
{
    Matrix<Interval> J(this->result_size(),this->argument_size());
    Box w=x-this->centre();
    array< array<I> > powers=this->_powers(w);

    for(uint j=0; j!=this->argument_size(); ++j) {
        for(MultiIndex m(this->argument_size()); m.degree()<this->order(); ++m) {
            MultiIndex n=m;
            int c=n[j];
            if(c!=0) {
                n.decrement_index(j);
                I a=c;
                for(uint k=0; k!=this->argument_size(); ++k) {
                    a*=powers[k][n[k]];
                }
                for(uint i=0; i!=this->result_size(); ++i) {
                    I amim=a*this->_centre_derivatives[i][m];
                    J[i][j]+=amim;
                }
            }
        }
        for(MultiIndex m=MultiIndex::first(this->argument_size(),this->order()); m.degree()<=this->order(); ++m) {
            MultiIndex n=m;
            int c=n[j];
            if(c!=0) {
                n.decrement_index(j);
                I a=c;
                for(uint k=0; k!=this->argument_size(); ++k) {
                    a*=powers[k][n[k]];
                }
                for(uint i=0; i!=this->result_size(); ++i) {
                    J[i][j]+=a*this->_domain_derivatives[i][m];
                }
            }
        }
    }

    return J;
}
 */




TaylorFunction
inverse(const TaylorFunction& p, const Point& v)
{
    assert(p.result_size()==p.argument_size());
    assert(p.argument_size()==v.size());

    // The following are only to simplfy testing.
    assert(v==Point(v.size(),0));
    assert(p.evaluate(v)==Point(p.result_size(),0));
    typedef Interval I;

    Point c=midpoint(Box(p.evaluate(v)));
    Matrix<Float> J=p.jacobian(v);

    Matrix<Float> invJ=inverse(J);

    // FIXME: Need to re-solve for image of centre. What should initial set be? Different code needed for Rational?
    Box invf=v;

    // FIXME: Give correct initial conditions
    TaylorFunction result;

    for(MultiIndex m(p.result_size()); m.degree()<=2; ++m) {
        if(m.degree()==0) {
            for(uint i=0; i!=p.argument_size(); ++i) {
                result._variables[i][m]=v[i];
            }
        } else if(m.degree()==1) {
            for(uint i=0; i!=p.argument_size(); ++i) {
                result._variables[i][m]=invJ[i][m.position()-1];
            }
        } else {
            // FIXME: Add code for higher indices
        }
    }
    return result;
}



inline TaylorFunction prod(const Matrix<Interval>& A, const TaylorFunction& f) {
    return operator*(A,f);
}

TaylorFunction
implicit(const TaylorFunction& f)
{
    ARIADNE_ASSERT(f.result_size()<=f.argument_size());
    // Solve the equation f(x,h(x))=0
    // Use D1f + D2f Dh = 0, so Dh=-D2f^-1 D1f
    uint rs=f.result_size();
    uint fas=f.argument_size();
    uint has=fas-rs;

    Vector<Interval> h_domain=project(f.domain(),range(0u,has));
    Vector<Interval> h_range=project(f.domain(),range(has,fas));
    TaylorFunction id=TaylorFunction::identity(h_domain);
    TaylorFunction h=TaylorFunction::constant(h_domain,h_range);
    //std::cerr<<"\nid="<<id<<"\nh0="<<h0<<"\n";

    for(uint k=0; k!=10; ++k) {
        Vector<Interval> ih_range=join(h.domain(),h.range());
        Matrix<Interval> Df=f.jacobian(ih_range);
        //std::cerr<<"  Df="<<Df<<std::endl;
        Matrix<Interval> D2f=project(Df,range(0,rs),range(has,fas));
        //std::cerr<<"  D2f="<<J<<std::endl;
        Matrix<Interval> D2finv=inverse(D2f);

        for(uint i=0; i!=rs; ++i) {
            const_cast<TaylorVariable&>(h.variables()[i]).set_error(0);
        }
        TaylorFunction idh=join(id,h);
        TaylorFunction fidxhx=compose(f,idh);
        //std::cerr<<"  f(x,h(x))="<<fh<<std::endl;
        TaylorFunction dh=prod(D2finv,fidxhx);
        //std::cerr<<"  dh="<<dh<<std::endl;
        h=h-dh;
    }
    //std::cerr<<"\n  f="<<f<<"\n  h[0]="<<h0<<"\n  h[1]="<<h1<<"\n\n";
    ARIADNE_ASSERT(h.result_size()==f.result_size());
    ARIADNE_ASSERT(h.argument_size()+h.result_size()==f.argument_size());
    return h;

}





std::pair< Interval, Vector<Interval> >
bounds(Vector<TaylorVariable> const& vfm,
       Vector<Interval> const& d,
       Float const& hmax,
       Vector<Interval> dmax)
{
    ARIADNE_ASSERT(vfm.result_size()==vfm.argument_size());
    ARIADNE_ASSERT(vfm.result_size()==d.size());
    ARIADNE_ASSERT(vfm.result_size()==dmax.size());
    // Try to find a time h and a set b such that inside(r+Interval<R>(0,h)*vf(b),b) holds

    // Set up constants of the method.
    // TODO: Better estimates of constants
    const Float INITIAL_MULTIPLIER=2;
    const Float MULTIPLIER=1.125;
    //const Float BOX_RADIUS_MULTIPLIER=1.03125;
    const uint EXPANSION_STEPS=8;
    const uint REDUCTION_STEPS=8;
    const uint REFINEMENT_STEPS=4;

    Vector<Interval> delta=d-midpoint(d);

    Float h=hmax;
    Float hmin=hmax/(1<<REDUCTION_STEPS);
    bool success=false;
    Vector<Interval> b,nb,df;
    Interval ih(0,h);
    while(!success) {
        ARIADNE_ASSERT(h>hmin);
        df=evaluate(vfm,d);
        b=d+(INITIAL_MULTIPLIER*ih)*df+delta;
        for(uint i=0; i!=EXPANSION_STEPS; ++i) {
            df=evaluate(vfm,b);
            nb=d+ih*df;
            if(subset(nb,b)) {
                success=true;
                break;
            } else {
                b=d+MULTIPLIER*ih*df+delta;
            }
        }
        if(!success) {
            h/=2;
            ih=Interval(0,h);
        }
    }

    ARIADNE_ASSERT(subset(nb,b));

    Vector<Interval> vfb;
    vfb=evaluate(vfm,b);

    for(uint i=0; i!=REFINEMENT_STEPS; ++i) {
        b=nb;
        vfb=evaluate(vfm,b);
        nb=d+ih*vfb;
        ARIADNE_ASSERT_MSG(subset(nb,b),std::setprecision(20)<<"refinement "<<i<<": "<<nb<<" is not a inside of "<<b);
    }

    // Check result of operation
    // We use "possibly" here since the bound may touch
    ARIADNE_ASSERT(subset(nb,b));

    ARIADNE_ASSERT(b.size()==vfm.size());
    return std::make_pair(h,nb);
}

Vector<TaylorVariable>
flow(const Vector<TaylorVariable>& vf, const Vector<Interval>& d, const Interval& h, const Vector<Interval>& b)
{
    ARIADNE_ASSERT(vf.result_size()==vf.argument_size());
    ARIADNE_ASSERT(vf.size()==d.size());
    ARIADNE_ASSERT(vf.size()==b.size());

    ARIADNE_ASSERT(h.l==0.0 || h.l==-h.u);

    typedef Interval I;
    typedef Float A;
    uint n=vf.size();
    uint so=1u;
    uint to=6u;

    for(uint i=0; i!=n; ++i) { const_cast<TaylorVariable&>(vf[i]).clobber(so,to-1); }

    Interval hh(-h.u,h.u);

    Vector<TaylorVariable> y(n,TaylorVariable(n+1));
    Vector<TaylorVariable> yp(n,TaylorVariable(n+1));
    Vector<TaylorVariable> yz(n,TaylorVariable(n+1));
    for(uint i=0; i!=n; ++i) {
        yz[i].set_gradient(i,1.0);
        yz[i]*=sub_ivl(d[i].u/2,d[i].l/2);
        yz[i]+=add_ivl(d[i].l/2,d[i].u/2);
    }



    for(uint i=0; i!=n; ++i) {
        y[i]=yz[i];
    }


    //std::cerr << "\ny[0]=" << y << std::endl << std::endl;
    for(uint j=0; j!=to; ++j) {
        yp=compose(vf,b,y);
        //std::cerr << "yp["<<j+1<<"]=" << yp << std::endl;
        //for(uint i=0; i!=n; ++i) { yp[i].clobber(so,to-1); }
        //std::cerr << "yp["<<j+1<<"]=" << yp << std::endl;
        for(uint i=0; i!=n; ++i) {
            y[i]=antiderivative(yp[i],hh,n);
            y[i]+=yz[i];
        }
        //std::cerr << "y["<<j+1<<"]=" << y << std::endl << std::endl;
    }

    if(h.l==0) {
        Vector<Interval> d(n+1,Interval(-1,+1)); d[n]=Interval(0,1);
        Vector<TaylorVariable> ht=TaylorVariable::scaling(d);
        y=compose(y,Vector<Interval>(n+1,Interval(-1,1)),ht);
    }

    //for(uint i=0; i!=n; ++i) { y[i].clobber(so,to); }
    for(uint i=0; i!=n; ++i) { y[i].sweep(0.0); }

    ARIADNE_ASSERT(y.result_size()==vf.result_size());
    ARIADNE_ASSERT(y.argument_size()==vf.argument_size()+1);
    return y;
}



TaylorFunction
flow(const TaylorFunction& p, const Vector<Interval>& domain, const Interval& time)
{
    return TaylorFunction(flow(p.variables(),domain,time,p.domain()));
}


bool
refines(const TaylorFunction& f, const TaylorFunction& g)
{
    if(f.domain()==g.domain()) {
        return refines(f.variables(),g.variables());
    } else if(subset(g.domain(),f.domain())) {
        return refines(restrict(f,g.domain()).variables(),g.variables());
    } else {
        return false;
    }
}


array< array<Interval> >
TaylorFunction::_powers(const Vector<Interval>& v) const
{
    uint order=21;
    array< array<I> > powers(this->argument_size(), array<I>(order));
    for(uint i=0; i!=this->argument_size(); ++i) {
        powers[i][0]=1;
        if(order>=1) {
            powers[i][1]=v(i);
            if(order>=2) {
                powers[i][2]=pow(v(i),2);
                for(uint j=3; j<=order; ++j) {
                    powers[i][j]=powers[i][2]*powers[i][j-2];
                }
            }
        }
    }
    return powers;
}



std::ostream&
TaylorFunction::write(std::ostream& os) const
{
    os << "TaylorFunction( "<<this->_domain<<" , ";
    for(uint i=0; i!=this->result_size(); ++i) {
        os << (i==0?'[':',')<<this->_variables[i].expansion()<<","<<this->_variables[i].error();
    }
    return os << "] )";
    /*
    os << "TaylorFunction(\n";
    os << "  domain=" << this->domain() << ",\n" << std::flush;
    os << "  range=" << this->range() << ",\n" << std::flush;
    os << "  model=" << this->_variables << "\n" << std::flush;
    os << ")\n";
    return os;
    */
}



std::ostream&
operator<<(std::ostream& os, const TaylorFunction& p)
{
    return p.write(os);
}



/*
latexstream&
operator<<(Output::latexstream& texs, const TaylorFunction& p)
{
    using namespace Function;
    texs << "%TaylorFunction\n";
    texs << "\\ensuremath{\n";
    texs << "\\left( \\begin{array}{c}\n";
    char var='x';
    for(uint i=0; i!=p.result_size(); ++i) {
        bool first = true;
        if(i!=0) { texs << "\\\\"; }
        for(MultiIndex j(p.argument_size()); j.degree()<=p.order(); ++j) {
            const Interval& a=p.centre_derivatives()[i][j];
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
*/



TaylorExpression::TaylorExpression()
    : _domain(), _model() { }

TaylorExpression::TaylorExpression(const Vector<Interval>& domain)
    : _domain(domain), _model(domain.size()) { }

TaylorExpression::TaylorExpression(const Vector<Interval>& domain, const TaylorVariable& model)
    : _domain(domain), _model(model) { }

TaylorExpression TaylorExpression::variable(const Vector<Interval>& domain, Nat i) {
    TaylorExpression t(domain); t._model=TaylorVariable::scaling(domain,i); return t; }

Vector<TaylorExpression> TaylorExpression::variables(const Vector<Interval>& domain) {
    Vector<TaylorExpression> r(domain.size(),TaylorExpression(domain));
    for(uint i=0; i!=r.size(); ++i) { r[i]=variable(domain,i); }
    return r; }

const Vector<Interval>& TaylorExpression::domain() const {
    return  this->_domain; }
const TaylorVariable& TaylorExpression::model() const {
    return  this->_model; }

TaylorExpression operator+(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c+t._model); }
TaylorExpression operator-(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c-t._model); }
TaylorExpression operator*(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c*t._model); }
TaylorExpression operator/(const Float& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c/t._model); }

TaylorExpression operator+(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c+t._model); }
TaylorExpression operator-(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c-t._model); }
TaylorExpression operator*(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c*t._model); }
TaylorExpression operator/(const Interval& c, const TaylorExpression& t) {
    return TaylorExpression(t._domain,c/t._model); }

TaylorExpression operator+(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model+c); }
TaylorExpression operator-(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model-c); }
TaylorExpression operator*(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model*c); }
TaylorExpression operator/(const TaylorExpression& t, const Float& c) {
    return TaylorExpression(t._domain,t._model/c); }

TaylorExpression operator+(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model+c); }
TaylorExpression operator-(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model-c); }
TaylorExpression operator*(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model*c); }
TaylorExpression operator/(const TaylorExpression& t, const Interval& c) {
    return TaylorExpression(t._domain,t._model/c); }

TaylorExpression operator+(const TaylorExpression& t) {
    return TaylorExpression(t._domain,+t._model);  }
TaylorExpression operator-(const TaylorExpression& t) {
    return TaylorExpression(t._domain,-t._model);  }
TaylorExpression operator+(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model+t2._model); }
TaylorExpression operator-(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model-t2._model); }
TaylorExpression operator*(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model*t2._model); }
TaylorExpression operator/(const TaylorExpression& t1, const TaylorExpression& t2) {
    ARIADNE_ASSERT(t1._domain==t2._domain); return TaylorExpression(t1._domain,t1._model/t2._model); }

TaylorFunction operator*(const TaylorExpression& t, const Vector<Float>& v) {
    return TaylorFunction(Vector<TaylorVariable>(t._model*v)); }
TaylorFunction operator*(const Vector<Float>& v, const TaylorExpression& t) {
    return t*v; };

std::ostream& operator<<(std::ostream& os, const TaylorExpression& t) {
    return os<<"("<<t._domain<<","<<t._model<<")"; }

} // namespace Ariadne
