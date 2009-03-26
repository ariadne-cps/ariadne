/***************************************************************************
 *            taylor_variable.cc
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
 *  GNU Library General Public License for more detai1ls.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#include <iomanip>

#include "config.h"
#include "rounding.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "expansion.h"
#include "differential.h"
#include "taylor_model.h"
#include "taylor_variable.h"
#include "exceptions.h"

namespace Ariadne {



TaylorVariable::TaylorVariable()
    : _domain(0), _model(0)
{ }

TaylorVariable::TaylorVariable(const DomainType& d)
    : _domain(d), _model(d.size())
{
}

TaylorVariable::TaylorVariable(const DomainType& d, const TaylorModel& m)
    : _domain(d), _model(m)
{
}

TaylorVariable::TaylorVariable(const DomainType& d, const ExpansionType& f, const ErrorType& e)
    : _domain(d), _model(f,e)
{
    ARIADNE_ASSERT(d.size()==f.argument_size());
}



void _set_scaling(TaylorVariable& x, const Interval& ivl, uint j)
{
    rounding_mode_t rounding_mode=get_rounding_mode();
    set_rounding_mode(upward);
    const double& l=ivl.lower();
    const double& u=ivl.upper();
    volatile double pc=u; pc+=l;
    volatile double nc=-u; nc-=l;
    volatile double pg=u; pg-=l;
    volatile double ng=l; ng-=u;
    x.error()=(pc+nc+pg+ng)/4;
    set_rounding_mode(to_nearest);
    MultiIndex a(x.argument_size());
    x.expansion().append(a,(l+u)/2);
    ++a[j];
    x.expansion().append(a,(l+u)/2);
    set_rounding_mode(rounding_mode);
}

TaylorVariable TaylorVariable::zero(const DomainType& d)
{
    TaylorVariable x(d);
    x.set_value(0.0);
    return x;
}

TaylorVariable TaylorVariable::constant(const DomainType& d, const Float& c)
{
    TaylorVariable x(d);
    x.set_value(c);
    return x;
}

TaylorVariable TaylorVariable::variable(const DomainType& d, unsigned int j)
{
    TaylorVariable x(d);
    _set_scaling(x,d[j],j);
    return x;
}

TaylorVariable TaylorVariable::scaling(const Interval& d, const Interval& r)
{
    TaylorVariable x(Vector<Interval>(1u,d));
    _set_scaling(x,r,0);
    return x;
}

TaylorVariable TaylorVariable::scaling(const Vector<Interval>& d, const Interval&, uint j)
{
    ARIADNE_ASSERT(j<d.size());
    TaylorVariable x(d);
    _set_scaling(x,d[j],j);
    return x;
}

Vector<TaylorVariable> TaylorVariable::constants(const Vector<Interval>& d, const Vector<Float>& c)
{
    Vector<TaylorVariable> x(c.size());
    for(uint i=0; i!=c.size(); ++i) {
        x[i]=TaylorVariable::constant(d,c[i]);
    }
    return x;
}

Vector<TaylorVariable> TaylorVariable::constants(const Vector<Interval>& d, const Vector<Interval>& c)
{
    Vector<TaylorVariable> x(c.size(),TaylorVariable(d));
    for(uint i=0; i!=c.size(); ++i) {
        x[i]=c[i];
    }
    return x;
}

Vector<TaylorVariable> TaylorVariable::variables(const Vector<Interval>& d)
{
    Vector<TaylorVariable> x(d.size());
    for(uint i=0; i!=d.size(); ++i) {
        x[i]=TaylorVariable::variable(d,i);
    }
    return x;
}

Vector<TaylorVariable> TaylorVariable::scalings(const Vector<Interval>& d, const Vector<Interval>& r)
{
    ARIADNE_ASSERT(d.size()==r.size());
    Vector<TaylorVariable> x(r.size(),TaylorVariable(d));
    for(uint i=0; i!=d.size(); ++i) {
        _set_scaling(x[i],r[i],i);
    }
    return r;
}


bool TaylorVariable::operator==(const TaylorVariable& tv) const
{
    return this->_domain==tv._domain && this->_model==tv._model;
}


TaylorVariable& operator+=(TaylorVariable& x, const TaylorVariable& y) {
    ARIADNE_ASSERT(subset(x.domain(),y.domain()));
    if(x.domain()==y.domain()) { x._model+=y._model; }
    else { x._model+=restrict(y,x.domain())._model; }
    return x;
}

TaylorVariable& operator-=(TaylorVariable& x, const TaylorVariable& y) {
    ARIADNE_ASSERT(subset(x.domain(),y.domain()));
    if(x.domain()==y.domain()) { x._model-=y._model; }
    else { x._model-=restrict(y,x.domain())._model; }
    return x;
}


TaylorVariable operator+(const TaylorVariable& x1, const TaylorVariable& x2) {
    if(x1._domain==x2._domain) {
        return TaylorVariable(x1._domain,x1._model+x2._model); }
    else {
        TaylorVariable::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorVariable(domain,restrict(x1,domain)._model+restrict(x2,domain)._model);}
}

TaylorVariable operator-(const TaylorVariable& x1, const TaylorVariable& x2) {
    if(x1._domain==x2._domain) {
        return TaylorVariable(x1._domain,x1._model-x2._model); }
    else {
        TaylorVariable::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorVariable(domain,restrict(x1,domain)._model-restrict(x2,domain)._model);}
}

TaylorVariable operator*(const TaylorVariable& x1, const TaylorVariable& x2) {
    if(x1._domain==x2._domain) {
        return TaylorVariable(x1._domain,x1._model*x2._model); }
    else {
        TaylorVariable::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorVariable(domain,restrict(x1,domain)._model*restrict(x2,domain)._model);}
}

TaylorVariable operator/(const TaylorVariable& x1, const TaylorVariable& x2) {
    if(x1._domain==x2._domain) {
        return TaylorVariable(x1._domain,x1._model/x2._model); }
    else {
        TaylorVariable::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorVariable(domain,restrict(x1,domain)._model/restrict(x2,domain)._model);}
}





TaylorVariable max(const TaylorVariable& x1, const TaylorVariable& x2) {
    if(x1._domain==x2._domain) {
        return TaylorVariable(x1._domain,max(x1._model,x2._model)); }
    else {
        TaylorVariable::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorVariable(domain,max(restrict(x1,domain)._model,restrict(x2,domain)._model));}
}

TaylorVariable min(const TaylorVariable& x1, const TaylorVariable& x2) {
    if(x1._domain==x2._domain) {
        return TaylorVariable(x1._domain,min(x1._model,x2._model)); }
    else {
        TaylorVariable::DomainType domain=intersection(x1._domain,x2._domain);
        return TaylorVariable(domain,min(restrict(x1,domain)._model,restrict(x2,domain)._model));}
}





Interval
TaylorVariable::evaluate(const Vector<Float>& v) const
{
    ARIADNE_NOT_IMPLEMENTED;
}

Interval
TaylorVariable::evaluate(const Vector<Interval>& v) const
{
    ARIADNE_NOT_IMPLEMENTED;
}



TaylorVariable restrict(const TaylorVariable& tv, const Vector<Interval>& d) {
    ARIADNE_NOT_IMPLEMENTED;
}

pair<TaylorVariable,TaylorVariable>
split(const TaylorVariable& tv, uint j)
{
    ARIADNE_NOT_IMPLEMENTED;
}

bool refines(const TaylorVariable& tv1, const TaylorVariable& tv2)
{
    if(tv1.domain()==tv2.domain()) { return refines(tv1.model(),tv2.model()); } if(subset(tv2.domain(),tv1.domain())) { return refines(restrict(tv1,tv2.domain()).model(),tv2.model()); }
    else { return false; }
}

Vector<TaylorVariable> compose(const Vector<TaylorVariable>& x, const Vector<TaylorVariable>& y) {
    ARIADNE_NOT_IMPLEMENTED; }
TaylorVariable compose(const TaylorVariable& x, const Vector<TaylorVariable>& y) {
    ARIADNE_NOT_IMPLEMENTED; }

Vector<TaylorVariable>
prod(const Matrix<Interval>& A,
     const Vector<TaylorVariable>& x)
{
    ARIADNE_ASSERT(x.size()>0);
    ARIADNE_ASSERT(A.column_size()==x.size());
    for(uint i=1; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].argument_size()==x[0].argument_size()); }

    Vector<TaylorVariable> r(A.row_size(),TaylorVariable(x[0].domain()));
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint j=0; j!=A.column_size(); ++j) {
            r[i]+=A[i][j]*x[j];
        }
    }
    return r;
}



Vector<TaylorVariable>
implicit(const Vector<TaylorVariable>& f)
{
    ARIADNE_ASSERT(f.result_size()<=f.argument_size());
    // Solve the equation f(x,h(x))=0
    // Use D1f + D2f Dh = 0, so Dh=-D2f^-1 D1f
    uint rs=f.result_size();
    uint fas=f.argument_size();
    uint has=fas-rs;

    Vector<Interval> h_domain=project(f.domain(),range(0u,has));
    Vector<Interval> h_range=project(f.domain(),range(has,fas));
    Vector<TaylorVariable> id=TaylorVariable::variables(h_domain);
    Vector<TaylorVariable> h=TaylorVariable::constants(h_domain,h_range);
    //std::cerr<<"\nid="<<id<<"\nh0="<<h0<<"\n";

    for(uint k=0; k!=10; ++k) {
        Vector<Interval> ih_range=join(h_domain,h_range);
        Matrix<Interval> Df=f.jacobian(ih_range);
        //std::cerr<<"  Df="<<Df<<std::endl;
        Matrix<Interval> D2f=project(Df,range(0,rs),range(has,fas));
        //std::cerr<<"  D2f="<<J<<std::endl;
        Matrix<Interval> D2finv=inverse(D2f);

        for(uint i=0; i!=rs; ++i) {
            h[i].set_error(0);
        }
        Vector<TaylorVariable> idh=join(id,h);
        Vector<TaylorVariable> fidxhx=compose(f,idh);
        //std::cerr<<"  f(x,h(x))="<<fh<<std::endl;
        Vector<TaylorVariable> dh=prod(D2finv,fidxhx);
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
flow(const Vector<TaylorVariable>& vf, const Vector<Interval>& d, const Interval& h)
{
    ARIADNE_ASSERT(vf.result_size()==vf.argument_size());
    ARIADNE_ASSERT(vf.size()==d.size());

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
        yz[i][MultiIndex::unit(n,i)]=1.0;
        yz[i]*=sub_ivl(d[i].u/2,d[i].l/2);
        yz[i]+=add_ivl(d[i].l/2,d[i].u/2);
    }



    for(uint i=0; i!=n; ++i) {
        y[i]=yz[i];
    }


    //std::cerr << "\ny[0]=" << y << std::endl << std::endl;
    for(uint j=0; j!=to; ++j) {
        yp=compose(vf,y);
        //std::cerr << "yp["<<j+1<<"]=" << yp << std::endl;
        //for(uint i=0; i!=n; ++i) { yp[i].clobber(so,to-1); }
        //std::cerr << "yp["<<j+1<<"]=" << yp << std::endl;
        for(uint i=0; i!=n; ++i) {
            y[i]=antiderivative(yp[i],n);
            y[i]+=yz[i];
        }
        //std::cerr << "y["<<j+1<<"]=" << y << std::endl << std::endl;
    }

    if(h.l==0) {
        Vector<Interval> dom(n+1,Interval(-1,+1));
        Vector<Interval> rng(n+1,Interval(-1,+1)); rng[n]=Interval(0,1);
        Vector<TaylorVariable> ht=TaylorVariable::scalings(dom,rng);
        y=compose(y,ht);
    }

    //for(uint i=0; i!=n; ++i) { y[i].clobber(so,to); }
    for(uint i=0; i!=n; ++i) { y[i].sweep(0.0); }

    ARIADNE_ASSERT(y.result_size()==vf.result_size());
    ARIADNE_ASSERT(y.argument_size()==vf.argument_size()+1);
    return y;
}


std::ostream&
operator<<(std::ostream& os, const TaylorVariable& tv) {
    return os << "TaylorVariable(" << tv.domain() << "," << tv.expansion() << "," << tv.error() << ")";
}



Vector< Expansion<Float> >
Vector<TaylorVariable>::expansion() const
{
    Vector< Expansion<Float> > r(this->size());
    for(uint i=0; i!=this->size(); ++i) {
        r[i]=(*this)[i].expansion();
    }
    return r;
}

Vector<Float>
Vector<TaylorVariable>::error() const
{
    Vector<Float> r(this->size());
    for(uint i=0; i!=this->size(); ++i) {
        r[i]=(*this)[i].error();
    }
    return r;
}

Vector<Float>
Vector<TaylorVariable>::value() const
{
    Vector<Float> r(this->size());
    for(uint i=0; i!=this->size(); ++i) {
        r[i]=(*this)[i].value();
    }
    return r;
}

Vector<Interval>
Vector<TaylorVariable>::range() const
{
    Vector<Interval> r(this->size());
    for(uint i=0; i!=this->size(); ++i) {
        r[i]=(*this)[i].range();
    }
    return r;
}

Matrix<Interval>
Vector<TaylorVariable>::jacobian(const Vector<Interval>& x) const
{
    Vector< Differential<Interval> > s(this->argument_size(),this->argument_size(),1u);
    for(uint j=0; j!=this->argument_size(); ++j) {
        Interval dj=this->domain()[j];
        s[j].set_value((x[j]-dj.midpoint())/dj.radius());
        s[j].set_gradient(j,1/dj.radius());
    }
    Vector< Expansion<Float> > p=this->expansion();
    Vector< Differential<Interval> > d=Ariadne::evaluate(p,s);
    //std::cerr<<"  x="<<x<<"\n  p="<<p<<"\n"<<"  s="<<s<<"\n  p.s="<<d<<"\n  J="<<d.jacobian()<<"\n"<<std::endl;
    return d.jacobian();
}


void Vector<TaylorVariable>::check() const
{
    for(uint i=0; i!=this->size(); ++i) {
        ARIADNE_ASSERT((*this)[0].argument_size()==(*this)[i].argument_size());
    }
    for(uint i=0; i!=this->size(); ++i) {
        ARIADNE_ASSERT((*this)[0].domain()==(*this)[i].domain());
    }
}


} //namespace Ariadne


