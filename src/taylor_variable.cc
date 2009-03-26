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
#include "taylor_variable.h"
#include "exceptions.h"

namespace Ariadne {


const double TaylorVariable::em=2.2204460492503131e-16;
const double TaylorVariable::ec=em/2;

template<class X> X Expansion<X>::_zero = 0.0;

double TaylorVariable::_default_sweep_threshold=1e-18;
uint TaylorVariable::_default_maximum_degree=16;



TaylorVariable::TaylorVariable()
    : _domain(0), _expansion(0), _error(0),
      _sweep_threshold(_default_sweep_threshold),
      _maximum_degree(_default_maximum_degree),
      _maximum_index(0,_default_maximum_degree)
{ }

TaylorVariable::TaylorVariable(const DomainType& d)
    : _domain(d), _expansion(d.size()), _error(0),
      _sweep_threshold(_default_sweep_threshold),
      _maximum_degree(_default_maximum_degree),
      _maximum_index(d.size(),_default_maximum_degree)
{
}

TaylorVariable::TaylorVariable(const DomainType& d, const ExpansionType& f, const ErrorType& e)
    : _domain(d), _expansion(f), _error(e),
      _sweep_threshold(_default_sweep_threshold),
      _maximum_degree(_default_maximum_degree),
      _maximum_index(d.size(),_default_maximum_degree)
{
    ARIADNE_ASSERT(f.empty() || d.size()==f.argument_size());
    ARIADNE_ASSERT(this->_error>=0);
}

/*
TaylorVariable::TaylorVariable(uint as, uint deg, const double* ptr, const double& err)
    : _expansion(as), _error(err),
      _sweep_threshold(_default_sweep_threshold),
      _maximum_degree(_default_maximum_degree),
      _maximum_index(as,_default_maximum_degree)
{
    for(MultiIndex j(as); j.degree()<=deg; ++j) {
        if(*ptr!=0.0 || j.degree()<=1) { this->_expansion.append(j,*ptr); }
        ++ptr;
    }
    this->_expansion.cleanup();
    if(err<0) { std::cerr<<err<<std::endl; } ARIADNE_ASSERT(this->_error>=0);
}

TaylorVariable::TaylorVariable(uint as, uint deg, double d0, ...)
    : _expansion(as), _error(0),
      _sweep_threshold(_default_sweep_threshold),
      _maximum_degree(_default_maximum_degree),
      _maximum_index(as,_default_maximum_degree)
{
    double x=d0;
    va_list args;
    va_start(args,d0);
    for(MultiIndex j(as); j.degree()<=deg; ++j) {
        if(x!=0.0 || j.degree()<=1) { this->_expansion.append(j,x); }
        x=va_arg(args,double);
    }
    this->_error=x;
    va_end(args);
    this->_expansion.cleanup();
    ARIADNE_ASSERT(this->_error>=0);
}

TaylorVariable::TaylorVariable(uint as, uint nnz, uint a0, ...)
    : _expansion(as), _error(0),
      _sweep_threshold(_default_sweep_threshold),
      _maximum_degree(_default_maximum_degree),
      _maximum_index(as,_default_maximum_degree)
{
    MultiIndex a(as);
    uint aj=0u;
    double x=0.0;
    va_list args;
    va_start(args,a0);
    for(uint k=0; k!=nnz; ++k) {
        for(uint j=0; j!=as; ++j) {
            if(j!=0 || k!=0) { aj=va_arg(args,uint); }
            a[j]=aj;
        }
        x=va_arg(args,double);
        this->_expansion[a]=x;
    }
    double e=va_arg(args,double);
    this->error()=e;
    va_end(args);
    this->_expansion.cleanup();
    ARIADNE_ASSERT(this->_error>=0);
}
*/

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





TaylorVariable&
TaylorVariable::operator=(const Float& c)
{
    this->_expansion.clear();
    this->_expansion[MultiIndex::zero(this->argument_size())]=c;
    this->_error=0;
    return *this;
}

TaylorVariable&
TaylorVariable::operator=(const Interval& c)
{
    this->_expansion.clear();
    this->_expansion[MultiIndex::zero(this->argument_size())]=c.midpoint();
    this->_error=c.radius();
    return *this;
}

bool TaylorVariable::operator==(const TaylorVariable& tv) const
{
    return this->_domain==tv._domain && this->_expansion==tv._expansion && this->_error == tv._error;
}


void scal(TaylorVariable& r, const Float& c)
{
    //std::cerr<<"TaylorVariable::scal(Float c) c="<<c<<std::endl;
    Float& re=r.error();

    // Shortcuts for special cases
    if(c==0.0) {
        r.expansion().clear();
        re=0;
        return;
    }
    if(c==1.0) {
        return;
    }
    if(c==0.5 || c==2.0 || c==-2.0 || c==-1.0 || c==-0.5) {
        // Operation can be performed exactly
        for(TaylorVariable::iterator iter=r.begin(); iter!=r.end(); ++iter) {
            iter->second*=c;
        }
        re*=abs(c);
        return;
    }

    // General case with error analysis
    set_rounding_upward();
    volatile Float te=0; // Twice the maximum accumulated error
    for(TaylorVariable::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        volatile Float u=riter->second*c;
        volatile Float t=-riter->second;
        volatile Float ml=t*c;
        te+=(u+ml);
    }
    re*=abs(c);
    re+=te/2;

    set_rounding_to_nearest();
    for(TaylorVariable::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->second*=c;
    }
    return;
}



void scal(TaylorVariable& r, const Interval& c)
{
    //std::cerr<<"TaylorVariable::scal(Interval c) c="<<c<<std::endl;
    Float& re=r.error();
    set_rounding_upward();
    volatile Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    for(TaylorVariable::const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        const Float& rv=riter->second;
        volatile Float mrv=-rv;
        if(rv>=0) {
            u=rv*c.u;
            ml=mrv*c.l;
        } else {
            u=rv*c.l;
            ml=mrv*c.u;
        }
        te+=(u+ml);
    }
    re*=mag(c);
    re+=te/2;

    set_rounding_to_nearest();
    Float m=(c.u+c.l)/2;
    for(TaylorVariable::iterator riter=r.begin(); riter!=r.end(); ++riter) {
        riter->second*=m;
    }
    return;
}



void acc(TaylorVariable& r, const Float& c)
{
    // Compute self+=c
    Float& rv=r.value();
    Float& re=r.error();
    set_rounding_upward();
    volatile Float rvu=rv+c;
    volatile Float mrvl=(-rv)-c;
    //std::cerr<<"re="<<re.u<<" ";
    re+=(rvu+mrvl)/2;
    //std::cerr<<"nre="<<re.u<<"\n";
    set_rounding_to_nearest();
    rv+=c;
    return;
}



void acc(TaylorVariable& r, const Interval& c)
{
    // Compute self+=c

    Float& rv=r.value();
    Float& re=r.error();
    set_rounding_upward();
    volatile Float rvu=rv+c.u;
    volatile Float mrvl=(-rv)-c.l;
    re+=(rvu+mrvl)/2;
    set_rounding_to_nearest();
    volatile Float m=(c.u+c.l)/2;
    rv+=m;
    return;
}


void acc1(TaylorVariable& x, const TaylorVariable& y)
{
    // Compute self+=y
    Float& xe=x.error();

    set_rounding_upward();
    Float te=0;
    for(TaylorVariable::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        const Float& yv=yiter->second;;
        Float& xv=x[yiter->first];
        if(xv!=0) {
            volatile Float u=xv+yv;
            volatile Float t=-xv;
            volatile Float ml=t-yv;
            te+=(u+ml);
            //std::cerr<<" xv="<<xv<<" yv="<<yv<<" u="<<u<<" ml="<<ml<<" d="<<(u+ml)<<" te="<<te<<"\n";
        }
    }
    xe+=te/2;
    xe+=y.error();

    set_rounding_to_nearest();
    for(TaylorVariable::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        x[yiter->first]+=yiter->second;
    }

    return;
}


void acc2(TaylorVariable& x, const TaylorVariable& y)
{
    // Compute x+=y
    TaylorVariable r(x.argument_size());

    set_rounding_upward();
    Float te=0.0;
    TaylorVariable::const_iterator xiter=x.begin();
    TaylorVariable::const_iterator yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->first<yiter->first) {
            ++xiter;
        } else if(yiter->first<xiter->first) {
            ++yiter;
        } else {
            const Float& xv=xiter->second;
            const Float& yv=yiter->second;
            volatile Float u=xv+yv;
            volatile Float t=-xv;
            volatile Float ml=t-yv;
            te+=(u+ml);
            ++xiter; ++yiter;
        }
    }
    r.error()=x.error();
    r.error()+=y.error();
    r.error()+=(te/2);

    set_rounding_to_nearest();
    xiter=x.begin();
    yiter=y.begin();
    while(xiter!=x.end() && yiter!=y.end()) {
        if(xiter->first<yiter->first) {
            r.expansion().append(xiter->first,xiter->second);
            ++xiter;
        } else if(yiter->first<xiter->first) {
            r.expansion().append(yiter->first,yiter->second);
            ++yiter;
        } else {
            Float c=xiter->second+yiter->second;
            if(c!=0.0) { r.expansion().append(xiter->first,xiter->second+yiter->second); }
            ++xiter; ++yiter;
        }
    }
    while(xiter!=x.end()) {
        r.expansion().append(xiter->first,xiter->second);
        ++xiter;
    }
    while(yiter!=y.end()) {
        r.expansion().append(yiter->first,yiter->second);
        ++yiter;
    }

    x.expansion().swap(r.expansion());
    x.error()=r.error();

    return;
}


inline void acc(TaylorVariable& r, const TaylorVariable& x)
{
    acc2(r,x);
}


struct Ivl { double u; double ml; };

void mulacc1(TaylorVariable& r, const TaylorVariable& x, const TaylorVariable& y)
{
    // Compute self+=x*y
    typedef TaylorVariable::const_iterator const_iterator;
    typedef std::map<MultiIndex,Ivl>::const_iterator ivl_const_iterator;
    Float& re=r.error();
    std::map<MultiIndex,Ivl> z;

    set_rounding_upward();
    for(const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        for(const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const Float& xv=xiter->second;;
            const Float& yv=yiter->second;;
            Ivl& zv=z[xiter->first+yiter->first];
            zv.u+=xv*yv;
            volatile double t=-xv;
            zv.ml+=t*yv;
        }
    }

    for(const_iterator riter=r.begin(); riter!=r.end(); ++riter) {
        Ivl& zv=z[riter->first];
        const Float& rv=riter->second;
        zv.u+=rv; zv.ml-=rv;
    }

    volatile Float te=0;
    for(ivl_const_iterator ziter=z.begin(); ziter!=z.end(); ++ziter) {
        const Ivl& zv=ziter->second;
        te+=(zv.u+zv.ml);
    }
    te/=2;

    Float xs=0;
    for(const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->second);
    }

    Float ys=0;
    for(const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->second);
    }

    const Float& xe=x.error();
    const Float& ye=y.error();

    re+=xs*ye+ys*xe+te+xe*ye;

    set_rounding_to_nearest();
    for(ivl_const_iterator ziter=z.begin(); ziter!=z.end(); ++ziter) {
        const Ivl& zv=ziter->second;
        r[ziter->first]=(zv.u-zv.ml)/2;
    }

    return;
}


void mulacc2(TaylorVariable& r, const TaylorVariable& x, const TaylorVariable& y)
{
    ARIADNE_ASSERT(r.begin()->first.degree()==0);
    TaylorVariable t(x.domain());
    for(TaylorVariable::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        set_rounding_upward();
        double te=0.0;
        for(TaylorVariable::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            const Float& xv=xiter->second;
            const Float& yv=yiter->second;
            volatile Float u=xv*yv;
            volatile Float ml=-xv; ml=ml*yv;
            te+=(u+ml);
        }
        t.error()=te/2;
        set_rounding_to_nearest();
        for(TaylorVariable::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            t.expansion().append(xiter->first,yiter->first,xiter->second*yiter->second);
        }
        r+=t;
        //std::cerr<<"  t="<<t<<"\n  r="<<r<<std::endl;
        t.expansion().clear();
        t.error()=0.0;
    }

    set_rounding_upward();
    Float xs=0;
    for(TaylorVariable::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        xs+=abs(xiter->second);
    }

    Float ys=0;
    for(TaylorVariable::const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
        ys+=abs(yiter->second);
    }

    Float& re=r.error();
    const Float& xe=x.error();
    const Float& ye=y.error();
    re+=xs*ye+ys*xe+xe*ye;

    r.clean();

    set_rounding_to_nearest();
    return;
}


inline void mulacc(TaylorVariable& r, const TaylorVariable& x, const TaylorVariable& y)
{
    mulacc2(r,x,y);
}


TaylorVariable&
TaylorVariable::clean()
{
    return this->clean(this->_maximum_degree,this->_sweep_threshold);
}

TaylorVariable&
TaylorVariable::clean(uint deg, double eps)
{
    TaylorVariable r(this->argument_size());
    r.error()=this->error();

    set_rounding_upward();
    for(TaylorVariable::const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->first.degree()==0 || (iter->first.degree() <= deg && abs(iter->second)>=eps)) {
            r.expansion().append(iter->first,iter->second);
        } else {
            r.error()+=abs(iter->second);
        }
    }
    set_rounding_to_nearest();

    this->expansion().swap(r.expansion());
    this->error()=r.error();

    return *this;
}

TaylorVariable&
TaylorVariable::sweep()
{
    return this->sweep(this->_sweep_threshold);
}

TaylorVariable&
TaylorVariable::sweep(double m)
{
    this->clean(255,m);
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        Float a=abs(iter->second);
        if(abs(iter->second)<=m) {
            this->_error+=a;
           _expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    return *this;
}

TaylorVariable&
TaylorVariable::truncate()
{
    return this->truncate(this->_maximum_degree);
}

TaylorVariable&
TaylorVariable::truncate(uint d)
{
    this->clean(d,0.0);
    set_rounding_upward();
    Float e=0;
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        if(iter->first.degree()>d) {
            e+=abs(iter->second);
            this->_expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    this->_error+=e;
    set_rounding_to_nearest();
    return *this;
}

TaylorVariable&
TaylorVariable::truncate(const MultiIndex& a)
{
    ARIADNE_ASSERT(a.size()==this->argument_size());
    set_rounding_upward();
    Float e=0;
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        bool erase=false;
        const MultiIndex& b=iter->first;
        for(uint i=0; i!=a.size(); ++i) {
            if(a[i] < b[i]) {
                erase=true;
                break;
            }
        }
        if(erase) {
            e+=abs(iter->second);
            this->_expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    this->_error+=e;
    set_rounding_to_nearest();
    return *this;
}

TaylorVariable&
TaylorVariable::truncate(const MultiIndexBound& b)
{
    ARIADNE_ASSERT(b.size()==this->argument_size());
    set_rounding_upward();
    Float e=0;
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        if(!(iter->first<=b)) {
            e+=abs(iter->second);
            this->_expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    this->_error+=e;
    set_rounding_to_nearest();
    return *this;
}

TaylorVariable&
TaylorVariable::clobber()
{
    this->_error=0;
    return *this;
}

TaylorVariable&
TaylorVariable::clobber(uint o)
{
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        if(iter->first.degree()>o) {
           _expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    this->_error=0;
    return *this;
}

TaylorVariable&
TaylorVariable::clobber(uint so, uint to)
{
    uint n=this->argument_size()-1;
    for(iterator iter=this->_expansion.begin(); iter!=this->end(); ) {
        const MultiIndex& a=iter->first;
        if(a[n]>to || a.degree()>so+a[n]) {
            this->_expansion.erase(iter++);
        } else {
            ++iter;
        }
    }
    this->_error=0;
    return *this;
}



TaylorVariable&
operator+=(TaylorVariable& x, const TaylorVariable& y)
{
    if(x.argument_size()==0) {
        x=TaylorVariable::constant(y.domain(),x.value());
    }
    ARIADNE_ASSERT_EQUAL(x.domain(),y.domain());
    acc(x,y); return x;

}

TaylorVariable&
operator-=(TaylorVariable& x, const TaylorVariable& y)
{
    ARIADNE_ASSERT_EQUAL(x.domain(),y.domain());
    acc(x,neg(y)); return x;
}


TaylorVariable&
operator+=(TaylorVariable& x, const Float& c)
{
    acc(x,c); return x;
}

TaylorVariable&
operator-=(TaylorVariable& x, const Float& c)
{
    acc(x,-c); return x;
}


TaylorVariable&
operator+=(TaylorVariable& x, const Interval& c)
{
    acc(x,c); return x;
}

TaylorVariable&
operator-=(TaylorVariable& x, const Interval& c)
{
    acc(x,-c); return x;
}

TaylorVariable&
operator*=(TaylorVariable& x, const Float& c)
{
    scal(x,c); return x;
}

TaylorVariable&
operator*=(TaylorVariable& x, const Interval& c)
{
    scal(x,c); return x;
}


TaylorVariable&
operator/=(TaylorVariable& x, const Float& c)
{
    if(c==0) {
        ARIADNE_THROW(DivideByZeroException,"operator/=(TaylorVariable x,Float c)","x="<<x<<" c="<<c);
    }
    scal(x,Interval(1.0)/c); return x;
}


TaylorVariable&
operator/=(TaylorVariable& x, const Interval& c)
{
    if(c.upper()>=0 && c.lower()<=0) {
        ARIADNE_THROW(DivideByZeroException,"operator/=(TaylorVariable x,Interval c)","x="<<x<<" c="<<c);
    }
    scal(x,1.0/c); return x;
}




TaylorVariable
operator+(const TaylorVariable& x) {
    return x;
}

TaylorVariable
operator-(const TaylorVariable& x) {
    return neg(x);
}


TaylorVariable
operator+(const TaylorVariable& x, const TaylorVariable& y) {
    ARIADNE_ASSERT(x.domain()==y.domain());
    TaylorVariable r(x); acc(r,y); return r;
}

TaylorVariable
operator-(const TaylorVariable& x, const TaylorVariable& y) {
    ARIADNE_ASSERT(x.domain()==y.domain());
    TaylorVariable r=neg(y); acc(r,x); return r;
}

TaylorVariable
operator*(const TaylorVariable& x, const TaylorVariable& y) {
    ARIADNE_ASSERT(x.domain()==y.domain());
    TaylorVariable r(x.domain()); r.set_value(0.0); mulacc(r,x,y); return r;
}

TaylorVariable
operator/(const TaylorVariable& x, const TaylorVariable& y) {
    ARIADNE_ASSERT(x.domain()==y.domain());
    TaylorVariable r(x.domain()); r.set_value(0.0); mulacc(r,x,rec(y)); return r;
}



TaylorVariable
operator+(const TaylorVariable& x, const Float& c) {
    TaylorVariable r(x); acc(r,c); return r;
}

TaylorVariable
operator-(const TaylorVariable& x, const Float& c) {
    TaylorVariable r(x); acc(r,-c); return r;
}

TaylorVariable
operator*(const TaylorVariable& x, const Float& c) {
    TaylorVariable r(x); scal(r,c); return r;
}

TaylorVariable
operator/(const TaylorVariable& x, const Float& c) {
    TaylorVariable r(x); scal(r,Interval(1)/c); return r;
}

TaylorVariable
operator+(const Float& c, const TaylorVariable& x) {
    TaylorVariable r(x); acc(r,c); return r;
}

TaylorVariable
operator-(const Float& c, const TaylorVariable& x) {
    TaylorVariable r=neg(x); acc(r,c); return r;
}

TaylorVariable
operator*(const Float& c, const TaylorVariable& x) {
    TaylorVariable r(x); scal(r,c); return r;
}

TaylorVariable
operator/(const Float& c, const TaylorVariable& x) {
    TaylorVariable r(x); scal(r,Interval(1)/c); return r;
}



TaylorVariable
operator+(const TaylorVariable& x, const Interval& c) {
    TaylorVariable r(x); acc(r,c); return r;
}

TaylorVariable
operator-(const TaylorVariable& x, const Interval& c) {
    TaylorVariable r(x); acc(r,-c); return r;
}

TaylorVariable
operator*(const TaylorVariable& x, const Interval& c) {
    TaylorVariable r(x); scal(r,c); return r;
}

TaylorVariable
operator/(const TaylorVariable& x, const Interval& c) {
    TaylorVariable r(x); scal(r,1/c); return r;
}

TaylorVariable
operator+(const Interval& c, const TaylorVariable& x) {
    TaylorVariable r(x); acc(r,c); return r;
}

TaylorVariable
operator-(const Interval& c, const TaylorVariable& x) {
    TaylorVariable r=neg(x); acc(r,c); return r;
}

TaylorVariable
operator*(const Interval& c, const TaylorVariable& x) {
    TaylorVariable r(x); scal(r,c); return r;
}

TaylorVariable
operator/(const Interval& c, const TaylorVariable& x) {
    TaylorVariable r=rec(x); scal(r,c); return r;
}





TaylorVariable max(const TaylorVariable& x, const TaylorVariable& y) {
    Interval xr=x.range();
    Interval yr=y.range();
    if(xr.lower()>=yr.upper()) {
        return x;
    } else if(yr.lower()>=xr.upper()) {
        return y;
    } else {
        TaylorVariable z(x.argument_size());
        z.value()=max(x.value(),y.value());
        z.error()=(max(xr.u,yr.u)-max(x.value(),y.value()));
        return z;
    }
}





TaylorVariable min(const TaylorVariable& x, const TaylorVariable& y) {
    return -max(-x,-y);
}

TaylorVariable abs(const TaylorVariable& x) {
    Interval xr=x.range();
    if(xr.lower()>=0.0) {
        return x;
    } else if(xr.upper()<=0.0) {
        return -x;
    } else {
        TaylorVariable z(x.argument_size());
        Float xv=x.value();
        z.value()=(abs(xv));
        z.error()=(abs(xr.u)-abs(xv));
        return z;
    }

}

Interval _sum(const TaylorVariable& x) {
    typedef TaylorVariable::const_iterator const_iterator;
    Interval r=0;
    for(const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        r+=xiter->second;
    }
    return r;
}

/*
TaylorVariable add(const TaylorVariable& x, const TaylorVariable& y) {
    TaylorVariable r=x;
    r+=y;
    return r;
}

TaylorVariable mul(const TaylorVariable& x, const TaylorVariable& y) {
    ARIADNE_ASSERT(x.argument_size()==y.argument_size());
    typedef TaylorVariable::const_iterator const_iterator;
    TaylorVariable r(x.argument_size());
    Interval& e=r.error();
    for(const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        for(const_iterator yiter=y.begin(); yiter!=y.end(); ++yiter) {
            acc(e,r[xiter->first+yiter->first],xiter->second,yiter->second);
        }
    }
    e += x.error() * _sum(y) + _sum(x) * y.error() + x.error() * y.error();
    return r;
}
*/

TaylorVariable neg(const TaylorVariable& x) {
    TaylorVariable r(x.domain());
    for(TaylorVariable::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        r.expansion().append(xiter->first,-xiter->second);
    }
    r.error()=x.error();
    return r;
}

Interval
TaylorVariable::range() const {
    Interval r(-this->error(),+this->error());
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->first.degree()==0) {
            r+=iter->second;
        } else {
            r+=abs(iter->second)*Interval(-1,1);
        }
    }
    return r;

    /* FIXME: The following code does not work with optimisation turned on using gcc */
    set_rounding_mode(upward);
    volatile Float t=this->error();
    volatile Float v=0.0;
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        if(iter->first.degree()==0) {
            v=iter->second;
        } else {
            t+=abs(iter->second);
        }
    }
    set_rounding_mode(to_nearest);
    return v+Interval(-t,t);
}


Interval
TaylorVariable::evaluate(const Vector<Float>& v) const
{
    return this->evaluate(Vector<Interval>(v));
}

Interval
TaylorVariable::evaluate(const Vector<Interval>& v) const
{
    ARIADNE_ASSERT(this->argument_size()==v.size());
    ARIADNE_ASSERT_MSG(subset(v,this->domain()),"x="<<v<<", d="<<this->domain());
    Interval r=this->error();
    for(const_iterator iter=this->begin(); iter!=this->end(); ++iter) {
        Interval t(iter->second);
        for(uint j=0; j!=iter->first.size(); ++j) {
            t*=pow(v[j],iter->first[j]);
        }
        r+=t;
    }
    return r;
}

template<class X> class Series;
typedef Series<Interval>(*series_function_pointer)(uint,const Interval&);

struct TaylorSeries {
    typedef Series<Interval>(*series_function_pointer)(uint,const Interval&);
    TaylorSeries(uint d) : expansion(d+1), error(0) { }
    TaylorSeries(uint degree, series_function_pointer function,
                 const Float& centre, const Interval& domain);
    uint degree() const { return expansion.size()-1; }
    Float& operator[](uint i) { return expansion[i]; }
    array<Float> expansion;
    Interval error;
    void sweep(Float e) {
        for(uint i=0; i<=degree(); ++i) {
            if(abs(expansion[i])<=e) {
                error+=expansion[i]*Interval(-1,1);
                expansion[i]=0; } } }
};


TaylorVariable
_compose(const TaylorSeries& ts, const TaylorVariable& tv, Float eps)
{
    //std::cerr<<"_compose(TaylorSeries,TaylorVariable,Error)\n";
    //std::cerr<<"\n  ts="<<ts<<"\n  tv="<<tv<<"\n";
    Float& vref=const_cast<Float&>(tv.value());
    Float vtmp=vref;
    vref=0.0;
    TaylorVariable r(tv);
    r+=ts.expansion[ts.expansion.size()-1];
    for(uint i=1; i!=ts.expansion.size(); ++i) {
        //std::cerr<<"    r="<<r<<std::endl;
        r=r*tv;
        r+=ts.expansion[ts.expansion.size()-i-1];
        r.sweep(eps);
    }
    //std::cerr<<"    r="<<r<<std::endl;
    r+=ts.error;
    //std::cerr<<"    r="<<r<<std::endl;
    vref=vtmp;
    return r;
}

TaylorVariable
compose(const TaylorSeries& ts, const TaylorVariable& tv)
{
    return _compose(ts,tv,TaylorVariable::ec);
}


// Compose using the Taylor formula directly. The final term is the Taylor series computed
// over the range of the series. This method tends to suffer from blow-up of the
// truncation error
TaylorVariable
_compose1(const series_function_pointer& fn, const TaylorVariable& tv, Float eps)
{
    static const uint DEGREE=18;
    static const Float TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    Float c=tv.value();
    Interval r=tv.range();
    Series<Interval> centre_series=fn(d,Interval(c));
    Series<Interval> range_series=fn(d,r);

    Float truncation_error_estimate=mag(range_series[d])*pow(mag(r-c),d);
    if(truncation_error_estimate>TRUNCATION_ERROR) {
        std::cerr<<"Warning: Truncation error estimate "<<truncation_error_estimate
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n";
    }

    TaylorVariable x=tv-c;
    TaylorVariable res(tv.domain());
    res+=range_series[d];
    for(uint i=0; i!=d; ++i) {
        //std::cerr<<"i="<<i<<" r="<<res<<"\n";
        res=centre_series[d-i-1]+x*res;
        res.sweep(eps);
    }
    //std::cerr<<"i="<<d<<" r="<<res<<"\n";
    return res;
}

// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. The radius of convergence of this method is still quite low,
// typically only half of the radius of convergence of the power series itself
TaylorVariable
_compose2(const series_function_pointer& fn, const TaylorVariable& tv, Float eps)
{
    static const uint DEGREE=20;
    static const Float TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    Float c=tv.value();
    Interval r=tv.range();
    Series<Interval> centre_series=fn(d,Interval(c));
    Series<Interval> range_series=fn(d,r);

    //std::cerr<<"c="<<c<<" r="<<r<<" r-c="<<r-c<<" e="<<mag(r-c)<<"\n";
    //std::cerr<<"cs[d]="<<centre_series[d]<<" rs[d]="<<range_series[d]<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";
    Float truncation_error=mag(range_series[d]-centre_series[d])*pow(mag(r-c),d);
    //std::cerr<<"te="<<truncation_error<<"\n";
    if(truncation_error>TRUNCATION_ERROR) {
        std::cerr<<"Warning: Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n";
    }

    TaylorVariable x=tv-c;
    TaylorVariable res(tv.domain());
    res+=centre_series[d];
    for(uint i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        res.sweep(eps);
    }
    res+=truncation_error*Interval(-1,1);
    return res;
}


// Compose using the Taylor formula with a constant truncation error. This method
// is usually better than _compose1 since there is no blow-up of the trunction
// error. This method is better than _compose2 since the truncation error is
// assumed at the ends of the intervals
TaylorVariable
_compose3(const series_function_pointer& fn, const TaylorVariable& tv, Float eps)
{
    static const uint DEGREE=20;
    static const Float TRUNCATION_ERROR=1e-8;
    uint d=DEGREE;
    Float c=tv.value();
    Interval r=tv.range();
    Series<Interval> centre_series=fn(d,Interval(c));
    Series<Interval> range_series=fn(d,r);

    //std::cerr<<"c="<<c<<" r="<<r<<" r-c="<<r-c<<" e="<<mag(r-c)<<"\n";
    //std::cerr<<"cs[d]="<<centre_series[d]<<" rs[d]="<<range_series[d]<<"\n";
    //std::cerr<<"cs="<<centre_series<<"\nrs="<<range_series<<"\n";
    Interval se=range_series[d]-centre_series[d];
    Interval e=r-c;
    Interval p=pow(e,d-1); p.l*=-e.l; p.u*=e.u;
    //std::cerr<<"se="<<se<<" e="<<e<<" p="<<p<<std::endl;
    // FIXME: Here we assume the dth derivative of f is monotone increasing
    Float truncation_error=max(se.l*p.l,se.u*p.u);
    //std::cerr<<"te="<<truncation_error<<"\n";
    if(truncation_error>TRUNCATION_ERROR) {
        std::cerr<<"Warning: Truncation error estimate "<<truncation_error
                 <<" is greater than maximum allowable truncation error "<<TRUNCATION_ERROR<<"\n";
    }

    TaylorVariable x=tv;
    TaylorVariable res(tv.domain());
    res+=centre_series[d];
    for(uint i=0; i!=d; ++i) {
        res=centre_series[d-i-1]+x*res;
        res.sweep(eps);
    }
    res+=truncation_error*Interval(-1,1);
    return res;
}


TaylorVariable
_compose(const series_function_pointer& fn, const TaylorVariable& tv, Float eps) {
    //std::cerr<<"_compose(SeriesFunction,TaylorVariable,Error)\n";
    return _compose3(fn,tv,eps);
}


TaylorVariable sqr(const TaylorVariable& x) {
    TaylorVariable r=x*x;
    return r;
}

TaylorVariable pow(const TaylorVariable& x, int n) {
    TaylorVariable r(x.domain()); r+=1;
    TaylorVariable p(x);
    while(n) {
        if(n%2) { r=r*p; }
        p=sqr(p);
        n/=2;
    }
    return r;
}

TaylorVariable sqrt(const TaylorVariable& x) {
    //std::cerr<<"rec(TaylorVariable)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], scale by constant a such that rl/a=1-d; ru/a=1+d
    Interval r=x.range();
    assert(r.l>0);
    Float a=(r.l+r.u)/2;
    set_rounding_upward();
    Float eps=(r.u-r.l)/(r.u+r.l);
    set_rounding_to_nearest();
    assert(eps<1);
    uint d=uint(log((1-eps)*x._sweep_threshold)/log(eps)+1);
    //std::cerr<<"x="<<x<<std::endl;
    //std::cerr<<"x/a="<<x/a<<" a="<<a<<std::endl;
    TaylorVariable y=(x/a)-1.0;
    //std::cerr<<"y="<<y<<std::endl;
    TaylorVariable z(x.domain());
    Series<Interval> sqrt_series=Series<Interval>::sqrt(d,Interval(1));
    //std::cerr<<"sqrt_series="<<sqrt_series<<std::endl;
    //std::cerr<<"y="<<y<<std::endl;
    z+=sqrt_series[d-1];
    for(uint i=0; i!=d; ++i) {
        z=sqrt_series[d-i-1] + z * y;
        z.sweep(); z.truncate();
        //std::cerr<<"z="<<z<<std::endl;
    }
    Float trunc_err=pow(eps,d)/(1-eps)*mag(sqrt_series[d]);
    //std::cerr<<"te="<<trunc_err<<" te*[-1,+1]="<<trunc_err*Interval(-1,1)<<std::endl;
    z.error()+=trunc_err;
    //std::cerr<<"z="<<z<<std::endl;
    Interval sqrta=sqrt(Interval(a));
    //std::cerr<<"sqrt(a)="<<sqrta<<std::endl;
    z*=sqrt(Interval(a));
    //std::cerr<<"z="<<z<<std::endl;
    return z;
}

TaylorVariable rec(const TaylorVariable& x) {
    //std::cerr<<"rec(TaylorVariable)\n";
    // Use a special routine to minimise errors
    // Given range [rl,ru], scale by constant a such that rl/a=1-d; ru/a=1+d
    Interval r=x.range();
    if(r.upper()>=0 && r.lower()<=0) {
        ARIADNE_THROW(DivideByZeroException,"rec(TaylorVariable x)","x="<<x);
    }
    Float a=(r.l+r.u)/2;
    set_rounding_upward();
    Float eps=abs((r.u-r.l)/(r.u+r.l));
    set_rounding_to_nearest();
    assert(eps<1);
    uint d=uint(log((1-eps)*x._sweep_threshold)/log(eps))+1;
    //std::cerr<<"  x="<<x<<"\n";
    //std::cerr<<"  r="<<r<<"\n";
    //std::cerr<<"  a="<<a<<"\n";
    TaylorVariable y=1-(x/a);
    //std::cerr<<"  y="<<y<<"\n";
    TaylorVariable z(x.domain());
    z+=Float(d%2?-1:+1);
    for(uint i=0; i!=d; ++i) {
        z=1.0 + z * y;
        z.sweep(); z.truncate();
        //std::cerr<<"z.sweep_threshold="<<z.sweep_threshold()<<"\n";
    }
    //std::cerr<<"  z="<<z<<"\n";
    Float te=pow(eps,d)/(1-eps);
    //std::cerr<<"  te="<<te<<"\n";
    set_rounding_upward();
    Float nze=te+z.error();
    //std::cerr<<"  nze="<<nze<<"\n";

    set_rounding_to_nearest();
    z.set_error(nze);
    //std::cerr<<"  z="<<z<<"\n";
    z/=a;
    //std::cerr<<"  z="<<z<<"\n";
    return z;
}

TaylorVariable log(const TaylorVariable& x) {
    // Use a special routine to minimise errors
    // Given range [rl,ru], scale by constant a such that rl/a=1-d; ru/a=1+d
    Interval r=x.range();
    assert(r.l>0);
    Float a=(r.l+r.u)/2;
    set_rounding_upward();
    Float eps=(r.u-r.l)/(r.u+r.l);
    set_rounding_to_nearest();
    assert(eps<1);
    uint d=uint(log((1-eps)*x._sweep_threshold)/log(eps)+1);
    TaylorVariable y=x/a-1;
    TaylorVariable z(x.domain());
    z+=Float(d%2?-1:+1)/d;
    for(uint i=1; i!=d; ++i) {
        z=Float((d-i)%2?+1:-1)/(d-i) + z * y;
        z.sweep(); z.truncate();
    }
    z=z*y;
    z.sweep();
    Float trunc_err=pow(eps,d)/(1-eps)/d;
    return z+log(Interval(a))+trunc_err*Interval(-1,1);
}

TaylorVariable exp(const TaylorVariable& x) {
    static const uint DEG=18;
    return _compose(&Series<Interval>::exp,x,x._sweep_threshold);
    return compose(TaylorSeries(DEG,&Series<Interval>::exp,
                                x.value(),x.range()),x);
}

TaylorVariable sin(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::sin,
                                x.value(),x.range()),x);
}

TaylorVariable cos(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::cos,
                                x.value(),x.range()),x);
}

TaylorVariable tan(const TaylorVariable& x) {
    return sin(x)*rec(cos(x));
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::tan,
                                x.value(),x.range()),x);
}

TaylorVariable asin(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::asin,
                                x.value(),x.range()),x);
}

TaylorVariable acos(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::acos,
                                x.value(),x.range()),x);
}

TaylorVariable atan(const TaylorVariable& x) {
    static const uint DEG=18;
    return compose(TaylorSeries(DEG,&Series<Interval>::atan,
                                x.value(),x.range()),x);
}


inline int pow2(uint k) { return 1<<k; }
inline int powm1(uint k) { return (k%2) ? -1 : +1; }



Interval
evaluate(const TaylorVariable& tv, const Vector<Interval>& x)
{
    return tv.evaluate(x);
}

Vector<Interval>
evaluate(const Vector<TaylorVariable>& tv, const Vector<Interval>& x)
{
    Vector<Interval> r(tv.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=evaluate(tv[i],x);
    }
    return r;
}


TaylorVariable antiderivative(const TaylorVariable& x, uint k) {
    ARIADNE_ASSERT(k<x.argument_size());
    // General case with error analysis
    Interval dk=x.domain()[k];
    TaylorVariable r(x.domain());
    const Float& xe=x.error();
    Float& re=r.error();
    re=xe;
    set_rounding_upward();
    volatile Float dkru=(dk.u-dk.l)/2;
    volatile Float dkrl=-dk.u; dkrl+=dk.l; dkrl/=(-2);
    volatile Float u,ml;
    Float te=0; // Twice the maximum accumulated error
    for(TaylorVariable::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        const uint c=xiter->first[k]+1;
        const Float& xv=xiter->second;
        volatile Float mxv=-xv;
        if(xv>=0) {
            u=xv*dkru/c;
            ml=mxv*dkrl/c;
        } else {
            u=xv*dkru/c;
            ml=mxv*dkrl/c;
        }
        te+=(u+ml);
    }
    re+=te/2;

    set_rounding_to_nearest();
    Float dkr=(dk.u-dk.l)/2;
    MultiIndex rindex(r.argument_size());
    r[rindex]=0.0;
    for(TaylorVariable::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        const uint c=xiter->first[k]+1;
        rindex=xiter->first;
        rindex[k]=c;
        r[rindex]=xiter->second*dkr/c;
    }

    return r;
}

Vector<TaylorVariable> antiderivative(const Vector<TaylorVariable>& x, uint k) {
    ARIADNE_ASSERT(k<x.argument_size());
    Vector<TaylorVariable> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=antiderivative(x[i],k);
    }
    return r;
}


pair<TaylorVariable,TaylorVariable>
split(const TaylorVariable& tv, uint j)
{
    ARIADNE_ASSERT(j<tv.argument_size());
    uint as=tv.argument_size();

    MultiIndex index(as);
    Float value;

    Vector<Interval> d1=tv.domain();
    d1[j]=Interval(d1[j].lower(),d1[j].midpoint());
    TaylorVariable r1(d1);
    for(TaylorVariable::const_iterator iter=tv.begin();
        iter!=tv.end(); ++iter)
    {
        index=iter->first;
        value=iter->second;
        uint k=index[j];
        for(uint l=0; l<=k; ++l) {
            index[j]=l;
            r1[index] += bin(k,l) * value / pow2(k);
        }
    }

    Vector<Interval> d2=tv.domain();
    d2[j]=Interval(d2[j].midpoint(),d2[j].upper());
    TaylorVariable r2(d2);
    for(TaylorVariable::const_iterator iter=tv.begin();
        iter!=tv.end(); ++iter)
    {
        index=iter->first;
        value=iter->second;
        uint k=index[j];
        for(uint l=0; l<=k; ++l) {
            index[j]=l;
            // Need brackets in expression below to avoid converting negative signed
            // integer to unsigned
            r2[index] += powm1(k-l) * (bin(k,l) * value / pow2(k));
        }
    }
    // FIXME: Add roundoff errors when computing new expansions

    r1.error()=tv.error();
    r2.error()=tv.error();

    return make_pair(r1,r2);
}


pair< Vector<TaylorVariable>, Vector<TaylorVariable> >
split(const Vector<TaylorVariable>& tv, uint j)
{
    ARIADNE_ASSERT(j<tv.argument_size());
    Vector<TaylorVariable> r1(tv.size());
    Vector<TaylorVariable> r2(tv.size());
    for(uint i=0; i!=tv.size(); ++i) {
        make_lpair(r1[i],r2[i])=split(tv[i],j);
    }
    return make_pair(r1,r2);
}

TaylorVariable
unscale(const TaylorVariable& tv, const Interval& ivl)
{
    // Scale tv so that the interval ivl maps into [-1,1]
    // The result is given by  (tv-c)*s where c is the centre
    // and s the reciprocal of the radius of ivl
    const Float& l=ivl.l;
    const Float& u=ivl.u;

    TaylorVariable r=tv;
    Interval c=Interval(l/2)+Interval(u/2);
    Interval s=2/(Interval(u)-Interval(l));
    r-=c;
    r*=s;

    return r;
}

Vector<TaylorVariable>
unscale(const Vector<TaylorVariable>& tvs, const Vector<Interval>& ivls)
{
    Vector<TaylorVariable> r(tvs.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=unscale(tvs[i],ivls[i]);
    }
    return r;
}


TaylorVariable
scale(const TaylorVariable& tv, const Interval& ivl)
{
    // Scale tv so that the interval [-1,1] maps into ivl
    // The result is given by  (tv*s)+c where c is the centre
    // and r the radius of ivl
    const Float& l=ivl.l;
    const Float& u=ivl.u;

    TaylorVariable r=tv;
    Interval c=add_ivl(l/2,u/2);
    Interval s=sub_ivl(u/2,l/2);
    r*=s;
    r+=c;

    return r;
}

Vector<TaylorVariable>
scale(const Vector<TaylorVariable>& tvs, const Vector<Interval>& ivls)
{
    Vector<TaylorVariable> r(tvs.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=scale(tvs[i],ivls[i]);
    }
    return r;
}


Vector<TaylorVariable>
combine(const Vector<TaylorVariable>& tv1, const Vector<TaylorVariable>& tv2)
{
    const uint rs1=tv1.result_size();
    const uint rs2=tv2.result_size();
    const uint as1=tv1.argument_size();
    const uint as2=tv2.argument_size();
    Vector<Interval> d=join(tv1.domain(),tv2.domain());
    Vector< Expansion<Float> > f(rs1+rs2);
    for(uint i=0; i!=rs1; ++i) { f[i]=tv1[i].expansion().embed(as1+as2,0u); }
    for(uint i=0; i!=rs2; ++i) { f[rs1+i]=tv2[i].expansion().embed(as1+as2,as1); }
    Vector<Float> e(rs1+rs2);
    for(uint i=0; i!=rs1; ++i) { e[i]=tv1[i].error(); }
    for(uint i=0; i!=rs2; ++i) { e[rs1+i]=tv2[i].error(); }
    return Vector<TaylorVariable>(d,f,e);
}

Vector<TaylorVariable>
combine(const Vector<TaylorVariable>& tv1, const TaylorVariable& tv2)
{
    return combine(tv1,Vector<TaylorVariable>(1u,tv2));
}


Vector<TaylorVariable>
embed(const Vector<TaylorVariable>& tv1, const Vector<Interval>& dom2)
{
    const uint rs=tv1.result_size();
    const uint as1=tv1.argument_size();
    const uint as2=dom2.size();
    Vector<Interval> d=join(tv1.domain(),dom2);
    Vector< Expansion<Float> > f(rs);
    for(uint i=0; i!=rs; ++i) { f[i]=tv1[i].expansion().embed(as1+as2,0u); }
    Vector<Float> e(rs);
    for(uint i=0; i!=rs; ++i) { e[i]=tv1[i].error(); }
    return Vector<TaylorVariable>(d,f,e);
}

Vector<TaylorVariable>
embed(const Vector<TaylorVariable>& tv1, const Interval& dom2)
{
    return embed(tv1,Vector<Interval>(1,dom2));
}

TaylorVariable
embed(const TaylorVariable& tv1, const Interval& dom2)
{
    const uint as=tv1.argument_size();
    return TaylorVariable(join(tv1.domain(),dom2),tv1.expansion().embed(as+1,as),tv1.error());
}

Vector<TaylorVariable>
embed(const Vector<Interval>& dom1, const Vector<TaylorVariable>& tv2)
{
    const uint rs=tv2.result_size();
    const uint as1=dom1.size();
    const uint as2=tv2.argument_size();
    Vector<Interval> d=join(dom1,tv2.domain());
    Vector< Expansion<Float> > f(rs);
    for(uint i=0; i!=rs; ++i) { f[i]=tv2[i].expansion().embed(as1+as2,as1); }
    Vector<Float> e(rs);
    for(uint i=0; i!=rs; ++i) { e[i]=tv2[i].error(); }
    return Vector<TaylorVariable>(d,f,e);
}

TaylorVariable
embed(const Vector<Interval>& dom1, const TaylorVariable& tv2)
{
    const uint as1=dom1.size();
    const uint as2=tv2.argument_size();
    return TaylorVariable(join(dom1,tv2.domain()),tv2.expansion().embed(as1+as2,as1),tv2.error());
}

bool
refines(const TaylorVariable& tv1, const TaylorVariable& tv2)
{
    ARIADNE_ASSERT(tv1.argument_size()==tv2.argument_size());
    TaylorVariable d=tv2;
    d.error()=0.0;
    d-=tv1;

    set_rounding_upward();
    Float e=d.error();
    for(TaylorVariable::const_iterator iter=d.begin(); iter!=d.end(); ++iter) {
        e+=abs(iter->second);
    }
    set_rounding_to_nearest();
    return e <= tv2.error();
}

bool
refines(const Vector<TaylorVariable>& tv1, const Vector<TaylorVariable>& tv2)
{
    ARIADNE_ASSERT(tv1.size()==tv2.size());
    for(uint i=0; i!=tv1.size(); ++i) {
        if(!refines(tv1[i],tv2[i])) { return false; }
    }
    return true;
}


void
_compose1(Vector<TaylorVariable>& r,
          const Vector<TaylorVariable>& x,
          const Vector<TaylorVariable>& ys)
{
    //std::cerr<<"x ="<<x<<"\n";
    //std::cerr<<"ys="<<ys<<"\n";
    // Brute-force, unoptimised approach
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=TaylorVariable::constant(ys.domain(),0.0);
        r[i].set_error(x[i].error());
        for(TaylorVariable::const_iterator iter=x[i].begin(); iter!=x[i].end(); ++iter) {
            TaylorVariable t=TaylorVariable::constant(ys.domain(),iter->second);
            for(uint j=0; j!=iter->first.size(); ++j) {
                TaylorVariable p=pow(ys[j],iter->first[j]);
                t=t*p;
            }
            r[i]+=t;
            //std::cerr<<"a="<<iter->first<<" c="<<iter->second<<" t="<<t<<" r="<<r<<"\n";
        }
    }
}


Vector<TaylorVariable>
unchecked_compose(const Vector<TaylorVariable>& x,
                  const Vector<TaylorVariable>& y)
{
   Vector<Interval> bx=x[0].domain();

    Vector<TaylorVariable> ys(y.size());
    for(uint i=0; i!=y.size(); ++i) {
        ys[i]=unscale(y[i],bx[i]);
    }

    Vector<TaylorVariable> r(x.size());
    _compose1(r,x,ys);

    return r;
}

TaylorVariable
unchecked_compose(const TaylorVariable& x,
                  const Vector<TaylorVariable>& y)
{
   return unchecked_compose(Vector<TaylorVariable>(1u,x),y)[0];
}


Vector<TaylorVariable>
compose(const Vector<TaylorVariable>& x,
        const Vector<TaylorVariable>& y)
{
    ARIADNE_ASSERT(x.size()>0);
    ARIADNE_ASSERT(y.size()>0);
    ARIADNE_ASSERT(y.size()==x[0].domain().size());

    for(uint i=1; i!=x.size(); ++i) { ARIADNE_ASSERT(x[i].domain()==x[0].domain()); }
    for(uint i=1; i!=x.size(); ++i) { ARIADNE_ASSERT(y[i].domain()==y[0].domain()); }

    ARIADNE_ASSERT(subset(y.range(),x.domain()));

    Vector<TaylorVariable> r=unchecked_compose(x,y);

    ARIADNE_ASSERT(r.result_size()==x.result_size());
    ARIADNE_ASSERT(r.argument_size()==y.argument_size());

    return r;
}


TaylorVariable
compose(const TaylorVariable& x,
        const Vector<TaylorVariable>& y)
{
    return compose(Vector<TaylorVariable>(1u,x),y)[0];
}

/*

Vector<TaylorVariable>
prod(const Matrix<Float>& A,
     const Vector<TaylorVariable>& x)
{
    ARIADNE_ASSERT(A.column_size()==x.size());

    Vector<TaylorVariable> r=TaylorVariable::zeroes(A.row_size(),x.argument_size());
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint j=0; j!=A.column_size(); ++j) {
            r[i]+=A[i][j]*x[j];
        }
    }
    return r;
}


Vector<TaylorVariable>
prod(const Matrix<Interval>& A,
     const Vector<TaylorVariable>& x)
{
    ARIADNE_ASSERT(A.column_size()==x.size());

    Vector<TaylorVariable> r=TaylorVariable::zeroes(A.row_size(),x.argument_size());
    for(uint i=0; i!=A.row_size(); ++i) {
        for(uint j=0; j!=A.column_size(); ++j) {
            r[i]+=A[i][j]*x[j];
        }
    }
    return r;
}

*/

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
        Vector<Interval> ih_range=join(h.domain(),h.range());
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
        yz[i].set_gradient(i,1.0);
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




Differential<Float> differential(const TaylorVariable& x) {
    Differential<Float> sd(x.argument_size(),x.degree());
    for(TaylorVariable::const_iterator xiter=x.begin(); xiter!=x.end(); ++xiter) {
        sd[xiter->first]=xiter->second;
    }
    return sd;
}





std::string
TaylorVariable::str() const
{
    std::stringstream ss;
    for(TaylorVariable::const_iterator iter=this->_expansion.begin();
        iter!=this->_expansion.end(); ++iter)
    {
        const MultiIndex& j=iter->first;
        const Float& c=iter->second;
        ss << j << ": " << c << "\n";
    }
    ss<<"error: "<<this->_error<<"\n";
    return ss.str();
}


std::ostream&
operator<<(std::ostream& os, const TaylorVariable& tv) {
    //    return os << "TaylorVariable( expansion=" << tv.expansion() << ", error=" << tv.error()
    //              << ", range="<<tv.range()<<" )";
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


