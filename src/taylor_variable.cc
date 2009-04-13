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
#include "expression_interface.h"
#include "exceptions.h"

namespace Ariadne {


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

TaylorVariable::TaylorVariable(const DomainType& d, const ExpressionInterface& f)
    : _domain(d), _model(f.argument_size())
{
    ARIADNE_ASSERT(d.size()==f.argument_size());
    Vector<TaylorModel> x=TaylorModel::scalings(d);
    this->_model=f.evaluate(x);
}


TaylorVariable TaylorVariable::constant(const Vector<Interval>& d, const Float& c)
{
    return TaylorVariable(d,TaylorModel::constant(d.size(),c));
}

TaylorVariable TaylorVariable::constant(const Vector<Interval>& d, const Interval& c)
{
    return TaylorVariable(d,TaylorModel::constant(d.size(),c));
}

TaylorVariable TaylorVariable::variable(const Vector<Interval>& d, uint j)
{
    return TaylorVariable(d,TaylorModel::scaling(d.size(),j,d[j]));
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
TaylorVariable::evaluate(const Vector<Float>& x) const
{
    return this->evaluate(Vector<Interval>(x));
}

Interval
TaylorVariable::evaluate(const Vector<Interval>& x) const
{
    Vector<Interval> sx=Ariadne::evaluate(TaylorModel::scalings(this->_domain),x);
    return Ariadne::evaluate(this->_model,sx);
}



TaylorVariable restrict(const TaylorVariable& tv, const Vector<Interval>& d) {
    ARIADNE_ASSERT(subset(d,tv.domain()));
    if(d==tv.domain()) { return tv; }
    Vector<TaylorModel> s=TaylorModel::rescalings(tv.domain(),d);
    return TaylorVariable(d,compose(tv._model,s));
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

Matrix<Interval>
jacobian(const Vector<TaylorVariable>& tv, const Vector<Interval>& x);

Vector<TaylorVariable>
implicit(const Vector<TaylorVariable>& f)
{
    ARIADNE_ASSERT(f.size()>0);
    ARIADNE_ASSERT(f.size()<=f[0].argument_size());
    // Solve the equation f(x,h(x))=0
    // Use D1f + D2f Dh = 0, so Dh=-D2f^-1 D1f
    uint rs=f.size();
    uint fas=f[0].argument_size();
    uint has=fas-rs;

    Vector<Interval> h_domain=project(f[0].domain(),range(0u,has));
    Vector<Interval> h_range=project(f[0].domain(),range(has,fas));
    Vector<TaylorVariable> id=TaylorVariable::variables(h_domain);
    Vector<TaylorVariable> h=TaylorVariable::constants(h_domain,h_range);
    //std::cerr<<"\nid="<<id<<"\nh0="<<h0<<"\n";

    for(uint k=0; k!=10; ++k) {
        Vector<Interval> ih_range=join(h_domain,h_range);
        Matrix<Interval> Df=jacobian(f,ih_range);
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
    ARIADNE_ASSERT(h.size()==f.size());
    ARIADNE_ASSERT(h[0].argument_size()+h.size()==f[0].argument_size());
    return h;

}

TaylorVariable implicit(const TaylorVariable& f) {
    Vector<Interval> h_domain=project(f.domain(),range(0u,f.argument_size()-1u));
    Interval h_codomain=f.domain()[f.argument_size()-1u];
    TaylorModel h_model=implicit(f.model());
    ARIADNE_ASSERT(h_model.argument_size()+1==f.model().argument_size());
    TaylorModel hrs_model=h_model.rescale(Interval(-1,+1),h_codomain);
    ARIADNE_ASSERT(hrs_model.argument_size()+1==f.model().argument_size());
    return TaylorVariable(h_domain,hrs_model);
}



std::ostream&
operator<<(std::ostream& os, const TaylorVariable& tv) {
    return os << "TaylorVariable(" << tv.domain() << "," << tv.expansion() << "," << tv.error() << ")";
}





bool
check(const Vector<TaylorVariable>& tv)
{
    for(uint i=0; i!=tv.size(); ++i) {
        if(tv[0].domain()!=tv[i].domain()) { return false; }
    }
    return true;
}

Vector< Expansion<Float> >
expansion(const Vector<TaylorVariable>& x)
{
    Vector< Expansion<Float> > r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].expansion();
    }
    return r;
}

Vector<Float>
error(const Vector<TaylorVariable>& x)
{
    Vector<Float> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].error();
    }
    return r;
}

Vector<Float>
value(const Vector<TaylorVariable>& x)
{
    Vector<Float> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].value();
    }
    return r;
}

Vector<Interval>
ranges(const Vector<TaylorVariable>& x)
{
    Vector<Interval> r(x.size());
    for(uint i=0; i!=x.size(); ++i) {
        r[i]=x[i].range();
    }
    return r;
}


Vector<Interval>
evaluate(const Vector<TaylorVariable>& tv, const Vector<Interval>& x)
{
    Vector<Interval> r(tv.size());
    for(uint i=0; i!=tv.size(); ++i) {
        r[i]=evaluate(tv[i],x);
    }
    return r;
}

Matrix<Interval>
jacobian(const Vector<TaylorVariable>& tv, const Vector<Interval>& x)
{
    ARIADNE_ASSERT(check(tv));
    const Vector<Interval>& dom=tv[0].domain();
    const uint n=dom.size();
    Vector< Differential<Interval> > s(n,n,1u);
    for(uint j=0; j!=n; ++j) {
        Interval dj=dom[j];
        s[j].set_value((x[j]-dj.midpoint())/dj.radius());
        s[j].set_gradient(j,1/dj.radius());
    }
    Vector< Expansion<Float> > p=expansion(tv);
    Vector< Differential<Interval> > d=evaluate(p,s);
    //std::cerr<<"  x="<<x<<"\n  p="<<p<<"\n"<<"  s="<<s<<"\n  p.s="<<d<<"\n  J="<<d.jacobian()<<"\n"<<std::endl;
    return d.jacobian();
}





} //namespace Ariadne


