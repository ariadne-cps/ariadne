/***************************************************************************
 *            affine_model.cc
 *
 *  Copyright 2009  Pieter Collins
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

#include "config.h"

#include "numeric.h"
#include "vector.h"
#include "function.h"
#include "taylor_model.h"
#include "affine_model.h"

#include "affine.h"
#include "taylor_function.h"
#include <include/vector.h>

namespace Ariadne {

AffineModel<Interval> operator+(const AffineModel<Interval>& a1, const AffineModel<Interval>& a2) {
    ARIADNE_ASSERT_MSG(a1.argument_size()==a2.argument_size(),"a1="<<a1<<" a2="<<a2);
    Nat n=a1.argument_size();
    AffineModel<Interval> r(n);
    r=a1.value()+a2.value();
    for(uint i=0; i!=n; ++i) {
        r[i]=a1[i]+a2[i];
    }

    set_rounding_upward();

    Float te=0.0;
    for(uint j=0; j!=n; ++j) {
        Float mrjl = (-a1.gradient(j))-a2.gradient(j);
        Float  rju = ( a1.gradient(j))+a2.gradient(j);
        te+=(rju+mrjl);
    }
    Float mrl = (-a1.value())-a2.value();
    Float  ru = ( a1.value())+a2.value();
    te += (ru+mrl);

    Float re=0.0;
    r.set_error( te/2 + (a1.error()+a2.error()) );

    set_rounding_to_nearest();

    return r;
}

AffineModel<Interval> operator+(const Interval& c, const AffineModel<Interval>& a) {
    AffineModel<Interval> r=a;
    Float cm=midpoint(c);
    r.set_value( cm + a.value() );

    set_rounding_upward();

    Float mrl = (-a.value())-cm;
    Float  ru = ( a.value())+cm;
    Float te = (ru+mrl)/2;

    r.set_error( a.error() + max(c.upper()-cm,cm-c.lower()) + te);

    set_rounding_to_nearest();

    return r;
}

AffineModel<Interval> operator+(const AffineModel<Interval>& a, const Interval& c) {
    return c+a;
}

AffineModel<Interval> operator*(const Interval& c, const AffineModel<Interval>& a) {
    Nat n=a.argument_size();
    Float cm=midpoint(c);
    AffineModel<Interval> r(n);
    r=a.value()*cm;
    for(uint i=0; i!=n; ++i) {
        r[i]=a[i]*cm;
    }

    set_rounding_upward();

    Float te=0.0;
    for(uint j=0; j!=n; ++j) {
        Float mca=(-cm)*a.gradient(j);
        Float ca= cm*a.gradient(j);
        te+=(ca+mca);
    }
    Float mca=(-cm)*a.value();
    Float ca= cm*a.value();

    Float re=0.0;
    if(c.lower()!=c.upper()) {
        Float ce=max(c.upper()-cm,cm-c.lower());
        for(uint j=0; j!=n; ++j) {
            re+=abs(a.gradient(j)*ce);
        }
    }

    r.set_error(abs(cm)*a.error() + ((ca+mca) + te)/2 + re);

    set_rounding_to_nearest();

    return r;
}

AffineModel<Interval> operator*(const AffineModel<Interval>& a, const Interval& c) {
    return c*a;
}

IntervalAffineModel affine_model(const IntervalAffine& a) {
    IntervalAffineModel am(a.argument_size());
    am = midpoint(a.value());
    for(uint j=0; j!=a.argument_size(); ++j) {
        am[j] = midpoint(a[j]);
    }
    set_rounding_upward();
    Float e = 0.0;
    for(uint j=0; j!=a.argument_size(); ++j) {
        e += max(a.gradient(j).upper()-am.gradient(j),am.gradient(j)-a.gradient(j).lower());
    }
    e += max(a.value().upper()-am.value(),am.value()-a.value().lower());
    am.set_error(e);
    set_rounding_to_nearest();
    return am;
}

AffineModel<Interval> affine_model(const TaylorModel<Interval>& taylor_model) {
    AffineModel<Interval> affine_model(taylor_model.argument_size());

    rounding_mode_t rnd=get_rounding_mode();
    set_rounding_upward();
    for(TaylorModel<Interval>::const_iterator iter=taylor_model.begin(); iter!=taylor_model.end(); ++iter) {
        if(iter->key().degree()>=2) {
            affine_model.set_error(abs(iter->data()+affine_model.error()));
        } else if(iter->key().degree()==1) {
            for(uint i=0; i!=taylor_model.argument_size(); ++i) {
                if(iter->key()[i]==1) {
                    affine_model.set_gradient(i,iter->data());
                    break;
                }
            }
        } else {
            affine_model.set_value(iter->data());
        }
    }
    affine_model.set_error(taylor_model.error()+affine_model.error());
    set_rounding_mode(rnd);

    return affine_model;
}

Vector< AffineModel<Interval> > affine_models(const Vector< TaylorModel<Interval> >& models)
{
    Vector< AffineModel<Interval> > result(models.size(),AffineModel<Interval>(models.size()==0?0u:models[0].argument_size()));
    for(uint i=0; i!=result.size(); ++i) { result[i]=affine_model(models[i]); }
    return result;
}

AffineModel<Interval> affine_model(const IntervalVector& domain, const ValidatedScalarFunction& function)
{
    return affine_model(ScalarTaylorFunction(domain,function,AffineSweeper()).model());
}

Vector< AffineModel<Interval> > affine_models(const IntervalVector& domain, const ValidatedVectorFunction& function)
{
    return affine_models(VectorTaylorFunction(domain,function,AffineSweeper()).models());
}

std::ostream& operator<<(std::ostream& os, const AffineModel<Interval>& f)
{
    os << f.value();
    for(uint j=0; j!=f.argument_size(); ++j) {
        if(f.gradient(j)>0.0) { os << "+"; }
        if(f.gradient(j)!=0.0) { os << f.gradient(j) << "*x" << j; }
    }
    return os << "+/-" << f.error();

}

} //namespace Ariadne


