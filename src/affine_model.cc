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

ValidatedAffineModel operator+(const ValidatedAffineModel& a1, const ValidatedAffineModel& a2) {
    ARIADNE_ASSERT_MSG(a1.argument_size()==a2.argument_size(),"a1="<<a1<<" a2="<<a2);
    Nat n=a1.argument_size();
    ValidatedAffineModel r(n);
    r=a1.value()+a2.value();
    for(uint i=0; i!=n; ++i) {
        r[i]=a1[i]+a2[i];
    }

    set_rounding_upward();

    RawNumberType te=0.0;
    for(uint j=0; j!=n; ++j) {
        RawNumberType mrjl = (-a1.gradient(j))-a2.gradient(j);
        RawNumberType  rju = ( a1.gradient(j))+a2.gradient(j);
        te+=(rju+mrjl);
    }
    RawNumberType mrl = (-a1.value())-a2.value();
    RawNumberType  ru = ( a1.value())+a2.value();
    te += (ru+mrl);

    RawNumberType re=0.0;
    r.set_error( te/2 + (a1.error()+a2.error()) );

    set_rounding_to_nearest();

    return r;
}

ValidatedAffineModel operator+(const ValidatedNumberType& c, const ValidatedAffineModel& a) {
    ValidatedAffineModel r=a;
    RawNumberType cm=midpoint(c);
    r.set_value( cm + a.value() );

    set_rounding_upward();

    RawNumberType mrl = (-a.value())-cm;
    RawNumberType  ru = ( a.value())+cm;
    RawNumberType te = (ru+mrl)/2;

    r.set_error( a.error() + max(c.upper()-cm,cm-c.lower()) + te);

    set_rounding_to_nearest();

    return r;
}

ValidatedAffineModel operator+(const ValidatedAffineModel& a, const ValidatedNumberType& c) {
    return c+a;
}

ValidatedAffineModel operator*(const ValidatedNumberType& c, const ValidatedAffineModel& a) {
    Nat n=a.argument_size();
    RawNumberType cm=midpoint(c);
    ValidatedAffineModel r(n);
    r=a.value()*cm;
    for(uint i=0; i!=n; ++i) {
        r[i]=a[i]*cm;
    }

    set_rounding_upward();

    RawNumberType te=0.0;
    for(uint j=0; j!=n; ++j) {
        RawNumberType mca=(-cm)*a.gradient(j);
        RawNumberType ca= cm*a.gradient(j);
        te+=(ca+mca);
    }
    RawNumberType mca=(-cm)*a.value();
    RawNumberType ca= cm*a.value();

    RawNumberType re=0.0;
    if(c.lower()!=c.upper()) {
        RawNumberType ce=max(c.upper()-cm,cm-c.lower());
        for(uint j=0; j!=n; ++j) {
            re+=abs(a.gradient(j)*ce);
        }
    }

    r.set_error(abs(cm)*a.error() + ((ca+mca) + te)/2 + re);

    set_rounding_to_nearest();

    return r;
}

ValidatedAffineModel operator*(const ValidatedAffineModel& a, const ValidatedNumberType& c) {
    return c*a;
}

ValidatedAffineModel affine_model(const ValidatedAffine& a) {
    ValidatedAffineModel am(a.argument_size());
    am = midpoint(a.value());
    for(uint j=0; j!=a.argument_size(); ++j) {
        am[j] = midpoint(a[j]);
    }
    set_rounding_upward();
    RawNumberType e = 0.0;
    for(uint j=0; j!=a.argument_size(); ++j) {
        e += max(a.gradient(j).upper()-am.gradient(j),am.gradient(j)-a.gradient(j).lower());
    }
    e += max(a.value().upper()-am.value(),am.value()-a.value().lower());
    am.set_error(e);
    set_rounding_to_nearest();
    return am;
}

ValidatedAffineModel affine_model(const ValidatedTaylorModel& taylor_model) {
    ValidatedAffineModel affine_model(taylor_model.argument_size());

    rounding_mode_t rnd=get_rounding_mode();
    set_rounding_upward();
    for(ValidatedTaylorModel::const_iterator iter=taylor_model.begin(); iter!=taylor_model.end(); ++iter) {
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

Vector< ValidatedAffineModel > affine_models(const Vector< ValidatedTaylorModel >& models)
{
    Vector< ValidatedAffineModel > result(models.size(),ValidatedAffineModel(models.size()==0?0u:models[0].argument_size()));
    for(uint i=0; i!=result.size(); ++i) { result[i]=affine_model(models[i]); }
    return result;
}

ValidatedAffineModel affine_model(const Box& domain, const ValidatedScalarFunction& function)
{
    return affine_model(ScalarTaylorFunction(domain,function,AffineSweeper()).model());
}

Vector< ValidatedAffineModel > affine_models(const Box& domain, const ValidatedVectorFunction& function)
{
    return affine_models(VectorTaylorFunction(domain,function,AffineSweeper()).models());
}

std::ostream& operator<<(std::ostream& os, const ValidatedAffineModel& f)
{
    os << f.value();
    for(uint j=0; j!=f.argument_size(); ++j) {
        if(f.gradient(j)>0.0) { os << "+"; }
        if(f.gradient(j)!=0.0) { os << f.gradient(j) << "*x" << j; }
    }
    return os << "+/-" << f.error();

}

} //namespace Ariadne


