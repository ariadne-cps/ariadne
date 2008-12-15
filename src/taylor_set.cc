/***************************************************************************
 *            taylor_set.cc
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
 
#include "macros.h"
#include "exceptions.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "multi_index.h"
#include "sparse_differential.h"
#include "function.h"
#include "approximate_taylor_model.h"

#include "taylor_variable.h"
#include "taylor_set.h"

#include "zonotope.h"
#include "polytope.h"
#include "graphics.h"

namespace Ariadne {

TaylorSet::TaylorSet(uint d) 
    : _variables(d)
{
}


TaylorSet::TaylorSet(const ApproximateTaylorModel& atm) 
    : _variables(atm.result_size())
{
    for(size_t i=0; i!=atm.result_size(); ++i) {
        this->_variables[i]=TaylorVariable(atm.expansion()[i],Interval(0));
    }
}







tribool 
TaylorSet::disjoint(const Vector<Interval>& bx) const
{
    return Ariadne::disjoint(this->bounding_box(),bx) || indeterminate;
}


tribool 
TaylorSet::overlaps(const Vector<Interval>& bx) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


tribool
TaylorSet::subset(const Vector<Interval>& bx) const
{
    ARIADNE_NOT_IMPLEMENTED;
}


Vector<Interval>
TaylorSet::bounding_box() const
{
    Vector<Interval> r(this->dimension());
    for(uint i=0; i!=this->dimension(); ++i) {
        r[i]=(*this)[i].range();
    }
    return r;
}


Zonotope 
zonotope(const TaylorSet& ts)
{
    uint d=ts.dimension();
    uint ng=ts.number_of_generators();
    Vector<Float> c(d);
    Matrix<Float> G(d,ng);
    Vector<Float> e(d);
    for(uint i=0; i!=d; ++i) {
        if(ts[i].error()==0) {
            c[i]=ts[i].expansion().value();
            e[i]=0;
        } else {
            Interval ce=ts[i].expansion().value()+ts[i].error();
            c[i]=midpoint(ce);
            e[i]=radius(ce);
        }
        for(uint j=0; j!=ng; ++j) {
            G[i][j]=ts[i].expansion().gradient(j);
        }
    }

    for(uint i=0; i!=d; ++i) {
        const SparseDifferential<Float>& sd=ts[i].expansion();
        for(SparseDifferential<Float>::const_iterator iter=sd.begin(); iter!=sd.end(); ++iter) {
            if(iter->first.degree()>=2) {
                e[i]=add_up(e[i],abs(iter->second));
            }
        }
    }
    return Zonotope(c,G,e);
}


std::ostream&
TaylorSet::write(std::ostream& os) const
{
    os << "TaylorSet(\n";
    os << "  dimension=" << this->dimension() << ",\n" << std::flush;
    os << "  variables=" << this->_variables << ",\n" << std::flush;
    os << ")\n";
    return os;
}



void draw(Figure& fig, const TaylorSet& ts) {
    draw(fig,zonotope(ts));
}

} // namespace Ariadne
 
