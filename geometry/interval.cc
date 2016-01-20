/***************************************************************************
 *            interval.cc
 *
 *  Copyright 2013-14  Pieter Collins
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

#include "interval.h"

namespace Ariadne {

Interval<Float64Value> widen_domain(Interval<Float64UpperBound> const& ivl) {
    volatile float min=std::numeric_limits<float>::min();
    volatile double l=-ivl.lower().get_d();
    volatile double u=ivl.upper().get_d();
    volatile double neg_l=-l;
    volatile float neg_rl=neg_l;
    volatile float ru=u;
    if(l==u) { neg_rl+=min; ru+=min; }
    if(neg_rl<neg_l) { neg_rl+=min; }
    if(ru<u) { ru+=min; }
    return Interval<Float64Value>(-neg_rl,ru);
}

InputStream&
operator>>(InputStream& is, Interval<Float64Value>& ivl)
{
    Float64 l,u;
    char cl,cm,cr;
    is >> cl >> l >> cm >> u >> cr;
    ARIADNE_ASSERT(is);
    ARIADNE_ASSERT(cl=='[' || cl=='(');
    ARIADNE_ASSERT(cm==':' || cm==',' || cm==';');
    ARIADNE_ASSERT(cr==']' || cr==')');
    ivl.set(Float64Value(l),Float64Value(u));
    return is;
}

} // namespace Ariadne
