/***************************************************************************
 *            vector.cc
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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

#include "numeric/numeric.h"
#include "config.h"

#include "utility/macros.h"
#include "algebra/vector.h"

namespace Ariadne {

Vector<ExactFloat>const& make_exact(const Vector<ApproximateFloat>& av) {
    return reinterpret_cast<Vector<ExactFloat>const&>(av);
}

Vector<ValidatedFloat> make_bounds(const Vector<ErrorFloat>& ev) {
    Vector<ValidatedFloat> r(ev.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=make_bounds(ev[i]);
    }
    return r;
}

UpperFloat sup_error(const Vector<ValidatedFloat>& x) {
    UpperFloat e(0);
    for(uint i=0; i!=x.size(); ++i) {
        e=max(e,x[i].error());
    }
    return e;
}

Vector<ExactFloat> midpoint(const Vector<ValidatedFloat>& x) {
    Vector<ExactFloat> r(x.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=midpoint(x[i]);
    }
    return r;
}

bool models(const Vector<ValidatedFloat>& x1, const Vector<ExactFloat>& x2) {
    assert(x1.size()==x2.size());
    for(uint i=0; i!=x1.size(); ++i) {
        if(!models(x1[i],x2[i])) { return false; }
    }
    return true;
}

bool consistent(const Vector<ValidatedFloat>& x1, const Vector<ValidatedFloat>& x2) {
    assert(x1.size()==x2.size());
    for(uint i=0; i!=x1.size(); ++i) {
        if(!consistent(x1[i],x2[i])) { return false; }
    }
    return true;
}

bool inconsistent(const Vector<ValidatedFloat>& x1, const Vector<ValidatedFloat>& x2) {
    return !consistent(x1,x2);
}

bool refines(const Vector<ValidatedFloat>& x1, const Vector<ValidatedFloat>& x2) {
    assert(x1.size()==x2.size());
    for(uint i=0; i!=x1.size(); ++i) {
        if(!refines(x1[i],x2[i])) { return false; }
    }
    return true;
}

Vector<ValidatedFloat> refinement(const Vector<ValidatedFloat>& x1, const Vector<ValidatedFloat>& x2) {
    assert(x1.size()==x2.size());
    Vector<ValidatedFloat> r(x1.size());
    for(uint i=0; i!=r.size(); ++i) {
        r[i]=refinement(x1[i],x2[i]);
    }
    return r;
}


} // namespace Ariadne
