/***************************************************************************
 *            propagator.cc
 *
 *  Copyright  2009  Pieter Collins
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

#include "numeric.h"
#include "propagator.h"
#include "function.h"

namespace Ariadne {

template class Propagator<Interval>;

void reduce(Vector<Interval>& x, Propagator<Interval>& p) {
    p.assign(x);
    p.propagate(x);
}

void reduce(const ScalarFunction& f, Vector<Interval>& x) {
    Propagator<Interval> r=f.propagator();
    r.assign(x);
    r.propagate(x);
}


} // namespace Ariadne
