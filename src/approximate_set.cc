/***************************************************************************
 *            approximate_set.cc
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
 

#include <iosfwd>

#include "tribool.h"
#include "approximate_set.h"

namespace Ariadne {

Grid g(2);
HybridGrid hg();

void test() {
    InnerApproximation inna(g);
    LowerApproximation lowa(g);
    OuterApproximation outa(g);
    MetricApproximation loca(g);

    HybridInnerApproximation hinna();
    HybridLowerApproximation hlowa();
    HybridOuterApproximation houta();
    HybridMetricApproximation hloca();


}

}
