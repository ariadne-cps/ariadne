/***************************************************************************
 *            rational.cc
 *
 *  Copyright 2008-10  Pieter Collins
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



#include "config.h"

#include "utility/macros.h"
#include "numeric/rational.h"

namespace Ariadne {


Rational sqr(const Rational& q) {
    return q*q;
}

Rational pow(const Rational& q, uint n) {
    if(n==0) { return 1; }
    Rational r=1; Rational p=q; uint m=n;
    while(m>=1) { if(m%2) { r*=p; } m/=2; p=p*p; }
    return r;
}

Rational pow(const Rational& q, int n) {
    if(n>=0) { return pow(q,uint(n)); }
    else { return pow(1/q,uint(-n)); }
}


} // namespace Ariadne

