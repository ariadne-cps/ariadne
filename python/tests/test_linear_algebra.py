#!/usr/bin/python

##############################################################################
#            test_linear_algebra.py
#
#  Copyright 2007  Pieter Collins <Pieter.Collins@cwi.nl>
##############################################################################

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR Aa PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

from ariadne import *

def test_linear_algebra():
    Precision=DoublePrecision
    ApproximateScalar=FloatDPApproximation
    ValidatedScalar=FloatDPBounds
    ExactScalar=FloatDPValue
    ApproximateVector=FloatDPApproximationVector
    ValidatedVector=FloatDPBoundsVector
    ApproximateMatrix=FloatDPApproximationMatrix
    ValidatedMatrix=FloatDPBoundsMatrix

    n=2
    d=2.125
    pr=Precision()
    one=ExactScalar(1.000,pr)
    xa=ApproximateScalar(2.125,pr)
    xb=ValidatedScalar(2.00,2.25,pr)
    va=ApproximateVector(2)
    va=ApproximateVector(va)
    va=ApproximateVector([1.125,2.125])
    va=ApproximateVector([1.125,x])
    vb=BoundsVector(2)
    vb=BoundsVector(va)
    vb=BoundsVector([1.125,2.125])
    vb=BoundsVector([{2:3},{3:4}])
    vb=BoundsVector(vb)
    #cv=Covector([1.125,x])
    #icv=IntervalCovector([1.125,x])
    Aa=ApproximateMatrix(2,2)
    Aa=ApproximateMatrix(Aa)
    Aa=ApproximateMatrix([[x,1],[1.0,one]])
    #Ab=BoundsMatrix([["[1.875,2.125]","[0.875,1.125]"],["[0.875,1.125]","[0.875,1.125]"]])
    Ab=BoundsMatrix([[{1.875:2.125},{0.875:1.125}],[{0.875:1.125},{0.875:1.125}]])

    (-va,-vb)
    (va+va,va+vb,vb+va,vb+vb)
    (va-va,va-vb,vb-va,vb-vb)
    (x*va,x*vb,ix*va,ix*vb)
    (va*x,vb*x,va*ix,vb*ix)
    (va/x,vb/x,va/ix,vb/ix)

    (n*va,n*vb,va*n,vb*n,va/n,vb/n)
    (d*va,d*vb,va*d,vb*d,va/d,vb/d)

    #(-cv,-icv)
    #(cv+cv,cv+icv,icv+cv,icv+icv)
    #(cv-cv,cv-icv,icv-cv,icv-icv)
    #(x*cv,x*icv,ix*cv,ix*icv)
    #(cv*x,icv*x,cv*ix,icv*ix)
    #(cv/x,icv/x,cv/ix,icv/ix)

    #(n*cv,n*icv,cv*n,icv*n,cv/n,icv/n)
    #(d*cv,d*icv,cv*d,icv*d,cv/d,icv/d)

    (-Aa,-Ab)
    (Aa+Aa,Aa+Ab,Ab+Aa,Ab+Ab)
    (Aa-Aa,Aa-Ab,Ab-Aa,Ab-Ab)
    (x*Aa,x*Ab,ix*Aa,ix*Ab)
    (Aa*x,Ab*x,Aa*ix,Ab*ix)
    (Aa/x,Ab/x,Aa/ix,Ab/ix)

    (n*Aa,n*Ab,Aa*n,Ab*n,Aa/n,Ab/n)
    (d*Aa,d*Ab,Aa*d,Ab*d,Aa/d,Ab/d)

    #(cv*va,icv*va,cv*vb,icv*vb)
    (Aa*va,Ab*va,Aa*vb,Ab*vb)
    #(cv*Aa,icv*Aa,cv*Ab,icv*Ab)
    (Aa*Aa,Ab*Aa,Aa*Ab,Ab*Ab)
