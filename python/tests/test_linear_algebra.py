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

print(dir(FloatDPBoundsVector))

def exact(x): return Dyadic(ExactDouble(x))

Precision=DoublePrecision
ApproximateScalar=FloatDPApproximation
ValidatedScalar=FloatDPBounds
ExactScalar=FloatDPValue
ApproximateVector=FloatDPApproximationVector
ValidatedVector=FloatDPBoundsVector
ApproximateCovector=FloatDPApproximationCovector
ValidatedCovector=FloatDPBoundsCovector
ApproximateMatrix=FloatDPApproximationMatrix
ValidatedMatrix=FloatDPBoundsMatrix

n=2
d=2.125
pr=Precision()
one=ExactScalar(exact(1.000),pr)
xa=ApproximateScalar(2.125,pr)
xb=ValidatedScalar(exact(2.00),exact(2.25),pr)

def test_vector():
    va=ApproximateVector(2,pr)
    va=ApproximateVector(va)
    va=ApproximateVector([1.125,2.125],pr)
    va=ApproximateVector([1.125,xa],pr)
    vb=ValidatedVector(2,pr)
    vb=ValidatedVector([exact(1.125),exact(2.125)],pr)
    vb=ValidatedVector([{2:3},{3:4}],pr)
    vb=ValidatedVector(vb)

    (+va,-va,va+va,va-va,xa*va,va*xa,va/xa)
    (+vb,-vb,vb+vb,vb-vb,xb*vb,vb*xb,vb/xb)
    
    # Mixed vector operations
    # Disallowed mixed ValidatedVector - ApproximateScalar operations
    (va+vb,va-vb,va*xb,va/xb)
    (vb+va,vb-va)
    
    # NOTE: No GenericScalar - ConcreteVector operations
#    (va*n,va/n)
#    (va*d,va/d)
#    (vb*n,vb/n)
 
 
def test_covector():
    ua=ApproximateCovector([1.125,xa],pr)
    ua=ApproximateCovector([1.125,2.125],pr)
    ub=ValidatedCovector([exact(1.125),xb],pr)
    ub=ValidatedCovector([1,2],pr)

    (+ua,-ua,ua+ua,ua-ua,xa*ua,ua*xa,ua/xa)
    (+ub,-ub,ub+ub,ub-ub,xb*ub,ub*xb,ub/xb)

    (ua+ub,ua-ub,ua*xb,ua/xb)
    (ub+ua,ub-ua)
    
#    (ua*n,ua/n)
#    (ua*d,ua/d)
#    (ub*n,ub/n)


def test_matrix():
    va=ApproximateVector([1,2],pr)
    vb=ValidatedVector([1,2],pr)
    ua=ApproximateCovector([1,2],pr)
    ub=ValidatedCovector([1,2],pr)

    Aa=ApproximateMatrix(2,2,pr)
    Aa=ApproximateMatrix([[xa,1],[1.0,one]],pr)
    Aa=ApproximateMatrix([[2.125,1],[1.0,1]],pr)
    Aa=ApproximateMatrix(Aa)
    Ab=ValidatedMatrix([[{exact(1.875):exact(2.125)},{exact(0.875):exact(1.125)}],[{exact(0.875):exact(1.125)},{exact(0.875):exact(1.125)}]],pr)

    (+Aa,-Aa,Aa+Aa,Aa-Aa,Aa*Aa,xa*Aa,Aa*xa,Aa/xa,Aa*va,ua*Aa)
    (+Ab,-Ab,Ab+Ab,Ab-Ab,Ab*Ab,xb*Ab,Ab*xb,Ab/xb,Ab*vb,ub*Ab)

    # NOTE: No ValidatedMatrix - ApproximateScalar operations
    (Aa+Ab,Aa-Ab,Aa*Ab,Aa*xb,Aa/xb,Aa*vb) #ua*Ab
    (Ab+Aa,Ab-Aa,Ab*Aa,xb*Aa,ub*Aa) #Ab*va

#    (Aa*n,Aa/n)
#    (Aa*d,Aa/d)
#    (Ab*n,Ab/n)

