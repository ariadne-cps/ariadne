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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

from ariadne import *

n=2
d=2.125
one=Float(1.000)
x=Float(2.125)
ix=Interval(2.00,2.25)
v=Vector(2)
v=Vector(v)
v=Vector([1.125,2.125])
v=Vector([1.125,x])
iv=IntervalVector(2)
iv=IntervalVector(v)
iv=IntervalVector([(2,3),["3",4]])
iv=IntervalVector(iv)
cv=Covector([1.125,x])
icv=IntervalCovector([1.125,x])
A=Matrix(2,2)
A=Matrix(A)
A=Matrix([[x,1],[1.0,one]])
iA=IntervalMatrix([["[1.875,2.125]","[0.875,1.125]"],["[0.875,1.125]","[0.875,1.125]"]])

(-v,-iv)
(v+v,v+iv,iv+v,iv+iv)
(v-v,v-iv,iv-v,iv-iv)
(x*v,x*iv,ix*v,ix*iv)
(v*x,iv*x,v*ix,iv*ix)
(v/x,iv/x,v/ix,iv/ix)

(n*v,n*iv,v*n,iv*n,v/n,iv/n)
(d*v,d*iv,v*d,iv*d,v/d,iv/d)

(-cv,-icv)
(cv+cv,cv+icv,icv+cv,icv+icv)
(cv-cv,cv-icv,icv-cv,icv-icv)
(x*cv,x*icv,ix*cv,ix*icv)
(cv*x,icv*x,cv*ix,icv*ix)
(cv/x,icv/x,cv/ix,icv/ix)

(n*cv,n*icv,cv*n,icv*n,cv/n,icv/n)
(d*cv,d*icv,cv*d,icv*d,cv/d,icv/d)

(-A,-iA)
(A+A,A+iA,iA+A,iA+iA)
(A-A,A-iA,iA-A,iA-iA)
(x*A,x*iA,ix*A,ix*iA)
(A*x,iA*x,A*ix,iA*ix)
(A/x,iA/x,A/ix,iA/ix)

(n*A,n*iA,A*n,iA*n,A/n,iA/n)
(d*A,d*iA,A*d,iA*d,A/d,iA/d)

(cv*v,icv*v,cv*iv,icv*iv)
(A*v,iA*v,A*iv,iA*iv)
(cv*A,icv*A,cv*iA,icv*iA)
(A*A,iA*A,A*iA,iA*iA)
