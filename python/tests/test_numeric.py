#!/usr/bin/python

##############################################################################
#            test_numeric.py
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
Float=float
n=2
d=2.125
z=Integer(3)
q=Rational(17,8)
x=Float(2.25)
ix=Interval()
ix=Interval(5)
ix=Interval(2.25)
ix=Interval(2.00,2.25)
ix=Interval({2.00:2.25})

ax=Float(2.25)
ix=Interval(2.25)
sx=Real(2.25)
ex=ExactFloat(2.25)

#(z+z,z-z,z*z)
#(z+n,z-n,z*n)
#(n+z,n-z,n*z)

#(q+q,q-q,q*q,q/q)
#(q+n,q-n,q*n,q/n)
#(n+q,n-q,n*q,n/q)
#(q+z,q-z,q*z,q/z)
#(z+q,z-q,z*q,z/q)
#(q+d,q-d,q*d,q/d)
#(d+q,d-q,d*q,d/q)

(ax+ax,ax-ax,ax*ax,ax/ax)
(ax+n,ax-n,ax*n,ax/n)
(n+ax,n-ax,n*ax,n/ax)
(ax+d,ax-d,ax*d,ax/d)
(d+ax,d-ax,d*ax,d/ax)
(ax+ix,ax-ix,ax*ix,ax/ix)
#(ix+ax,ix-ax,ix*ax,ix/ax)
(ax+sx,ax-sx,ax*sx,ax/sx)
#(sx+ax,sx-ax,sx*ax,sx/ax)
(ax+ex,ax-ex,ax*ex,ax/ex)
#(ex+ax,ex-ax,ex*ax,ex/ax)


(ix+ix,ix-ix,ix*ix,ix/ix)
(ix+n,ix-n,ix*n,ix/n)
(n+ix,n-ix,n*ix,n/ix)
(ix+d,ix-d,ix*d,ix/d)
(d+ix,d-ix,d*ix,d/ix)
(ix+sx,ix-sx,ix*sx,ix/sx)
#(sx+ix,sx-ix,sx*ix,sx/ix)
(ix+ex,ix-ex,ix*ex,ix/ex)
#(ex+ix,ex-ix,ex*ix,ex/ix)

(sx+sx,sx-sx,sx*sx,sx/sx)
(sx+n,sx-n,sx*n,sx/n)
(n+sx,n-sx,n*sx,n/sx)
(sx+d,sx-d,sx*d,sx/d)
(d+sx,d-sx,d*sx,d/sx)
(sx+ex,sx-ex,sx*ex,sx/ex)
#(ex+sx,ex-sx,ex*sx,ex/sx)

#(z==z,z!=z,z<z,z<=z,z>z,z>=z)
#(z==n,z!=n,z<n,z<=n,z>n,z>=n)
#(n==z,n!=z,n<z,n<=z,n>z,n>=z)

#(q==q,q!=q,q<q,q<=q,q>q,q>=q)
#(q==n,q!=n,q<n,q<=n,q>n,q>=n)
#(n==q,n!=q,n<q,n<=q,n>q,n>=q)
#(q==z,q!=z,q<z,q<=z,q>z,q>=z)
#(z==q,z!=q,z<q,z<=q,z>q,z>=q)
#(q==d,q!=d,q<d,q<=d,q>d,q>=d)
#(d==q,d!=q,d<q,d<=q,d>q,d>=q)

(ax==ax,ax!=ax,ax<ax,ax<=ax,ax>ax,ax>=ax)
(ax==n,ax!=n,ax<n,ax<=n,ax>n,ax>=n)
(n==ax,n!=ax,n<ax,n<=ax,n>ax,n>=ax)
(ax==d,ax!=d,ax<d,ax<=d,ax>d,ax>=d)
(d==ax,d!=ax,d<ax,d<=ax,d>ax,d>=ax)

