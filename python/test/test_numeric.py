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
x=Float(2.25)
    
n=2
d=2.125
z=Integer(3)
x=Float(2.25)
ix=FuzzyFloat(2.00,2.25)
q=Rational(17,8)

(z+z,z-z,z*z)
(z+n,z-n,z*n)
(n+z,n-z,n*z)

(x+x,x-x,x*x,x/x)
(x+n,x-n,x*n,x/n)
(n+x,n-x,n*x,n/x)
(x+d,x-d,x*d,x/d)
(d+x,d-x,d*x,d/x)

(q+q,q-q,q*q,q/q)
(q+n,q-n,q*n,q/n)
(n+q,n-q,n*q,n/q)
(q+z,q-z,q*z,q/z)
(z+q,z-q,z*q,z/q)
(q+d,q-d,q*d,q/d)
(d+q,d-q,d*q,d/q)
(q+x,q-x,q*x,q/x)
(x+q,x-q,x*q,x/q)

(ix+ix,ix-ix,ix*ix,ix/ix)
(ix+n,ix-n,ix*n,ix/n)
(n+ix,n-ix,n*ix,n/ix)
(ix+d,ix-d,ix*d,ix/d)
(d+ix,d-ix,d*ix,d/ix)
(ix+x,ix-x,ix*x,ix/x)
(x+ix,x-ix,x*ix,x/ix)

(z==z,z!=z,z<z,z<=z,z>z,z>=z)
(z==n,z!=n,z<n,z<=n,z>n,z>=n)
(n==z,n!=z,n<z,n<=z,n>z,n>=z)

(x==x,x!=x,x<x,x<=x,x>x,x>=x)
(x==n,x!=n,x<n,x<=n,x>n,x>=n)
(n==x,n!=x,n<x,n<=x,n>x,n>=x)
(x==d,x!=d,x<d,x<=d,x>d,x>=d)
(d==x,d!=x,d<x,d<=x,d>x,d>=x)

(q==q,q!=q,q<q,q<=q,q>q,q>=q)
(q==n,q!=n,q<n,q<=n,q>n,q>=n)
(n==q,n!=q,n<q,n<=q,n>q,n>=q)
(q==z,q!=z,q<z,q<=z,q>z,q>=z)
(z==q,z!=q,z<q,z<=q,z>q,z>=q)
(q==d,q!=d,q<d,q<=d,q>d,q>=d)
(d==q,d!=q,d<q,d<=q,d>q,d>=q)
(q==x,q!=x,q<x,q<=x,q>x,q>=x)
(x==q,x!=q,x<q,x<=q,x>q,x>=q)
