#!/usr/bin/python

##############################################################################
#            test_import.py
#
#  Copyright 2006  Pieter Collins <Pieter.Collins@cwi.nl>
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

from numeric import *
print dir()

n=2
d=2.125
z=Integer(3)
x=Float(2.125)
ix=Interval(2.00,2.25)
q=Rational(17,8)

print (z+z,z-z,z*z)
print (z+n,z-n,z*n)
print (n+z,n-z,n*z)
print 
print (x+x,x-x,x*x,x/x)
print (x+n,x-n,x*n,x/n)
print (n+x,n-x,n*x,n/x)
print (x+d,x-d,x*d,x/d)
print (d+x,d-x,d*x,d/x)
print 
print (q+q,q-q,q*q,q/q)
print (q+n,q-n,q*n,q/n)
print (n+q,n-q,n*q,n/q)
print (q+z,q-z,q*z,q/z)
print (z+q,z-q,z*q,z/q)
print (q+d,q-d,q*d,q/d)
print (d+q,d-q,d*q,d/q)
print (q+x,q-x,q*x,q/x)
print (x+q,x-q,x*q,x/q)
print 
print (ix+ix,ix-ix,ix*ix,ix/ix)
print (ix+n,ix-n,ix*n,ix/n)
print (n+ix,n-ix,n*ix,n/ix)
print (ix+d,ix-d,ix*d,ix/d)
print (d+ix,d-ix,d*ix,d/ix)
print (ix+x,ix-x,ix*x,ix/x)
print (x+ix,x-ix,x*ix,x/ix)
print 
print (z==z,z!=z,z<z,z<=z,z>z,z>=z)
print (z==n,z!=n,z<n,z<=n,z>n,z>=n)
print (n==z,n!=z,n<z,n<=z,n>z,n>=z)
print 
print (x==x,x!=x,x<x,x<=x,x>x,x>=x)
print (x==n,x!=n,x<n,x<=n,x>n,x>=n)
print (n==x,n!=x,n<x,n<=x,n>x,n>=x)
print (x==d,x!=d,x<d,x<=d,x>d,x>=d)
print (d==x,d!=x,d<x,d<=x,d>x,d>=x)
print 
print (q==q,q!=q,q<q,q<=q,q>q,q>=q)
print (q==n,q!=n,q<n,q<=n,q>n,q>=n)
print (n==q,n!=q,n<q,n<=q,n>q,n>=q)
print (q==z,q!=z,q<z,q<=z,q>z,q>=z)
print (z==q,z!=q,z<q,z<=q,z>q,z>=q)
print (q==d,q!=d,q<d,q<=d,q>d,q>=d)
print (d==q,d!=q,d<q,d<=q,d>q,d>=q)
print (q==x,q!=x,q<x,q<=x,q>x,q>=x)
print (x==q,x!=q,x<q,x<=q,x>q,x>=q)
