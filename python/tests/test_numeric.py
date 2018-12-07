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

def test_regression():
    # New regression tests
    dp = DoublePrecision()
    mp = MultiplePrecision(52)
    # Questions: Should the following conversions from double be allowed?
    ExactDouble(1.3) # Can we construct ExactDouble when input probably not exact?
    Dyadic(1.375)
    FloatDP(1.375,dp)
    FloatMP(1.375,mp)
    #FloatDPValue(1.375,dp)
    #FloatMPValue(1.375,mp)

    z=Integer(1)
    w=Dyadic(1)
    q=Rational(1)

    # Check application of overloaded operators, especially on builtin types, returns the correct value
    assert(type(hlf(1))==Dyadic)
    print("type(rec(1)):",type(rec(1)))
    print("type(rec(z)):",type(rec(z)))

    assert(type(rec(1))==Rational)
    assert(type(rec(z))==Rational)
    assert(type(rec(z))==Rational)
    print("type(exp(1)):",type(exp(1)))
    assert(type(exp(1))==Real)
    assert(type(exp(z))==Real)

    ValidatedNumber({1:2})

    0*Dyadic(0)

    # Check construction from builtin integers
    assert(type(FloatDPBounds(1,2,dp))==FloatDPBounds);

    # Check construction from builtin floats
    assert(type(FloatDPBounds(1,2,dp))==FloatDPBounds);

    assert(type(FloatDPBounds({1:2},dp))==FloatDPBounds);

def check_arithmetic(x1,x2,r,qr=None):
    if qr==None: qr=r;
    assert(type(x1+x2)==type(r));
    assert(type(x1-x2)==type(r));
    assert(type(x1*x2)==type(r));
    assert(type(x1/x2)==type(qr));
    assert(type(x2+x1)==type(r));
    assert(type(x2-x1)==type(r));
    assert(type(x2*x1)==type(r));
    assert(type(x2/x1)==type(qr));
    assert(type(max(x1,x2))==type(r))
    assert(type(min(x1,x2))==type(r))
    assert(type(max(x2,x1))==type(r))
    assert(type(min(x2,x1))==type(r))

    assert(type(rec(x2))==type(qr))


def types(lst):
    return [type(val) for val in lst]

def test_algebraic():
    n=1
    z=Integer(n)
    w=Dyadic(z)
    d=Decimal(w)
    q=Rational(d)
    r=Real(q)

    print(types([n,w,max(n,w)]))

    check_arithmetic(n,z,z,q)
    check_arithmetic(n,w,w,q)
    check_arithmetic(n,q,q)
    check_arithmetic(n,r,r)
    check_arithmetic(z,z,z,q)
    check_arithmetic(z,w,w,q)
    check_arithmetic(z,q,q)
    check_arithmetic(z,r,r)
    check_arithmetic(w,w,w,q)
    check_arithmetic(w,q,q)
    check_arithmetic(w,r,r)
    check_arithmetic(q,q,q)
    check_arithmetic(q,r,r)
    check_arithmetic(r,r,r)



def test():
    test_algebraic()

    def exact(x): return Dyadic(FloatDPValue(FloatDP(x,DoublePrecision())))

    n=2
    d=2.125
    z=Integer(3)
    w=Dyadic(9,2)
    q=Rational(9,4)
    r=Real(q)

    dp=DoublePrecision()
    x=FloatDPApproximation(2.25,dp)
    bx=FloatDPBounds(dp)
    bx=FloatDPBounds(5,dp)
    bx=FloatDPBounds(exact(2.25),dp)
    bx=FloatDPBounds(exact(2.00),exact(2.25),dp)
#    bx=FloatDPBounds({exact(2.00):exact(2.25)},dp)


    ax=FloatDPApproximation(2.25,dp)
    bx=FloatDPBounds(exact(2.25),dp)
    vx=FloatDPValue(2.25)

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
#    (n+ax,n-ax,n*ax,n/ax)
    (ax+d,ax-d,ax*d,ax/d)
#    (d+ax,d-ax,d*ax,d/ax)
    (ax+bx,ax-bx,ax*bx,ax/bx)
    #(bx+ax,bx-ax,bx*ax,bx/ax)
#    (ax+r,ax-r,ax*r,ax/r)
    #(r+ax,r-ax,r*ax,r/ax)
    (ax+vx,ax-vx,ax*vx,ax/vx)
    #(vx+ax,vx-ax,vx*ax,vx/ax)


    (bx+bx,bx-bx,bx*bx,bx/bx)
    (bx+n,bx-n,bx*n,bx/n)
#    (n+bx,n-bx,n*bx,n/bx)
#    (bx+r,bx-r,bx*r,bx/r)
    #(r+bx,r-bx,r*bx,r/bx)
    (bx+vx,bx-vx,bx*vx,bx/vx)
    #(vx+bx,vx-bx,vx*bx,vx/bx)

    (r+r,r-r,r*r,r/r)
    (r+n,r-n,r*n,r/n)
#    (n+r,n-r,n*r,n/r)
#    (r+vx,r-vx,r*vx,r/vx)
    #(vx+r,vx-r,vx*r,vx/r)

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
#    (ax==n,ax!=n,ax<n,ax<=n,ax>n,ax>=n)
#    (n==ax,n!=ax,n<ax,n<=ax,n>ax,n>=ax)
#    (ax==d,ax!=d,ax<d,ax<=d,ax>d,ax>=d)
#    (d==ax,d!=ax,d<ax,d<=ax,d>ax,d>=ax)
