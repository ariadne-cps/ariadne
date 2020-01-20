#!/usr/bin/python

# -*- coding: utf-8 -*-

##############################################################################
#            tutorial.py
#
#  Copyright 2009-18  Pieter Collins <Pieter.Collins@cwi.nl>
##############################################################################

# This file is part of Ariadne.

# Ariadne is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Ariadne is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Ariadne. If not, see <https://www.gnu.org/licenses/>.

# Import all classes in the ariadne module
from ariadne import *

# Print a list of all available classes and functions
print(dir())

# Output a
###############################################################################
## [Numeric demonstration]

# Number classes and their constructors

# Create an integer
z=Integer(5)
print("z:",z,type(z))

# Create a dyadic; can convert from Integer
w=Dyadic(z)
w=Dyadic(5)
w=Dyadic(11,3)
w=11/two**3
print("w:",w,type(w))

# Create a decimal number; can convert from Dyadic
g=Decimal(w)
g=Decimal(9.81)
g=Decimal("9.81")
print("g:",g,type(g))

# Create a rational; can convert from Dyadic, Decimal
q=Rational(w);
q=Rational(5);
q=Rational(11,8);
print("q:",q,type(q))

# Create a real number
r=Real(q)
r=Real(sin(q))
print("r:",r,type(r))


# Store a double-precision floating-point number
d=ExactDouble(1.375)
d=exact(1.375)
print("d:",d,type(d))
# Can convert an ExactDouble to a Dyadic number.
w=Dyadic(d)

# Define shorthands for defining Ariadne values from input
def dec(x): return Decimal(x)
def ex(x): return ExactDouble(x)
def dy(x): return Dyadic(ExactDouble(x))


# Specify precisions of floating-point number types
dpr=DoublePrecision()
dpr=double_precision
mpr=MultiplePrecision(128)

# Create a raw double-precision number
x=FloatDP(Dyadic(exact(1.75)),dpr)
# Create a raw multiple-precision number
x=FloatMP(Dyadic(exact(1.75)),mpr)

# Create double-precision bounds for a value
x=FloatDPBounds(Decimal(1.2),dpr) # Creates the interval [1.19999999999999996:1.20000000000000018]
print("FloatDPBounds(1.2):",x)

# Create double-precision bounds for a range of values
x=FloatDPBounds(Rational(11,10),Rational(14,10),dpr) # Creates the interval [1.09999999999999987:1.40000000000000013]
print("FloatDPBounds(11/10,14/10,dpr):",x)

# Create multiple-precision bounds for a value
x=FloatMPBounds(dec("1.2"),mpr) # Creates the interval [1.19999999999999996:19999999999999996]
print("FloatMPBounds(1.2,mpr):",x)

# Create multiple-precision bounds for a range of values
x=FloatMPBounds(dec(1.5),dec(2.25),mpr) # Creates the interval [1.5,2.25]
print("FloatMPBounds(1.5,2.25,mpr):",x)

# Create multiple-precision bounds for a range of values
x=FloatMPBounds(Rational(11,10),Rational(14,10),mpr) # Creates the interval [1.10000000000000009:1.39999999999999991]
print("FloatMPBounds(11/10,14/10,mpr):",x)

# Create a double-precision approximation
xa=FloatDPApproximation(1.23,dpr)
print("FloatDPApproximation(1.23,dpr):",xa)

xa=FloatMPApproximation(1.23,mpr)
print("FloatMPApproximation(1.23,dpr):",xa)
xa=FloatMPApproximation(Decimal("1.23"),mpr)
print("FloatMPApproximation(Decimal(\"1.23\"),dpr):",xa)
print("\n")

## [Numeric demonstration]
###############################################################################
## [Linear Algebra demonstration]

# Create an interval vector
b=FloatDPBoundsVector([1,{2:3},4],dpr)
print("b:",b)

Aq=RationalMatrix([[1,2,4],[3,1,2],[0,0,1]])
print("Aq:",Aq)

# Create an interval matrix
A=FloatDPBoundsMatrix([[1,2,4],[3,1,2],[0,0,1]],dpr)

A=FloatDPBoundsMatrix([[1,{2:3},4],[dy(1.5),dec(1.1),2],[0,0,1]],dpr)

print("A:",A)

# Solve the linear equation Ax=b
x=solve(A,b)
print("A\\b:",x)

print("\n")

## [Linear Algebra demonstration]
###############################################################################
## [Function demonstration]

# Create a user-defined scalar-valued function
argument_size=2
x=[EffectiveScalarMultivariateFunction.coordinate(argument_size,index) for index in range(0,argument_size)]
f=sqrt(sqr(x[0])+sqr(x[1]))
print("f:",f)
a=FloatDPApproximationVector([4,3],dpr)
print("f(a):",f(a))
b=FloatDPApproximationVector([4,3],dpr)
print("f(b):",f(b))

result_size=2
argument_size=3
x=[EffectiveScalarMultivariateFunction.coordinate(argument_size,index) for index in range(0,argument_size)]
f=EffectiveVectorMultivariateFunction([sqrt(sqr(x[0])+sqr(x[1])),x[1]+x[2]])
print("f:",f)
a=FloatDPApproximationVector([4,3,0],dpr)
print("f(a):",f(a))
b=FloatDPBoundsVector([4,3,0],dpr)
print("evaluate(f,b):",evaluate(f,b))

p0=FloatDPApproximationMultivariatePolynomial.coordinate(2,0);
p1=FloatDPApproximationMultivariatePolynomial.coordinate(2,1);
print("p1:",p1)
a3=FloatDPApproximation(3,dpr)
a5=FloatDPApproximation(5,dpr)
print(p0*(p0+p1*a3)+a5)

print("\n")

## [Function demonstration]
###############################################################################
## [Calculus demonstration]

# Create a box to act as the domain of a Taylor function
dom=ExactBoxType([{4:7},{1:6},{-1:+1}])
print("dom:",dom)

# Create a sweeper to control the accuracy of a Taylor function
swp=ThresholdSweeper(dpr,1e-8)
#swp=GradedSweeper(dpr,6)
print("swp:",swp)

# Create the scalar Taylor model representing the function x1 on the domain dom
t1=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom,1,swp)
print("t1:",t1,"\n")

# Create all Taylor variables on the domain dom
t=[ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom,i,swp) for i in range(0,3)]
print("t[0]:",t[0])
print("t[1]:",t[1])
print("t[2]:",t[2])
print()

# Make shorthands for the variable names
tx=t[0]
ty=t[1]
tz=t[2]
print("tx:",tx,"\nty:",ty,"\ntz:",tz,"\n")

# Make a shorthand for constructing Taylor expressions
def T(dom,j):
    return ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom,j,swp)

#Create a ValidatedScalarMultivariateTaylorFunctionModelDP from a Polynomial
p=EffectiveScalarMultivariateFunction.coordinate(3,0)
tp=ValidatedScalarMultivariateTaylorFunctionModelDP(dom,p,swp)
print("p:",p,"\ntp: ",tp,"\n")

# The domain D of tx.
tx.domain()
# A not-necessarily tight over-approximation to p(D)+/-e.
tx.codomain()
# An over-approximation to p(D)+/-e.
tx.range()
### Convert to a polynomial expression.
##tx.polynomial()

# Arithmetic on Taylor models
+tx; -tx; tx+ty; tx-ty; tx*ty; tx/ty;

# Define some constants
n=int(3); c=FloatDPValue(Dyadic(exact(1.5)),dpr); b=FloatDPBounds(Dyadic(exact(1.25)),Dyadic(exact(1.75)),dpr);
# Mixed arithmetic
tx+c; tx-c; tx*c; tx/c;
    #c+tx; c-tx; c*tx; c/tx;

tx+b; tx-b; tx*b; tx/b;
    #b+tx; b-tx; b*tx; b/tx;

# Inplace operations
tx+=ty; tx-=ty;
tx+=c; tx-=c; tx*=c; tx/=c;
tx+=b; tx-=b; tx*=b; tx/=b;

# Reset tx,ty
tx=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom,0,swp);
ty=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom,1,swp);

# Comparison functions
min(tx,ty); max(tx,ty); abs(tx);

print("tx:",tx,type(tx))
# Univariate arithmetical functions
neg(tx);
rec(tx);
sqr(tx);
pow(tx,n);

# Algebraic and transcendental functions
sqrt(tx);
exp(tx); log(tx);
sin(tx), cos(tx), tan(tx/2)
print()

# Non-arithmetic functions

# Restrict to a subdomain
d1=ExactBoxType([{-1:1},{-1:1}])
d2=ExactBoxType([{dy(-0.125):dy(0.125)},{dy(0.5):dy(0.75)}])
w=T(d1,0)*T(d1,1)
rw=restrict(w,d2)
print("restrict(w,d2):",rw)

# Embed the domain of tx in a space of higher dimension
dom=ExactBoxType([{-1:+1}])
    #etx=embed(tx,dom)
    #print("embed(tx,dom):",etx)

#Join tx and ty into a TaylorFunction
tf=join(tx,ty)
print("join(tx,ty):",tf)

tg=combine(tx,ty)
print("combine(tx,ty):",tg)
print

# Function composition
dom=ExactBoxType([{4:7},{1:6},{-1:1}])
th=ValidatedVectorMultivariateTaylorFunctionModelDP.identity(dom,swp)
cd=th.codomain()

##f=VectorFunction.affine(FloatMatrix([[2,1,0],[1,1,1]]),FloatVector([1,1]))
f=EffectiveVectorMultivariateFunction(2,3)
g=EffectiveScalarMultivariateFunction.coordinate(3,0)

tg=ValidatedScalarMultivariateTaylorFunctionModelDP(cd,g,swp)
tf=ValidatedVectorMultivariateTaylorFunctionModelDP(cd,f,swp)

# Compose an scalar multivariate function and a vector Taylor function
compose(g,th)
print("compose(g,th):",compose(g,th))

# Compose a vector multivariate function and a vector Taylor function
compose(f,th)
print("compose(f,th):",compose(f,th))

# Compose two Taylor functions
assert(subset(th.codomain(),tf.domain()))
tr=compose(tf,th)
print("compose(tf,th):",compose(tf,th))

# Compose a Taylor expression and a Taylor function
assert(subset(th.codomain(),tg.domain()))
tr=compose(tg,th)
print("compose(tg,th):",compose(tg,th))
print()



# Solution of parameterised algebraic equations

# Compute the solution h to the vector equation f(x,h(x))=0
dom=ExactBoxType([{-1:+1},{-1:+1},{-1:+1}])
x=EffectiveScalarMultivariateFunction.coordinate(3,0)
y0=EffectiveScalarMultivariateFunction.coordinate(3,1)
y1=EffectiveScalarMultivariateFunction.coordinate(3,2)
f=join(x+4*y0+y1,y0+y1)
slv=IntervalNewtonSolver(1e-8,6)
print("f:",f)
h=slv.implicit(f,ExactBoxType([{-1:+1}]),ExactBoxType([{-1:+1},{-1:+1}]))
print("implicit(f):",h)

# Compute the solution h to the scalar equation g(x,h(x))=0
# with f(x,y)=4+x-y^2, so y=sqrt(4+x)
dom=ExactBoxType([{-1:+1},{-1:+1}])
x=EffectiveScalarMultivariateFunction.coordinate(2,0)
y=EffectiveScalarMultivariateFunction.coordinate(2,1)
g=x-4*y+y*y
print("g:",g)
h=slv.implicit(g, ExactBoxType([{-1:+1}]),ExactIntervalType({-1:+1}))
print("implicit(g):",h)
print()

# Differentiation, integration and differential equations
x0=T(dom,0)
x1=T(dom,1)
f=1+2*x0+3*x0*x1+x0*x0
print("f:",f)

# Compute the derivative of f with respect to x[j], assuming that the error is constant
dtf0=derivative(f,0)
dtf1=derivative(f,1)
print("derivative(f,0):",dtf0)
print("derivative(f,1):",dtf1)

# Compute an antiderivative of f with respect to x[j]
# The value is zero at the midpoint of the domain of x[j]
Itf0=antiderivative(f,0)
Itf1=antiderivative(f,1)
print("antiderivative(f,0):",Itf0)
print("antiderivative(f,1):",Itf1)
print()

# Compute the flow of the Taylor function f starting in the domain dom for time interval [-h,+h]
dom=ExactBoxType([{-1:+1}])
bbx=ExactBoxType([{-4:+4}])
h=1/two
o=8 # Temporal order
f=ValidatedVectorMultivariateTaylorFunctionModelDP.identity(bbx,swp)
f=ValidatedVectorMultivariateFunction.identity(1)
integrator=GradedTaylorSeriesIntegrator(1e-8)

phis=integrator.flow(f,dom,h)
assert(len(phis)==1)
phi=phis[0]
print("phi:",phi,type(phi))
print()

# Compute two time steps of the flow of the Taylor function f starting in domain D for the interval [h,2h]
phi0=phi
print("phi.domain():",phi.domain(),", h:",h)
phi0h=partial_evaluate(phi,1,FloatDPBounds(h,dpr))
dom1=phi0h.codomain()
phi=integrator.flow(f,dom1,h)[0]
print("phi:",phi)
tr=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate([{0:2*h}],0,swp)-h
tr=ValidatedScalarMultivariateFunctionModelDP(tr)
phi1=compose(phi,combine(phi0h,tr))
print("phi0:",phi0)
print("phi1:",phi1)
print()



## Contractors

f1=x0+x0*x0/2
f2=x0+FloatDPBounds(-1,1,dpr)
print("f1:",f1)
print("f2:",f1)
b=refines(f1,f2)
print("refines(f1,f2):",b)
b=inconsistent(f1,f2)
print("inconsistent(f1,f2):",b)
g=refinement(f1,f2)
print("refinement(f1,f2):",g)
print()

# Set up algebraic equation
dom1=ExactBoxType([{-1:1}])
dom2=ExactIntervalType(-1,1)
dom=product(dom1,dom2)
x=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom,0,swp)
y=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom,1,swp)
f=x-4*y+y*y
i=ValidatedVectorMultivariateTaylorFunctionModelDP.identity(dom1,swp)
h=ValidatedScalarMultivariateTaylorFunctionModelDP.constant(dom1,FloatDPBounds(-1,+1,dpr),swp)

# Newton contractor to solve f(x,h(x))=0 for scalar f
df1 = derivative(f,1).range()
df1 = FloatDPBounds(df1.lower(),df1.upper())
h=refinement(compose(f,join(i,h))/df1,h)
print("h:",phi)
print()


# Set up differential equation
bbx=ExactBoxType([{-1:1},{-1:1}])
o=ValidatedScalarMultivariateTaylorFunctionModelDP.constant(bbx,1,swp)
x=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(bbx,0,swp)
y=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(bbx,1,swp)
f=join(o,x) # [dot(x),dot(y)]=[1,x]
dom=ExactBoxType([{0:dy(0.125)},{0:dy(0.125)}])
h=Dyadic(exact(0.5))
dom0=product(dom,ExactIntervalType(-h,+h))
x0=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom0,0,swp)
y0=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom0,1,swp)
t=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom0,2,swp)
phi=join(x0,y0)

# Picard operator to solve dot(phi)(x,t) = f(phi(x,t))
phi=antiderivative(compose(f,phi),2)
print("phi:",phi)
print

## [Calculus demonstration]
###############################################################################
