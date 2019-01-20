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
pr = DoublePrecision()


# Print a list of all available classes and functions
print(dir(),"\n")

###############################################################################
# Numeric submodule
###############################################################################

# Number classes and their constructors

# Create a double-precision number from an exact value
d=ExactDouble(1.375)
d=exact(1.375)
print("d:",d,type(d))

# Create an integer
z=Integer(5)
print("z:",z,type(z))

# Create a dyadic; can convert from ExactDouble and Integer
w=Dyadic(d)
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
print("r:",r,type(r))


# Define shorthands for defining Ariadne values from input
def dec(x): return Decimal(x)
def ex(x): return ExactDouble(x)
def dy(x): return Dyadic(ExactDouble(x))

# Specify precisions
mpr=MultiplePrecision(128)
dpr=DoublePrecision()


# Create double-precision bounds for a value
x=FloatDPBounds(dec(1.2),dpr) # Creates the interval [1.19999999999999996:1.20000000000000018]
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

###############################################################################
# Linear Algebra submodule
###############################################################################

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

###############################################################################
# Function submodule
###############################################################################

print(dir(EffectiveScalarMultivariateFunction))


# Create a user-defined scalar-valued function
argument_size=2
x=[EffectiveScalarMultivariateFunction.coordinate(argument_size,index) for index in range(0,argument_size)]
f=sqrt(sqr(x[0])+sqr(x[1]))
print(f)
print(f(FloatDPApproximationVector([4,3],dpr)))
print(f(FloatDPBoundsVector([4,3],dpr)))

result_size=2
argument_size=3
x=[EffectiveScalarMultivariateFunction.coordinate(argument_size,index) for index in range(0,argument_size)]
f=EffectiveVectorMultivariateFunction([sqrt(sqr(x[0])+sqr(x[1])),x[1]+x[2]])
print(f)
print(f(FloatDPApproximationVector([4,3,0],dpr)))
print(evaluate(f,(FloatDPBoundsVector([4,3,0],dpr))))



###############################################################################
# Calculus submodule
###############################################################################

# Create a box to act as the domain of a Taylor function
d=ExactBox([{4:7},{1:6},{-1:+1}])
print("d:",d)

# Create a sweeper to control the accuracy of a Taylor function
swp=ThresholdSweeper(pr,1e-8)
#swp=GradedSweeper(pr,6)
print("swp:",swp)

# Create the scalar Taylor model representing the function x1 on the domain d
t1=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(d,1,swp)
print("t1:",t1,"\n")

# Create all Taylor variables on the domain dom
t=[ValidatedScalarMultivariateTaylorFunctionModel.coordinate(d,i,swp) for i in range(0,3)]
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
def t(d,j):
    return ValidatedScalarMultivariateTaylorFunctionModel.coordinate(d,j,swp)

#Create a ValidatedScalarMultivariateTaylorFunctionModel from a Polynomial
p=EffectiveScalarMultivariateFunction.coordinate(3,0)
tp=ValidatedScalarMultivariateTaylorFunctionModel(d,p,swp)
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
n=int(3); c=FloatDPValue(Dyadic(1.5),pr); b=FloatDPBounds(Dyadic(1.25),Dyadic(1.75),pr)

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
tx=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(d,0,swp);
ty=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(d,1,swp);

# Comparison functions
min(tx,ty); max(tx,ty); abs(tx);

print("tx:",tx,type(tx))
# Univariate arithmetical functions
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
d1=ExactBox([{-1:1},{-1:1}])
d2=ExactBox([{dy(-0.125):dy(0.125)},{dy(0.5):dy(0.75)}])
w=t(d1,0)*t(d1,1)
rw=restrict(w,d2)
print("restrict(w,d2):",rw)

# Embed the domain of x in a space of higher dimension
d=ExactBox([{-1:+1}])
    #ex=embed(x,d)
    #print("embed(x,d):",ex)

#Join tx and ty into a TaylorFunction
tf=join(tx,ty)
print("join(tx,ty):",tf)

tg=combine(tx,ty)
print("combine(tx,ty):",tg)
print

# Function composition
d=ExactBox([{4:7},{1:6},{-1:1}])
th=ValidatedVectorMultivariateTaylorFunctionModel.identity(d,swp)
cd=th.codomain()

##f=VectorFunction.affine(FloatMatrix([[2,1,0],[1,1,1]]),FloatVector([1,1]))
f=EffectiveVectorMultivariateFunction(2,3)
g=EffectiveScalarMultivariateFunction.coordinate(3,0)

tg=ValidatedScalarMultivariateTaylorFunctionModel(cd,g,swp)
tf=ValidatedVectorMultivariateTaylorFunctionModel(cd,f,swp)

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
d=ExactBox([{-1:+1},{-1:+1},{-1:+1}])
x=EffectiveScalarMultivariateFunction.coordinate(3,0)
y0=EffectiveScalarMultivariateFunction.coordinate(3,1)
y1=EffectiveScalarMultivariateFunction.coordinate(3,2)
f=join(x+4*y0+y1,y0+y1)
slv=IntervalNewtonSolver(1e-8,6)
print("f:",f)
h=slv.implicit(f,ExactBox([{-1:+1}]),ExactBox([{-1:+1},{-1:+1}]))
print("implicit(f):",h)

# Compute the solution h to the scalar equation g(x,h(x))=0
# with f(x,y)=4+x-y^2, so y=sqrt(4+x)
d=ExactBox([{-1:+1},{-1:+1}])
x=EffectiveScalarMultivariateFunction.coordinate(2,0)
y=EffectiveScalarMultivariateFunction.coordinate(2,1)
g=x-4*y+y*y
print("g:",g)
h=slv.implicit(g, ExactBox([{-1:+1}]),ExactInterval({-1:+1}))
print("implicit(g):",h)
print()

# Differentiation, integration and differential equations
x=t(d,0)
y=t(d,1)
f=1+2*x+3*x*y+x*x
print("f:",f)

# Compute the derivative of f with respect to x[j], assuming that the error is constant
df0=derivative(f,0)
df1=derivative(f,1)
print("derivative(f,0):",df0)
print("derivative(f,1):",df1)

# Compute an antiderivative of f with respect to x[j]
# The value is zero at the midpoint of the domain of x[j]
if0=antiderivative(f,0)
if1=antiderivative(f,1)
print("antiderivative(f,0):",if0)
print("antiderivative(f,1):",if1)
print()

# Compute the flow of the Taylor function f starting in the domain D for time interval [-h,+h]
b=ExactBox([{-4:4}])
d=ExactBox([{-1:+1}])
h=1/two
o=8 # Temporal order
f=ValidatedVectorMultivariateTaylorFunctionModel.identity(b,swp)
f=ValidatedVectorMultivariateFunction.identity(1)
integrator=GradedTaylorSeriesIntegrator(1e-8)
phis=integrator.flow(f,d,h)
phi=phis[0]
print("phi:",phi,type(phi))
print()

# Compute two time steps of the flow of the Taylor function f starting in domain D for the interval [h,2h]
phi0=phi
print("phi.domain():",phi.domain(),", h:",h)
phi0h=partial_evaluate(phi,1,FloatDPBounds(h,pr))
dd=phi0h.codomain()
phi=integrator.flow(f,dd,h)[0]
tr=ValidatedScalarMultivariateTaylorFunctionModel.coordinate([{0:2*h}],0,swp)-h
tr=ValidatedScalarMultivariateFunctionModel(tr)
print(type(phi0),type(tr))
phi1=compose(phi,combine(phi0h,tr))
print("phi0:",phi0)
print("phi1:",phi1)
print()



## Contractors

x1=x+x*x/2
x2=x+FloatDPBounds(-1,1,pr)
print("x1:",x1)
print("x2:",x1)
b=refines(x1,x2)
print("refines(x1,x2):",b)
b=inconsistent(x1,x2)
print("inconsistent(x1,x2):",b)
y=refinement(x1,x2)
print("refinement(x1,x2):",y)
print()

# Set up algebraic equation
d1=ExactBox([{-1:1}])
d2=ExactInterval(-1,1)
d=product(d1,d2)
x=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(d,0,swp)
y=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(d,1,swp)
f=x-4*y+y*y
i=ValidatedVectorMultivariateTaylorFunctionModel.identity(d1,swp)
h=ValidatedScalarMultivariateTaylorFunctionModel.constant(d1,FloatDPBounds(-1,+1,pr),swp)

# Newton contractor to solve f(x,h(x))=0 for scalar f
df1 = derivative(f,1).range()
df1 = FloatDPBounds(df1.lower(),df1.upper())
h=refinement(compose(f,join(i,h))/df1,h)
print("h:",phi)
print()


# Set up differential equation
b=ExactBox([{-1:1},{-1:1}])
o=ValidatedScalarMultivariateTaylorFunctionModel.constant(b,1,swp)
x=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(b,0,swp)
y=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(b,1,swp)
f=join(o,x) # [dot(x),dot(y)]=[1,x]
d=ExactBox([{0:dy(0.125)},{0:dy(0.125)}])
h=Dyadic(exact(0.5))
d0=product(d,ExactInterval(-h,+h))
x0=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(d0,0,swp)
y0=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(d0,1,swp)
t=ValidatedScalarMultivariateTaylorFunctionModel.coordinate(d0,2,swp)
phi=join(x0,y0)

# Picard operator to solve dot(phi)(x,t) = f(phi(x,t))
phi=antiderivative(compose(f,phi),2)
print("phi:",phi)
print

# End
