#!/usr/bin/python

# -*- coding: utf-8 -*-

##############################################################################
#            tutorial.py
#
#  Copyright 2009  Pieter Collins <Pieter.Collins@cwi.nl>
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
print dir(),"\n"

###############################################################################
# Numeric submodule
###############################################################################

# Create an integer
z=Integer(8)
print z

# Create a rational
q=Rational(8,6)
print q

# Create an exact float (dyadic)
x=Dyadic(1.25)
print x

# Create a decimal number
g=Decimal(9.81)
g=Decimal("9.81")
print g

# Create a real number
r=sin(1)
print r

# Create an interval
i=Interval(1.5,2.25) # Creates the interval [1.5,2.25]
print "Interval(1.5,2.25):",i

# Create an interval
i=Interval(1.2) # Creates the interval [1.19999999999999996:19999999999999996]
print "Interval(1.2):",i
#i=Interval("1.2") # Creates the interval [1.19999999999999996:19999999999999996]
#print "Interval(\"1.2\"):",i

# Create an interval
i=Interval(1.1,1.4) # Creates the interval [1.10000000000000009:1.39999999999999991]
print "Interval(1.1,1.4):",i

# Create an interval
i=Interval(Rational(12,10)) # Creates the interval [1.19999999999999996:1.20000000000000018]
print "Interval(Rational(12,10)):",i

# Create an interval
i=Interval(Rational(11,10),Rational(14,10)) # Creates the interval [1.09999999999999987:1.40000000000000013]
print "Interval(Rational(11,10),Rational(14,10)):",i

# Create an approximate float
x=Float(1.23)
print x

###############################################################################
# Linear Algebra submodule
###############################################################################

# Create an interval vector
b=IntervalVector([1,{2:3},4])
print "b:",b

# Create an interval matrix
A=IntervalMatrix([[1,{2:3},4],[{1.5:1.5},Decimal("1.1"),2],[0,0,1]])
print "A:",A

# Solve the linear equation Ax=b
solve(A,b)

print "\n\n"

###############################################################################
# Function submodule
###############################################################################

print dir(RealScalarFunction)

# Create a user-defined scalar-valued function
argument_size=2
x=[RealScalarFunction.coordinate(argument_size,index) for index in range(0,argument_size)]
f=sqrt(sqr(x[0])+sqr(x[1]))
print f
print f(FloatVector([4,3]))
print f(IntervalVector([4,3]))

result_size=2
argument_size=3
x=[RealScalarFunction.coordinate(argument_size,index) for index in range(0,argument_size)]
f=RealVectorFunction([sqrt(sqr(x[0])+sqr(x[1])),x[1]+x[2]])
print f
print f(FloatVector([4,3,0]))
print evaluate(f,(IntervalVector([4,3,0])))



###############################################################################
# Calculus submodule
###############################################################################


# Create a box to act as the domain of a Taylor function
d=IntervalVector([{4:7},{1:6},{-1:+1}])
print "d:",d

# Create a sweeper to control the accuracy of a Taylor function
swp=ThresholdSweeper(1e-8)
print "swp:",swp

# Create the scalar Taylor model representing the function x1 on the domain d
t=ScalarTaylorFunction.coordinate(d,1,swp)
print "t:",t,"\n"

# Create all Taylor variables on the domain dom
tv=[ScalarTaylorFunction.coordinate(d,i,swp) for i in range(0,3)]
print "tv[0]:",tv[0]
print "tv[1]:",tv[1]
print "tv[2]:",tv[2]
print

# Make shorthands for the variable names
x=tv[0]
y=tv[1]
z=tv[2]
print "x:",x,"\ny:",y,"\nz:",z,"\n"

# Make a shorthand for constructing Taylor expressions
def t(d,j):
    return ScalarTaylorFunction.coordinate(d,j,swp)

#Create a ScalarTaylorFunction from a Polynomial
p=RealScalarFunction.coordinate(3,0)
tp=ScalarTaylorFunction(d,p,swp)
print "p:",p,"\ntp: ",tp,"\n"

# The domain D of x.
x.domain()
# A not-necessarily tight over-approximation to p(D)+/-e.
x.codomain()
# An over-approximation to p(D)+/-e.
x.range()
### Convert to a polynomial expression.
##x.polynomial()

# Arithmetic on Taylor models
+x; -x; x+y; x-y; x*y; x/y;

# Define some constants
n=int(-3); c=Float(1.5); i=Interval(1.25,1.75)

# Mixed arithmetic
x+c; x-c; x*c; x/c;
c+x; c-x; c*x; c/x;

x+i; x-i; x*i; x/i;
i+x; i-x; i*x; i/x;

# Inplace operations
x+=y; x-=y;
x+=c; x-=c; x*=c; x/=c;
x+=i; x-=i; x*=i; x/=i;

# Reset x,y
x=ScalarTaylorFunction.coordinate(d,0,swp);
y=ScalarTaylorFunction.coordinate(d,1,swp);

# Comparison functions
min(x,y); max(x,y); abs(x);

# Compound arithmetical functions
neg(x); rec(x); sqr(x); pow(x,n);

# Algebraic and transcendental functions
sqrt(x);
exp(x); #log(x);
sin(x), cos(x), tan(x/100)
print

# Non-arithmetic functions

# Restrict to a subdomain
d1=Box([[-1,1],[-1,1]])
d2=Box([[-0.125,0.125],[0.5,0.75]])
w=t(d1,0)*t(d1,1)
rw=restrict(w,d2)
print "restrict(w,d2):",rw

# Embed the domain of x in a space of higher dimension
d=Box([[-1,+1]])
ex=embed(x,d)
print "embed(x,d):",ex

#Join x and y into a TaylorFunction
f=join(x,y)
print "join(x,y):",f

g=combine(x,y)
print "combine(x,y):",g
print

# Function composition
d=IntervalVector([[4,7],[1,6],[-1,1]])
th=VectorTaylorFunction.identity(d,swp)
cd=th.codomain()

##f=VectorFunction.affine(FloatMatrix([[2,1,0],[1,1,1]]),FloatVector([1,1]))
f=RealVectorFunction(2,3)
g=RealScalarFunction.coordinate(3,0)

tg=ScalarTaylorFunction(cd,g,swp)
tf=VectorTaylorFunction(cd,f,swp)

# Compose an expression and a Taylor function
compose(g,th)
print "compose(g,th):",compose(g,th)

# Compose a function and a Taylor function
compose(f,th)
print "compose(f,th):",compose(f,th)

# Compose two Taylor functions
assert(subset(th.codomain(),tf.domain()))
tr=compose(tf,th)
print "compose(tf,th):",compose(tf,th)

# Compose a Taylor expression and a Taylor function
assert(subset(th.codomain(),tg.domain()))
tr=compose(tg,th)
print "compose(tg,th):",compose(tg,th)
print

# Solution of parameterised algebraic equations

# Compute the solution h to the vector equation f(x,h(x))=0
d=Box([[-1,1],[-1,1],[-1,1]])
a=t(d,0)
x=t(d,1)
y=t(d,2)
f=join(a+4*x+y,x+y)
slv=IntervalNewtonSolver(1e-8,6,)
print "f:",f
h=slv.implicit(f,Box([[-1,1],[-1,1]]),Box([[-1,1]]))
print "implicit(f):",h

# Compute the solution h to the scalar equation g(x,h(x))=0
# with f(x,y)=4+x-y^2, so y=sqrt(4+x)
d=Box([[-1,1],[-1,1]])
x=t(d,0)
y=t(d,1)
g=x-4*y+y*y
print "g:",g
h=implicit(g)
print "implicit(g):",h
print

# Differentiation, integration and differential equations
f=1+2*x+3*x*y+x*x
print "f:",f

# Compute the derivative of f with respect to x[j], assuming that the error is constant
df0=derivative(f,0)
df1=derivative(f,1)
print "derivative(f,0):",df0
print "derivative(f,1):",df1

# Compute an antiderivative of f with respect to x[j]
# The value is zero at the midpoint of the domain of x[j]
if0=antiderivative(f,0)
if1=antiderivative(f,1)
print "antiderivative(f,0):",if0
print "antiderivative(f,1):",if1

# Compute the flow of the Taylor function f starting in the domain D for time interval [-h,+h]
b=Box([[-4,4]])
d=Box([[-1,1]])
h=0.5
o=8 # Temporal order
f=VectorTaylorFunction.identity(b)
phi=flow(f,d,h,o)
print "phi:",phi
print

# Compute two time steps of the flow of the Taylor function f starting in domain D for the interval [h,2h]
phi0=phi
phi0h=partial_evaluate(phi,1,h)
dd=phi0h.range()
phi=flow(f,dd,h,o)
tr=ScalarTaylorFunction.variable([{0:2*h}],0)-h
phi1=compose(phi,combine(phi0h,tr))
print "phi0:",phi0
print "phi1:",phi1
print
sys.exit()



## Contractors

x1=x+x*x/2
x2=x+Interval(-1,1)
print "x1:",x1
print "x2:",x1

b=refines(x1,x2)
print "refines(x1,x2):",b
b=disjoint(x1,x2)
print "disjoint(x1,x2):",b

y=intersection(x1,x2)
print "intersection(x1,x2):",y
print

# Set up algebraic equation
d1=Box([[-1,1]])
d2=Interval(-1,1)
d=join(d1,d2)
x=ScalarTaylorFunction.variable(d,0)
y=ScalarTaylorFunction.variable(d,1)
f=x-4*y+y*y
i=VectorTaylorFunction.identity(d1)
h=ScalarTaylorFunction.constant(d1,Interval(-1,+1))

# Newton contractor to solve f(x,h(x))=0 for scalar f
h=intersection(rec(derivative(f,1).range())*compose(f,join(i,h)),h)
print "h:",phi
print


# Set up differential equation
b=Box([[-1,1],[-1,1]])
o=ScalarTaylorFunction.constant(b,1)
x=ScalarTaylorFunction.variable(b,0)
y=ScalarTaylorFunction.variable(b,1)
f=join(o,x) # [dot(x),dot(y)]=[1,x]
d=Box([[0,0.125],[0,0.125]])
h=Float(0.5)
d0=join(d,Interval(-h,+h))
x0=ScalarTaylorFunction.variable(d0,0)
y0=ScalarTaylorFunction.variable(d0,1)
t=ScalarTaylorFunction.variable(d0,2)
phi=join(x0,y0)

# Picard operator to solve dot(phi)(x,t) = f(phi(x,t))
phi=antiderivative(compose(f,phi),2)
print "phi:",phi
print



# End
