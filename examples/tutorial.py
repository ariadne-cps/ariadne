#!/usr/bin/python

##############################################################################
#            tutorial.py
#
#  Copyright 2009  Pieter Collins <Pieter.Collins@cwi.nl>
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

# Import all classes in the ariadne module
from ariadne import *

# Set Ariadne name of builtin float type
Float=float

def box(x): return IVector(x)

# Print a list of all available classes and functions
print dir(),"\n"

# Print a list of all available classes and functions
print dir()

# Create an integer
z=Integer(8)
print z

# Create a rational
q=Rational(8,6)
print q

# Create a float
x=Float(1.25)
print x

# Create an interval
i=Interval(1.5,2.25) # Creates the interval [1.5,2.25]
print "i:",i

# Create an interval
i=Interval(1.2) # Creates the interval [1.19999999999999996:19999999999999996]
print "i:",i

# Create an interval
i=Interval(1.1,1.4) # Creates the interval [1.10000000000000009:1.39999999999999991]

print "i:",i

# Create an interval
i=Interval(Rational(12,10)) # Creates the interval [1.19999999999999996:1.20000000000000018]
print "i:",i

# Create an interval
i=Interval(Rational(11,10),Rational(14,10)) # Creates the interval [1.09999999999999987:1.40000000000000013]
print "i:",i


# Create an interval vector
v=IVector([1,[2,3],4])
print "v:",v

# Create an interval matrix
A=IMatrix([[1,[2,3],4],[1.5,1.1,2]])
print "A:",A

# Solve the linear equation Ax=b
solve(A,b)


print "\n\n\n\n\n"



argument_size=2
function=lambda x:sqrt(sqr(x[0])+sqr(x[1]))
f=UserExpression(argument_size,function)




# Create a Polynomial expression in three unknowns, with value \f$x_0\f$.
#px=PolynomialExpression.variable(3,0)

# Make a shorthand for constructing Polynomial expressions
def p(n,j):
    return PolynomialExpression.variable(n,j)

# Arithmetic for Polynomial expressions
#+p; -p; p+q; p-q; p*q;

#p+c; p-c; p*c; p/c;
#c+p; c-p; c*p;

#p+i; p-i; p*i; p/i;
#i+p; i-p; i*p;

#p+=q; p-=q;
#p+=c; p-=c; p*=c; p/=c;
#p+=i; p-=i; p*=i; p/=i;

print "\n\n\n\n\n"





# Create a box to act as the domain of a Taylor expression
d=IVector([[4,7],[1,6],[-1,1]])
print "d:",d

# Create the TaylorExpression representing the function x1 on the domain d
t=TaylorExpression.variable(d,1)
print "t:",t

# Create all TaylorExpression variables on the domain dom
tv=TaylorExpression.variables(d)
print "tv[0]:",tv[0]
print "tv[1]:",tv[1]
print "tv[2]:",tv[2]

# Make shorthands for the variable names
x=tv[0]
y=tv[1]
z=tv[2]
print "x:",x,"\ny:",y,"\nz:",z

# Make a shorthand for constructing Taylor expressions
def t(d,j):
    return TaylorExpression.variable(d,j)

#Create a TaylorExpression from a PolynomialExpression
t=TaylorExpression(p,D)

# The domain D of x.
x.domain()
# A not-necessarily tight over-approximation to p(D)+/-e.
x.codomain()
# An over-approximation to p(D)+/-e.
x.range()
# The function m(x) defined by m(s(x)) = p(x), where s is the scaling function from the unit box to the domain D.
x.model()
# Convert to a polynomial expression.
x.polynomial()

# Arithmetic on Taylor expressions
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

# Reset x
x=v[0]; 
x=v[1]; 

# Comparison functions
min(x,y); max(x,y); abs(x);

# Compound arithmetical functions
neg(x); rec(x); sqr(x); pow(x,n);

# Algebraic and transcendental functions
sqrt(x); 
exp(x); log(x);
sin(x), cos(x), tan(x/100)

x=embed(x,box([-1,1]))

#Join x and y into a TaylorFunction
f=join(x,y)
print "join(x,y):",f

x=combine(x,y)
print x.domain()

# Compose a TaylorExpression and a TaylorFunction
print compose(x,g)

# Compose two TaylorFunctions and a TaylorFunction
compose(join(x,y),g)

# Compose an polynomial expression and a TaylorFunction
print compose(p,f)


gd=box([[0,10],[1,7]])
print f.codomain(),gd


g=t(gd,0)*t(gd,1)
print compose(g,f)

assert( subset(f.codomain(),gd) )

print compose(p(2,0)*p(2,1),join(t(3,0),t(3,1)))

# Create a TaylorFunction equal to the identity on dom
f=TaylorFunction.identity(d)
print "f:",f

implicit(f)

implicit(x)

antiderivative(f,j)

derivative(midpoint(f),j)

derivative(midpoint(f),[j0,j1])

flow(f,D,h,B)
flow(tf,D,h)

## Contractors

refines(x0,x1)
subset(x0,x1)
intersects(x0,x1)


intersection(x0,x1)


# Newton contractor to solve f(x,h(x))=0
f=(4+tx)*tx
h=inverse(derivative2(f))*compose(f,join(i,h))

# Picard operator to solve dot(phi)(x,t) = f(phi(x,t))
f=tx
phi=antiderivative(compose(f,phi),i)



# End
