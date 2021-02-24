#!/usr/bin/python

##############################################################################
#            modules_tutorial.py
#
#  Copyright 2009-21  Pieter Collins <pieter.collins@maastrichtuniversity.nl>
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


def numeric_demonstration():
    #! [Numeric demonstration]

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
    dp=DoublePrecision()
    dp=double_precision
    mp=MultiplePrecision(128)

    # Create a raw double-precision number
    xdp=FloatDP(exact(1.75),dp)
    # Create a raw multiple-precision number
    xmp=FloatMP(exact(1.75),mp)

    # Create double-precision bounds for a value
    xdpb=FloatDPBounds(Decimal(1.2),dp) # Creates the interval [1.19999999999999996:1.20000000000000018]
    print("FloatDPBounds(1.2):",xdpb)

    # Create double-precision bounds for a range of values
    xdpb=FloatDPBounds(Rational(11,10),Rational(14,10),dp) # Creates the interval [1.09999999999999987:1.40000000000000013]
    print("FloatDPBounds(11/10,14/10,dp):",xdpb)

    # Create multiple-precision bounds for a value
    xmpb=FloatMPBounds(dec("1.2"),mp) # Creates the interval [1.19999999999999996:19999999999999996]
    print("FloatMPBounds(1.2,mp):",xmpb)

    # Create multiple-precision bounds for a range of values
    xmpb=FloatMPBounds(dec(1.5),dec(2.25),mp) # Creates the interval [1.5,2.25]
    print("FloatMPBounds(1.5,2.25,mp):",xmpb)

    # Create multiple-precision bounds for a range of values
    xmpb=FloatMPBounds(Rational(11,10),Rational(14,10),mp) # Creates the interval [1.10000000000000009:1.39999999999999991]
    print("FloatMPBounds(11/10,14/10,mp):",xmpb)

    # Create a double-precision approximation
    xdpa=FloatDPApproximation(1.23,dp)
    print("FloatDPApproximation(1.23,dp):",xdpa)

    xmpa=FloatMPApproximation(1.23,mp)
    print("FloatMPApproximation(1.23,dp):",xmpa)
    xmpa=FloatMPApproximation(Decimal("1.23"),mp)
    print("FloatMPApproximation(Decimal(\"1.23\"),dp):",xmpa)
    #! [Numeric demonstration]


def linear_algebra_demonstration():
    #! [Linear Algebra demonstration]

    # Create an interval vector
    b=FloatDPBoundsVector([1,{2:3},{exact(3.875):exact(4.125)}],dp)
    print("b:",b)

    Aq=RationalMatrix([[1,2,4],[3,1,2],[0,0,1]])
    print("Aq:",Aq)

    # Create an interval matrix
    A=FloatDPBoundsMatrix([[1,2,4],[3,exact(1.5),2],[0,0,1]],dp)
    print("A:",A)

    # Solve the linear equation Ax=b
    x=solve(A,b)
    print("A\\b:",x)
  #! [Linear Algebra demonstration]


def function_demonstration():
    #! [Function demonstration]

    # Create a user-defined scalar-valued function
    x=[EffectiveScalarMultivariateFunction.coordinate(2,index) for index in range(0,2)]
    sf=sqrt(sqr(x[0])+sqr(x[1]))
    print("sf:",sf)
    a=FloatDPApproximationVector([4,3],dp)
    sfa=sf(a)
    print("sf(a):",sfa)
    b=FloatDPBoundsVector([4,3],dp)
    sfb=sf(b)
    print("sf(b):",sfb)

    x=EffectiveVectorMultivariateFunction.identity(3)
    vf=EffectiveVectorMultivariateFunction([sqrt(sqr(x[0])+sqr(x[1])),x[1]+x[2]])
    print("vf:",vf)
    a=FloatDPApproximationVector([4,3,0],dp)
    vfa=vf(a)
    print("vf(a):",vfa)
    b=FloatDPBoundsVector([4,3,0],dp)
    vfb=evaluate(vf,b)
    print("evaluate(vf,b):",vfb)

    dsf0=derivative(sf,0);
    print("derivative(sf,0):",dsf0);

    p0=FloatDPApproximationMultivariatePolynomial.coordinate(2,0,dp);
    p1=FloatDPApproximationMultivariatePolynomial.coordinate(2,1,dp);
    print("p1:",p1)
    a3=FloatDPApproximation(3,dp)
    a5=FloatDPApproximation(5,dp)
    q=p0*(p0+p1*a3)+a5
    print("q:",q)
    #! [Function demonstration]


def calculus_demonstration():
    #! [Calculus demonstration]

    # Create a box to act as the domain of a Taylor function
    dom=BoxDomainType([{4:7},{1:6},{-1:+1}])
    print("dom:",dom)

    # Create a sweeper to control the accuracy of a Taylor function
    swp=ThresholdSweeperDP(dp,1e-8)
    #swp=GradedSweeper(dp,6)
    print("swp:",swp)

    # Create the scalar Taylor model representing a constant function with given value on the domain dom
    tc=ValidatedScalarMultivariateTaylorFunctionModelDP.constant(dom,dec(4.2),swp)
    print("tc:",tc,"\n")

    # Create the scalar Taylor model representing the function x1 on the domain dom
    t1=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom,1,swp)
    print("t1:",t1,"\n")

    # Create all Taylor variables on the domain dom
    t=[ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom,i,swp) for i in range(0,3)]
    print("t[0]:",t[0])
    print("t[1]:",t[1])
    print("t[2]:",t[2])

    # Make shorthands for the variable names
    tx=t[0]
    ty=t[1]
    tz=t[2]
    print("tx:",tx)
    print("ty:",ty)
    print("tz:",tz)

    #Create a ValidatedScalarMultivariateTaylorFunctionModelDP from a Function
    p=EffectiveScalarMultivariateFunction.coordinate(3,0)
    tp=ValidatedScalarMultivariateTaylorFunctionModelDP(dom,p,swp)
    print("p:",p)
    print("tp:",tp)

    # The domain D of tx.
    print("tx.domain():",tx.domain())
    # A not-necessarily tight over-approximation to p(D)+/-e.
    print("tx.codomain():",tx.codomain())
    # An over-approximation to p(D)+/-e.
    print("tx.range():",tx.range())

    # Define some constants
    n=int(3)
    c=FloatDPValue(exact(0.75),dp)
    b=FloatDPBounds(exact(0.625),exact(0.875),dp);

    # Arithmetic on Taylor models
    +tx; -tx; tx+ty; tx-ty; tx*ty; tx/ty;

    # Mixed arithmetic
    tx+c; tx-c; tx*c; tx/c;
    c+tx; c-tx; c*tx; c/tx;

    tx+b; tx-b; tx*b; tx/b;
    b+tx; b-tx; b*tx; b/tx;

    # Lattice functions
    min(tx,ty); max(tx,ty); abs(tx);

    # Univariate arithmetical functions
    neg(tx);
    rec(tx);
    sqr(tx);
    pow(tx,n);

    # Algebraic and transcendental functions
    sqrt(tx);
    exp(tx); log(tx);
    sin(tx), cos(tx), tan(tx/2)

    # Inplace operations
    tx+=ty; tx-=ty;
    tx+=c; tx-=c; tx*=c; tx/=c;
    tx+=b; tx-=b; tx*=b; tx/=b;

    # Set coordinate functions tx0, tx1
    dom=BoxDomainType([{-1:+1},{-1:+1}])
    tx0=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom,0,swp);
    tx1=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate(dom,1,swp);

    tf=ValidatedScalarMultivariateTaylorFunctionModelDP(1+2*tx0+3*tx0*tx1+tx0*tx0)
    # Compute an antiderivative of f with respect to x[0] whose value is zero at the midpoint of the domain of x[0]
    Itf0=antiderivative(tf,0)
    # Compute the antiderivative of f with respect to x[1] whose value is zero when x[1]=c
    Itf1=antiderivative(tf,1,c)
    print("antiderivative(tf,0):",Itf0)
    print("antiderivative(tf,1):",Itf1)

    # Restrict to a subdomain
    subdom=BoxDomainType([{exact(-0.125):exact(0.125)},{exact(0.5):exact(0.75)}])
    w=tx0*tx1
    rw=restriction(w,subdom)
    print("restriction(w,subdom):",rw)

    #Join tf and tg into a VectorTaylorFunction x:->[tf(x),tg(x)]
    tg=ValidatedVectorMultivariateTaylorFunctionModelDP([exp(tx0)*cos(tx1)])
    tfg=join(tf,tg)
    print("join(tf,tg):",tfg)

    # Combine tf and tg into a VectorTaylorFunction (x,y):->[tf(x),tg(y)]
    tfg=combine(tf,tg)
    print("combine(tf,tg):",tfg)

    # Function composition
    dom=BoxDomainType([{4:7},{1:6},{-1:1}])
    th=ValidatedVectorMultivariateTaylorFunctionModelDP.identity(dom,swp)
    cd=th.codomain()

    f=EffectiveScalarMultivariateFunction.coordinate(3,0)
    g=EffectiveVectorMultivariateFunction.zeros(2,3)

    tf=ValidatedScalarMultivariateTaylorFunctionModelDP(cd,f,swp)
    tg=ValidatedVectorMultivariateTaylorFunctionModelDP(cd,g,swp)

    # Compose an scalar multivariate function and a vector Taylor function
    fth=compose(f,th)
    print("compose(f,th):",fth)

    # Compose a vector multivariate function and a vector Taylor function
    gth=compose(g,th)
    print("compose(g,th):",gth)

    # Compose a scalar and a vector Taylor function
    assert(subset(th.codomain(),tf.domain()))
    tfh=compose(tf,th)
    print("compose(tf,th):",tfh)

    # Compose two vector Taylor functions
    assert(subset(th.codomain(),tg.domain()))
    tgh=compose(tg,th)
    print("compose(tg,th):",tgh)


    # Consistency and refinement checking
    dom=BoxDomainType([{-1:+1},{-1:+1}]);
    tf1=tx0+tx0*tx0/2
    tf2=tx0+FloatDPBounds(-1,1,dp)
    print("tf1:",tf1)
    print("tf2:",tf1)
    ic=inconsistent(tf1,tf2)
    print("inconsistent(tf1,tf2):",ic)
    rf=refines(tf1,tf2)
    print("refines(tf1,tf2):",rf)
    tf=refinement(tf1,tf2)
    print("refinement(tf1,tf2):",tf)
    #! [Calculus demonstration]



def algebraic_solver_demonstration():
    #! [Algebraic Solver demonstration]

    # Compute the solution h to the vector equation f(x,h(x))=0
    dom=BoxDomainType([{-1:+1},{-1:+1},{-1:+1}])
    x=EffectiveScalarMultivariateFunction.coordinate(3,0)
    y0=EffectiveScalarMultivariateFunction.coordinate(3,1)
    y1=EffectiveScalarMultivariateFunction.coordinate(3,2)
    f=join(x+4*y0+y1,y0+y1)
    slv=IntervalNewtonSolver(1e-8,6)
    print("f:",f)
    hf=slv.implicit(f,BoxDomainType([{-1:+1}]),BoxDomainType([{-1:+1},{-1:+1}]))
    print("implicit(f):",hf)

    # Compute the solution h to the scalar equation g(x,h(x))=0
    # with f(x,y)=4+x-y^2, so y=sqrt(4+x)
    dom=BoxDomainType([{-1:+1},{-1:+1}])
    x=EffectiveScalarMultivariateFunction.coordinate(2,0)
    y=EffectiveScalarMultivariateFunction.coordinate(2,1)
    g=x-4*y+y*y
    print("g:",g)
    hg=slv.implicit(g, BoxDomainType([{-1:+1}]),IntervalDomainType({-1:+1}))
    print("implicit(g):",hg)
    #! [Algebraic Solver demonstration]


def differential_solver_demonstration():
    #! [Differential Solver demonstration]

    # Compute the flow of the Taylor function f starting in the domain dom for time interval [-h,+h]
    dom=BoxDomainType([{-1:+1}])
    bbx=BoxDomainType([{-4:+4}])
    h=1/two
    f=ValidatedVectorMultivariateFunction.identity(1)
    integrator=GradedTaylorSeriesIntegrator(1e-8)

    print(f,dom,h)
    phis=integrator.flow(f,dom,h)
    assert(len(phis)==1)
    phi=phis[0]
    print("phi:",phi,type(phi))

    # Compute two time steps of the flow of the Taylor function f starting in domain D for the interval [h,2h]
    phi0=phi
    print("phi.domain():",phi.domain())
    print("h:",h)
    phi0h=partial_evaluate(phi,1,FloatDPBounds(h,dp))
    dom1=phi0h.codomain()
    phi=integrator.flow(f,dom1,h)[0]
    print("phi:",phi)
    swp=GradedSweeperDP(dp,6);
    tr=ValidatedScalarMultivariateTaylorFunctionModelDP.coordinate([{0:2*h}],0,swp)-h
    tr=ValidatedScalarMultivariateFunctionModelDP(tr)
    phi1=compose(phi,combine(phi0h,tr))
    print("phi0:",phi0)
    print("phi1:",phi1)
    #! [Differential Solver demonstration]


if __name__=='__main__':
    numeric_demonstration()
    print()
    linear_algebra_demonstration()
    print()
    function_demonstration()
    print()
    calculus_demonstration()
    print()
    algebraic_solver_demonstration()
    print()
    differential_solver_demonstration()
    print()

