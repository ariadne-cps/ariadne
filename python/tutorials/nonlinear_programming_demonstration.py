#!/usr/bin/python3

# nonlinear_programming_demonstration.py
# Copyright  2023  Pieter Collins
#
# This file is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This file is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTAXILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file. If not, see <https://www.gnu.org/licenses/>.


from pyariadne import *

# Set logging verbosity
Logger.instance().configuration().set_verbosity(0)

# Available approximate optimisers
opt = InteriorPointOptimiser()
opt = InfeasibleInteriorPointOptimiser()
opt = SplitInfeasibleInteriorPointOptimiser()


# Constructors for numbers and vectors
y=ValidatedNumber({"0.625":"0.7"})
print("y:",y)
x=FloatDPBounds({"0.625":"0.7"},dp)
print("x:",x)
v=Vector[FloatDPBounds]([{"0.625":"0.7"},("-1.5","-1.25")],dp)
print("v:",v)

#Constructors for polynomials
p=MultivariatePolynomial[FloatDPBounds].coordinates(2,dp)
q=2*p[0]*(p[0]-p[1]-3)
print("q:",q)
print("q**2:",q**2)
q=MultivariatePolynomial[FloatDPBounds]({ (1,1):-2, (2,0):2, (1,0):-6},dp)
print("q:",q)
print("q[(1,1)]:",q[(1,1)])
a=MultiIndex((1,1))
a=MultiIndex([1,1])
print("a:",a)
print("q[a]:",q[a])
print("q.coefficient(a):",q.coefficient(a))
print("q.terms():",[(term.index(),term.coefficient()) for term in q])
v=Vector[FloatDPBounds]([3,5],dp)
print("q(v):",q(v))
print("evaluate(q,v):",evaluate(q,v))
print("compose(q,p):",compose(q,p))
print("derivative(q,0):",derivative(q,0))
print("derivative(q,1):",derivative(q,1))
print("antiderivative(q,0):",antiderivative(q,0))
print("antiderivative(q,1):",antiderivative(q,1))

def to_function(p):
    id=ValidatedVectorMultivariateFunction.identity(2)
    first_term=True
    for term in q:
        a=term.index()
        c=term.coefficient()
        t=c * ( (id[0]**a[0]) * (id[1]**a[1]) )
        if first_term:
            first_term=False
            f=t
        else:
            f=f+t
    return f
print("to_function(q):",to_function(q))
print()


# Construct functions and Boxes
g=Function(2, lambda x : [sqr(x[0])+2*sqr(x[1])-1])
print("g:",g,type(g))
f=Function(2, lambda x :  2*x[0]+3*x[1])
print("f:",f,type(f))
D=ExactBoxType([{-2:2},(-1,1)])
print("D:",D,type(D))
C=ExactBoxType([{0:0}])
C=ExactBoxType([{"-0.0625":"0.0625"}])
print("C:",C,type(C))

opt = InfeasibleInteriorPointOptimiser()
# Simple feasible case
D1=ExactBoxType([{-1:2},(-1,1)])
fp1 = ValidatedFeasibilityProblem(D1,g,C)
print("fp1:",fp1)
print("opt.feasible(D1,g,C):",opt.feasible(fp1))
# Construct using constraints
g0=g[0]
C0l=ExactNumber(C[0].lower_bound())
C0u=ExactNumber(C[0].upper_bound())
gC0=C0l<=(g0<=C0u)
print("gC0:",gC0,type(gC0))
fp1 = ValidatedFeasibilityProblem(D1,[C0l<=(g0<=C0u)])
fp1 = ValidatedFeasibilityProblem(D1,[(C0l<=g0)<=C0u])
# Note that we cannot construct a constraint using Cl <= g <= Cu, need (Cl <= g) <= Cu or Cl <= (g <= Cu)
# We can also construct constraints C0l<=g0, g0>=C0u, g0<=C0u, g0=C0
print("fp1:",fp1)
print("opt.feasible(D1,g,C):",opt.feasible(fp1))

# Feasible case
D2=ExactBoxType([{"0.75":2},("0.375",1)])
print("opt.feasible(D2,g,C):",opt.feasible(ValidatedFeasibilityProblem(D2,g,C)))

# Boundary case; only feasibly point is (0.75,0.5)
D3=ExactBoxType([{"0.75":2},("0.5",1)])
print("opt.feasible(D3,g,C):",opt.feasible(ValidatedFeasibilityProblem(D3,g,C)))

# infeasible case
D4=ExactBoxType([{"0.5":2},("0.75",1)])
print("opt.feasible(D4,g,C):",opt.feasible(ValidatedFeasibilityProblem(D4,g,C)))

print()


# Check feasibility
fc = FeasibilityChecker()
X=Vector[ValidatedNumber]([{"0.625":"0.75"},{"0.5":"0.625"}])
print("X:",X)
print("fc.validate_feasibility(g,X):",fc.validate_feasibility(g,X))
print()


# Solve an optimisation problem
print("opt:",opt)
p = ApproximateOptimisationProblem(f,D,g,C)
print("p:",p)
# Initial values of primal and dual variables
x0=Vector[FloatDPApproximation]([".625","0.7"],dp)
y0=Vector[FloatDPApproximation]([".125"],dp)
print("x0:",x0," y0:",y0)

# Examine data used for a step of the solver
step_data = opt.initial_step_data(p)
print("step_data",step_data)

# Compute approximate local optimum
vxy_opt = opt.minimise_hotstarted(p,x0,y0)
print("optimal_vxy:",vxy_opt,type(vxy_opt))
v_opt = vxy_opt.value()
x_opt = vxy_opt.primal()
y_opt = vxy_opt.dual()
xa=Vector[FloatDPApproximation](x_opt,dp)
print("f(x):",f(xa),"  g(x):",g(xa))
print()

# Rigorous solvers
opt=InfeasibleKarushKuhnTuckerOptimiser()
p = ValidatedOptimisationProblem(f,D,g,C)
print("opt:",opt)
print("p:",p)
vxy_opt = opt.minimise(p)
print("opt.minimise(p):",vxy_opt,type(vxy_opt))
print()



# Interval Newton approach for unconstrained optimisation

print("Hand-coded interval Newton solver")
f=Function(2, lambda x : (((x[0]+3)*x[0]-2)*x[0]+1)*x[0]-5*x[0]*x[1]+3*x[1]*x[1])
g=derivatives(f)
print("f:",f)
print("f':",g)

tolerance=Decimal("0.000001")

def error(X):
    e=X[0].error()
    for i in range(len(X)):
        e=max(e,X[i].error())
    # return e
    return Dyadic(cast_exact(FloatDPApproximation(e))) # Conversion should be unnecessary here


B=ExactBoxType([("0.5","0.75"),("0.25","0.5")])
print("B",B,type(B))
X=cast_singleton(B)
is_refinement=False
while possibly(error(X)>tolerance):
    x=cast_exact(X) # The midpoint of X
    ddfX=hessian(f,X)
    dfx=transpose(gradient(f,x))
    try:
        dX=lu_solve(ddfX,dfx)
    except:
        print("Cannot verify critical point: singular matrix")
        break
    nX=x-dX
    print("X:",X," x:",x," nX",nX)
    if inconsistent(nX,X):
        print("No critical point in B")
        break
    elif (refines(X,nX)):
        print("No refinement possible: inconclusive")
        break
    elif (refines(nX,X)):
        print("Refinement found! Unique critical point in B.")
        is_refinement=True
        X=nX
    else:
        X=refinement(X,nX)
if is_refinement:
    print("Verified unique critical point in B at",X)
else:
    print("Could not verifiy critical point in B")
print()

# Solve using a builtin Ariadne solver
tolerance=1e-6
max_steps=32
algslv=IntervalNewtonSolver(tolerance, max_steps)

B=ExactBoxType([("0.0","1.0"),("0.0","1.0")])
bxs=algslv.solve_all(g,B)
print(bxs)

