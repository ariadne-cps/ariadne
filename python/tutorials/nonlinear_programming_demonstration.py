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


g=Function(2, lambda x : [sqr(x[0])+2*sqr(x[1])-1])

f=Function(2, lambda x :  2*x[0]+3*x[1])
D=BoxDomainType([{-2:2},(-1,1)])
C=BoxDomainType([{0:0}])
C=BoxDomainType([{"-0.0625":"0.0625"}])

opt = InfeasibleInteriorPointOptimiser()
# Simple feasible case
D1=BoxDomainType([{-1:2},(-1,1)])
fp1 = ValidatedFeasibilityProblem(D1,g,C)
print("fp1:",fp1)
print("opt.feasible(D1,g,C):",opt.feasible(fp1))

# Feasible case
D2=BoxDomainType([{"0.75":2},("0.375",1)])
print("opt.feasible(D2,g,C):",opt.feasible(ValidatedFeasibilityProblem(D2,g,C)))

# Boundary case; only feasibly point is (0.75,0.5)
D3=BoxDomainType([{"0.75":2},("0.5",1)])
print("opt.feasible(D3,g,C):",opt.feasible(ValidatedFeasibilityProblem(D3,g,C)))

# infeasible case
D4=BoxDomainType([{"0.5":2},("0.75",1)])
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
v_opt = vxy_opt[0]
x_opt = vxy_opt[1]
y_opt = vxy_opt[2]
print("optimal v:",v_opt," x:",x_opt," y:",y_opt)
xa=Vector[FloatDPApproximation](x_opt,dp)
print("f(x):",f(xa),"  g(x):",g(xa))
print()

# Rigorous solvers
opt=InfeasibleKarushKuhnTuckerOptimiser()
p = ValidatedOptimisationProblem(f,D,g,C)
vxy_opt = opt.minimise(p)
print(vxy_opt)



# Interval Newton approach for unconstrained optimisation

f=Function(2, lambda x : (((x[0]+3)*x[0]-2)*x[0]+1)*x[0]-5*x[0]*x[1]+3*x[1]*x[1])
print(f)
g=derivatives(f)
print(g)

tolerance=Decimal("0.000001")

def error(X):
    e=X[0].error()
    for i in range(len(X)):
        e=max(e,X[i].error())
    # return e
    return Dyadic(cast_exact(FloatDPApproximation(e))) # Conversion should be unnecessary here


B=BoxDomainType([("0.5","0.75"),("0.25","0.5")])
print("B",B,type(B))
X=cast_singleton(B)
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
        X=nX
    else:
        X=refinement(X,nX)
print()

# Solve using a builtin Ariadne solver
tolerance=1e-6
max_steps=32
algslv=IntervalNewtonSolver(tolerance, max_steps)

B=BoxDomainType([("0.0","1.0"),("0.0","1.0")])
bxs=algslv.solve_all(g,B)
print(bxs)

