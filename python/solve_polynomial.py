#!/usr/bin/python

##############################################################################
#            solve_polynomial.py
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

"""A module for finding all positive solutions of a system of polynomial equations using interval solvers.

Example usage:
from solve_polynomial import *          # imports module and Ariadne supporting types
[x,y]=variables(2)                      # sets x,y to be the coordinate variables in R2
p=PolynomialFunction([3-x-y*y,1.5-x])   # defines a polynomial function p
e=1e-12                                 # sets the error bound e for the solver
s=positive_solutions(p,e)               # computes the set of positive solutions s
for y in s: print y,p(y)                # prints the solutions and evaluates p at them

"""

from ariadne import Box,Polynomial,PolynomialFunction,TaylorFunction,KrawczykSolver,variables
from sys import exit
import __builtin__

Polynomial.__repr__=Polynomial.__str__
TaylorFunction.__repr__=TaylorFunction.__str__
#ExpressionInterface.__repr__=ExpressionInterface.__str__

def min_zero(n):
    """Compute the minumum binary digit of n which is equal to zero."""
    k=0
    while n & __builtin__.pow(2,k):
        k+=1
    return k


def flip_coefficients(p,k):
    """The polynomial computed from p by setting the power of x[k] in each term to d[k]-a[k], where d[k] is the degree of p in variable x[k], and a[k] is the degree of x[k] in the current term.

    The roots of flip_coefficients(p,k) are the same as those of p, except with x[k] replaced by 1/x[k].
    """
    assert(p.argument_size()>k)
    q=Polynomial(p.argument_size())

    # Compute the degree of p in coefficient x[k]
    dk=0
    for m in p:
        dk=__builtin__.max(dk,m.key()[k])

    # Compute the new polynomial q
    for m in p:
        a=m.key()
        c=m.data()
        a[k]=dk-a[k]
        q.insert(a,c)

    return q


def flip_point(x,j):
    """Sets x[i] to 1/x[i] whenever the ith bit of j is equal to 1."""
    for i in range(0,len(x)):
        if j & __builtin__.pow(2,i):
            x[i]=1.0/x[i]


def flips(n):
    i=0
    j=0
    while i!=__builtin__.pow(2,n)-1:
        k=min_zero(i)
        j=j ^ __builtin__.pow(2,k)
        i+=1
        print i,j


def positive_solutions(f,e):
    """Finds all positive solutions of the polynomial function f to accuracy e using an interval solver."""
    solver=KrawczykSolver(1e-8,12)

    pf=PolynomialFunction(f)
    res_sz=pf.result_size()
    arg_sz=pf.argument_size()

    i=0 # A counter which controls the iteration
    q=0 # The current actual quadrant
    unit_box=Box([{0:1}]*res_sz)
    solutions=[]
    new_solutions=solver.solve(pf,unit_box)
    #print j,pf,new_solutions,new_solutions
    solutions+=new_solutions
    while i!=__builtin__.pow(2,arg_sz)-1:
        k=min_zero(i)
        q=q ^ __builtin__.pow(2,k)
        i+=1
        for l in range(0,res_sz):
            pf[l]=flip_coefficients(pf[l],k)
        new_solutions=solver.solve(pf,unit_box)
        #print j,pf,new_solutions,
        for solution in new_solutions:
            flip_point(solution,q)
        #print new_solutions
        solutions+=new_solutions

    return solutions


def parameterised_solution(f,ia,ix,e):
    """Finds all positive solutions of the equation f(a,x)=0 in the parameter range ia within ix to error e using an interval solver."""
    solver=KrawczykSolver(e,12)
    return solver.implicit(f,ia,ix)



if __name__=="__main__":
    print __doc__

    [a,x]=variables(2)
    f=PolynomialFunction([1-5*x-a*x*x])
    pd=Box([{-0.25:0.25}])
    xd=Box([{0:1}])
    e=1e-6
    solver=KrawczykSolver(e,4)
    print solver.implicit(f,pd,xd)
    #print positive_parameterised_solutions(f,pd,e)
