#!/usr/bin/python
from ariadne import FloatMatrix, FloatVector, IntervalVector,InteriorPointSolver,transpose,inverse,midpoint
#from scipy import matrix as Matrix
#from scipy import vector as Vector
Vector=FloatVector
Matrix=FloatMatrix

def dot(v1,v2):
    s=0.0
    for i in range(len(v1)):
        s=s+v1[i]*v2[i]
    return s

def emin(v):
    r=v[0]
    for i in range(len(v)):
        r=min(r,v[0])
    return r

def esub(v,s):
    return Vector( [v[i]-s for i in range(len(v)) ] )

def emul(v1,v2):
    return Vector( [v1[i]*v2[i] for i in range(len(v1)) ] )

def ediv(v1,v2):
    return Vector( [v1[i]/v2[i] for i in range(len(v1)) ] )

def erec(v):
    return Vector( [1.0/v[i] for i in range(len(v)) ] )

def adat(A,D):
    S=FloatMatrix(A.row_size(),A.row_size())
    for i in range(A.row_size()):
        for j in range(A.row_size()):
            S[i,j]=0.0;
            for k in range(A.column_size()):
                S[i,j]+=A[i,k]*D[k]*A[j,k]
    return S

def damul(D,A):
    R=A*0.0
    for i in range(A.row_size()):
        for j in range(A.column_size()):
            R[i,j]=D[i]*A[i,j]
    return S

def admul(A,D):
    R=A*0.0
    for i in range(A.row_size()):
        for j in range(A.column_size()):
            R[i,j]=A[i,j]*D[j]
    return S


# The standard linear programming problem max cx; Ax = b, x>=0.
class StandardLinearProgrammingProblem:
    def __init__(self, A,b,c):
        self.A=A
        self.b=b
        self.c=c
        assert(A.row_size()==b.size())
        assert(A.column_size()==c.size())
        assert(A.column_size()>A.row_size()>0)


# The constrained feasibility problem Ax <= c, l<=x<=u
class FeasibilityProblem:
    def __init__(self, A,c,l,u):
        self.A=A
        self.c=c
        self.l=l
        self.u=u
        assert(A.row_size()==c.size())
        assert(A.column_size()==l.size())
        assert(A.column_size()==u.size())


class PyInteriorPointSolver:
    def solve(self, problem):
        assert(isinstance(problem,StandardLinearProgrammingProblem))
        self.solve_lp_problem(problem.A,problem.b,problem.c)

    def lp_step(self,A,b,c,x,y,z):
        m=len(b); n=len(c)
        gamma=1.0/1024
        sigma=1.0/8

        mu=dot(x,z)/n

        rb=A*x-b
        rc=At*y+z-c
        rs=esub(emul(x,z),sigma*mu)

        ry=A*(ediv(rs-emul(x,rc),z))-rb

        S=adat(A,ediv(x,z))
        Sinv=inverse(S)

        dy=Sinv*ry
        dz=-rc-dy*A
        dx=-ediv((rs+emul(x,dz)),z)

        allpositive=False
        a=1.0/scale
        while not allpositive:
            a=a*scale
            nx=x+a*dx
            nz=z+a*dz
            allpositive = emin(nx)>0.0 and emin(nz)>0.0 and emin(emul(nx,nz))>gamma*mu
        ny=midpoint(y+a*dy)

        x=nx; y=ny; z=nz


    def solve_lp_problem(self,A,b,c,x,y,z):
        max_error=1e-8
        max_steps=24
        # Primal variables x, dual variables y, slack variables z
        # Ax=b yA+z=c
        cx=dot(c,x)
        yb=dot(y,b)

        print "A:",A," b:",b," c:",c
        print "x:",x," y:",y," z:",z
        print "cx:",cx,"yb:",yb

        steps=0
        while(steps<max_steps and cx-yb<max_error):
            lp_step(A,b,c,x,y,z)
            cx=dot(c,x)
            yb=dot(y,b)
            steps+=1
        return (Interval(yb,cx),x,y)

    def check_feasibility(self,A,b,c,x,y):
        iA=IntervalMatrix(A)
        iAT=transpose(iA)
        ib=IntervalVector(b)
        ic=IntervalVector(c)
        ix=IntervalVector(x)
        iy=IntervalVector(y)
        iS=iA*iAT
        iSinv=inverse(iS)
        ix+=iAT*(iSinv*(ib-iA*ix))
        iz=ic-iAT*iy
        return emin(x).lower()>0.0 and emin(z).lower()>0.0

if __name__=='__main__':
    A=Matrix([[1,3,1,0],[-2,1,0,1]])
    b=Vector([2,1])
    c=Vector([-1,-1,0,0])
    x0=Vector([0.125,0.125,1.5,1.125])
    y0=Vector([-10.0,-2.5])
    z0=midpoint(c-transpose(A)*y0)
    (x,y,z)=InteriorPointSolver().optimize(A,b,c,x0,y0,z0)

    print x,y,z
