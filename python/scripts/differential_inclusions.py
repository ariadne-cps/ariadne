#!/usr/bin/python
#EVALUATING AT BOTH PARAMETERS
# Import all classes in the ariadne module
from numpy import *
import sys
import time
from ariadne import *

pi=pi()

def join_list(l):
    r=l[0]
    for i in range(1,len(l)):
        r=join(r,l[i])
    return r

def combine_list(l):
    r=l[0]
    for i in range(1,len(l)):
        r=combine(r,l[i])
    return r


def lohner_approximation(f):
    n=f.result_size()
    models=f.models()
    b=FloatVector(n)
    e=FloatVector(n)
    A=FloatMatrix(n,models[0].argument_size())
    for i in range(n):
        b[i]=models[i].value()
        for j in range(models[0].argument_size()):
            A[i,j]=models[i].gradient(j)
        e[i]=models[i].error()
    z=Zonotope(b,A,e)
    #print z.error()
    z=orthogonal_approximation(z)
    #print "z:",z
    p=polytope(z)
    #print "p:",p
    return p

    
class DifferentialInclusionIntegratorBase:
    def __init__(self, sweeper=ThresholdSweeper(1e-8), step_size=0.125, simplification_steps=128, number_of_variables_to_keep=24):
        self.sweeper=sweeper
        self.step_size=Float(step_size)
        self.number_of_steps_between_simplifications=simplification_steps
        self.number_of_variables_to_keep=number_of_variables_to_keep

    def expand_errors(self,Phi):
        domain=Phi.domain()
        errors=Phi.errors()*Interval(-1,+1)
        print "Uniform errors:",errors
        for i in range(Phi.result_size()): Phi[i].set_error(0)
        error_function=VectorTaylorFunction.identity(errors,self.sweeper)
        return embed(Phi,errors)+embed(domain,error_function)

    def solve(self,f,V,X0,T):
        """Solve the differential inclusion dot(x) in f(x)+V for x(0) in X0 up to time T."""
        print "\nf:",f,"\nV:",V,"\nX0:",X0,"\nT:",T

        # Ensure all arguments have the correct type
        f=VectorFunction(f)
        V=IntervalVector(V)
        X0=Box(X0)
        T=Real(T)

        n=X0.size()
        assert(f.result_size()==n)
        assert(f.argument_size()==n)
        assert(V.size()==n)

        h=self.step_size

        evolve_function = VectorTaylorFunction.identity(X0,self.sweeper)
        t=Float(0.0)

        evolve_sets = [IntervalConstrainedImageSet(evolve_function.domain(),evolve_function.function())]
        reach_sets = []

        step = 0
        while t<Float(T):
            print "\nt:",t,"h:",h
            if t+h>Float(T): h=Float(T)-t

            D = evolve_function.codomain()
            B=self.compute_bound(f,V,D,h)

            Phi = self.compute_step(f,V,D,h,B)
            p = Phi.argument_size()-(2*n+1)
            
            E=evolve_function
            W=VectorTaylorFunction.identity(Phi.domain()[n:n+p],self.sweeper)
            U=VectorTaylorFunction.identity(Phi.domain()[n+p:2*n+p],self.sweeper)
            H=VectorTaylorFunction.identity([Interval(0.0,h)],self.sweeper)

            Psi = partial_evaluate(Phi,2*n+p,Interval(h))
            Phi.set_sweeper(ThresholdSweeper(4e-3))
            Phi.sweep()

            reach_function=compose(Phi,combine_list([E,W,U,H]))
            reach_function=self.expand_errors(reach_function)
            evolve_function=compose(Psi,combine_list([E,W,U]))

            step+=1
            if step%self.number_of_steps_between_simplifications==0:
                evolve_function=self.simplify(evolve_function)

            print "evolve_function.range():",evolve_function.range()
            t=t+h
            reach_sets.append(lohner_approximation(reach_function))
            #reach_sets.append(IntervalConstrainedImageSet(reach_function.domain(),reach_function.function()))

            print evolve_function
            evolve_sets.append(lohner_approximation(evolve_function))
            e=evolve_function.errors()
            #evolve_sets.append(IntervalConstrainedImageSet(evolve_function.domain(),evolve_function.function()))

        return (evolve_sets,reach_sets)

    def simplify(self,phi):
        assert isinstance(phi,VectorTaylorFunction)

        m=phi.argument_size()
        n=phi.result_size()
        #Compute effect of error terms, but not of original variables
        C=[ [0.0]*n for i in range(m) ]
        e=[]
        for i in range(n):
            p=phi[i].model().expansion()
            for (a,c) in p.items():
                for j in range(m):
                    if a[j]!=0:
                        C[j][i] = C[j][i]+abs(c)

        for i in range(n):
            e.append(phi[i].error())
        print "C",C,"\ne",e

        Ce = [0]*m
        for j in range(m):
            for i in range(n):
                Ce[j] += C[j][i]
        print "\nCe:",Ce
        SCe=list(Ce)
        SCe.sort()
        print "\nSortCe:",SCe
        keep_indices = []
        remove_indices = []
        number_of_variables_to_keep=self.number_of_variables_to_keep
        if m<self.number_of_variables_to_keep: number_of_variables_to_keep=m
        threshold = (SCe[-number_of_variables_to_keep]+SCe[1-number_of_variables_to_keep])/2
        for j in range(m):
            if Ce[j] < threshold:
                remove_indices.append(j)
            else:
                keep_indices.append(j)
        print keep_indices
        print remove_indices

        old_domain=phi.domain()
        new_domain=join_list([ old_domain[j] for j in keep_indices ])
        projection=VectorTaylorFunction(m,new_domain,self.sweeper)
        for i in range(new_domain.size()): projection[keep_indices[i]]=ScalarTaylorFunction.coordinate(new_domain,i,self.sweeper)
        for i in range(len(remove_indices)): projection[remove_indices[i]]=ScalarTaylorFunction.constant(new_domain,old_domain[remove_indices[i]],self.sweeper)
        phi=compose(phi,projection)

        phi = self.expand_errors(phi)
        return phi

    def compute_bound(self,f,V,D,h):
        """Compute a bound B for the differential inclusion dot(x) in f(x)+V for x(0) in D for step size h"""
        print "D:",D,
        B = D + Interval(0,2*h)*(f(D)+V)

        if subset(D+h*(f(B)+V),B):
            for i in range(4):
                B=D+h*(f(B)+V)
        else:
            for i in range(16):
                B=D+h*(f(B)+V)

        print "B:",B
        return B



class DifferentialInclusionIntegrator3rdOrder(DifferentialInclusionIntegratorBase):
    def compute_norms(self,f,B):
        """Compute the norms K=|f(B)|, L=|Df(B)|, H=|D2f(B)| and LN=l(Df(B))"""
        # Estimate error terms
        Df=f.differentials(B,2)
        print "  Df=",Df
        (K,L,H,LN)=(Float(0),Float(0),Float(0),Float(0))
        for i in range(f.result_size()):
            Dfi=Df[i].expansion()
            (Ki,Li,Hi,LNi)=(Float(0),Float(0),Float(0),Float(0))
            for (a,c) in Dfi.items():
                if a.degree()==0:
                    Ki += mag(c)
                elif a.degree()==1:
                    Li += mag(c)
                    if a[i]==1: LNi += c.upper()
                    else: LNi += mag(c)
                else:
                    assert(a.degree()==2)
                    Hi += mag(c)
            K=max(K,Ki); L=max(L,Li); H=max(H,Hi); LN=max(LN,LNi)
            return (K,L,H,LN)

    def compute_step(self,f,V,D,h,B):
        """Compute a time-step for the differential inclusion dot(x) in f(x)+V for x(0) in D assuming bound B"""

        n=D.size()
        (K,L,H,LN)=self.compute_norms(f,B)
        KV=mag(norm(V))
        if LN>0.0: eLN=(exp(LN*h)-1)/(LN*h)
        else: eLN = Float(0.0)
        
        e = (7.0/48*KV*H*(K+KV)+7.0/8*KV*(L*L+H*(K+5*KV/2))*eLN)/(1-h*L/2)*pow(h,3)
        print "e:",e,"h:",h,"e/h^3:",e/pow(h,3)
        print "  K:",K,"KV:",KV,"L:",L,"LN:",LN,"eLN:",eLN,"H:",H,"e:",e

        E=e*IntervalVector(n,Interval(-1,+1))

        # Set up estimate for differential equation
        # Use affine wi=ai0+(t-h/2)ai1/h
        #   with |ai0|<=Vi, |ai1|<=3Vi, |wi(t)|<=5Vi/2, and |w'(t)|<=3Vi/h.
        # Alternatively use step inputs wi=ai,0 for t<h/2 and wi=ai,1 for t>h/2
        #   with |ai0|,|ai1|<= 2Vi , and |wi(t)|<=2Vi.

        # The flow is a function of n state variables, 2n parameter variables and 1 time variable
        # Assume for now noise in all variables
        # We have dot[xi](t)=fi(x(t))+wi(t)
        swp=self.sweeper
        DVH=join_list([D,V,3*V,E,Interval(-h,+h)])
        #print "DVH["+str(DVH.size())+"]:",DVH
        z=ScalarTaylorFunction(DVH,swp)
        xat=VectorTaylorFunction.identity(DVH,swp)
        w=VectorTaylorFunction(DVH,VectorFunction(n,4*n+1),swp)
        for i in range(n): w[i]=xat[n+i]+xat[2*n+i]*(xat[4*n]-h/2)/h
        #print "w:",w

        x0=VectorTaylorFunction(n,DVH,swp)
        for i in range(n): x0[i]=xat[i]
        #print "x0:",x0

        phi=VectorTaylorFunction(n,DVH,swp)
        for i in range(n): phi[i]=ScalarTaylorFunction(DVH,swp)+B[i]
        #print "phi0:",phi

        for i in range(6):
            phi=antiderivative(compose(f,phi)+w,4*n)+x0
            #print i,phi.errors()
        print (derivative(phi,4*n)-(compose(f,phi)+w)).range()

        for i in range(n):
            phi[i]=phi[i]+xat[3*n+i]

        return phi



class DifferentialInclusionIntegrator2ndOrder(DifferentialInclusionIntegratorBase):
    def compute_norms(self,f,B):
        """Compute the norms K=|f(B)|, L=|Df(B)|, and LN=l(Df(B))"""
        # Estimate error terms
        Df=f.differentials(B,1)
        print "  Df=",Df
        (K,L,LN)=(Float(0),Float(0),Float(0))
        for i in range(f.result_size()):
            Dfi=Df[i].expansion()
            (Ki,Li,LNi)=(Float(0),Float(0),Float(0))
            for (a,c) in Dfi.items():
                if a.degree()==0:
                    Ki += mag(c)
                elif a.degree()==1:
                    Li += mag(c)
                    if a[i]==1: LNi += c.upper()
                    else: LNi += mag(c)
            K=max(K,Ki); L=max(L,Li); LN=max(LN,LNi)
            return (K,L,LN)

    def compute_step(self,f,V,D,h,B):
        """Compute a time-step for the differential inclusion dot(x) in f(x)+V for x(0) in D assuming bound B"""

        n=D.size()
        (K,L,LN)=self.compute_norms(f,B)
        KV=mag(norm(V))
        if LN>0.0: eLN=(exp(LN*h)-1)/(LN*h)
        else: eLN = Float(1.0)
        e = pow(h,2)*(2*KV*L*eLN)
        print "e:",e,"h:",h,"e/h^2:",e/pow(h,2)

        print "  K:",K,"KV:",KV,"L:",L,"LN:",LN,"eLN:",eLN,"e:",e

        E=e*IntervalVector(n,Interval(-1,+1))

        # Set up estimate for differential equation
        # Use affine wi=ai0+(t-h/2)ai1/h
        #   with |ai0|<=Vi, |ai1|<=3Vi, |wi(t)|<=5Vi/2, and |w'(t)|<=3Vi/h.
        # Alternatively use step inputs wi=ai,0 for t<h/2 and wi=ai,1 for t>h/2
        #   with |ai0|,|ai1|<= 2Vi , and |wi(t)|<=2Vi.

        # The flow is a function of n state variables, 2n parameter variables and 1 time variable
        # Assume for now noise in all variables
        # We have dot[xi](t)=fi(x(t))+wi(t)
        swp=self.sweeper
        DVH=join_list([D,V,E,Interval(-h,+h)])
        #print "DVH["+str(DVH.size())+"]:",DVH
        z=ScalarTaylorFunction(DVH,swp)
        xat=VectorTaylorFunction.identity(DVH,swp)
        w=VectorTaylorFunction(DVH,VectorFunction(n,3*n+1),swp)
        for i in range(n): w[i]=xat[n+i]
        #print "w:",w

        x0=VectorTaylorFunction(n,DVH,swp)
        for i in range(n): x0[i]=xat[i]
        #print "x0:",x0

        phi=VectorTaylorFunction(n,DVH,swp)
        for i in range(n): phi[i]=ScalarTaylorFunction(DVH,swp)+B[i]
        #print "phi0:",phi

        for i in range(6):
            phi=antiderivative(compose(f,phi)+w,3*n)+x0
            #print i,phi.errors()
        print (derivative(phi,3*n)-(compose(f,phi)+w)).range()

        for i in range(n):
            phi[i]=phi[i]+xat[2*n+i]

        return phi



def damped_harmonic(integrator,evolution_time,damping,noise,delta):
    T=evolution_time
    d=damping
    v=noise
    e=delta
    
    I=Interval(-1,+1)
    x=VectorFunction.identity(2)
    f=VectorFunction([-x[0]*d-x[1],x[0]-x[1]*d])
    V=IntervalVector([I*v[0],I*v[1]])
    X0=Box([[1.0-e,1.0+e],[-e,+e]])
    
    start_time=time.clock()
    (evolve_sets,reach_sets) = integrator.solve(f,V,X0,T)
    end_time=time.clock() - start_time
    print 'end_time=', end_time

    #print reach_sets

    #print [set.bounding_box() for set in reach_sets]
    #print [2*(box.bounding_box()).radius() for box in reach_sets]
   
    fig=Figure()
    fig.set_bounding_box(Box([[-2,1.5],[-1.3,2]]))
    fig.set_line_colour(0.0,0.0,0.0)
    fig.set_fill_colour(1.0,0.0,0.8)
    for set in reach_sets: fig.draw(set)
    fig.set_fill_colour(0.75,0.0,0.6)
    for set in evolve_sets: fig.draw(set)
    fig.write('damped_harmonic-test')



if __name__=='__main__':
    #damped_harmonic( evolution_time=2*pi, damping=1.0/4, noise=(0.0,0.1), step_size=8.0/32 )
    #damped_harmonic( evolution_time=2*pi, damping=0.0, noise=(0.0,0.1), delta=0.01, step_size=2*pi/50 )
    integrator = DifferentialInclusionIntegrator3rdOrder(ThresholdSweeper(1e-8),2*pi/32, 64, 32)

    damped_harmonic( integrator, evolution_time=2*pi, damping=0.0, noise=[0.0,0.1], delta=0.01),
#damped_harmonic( evolution_time=6.5, damping=1.0/16, noise=1.0/8, step_size=1.0/4 )
