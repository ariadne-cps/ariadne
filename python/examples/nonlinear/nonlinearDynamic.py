
from ariadne import *
#from ariadne.numeric import *

import sys

# Construct Python-based nonlinear function
print "Constructing and evaluating nonlinear integration system defined in Python"
def fun(x):
     
     yF   = 0.5
     VMAX = 0.5
     K    = 10
     A    = 0.03

     xdot = x[3]*cos(x[2]);
     ydot = x[3]*sin(x[2]);
     thetadot = -((x[1]-yF)*sin(x[2])*x[3]/x[2])-K*x[3]*x[2];	
 
     if x[3]<VMAX:
 	      vdot = A
     else:
        vdot = 0

     return [ xdot, ydot, thetadot, vdot ]

fun.result_size=4
fun.argument_size=4
nl_fun=AriadneFunction(fun)
print "nl_function =",nl_fun

# Compare the behaviour of fun e nl_fun
for i in range(11):
  x = Float(i*0.1)
  print "fun(",x,") = ",fun([x, x, x, x])
  #print "nl_fun(",x,") = ",nl_fun(Vector([x, x, x, x]))

#sys.exit(0)

# Define the corresponding VectorField
vector_field=VectorField(nl_fun)

# Construct initial set and bounding set
initial_set=RectangularSet([[0.0,0.00001],[0.0,0.90001],[0.10001,0.57001],[0.0001,0.45000]])
bounding_set=RectangularSet([[0,20],[-1,2],[-1.6,1.6],[0,0.6]])

# Construct integrator (use Kuhn integrator for nonlinear systems)
integrator=KuhnIntegrator(3,4)

# Construct evolver
parameters=EvolutionParameters()
parameters.set_grid_length(0.1)
parameters.set_bounding_domain_size(2.0)
parameters.set_maximum_step_size(0.125)

evolver=VectorFieldEvolver(parameters,integrator)

# Compute chain-reachable set
print "Computing chain-reachable set..."
#chain_reach_set = evolver.chain_reach(vector_field,initial_set, bounding_set.bounding_box())
chain_reach_set = evolver.upper_reach(vector_field,initial_set,Rational(0.5))
print "Found", chain_reach_set.size(), "cells in grid with", chain_reach_set.capacity(), "cells."

# Reduce chain-reachable set to partition tree set
#print "Computing tree representation of chain-reachable set..."
#chain_reach_tree_set = PartitionTreeSet(chain_reach_set)
#print "Reduced to", chain_reach_tree_set.size(),"cells " \
#    "in partition tree with", chain_reach_tree_set.capacity(),"cells."

# Export to postscript output
print "Exporting to postscript output...",
epsbb=Box([[-4.1,4.1],[-4.1,4.1],[-4.1,4.1],[-4.1,4.1]]) # eps bounding box
eps=EpsPlot()
eps.open("nonlinear_system-mask_set.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("green")
eps.write(chain_reach_set)
eps.set_fill_colour("blue")
eps.write(initial_set)
eps.close()

# eps.open("nonlinear_system-tree_set.eps",epsbb)
# eps.set_pen_colour("black")
# eps.set_fill_colour("green")
# eps.set_line_style(0)
# eps.write(chain_reach_tree_set)
# eps.set_line_style(1)
# eps.set_fill_style(0)
# eps.write(chain_reach_tree_set.partition_tree())
# eps.set_fill_style(1)
# eps.set_fill_colour("blue")
# eps.write(initial_set)
# eps.close()
print " done."
