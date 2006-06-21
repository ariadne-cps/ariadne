##############################################################################
#            system_evolution.py
#
#  Copyright 2006  Pieter Collins, Alberto Casagrande
#  Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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

from ariadne.base import *
from ariadne.numeric import *
from ariadne.evaluation import *
from ariadne.geometry import *
from ariadne.linear_algebra import *
from ariadne.system import *
import sys

def write_reached(eps,reached,color):
  eps.set_fill_colour(color)
  for bs in reached.reached_regions():
    eps.write(bs.continuous_set())

def try_to_grid_reset_from(H, d_node, c_set, time, time_step, reach_set, flowed_set, max_jump, verbatim):
  eps.write(bs.continuous_set())
  if (max_jump==0):
    return
  for trans in H.get_discrete_transitions_leaving(d_node):
    act=trans.activation()
    if (interiors_intersect(c_set, act)):
      if (max_jump>0):
        max_jump=max_jump-1
      if (verbatim == 'yes'):
        print 'Jumping over',trans.name()
      reset=trans.reset()
      dest_d_node=trans.destination()
      dest_c_set=reset(c_set)
      bounded_time_grid_reachability_with_reach_set(H, dest_d_node, dest_c_set, time, time_step, reach_set, flowed_set, max_jump, verbatim)
	

def bounded_time_grid_reachability_with_reach_set(H, d_node, c_set, time, time_step, reach_set, flowed_set, max_jump,verbatim):
  if not (d_node in H.discrete_nodes()):
    raise 'bounded_time_reachability(...): initial discrete node does not belong to the automaton'
  if (verbatim == 'yes'):
    print 'Flowing on',d_node.name()
  grid_c_set=reach_set.grid_over_approximation(d_node, c_set)
  inv=(reach_set.reached_region_on(d_node)).bounds_set()
  already_flowed_set=(flowed_set.reached_region_on(d_node)).continuous_set()
  new_set=difference(grid_c_set,already_flowed_set)
  last_reached=regular_intersection(new_set,inv)
  #lohner=C1LohnerIntegrator(0.125,0.5,0.0625)
  lohner=C1LohnerIntegrator(0.5,1.0,0.1)
  vf=d_node.dynamic()
  while ((time>=0) and (not last_reached.empty())):
    if (verbatim == 'yes'):
      print 'Time',time
    already_flowed_set=(flowed_set.reached_region_on(d_node)).continuous_set()
    new_set=difference(last_reached,already_flowed_set)
    last_reached=regular_intersection(new_set,inv)
    reached=lohner.reach(vf, last_reached, inv, Real(time_step))
    reach_set.join_over_approximation(d_node, reached)
    try_to_grid_reset_from(H, d_node, reached, time, time_step, reach_set, flowed_set, max_jump, verbatim)
    last_reached=lohner.integrate(vf, last_reached, inv, Real(time_step))
    flowed_set.join_under_approximation(d_node, last_reached)
    print (flowed_set.reached_region_on(d_node)).continuous_set()
    time=time-time_step

def bounded_time_grid_reachability(H, d_node, c_set, time, time_step, approx, max_jump, verbatim):
  reach_set=HybridGridReachSet(H, approx)
  flowed_set=HybridGridReachSet(H, approx)
  bounded_time_grid_reachability_with_reach_set(H, d_node, c_set, time, time_step, reach_set, flowed_set, max_jump, verbatim)
  return [reach_set, flowed_set]

def try_to_reset_from(H, d_node, c_set, time, time_step, reach_set, flowed_set, max_jump, eps_img, verbatim):
  if (max_jump<=0):
    return
  for trans in H.get_discrete_transitions_leaving(d_node):
    act=trans.activation()
    if (interiors_intersect(c_set, act)):
      if (max_jump>0):
        max_jump=max_jump-1
      if (verbatim == 'yes'):
        print 'Jumping over',trans.name()
	print 'Reseting...',
      reset=trans.reset()
      dest_d_node=trans.destination()
      dest_c_set=reset(c_set)
      grid_c_set=reach_set.grid_over_approximation(dest_d_node, dest_c_set)
      try:
       new_c_set=ZonotopeListSet(RectangleListSet(PartitionTreeSet(grid_c_set)))
      except:
       new_c_set=ZonotopeListSet(RectangleListSet(grid_c_set))
      #new_c_set=dest_c_set
      if (verbatim == 'yes'):
        print 'done'
      if not (len(eps_img)==0):
        for eps in eps_img:
          write_reached(eps,reach_set,"blue")
          write_reached(eps,flowed_set,"red")
      bounded_time_reachability_with_reach_set(H, dest_d_node, new_c_set, time, time_step, reach_set, flowed_set, max_jump, eps_img, verbatim)
      print 'Backtracing...'
	

def bounded_time_reachability_with_reach_set(H, d_node, c_set, time, time_step, reach_set, flowed_set, max_jump, eps_img, verbatim):
  if not (d_node in H.discrete_nodes()):
    raise 'bounded_time_reachability(...): initial discrete node does not belong to the automaton'
  if (verbatim == 'yes'):
    print 'Flowing on',d_node.name()
  last_reached=c_set
  inv=d_node.invariant()
  last_reached=touching_intersection(last_reached,inv)
  vf=d_node.dynamic()
  #lohner=C1LohnerIntegrator(0.125,0.5,0.0625)
  lohner=C1LohnerIntegrator(0.125,0.005,5)
  grid_c_set=flowed_set.grid_over_approximation(d_node, c_set)
  already_flowed_set=(flowed_set.reached_region_on(d_node)).continuous_set()
  new_set=difference(grid_c_set,already_flowed_set)
  new_flow=not (new_set.empty())
  flowed_set.join_under_approximation(d_node, last_reached)
  while ((time>=0) and (new_flow)):
    if (verbatim == 'yes'):
      print 'Time %.3f' % time
      print 'Computing reached region...'
    reached=lohner.reach(vf, last_reached, Real(time_step))
    if (verbatim == 'yes'):
      print 'done'
      print 'Memoizing reached region...'
    reached=touching_intersection(reached,inv)
    #reach_set.add_a_hybrid_region(d_node, reached)
    reach_set.join_over_approximation(d_node,reached)
    if (verbatim == 'yes'):
      print 'done'
    try_to_reset_from(H, d_node, reached, time, time_step, reach_set, flowed_set, max_jump, eps_img, verbatim)
    if (verbatim == 'yes'):
      print 'Integrating region...'
    last_reached=lohner.integrate(vf, last_reached, Real(time_step))
    last_reached=touching_intersection(last_reached,inv)
    if (verbatim == 'yes'):
      print 'done'
      print 'Memoizing flowed region...'
    grid_c_set=flowed_set.grid_over_approximation(d_node, last_reached)
    already_flowed_set=(flowed_set.reached_region_on(d_node)).continuous_set()
    new_set=difference(grid_c_set,already_flowed_set)
    new_flow=not (new_set.empty())
    flowed_set.join_under_approximation(d_node, last_reached)
    new_flow=interiors_intersect(last_reached,inv)
    if (verbatim == 'yes'):
      print 'done'
    time=time-time_step

def bounded_time_reachability(H, d_node, c_set, time, time_step, approx, max_jump, eps_img, verbatim):
  #reach_set=HybridReachSet()
  #approx=1.0/(2**approx)
  reach_set=HybridGridReachSet(H, approx)
  flowed_set=HybridGridReachSet(H, approx)
  bounded_time_reachability_with_reach_set(H, d_node, c_set, time, time_step, reach_set, flowed_set, max_jump, eps_img, verbatim)
  return [reach_set, flowed_set]
 
