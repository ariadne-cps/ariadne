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
from ariadne.evaluation import *
from ariadne.geometry import *
from ariadne.linear_algebra import *
from ariadne.system import *
import sys

def try_to_reset_from(H, d_node, c_set, time, time_step, reach_set):
  for trans in H.get_discrete_transitions_leaving(d_node):
    act=trans.activation()
    if (interiors_intersect(c_set, act)):
      reset=trans.reset()
      dest_d_node=trans.destination()
      dest_c_set=reset(c_set)
      bounded_time_reachability_with_reach_set(H, dest_d_node, dest_c_set, time, time_step, reach_set)
	

def bounded_time_reachability_with_reach_set(H, d_node, c_set, time, time_step, reach_set):
  if not (d_node in H.discrete_nodes()):
    raise 'bounded_time_reachability(...): initial discrete node does not belong to the automaton'
  last_reached=c_set
  inv=d_node.invariant()
  while ((time>=0) and (interiors_intersect(last_reached,inv))):
    reached=integrate(d_node.dynamic(),last_reached,time_step,time_step)
    reach_set.add_a_hybrid_region(d_node, reached)
    try_to_reset_from(H, d_node, reached, time, time_step, reach_set)
    last_reached=reached
    time=time-time_step

def bounded_time_reachability(H, d_node, c_set, time, time_step):
  reach_set=HybridReachSet()
  bounded_time_reachability_with_reach_set(H, d_node, c_set, time, time_step, reach_set)
  return reach_set.reached_region()

