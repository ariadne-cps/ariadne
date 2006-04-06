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
import sys

def bounded_time_reachability(H, d_node, c_set, time, time_step):
	if not (d_node in H.discrete_nodes()):
		raise 'bounded_time_reachability(...): initial discrete node does not belong to the automaton'
	last_reached=c_set
	reach_set=[last_reached]
	inv=d_node.invariant()
	while ((time>=0) and (interiors_intersect(last_reached,inv))):
		reached=integrate(d_node.dynamic(),last_reached,time_step,time_step)
		for trans in H.get_discrete_transitions_leaving(d_node):
			act=trans.activation()
			if (interiors_intersect(reached,act)):
				reset=trans.reset()
				dest_d_node=trans.destination()
				dest_c_set=reset(reached)
				new_reach=bounded_time_reachability(H, dest_d_node, dest_c_set, time, time_step)
				reach_set.extend(new_reach)
		
		reach_set.append(reached)
		last_reached=reached
		time=time-time_step
	return reach_set

