/***************************************************************************
 *            classes.h
 *
 *  Wed Sep  1 16:03:39 2004
 *  Copyright  2004  Alberto Casagrande
 *  casagrande@dimi.uniud.it
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#ifndef _CLASSES_H
#define _CLASSES_H

namespace Ariadne {
/* for the class documentation look into the other dot-h-files */

/* from include/arc.h */
class Arc;
class LeavingArc;

/* from include/automaton.h */
class Automaton;

/* from include/basic_maintainer.h */
class BasicMaintainer;

/* from include/basic_set.h */
class BasicSet;

/* from include/basic_set_list.h */
class BasicSetList;
class Ex_BSL;
class In_BSL;

/* from include/buffer.h */
class Buffer;

/* from include/cluster_list.h */
class Interval;

/* from include/location.h */
class Location;

/* from include/maintain.h */
class Maintainer;
class HMaintainer;
class MaintainSystem;

/* from include/map.h */
class Map;
class LinearMap;

/* from include/reset.h */
class Reset_System;

/* from include/set.h */
class ASet;
class OverSet;
class UnderSet;
class Set;
class HSet;
class RSet;

/* from include/solver.h */
class Solver;

/* from include/trace.h */
class LocationTrace;

/* from include/variable.h */
class Variable;

/* from include/vectorfield.h */
class VectorField;
class LinearVectorField;


/* TODO: implement the folowing classes */
class _Maintain_Changes;
}

#endif /* _CLASSES_H */

