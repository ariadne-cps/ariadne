/***************************************************************************
 *            system.h
 *
 *  Copyright 2009  Davide Bresolin
 *
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

#ifndef COCONUT_SYSTEM_H
#define COCONUT_SYSTEM_H

#include "ariadne.h"

using namespace Ariadne;

// Variables and objects that will be used to define the system
MonolithicHybridAutomaton automaton;
Box initial_box;
DiscreteState initial_state;
DiscreteState safe(1001);
DiscreteState unsafe(999);

// Evolution parameters
Grid grid;
float MAX_ENCL_RADIUS = 1.0; /// Maximum enclosure radius
float MAX_STEP_SIZE = 0.25; /// Maximum step size
float LOCK_TOGRID_TIME = 10.0;
int LOCK_TOGRID_STEPS = 1;
int MAX_GRID_DEPTH = 16;

// Graphic output parameters
int nvar = 2;
int xvar = 0;
int yvar = 1;
Box graphic_box(2, -1.0, 1.0, -1.0,1.0);

// Initialization function
void build_automaton();

#endif