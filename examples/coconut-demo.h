/***************************************************************************
 *            coconut-demo.h
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

#ifndef COCONUT_DEMO_H
#define COCONUT_DEMO_H

#include "ariadne.h"

using namespace Ariadne;

// Variables and objects that will be used to define the system 
extern HybridAutomaton automaton;
extern Box initial_box;
extern DiscreteState initial_state;
extern DiscreteState safe;
extern DiscreteState unsafe;

// Evolution parameters
extern Grid grid;
extern float MAX_ENCL_RADIUS; /// Maximum enclosure radius
extern float MAX_STEP_SIZE; /// Maximum step size
extern float LOCK_TOGRID_TIME;
extern int LOCK_TOGRID_STEPS;
extern int MAX_GRID_DEPTH;

// Graphic output parameters
extern int nvar;
extern int xvar;
extern int yvar;
extern Box graphic_box;

// Initialization function
void build_automaton();

#endif