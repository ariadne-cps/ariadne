/****************************************************************************
 *            arc.cpp
 *
 *  Wen Aug  4 17:44:20 2004
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
 *  Foundation, Inc., 59 Temple Place - Suite 330,
 */

#include "arc.h"
#include "location.h"
#include "set.h"
#include "map.h"

using namespace Ariadne;
	
LeavingArc::LeavingArc(Location *dest, ASet *act, Map *reset){

	this->destination=dest;
	this->activation=act;
	this->reset=reset;

}
		
LeavingArc::~LeavingArc() {
			
	delete this->activation;
	delete this->reset;

	// the destination location shouldn't be delete
	this->destination=NULL;

	// safer, but (probabily) useless stuff
	this->reset=NULL;
	this->activation=NULL;
}
	

Location* LeavingArc::get_destination() const {
	return this->destination;
}
	
const ASet* LeavingArc::get_activation() const {
	return this->activation;
}
 

Arc::Arc(Location *source, Location *dest, 
		ASet *act, Map *reset) : LeavingArc(dest,act,reset) {

	this->source=source;
}

Arc::~Arc() {
	delete this->activation;
	delete this->reset;		

	// useless, but safer code
	this->activation=NULL;
	this->reset=NULL;		
	this->source=NULL;
	this->destination=NULL;
}
	
Location* Arc::get_source() const {
	return this->source;
}

