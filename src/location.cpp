/***************************************************************************
 *            location.cpp
 *
 *  Thu Aug  5 16:30:32 2004
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

#include <string>	
#include <list>

#include "location.h"
#include "arc.h"
#include "set.h"
#include "vectorfield.h"

using namespace Ariadne;
	
Location::Location(const std::string &name, VectorField *field, ASet *inv){
			
	this->name=name;
	this->field=field;
	this->invariant=inv;

	this->id=this->smallest_unused_id;

	this->smallest_unused_id++;		
};
		
Location::~Location(){
	delete this->field;
	delete this->invariant;	

	this->field=NULL;
	this->invariant=NULL;

	//! delete the leaving arc list
	std::list<LeavingArc *>::iterator i;

	for ( i=(this->leaving_arcs).begin(); 
			i!=(this->leaving_arcs).end(); 
			i++ ) {
		
		delete (*i);
		(this->leaving_arcs).pop_front();
	}

	

}
	
const VectorField* Location::get_vectorfield() const{
	return this->field;
}
		
const std::string& Location::get_name() const{
	return this->name;
}
	
const LeavingArc* Location::add_a_leaving_arc(Location *dest, 
		ASet *act, Map *reset){
		
	//! create a new leaving arc 
	LeavingArc *l_arc= new LeavingArc(dest,act,reset);

	/*! push the new arc in front of the leaving
	 * arc list */
	(this->leaving_arcs).push_front(l_arc);

	return l_arc;
}
		
void Location::Init() {
	this->smallest_unused_id=0;
}
		
const ASet* Location::get_invariant() const {
	return this->invariant;
}

const std::list<LeavingArc *>& Location::get_leaving_arcs() const {
	return this->leaving_arcs;
}

bool Location::operator==(const Location &a) const{
	return (this->id==a.id);	
} 
		
bool Location::operator!=(const Location &a) const{
	return (this->id!=a.id);	
} 

