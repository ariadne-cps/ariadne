/***************************************************************************
 *            automaton.cpp
 *
 *  Wen Aug  4 17:46:31 2004
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
 
#include <string>
#include <list>

#include "automaton.h"
#include "set.h"
#include "location.h"
#include "variable.h"

using namespace Ariadne;

void Automaton::Init() {
	this->smallest_unused_id=0;
}

		
Automaton::Automaton(const std::string &name, 
		std::list<Location *> *locations, 
		std::list<HSet> *init,
		VecVariable *v_variable) {	

	this->name=name;
	this->locations=locations;
	this->i_regions=init;

	this->id=this->smallest_unused_id;

	this->smallest_unused_id++;

	this->v_variable=v_variable;
}
		
Automaton::Automaton(const std::string &name, std::list<HSet> *init) {	
	
	this->name=name;
	this->locations=new std::list<Location *>();
	this->i_regions=init;
	this->v_variable=NULL;

	this->id=this->smallest_unused_id;

	this->smallest_unused_id++;
}
		
Automaton::Automaton(const std::string &name) {
	
	this->name=name;
	this->locations=new std::list<Location *>();
	this->i_regions=new std::list<HSet>();
	this->v_variable=NULL;
			
	this->id=this->smallest_unused_id;

	this->smallest_unused_id++;
}
	
Automaton::~Automaton(){
	std::list<Location *>::iterator i;

	/* delete all the location list's elements 
	 * and clear the list itself */
	for ( i=(this->locations)->begin(); 
			i!=(this->locations)->end(); 
			i++ ) {
		
		delete (*i);
		(this->locations)->pop_front();
	}	
	
	delete this->locations;
	
	/* clear the initial region list */
    	(this->i_regions)->clear();
			
	delete this->i_regions;

	/* TODO: check this code */
	delete this->v_variable;
			
	/* the following code is useless, but could be safer
	 * to include the following code */
	this->locations=NULL;
	this->i_regions=NULL;
	this->v_variable=NULL;
			
}	

const std::string &Automaton::get_name(){
	return this->name;
}

bool Automaton::operator==(const Automaton &a) {
	return (this->id==a.id);	
} 
		
bool Automaton::operator!=(const Automaton &a) {
	return (this->id!=a.id);	
}

