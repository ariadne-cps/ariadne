/***************************************************************************
 *            variable.cpp
 *
 *  Fri Aug  6 12:08:31 2004
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

#include "variable.h"

using namespace Ariadne;
	
Variable::Variable(const std::string &name){
	this->name=name;
	this->id=this->smallest_unused_id;
	this->smallest_unused_id++;

	this->type=REAL_VAR;
}

Variable::Variable(const std::string &name, VariableID id){
	this->name=name;
	this->id=id;
			
	this->type=REAL_VAR;
}
		
Variable::Variable(const std::string &name, VariableType type){
	this->name=name;
	this->type=type;
}
		

Variable::~Variable(){}

bool Variable::has_id(const VariableID id) const{
	return (this->id==id);
}

	
bool Variable::has_name(const std::string &name) const{
	return (this->name==name);
}

void Variable::set(const std::string &name, VariableID id){
	this->name=name;
	this->id=id;
}

