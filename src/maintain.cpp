/***************************************************************************
 *            maintain.cpp
 *
 *  Fri Aug 06 10:45:12 2004
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
 
#include "maintain.h"
#include "basic_maintainer.h"

using namespace Ariadne;

Maintainer::Maintainer(BasicMaintainer* over, BasicMaintainer *under){
	this->over= over;
	this->under= under;
}

Maintainer::Maintainer(const Maintainer &templ) {
	this->over= (templ.over)->copy();
	this->under=(templ.under)->copy();
}
		
Maintainer::~Maintainer() {
	delete (this->over);	
	delete (this->under);

	// useless, but safer
	this->over=NULL;
	this->under=NULL;
}

