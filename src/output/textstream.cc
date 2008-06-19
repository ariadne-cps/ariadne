/****************************************************************************
 *            textstream.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins, Davide Bresolin
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl, bresolin@sci.univr.it
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

#include "base/stlio.h"
#include "output/textstream.h"

namespace Ariadne { 


textstream::~textstream()
{
}

textstream::textstream()
  : _os_ptr(&std::cout)
{
}

textstream::textstream(std::ostream& os)
  : _os_ptr(&os)
{
}

void
textstream::redirect(std::ostream& os)
{
  this->_os_ptr=&os;
}

void 
textstream::writenl() 
{
  (*this) << std::endl;
}


textfstream::textfstream()
  : textstream(), _ofs_ptr(new std::ofstream())
{
  this->textstream::redirect(*_ofs_ptr);
}

textfstream::~textfstream()
{
  this->close(); 
  delete this->_ofs_ptr;
}

void
textfstream::open(const char* fn)
{
  this->_ofs_ptr->open(fn);
}


void 
textfstream::close() 
{
  this->_ofs_ptr->close();
}




}
