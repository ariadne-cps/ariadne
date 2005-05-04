/***************************************************************************
 *            exception.h
 *
 *  2 May 2005
 *  Copyright  2005  Pieter Collins, Alberto Casagrande
 *  Email  Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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
 
/*! \file exception.h
 *  \brief Exceptions, error handling and assertions.
 */

#ifndef _EXCEPTION_H
#define _EXCEPTION_H

#include <stdexcept>
#include <iosfwd>

namespace Ariadne {


    class invalid_input : public std::runtime_error {
      public:
	invalid_input(const std::string& str) : std::runtime_error(str) { }
    };
    
}

#endif /* _EXCEPTION_H */
