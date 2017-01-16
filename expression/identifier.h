/***************************************************************************
 *            identifier.h
 *
 *  Copyright 2008-16  Pieter Collins
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


/*! \file identifier.h
 *  \brief Strings used as names for constants and variables.
 */

#ifndef ARIADNE_IDENTIFIER_H
#define ARIADNE_IDENTIFIER_H

#include "utility/string.h"

namespace Ariadne {

//! \ingroup ExpressionModule
//! \brief A class representing the name of a variable.
//! \details A proxy for a standard string; used to distinguish a string used as a variable name from a value.
//! \sa Variable
class Identifier : public String
{
  public:
    Identifier() : String() { }
    Identifier(const char* cstr) : String(cstr) { }
    //! \brief Construct an identifier from a standard string.
    Identifier(const std::string& str) : String(str) { }
};

} // namespace Ariadne

#endif /* ARIADNE_IDENTIFIER_H */
