/***************************************************************************
 *            macros/precondition.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_MACROS_PRECONDITION_H
#define ARIADNE_MACROS_PRECONDITION_H

#include <iostream>
#include <exception>

#include "macros/throw.h"

/*! \brief Evaluates \a expression in a boolean context and checks if the result is \a true. 
 *  If the result is false, throws an exception. */
#define ARIADNE_PRECONDITION(expression) \
{ \
  bool result = (expression); \
  if(!result) {                                                    \
    std::cerr << __FILE__ << ":" << __LINE__ << ": " << __FUNCTION__ << ": Precondition `" << #expression << "' failed.\n" << std::endl; \
    ARIADNE_THROW(std::runtime_error,__PRETTY_FUNCTION__,": Precondition `" << #expression << "' failed.") \
  } \
} \

#endif // ARIADNE_MACROS_PRECONDITION_H
