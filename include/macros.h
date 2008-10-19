/***************************************************************************
 *            macros.h
 *
 *  Copyright 2008  Pieter Collins
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
 
/*! \file macros.h
 *  \brief Commonly used macros.
 */

#ifndef ARIADNE_MACROS_H
#define ARIADNE_MACROS_H

#include <sstream>
#include <stdexcept>

#include "exceptions.h"

#define ARIADNE_THROW(except,func,msg)          \
{ \
  std::stringstream ss; \
  ss << #except " in " << func << " " << msg;    \
  throw except(ss.str()); \
} \

#define ARIADNE_ASSERT(expression) \
{ \
  bool result = (expression); \
  if(!result) { \
    ARIADNE_THROW(std::runtime_error,__FILE__<<":"<<__LINE__<<": "<<__PRETTY_FUNCTION__,"Assertion `" << #expression << "' failed.\n"); \
  } \
} \

#define ARIADNE_LOG(level,msg)                  \
  if(verbosity >= level) { std::clog << msg << std::flush; }

#define ARIADNE_NOT_IMPLEMENTED                        \
  throw NotImplemented(__PRETTY_FUNCTION__); 


#endif // ARIADNE_MACROS_H
