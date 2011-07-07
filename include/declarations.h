/***************************************************************************
 *            declarations.h
 *
 *  Copyright 2011  Pieter Collins
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

/*! \file declarations.h
 *  \brief Forward declarations of types and classes.
 */
#ifndef ARIADNE_DECLARATIONS_H
#define ARIADNE_DECLARATIONS_H

#include <iosfwd>
#include <boost/logic/tribool.hpp>

namespace Ariadne {

//! Internal name for output stream.
typedef std::ostream OutputStream;

//! Internal name for builtin boolean type.
typedef bool Bool;
//! Internal name for three-valued logical type.
typedef boost::tribool Tribool;
//! Internal name for string objects.
typedef std::string String;
//! Internal name for builtin unsigned integers.
typedef unsigned int Nat;
//! Internal name for builtin integers.
typedef int Int;

} // namespace Ariadne

#endif
