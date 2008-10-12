/***************************************************************************
 *            tribool.h
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
 
/*! \file tribool.h
 *  \brief Three-valued logic variable.
 */

#ifndef ARIADNE_TRIBOOL_H
#define ARIADNE_TRIBOOL_H

#include <boost/logic/tribool.hpp>
#include <boost/logic/tribool_io.hpp>

using boost::logic::tribool;
using boost::logic::indeterminate;

inline bool definitely(tribool tb) { return tb; }
inline bool possibly(tribool tb) { return tb || indeterminate(tb); }

inline tribool operator^(tribool tb1, tribool tb2) { return (tb1&&!tb2)||(!tb1&&tb2); }

#endif // ARIADNE_TRIBOOL_H
