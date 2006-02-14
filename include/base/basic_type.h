/***************************************************************************
 *            basic_type.h
 *
 *  31 January 2006
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

/*! \file basic_type.h
 *  \brief Basic type definitions.
 */

#ifndef _ARIADNE_BASIC_TYPE_H
#define _ARIADNE_BASIC_TYPE_H

#include <base/array.h>
#include <base/sequence.h>

namespace Ariadne {
    
    typedef unsigned short dimension_type;
    typedef size_t size_type;
    typedef int index_type;

    //typedef std::vector<bool> BooleanArray;
    typedef array<bool> BooleanArray;
    typedef sequence<dimension_type> SubdivisionSequence;

}

#endif /* _ARIADNE_BASIC_TYPE_H */
