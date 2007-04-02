/***************************************************************************
 *            array.h
 *
 *  4 October 2004
 *  Copyright  2004  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or5
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
 
/*! \file array.decl.h
 *  \brief Forward declarations for STL style arrays.
 */

#ifndef ARIADNE_ARRAY_DECL_H
#define ARIADNE_ARRAY_DECL_H

namespace Ariadne {
  namespace Base {
    
    // Need to give default size in first declaration.
    template<class T, unsigned short int N=0> class array;
    
  }
}

  
#endif /* ARIADNE_ARRAY_DECL_H */
