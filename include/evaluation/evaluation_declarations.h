/***************************************************************************
 *            evaluation_declarations.h
 *
 *  Copyright 6 February 2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file evaluation_declarations.h
 *  \brief Forward declarations for the Evaluation module.
 */

#ifndef _ARIADNE_EVALUATION_DECLARATIONS_H
#define _ARIADNE_EVALUATION_DECLARATIONS_H

namespace Ariadne {
  namespace Evaluation {

    template <typename R> class Map;
    template <typename R> class AffineMap;

    template <typename R> class VectorField;
    template <typename R> class AffineVectorField;

    template <typename R> class HenonMap;
    template <typename R> class LorenzSystem;
      
  }
}

#endif /* _ARIADNE_EVALUATION_DECLARATIONS_H */
