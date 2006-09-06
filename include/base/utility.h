/***************************************************************************
 *            utilty.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file utility.h
 *  \brief Miscellaneous, general-purpose functions.
 */

#ifndef _ARIADNE_UTILITY_H
#define _ARIADNE_UTILITY_H

#include <iosfwd>
#include "../declarations.h"

namespace Ariadne {
  namespace Base {
    
    /*! \brief Convert an element \a x of type \a Arg to type \a Res. */
    //template<typename Res, typename Arg> inline Res convert_to(const Arg& x) { return Res(x); }

    /*! \brief The name of class T. */
    //template<typename T> inline std::string name();
    
  }    
}

#endif /* _ARIADNE_UTILITY_H */
