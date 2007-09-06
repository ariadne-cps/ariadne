/***************************************************************************
 *            pointer.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

/*! \file pointer.h
 *  \brief Smart pointers.
 */

#ifndef ARIADNE_POINTER_H
#define ARIADNE_POINTER_H

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

namespace Ariadne {
  namespace Base {
    using boost::shared_ptr;
    using boost::scoped_ptr;
  }
}

#ifdef DOXYGEN
namespace Ariadne {
  /*! \brief Reference-counted shared pointer. Initialise with a dynamically-allocated object; destruction is performed automatically. */
  template<class T> class shared_ptr { };
  /*! \brief Scope smart pointer. Initialise with a dynamically-allocated object; destruction is performed automatically when object goes out of scope. Cannot be used in containers. */
  template<class T> class scoped_ptr { };
}
#endif

#endif /* ARIADNE_POINTER_H */
