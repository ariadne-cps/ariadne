/***************************************************************************
 *            logging.h
 *
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
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

#ifndef _ARIADNE_LOGGING_H
#define _ARIADNE_LOGGING_H

#include <iostream>

namespace Ariadne {

#ifdef DEBUG
  static const int default_verbosity=1; 
#else
  static const int default_verbosity=0; 
#endif

#define ARIADNE_SET_VERBOSITY(where,level) \
  { where::verbosity=level; }

#define ARIADNE_LOG(level,msg) \
  if(verbosity > 
namespace Numeric { static int verbosity=default_verbosity; }
namespace LinearAlgebra { static int verbosity=default_verbosity; }
namespace Combinatoric { static int verbosity=default_verbosity; }
namespace Geometry { static int verbosity=default_verbosity; }
namespace System { static int verbosity=default_verbosity; }
namespace Evaluation { static int verbosity=default_verbosity; }
  
static int maximum_verbosity() {
  return std::max(
          std::max(Numeric::verbosity,LinearAlgebra::verbosity),
          std::max(
            std::max(Combinatoric::verbosity,Geometry::verbosity),
            std::max(System::verbosity,Evaluation::verbosity)
          )
        );
  

}




} // namespace Ariadne

#endif // _ARIADNE_LOGGING_H
