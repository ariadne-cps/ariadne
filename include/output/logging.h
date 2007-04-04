/***************************************************************************
 *            output/logging.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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
#include <fstream>

#define ARIADNE_LOG(level,msg) \
  if(verbosity > level) { std::clog << msg; }

namespace Ariadne {
  namespace Output {

    // Global log output file
    extern std::ofstream log_file_stream;
  
    /*! \brief Redirect logging output to file \a filename. */
    void redirect_log(const char* filename);

  }

  namespace Numeric { extern int verbosity; }
  namespace LinearAlgebra { extern int verbosity; }
  namespace Combinatoric { extern int verbosity; }
  namespace Geometry { extern int verbosity; }
  namespace System { extern int verbosity; }
  namespace Evaluation { extern int verbosity; }
}

#endif // _ARIADNE_LOGGING_H
