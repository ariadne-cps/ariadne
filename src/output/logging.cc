/***************************************************************************
 *            logging.cc
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

#include "output/logging.h"

namespace Ariadne {
  namespace Output {
    std::ofstream log_file_stream;

    void redirect_log(const char* filename) {
      if(log_file_stream.is_open()) {
        log_file_stream.close();
      }
      log_file_stream.open(filename);
      std::clog.rdbuf( log_file_stream.rdbuf() );
    }

  }

  namespace Numeric { int verbosity=0; } 
  namespace LinearAlgebra { int verbosity=0; } 
  namespace Combinatoric { int verbosity=0; } 
  namespace Geometry { int verbosity=0; } 
  namespace System { int verbosity=0; } 
  namespace Evaluation { int verbosity=0; } 
}
