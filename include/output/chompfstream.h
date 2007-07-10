/****************************************************************************
 *            chompfstream.h
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_CHOMPFSTREAM_H
#define ARIADNE_CHOMPFSTREAM_H

/*! \file chompfstream.h
 *  \brief Output to CHomP computational homology programs
 */
 
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

#include "../combinatoric/lattice_set.h"
#include "../combinatoric/lattice_map.h"

namespace Ariadne {
  namespace Output {
    
    class chompfstream {
     public:
      chompfstream();
      chompfstream(const char* fn);
      ~chompfstream();

      void open(const char* fn);
      void close();
     private:
      friend chompfstream& operator<<(chompfstream&, const char*);
      friend chompfstream& operator<<(chompfstream&, const Combinatoric::LatticeCell&);
      friend chompfstream& operator<<(chompfstream&, const Combinatoric::LatticeMaskSet&);
      friend chompfstream& operator<<(chompfstream&, const Combinatoric::LatticeMultiMap&);
     private:
      std::ofstream _ofs;
    };

  }
}

#endif /* ARIADNE_CHOMPFSTREAM_H */


