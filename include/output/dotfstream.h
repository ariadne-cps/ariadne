/****************************************************************************
 *            dotfstream.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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

#ifndef ARIADNE_DOTFSTREAM_H
#define ARIADNE_DOTFSTREAM_H

/*! \file dotfstream.h
 *  \brief Dot graph output
 */
 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>

#include "../system/hybrid_automaton.h"

namespace Ariadne {
  namespace Output {
 
    
    class dotfstream
      : private std::ofstream 
    {
     public:
      dotfstream() : std::ofstream() { }
      ~dotfstream() { this->close(); }
      
      void open(const char* fn) { this->std::ofstream::open(fn); }
      void close() { this->std::ofstream::close(); }
    };

    
    template<class R>
    dotfstream& 
    operator<<(dotfstream& dots, const System::HybridAutomaton<R>& ha) 
    {
      typedef typename System::HybridAutomaton<R>::discrete_transition_iterator transition_iterator;
      
      std::string f_name=ha.name();
      std::ostream& os=dots;

      size_t arc_number=0;
      
      os << "digraph \""<< ha.name()<<"\" {" << std::endl
         << " rankdir=LR; "<< std::endl
         << " node [shape = circle]; "<< std::endl;
      
      for (transition_iterator iter=ha.transitions().begin(); iter!=ha.transitions().end(); ++iter) {
        System::DiscreteTransition<R>& dt=*iter;

        os << "\"" <<  dt.source().id() << "\" -> \"" 
           << dt.destination().id()
           << "\" [ label=\"a_" << dt.id()++ << "\" ]; " << std::endl;
      }      
      
      os << "}" << std::endl;

      return dots;
    }

  }
}


#endif /* ARIADNE_DOTFSTREAM_H */
