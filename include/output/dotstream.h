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

#ifndef ARIADNE_DOTSTREAM_H
#define ARIADNE_DOTSTREAM_H

/*! \file dotstream.h
 *  \brief Dot graph output
 */
 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstring>

#include "../system/set_based_hybrid_automaton.h"
#include "../system/constraint_based_hybrid_automaton.h"

namespace Ariadne {
  namespace Output {
 
    
    /*!\brief A stream for output to the graph drawing package Dot. */
    class dotstream
    {
     public:
      dotstream() : _os_ptr(&std::cout) { }
      dotstream(std::ostream& os) : _os_ptr(&os) { }
      void redirect(std::ostream& os) { this->_os_ptr=&os; }
      std::ostream& ostream() { return *this->_os_ptr; }
      ~dotstream() { }
     private:
      std::ostream* _os_ptr;
    };

    class dotfstream
      : private dotstream 
    {
     public:
      dotfstream() : dotstream(), _ofs_ptr(new std::ofstream()) { this->dotstream::redirect(*this->_ofs_ptr); }
      ~dotfstream() { this->close(); delete this->_ofs_ptr; }
      
      void open(const char* fn) { this->_ofs_ptr->open(fn); }
      void close() { this->_ofs_ptr->close(); }
     private:
      std::ofstream* _ofs_ptr;
    };

    
    template<class R>
    dotstream& 
    operator<<(dotstream& dots, const System::ConstraintBasedHybridAutomaton<R>& ha) 
    {
      typedef typename System::ConstraintBasedHybridAutomaton<R>::discrete_transition_const_iterator transition_iterator;
      
      std::string f_name=ha.name();
      std::ostream& os=dots;

      size_t arc_number=0;
      
      os << "digraph \""<< ha.name()<<"\" {" << std::endl
         << " rankdir=LR; "<< std::endl
         << " node [shape = circle]; "<< std::endl;
      
      for (transition_iterator iter=ha.transitions().begin(); iter!=ha.transitions().end(); ++iter) {
        System::ConstraintBasedDiscreteTransition<R>& dt=*iter;

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
