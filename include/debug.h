/***************************************************************************
 *            debug.h
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

#ifndef _ARIADNE_DEBUG_H
#define _ARIADNE_DEBUG_H

#include <iostream>

namespace Ariadne {

#ifdef DEBUG
  static const int default_debug_level=1; 
#else
  static const int default_debug_level=0; 
#endif

  class dbgstream;
  template<class T> dbgstream& operator<<(dbgstream& dbgs, const T& t);

  class dbgstream : public std::ostream
  {
   public:
    dbgstream(std::ostream& os, int debug_level=default_debug_level) : _stream(os), _debug_level(debug_level) { }
    dbgstream(int debug_level=default_debug_level) : _stream(std::cerr), _debug_level(debug_level) { }
   private:
    template<class T> friend dbgstream& operator<<(dbgstream& dbgs, const T& t);
   private:
    std::ostream& _stream;
    int _debug_level;
  };


  template<class T>
  inline
  dbgstream& 
  operator<<(dbgstream& dbgs, const T& t) 
  {
    if(dbgs._debug_level>0) { dbgs._stream << t; }
    return dbgs;
  }

} // namespace Ariadne

#endif // _ARIADNE_DEBUG_H
