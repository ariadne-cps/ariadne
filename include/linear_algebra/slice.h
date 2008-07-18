/***************************************************************************
 *            slice.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file slice.h
 *  \brief Slicing.
  */

#ifndef ARIADNE_SLICE_H
#define ARIADNE_SLICE_H 

#include "base/types.h"

namespace Ariadne {
  
    
    class Slice 
    {
     public:
      Slice(size_type start, size_type size, size_type stride=1u)
        : _start(start), _size(size), _stride(stride) { }
      size_type start() const { return this->_start; }
      size_type stop() const { return this->_start+this->_size*this->_stride; }
      size_type size() const { return this->_size; }
      size_type stride() const { return this->_stride; }
     private:
      size_type _start;
      size_type _size;
      size_type _stride;
    };

    inline Slice slice(size_type start, size_type size, size_type stride=1u) {
      return Slice(start,size,stride); }

    inline Slice range(size_type start, size_type stop) {
      return Slice(start,stop-start,1u); }

    class Range 
    {
     public:
      Range(size_type start, size_type stop)
        : _start(start), _size(stop-start) { }
      size_type start() const { return this->_start; }
      size_type stop() const { return this->_start+this->_size; }
      size_type size() const { return this->_size; }
      size_type stride() const { return 1u; }
     private:
      size_type _start;
      size_type _size;
    };



} // namespace Ariadne


#endif /* ARIADNE_SLICE_H */
