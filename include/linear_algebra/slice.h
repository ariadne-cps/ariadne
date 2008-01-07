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
 
/*! \file vector.h
 *  \brief Vector types and vector operations.
  */

#ifndef ARIADNE_SLICE_H
#define ARIADNE_SLICE_H 

namespace Ariadne {
  namespace LinearAlgebra {
    
    class Slice 
    {
     public:
      Slice(size_type start, size_type size, size_type stride=1u)
        : _start(start), _size(size), _stride(stride) { }
      size_type start() const { return this->_start; }
      size_type size() const { return this->_size; }
      size_type stride() const { return this->_stride; }
     private:
      size_type _start;
      size_type _size;
      size_type _stride;
    };

    inline Slice slice(size_type start, size_type size, size_type stride=1u) {
      return Slice(start,size,stride); }

    inline Slice range(size_type start, size_type finish, size_type stride=1u) {
      return Slice(start,(finish-start-1)/stride+1,stride); }

  }
}


#endif /* ARIADNE_SLICE_H */
