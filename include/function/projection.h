/***************************************************************************
 *            projection.h
 *
 *  Copyright  2008  Pieter Collins
 *
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
 
/*! \file projection.h
 *  \brief Projections.
  */

#ifndef ARIADNE_PROJECTION_H
#define ARIADNE_PROJECTION_H 

#include "base/types.h"

namespace Ariadne {
  
    
    class Projection 
    {
     public:
      Projection(size_type result_size, size_type argument_size, const size_type* coordinates)
        : _argument_size(argument_size), _coordinates(coordinates,coordinates+result_size) { }
      static Projection slice(size_type result_size, size_type argument_size,
                              size_type start, size_type stride=1u) {
        ARIADNE_ASSERT(start+(result_size-1)*stride<argument_size);
        Projection result(result_size, size_type argument_size);
        for(uint i=0; i!=result_size; ++i) { 
          result._coordinates[i]=start+i*stride; }
        return result;
      }
      static Projection range(size_type result_size, size_type argument_size,
                              size_type start) {
        return slice(result_size,argument_size,start,1u); }
      size_type result_size() const { return this->_result_size; }
      size_type argument_size() const { return this->_coordinates.size(); }
      const size_type& operator[](size_type i) const { return this->_coordinates[i]; }
     private:
      Projection(size_type result_size, size_type argument_size)
        : _argument_size(argument_size), _coordinates(result_size) { }
     private:
      size_type _argument_size;
      array<size_type> _coordinates;
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
